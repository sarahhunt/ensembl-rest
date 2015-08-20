=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package EnsEMBL::REST::Model::ga4gh::variants;

use Moose;
extends 'Catalyst::Model';

use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;

use Data::Dumper;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');

our $DEBUG = 0;

=head2

Now mapping ensembl source        => GA4GH dataset
                    VCFcollection => variationSet

08/2015 Master version

=cut

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub fetch_gavariant {

  my ($self, $data ) = @_; 

  ## get the VCF collection object for the required set
  $data->{vcf_collection} = $self->get_VCFcollection($data->{variantSetId});

  ## format sample names if filtering by sample required
  $data->{req_samples} = $self->check_sample_info($data->{callSetIds}) if defined $data->{callSetIds}->[0] ;

  ## get genotype data - pull out set by region & filter on variant name if supplied
  return $self->fetch_by_region($data);

}

## get VCFCollections object for required VariantSet
sub get_VCFcollection{

  my ($self, $variantSetId ) = @_;


  $Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $self->{ga_config};
  my $vca = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new();

 
  my $vcf_coll = $vca->fetch_by_id($variantSetId);

  $self->context()->go( 'ReturnError', 'custom', [ " Failed to find the specified variantSetId"])
    unless defined $vcf_coll; 

  return $vcf_coll;
}
  

## format sample names if filtering by sample required
sub check_sample_info{

  my ($self, $callSetIds ) = @_;

  my %req_samples; 

  foreach my $sample ( @{$callSetIds} ){
    $req_samples{$sample} = 1;
  }

  return \%req_samples;
}


sub fetch_by_region{

  my $self = shift;
  my $data = shift;

  ## if this is a new request set first token to start of region 
  $data->{pageToken} = $data->{start} unless exists  $data->{pageToken} && $data->{pageToken} =~/\d+/;


  ## get the next range of variation data for the token - should return slightly more than required
  my ($var_info, $next_token)  = $self->get_next_by_token($data);


  warn "No variants are available for this region\n"   
     unless defined $var_info && scalar(@{$var_info}) >0 && $DEBUG ==1;


  return ({ "variants"      => $var_info,
            "nextPageToken" => $next_token
          });
  
}


## extract genotypes & apply filtering
## input: array of strings - format : 'NA10000:0|1:44:23'
## output: array of individual genotype hashes
sub sort_genotypes {
  my ($self, $parser, $data, $is_remapped) = @_;


  my $geno_strings  = $parser->get_samples_info();

  my @genotypes;

  foreach my $sample(keys %{$geno_strings}){

    ## filter by individual if required
    next if defined $data->{callSetIds}->[0] && 
      !defined $data->{req_samples}->{$sample};

    my $gen_hash;
    $gen_hash->{callSetId}    = $sample;
    $gen_hash->{callSetName}  = $sample;

    my @g = split/\||\//, $geno_strings->{$sample}->{GT};
    foreach my $g(@g){
      ## force genotype to be numeric
      push @{$gen_hash->{genotype}}, numeric($g);
    } 

    ## place holders
    $gen_hash->{phaseset}           = '';
    if( $is_remapped ){
      ## minimal data relevant if remapped
      $gen_hash->{genotypeLikelihood} = [];
      $gen_hash->{info}               = {};
    }
    else{

      foreach my $inf (keys %{$geno_strings->{$sample}} ){
        next if $inf eq 'GT';
        if( $inf eq 'GL'){
          my @gl = split/\,/, $geno_strings->{$sample}->{GL};
          foreach my $gl (@gl){
            push @{$gen_hash->{genotypeLikelihood}}, numeric($gl);
          }
        }
        else{
          if( $geno_strings->{$sample}->{$inf} =~/\,/){
            @{$gen_hash->{info}->{$inf}} = split/\,/, $geno_strings->{$sample}->{$inf};
          }
          else{ 
            $gen_hash->{info}->{$inf} = [$geno_strings->{$sample}->{$inf}];
          }
        }
      }
    }

    push @genotypes, $gen_hash;

  }
  return \@genotypes;
}



## extract a batch of results for request and hold new start pos in token string
sub get_next_by_token{

  my ($self, $data) = @_;


  ## look up filename from vcf collection object
  my $file  =  $data->{vcf_collection}->filename_template(); 
  $file =~ s/\#\#\#CHR\#\#\#/$data->{referenceName}/;
  $file = $self->{geno_dir} .'/'. $file;

 # return these ordered by position for simple pagination
  my @var;
  my $nextToken;

  ## exits here if unsupported chromosome requested
  return (\@var, $nextToken) unless -e $file;

  my $parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open( $file ) || die "Failed to get parser : $!\n";
  $parser->seek($data->{referenceName}, $data->{pageToken}, $data->{end});


  my $n = 0;

  while($n < ($data->{pageSize} + 1)){
    
    my $got_something = $parser->next();
    last if $got_something ==0;

    my $name = $parser->get_IDs->[0];
#    next unless defined $name ;


   if ($n == $data->{pageSize} ){
      ## batch complete 
      ## save next position for new page token
      $nextToken = $parser->get_raw_start;
      last;
    }


    ## add filter for variant name if required
    next if defined $data->{variantName} && $data->{variantName} =~/\w+/ &&  $data->{variantName} ne $name;

    next if $name=~ /esv/; ##skip these for now

    ## format array of genotypes
    my $genotype_calls = $self->sort_genotypes($parser, $data, $data->{vcf_collection}->{is_remapped});


    ## check there are genotypes to return
    next unless defined $genotype_calls;

    $n++;
    my $variation_hash;

    $variation_hash->{variantSetId}    = $data->{variantSetId};
    $variation_hash->{calls}           = $genotype_calls;

    $variation_hash->{names}           = [ $name ];
    $variation_hash->{id}              = $name;
    $variation_hash->{referenceBases}  = $parser->get_reference;
    $variation_hash->{alternateBases}  = \@{$parser->get_alternatives};
    $variation_hash->{referenceName}   = $parser->get_seqname ;

    ## position is zero-based + closed start of interval 
    $variation_hash->{start}           = numeric($parser->get_raw_start - 1);
    ## open end of interval
    $variation_hash->{end}             = numeric($parser->get_raw_end);

    $variation_hash->{created}         = $data->{vcf_collection}->created || undef;
    $variation_hash->{updated}         = $data->{vcf_collection}->updated || undef; 


    ## What can be trusted if variants re-mapped but not re-called? Start with: AC, AF, AN
    my $var_info = $parser->get_info();
    if( $data->{vcf_collection}->{is_remapped} ){
      $variation_hash->{info} = { AC => [$var_info->{AC}],
                                  AN => [$var_info->{AN}],
                                  AF => [$var_info->{AF}]
                                };
    }
    else{
      foreach my $k (keys %{$var_info}){
        $variation_hash->{info}->{$k}  =  [$var_info->{$k}];
      }
    }

     push @var, $variation_hash;

   
    ## exit if single required variant already found
    last if defined $data->{variantName} && $data->{variantName} =~/\w+/;
  }


  $parser->close();

  return (\@var, $nextToken) ;

}


=head2 getVariant

  Gets a Variant by ID.
  GET /variants/{id} will return a JSON version of Variant.

=cut

sub getVariant{

  my ($self, $id ) = @_; 

  my $c = $self->context();
  my $species = "homo_sapiens"; 
  my $varfeat;
 
  my $va  = $c->model('Registry')->get_adaptor($species, 'Variation', 'Variation');
  my $vfa = $c->model('Registry')->get_adaptor($species, 'Variation', 'VariationFeature');

  $vfa->db->include_failed_variations(0); ## don't extract multi-mapping variants
 
  my $var = $va->fetch_by_name($id);
  my $vf  = $vfa->fetch_all_by_Variation($var) if defined $var;  

  if (defined $vf->[0]) {
    $varfeat = $vf->[0];
    ## retrun 1KG data by default if available 
    my $var_info = $self->getSingleCallSets($vf->[0], $id);
 
    return ({ "variants"      => [$var_info]}) if exists $var_info->[0]->{name} ;   
  }
  elsif($id =~/\w+\:c\.\w+|\w+\:g\.\w+|\w+\:p\.\w+/ ){
    ## try to look up as HGVS
    eval { $varfeat = $vfa->fetch_by_hgvs_notation( $id ) };
  }
  
  $c->go( 'ReturnError', 'custom', [ " No variants are available with this id"])
     unless defined $varfeat;  

  ## return basic location info if available

   my $variation_hash;

   my @als = split/\//, $varfeat->allele_string();
   shift @als; ## remove reference allele
   $variation_hash->{name}            = $varfeat->variation_name();
   $variation_hash->{id}              = $id;
   $variation_hash->{referenceBases}  = $varfeat->ref_allele_string();
   $variation_hash->{alternateBases}  = \@als;
   $variation_hash->{referenceName}   = $varfeat->seq_region_name();

   ## position is zero-based + closed start of interval 
   $variation_hash->{start}           = $varfeat->seq_region_start() -1;
   ## open end of interval
   $variation_hash->{end}             = $varfeat->seq_region_end();
   $variation_hash->{created}         = 'null';
   $variation_hash->{updated}         = 'null';

   return $variation_hash;
}

=head2 getSingleCallSets

look up a default set of genotypes if queried by id

=cut
sub getSingleCallSets{

 my ($self, $varfeat, $id ) = @_;

  my $data;
  $data->{referenceName} = $varfeat->seq_region_name();
  $data->{start}         = $varfeat->seq_region_start() -1;
  $data->{end}           = $varfeat->seq_region_end();
  $data->{variantsetId}  = 1; ## return 1KG by default
  $data->{variantName}   = $id; ## check for supplied name not current database mane


  ## load VCFcollections object for variantSet 
  $data->{vcf_collection} = $self->get_VCFcollection($data->{variantsetId});

  ## create fake token -what should really be returned for get??
  $data->{pageSize} = 1;
  $data->{pageToken} = $data->{start};

  my ($var_info, $next_ds) = $self->get_next_by_token($data);

  ## exit if none found
  $self->context()->go( 'ReturnError', 'custom', [ " No variants are available for this region" ] )
     unless defined $var_info ;

  return $var_info;

}


=head
sub getGenotypesForSingleVariant{

  my ($self, $data, $vf ) = @_;

  my @datasetsreq = sort(keys %{$data->{files}});

  foreach my $ds(@datasetsreq){
  
    my $file = $self->{dir} .'/'. $data->{files}->{$current_ds};

    my $parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open( $file ) || die "Failed to get parser : $!\n";
    $parser->seek($vf->seq_region_name() ,$vf->seq_region_start(), $vf->seq_region_end());

=cut



sub numeric{
  my $string = shift;
  return $string * 1;
}


with 'EnsEMBL::REST::Role::Content';

__PACKAGE__->meta->make_immutable;

1;
