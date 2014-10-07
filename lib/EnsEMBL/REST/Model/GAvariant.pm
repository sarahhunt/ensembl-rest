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

package EnsEMBL::REST::Model::GAvariant;
## version using raw IO mods
use Moose;
extends 'Catalyst::Model';
use Data::Dumper;
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;


with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');
our $config_file = "/home/vagrant/src/ensembl-rest/ga_vcf_config.json"; 


sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub fetch_gavariant {

  my ($self, $data ) = @_; 

  ## load callSet to variantSet link 
  $data = $self->get_set_info($data);

  ## exclude samples outside specified variantSet and variantSet with no required samples
  $data = $self->check_sample_info($data) if defined $data->{callSetIds}->[0] ;

  ## set up files required
  $data->{files} = $self->get_req_file($data);

  ## get variation data - pull out set by region & filter on variant name if supplied
  return $self->fetch_by_region($data);

}

## link sample names to variantSet for later filtering
sub get_set_info{

  my ($self, $data ) = @_;

  ## extract required variantSets if supplied
  if(defined $data->{variantSetIds}->[0] && $data->{variantSetIds}->[0] >0){  
    foreach my $set (@{$data->{variantSetIds}}){
      $data->{required_set}->{$set} = 1;
    }
  }

  ## read config to get call set to variant set link
  open my $cf, $config_file ||
    $self->context()->go( 'ReturnError', 'custom', [ " Failed to find config to extract set ids variantSets"]);

  local $/ = undef;
  my $json_string = <$cf>;
  close $cf;

  my $config = JSON->new->decode($json_string) ||  
    $self->context()->go( 'ReturnError', 'custom', [ " Failed to parse config for variantSets"]); 

  foreach my $dataSet( @{$config->{collections}} ) {

    ## loop over callSets
    foreach my $callset_id(keys %{$dataSet->{individual_populations}} ){
 
     my $variantSet_id = $dataSet->{individual_populations}->{$callset_id}->[0];

      ## limit by variant set if required
      next if defined $data->{variantSetIds}->[0] && ! defined $data->{required_set}->{$variantSet_id} ;
 
      ## save sample to variantSet link
      $data->{sample2set}->{$callset_id} = $dataSet->{individual_populations}->{$callset_id}->[0];
    }
  }
  return $data;
}
  


## format sample names if filtering by sample required
## variantSet limitation takes presidence 
sub check_sample_info{

  my ($self, $data ) = @_;

  my %req_samples; ## store sample <-> variantSet link
  my %req_sets;    ## if samples are specified, don't look in sets not containing them

  foreach my $sample ( @{ $data->{callSetIds} } ){
    if (defined  $data->{sample2set}->{$sample}){ ## only set if set required or no set limit 
      $req_samples{$sample} = 1;
      $req_sets{ $data->{sample2set}->{$sample} } = 1;
    }
  }
 
  ## exit if the samples & sets are incompatible
  $self->context()->go( 'ReturnError', 'custom', [ " The specified callSets are not available in the specified variantSets"])
    unless scalar(keys %req_samples) >0;

  $data->{samples} = \%req_samples;

  ## reset variantSets to only those with samples
  $data->{required_set} = \%req_sets;
  
  return $data;
}

## place holder
sub get_req_file{

  my $self = shift;
  my $data = shift;

  my @files; 

  if(defined $data->{variantSetIds}->[0]){
    if($data->{variantSetIds}->[0] == 65){
      my $file =  "/home/vagrant/Genotypes/Illumina_platinum/NA12878_S1r.chr" . $data->{referenceName} . ".vcf.gz";
      push @files, $file;
    }
    else{
      my $file = "/home/vagrant/Genotypes/1KG_data/ALL.chr" . $data->{referenceName} . ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz";
      push @files, $file;
    }
  }
  else{
    my $file1 =  "/home/vagrant/Genotypes/Illumina_platinum/NA12878_S1r.chr" . $data->{referenceName} . ".vcf.gz";
    push @files, $file1;
    my $file2 = "/home/vagrant/Genotypes/1KG_data/ALL.chr" . $data->{referenceName} . ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz";
    push @files, $file2;
  } 

  return \@files;

}

sub fetch_by_region{

  my $self = shift;
  my $data = shift;

  my $c = $self->context();

  ## if this is a new request set first token to start of region and set 0  
  $data->{pageToken} = $data->{start} . "_0" unless exists  $data->{pageToken} && $data->{pageToken} =~/\d+/;

  ## get the next range of variation data for the token - should return slightly more than required
  my $var_info = $self->get_next_by_token($data);

  ## exit if none found
  $c->go( 'ReturnError', 'custom', [ " No variants are available for this search token: $data->{pageToken}, batch: $data->{pageSize}" ] )  
     unless scalar(@{$var_info}) >0;
  
  ## get the correct number of GAvariants - returning by variantSet means 1var => many GAvar
  my @var_response;
  my $gavariant_count = 0;

  my $next_token;
  ## last reported position & setid
  my ($last_pos, $last_set ) = split/\_/,  $data->{pageToken} ;

  foreach my $ga_var( @{$var_info}){
   
    last if $gavariant_count == $data->{pageSize};
     
    ## last batch may have ended mid-variant - skip any which have been returned
    next if $ga_var->{start} == $last_pos &&  $ga_var->{variantSetId} <= $last_set;

    ## save var and count total
    $gavariant_count++;
    push @var_response, $ga_var;
      
    if($gavariant_count == $data->{pageSize}){
      ## set next token to position and set of last reported result  
      $next_token =  $ga_var->{start} ."_". $ga_var->{variantSetId} ;
      last;
    }
  }
  
  print  localtime() . " responding\n";

  return ({ "variants"      => \@var_response,
            "nextPageToken" => $next_token
          });
  
}


## extract genotypes & apply filtering
## input: array of strings - format : 'NA10000:0|1:44:23'
## output: hash with key: variantSetID and value: array of individual genotype hashes
sub sort_genotypes {
  my ($self, $geno_strings, $data) = @_;

  my %genotypes;

  foreach my $geno_string(@{$geno_strings}){

    my ($sample, $call, $qual) = split(/\:/,$geno_string,3);

    ## filter by individual if required
    next if defined $data->{callSetIds}->[0] && 
      !defined $data->{samples}->{$sample};

    my $gen_hash;
    $gen_hash->{callSetId}    = $sample;
    $gen_hash->{callSetName}  = $sample;
    @{$gen_hash->{genotype}}  = split/\|/, $call;
   # @{$gen_hash->{genotypeLikelihood}}  = split/\|/, $qual if defined $qual;

    ## store genotypes by variationSetId
    push @{$genotypes{ $data->{sample2set}->{$sample}}}, $gen_hash
      if defined $data->{sample2set}->{$sample};
  }
  return \%genotypes;
}


=head

## extract all SO terms for the variants consequences & return string
sub Consequences{ 

  my $self = shift;
  my $vf   = shift;

  my %cons_list;

  my $consequences = $vf->get_all_OverlapConsequences();
  
  return undef unless defined $consequences->[0];

  foreach my $cons(@{$consequences}){
     $cons_list{ $cons->SO_term() } = 1;
   }
   my $cons_string = join(",", keys %cons_list);  

  return $cons_string;

}

=cut

## extract a batch of results for request and hold new start pos in token string
sub get_next_by_token{

  my ($self, $data) = @_;
 
  ## These are the last seq start and set reported
  my ($batch_start, $set_start) = (split/\_/,$data->{pageToken});

  ## get one extra to see if it is worth sending new page token (only one set may be requested)
  my $limit = $data->{pageSize} + 1;

  ## reduce look up region to a size large enough to hold the requested number of variants??
  my $batch_end =  $data->{end};

  my $parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($data->{files}->[0]) || die "Failed to get parser : $!\n";
  $parser->seek($data->{referenceName},$batch_start,$batch_end);

  ## return these ordered by position & set id to allow pagination
  my @var;

  my $n = 0;


  while($n < $limit){
    
    $parser->next();
    my $name = $parser->get_IDs->[0];

    ## add filter for variant name if require
    next if defined $data->{variantName} && $data->{variantName} =~/\w+/ &&  $data->{variantName} ne $name;

    my $raw_genotypes  = $parser->get_raw_individuals_info();
    ## extract arrays of genotypes by variantSet
    my $genotype_calls = $self->sort_genotypes($raw_genotypes, $data);

    my @sets_to_return; ## sort here for paging
    if(defined $data->{variantSetIds}->[0]){
      @sets_to_return = sort @{$data->{variantSetIds}};
    }
    else{
      @sets_to_return = sort (keys %{$genotype_calls});
    }
    ## loop over sets to divide up genotypes
    foreach my $set_required (@sets_to_return){

      ## check there are genotypes to return
      next unless exists $genotype_calls->{$set_required}->[0];

      $n++;
      my $variation_hash;

      $variation_hash->{variantSetId}    = $set_required;
      $variation_hash->{calls}           = $genotype_calls->{$set_required};

      $variation_hash->{name}            = $name;
      $variation_hash->{id}              = $name;
      $variation_hash->{referenceBases}  = $parser->get_reference;
      $variation_hash->{alternateBases}  = \@{$parser->get_alternatives};
      $variation_hash->{referenceName}   = $parser->get_seqname ;

      ## position is zero-based + closed start of interval 
      $variation_hash->{start}           = $parser->get_raw_start -1;
      ## open end of interval
      $variation_hash->{end}             = $parser->get_raw_end;

      push @var, $variation_hash;
    }
   
    ## exit if single required variant already found
    last if defined $data->{variantName} && $data->{variantName} =~/\w+/;
  }

  ## this should not happen
  $self->context()->go('ReturnError', 'custom', ["No data found in the required region"]) if $n ==0;
  
  return \@var;

}



with 'EnsEMBL::REST::Role::Content';

__PACKAGE__->meta->make_immutable;

1;
