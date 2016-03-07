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

package EnsEMBL::REST::Model::ga4gh::genotypephenotype;

## TO DO
##     Paging
##     Ontology storage/ look ups
##     IDs
##     sets
##     timestamps

use Moose;
extends 'Catalyst::Model';
use Catalyst::Exception;
use Scalar::Util qw/weaken/;
use Data::Dumper;

with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');

our $species = 'homo_sapiens';

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  weaken($c);
  return $self->new({ context => $c, %$self, @args });

}

## can specify by feature and/ or phenotype and/or evidence

sub fetch_g2p_results {

  my ($self, $data ) = @_;

  my $results;

  if( defined $data->{feature} && $data->{feature} ne ''){
    $results = $self->fetch_by_feature($data);
  }
  elsif( defined $data->{phenotype} && $data->{phenotype} ne ''){
    $results = $self->fetch_by_phenotype($data);
  } 
  elsif( defined $data->{evidence} && $data->{evidence} ne ''){
    $results = $self->fetch_by_evidence($data);
  } 
  else{
    ## shouldn't get here
    $self->context()->go( 'ReturnError', 'custom', [ "Nothing to search on "]);
  }

  my ($assoc, $nextPageToken) = $self->format_results($results, $data->{pageSize} );

  return{ associations  => $assoc,
          nextPageToken => $nextPageToken};
}


sub fetch_by_feature{

  my $self = shift;
  my $data = shift; 

  my $pfa = $self->context()->model('Registry')->get_adaptor($species, 'variation', 'PhenotypeFeature');

  my $pheno_feat_out;

  my @pfs = @{$pfa->fetch_all_by_object_id( $data->{feature} )};

  foreach my $pf (@pfs){
    #next if defined $data->{phenotype}->[0] && $pf->phenotype()->description() !~ /$data->{phenotype}->[0]/;
    push @$pheno_feat_out, $pf;
  }

  return $pheno_feat_out;
}


sub fetch_by_phenotype{

  my $self = shift;
  my $data = shift;

  ## read temp file for description <-> URI link pending db update
  my $terms = read_ontology_file("ontology");
  my $desc  = $terms->{$data->{phenotype}} ;
  Catalyst::Exception->throw("No pheno for $data->{phenotype}") unless defined $terms->{$data->{phenotype}} ;

  ## to database with descriptions
  my $pfa = $self->context()->model('Registry')->get_adaptor($species, 'variation', 'PhenotypeFeature');

  my $pheno_feat_out;

  my $count=0; ## only pass enough for page size - could be many

  my $pfs = $pfa->fetch_all_by_phenotype_description_source_name( $desc, 'NHGRI-EBI GWAS catalog' );

  foreach my $pf (@{$pfs}){
    last if defined $data->{pageSize} && $data->{pageSize}=~/\d+/ && $count== $data->{pageSize};

    ## TODO filter by evidence

    push @{$pheno_feat_out}, $pf;
    $count++;
  }

  return $pheno_feat_out;

}

## currently only recommending PMID look up 
## move to pheno feat attrib?? -  would work for ClinVar too
## VERY SLOW
sub fetch_by_evidence{

  my $self = shift;
  my $data = shift;

  my $pfa = $self->context()->model('Registry')->get_adaptor($species, 'variation', 'PhenotypeFeature');
  my $sta = $self->context()->model('Registry')->get_adaptor($species, 'variation', 'Study');

  my $pheno_feat_out;

  my $studies = $sta->fetch_all_by_external_reference( $data->{evidence} );

  unless ( defined $studies->[0] ){
    warn "No study found for $data->{evidence}\n";

  }

  foreach my $study (@{$studies}){
    my $pfs = $pfa->fetch_all_by_Study( $study);
    foreach my $pf (@{$pfs}){

      push @{$pheno_feat_out}, $pf;
    }
  }

  return $pheno_feat_out;

}

## create a set of FeaturePhenotypeAssociation's from a set of ensembl PhenotypeFeatures
sub format_results{

  my $self     = shift;
  my $pfs      = shift;
  my $pageSize = shift;

  my $nextPageToken;

  ## hack 
  my $ontol_info = read_ontology_file('description');

  my @assocs;

  my $count =0;

  foreach my $pf( @{$pfs} ){

    if (defined $pageSize &&  $pageSize eq $count){
      $nextPageToken = $pf->object_id();
      last;
    }
    my $assoc;
    $assoc->{id} = 'placeholder';

    my $feat = $self->format_feature( $pf );

    push @{$assoc->{features}}, $feat;

    ## Evidence - need ontology terms. How to store for look up??
    $assoc->{evidence} = [];

    ## PMID held here:
    my $ext_ref = $pf->study->external_reference() 
      if defined $pf->study() && defined $pf->study->external_reference(); 

    my $attribs = $pf->get_all_attributes();

    push @{$assoc->{evidence}}, { evidenceType => { sourceName    => 'IAO',
                                                    id            => "http://purl.obolibrary.org/obo/IAO_0000311",
                                                    term          => 'publication',
                                                    sourceVersion => undef},
                                  description  => $ext_ref 
                                }    if defined $ext_ref && $ext_ref =~/PMID/;

    push @{$assoc->{evidence}}, { evidenceType => { sourceName    => 'OBI',
                                                    id            => "http://purl.obolibrary.org/obo/OBI_0000175",
                                                    term          => 'p-value',
                                                    sourceVersion => undef},
                                  description  => $attribs->{p_value}
                                }    if defined $attribs->{p_value}; 


    push @{$assoc->{evidence}}, { evidenceType => { sourceName    => 'OBCS',
                                                    id            => "http://purl.obolibrary.org/obo/OBCS_0000085",
                                                    term          => 'standardized coefficient',
                                                    sourceVersion => undef},
                                  description  => $attribs->{beta_coefficient}
                                }    if defined $attribs->{beta_coefficient};

    push @{$assoc->{evidence}}, { evidenceType => { sourceName  => 'OBCS',
                                                    id          => "http://purl.obolibrary.org/obo/OBCS_0000054",
                                                    term        => 'odds ratio',
                                                    sourceVersion => undef},
                                  description  => $attribs->{odds_ratio}
                                }    if defined $attribs->{odds_ratio};


    ### THIS ISN'T IN THE SCHEMA
    $assoc->{info}  = {};
    $assoc->{info}->{risk_allele}     = $attribs->{risk_allele}     if defined $attribs->{risk_allele};
    $assoc->{info}->{associated_gene} = $attribs->{associated_gene} if defined $attribs->{associated_gene};

    ## for DDG2P
    $assoc->{info}->{inheritance_type} = $attribs->{inheritance_type}     if defined $attribs->{inheritance_type};    

    ## PhenotypeInstance record
    my $phenotype_id = $pf->phenotype()->dbID;
    $assoc->{phenotype}->{id} = $phenotype_id; ##temp
    ### look up from file for now ?? KEY off DB id??
    ## OntologyTerm record
    my $ensembl_pheno_desc = $pf->phenotype->description();
    my $lc_ensembl_pheno_desc = "\L$ensembl_pheno_desc";
    print "No ontol term for $ensembl_pheno_desc\n" unless defined $ontol_info->{$lc_ensembl_pheno_desc}->{id};
 
    my @a = split/\//, $ontol_info->{$lc_ensembl_pheno_desc}->{id};
    my $term = pop(@a);
    my $ontol_source = (split/\_/, $term)[0];

    $assoc->{phenotype}->{type}->{sourceName}    = $ontol_source;
    $assoc->{phenotype}->{type}->{id}            = $ontol_info->{ $lc_ensembl_pheno_desc }->{id};
    $assoc->{phenotype}->{type}->{term}          = $ontol_info->{ $lc_ensembl_pheno_desc }->{term};
    $assoc->{phenotype}->{type}->{sourceVersion} = "";

    ## hack to see how genes look
    if( $pf->source_name =~/Orphanet/){
      $assoc->{phenotype}->{type}->{sourceName}    = 'Orphanet';
      $assoc->{phenotype}->{type}->{id}            = 'http://www.orpha.net/ORDO/Orphanet_' . $attribs->{external_id};
      $assoc->{phenotype}->{type}->{term}          = $lc_ensembl_pheno_desc ;
    }

   ## bail if Phenotype Ontology term not available
   next unless defined $assoc->{phenotype}->{type}->{id}  ;

    $assoc->{phenotype}->{qualifier}   = '';
    $assoc->{phenotype}->{ageOfOnset}  = '';


    
    ## include source description - what makes sense here?
    $assoc->{description} = $pf->source_name ." ". $pf->phenotype()->description(); ;
    $assoc->{environmentalContexts} =  [ ];

    push @assocs, $assoc;

    $count++;
 }

  return \@assocs, $nextPageToken;
}

## move centrally ? ##
## in Ensembl API or EnsemblModel?
sub format_feature{

  my $self    = shift;
  my $pf      = shift;

  my $formatted;

  $formatted->{id}            = $pf->object_id();
  $formatted->{parentIds}     =  [];
  $formatted->{featureSetId}  =  $pf->source_name() . ":".  $pf->source_version() ;
  $formatted->{featureType}   =  feature_type($pf->type());  ## TODO ONTOL look up
  $formatted->{attributes}    = {} ;
  ## how to handle patches??
  $formatted->{referenceName} =  $pf->seq_region_name();
  $formatted->{start}         =  $pf->seq_region_start();
  $formatted->{end}           =  $pf->seq_region_end();

  return $formatted;
}


##temp hack - sort db lookup
sub feature_type{

  my $type = shift;

  return {  sourceName     => "Sequence Ontology",
            term           => "sequence_variant", 
            sourceVersion  => "",
            id             => "SO:0001060"
         } if $type  eq 'Variation';

  return {  sourceName     => "Sequence Ontology",
            term           => "gene",                  
            sourceVersion  => "",
            id             => "SO:0000704"
         } if $type eq 'Gene';

}

### FIX THIS
## temp - build lookup hash of desc terms prior to database'ing
sub read_ontology_file{

  my $type = shift;  ## what should be key?

  my %ontol_info;
  open my $lookup, "/home/vagrant/src/ensembl-rest/gwas_uri_term_desc.txt" || die "Failed to open temp_pheno_ontol.txt :$!\n";
  while(<$lookup>){
    chomp;

    my ($id, $term, $desc) = split/\t/; 
    if($type eq "ontology"){
      $ontol_info{$id} = $term;
    }
    else{
      $ontol_info{"\L$desc"}{id}   = $id;
      $ontol_info{"\L$desc"}{term} = $term;
    }
  }

  return \%ontol_info;
}



