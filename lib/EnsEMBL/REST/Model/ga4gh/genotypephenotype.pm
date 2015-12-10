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
use Data::Dumper;

with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');

our $species = 'homo_sapiens';

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });

}

## can specify by feature and/ or phenotype and/or evidence

sub fetch_g2p_results {

  my ($self, $data ) = @_;

  my $results;

  if( defined $data->{feature}->[0]){
    $results = $self->fetch_by_feature($data);
  }
  elsif( defined $data->{phenotype}->[0]){
    $results = $self->fetch_by_phenotype($data);
  } 
  elsif( defined $data->{evidence}->[0]){
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

  foreach my $feature ( @{$data->{feature}} ){
    my @pfs = @{$pfa->fetch_all_by_object_id( $feature )};

    foreach my $pf (@pfs){
      #next if defined $data->{phenotype}->[0] && $pf->phenotype()->description() !~ /$data->{phenotype}->[0]/;
      push @$pheno_feat_out, $pf;
    }
  }

  return $pheno_feat_out;
}


sub fetch_by_phenotype{

  my $self = shift;
  my $data = shift;

  ## temp file for description <-> URI link pending db update

  ## many descriptions => one id
  my %required;
  foreach my $uri ( @{$data->{phenotype}} ){
    $required{$uri} = 1;
  }

  my %uri;  ##content: uri\tdesc
  open my $temp_lookup, "/home/vagrant/src/ensembl-rest/gwas_uri_desc.txt " || die "Failed to file gwas_uri_desc.txt  for pheno URL :$!\n"; 
  while(<$temp_lookup>){
    chomp;
    my @a = split/\t/;
    $uri{$a[1]} = $a[0] if $required{$a[0]} ; 
  }

  ## to database with descriptions
  my $pfa = $self->context()->model('Registry')->get_adaptor($species, 'variation', 'PhenotypeFeature');

  my $pheno_feat_out;

  my $count=0; ## only pass enough for page size - could be many

  foreach my $desc (sort keys %uri ){
    warn "Looking for pheno : $desc\n";
    my $pfs = $pfa->fetch_all_by_phenotype_description_source_name( $desc );
  
   foreach my $pf (@{$pfs}){
       last if defined $data->{pageSize} && $data->{pageSize}=~/\d+/ && $count== $data->{pageSize};

      ## TODO filter by evidence

      push @{$pheno_feat_out}, $pf;
      $count++;
    }
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


  foreach my $pmid (@{$data->{evidence}}){

    my $studies = $sta->fetch_all_by_external_reference( $pmid );
print Dumper $studies;
    unless ( defined $studies->[0] ){
      warn "No study found for $pmid\n";
      next;
    }

    foreach my $study (@{$studies}){
      my $pfs = $pfa->fetch_all_by_Study( $study);
      foreach my $pf (@{$pfs}){

        push @{$pheno_feat_out}, $pf;
      }
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
  my $ontol_info = read_ontology_file();

  my @assocs;

  my $count;
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

    push @{$assoc->{evidence}}, { evidenceType => { ontologySourceName  => 'IAO',
                                                    ontologySourceID    => "http://purl.obolibrary.org/obo/IAO_0000311",
                                                    ontologySourceVersion => undef},
                                  description  => $ext_ref 
                                }    if defined $ext_ref && $ext_ref =~/PMID/;

    push @{$assoc->{evidence}}, { evidenceType => { ontologySourceName  => 'OBI',
                                                    ontologySourceID    => "http://purl.obolibrary.org/obo/OBI_0001442",
                                                    ontologySourceVersion => undef},
                                  description  => $pf->p_value()
                                }    if defined $pf->p_value(); 


    push @{$assoc->{evidence}}, { evidenceType => { ontologySourceName  => 'OBCS',
                                                    ontologySourceID    => "http://purl.obolibrary.org/obo/OBCS_0000085",
                                                    ontologySourceVersion => undef},
                                  description  => $pf->beta_coefficient()
                                }    if defined $pf->beta_coefficient();

    push @{$assoc->{evidence}}, { evidenceType => { ontologySourceName  => 'OBCS',
                                                    ontologySourceID    => "http://purl.obolibrary.org/obo/OBCS_0000054",
                                                    ontologySourceVersion => undef},
                                  description  => $pf->odds_ratio()
                                }    if defined $pf->odds_ratio();

    push @{$assoc->{evidence}}, { evidenceType => { ontologySourceName  => undef,
                                                    ontologySourceID    => "? risk allele",
                                                    ontologySourceVersion => undef},
                                  description  => $pf->risk_allele()
                                 }   if defined $pf->risk_allele();



    ## PhenotypeInstance record
    my $phenotype_id = $pf->phenotype()->dbID;
    $assoc->{phenotype}->{id} = $phenotype_id; ##temp
    ### look up from file for now ?? KEY off DB id??
    ## OntologyTerm record
    my $ensembl_pheno_desc = $pf->phenotype->description();
    my $lc_ensembl_pheno_desc = "\L$ensembl_pheno_desc";
    print "No ontol term for $ensembl_pheno_desc\n" unless defined $ontol_info->{$lc_ensembl_pheno_desc};
 
    my @a = split/\//, $ontol_info->{$lc_ensembl_pheno_desc};
    my $term = pop(@a);
    my $ontol_source = (split/\_/, $term)[0];

    $assoc->{phenotype}->{type}->{ontologySourceName}    = $ontol_source;
    $assoc->{phenotype}->{type}->{ontologySourceID}      = $ontol_info->{ $lc_ensembl_pheno_desc };
    $assoc->{phenotype}->{type}->{ontologySourceVersion} = "";

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


##temp hack
sub feature_type{

  my $type = shift;
  return {  source => "SO",
            name   => "sequence_variant",
            id     => "SO:0001060" } if $type =~/Variation/;
}


## temp - build lookup hash of ontology terms prior to database'ing
sub read_ontology_file{

  my %ontol_info;
  open my $lookup, "/home/vagrant/src/ensembl-rest/gwas_uri_desc.txt" || die "Failed to open temp_pheno_ontol.txt :$!\n";
  while(<$lookup>){
    chomp;

    my ( $ontol, $desc) = split/\t/; 
    $ontol_info{"\L$desc"} = $ontol;
  }

  return \%ontol_info;
}



