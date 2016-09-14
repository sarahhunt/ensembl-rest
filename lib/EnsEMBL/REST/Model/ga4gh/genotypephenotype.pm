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
##     IDs
##     timestamps  -  schema does not support

use Moose;
extends 'Catalyst::Model';
use Catalyst::Exception;
use Scalar::Util qw/weaken/;
use Data::Dumper;
use Digest::MD5 qw(md5_hex);
use EnsEMBL::REST::Model::ga4gh::ga4gh_utils;

with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro', weak_ref => 1);


## hard coded supported ources for now.
has 'supported_sets' => ( isa => 'HashRef', is => 'ro', lazy => 1, default => sub {
  return {
    map { md5_hex($_) => $_ } ("ClinVar", "Orphanet", "DDG2P", "NHGRI-EBI GWAS catalog")
  };
});



sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  weaken($c);
  return $self->new({ context => $c, %$self, @args });

}

## can specify by feature_id and/ or phenotype_id and/or evidence

sub fetch_g2p_results {

  my $self = shift;
  my $data = shift;

  my $supported_sets = $self->supported_sets();

  return { associations    => [],
           nextPageToken => undef} unless $supported_sets->{ $data->{phenotypeAssociationSetId} };

  $data->{source} = $supported_sets->{ $data->{phenotypeAssociationSetId} };

  my $results;

  if( $data->{featureIds} && defined $data->{featureIds}->[0] ){
    $results = $self->fetch_by_feature($data);
  }
  elsif( $data->{phenotypeIds} && defined $data->{phenotypeIds}->[0]){
    $results = $self->fetch_by_phenotype($data);
  }
  else{
    ## shouldn't get here - either feature or phenotype mandatory
    $self->context()->go( 'ReturnError', 'custom', [ "Nothing to search on "]);
  }

  my ($assoc, $nextPageToken) = $self->format_results($results, $data );

  return{ associations    => $assoc,
          nextPageToken => $nextPageToken};
}

## feature id is currently stable id eg ENSG00000176515.1 or prefixed extternal eg 85:rs6920220
sub fetch_by_feature{

  my $self = shift;
  my $data = shift; 

  my $required_phenos   = $self->list_to_hash($data->{phenotypeIds});
  my $required_evidence = $self->list_to_hash($data->{evidence});

  my $pfa = $self->context()->model('Registry')->get_adaptor('homo_sapiens', 'variation', 'PhenotypeFeature');

  my $pheno_feat_out;
  my @pfs;

  ## extract all phenotype features for the required features
  foreach my $feature_id(@{$data->{featureIds}}){
    my $stable_id = $feature_id;
    $stable_id =~ s/\.\d+$//;  ## for genes
    $stable_id =~ s/^\d+\://;  ## for variants

    push @pfs, @{$pfa->fetch_all_by_object_id( $stable_id )};
  }

  ## filter by phenotype & evidence
  foreach my $pf (@pfs){
    ## filter by phenotype
    next if  scalar(keys $required_phenos) > 0 &&  $required_phenos->{ $pf->phenotype()->dbID() };

## ************************** API update required- extract by source **************************
    next unless $pf->source_name() ne $data->{source};

    ## filter by evidence
    ## next if defined $required_evidence && ! $required_evidence{}
    push @$pheno_feat_out, $pf;

  }
  return $pheno_feat_out;
}


sub fetch_by_phenotype{

  my $self = shift;
  my $data = shift;


  my $pfa = $self->context()->model('Registry')->get_adaptor('homo_sapiens', 'Variation', 'PhenotypeFeature');

  my @pfs;
  foreach my $phenotype_id (@{$data->{phenotypeIds}}){
    push @pfs, @{$pfa->fetch_all_by_phenotype_id_source_name( $phenotype_id, $data->{source} ) };
  }

  my $pheno_feat_out;

  foreach my $pf (@pfs){

    ## TODO filter by evidence

    push @{$pheno_feat_out}, $pf;

  }

  return $pheno_feat_out;

}
=head

## This is now removed from schema - would be good to re-instate
## move to pheno feat attrib?? -  would work for ClinVar too
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

=cut


## create a set of FeaturePhenotypeAssociation's from a set of ensembl PhenotypeFeatures
sub format_results{

  my $self  = shift;
  my $pfs   = shift;
  my $data  = shift;

  my $nextPageToken;
  my @assocs = ();
  my $count  = 0;
  my %current_versions;

  my $start = (defined $data->{pageToken} ? 0 : 1);

  my $ont_ad    = $self->context->model('Registry')->get_adaptor('Multi', 'Ontology', 'OntologyTerm');
  my $gene_ad   = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'core', 'gene');

  my $release_version = $self->context->model('ga4gh::ga4gh_utils')->release_version();

  foreach my $pf( @{$pfs} ){  ## sort on seq_region + pos??

    ## paging
    $start = 1 if  $start == 0 && $pf->dbID() eq $data->{pageToken};
    next unless $start == 1;

    if (defined $data->{pageSize} &&  $data->{pageSize} eq $count){
      $nextPageToken = $pf->dbID();
      last;
    }

    ## feature ids used in GA4GH endpoints are more complex than those stored in phenotype feature
    unless ($current_versions{ $pf->object_id()} ){
      if ( $pf->type() eq 'Gene'){
        my $gene = $gene_ad->fetch_by_stable_id($pf->object_id() );
        $current_versions{ $pf->object_id()} = $gene->stable_id() .'.'. $gene->version();      
      }
      else{
        $current_versions{ $pf->object_id()} = $release_version .':'. $pf->object_id();
      }
    }
    my $assoc = { id                           => $pf->dbID(),
                  phenotypeAssociationSetId => $data->{phenotype_association_set_id},
                  featureIds                  => [ $current_versions{$pf->object_id()}],
                  environmentalContexts       => []
                }; 

    ($assoc->{evidence}, $assoc->{info}) = $self->format_attribs($pf);


    ## schema only suports one ontology term per phenotype
    ## add formatting to perl API when stable?
    my $ontology_accessions = $pf->phenotype()->ontology_accessions() ;
    my $ont_term            = $ont_ad->fetch_by_accession( $ontology_accessions->[0] );
    my $formatted_term = {};
    if (defined $ont_term ){
      $formatted_term = { sourceName    => $ont_term->ontology(),
                          id            => $ont_term->accession(),
                          term          => $ont_term->name(),
                          sourceVersion => $ont_term->ontology_version() 
                        };
    }

    ## PhenotypeInstance record
    $assoc->{phenotype} = { id           => $pf->phenotype()->dbID(),
                            type         => $formatted_term ,
                            qualifier    =>[],
                            ageOfOnset => undef,
                            description  =>  $pf->phenotype()->description(),
                            info         => {}
                          };

   ## TODO do this on extraction
   next if 

    push @assocs, $assoc;

    $count++;
 }

  return( \@assocs, $nextPageToken);
}


## convert attribs to evidence
sub format_attribs{

  my $self = shift;
  my $pf   = shift;

  ## PMID held here rather than as attrib:
  my $ext_ref = $pf->study->external_reference()
    if defined $pf->study() && defined $pf->study->external_reference();

  my $attribs = $pf->get_all_attributes();


  my @evidence = ();

  push @evidence,  $self->format_evidence( 'IAO',
                                           'http://purl.obolibrary.org/obo/IAO_0000311',
                                           'publication',
                                            $ext_ref )  if defined $ext_ref && $ext_ref =~/PMID/;

  push @evidence,  $self->format_evidence( 'OBI',
                                           'http://purl.obolibrary.org/obo/OBI_0000175',
                                           'p-value',
                                            $attribs->{p_value} )   if defined $attribs->{p_value};


  push @evidence,  $self->format_evidence( 'OBCS',
                                           'http://purl.obolibrary.org/obo/OBCS_0000085',
                                           'standardized coefficient',
                                            $attribs->{beta_coefficient} )   if defined $attribs->{beta_coefficient};

  push @evidence, $self->format_evidence( 'OBCS',
                                          'http://purl.obolibrary.org/obo/OBCS_0000054',
                                          'odds ratio',
                                           $attribs->{odds_ratio} )   if defined $attribs->{odds_ratio};


  my $info  = {};
  $info->{risk_allele}     = $attribs->{risk_allele}       if defined $attribs->{risk_allele} && $attribs->{risk_allele} =~ /\D+/;
  $info->{associated_gene} = $attribs->{associated_gene}   if defined $attribs->{associated_gene};

  ## for DDG2P
  $info->{inheritance_type} = $attribs->{inheritance_type} if defined $attribs->{inheritance_type};

  return (\@evidence, $info);

}
sub format_evidence{

  my @vals = @_;

  return { evidenceType => { sourceName  => $vals[1],
                             id           => $vals[2],
                             term         => $vals[3],
                             sourceVersion => undef },
           description  => $vals[4]
         };   
}

sub list_to_hash{

  my $self = shift;
  my $data = shift;

  return { map {$_ => 1 } @{$data} };

}


1;
