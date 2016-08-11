=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute 
Copyright [2016] EMBL-European Bioinformatics Institute
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

package EnsEMBL::REST::Model::ga4gh::phenotypes;

use Moose;
extends 'Catalyst::Model';
use Catalyst::Exception;
use Digest::MD5 qw(md5_hex);
use Scalar::Util qw/weaken/;
use Data::Dumper;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro', weak_ref => 1);

## hard coded sources for now.
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

## POST entry point
sub searchPhenotypes {
  
  my $self   = shift;
  my $data   = shift;


  ## check the set is supported
  my $supported_sets = $self->supported_sets();
  return { phenotypeAssociationSets   => [],
           nextPageToken => undef
         } unless $supported_sets->{ $data->{phenotype_association_set_id} };

  ## short term fix
  my $statement = (qq[ select distinct p.phenotype_id, p.description, poa.accession
                       from phenotype p, phenotype_ontology_accession poa, phenotype_feature pf, source s
                       where s.name =?
                       and pf.source_id = s.source_id
                       and p.phenotype_id = pf.phenotype_id
                       and p.phenotype_id = poa.phenotype_id
                       and poa.accession not like 'HP%'
                      ]);

  my $vardb_ad  = $self->context->model('Registry')->get_DBAdaptor('homo_sapiens', 'variation');
  my $ont_ad    = $self->context->model('Registry')->get_adaptor('Multi', 'Ontology', 'OntologyTerm');

  if( $data->{description}){
    ## look up ontology term by name or synonym with like '% %'
    my $ont_terms = $ont_ad->fetch_all_by_name( $data->{description} );

    return { phenotypeAssociationSets => [],
             nextPageToken            => undef
           } unless defined $ont_terms ;

    $statement .= " and poa.accession in (";
    foreach my $ont_term (@{$ont_terms }){
      $statement .= "'" . $ont_term->accession() ."', ";
    }
    $statement =~ s/,\s+$/ ) /;
  }  


  ## lots of optional query attributes
  $statement .= " and p.phenotype_id = $data->{id} "                      if $data->{id};
  $statement .= " and poa.accession = $data->{type} "                     if $data->{type} ;
  $statement .= " limit 100"; ## temp


  ## extract required meta data from core db
  my $pheno_ext_sth = $vardb_ad->dbc->db_handle->prepare( $statement) ||die;
  $pheno_ext_sth->execute( $supported_sets->{ $data->{phenotype_association_set_id}})   ||die;

  my $phenotypes = $pheno_ext_sth->fetchall_arrayref();

  my $n = 0;
  my $nextPageToken;

  my @phenos;
  my %done;  ## only supporting one OntologyTerm

  foreach my $l(@{$phenotypes}){

    next if $done{$l->[0]};
    $done{$l->[0]} = 1;

    my $ont   = $ont_ad->fetch_by_accession( $l->[2] ); 

    push @phenos, { id           => $l->[0],
                    type         => { id            => $ont->accession(),
                                      term          => $ont->name(), 
                                      sourceName    => $ont->ontology(),
                                      sourceVersion => undef 
                                    }, 
                        qualifier    => undef,
                        age_of_onset => undef,
                        description  => $l->[1],
                        info         => {}
                      };

    last if defined $data->{page_size} && $n ==  $data->{page_size} ; 
  }

  return { phenotypes     => \@phenos,
           nextPageToken  => $nextPageToken
         }; 
}

## GET entry point not actaully in the spec
sub getPhenotypeAssociationSet{

  my ($self, $id ) = @_; 

  my $data = { id       => $id,
               pageSize => 1
              };


  return ;

}



1;
