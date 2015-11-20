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

package EnsEMBL::REST::Model::ga4gh::annotationSets;

use Moose;
extends 'Catalyst::Model';
use Data::Dumper;
use Bio::EnsEMBL::Variation::Utils::VEP qw/read_cache_info get_version_data /;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');

our $species = 'homo_sapiens';

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

## TO DO
##  handle dataset
##  previous releases as annotation sets

## take version info from database or VEP cache?
##  - start with db
##  - access to historic versions needed

sub fetch_annotationSet {
  
  my $self   = shift;

  my $c = $self->context();
  
  my $core_ad = $c->model('Registry')->get_DBAdaptor($species, 'Core',    );
  my $var_ad  = $c->model('Registry')->get_DBAdaptor($species, 'variation');


#  my $var_meta = $var_ad->get_MetaContainer();
#  my $version = $var_meta->schema_version();

  my $annotationSet;

  ## variation required sources
  my $source_ad = $c->model('Registry')->get_adaptor($species, 'variation', 'Source');
  foreach my $name (qw[ dbSNP ClinVar]){
    my $source   = $source_ad->fetch_by_name($name);
    next unless defined $source;
    $annotationSet->{analysis}->{info}->{ $name . "_version" } = $source->version(); 
  }


  ## variation required meta
  my $meta_ext_sth = $var_ad->dbc->db_handle->prepare(qq[ select meta_key, meta_value from meta]);
  $meta_ext_sth->execute();
  my $stuff = $meta_ext_sth->fetchall_arrayref();

  my %meta;
  foreach my $l(@{$stuff}){
    $meta{$l->[0]} = $l->[1] if defined $l->[1];
  }

  $annotationSet->{analysis} =  { 'name'        => 'Ensembl',
                                  'created'    =>   $meta{"tv.timestamp"}
                                  };

  ## Incl. Var Set name
  $annotationSet->{id} = 'Ensembl_' . $meta{schema_version};
  $annotationSet->{analysis}->{info}->{Ensembl_version}  = $meta{schema_version};

  foreach my $v_attrib (qw [polyphen_version sift_version sift_protein_db 1000genomes_version ]){
     $annotationSet->{analysis}->{info}->{$v_attrib} = $meta{$v_attrib} if defined  $meta{$v_attrib};
  }


  ## core required meta
  my $cmeta_ext_sth = $core_ad->dbc->db_handle->prepare(qq[ select meta_key, meta_value from meta]);
  $cmeta_ext_sth->execute();
  my $core_meta = $cmeta_ext_sth->fetchall_arrayref();

  my %cmeta;
  foreach my $l(@{$core_meta}){
    $cmeta{$l->[0]} = $l->[1];
  }

  foreach my $c_attrib (qw[genebuild.id assembly.name assembly.accession gencode.version assembly.long_name genebuild.last_geneset_update genebuild.havana_datafreeze_date]){
     $annotationSet->{analysis}->{info}->{$c_attrib} = $cmeta{$c_attrib} if defined  $cmeta{$c_attrib};
  }

  return { annotationSets => [$annotationSet]}; 

}

## not useful without archive data
sub getAnnotationSet{

  my ($self, $id ) = @_; 

  my $c = $self->context();
  
  my $var_ad  = $c->model('Registry')->get_DBAdaptor($species, 'variation');
  my $var_meta = $var_ad->get_MetaContainer();
  my $version = $var_meta->schema_version();

  my $current = "Ensembl_" . $version; 

  ## exit if not current
  $self->context()->go( 'ReturnError', 'custom', [ " No data available for set $id" ] )
    unless $id =~/$current/i || $id eq 'Ensembl';

  return $self->fetch_annotationSet();

}


1;
