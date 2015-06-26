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

package EnsEMBL::REST::Model::ga4gh::referenceSets;

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

### using assembly name as id pending core support for sequence md5's 

sub fetch_referenceSet {
  
  my $self = shift;

  my $c = $self->context();
  
  my $post_data = $c->req->data;

  $c->log->debug(Dumper $post_data);

  $c->go( 'ReturnError', 'custom', [ ' Error - search by md5sum not currently supported'])
    if exists $post_data->{md5checksums} ;

  
  $c->go( 'ReturnError', 'custom', [ ' Error - search by accession not currently supported'])
    if exists $post_data->{accessions} ;

  return $self->getReferenceSet( $post_data->{assemblyId} );

}


sub getReferenceSet{

  my ($self,  $get_id ) = @_; 

  my $c = $self->context();

  my $core_ad = $c->model('Registry')->get_DBAdaptor($species, 'Core',    );
  my $cmeta_ext_sth = $core_ad->dbc->db_handle->prepare(qq[ select meta_key, meta_value from meta]);
  $cmeta_ext_sth->execute();
  my $core_meta = $cmeta_ext_sth->fetchall_arrayref();

  my %meta;
  foreach my $l(@{$core_meta}){
    $meta{$l->[0]} = $l->[1];
  }

  ## exit if not current
  $self->context()->go( 'ReturnError', 'custom', [ " No data available for this reference set: $get_id" ] )
    if defined $get_id &&  $get_id !~/$meta{"assembly.name"}/i;

  my $referenceSet;
  $referenceSet->{id}           = $meta{"assembly.name"};
  $referenceSet->{referenceIds} = [];
  $referenceSet->{md5checksum}  = 'md5';
  $referenceSet->{ncbiTaxonId}  = $meta{"species.alias"};
  $referenceSet->{description}  = $meta{"assembly.longname"};
  $referenceSet->{assemblyId}   = $meta{"assembly.name"};
  $referenceSet->{sourceURI}    = 'ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/dna/';  ##FIX!
  $referenceSet->{sourceAccessions} =[ $meta{"assembly.accession"}];
  $referenceSet->{isDerived}    = 'false';


  return { referenceSets => [$referenceSet]};

}


1;
