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
our @ISA =('EnsEMBL::REST::Model::ga4gh');

use Moose;
extends 'Catalyst::Model';
use Data::Dumper;
use EnsEMBL::REST::Model::ga4gh::ga4gh_utils;

with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');
use EnsEMBL::REST::Model::ga4gh::ga4gh_utils;


sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

## POST entry point
sub searchReferenceSet {
  
  my $self = shift;

  #$c->log->debug(Dumper $elf->context()->req->data);

  my ( $referenceSets, $nextPageToken)  =  $self->fetchData( $self->context()->req->data );

  return ({ referenceSets => $referenceSets,
            nextPageToken => $nextPageToken });

}

## GET entry point
sub getReferenceSet {

  my $self = shift;
  my $id = shift;

  my $data = { id => $id};

  my ($referenceSets, $nextPageToken) =  $self->fetchData( $data );

  $self->context()->go( 'ReturnError', 'custom', ["ERROR: no data for ReferenceSet $id"])
    unless defined $referenceSets &&  ref($referenceSets) eq 'ARRAY' ;
 
  return ($referenceSets->[0]);

}


## send both post & get here as few sets to check
sub fetchData{

  my $self  = shift;
  my $data  = shift;

  my $c = $self->context();

  ## read config
  my $config = $self->context->model('ga4gh::ga4gh_utils')->read_sequence_config();
  my $referenceSets =  $config->{referenceSets};
  my $nextPageToken;

  ## return empty array if no sets available by this id (behaviour not fully specified)
  return ( [], $nextPageToken) unless defined $referenceSets &&  ref($referenceSets) eq 'ARRAY' ;


  my @referenceSets;

  $data->{pageToken} = 0 unless defined $data->{pageToken} && $data->{pageToken} ne ''; 

  my $count = 0;
  foreach( my $n = $data->{pageToken}; $n <  scalar @{$referenceSets}; $n++ ) {

    my $refset_hash = $referenceSets->[$n];

    ## filter if an attrib supplied
    next if defined $data->{id}          &&  $data->{id}          ne '' 
                                         &&  $data->{id}          ne $refset_hash->{id}; ##  GET

    next if defined $data->{md5checksum} &&  $data->{md5checksum} ne ''
                                         &&  $data->{md5checksum} ne $refset_hash->{md5};

    next if defined $data->{accession}   &&  $data->{accession}   ne ''
                                         &&  $data->{accession}   ne $refset_hash->{sourceAccessions}->[0]; ##FIX for other accessions

    next if defined $data->{assemblyId}  &&  $data->{assemblyId}  ne ''
                                         &&  $data->{assemblyId}  ne $refset_hash->{id};


    ## paging - only return requested page size
    if (defined $data->{pageSize} && $data->{pageSize} ne '' && $count == $data->{pageSize}){
      $nextPageToken = $n;
      last;
    }

    ## format
    my $referenceSet;
    $referenceSet->{id}           = $refset_hash->{id};
    $referenceSet->{name}         = $refset_hash->{name};
    $referenceSet->{md5checksum}  = $refset_hash->{md5};
    $referenceSet->{ncbiTaxonId}  = $refset_hash->{ncbiTaxonId};
    $referenceSet->{description}  = "Homo sapiens " . $refset_hash->{id};
    $referenceSet->{assemblyId}   = $refset_hash->{id};
    $referenceSet->{sourceURI}    = $refset_hash->{sourceURI}; 
    $referenceSet->{sourceAccessions} = $refset_hash->{sourceAccessions} ;
    $referenceSet->{isDerived}    = $refset_hash->{isDerived};

    push @referenceSets, $referenceSet;

    $count++;
  
  }

  return ( \@referenceSets, $nextPageToken );

}

1;
