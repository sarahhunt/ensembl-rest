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

### using md5 of concatenated sequence MD5s as id

sub searchReferenceSet {
  
  my $self = shift;

  my $c = $self->context();
  
  my $post_data = $c->req->data;

#  $c->log->debug(Dumper $post_data);

  my ( $referenceSets, $nextPageToken)  =  $self->fetchData( $post_data );
  return ({ referenceSets => $referenceSets,
            nextPageToken => $nextPageToken });

}

sub getReferenceSet {

  my $self = shift;
  my $id = shift;

  my $data;
  $data->{id} = $id;
  my ($referenceSets, $pageToken) =  $self->fetchData( $data );

  return ($referenceSets->[0]);

}


## send both post & get here as few sets to check

sub fetchData{

  my ($self,  $data ) = @_; 

  my $c = $self->context();

  ## read config
  my $referenceSets = $self->read_config();

  my @referenceSets;

  $data->{pageToken} = 0 unless defined $data->{pageToken}; 
  my $nextPageToken;
  my $count = 0;
  foreach( my $n = $data->{pageToken}; $n <  scalar @{$referenceSets}; $n++ ) {

    my $refset_hash = $referenceSets->[$n];

    ## filter if an attrib supplied
    next if defined $data->{id}          &&  $data->{id}          ne $refset_hash->{id}; ##  GET
    next if defined $data->{md5checksum} &&  $data->{md5checksum} ne $refset_hash->{md5};
    next if defined $data->{accession}   &&  $data->{accession}   ne $refset_hash->{accession}->[0]; ##FIX for other accessions
    next if defined $data->{assemblyId}  &&  $data->{assemblyId}  ne $refset_hash->{id};

    if (defined $data->{pageSize} && $count == $data->{pageSize}){
      $nextPageToken = $n;
      last;
    }

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

  return (\@referenceSets, $nextPageToken);

}


##read config from JSON file

sub read_config{

  my $self = shift;

  open IN, $self->{ga_reference_config} ||
    $self->context()->go( 'ReturnError', 'custom', ["ERROR: Could not read from config file $self->{ga_reference_config}"]);
  local $/ = undef;
  my $json_string = <IN>;
  close IN;


  my $config = JSON->new->decode($json_string) ||
    $self->context()->go( 'ReturnError', 'custom', ["ERROR: Failed to parse config file $self->{ga_reference_config}"]);

  $self->context()->go( 'ReturnError', 'custom', [ " No data available " ] )
    unless $config->{referenceSets} && scalar @{$config->{referenceSets}};

  return $config->{referenceSets};
}


=head don't take from db as compliance data is fake 

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
  $referenceSet->{isDerived}    = 'true';


  return { referenceSets => [$referenceSet]};

}
=cut




1;
