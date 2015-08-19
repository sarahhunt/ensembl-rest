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

package EnsEMBL::REST::Model::ga4gh::callSet;

use Moose;
extends 'Catalyst::Model';
use Data::Dumper;
use Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;

with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');


sub build_per_context_instance {
  my ($self, $c, @args) = @_;

  return $self->new({ context => $c, %$self, @args });
}


=head2 fetch_callSets

  POST request entry point

  ga4gh/callsets/search -d 

{ "variantSetId": 1,
 "name": '' ,
 "pageToken":  null,
 "pageSize": 10
}

=cut

sub fetch_callSets {

  my ($self, $data ) = @_;

  my ($callsets, $newPageToken ) = $self->fetch_batch($data);

  my $return_data = { callSets  => $callsets} ;
  $return_data->{pageToken} = $newPageToken if defined $newPageToken ;

  return $return_data;
}

=head2 fetch_batch

Read config, apply filtering and format records
Handle paging and return nextPageToken if required

=cut

sub fetch_batch{

  my $self = shift;
  my $data = shift;

  ## ind_id to start taken from page token - start from 0 if none supplied [!!put ids back]
  $data->{pageToken} = 0  if (! defined $data->{pageToken} || $data->{pageToken} eq "");
  my $next_ind_id   =  $data->{pageToken} ;

  my @callsets;
  my $n = 1;
  my $newPageToken; ## save id of next individual to start with


  $Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $self->{ga_config};
  my $vca = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new();

  my $count_ind = 0;## for paging [!!put ids back]
  foreach my $dataSet ( @{$vca->fetch_all} ){

    ## loop over callSets
    $dataSet->{sample_populations} = $dataSet->{_raw_populations};
    foreach my $callset_id( sort( keys %{$dataSet->{sample_populations}} ) ){

      ## stop if there are enough saved for the required batch size
      last if defined  $newPageToken ;


      ## filter by name if required
      next if defined $data->{name} && $callset_id !~ /$data->{name}/; 

      ## filter by id from GET request (these will be different to names eventually)
      next if defined $data->{req_callset} && $callset_id !~ /$data->{req_callset}/;
 
      ## filter by variant set if required
      next if defined $data->{variantSetId} && $data->{variantSetId} ne ''  && 
        $data->{variantSetId} ne  $dataSet->{sample_populations}->{$callset_id}->[0] ;


      ## paging
      $count_ind++;
      ## skip ind already reported
      next if $count_ind <$next_ind_id;
      $newPageToken = $count_ind + 1  if (defined  $data->{pageSize}  &&  $data->{pageSize} =~/\w+/ && $n == $data->{pageSize});


      ## save info
      my $callset;
      $callset->{sampleId}       = $callset_id;
      $callset->{id}             = $callset_id;
      $callset->{name}           = $callset_id;
      $callset->{variantSetIds}  = [$dataSet->{sample_populations}->{$callset_id}->[0]]; 
      $callset->{info}           = {"assembly_version" => [ $dataSet->assembly]};
      $callset->{created}        = $dataSet->created();
      $callset->{updated}        = $dataSet->updated();
      push @callsets, $callset;
      $n++; ## keep track of batch size

    }
  }

  return (\@callsets, $newPageToken);
}


=head2 getCallSet

  GET entry point - get a CallSet by ID.
  ga4gh/callset/{id}

=cut

sub get_callSet{

  my ($self, $id ) = @_; 

  my $c = $self->context();

  my $data = {req_callset => $id};

  ## extract required call set 
  my ($callSets, $newPageToken ) = $self->fetch_batch($data);

  $self->context()->go( 'ReturnError', 'custom', [ " Failed to find a callSet with id: $id"])
    unless defined $callSets && defined $callSets->[0];

  return $callSets->[0];
}


1;
