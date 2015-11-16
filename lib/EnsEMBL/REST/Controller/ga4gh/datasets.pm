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

package EnsEMBL::REST::Controller::ga4gh::datasets;

use Moose;
use namespace::autoclean;
use Try::Tiny;
use Data::Dumper;

require EnsEMBL::REST;
EnsEMBL::REST->turn_on_config_serialisers(__PACKAGE__);

=pod

POST /datasets/search -d { "pageSize": 2, "pageToken":3}

GET requests : /dataset/id

returns id & description

=cut

BEGIN {extends 'Catalyst::Controller::REST'; }


sub get_request_POST {
  my ( $self, $c ) = @_;

}

sub get_request: Chained('/') PathPart('ga4gh/datasets/search') ActionClass('REST')  {
  my ( $self, $c ) = @_;

  my $post_data = $c->req->data;

  ## set a maximum page size - not likely to be neccessary
  $post_data->{pageSize} =  100 if !defined $post_data->{pageSize} || $post_data->{pageSize} > 100; 

  my $dataset;

  try {
    $dataset = $c->model('ga4gh::datasets')->fetch_datasets($post_data);
  } catch {
    $c->go('ReturnError', 'from_ensembl', [$_]);
  };

  $self->status_ok($c, entity => $dataset);

}




sub id: Chained('/') PathPart('ga4gh/datasets') ActionClass('REST') {}

sub id_GET {
  my ($self, $c, $id) = @_;

  my $dataset;

  try {
    $dataset = $c->model('ga4gh::datasets')->getDataset($id);
  } catch {
    $c->go('ReturnError', 'from_ensembl', [qq{$_}]) if $_ =~ /STACK/;
    $c->go('ReturnError', 'custom', [qq{$_}]);
  };
  $self->status_ok($c, entity => $dataset);
}

__PACKAGE__->meta->make_immutable;

1;