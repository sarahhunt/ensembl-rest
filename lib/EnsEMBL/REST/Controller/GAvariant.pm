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

package EnsEMBL::REST::Controller::GAvariant;

use Moose;
use namespace::autoclean;
use Try::Tiny;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;
require EnsEMBL::REST;
EnsEMBL::REST->turn_on_config_serialisers(__PACKAGE__);

=pod

POST requests : /variants/

{ "variantSetIds": [1],
 "variantName": '' ,
 "callSetIds": [],
 "referenceName": 7,
 "start":  140419275,
 "end": 140429275,
 "pageToken":  null,
 "maxResults": 10
}

application/json

=cut

BEGIN {extends 'Catalyst::Controller::REST'; }


sub get_request_POST {
  my ( $self, $c ) = @_;

}


sub get_request: Chained('/') PathPart('variants') ActionClass('REST')  {
  my ( $self, $c ) = @_;
  my $post_data = $c->req->data;

  $c->log->debug(Dumper $post_data);

  ## required by spec, so check early
  $c->go( 'ReturnError', 'custom', [ ' Cannot find "referenceName" key in your request' ] ) 
    unless exists $post_data->{referenceName} ;

  $c->go( 'ReturnError', 'custom', [ ' Cannot find "start" key in your request'])         
    unless exists $post_data->{start};

  $c->go( 'ReturnError', 'custom', [ ' Cannot find "end" key in your request'])   
    unless exists $post_data->{end};

  $c->go( 'ReturnError', 'custom', [ ' Cannot find "variantSetIds" key in your request'])
    unless exists $post_data->{variantSetIds}->[0];


  $c->go( 'ReturnError', 'custom', [ ' key "maxResults" must be a positive number'])
    if exists $post_data->{maxResults} && $post_data->{maxResults} <1;


  my $gavariant;

  try {
    $gavariant = $c->model('GAvariant')->fetch_gavariant($post_data);
  } catch {
    $c->go('ReturnError', 'from_ensembl', [$_]);
  };

  $self->status_ok($c, entity => $gavariant);
}


__PACKAGE__->meta->make_immutable;

1;
