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

package EnsEMBL::REST::Controller::ga4gh::annotate;

use Moose;
use namespace::autoclean;
use Try::Tiny;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;
require EnsEMBL::REST;
EnsEMBL::REST->turn_on_config_serialisers(__PACKAGE__);

=pod

POST requests: /ga4gh/variants/annotate -d 
{ "variants:" [org.ga4gh.models.Variant ],
  "pageToken":  null,
  "pageSize": 10
}

=cut

BEGIN {extends 'Catalyst::Controller::REST'; }


sub annotateVariants: Chained('/') PathPart('ga4gh/annotate/variants') ActionClass('REST') {}


sub annotateVariants_POST{

  my ( $self, $c ) = @_;

  my $post_data = $c->req->data;

  $c->log->debug(Dumper $post_data);

  $c->go( 'ReturnError', 'custom', [ ' Cannot find "variants" key in your request' ] )
    unless exists $post_data->{variants} ;

  my $variantAnnotation;

  try {
    $variantAnnotation = $c->model('ga4gh::annotate')->annotateVariants($post_data);
  } catch {
    $c->go('ReturnError', 'from_ensembl', [$_]);
  };

  $self->status_ok($c, entity => $variantAnnotation);

}


__PACKAGE__->meta->make_immutable;

1;
