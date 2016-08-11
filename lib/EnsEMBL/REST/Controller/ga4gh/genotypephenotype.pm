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

package EnsEMBL::REST::Controller::ga4gh::genotypephenotype;

use Moose;
use namespace::autoclean;
use Try::Tiny;
use Data::Dumper;

require EnsEMBL::REST;
EnsEMBL::REST->turn_on_config_serialisers(__PACKAGE__);

=pod

POST requests : /genotypephenotype/search -d

{ "phenotype_association_set_id" : "ppp"
  "feature_ids": null,
  "phenotype_ids": null,
  "evidence":  null, 
  "page_token": null,
  "page_size": 10
}

application/json

=cut

BEGIN {extends 'Catalyst::Controller::REST'; }


sub SearchFeaturesRequest_POST {
  my ( $self, $c ) = @_;

}


sub SearchFeaturesRequest: Chained('/') PathPart('ga4gh/genotypephenotype/search') ActionClass('REST')  {
  my ( $self, $c ) = @_;
  my $post_data = $c->req->data;

  $c->log->debug(Dumper $post_data);

  ## check a set is supplied
  $c->go( 'ReturnError', 'custom', [ " A phenotype_association_set_id must be supplied "])
    unless defined $post_data->{phenotype_association_set_id};


  ## check there is something to look 
  $c->go( 'ReturnError', 'custom', [ " Feature or Phenotype ids must be supplied "])
    unless defined $post_data->{feature_ids}   ||
           defined $post_data->{phenotype_ids};


  ## set a default page size if not supplied or not a number
  $post_data->{page_size} = 50 unless (defined  $post_data->{page_size} &&  
                                      $post_data->{page_size} =~ /\d+/ &&
                                      $post_data->{page_size} >0  );

  my $g2p_results;

  try {
    $g2p_results = $c->model('ga4gh::genotypephenotype')->fetch_g2p_results($post_data);
  } catch {
    $c->go('ReturnError', 'from_ensembl', [$_]);
  };

  $self->status_ok($c, entity => $g2p_results);
}



__PACKAGE__->meta->make_immutable;

1;
