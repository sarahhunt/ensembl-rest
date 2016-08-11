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

package EnsEMBL::REST::Model::ga4gh::phenotypeAssociationSets;

use Moose;
extends 'Catalyst::Model';
use Catalyst::Exception;
use Digest::MD5 qw(md5_hex);
use Scalar::Util qw/weaken/;


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
sub searchPhenotypeAssociationSets {
  
  my $self   = shift;
  my $data   = shift;

  ## only supporting one dataset currently
  return { phenotype_association_sets   => [],
           next_page_token => undef
         } unless $data->{datasetId} =~ /Ensembl/;  ## use ClinVar etc as database or release id? Former more robust?


  my ($phenotypeAssociationSet, $nextPageToken)  = $self->fetch_sets($data);

  return { phenotype_association_sets   => $phenotypeAssociationSet,
           next_page_token              => $nextPageToken
         }; 
}

## GET entry point
sub getPhenotypeAssociationSet{

  my ($self, $id ) = @_; 

  my $data = { required_set => $id,
               pageSize    => 1
              };

  my ($sets, $nextPageToken) =  $self->fetch_sets($data); 
  return $sets->[0];

}

sub fetch_sets{

  my $self = shift;
  my $data = shift;

  my @sets;

  my $source_ad = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Variation', 'Source');

  my $n = 0;
  my $nextPageToken;
  ## only a few sources supported as phenotypeassociation sets
  my $supported_sets = $self->supported_sets();
  foreach my $id (sort keys %{$supported_sets}){

    ## filter for GET
    next if defined $data->{required_set} && $id ne $data->{required_set};

    my $source_ob = $source_ad->fetch_by_name( $supported_sets->{$id} );
    next unless $source_ob;

    if (defined $data->{pageSize} && $n == $data->{pageSize}){
      $nextPageToken = $id;
      last;
    }

    push @sets, { id         => $id,
                  name       => $supported_sets->{$id},
                  dataset_id => 'Ensembl',
                  info       => { version =>  $source_ob->version(),
                                   url    =>  $source_ob->url()
                                 }
                 };

    $n++;

  }
  return (\@sets, $nextPageToken);
}


1;
