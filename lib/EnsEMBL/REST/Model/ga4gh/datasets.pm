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

package EnsEMBL::REST::Model::ga4gh::datasets;

use Moose;
extends 'Catalyst::Model';

use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;
use Digest::MD5 qw(md5_hex);
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');


sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub fetch_datasets{
  my ($self, $data ) = @_;

  my @datasets;

  ## paging
  my $count = 0;
  my $next;      
  my $start = 1;
  $start = 0 if defined $data->{pageToken} && $data->{pageToken} ne ''; 

  my $collections = $self->sort_collections;
  foreach my $id(keys %{$collections}){

    $start = 1 if (defined $data->{pageToken} &&  $data->{pageToken} eq $id);

    ## skip if in last batch
    next if $start == 0 ;

    ## store start of next batch before exiting this one
    if ($count == $data->{pageSize}){
      $next = $id;
      last; 
    }

    ## store and increment counter
    my %done;
    my $data_source = $collections->{$id}->source_name() ;
    unless (defined $done{$data_source}){

      ## don't use source->dbId as may clash with other entries 
      my $dataset = { id => $id, description => $data_source };
      push @datasets, $dataset;
      $count++;
      $done{$data_source} = 1;
    }
  }

  $self->context()->go( 'ReturnError', 'custom', [ " No datasets is available"])
     unless scalar(@datasets) >0 ;

  return { datasets => \@datasets, nextPageToken => $next};


}

sub sort_collections{

  my ($self) = @_;

  ## read config
  $Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $self->{ga_config};
  my $vca = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new();

  my %collections;
  foreach my $collection(@{$vca->fetch_all} ) { 

    my $ga_id = md5_hex($collection->source_name());
    $collections{$ga_id} = $collection;
  }

  return \%collections;

}



sub getDataset{

  my ($self, $id ) = @_; 


  ## read config
  $Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $self->{ga_config};
  my $vca = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new();

  ## extract requested data
  my $dataset;  

  foreach my $collection(@{$vca->fetch_all} ) {

     my $ga_id = md5_hex($collection->source_name());

     next unless $ga_id eq $id;

     $dataset = { id => $ga_id, description => $collection->source_name() };
  }

  $self->context()->go( 'ReturnError', 'custom', [ " No dataset is available with this id"])
     unless defined $dataset;

  return $dataset;  

}


with 'EnsEMBL::REST::Role::Content';

__PACKAGE__->meta->make_immutable;

1;
 
