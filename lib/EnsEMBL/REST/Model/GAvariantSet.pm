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

package EnsEMBL::REST::Model::GAvariantSet;

use Moose;
extends 'Catalyst::Model';
use Data::Dumper;


with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');
our $config_file = "/home/vagrant/src/ensembl-rest/ga_vcf_config.json"; 

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub fetch_ga_variantSet {

  my ($self, $data ) = @_;

  ## format set ids if filtering by set required
  if(defined $data->{datasetIds}->[0]){

    my %req_dataset;
    foreach my $set ( @{$data->{datasetIds}} ){
      $req_dataset{$set} = 1; 
    }
    $data->{req_datasets} = \%req_dataset;
  } 


  ## extract required variant sets
  my $varsets = $self->fetch_sets($data);

  return ({"variantSets" => $varsets});


}




sub fetch_sets{

  my $self = shift;
  my $data = shift;

  ## ind_id to start from is appended to page token - start from 0 if none supplied
  $data->{pageToken} = 0  unless ( defined $data->{pageToken} && $data->{pageToken} ne "");
  my $next_set_id    = $data->{pageToken} ;


  ## read config
  open my $cf, $config_file ||
    $self->context()->go( 'ReturnError', 'custom', [ " Failed to find config to extract set ids variantSets"]);

  local $/ = undef;
  my $json_string = <$cf>;
  close $cf;

  my $config = JSON->new->decode($json_string) ||  
    $self->context()->go( 'ReturnError', 'custom', [ " Failed to parse config for variantSets"]); 


  ## extract requested data
  my %variantSets;  
  foreach my $hash(@{$config->{collections}}) {

    ## limit by data set if required
    next if defined  $data->{req_datasets} &&  ! defined $data->{req_datasets}->{ $hash->{datasetId} }; 

    ## save variantSets by dataset 
    foreach my $varset(keys %{$hash->{sets}}){

      $variantSets{$varset}{datasetId}  = $hash->{datasetId};
      $variantSets{$varset}{varsetdesc} = $hash->{sets}->{$varset};
    }
  }

  ## create response
  my @varsets;
  my $n = 0;
  my $newPageToken; ## save id of next variationSet to start with

  foreach my $varset_id (sort (keys %variantSets)){    
 
    ## paging
    next if $varset_id < $next_set_id;

    if ( $n == $data->{pageSize}){
      $newPageToken = $varset_id if defined $data->{pageSize} && $n == $data->{pageSize};
      last;
    }
      
    my $varset;
    $varset->{id}          = $varset_id;
    $varset->{datasetId}   = $variantSets{$varset_id}->{datasetId};
    ## add meta
    push @varsets, $varset;
    $n++;
  }

  ## check there is something to return
  $self->context()->go( 'ReturnError', 'custom', [ " Failed to find any variantSets for this dataset id"]) if $n ==0;
 
  push @varsets, {"pageToken" => $newPageToken } if defined $newPageToken ;

  return (\@varsets);
  
}



1;
