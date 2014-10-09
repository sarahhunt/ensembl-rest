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

package EnsEMBL::REST::Model::GAcallSet;

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

sub fetch_ga_callSet {

  my ($self, $data ) = @_;

#  print Dumper $data;

  ## format set ids if filtering by set required
  if(defined $data->{variantSetIds}->[0]){

    my %req_variantset;
    foreach my $set ( @{$data->{variantSetIds}} ){
      $req_variantset{$set} = 1; 
    }
    $data->{req_variantsets} = \%req_variantset;
  } 

  ## extract required sets
  return $self->fetch_callSets($data);

}


sub fetch_callSets{

  my $self = shift;
  my $data = shift;

  ## ind_id to start taken from page token - start from 0 if none supplied [!!put ids back]
  $data->{pageToken} = 0  if (! defined $data->{pageToken} || $data->{pageToken} eq "");
  my $next_ind_id   =  $data->{pageToken} ;

  my @callsets;
  my $n = 0;
  my $newPageToken; ## save id of next individual to start with


   ## read config
  open my $cf, $config_file ||
    $self->context()->go( 'ReturnError', 'custom', [ " Failed to find config to extract set ids variantSets"]);

  local $/ = undef;
  my $json_string = <$cf>;
  close $cf;

  my $config = JSON->new->decode($json_string) ||  
    $self->context()->go( 'ReturnError', 'custom', [ " Failed to parse config for variantSets"]); 

  my $count_ind = 0;## for paging [!!put ids back]
  foreach my $hash(@{$config->{collections}}) {

    ## loop over variantSets
    foreach my $varSetId( sort(keys %{$hash->{sets}}) ){
      ## limit by data set if required
      next if defined  $data->{req_variantsets} &&  ! defined $data->{req_variantsets}->{$varSetId}; 
    }

    ## loop over callSets
    foreach my $callset_id( sort( keys %{$hash->{individual_populations}} ) ){
      
      ## limit by variant set if required
      next if defined $data->{req_variantsets} && ! defined $data->{req_variantsets}->{ $hash->{individual_populations}->{$callset_id}->[0] } ;
 
      ## paging
      $count_ind++;
      next if $count_ind <$next_ind_id;
 
      if (defined  $data->{pageSize}  &&  $data->{pageSize} =~/\w+/ && $n == $data->{pageSize}){
        $newPageToken = $count_ind;
        last;
      }
      
      ## save info
      my $callset;
      $callset->{sampleId}       = $callset_id;
      $callset->{id}             = $callset_id;
      $callset->{name}           = $callset_id;
      $callset->{variantSetIds}  = [$hash->{individual_populations}->{$callset_id}->[0]]; 
      $callset->{info}           = { "assembly_version" => "GRCh37"};
      push @callsets, $callset;
      $n++;

    }
  }

 
  my $return_data = { "callSets"  => \@callsets};
  $return_data->{"pageToken"} = $newPageToken if defined $newPageToken ;

  return $return_data; 
  
}



1;
