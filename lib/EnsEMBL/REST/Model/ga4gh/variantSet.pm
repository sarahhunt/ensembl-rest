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

package EnsEMBL::REST::Model::ga4gh::variantSet;

use Moose;
extends 'Catalyst::Model';
use Data::Dumper;
use Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');


sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });

}

=head fetch_variantSets

  POST request entry point

=cut

sub fetch_variantSets {

  my ($self, $data ) = @_;

  ## extract required variant sets
  my ($variantSets, $newPageToken ) = $self->fetch_sets($data);

  my $ret = { variantSets => $variantSets};
  $ret->{pageToken} = $newPageToken  if defined $newPageToken ;

  return $ret;
}

=head fetch_sets

Read config, apply filtering and format records

=cut

sub fetch_sets{

  my $self = shift;
  my $data = shift;

  ## varset id to start from is the page token - start from 0 if none supplied
  my $next_set_id = 0;
  $next_set_id    = $data->{pageToken} if ( defined $data->{pageToken} && $data->{pageToken} ne "");



  ## read config
  $Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $self->{ga_config};
  my $vca = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new();

  ## extract requested data
  my %variantSets;  

  foreach my $dataSet(@{$vca->fetch_all} ) {
     $dataSet->use_db(0);
    ## limit by data set if required
    next unless !defined $data->{datasetId} || $dataSet->id() eq $data->{datasetId} ; 

    ## get info descriptions from one of the VCF files    
    my $meta = $self->get_info($dataSet);

    ## add summary of essential info for meta for the data set from the config
    foreach my $key ("assembly", "source_name", "source_url"){
      my %meta;
      $meta{key}   = $key;
      $meta{value} = $dataSet->$key;
      push @{$meta}, \%meta;
    }

    ## loop through and save all available variantsets
    foreach my $population (@{$dataSet->get_all_Populations}){

      my $varset = $population->dbID();

      ## limit by variant set if required (for GET)
      next unless !defined  $data->{req_variantset} || $data->{req_variantset} eq $varset;

      $variantSets{$varset}{datasetId}   = $dataSet->id();
      my @m = @{$meta};
      my %m  = ( "key"   => "set_name", 
                 "value" => $population->name()
               );
      push @m, \%m;
      $variantSets{$varset}{metadata}    = \@m;
    }
  }

  ## create response with correct number of variantSets
  my @varsets;
  my $n = 0;
  my $newPageToken; ## save id of next variantSet to start with

  foreach my $varset_id (sort (keys %variantSets)){    

    ## paging - skip if already returned
    next if $varset_id < $next_set_id;
   
    if (defined $data->{pageSize} &&  $n == $data->{pageSize}){
      ## set next token and stop storing if page size reached
      $newPageToken = $varset_id if defined $data->{pageSize} && $n == $data->{pageSize};
      last;
    }
    ## store variantSet info to return
    my $varset;
    $varset->{id} = $varset_id;
    $varset->{datasetId} = $variantSets{$varset_id}->{datasetId};
    $varset->{metadata}  = $variantSets{$varset_id}->{metadata};
    push @varsets, $varset;
    $n++;

  }

  ## check there is something to return
  $self->context()->go( 'ReturnError', 'custom', [ " Failed to find any variantSets for this query"]) if $n ==0;
 
  my $ret = { variantSets => \@varsets};
  $ret->{pageToken} = $newPageToken  if defined $newPageToken ;
 
  return (\@varsets, $newPageToken );
  
}

=head2 get_info

 get info descriptions from one of the VCF files   
 to report as meta data

=cut

sub get_info{

  my $self = shift;
  my $dataSet  = shift;

  my $vcf_file  = $dataSet->filename_template();
  $vcf_file =~ s/\#\#\#CHR\#\#\#/22/;
  $vcf_file = $self->{dir} .'/'. $vcf_file;

  my $parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open( $vcf_file ) || die "Failed to get parser : $!\n";

  my @meta;

  ## metadata for info & format 
  foreach my $type ("INFO", "FORMAT"){

    my $data = $parser->get_metadata_by_pragma($type);

    foreach my $d(@{$data}){
      my $out;
      $out->{key} = $type;
      foreach my $k (keys %$d){
        $out->{"\L$k"} = $d->{$k};
      }
      push @meta, $out;
    }
  }
 
  return \@meta;
}




=head2 getVariantSet

  Gets a VariantSet by ID.
  GET /variantsset/{id}

=cut

sub getVariantSet{

  my ($self, $id ) = @_; 

  my $c = $self->context();

  my $data = {req_variantset => $id};

  ## extract required variant set 
  my ($variantSets, $newPageToken ) = $self->fetch_sets($data);

  ## would exit earlier if no data
  return $variantSets->[0];
}

1;
