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

## switched sample list look-up to VCF file rather than config - will be slower

sub fetch_batch{

  my $self = shift;
  my $data = shift;

  ## ind_id to start taken from page token - start from 0 if none supplied [!!put ids back]
  $data->{pageToken} = 0  if (! defined $data->{pageToken} || $data->{pageToken} eq "");
  my $next_ind_id   =  $data->{pageToken} ;

  my @callsets;
  my $n = 1;
  my $newPageToken; ## save id of next individual to start with

  $ENV{ENSEMBL_VARIATION_VCF_ROOT_DIR} = $self->{geno_dir};
  $Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $self->{ga_config};
  my $vca = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new();
  
  my %collections; 
  foreach my $vcf_collection ( @{$vca->fetch_all} ){
    $collections{$vcf_collection->id()} = $vcf_collection;
  }

  my $count_ind = 0;## for paging
  foreach my $vcfc_id (sort sort_num keys %collections){
    my $vcf_collection = $collections{$vcfc_id};

    $vcf_collection->use_db(0);
    ## filter by variant set if required
    next if defined $data->{variantSetId} && $data->{variantSetId} ne ''  &&
      $data->{variantSetId} ne  $vcf_collection->id() ;

    ## loop over callSets
    my $samples = $vcf_collection->get_all_Samples(); ## returned sorted

    foreach my $sample (@{$samples}){ 

      my $sample_name = $sample->name();
      ## stop if there are enough saved for the required batch size
      last if defined  $newPageToken ;


      ## filter by name if required
      next if defined $data->{name} && $sample_name !~ /$data->{name}/; 

      ## filter by id from GET request (these will be different to names eventually)
      next if defined $data->{req_callset} && $sample_name !~ /$data->{req_callset}/;
 

      ## paging
      $count_ind++;
      ## skip ind already reported
      next if $count_ind <$next_ind_id;
      $newPageToken = $count_ind + 1  if (defined  $data->{pageSize}  &&  $data->{pageSize} =~/\w+/ && $n == $data->{pageSize});


      ## save info
      my $callset;
      $callset->{sampleId}       = $sample_name;
      $callset->{id}             = $sample_name;
      $callset->{name}           = $sample_name;
      $callset->{variantSetIds}  = [$vcf_collection->id()]; 
      $callset->{info}           = {"assembly_version" => [ $vcf_collection->assembly() ],
                                    "variantSetName"   => [ $vcf_collection->source_name()] };
      $callset->{created}        = $vcf_collection->created();
      $callset->{updated}        = $vcf_collection->updated();
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

sub sort_num{
  $a<=>$b;
}
1;
