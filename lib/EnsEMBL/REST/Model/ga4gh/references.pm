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

package EnsEMBL::REST::Model::ga4gh::references;

use Moose;
extends 'Catalyst::Model';
use Data::Dumper;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');

our $species = 'homo_sapiens';

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}



sub fetch_references {
  
  my $self = shift;

  my $c = $self->context();
  
  my $post_data = $c->req->data;

#  $c->log->debug(Dumper $post_data);

  $c->go( 'ReturnError', 'custom', [ ' Error - search by md5sum not currently supported'])
    if exists $post_data->{md5checksums} ;

  my %req;
  if (exists $post_data->{accessions}){
    foreach my $acc (@{$post_data->{accessions}}){
      $req{$acc}++;
    }
  }
  if (exists $post_data->{referenceNames}){
    foreach my $name (@{$post_data->{referenceNames}}){
      $req{$name}++;
    }
  }

  return $self->getReferences(\%req);

}



sub getReferences{

  my $self = shift;
  my $req  = shift;

  my @references;
  my $count = 0;
  my $nextToken;
  my $start = 0;

  ## fetch version info needed for ftp path once
  my ($ens_version, $assembly ) = $self->getVersion(); 

  my $c = $self->context();

  my $slice_ad = $c->model('Registry')->get_adaptor($species, 'Core', 'slice');

  my @slices = @{$slice_ad->fetch_all('toplevel')};

  foreach my $slice(@slices){

    my $name   = $slice->seq_region_name();   
    my $tmp_id = $slice->name();
 
    next if (scalar (keys %$req)) > 0  && ! $req->{$name} && ! $req->{$tmp_id};
 
    ## rough paging - is order reliable?
    $start = 1 if defined $c->req->data->{pageToken} && $tmp_id eq $c->req->data->{pageToken};
    next if defined $c->req->data->{pageToken} && $start ==0;

    $count++;
    if (defined $c->req->data->{pageSize} && $count > $c->req->data->{pageSize}){

      $nextToken = $tmp_id;
      last;
    }
   
    my $ref = $self->formatSlice($slice, $ens_version, $assembly);
    push @references,  $ref;
  }

  return { references => \@references, nextPageToken => $nextToken };

}


### get single reference by 'id'

sub getReference{

  my ($self,  $get_id ) = @_; 

  my $c = $self->context();

  my $slice_ad = $c->model('Registry')->get_adaptor($species, 'Core', 'slice');

  my $slice = $slice_ad->fetch_by_name($get_id);

  ## exit if not current
  $self->context()->go( 'ReturnError', 'custom', [ " No data available for this reference : $get_id" ] )
    unless defined $slice;

  my ($ens_version, $assembly ) = $self->getVersion();

  my $reference = $self->formatSlice($slice, $ens_version, $assembly);
  return { references => [$reference]};

}

sub formatSlice{

  my $self     = shift;
  my $slice    = shift;
  my $ens_ver  = shift;
  my $assembly = shift;
  my %ref;

  $ref{start}  = 0; ##FIX
  $ref{length} = $slice->length();
  $ref{id}     = $slice->display_id();
  $ref{name}   = $slice->seq_region_name();

  my @alternative_names = @{$slice->get_all_synonyms('RefSeq_genomic')} ;
  $ref{sourceAccessions}  = [ $alternative_names[0]->name()] if defined $alternative_names[0];
  warn "Error- 2 accessions\n" if defined $alternative_names[1];

  $ref{ncbiTaxonId} = 9609; 
  $ref{isDerived}   = 'true'; ##ambiguity codes to Ns
  $ref{isPrimary}   = 'true';

  ## Fix these
  $ref{sequenceId}  = $ref{id} ;
  $ref{md5checksum} ='';
  $ref{url} =  'ftp://ftp.ensembl.org/pub/release-'. $ens_ver .'/fasta/homo_sapiens/dna/Homo_sapiens.'. $assembly.'.dna.chromosome.' . $ref{name} . '.fa.gz' 
   if length($ref{name})<3;

  $ref{sourceDivergence} = '';

  return \%ref;

}

sub getVersion{

  my $self  = shift;

  my $core_ad     = $self->context->model('Registry')->get_DBAdaptor($species, 'core');
  my $core_meta   = $core_ad->get_MetaContainer();
  my $ens_version = $core_meta->schema_version();

  my ($highest_cs) = @{$core_ad->get_CoordSystemAdaptor->fetch_all()};
  my $assembly = $highest_cs->version();

  return ($ens_version, $assembly );

}

1;
