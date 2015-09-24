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

  my ($references, $nextPageToken);

  if ( $self->context->req->data->{referenceSetId} =~ /compliance/i){
      ($references, $nextPageToken) =  $self->get_fake_data();
  }
  else{
      ($references, $nextPageToken) =   $self->get_data();
  }

  return ({ references    => $references, 
            nextPageToken => $nextPageToken});
}


## read data from database & add MD5's from config
sub get_data{

  my $self = shift;

  my $c = $self->context();

  my $post_data = $c->req->data;

#  $c->log->debug(Dumper $post_data);

  my @references;
  my $count = 0;
  my $nextToken;
  $post_data->{pageToken} = 0 unless defined $post_data->{pageToken};


  ## read config to get look up from seq name to MD5
  my $seq_to_md5 = $self->get_md5_lookup();

  ## fetch version info needed for ftp path once
  my ($ens_version, $assembly ) = $self->getVersion(); 


  ## extract & loop through sequences
  my $slice_ad = $c->model('Registry')->get_adaptor($species, 'Core', 'slice');
  my @slices = @{$slice_ad->fetch_all('toplevel')};


  foreach (my $n = $post_data->{pageToken}; $n < scalar(@slices); $n++){

    my $name   = $slices[$n]->seq_region_name();   
    my $tmp_id = $slices[$n]->name(); 

    ## filter by any supplied attributes
    next unless defined $seq_to_md5->{$name}; ## patch sequence may not be in the set
    next if defined $post_data->{md5checksum} && $post_data->{md5checksum} ne $seq_to_md5->{$name};
    next if defined $post_data->{accession}   && $post_data->{accession} ne $slices[$n]->get_all_synonyms('RefSeq_genomic')->[0]->name(); 

    ## rough paging - assuming order reliable
    $count++;
    if (defined $post_data->{pageSize} && $count > $post_data->{pageSize}){
      $nextToken = $n;
      last;
    }
   
    my $ref = $self->formatSlice($slices[$n], $ens_version, $assembly, $seq_to_md5->{$name} ); 
    push @references,  $ref;
  }


  return (\@references, $nextToken ) ;

}


### get single reference by 'id'

sub getReference{

  my ($self,  $get_id ) = @_; 

  ## read config to get look up from MD5
  my $config = $self->read_config();

  my $seq;
  foreach my $referenceSet (@{$config->{referenceSets}}){
    foreach my $seqname( keys %{$referenceSet->{sequences}}){
# warn "Found $seqname => " . $referenceSet->{sequences}->{$seqname} . " Looking for $get_id \n";
      next if $referenceSet->{sequences}->{$seqname} ne $get_id ;

      $seq->{name}  = $seqname;
      $seq->{md5}   = $referenceSet->{sequences}->{$seqname};
      $seq->{refset_id}   = $referenceSet->{name};## track if compliance 
      $seq->{refset_name} = $referenceSet->{id};
    }
  }

  ## exit if MD5 not found
  $self->context()->go( 'ReturnError', 'custom', [ " No data available for this reference : $get_id" ] )
    unless defined $seq->{name};


  ## compliance suite uses fake data - not in db => handle separately
  if( $seq->{refset_name} =~ /compliance/){
    return $self->format_config_hash( $seq );
  }
  
  ### look up full info from database

  my $c = $self->context();

  my $slice_ad = $c->model('Registry')->get_adaptor($species, 'Core', 'slice');

  my $slice = $slice_ad->fetch_by_region('toplevel', $seq->{name});

  ## exit if not current
  $self->context()->go( 'ReturnError', 'custom', [ " No data available for this reference : $get_id" ] )
    unless defined $slice;

  my ($ens_version, $assembly ) = $self->getVersion();

  my $reference = $self->formatSlice($slice, $ens_version, $assembly, $get_id);
  return $reference;

}

sub formatSlice{

  my $self     = shift;
  my $slice    = shift;
  my $ens_ver  = shift;
  my $assembly = shift;
  my $md5      = shift;

  my %ref;

  $ref{start}  = 0; ##FIX
  $ref{length} = $slice->length();
  $ref{id}     = $md5;
  $ref{name}   = $slice->seq_region_name();

  my @alternative_names = @{$slice->get_all_synonyms('RefSeq_genomic')} ;
  $ref{sourceAccessions}  = [ $alternative_names[0]->name()] if defined $alternative_names[0];
  warn "Error- 2 accessions\n" if defined $alternative_names[1];

  $ref{ncbiTaxonId} = 9609; 
  $ref{isDerived}   = 'true'; ##ambiguity codes to Ns
  $ref{isPrimary}   = 'true';

  ## Fix these
  $ref{sequenceId}  = $ref{id} ;
  $ref{md5checksum} = $md5;
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

## TEMP: get md5's from config for get request 
sub get_md5_lookup{

  my $self = shift;
  my $set  = shift;

  my $config = $self->read_config();

  my $md5;

  foreach my $referenceSet (@{$config->{referenceSets}}){
    next if defined $set && $referenceSet->{id} ne $set;

    foreach my $seqname( keys %{$referenceSet->{sequences}}){
      $md5->{$seqname} = $referenceSet->{sequences}->{$seqname};
    }
  }
  return $md5;
}



##read config from JSON file to get MD5's
sub read_config{

  my $self = shift;
  my $set  = shift;

  open IN, $self->{ga_reference_config} ||
    $self->context()->go( 'ReturnError', 'custom', ["ERROR: Could not read from config file $self->{ga_reference_config}"]);
  local $/ = undef;
  my $json_string = <IN>;
  close IN;


  my $config = JSON->new->decode($json_string) ||
    $self->context()->go( 'ReturnError', 'custom', ["ERROR: Failed to parse config file $self->{ga_reference_config}"]);

  $self->context()->go( 'ReturnError', 'custom', [ " No data available " ] )
    unless $config->{referenceSets} && scalar @{$config->{referenceSets}};

  return $config;
}

## compliance data is a random slice => not in the db
## bolted on - should be handled better

sub get_fake_data{

  my $self = shift;

  my $c = $self->context();

  ## get array of sequences from config
  my $refseqs = $self->sequences_from_config(); 

  my $formatted_seqs;

  foreach my $refseq(@{$refseqs}){

    my $ref = $self->format_config_hash($refseq);
    push @{$formatted_seqs}, $ref;
  }
  return $formatted_seqs;
}

## used by POST & GET
sub format_config_hash {

  my $self   = shift;
  my $refseq = shift;
  my $ref;

  $ref->{id}          = $refseq->{md5};
  $ref->{name}        = $refseq->{name};
  $ref->{md5checksum} = $refseq->{md5};
  $ref->{sequenceId}  = $refseq->{md5} ; ## what is better here?

  ## Potential error
  $ref->{ncbiTaxonId} = 9609;
  $ref->{isDerived}   = 'true';
  $ref->{isPrimary}   = 'true';
  $ref->{url}         = 'https://github.com/ga4gh/compliance';

  ## Fix these
  $ref->{sourceAccessions}  = '';
  $ref->{start}  = 0;
  $ref->{length} = '';
  $ref->{sourceDivergence} = ''; 

  return $ref;
}

## compliance data not in db - get everything from config 
sub get_sequences_from_config{

  my $self = shift;
  my $set  = shift;

  my $config = $self->read_config();

  my $refseqs;

  foreach my $referenceSet (@{$config->{referenceSets}}){
    next if $referenceSet->{id} ne $set;

     foreach my $seqname( keys ${$referenceSet->{sequences}}){
      
      my $seq;
      $seq->{name}        = $referenceSet->{name};
      $seq->{md5}         = $referenceSet->{sequences}->{$seqname};
      $seq->{refset_id}   = $referenceSet->{name};## track if compliance 
      $seq->{refset_name} = $referenceSet->{id};
      push @{$refseqs}, $seq;
    }
  }
  return $refseqs;
}


1;
