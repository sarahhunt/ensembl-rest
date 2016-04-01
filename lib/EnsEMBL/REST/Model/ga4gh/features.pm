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

package EnsEMBL::REST::Model::ga4gh::features;

use Moose;
extends 'Catalyst::Model';

use Bio::DB::Fasta;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::Utils::VEP qw( read_cache_info );
use Scalar::Util qw/weaken/;
use Data::Dumper;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');


## What features to return?
## those required for Var Ann initially
has 'allowed_features' => ( isa => 'HashRef', is => 'ro', lazy => 1, default => sub {
  return {
    map { $_ => 1 } qw/gene transcript cds exon /
  };
});

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  weaken($c);
  return $self->new({ context => $c, %$self, @args });
}

## potential weirdness:
##  next token is: 
##     - numerically sorted transcript dbID in 'all' mode
##     - seq & pos start if in region mode

sub searchFeatures {

  my ($self, $data ) = @_; 

  ## TEMP - using genebuild version
  $data->{current_set} = $self->getSet();
  my $allowed_features = $self->allowed_features();


  ## check feature type is supported
  my %required_types;
  my $get_all_types = 1;
  if (exists $data->{featureTypes}->[0] ){
    $get_all_types = 0;

    foreach my $term (@{$data->{featureTypes}}){
      $required_types{$term} = 1 if $allowed_features->{$term} ==1;
    }
  }


  ## exit if unsupported feature requested
  Catalyst::Exception->throw(" Request for unsupported feature type" )
    unless $get_all_types == 1  ||  scalar keys(%required_types) >0  ; 


  ## get slice  
  my $sla = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Core', 'Slice');
 
  ## modify start if required for paging
  ## FIXME 1st set needs all overlapping; rest do not 
  my ($next_seq_start, $next_seq_pos ) = split/\_/, $data->{pageToken} if defined $data->{pageToken};

  $data->{start} = $next_seq_pos if defined $next_seq_pos;
  my $location   = $data->{referenceName} ."\:" . $data->{start} . "\-" . $data->{end};
  $data->{slice} = $sla->fetch_by_toplevel_location( $location );


  ## get features of required type 
  my @features;

  if( exists $required_types{transcript} || $get_all_types == 1 ){
    my $transcripts = $self->extractTranscriptsBySegment( $data );
    push @features, @{$transcripts};
  }
  if( exists $required_types{gene} || $get_all_types == 1 ){
    my $genes = $self->extractGenesBySegment( $data );
    push @features, @{$genes};
  }

  if( exists $required_types{exons} || $get_all_types == 1 ){
    my $exons = $self->extractExonsBySegment( $data );
    push @features, @{$exons};
  }


  ## sort & trim features
  my $sorted_features = sort_features(\@features);

  my @return_features;
  my $nextPageToken;
  my $feature_count = scalar @{$sorted_features};
  if(@{$sorted_features} > $data->{pageSize}){
    @return_features = splice(@{$sorted_features}, 0, $data->{pageSize});
    $nextPageToken = $sorted_features->[0]->{start};   ### FIX THIS - won't work for colocated features
  }
  else{
    @return_features = @{$sorted_features};
  }

  return ({ features      => \@return_features, 
            nextPageToken => $nextPageToken});

}

## extract data from db by region
sub extractTranscriptsBySegment{

  my ($self, $data ) = @_;

  my @features;
  my $count = 0;

  my $tra = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Core', 'Transcript');
  my $transcripts = $tra->fetch_all_by_Slice( $data->{slice} ) ;

  foreach my $tr (@{$transcripts}){

    ## if pageSize +1 reached, return
    last if $count == $data->{pageSize} +1 ;
 
    my $gafeat = $self->formatTranscript($tr, $data);
    push @features,  $gafeat;
    ## keep count for pageSize
    $count++;

  }
  return \@features;
}


sub extractGenesBySegment{

  my ($self, $data ) = @_;
  
  my @features;
  my $count = 0;
  
  my $ga = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Core', 'Gene');
  my $genes = $ga->fetch_all_by_Slice( $data->{slice} ) ;

  foreach my $gene (@{$genes}){

    ## if pageSize +1 reached, return
    last if $count == $data->{pageSize} +1 ;
 
    my $gafeat = $self->formatGene($gene, $data);
    push @features,  $gafeat;
    ## keep count for pageSize
    $count++;

  }
  return \@features;  
}



sub extractExonsBySegment{

  my ($self, $data ) = @_;
 
  my @features;
  my $count = 0;
  
  my $ea = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Core', 'Exon');
  my $exons = $ea->fetch_all_by_Slice( $data->{slice} ) ;

  foreach my $exon (@{$exons}){

    ## if pageSize +1 reached, return
    last if $count == $data->{pageSize} +1 ;

    my $gafeat = $self->formatExon($exon, $data);
    push @features,  $gafeat;
    ## keep count for pageSize
    $count++;

  }
  return \@features;
}


## support more feature types..
sub getFeature{

  my $self = shift;
  my $id   = shift;
  
  my $data;
  $data->{current_set} = $self->getSet();

  return $self->getTranscript($id, $data) if $id =~/ENST/;
  return $self->getGene($id, $data)       if $id =~/ENSG/;
  return $self->getProtein($id, $data)    if $id =~/ENSP/;
  return $self->getExon($id, $data)       if $id =~/ENSE/;

}

## look up transcript by id
## FIX: to use UUID

sub getTranscript{

  my $self = shift;
  my $id   = shift;
  my $data = shift;

  my $tra = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Core', 'Transcript');
  my $tr = $tra->fetch_by_stable_id( $id, $data );

  Catalyst::Exception->throw("  Cannot find transcript feature for id " . $id  )
    unless defined $tr ;

  return $self->formatTranscript($tr, $data);

}

## turn ensembl transcript into GA4GH transcript

sub formatTranscript{

  my $self = shift;
  my $tr   = shift;
  my $data = shift;

  my $feature = $self->formatFeature($tr, $data, 'transcript');

  $feature->{parentId}  = $tr->get_Gene()->stable_id();
  $feature->{childIds}  = [$tr->translation()->stable_id()] if defined $tr->translation();

=head
  my $strand;
  $tr->seq_region_strand() eq 1 ? $strand = 'POS_STRAND' 
                                : $strand = 'NEG_STRAND';

  my $protein_id = $tr->translation()->stable_id();

  my $gaFeature  = { id            => $tr->stable_id(),
                     parentId      => $tr->get_Gene()->stable_id(),
                     childIds      => [$protein_id],
                     featureSetId  => $data->{current_set},
                     referenceName => $tr->seq_region_name(),
                     start         => $tr->seq_region_start() - 1,
                     end           => $tr->seq_region_end(),
                     strand        => $strand
                    };


  ## look up ontology info if non cached
  $data->{ontol}->{transcript} = $self->fetchSO('transcript') 
    unless exists $data->{ontol}->{transcript};

  $gaFeature->{featureType} = $data->{ontol}->{transcript};

  

  ## what is interesting here?
  $gaFeature->{attributes} = { version => $tr->version(),
                               biotype => $tr->biotype(),
                               gene    => $tr->get_Gene()->display_id(),
                               created => $tr->created_date(),           
                               updated => $tr->modified_date(),
                               source  => $tr->source()
                             };

  $gaFeature->{attributes}->{external_name} = $tr->external_name() if defined $tr->external_name(); 
=cut
  return $feature;

}

sub getExon{

  my $self = shift;
  my $id   = shift;
  my $data = shift;

  my $ea   = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Core', 'Exon');
  my $exon = $ea->fetch_by_stable_id( $id, $data );

  Catalyst::Exception->throw(" Cannot find exon feature for id " . $id)
    unless defined $exon ;

  return $self->formatExon($exon, $data) ;

}

## turn ensembl gene into GA4GH transcript

sub formatExon{

  my $self = shift;
  my $exon = shift;
  my $data = shift;

  my $feature = $self->formatFeature($exon, $data, 'exon');

  ## set parent ids!!
  

  return $feature;

}


sub getGene{

  my $self = shift;
  my $id   = shift;
  my $data = shift;

  my $ga   = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Core', 'Gene');
  my $gene = $ga->fetch_by_stable_id( $id, $data );

  Catalyst::Exception->throw(" Cannot find gene feature for id " . $id)
    unless defined $gene ;

  return $self->formatGene($gene, $data) ;

}
 
## turn ensembl gene into GA4GH transcript

sub formatGene{

  my $self = shift;
  my $gene = shift;
  my $data = shift;

  my $feature = $self->formatFeature($gene, $data, 'gene');

  ## set child ids to transcripts
  my $childIds;
  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}){
    push @{$childIds}, $transcript->stable_id()
  }
  $feature->{childIds} = $childIds;

  return $feature;

}

sub formatFeature{

  my $self = shift;
  my $feat = shift;
  my $data = shift;
  my $type = shift;

  my $strand;
  $feat->seq_region_strand() eq 1 ? $strand = 'POS_STRAND' 
                                  : $strand = 'NEG_STRAND';



  my $gaFeature  = { id            => $feat->stable_id(),
                     parentId      => undef, 
                     childIds      => [],
                     featureSetId  => $data->{current_set},
                     referenceName => $feat->seq_region_name(),
                     start         => $feat->seq_region_start() - 1,
                     end           => $feat->seq_region_end(),
                     strand        => $strand
                    };


  ## look up ontology info if non cached
  $data->{ontol}->{$type} = $self->fetchSO($type) 
    unless exists $data->{ontol}->{$type};

  $gaFeature->{featureType} = $data->{ontol}->{$type};

  

  ## what is interesting here?
  $gaFeature->{attributes} = { version => $feat->version(),  
                               created => $feat->created_date(),
                               updated => $feat->modified_date()
                             };

  unless ($type =~/exon/){ 
    $gaFeature->{attributes}->{biotype} = $feat->biotype();
    $gaFeature->{attributes}->{status}  = $feat->status();
    $gaFeature->{attributes}->{source}  = $feat->source();

    $gaFeature->{attributes}->{external_name} = $feat->external_name() 
      if defined $feat->external_name();
  }

  return $gaFeature;

}

sub getProtein{

  my ($self, $id, $data) = @_;

  my $ta   = $self->context->model('Registry')->get_adaptor('homo_sapiens', 'Core', 'translation');
  my $translation = $ta->fetch_by_stable_id( $id, $data );

  Catalyst::Exception->throw(" Cannot find gene feature for id " . $id)
    unless defined $translation ;

  return $self->formatProtein($translation, $data) ;

}

## Not a great match.. handle differently?
sub formatProtein{

  my $self = shift;
  my $feat = shift;
  my $data = shift;

  my $gaFeature  = { id            => $feat->stable_id(),
                     parentId      => $feat->transcript()->stable_id(),
                     childIds      => [],
                     featureSetId  => $data->{current_set},
                     referenceName => undef,
                     start         => $feat->genomic_start() - 1,
                     end           => $feat->genomic_end(),
                     strand        => undef
                    };


  ## look up ontology info if non cached
  $data->{ontol}->{'polypeptide'} = $self->fetchSO('polypeptide')
    unless exists $data->{ontol}->{'polypeptide'};

  $gaFeature->{featureType} = $data->{ontol}->{'polypeptide'};


  ## what is interesting here?
  $gaFeature->{attributes} = { version => $feat->version(),
                               created => $feat->created_date(),
                               updated => $feat->modified_date()
                             };

  return $gaFeature;
}


## Look up onology information
## Return GA4GH OntologyTerm
sub fetchSO{

  my $self = shift;
  my $type = shift;

  my $onta = $self->context->model('Registry')->get_adaptor('Multi', 'Ontology', 'OntologyTerm');
  my $ont  = $onta->fetch_all_by_name($type); 

  my $ontologyTerm = { id            => $ont->[0]->accession(),
                       term          => $ont->[0]->name(), 
                       sourceName    => $ont->[0]->ontology(),
                       sourceVersion => undef
                     };

  return $ontologyTerm;
}

## create temp feature set name from current genebuild version
## replace with GA4GH id when format available
## TARK integration needed
sub getSet{

  my $self = shift;

  my $core_ad = $self->context->model('Registry')->get_DBAdaptor('homo_sapiens', 'Core'  );

  my $meta_ext_sth = $core_ad->dbc->db_handle->prepare(qq[ select meta_value from meta where meta_key = 'genebuild.id']);
  $meta_ext_sth->execute();
  my $genebuild_version = $meta_ext_sth->fetchall_arrayref();

  return $genebuild_version->[0]->[0];
}

sub sort_features{

  my $feat = shift;

  my @sorted = sort {  $a->{start} <=> $b->{start} } @{$feat};

  return \@sorted;
}

1;
