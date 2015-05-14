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

use Data::Dumper;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');

our $species = 'homo_sapiens';

## Switch to VEP cache for speed?
## What features to return?

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub searchFeatures {

  my ($self, $data ) = @_; 

  ## check feature type (ontology terms) supported
  ## initially default to transcripts only
  my $type_problem;
  if (exists $data->{features} ){
    foreach my $ontolterm (@{$data->{features}}){
      print "Requesting feature $ontolterm->{name}\n";
      $type_problem = $ontolterm->{name}  unless $ontolterm->{name} eq 'transcript';
    }
  }
  $self->context->go( 'ReturnError', 'custom', [ ' Request for unsupported feature type ' . $type_problem ] )
    if defined $type_problem ; 


  my $features;
  my $nextToken;

  if( exists $data->{range}){ ## is a specific path defined?
    foreach my $segment ( @{$data->{range}} ){

      ($features, $nextToken) = $self->extractFeaturesBySegment($data, $segment);
      last if defined $nextToken;
    }
  }
  else{
    ## default initial values
    my $next_seq_start = 1; 
    my $next_pos_start = 1;

    ## or take from token
    ($next_seq_start, $data->{from}) = split/\_/, $data->{pageToken} if defined $data->{pageToken};
     
    ## fudge a segment
    my $segment;
    ### Don't expect ids yet    $segment->{start}->{base}->{sequenceID} 
    $segment->{start}->{base}->{referenceName} = $next_seq_start;
    $segment->{start}->{base}->{position}      = $data->{from};
    $segment->{length}                         = 100000;   ## what is a sensible default slice length??
    ($features, $nextToken) = $self->extractFeaturesBySegment($data, $segment);
  }

  ## FIX: may not have a nextToken
  return ({ features     => $features, 
            nextPageToken => $nextToken}); 

}


## extract data from db
sub extractFeaturesBySegment{

  my ($self, $data, $segment ) = @_;

  ## db stuff
  my $sla = $self->context->model('Registry')->get_adaptor($species, 'Core', 'Slice');
  my $tra = $self->context->model('Registry')->get_adaptor($species, 'Core', 'Transcript');


  my @features;
  my $nextPageToken;

  my $end = $segment->{start}->{base}->{position} + $segment->{length};
  my $location = "$segment->{start}->{base}->{referenceName}\:$segment->{start}->{base}->{position}\-$end";
  my $slice = $sla->fetch_by_toplevel_location( $location );
  print "Using slice for location : $location\n";
  my $transcripts = $tra->fetch_all_by_Slice( $slice ) ;
  
  my $count = 0;
  foreach my $tr (@{$transcripts}){

    ## if pageSize reached save next position & return
    if ($count == $data->{pageSize}){
      print "Next trans should be : ". $tr->stable_id() ."\n";
      $nextPageToken = $tr->seq_region_name() . "_" . $tr->seq_region_start();
      last;
    }
    
    next if exists $data->{pageToken} && $data->{from} > $tr->seq_region_start();
 
    print "Transcript $count: ". $tr->stable_id() ."\n";
    my $gafeat = $self->formatTranscript($tr, $data);
    push @features,  $gafeat;
    ## keep count for pageSize
    $count++;

  }
  return (\@features, $nextPageToken);
}


## support more feature types..
sub getFeature{

  my $self = shift;
  my $id   = shift;

  return $self->getTranscript($id);
}

## look up transcript by id
## FIX: to use UUID

sub getTranscript{

  my $self = shift;
  my $id   = shift;

  my $tra = $self->context->model('Registry')->get_adaptor($species, 'Core', 'Transcript');
  my $tr = $tra->fetch_by_stable_id( $id );

  $self->context->go( 'ReturnError', 'custom', [ ' Cannot find transcript feature for id ' . $id ] )
    unless defined $tr ;

  return ({ features => [$self->formatTranscript($tr) ]});

}

## turn ensembl transcript into GA4GH transcript

sub formatTranscript{

  my $self = shift;
  my $tr   = shift;
  my $data = shift;

  my $gaFeature;
  $gaFeature->{id}            = $tr->stable_id();
  $gaFeature->{featureSetId}  = 'placeholder_Ensembl79';

  my $segment;
  ## FIX: check graph stuff
  $segment->{start}->{base}->{referenceName} = $tr->seq_region_name();
  $segment->{start}->{base}->{position}      = $tr->seq_region_start();
#  $segment->{start}->{strand} #not needed?
  $segment->{length}                         = $tr->seq_region_end() - $tr->seq_region_start();

  push @{$gaFeature->{path}}, $segment;

  ## look up ontology info if non supplied
  $data->{ontol}->{transcript} = $self->fetchSO('transcript') 
    unless exists $data->{ontol}->{transcript};

  $gaFeature->{featureType} = { id     =>  $data->{ontol}->{transcript}->{id},
                                name   =>  $data->{ontol}->{transcript}->{name},
                                source =>  $data->{ontol}->{transcript}->{source} };

  ## what is interesting here?
  $gaFeature->{attributes} = { version => $tr->version(),
                               biotype => $tr->biotype(),
                               gene    => $tr->get_Gene()->display_id(),
                               created => $tr->created_date(),           
                               updated => $tr->modified_date() 	
                             };

  return $gaFeature;

}

## Look up onology information
sub fetchSO{

  my $self = shift;
  my $type = shift;

  my $onta = $self->context->model('Registry')->get_adaptor('Multi', 'Ontology', 'OntologyTerm');
  my $ont  = $onta->fetch_all_by_name($type); 

  my $ontologyTerm;
  $ontologyTerm->{id}     = $ont->[0]->accession();
  $ontologyTerm->{source} = $ont->[0]->ontology();
  $ontologyTerm->{name}   = $ont->[0]->name();

  return $ontologyTerm;
}

1;
