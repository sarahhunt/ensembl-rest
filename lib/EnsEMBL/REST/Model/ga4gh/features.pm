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

## FIX: pagination over multiple regions

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

## potential weirdness:
##  next token is: 
##     - numerically sorted transcript dbID in 'all' mode
##     - seq & pos start if in region mode

sub searchFeatures {

  my ($self, $data ) = @_; 

  $data->{current_set} = $self->getSet();

print localtime () . " Starting\n "; 
  ## check feature type (ontology terms) is supported
  ## initially default to transcripts only
  my $type_problem;
  if (exists $data->{features} ){
    foreach my $ontolterm (@{$data->{features}}){
      $type_problem = $ontolterm->{name}  unless $ontolterm->{name} eq 'transcript';
    }
  }
  $self->context->go( 'ReturnError', 'custom', [ ' Request for unsupported feature type ' . $type_problem ] )
    if defined $type_problem ; 


  my $features;
  my $nextToken;

  if( exists $data->{range}){
    ## extract specific region(s) if requested

    foreach my $segment ( @{$data->{range}} ){
      #print "using requested region\n";
      my $rangefeatures;
      ($rangefeatures, $nextToken) = $self->extractFeaturesBySegment($data, $segment);
      push @{$features}, @{$rangefeatures};
      last if defined $nextToken;
    }
  }
  else{
    ## page through everything in the database
    ($features, $nextToken) = $self->extractAllFeatures($data);
  }
print localtime () . " Ending\n ";
  ## FIX: may not have a nextToken
  return ({ features     => $features, 
            nextPageToken => $nextToken});

}

sub extractAllFeatures{

  my $self = shift;
  my $data = shift;


  ## take from db in transcript_id order if all requested
  my $tra = $self->context->model('Registry')->get_adaptor($species, 'Core', 'Transcript');

  my $get_id_sth = $tra->prepare(qq[ select transcript_id 
                                     from transcript
                                     where transcript_id >? 
                                     order by transcript_id
                                     limit ?    
                                  ]);

  my $next_trans = $data->{pageToken} || 0;
  my $limit      = $data->{pageSize};
 

  $get_id_sth->execute($next_trans, $limit );        
  my $trans_dbID =  $get_id_sth->fetchall_arrayref();

  my @featurelist;
  foreach my $dbID (@{$trans_dbID}){ 
    push @featurelist, $dbID->[0];
  }

  ## save last db id to start from next time
  my $last = pop @{$trans_dbID};
  my $nextToken = $last->[0];

  my $transcripts = $tra->fetch_all_by_dbID_list(\@featurelist);

  my $features; 
  foreach my $transcript (@{$transcripts}){
    my $gafeature =  $self->formatTranscript($transcript, $data) ;
    push @{$features}, $gafeature;
  }

  return ($features,  $nextToken); 

}



## extract data from db by region
sub extractFeaturesBySegment{

  my ($self, $data, $segment ) = @_;

  my $sla = $self->context->model('Registry')->get_adaptor($species, 'Core', 'Slice');
  my $tra = $self->context->model('Registry')->get_adaptor($species, 'Core', 'Transcript');


  my @features;
  my $nextPageToken;
  my $count = 0;

  ##fix paging for multiple regions if required
  my ($next_seq_start, $next_seq_pos ) = split/\_/, $data->{pageToken} if defined $data->{pageToken};

  my $end = $segment->{start}->{base}->{position} + $segment->{length};
  my $location = "$segment->{start}->{base}->{referenceName}\:$segment->{start}->{base}->{position}\-$end";
  my $slice = $sla->fetch_by_toplevel_location( $location );
  #print "Using slice for location : $location\n";
  my $transcripts = $tra->fetch_all_by_Slice( $slice ) ;

  foreach my $tr (@{$transcripts}){

    ## if pageSize reached save next position & return
    if ($count == $data->{pageSize}){ 
      $nextPageToken = $tr->seq_region_name() . "_" . $tr->seq_region_start();
      last;
    }
    
    next if exists $data->{pageToken} &&  $next_seq_pos > $tr->seq_region_start();
 
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
  
  my $data;
  $data->{current_set} = $self->getSet();

  return $self->getTranscript($id, $data);
}

## look up transcript by id
## FIX: to use UUID

sub getTranscript{

  my $self = shift;
  my $id   = shift;
  my $data = shift;

  my $tra = $self->context->model('Registry')->get_adaptor($species, 'Core', 'Transcript');
  my $tr = $tra->fetch_by_stable_id( $id, $data );

  $self->context->go( 'ReturnError', 'custom', [ ' Cannot find transcript feature for id ' . $id ] )
    unless defined $tr ;

  return ({ features => [$self->formatTranscript($tr, $data) ]});

}

## turn ensembl transcript into GA4GH transcript

sub formatTranscript{

  my $self = shift;
  my $tr   = shift;
  my $data = shift;

  my $gaFeature;
  $gaFeature->{id}            = $tr->stable_id();
  $gaFeature->{featureSetId}  = $data->{current_set};

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
                               updated => $tr->modified_date(),
                               source  => $tr->source()
                             };

  $gaFeature->{attributes}->{external_name} = $tr->external_name() if defined $tr->external_name(); 

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

## create temp feature set name from curent db version
## replace with GA4GH id when format available
sub getSet{

  my $self = shift;

  my $var_ad   = $self->context->model('Registry')->get_DBAdaptor($species, 'variation');
  my $var_meta = $var_ad->get_MetaContainer();
  my $version  = $var_meta->schema_version();

  my $set = "Ensembl_" . $version; 

  return $set;
}

1;
