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

package EnsEMBL::REST::Model::ga4gh::variantannotations;

use Moose;
extends 'Catalyst::Model';

use Bio::DB::Fasta;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Time::HiRes qw(gettimeofday);
use Data::Dumper;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');

## speed up
##  - fasta look up
##  - order trans var to not call unnecessary methods
##  - Is GFF faster than db?

## proper ids (awaiting guidelines)

## handle genes as well as transcripts

## paging forever issue if rare consequence - seek by constraint?

## check effects are valid else page for ever
 
sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub searchVariantAnnotations {

  my ($self, $data ) = @_;

  ## set species
  $data->{species} = "homo_sapiens";

  ## format look up lists if any specified
  $data->{required_features} = $self->extractRequired( $data->{features}, 'features') if $data->{features}->[0];
  $data->{required_effects}  = $self->extractRequired( $data->{effects}, 'SO' )       if $data->{effects}->[0];
  #print "Input: "; print Dumper $data;
  return $self->searchVariantAnnotations_by_features( $data)
    if exists $data->{features}->[0];

  return $self->searchVariantAnnotations_by_region( $data)
    if exists $data->{start};

  warn "uncaught error in searchVariantAnnotations\n";

}

## probably a bad idea to allow multiple features
##  - more complex for server & client in paging 

sub searchVariantAnnotations_by_features {

  my ($self, $data ) = @_;
print "In searchVariantAnnotations_by_features with page token $data->{pageToken}\n";
  ## return values
  my @annotations;
  my $nextPageToken;

  my $c = $self->context();
  
  my $tra = $c->model('Registry')->get_adaptor($data->{species}, 'Core',      'Transcript');
  my $tva = $c->model('Registry')->get_adaptor($data->{species}, 'variation', 'TranscriptVariation');


  ## hacky paging
  my $running_total = 0;

  my ($current_trans, $next_pos);
  my $ok_trans = 1;
  if (exists $data->{pageToken}){
     ($current_trans, $next_pos) = split/\_/, $data->{pageToken}; 
     $ok_trans = 0 ;
  }

  ## extract one transcripts worth at once for paging
  ## should records be merged where transcripts overlap??
  foreach my $req_feat ( @{$data->{features}} ){
    print "Starting to look for feature $req_feat->{id} ok_trans is $ok_trans\n";
    $ok_trans = 1 if defined $current_trans && $req_feat->{id} eq $current_trans;
    next unless $ok_trans == 1;

    ## FIX: silent fail 
    next unless $req_feat->{featureType}->{name} eq "transcript";
      
    my $transcript = $tra->fetch_by_stable_id( $req_feat->{id} );  
    $c->go('ReturnError', 'custom', [" feature $req_feat->{id} not found"])
      if !$transcript;

    my $tvs;
    if(exists $data->{required_effects}){
      my @cons_terms = (keys %{$data->{required_effects} });
      $tvs = $tva->fetch_all_by_Transcripts_SO_terms([$transcript, \@cons_terms ]);	
    }
    else{
      $tvs = $tva->fetch_all_by_Transcripts([$transcript]);
    }

    ## create an annotation record at the variation_feature level
    foreach my $tv (@{$tvs}){

      my $pos = $tv->variation_feature->seq_region_start();
      next if defined $next_pos && $pos < $next_pos;
      last if defined $nextPageToken ;      

      my $var_ann;
      $var_ann->{variantId} = $tv->variation_feature->variation_name();
print "Checking annot for $var_ann->{variantId}\n";
      $var_ann->{annotationSetId} = 'Ensembl_79';
      $var_ann->{created} = 'FIXME_release_date';

      my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
      foreach my $tva(@{$tvas}) {
        if($running_total == $data->{pageSize}){
          $nextPageToken = $req_feat->{id} . "_" . $pos;
          last;
        } 
        my $ga_annotation  = $self->formatTVA( $tva, $tv ) ;
        push @{$var_ann->{transcriptEffects}}, $ga_annotation;
        $running_total++;
      }
     push @annotations, $var_ann;
    }
  }


  ## may have no annotations if a transcript & consequence specified
  $self->context->go('ReturnError', 'custom', [" No annotations found for this query "])
      if $running_total ==0;

  return ({"variant_annotation"  => \@annotations,
           "nextPageToken"       => $nextPageToken });
}



sub searchVariantAnnotations_by_region {

  my ($self, $data ) = @_; 

  my $c = $self->context(); 

  my $sla = $c->model('Registry')->get_adaptor($data->{species}, 'Core',      'Slice');

  my $start = $data->{start};
  $start = $data->{pageToken} if defined $data->{pageToken};  
  
  my $location = "$data->{referenceName}\:$start\-$data->{end}";
  my $slice = $sla->fetch_by_toplevel_location($location);
#  print "Looking by region : $location\n";
  $c->go('ReturnError', 'custom', ["sequence $location available for this assembly"])
   if !$slice;

  my ($annotations, $nextPageToken) =  $self->extractVFbySlice($data, $slice);

  return ({"variant_annotation"  => $annotations,
           "nextPageToken"       => $nextPageToken });
}

## get VF by slice and return appropriate number + next token 

sub extractVFbySlice{

  my $self  = shift;
  my $data  = shift;
  my $slice = shift;
  my $count = shift; ## handle multiple regions in one response. Is this a good idea?

  $count ||= 0;

  my @response;

  my $vfa = $self->context->model('Registry')->get_adaptor($data->{species}, 'Variation', 'VariationFeature');
  $vfa->db->include_failed_variations(0); ## don't extract multi-mapping variants

  my $vfs;
  if( exists $data->{required_effects}){
    my @cons_terms = (keys %{$data->{required_effects} });
print "Sending in effects : " ; print Dumper @cons_terms;
    #$vfs = $vfa->fetch_all_by_Slice_SO_terms($slice, \@cons_terms);
    my $constraint = " vf.consequence_types IN (";
    foreach my $cons(@cons_terms){ $constraint .= "\"$cons\",";}
    $constraint =~ s/\,$/\)/;
    $vfs = $vfa->fetch_all_by_Slice_constraint($slice, $constraint);
  }
  else{
     $vfs = $vfa->fetch_all_by_Slice($slice);
  }
  my $next_pos; ## save this for paging
  foreach my $vf(@{$vfs}){
    warn "seeking annot for " . $vf->variation_name() . "\n";
    ## use next variant location as page token
    if ($count == $data->{pageSize}){
      $next_pos = $vf->seq_region_start();
      last;
    }

    ## filter by variant name if required
    next if defined $data->{variantName} && $vf->{variation_name} ne  $data->{variantName};

    ## extract & format - may not get a result if filtering applied
    my $var_an = $self->fetchByVF($vf, $data);
    if (defined $var_an && $var_an ne ''){
      push @response, $var_an if defined $var_an;
      $count++;
    }
  }

  return ( \@response, $next_pos );

}

## extact and check annotation for a single variant

sub fetchByVF{

  my $self = shift;
  my $vf   = shift;
  my $data = shift;

  my $var_ann;
  $var_ann->{variantId} = $vf->variation_name();
  $var_ann->{annotationSetId} = 'Ensembl_79';
  $var_ann->{created} = 'FIXME_release_date';

  my $tvs =  $vf->get_all_TranscriptVariations();

  foreach my $tv (@{$tvs}){   


    ## check if a feature list was specified
    next if scalar @{$data->{features}}>0 && !exists $data->{required_features}->{ $tv->transcript()->stable_id()} ; 

    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
    foreach my $tva(@{$tvas}) {

      my $ga_annotation  = $self->formatTVA( $tva, $tv ) ;
      push @{$var_ann->{transcriptEffects}}, $ga_annotation;
   }

  ## add co-located
  }
  return $var_ann if exists $var_ann->{transcriptEffects}->[0];
}


sub formatTVA{

  my $self = shift;
  my $tva  = shift;
  my $tv   = shift;

  my $ga_annotation;

  $ga_annotation->{alternateBase}       = $tva->variation_feature_seq();

  ## do HGVS - only return if there is a value?
  $ga_annotation->{HGVSg} = $tva->hgvs_genomic(); 
  $ga_annotation->{HGVSc} = $tva->hgvs_transcript() || undef; 
  $ga_annotation->{HGVSp} = $tva->hgvs_protein()    || undef;

  ## get consequences & impact
  my $ocs = $tva->get_all_OverlapConsequences();

  foreach my $oc(@{$ocs}) {
    my $term = $oc->SO_term();
    my $acc  = $oc->SO_accession();
    my $ontolTerm = { id             => $acc,
                      name           => $term,
                      ontologySource => 'Sequence Ontology'      
                     };

    push @{$ga_annotation->{effects}}, $ontolTerm;
    $ga_annotation->{impact} = $oc->impact() ;
  }

  my $cdna_start = $tv->cdna_start();
  if( defined $cdna_start ){
    $ga_annotation->{cDNALocation}->{overlapStart}   =   $cdna_start  - 2;
    $ga_annotation->{cDNALocation}->{overlapEnd}     =   $tv->cdna_end() - 1;
  }
 
  my $codon = $tva->codon() ;
  if( defined $codon ){
    $ga_annotation->{cdsLocation}->{referenceSequence} = $tv->get_reference_TranscriptVariationAllele->codon();
    $ga_annotation->{cdsLocation}->{alternateSequence} = $codon;
    $ga_annotation->{cdsLocation}->{overlapStart}      = $tv->cds_start() -2;
    $ga_annotation->{cdsLocation}->{overlapEnd}        = $tv->cds_end() -1;
  }


  my $peptide =  $tva->peptide();
  if( defined $peptide ){
    $ga_annotation->{proteinLocation}->{referenceSequence} = $tv->get_reference_TranscriptVariationAllele->peptide();
    $ga_annotation->{proteinLocation}->{alternateSequence} = $peptide;
    $ga_annotation->{proteinLocation}->{overlapStart}      = $tv->translation_start()  - 2;
    $ga_annotation->{proteinLocation}->{overlapEnd}        = $tv->translation_end()  -1;

    ## Extract protein impact information for missense variants
    my $protein_impact = $self->protein_impact($tva); 
    $ga_annotation->{analysisResults} = $protein_impact if @{$protein_impact} >0; 
  }

  return $ga_annotation;
}



## extract & format sift and polyphen results
sub protein_impact{

  my $self = shift; 
  my $tva  = shift;

  my @analysisResults;

  my $sift_analysis;
  $sift_analysis->{analysisResult} = $tva->sift_prediction() || undef;
  $sift_analysis->{analysisScore}  = $tva->sift_score()      || undef;

  if (defined $sift_analysis->{analysisResult}){
    ## move to anotationset?  
    $sift_analysis->{analysis}     = { id          => 'placeholder',
                                       description => 'SIFT',
                                       software    => 'SIFT.5.2.2',
                                     };
    push @analysisResults , $sift_analysis;
  }    
      
  my $polyphen_analysis;
  $polyphen_analysis->{analysisResult} = $tva->polyphen_prediction() || undef;
  $polyphen_analysis->{analysisScore}  = $tva->polyphen_prediction() || undef;

  if (defined $polyphen_analysis->{analysisResult}){
    $polyphen_analysis->{analysis}     = { id          => 'placeholder',
                                           description => 'Polyphen',
                                           software    => 'Polyphen.2.2.2_r405',
                                          };
                                               
    push @analysisResults , $polyphen_analysis;
  }

  return \@analysisResults ;
}


## put required features ids/ SO accessions in a hash for look up

sub extractRequired{

  my $self     = shift;
  my $req_list = shift;
  my $type     = shift;
  my $req_hash;

  foreach my $required ( @{$req_list} ){
    $req_hash->{$required->{id}}   = 1 if $type eq 'feature';
    $req_hash->{$required->{name}} = 1 if $type eq 'SO';
  }
  return $req_hash;
}


1;
