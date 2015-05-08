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

## proper ids (awaiting guidlines)

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
  $data->{required_features} = $self->extractRequired( $data->{features} ) if $data->{features}->[0];
  $data->{required_effects}  = $self->extractRequired( $data->{effects} )  if $data->{effects}->[0];
  print Dumper $data;
  return $self->searchVariantAnnotations_by_features( $data)
    if exists $data->{features}->[0];

  return $self->searchVariantAnnotations_by_region( $data)
    if exists $data->{start};

  warn "uncaught error in searchVariantAnnotations\n";

}

sub searchVariantAnnotations_by_features {

  my ($self, $data ) = @_;

  my $c = $self->context();

  my $tra = $c->model('Registry')->get_adaptor($data->{species}, 'Core', 'Transcript');

#  foreach my $req_feat (@{$data->{features}}){
  my $transcript = $tra->fetch_by_stable_id( $data->{features}->[0] );
  
  $c->go('ReturnError', 'custom', [" feature $data->{features}->[0] not found"])
   if !$transcript;
   
  my $slice = $transcript->feature_Slice();
  $c->go('ReturnError', 'custom', ["location for $data->{features}->[0] not available for this assembly"])
   if !$slice;

  return $self->extractVFbySlice($data, $slice);

}

sub searchVariantAnnotations_by_region {

  my ($self, $data ) = @_; 

  my $c = $self->context(); 

  my $sla = $c->model('Registry')->get_adaptor($data->{species}, 'Core',      'Slice');

  my $start = $data->{start};
  $start = $data->{pageToken} if defined $data->{pageToken};  
  
  my $location = "$data->{referenceName}\:$start\-$data->{end}";
  my $slice = $sla->fetch_by_toplevel_location($location);
  
  $c->go('ReturnError', 'custom', ["sequence $location available for this assembly"])
   if !$slice;

  return $self->extractVFbySlice($data, $slice);

}

## get VF by slice and return appropriate number + next token 

sub extractVFbySlice{

  my $self  = shift;
  my $data  = shift;
  my $slice = shift;

  my @response;
  my $count = 0;

  my $vfa = $self->context->model('Registry')->get_adaptor($data->{species}, 'Variation', 'VariationFeature');
  $vfa->db->include_failed_variations(0); ## don't extract multi-mapping variants

  my $vfs = $vfa->fetch_all_by_Slice($slice);

  my $next_pos; ## save this for paging
  foreach my $vf(@{$vfs}){

    ## use next variant location as page token
    if ($count == $data->{pageSize}){
      $next_pos = $vf->seq_region_start();
      last;
    }

    ## filter by variant name if required
    next if defined $data->{variantName} && $vf->{variation_name} ne  $data->{variantName};

    ## extract & format - may not get a result if filtering applied
    my $var_an = $self->get_annotation($vf, $data);
    if (defined $var_an && $var_an ne ''){
      push @response, $var_an if defined $var_an;
      $count++;
    }
  }
  
  return ({"variant_annotation"  => \@response,
           "nextPageToken"       => $next_pos});

}

## extact and check annotation for a single variant

sub get_annotation{

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
      next if $tva->is_reference();

      my $ga_annotation;

      $ga_annotation->{alternateBase}       = $tva->variation_feature_seq();

      ## do HGVS - only return if there is a value?
      $ga_annotation->{HGVSg} = $tva->hgvs_genomic(); 
      $ga_annotation->{HGVSc} = $tva->hgvs_transcript() || undef; 
      $ga_annotation->{HGVSp} = $tva->hgvs_protein() || undef;

      ## get consequences & impact
      my $ocs = $tva->get_all_OverlapConsequences();

      my $keep_cons = 1;                                  ## by default report all consequences 
      $keep_cons = 0 if defined $data->{effects}->[0];    ## unless specific terms requested

      foreach my $oc(@{$ocs}) {
        my $cons = $oc->SO_term();

        $keep_cons = 1  if $data->{required_effects}->{ $cons } ;        
        push @{$ga_annotation->{effects}}, $cons;
        $ga_annotation->{impact} = $oc->impact() ;
      }
      ## skip if consequences specified but not found
      next unless $keep_cons ==1;

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

     push @{$var_ann->{transcriptEffects}}, $ga_annotation;
   }

  ## add co-located
  }
  return $var_ann if exists $var_ann->{transcriptEffects}->[0];
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


## put required features/effects in a hash for look up

sub extractRequired{

  my $self = shift;
  my $req_list = shift;

  my $req_hash;
  foreach my $required ( @{$req_list} ){
    $req_hash->{$required} = 1; 
  }
  return $req_hash;
}


1;
