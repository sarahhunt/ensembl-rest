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

our $species = 'homo_sapiens';
 
sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub searchVariantAnnotations {
 
  my ($self, $data ) = @_;

  ## format look up lists if any specified
  $data->{required_features} = $self->extractRequired( $data->{feature_ids}, 'features') if $data->{feature_ids}->[0];
  $data->{required_effects}  = $self->extractRequired( $data->{effects}, 'SO' )          if $data->{effects}->[0];

  if ( $data->{variantAnnotationSetId} =~/compliance/){
    return $self->searchVariantAnnotations_compliance($data);
  }
  else{
     return $self->searchVariantAnnotations_database($data);
  }
}

sub searchVariantAnnotations_database {

  my ($self, $data ) = @_;

  ## temp set id
  $data->{current_set} = $self->getSet();

  $self->context->go('ReturnError', 'custom', [" No annotations available for this set: " . $data->{variantAnnotationSetId} ])
    if defined $data->{variantAnnotationSetId}  && $data->{variantAnnotationSetId} ne $data->{current_set} && $data->{variantAnnotationSetId} ne 'Ensembl'; 


  ## loop over features if supplied
  return $self->searchVariantAnnotations_by_features( $data)
    if exists $data->{feature_ids}->[0];

  ## search by region otherwise
  return $self->searchVariantAnnotations_by_region( $data)
    if exists $data->{start};

  warn "uncaught error in searchVariantAnnotations\n";

}

## probably a bad idea to allow multiple features
##  - more complex for server & client in paging 

sub searchVariantAnnotations_by_features {

  my ($self, $data ) = @_;

  ## return values
  my @annotations;
  my $nextPageToken;

  my $c = $self->context();
  
  my $tra = $c->model('Registry')->get_adaptor($species, 'Core',      'Transcript');
  my $tva = $c->model('Registry')->get_adaptor($species, 'variation', 'TranscriptVariation');


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
  foreach my $req_feat ( @{$data->{feature_ids}} ){
#    print "Starting to look for feature $req_feat ok_trans is $ok_trans\n";
    $ok_trans = 1 if defined $current_trans && $req_feat eq $current_trans;
    next unless $ok_trans == 1;

    ## may have no annotations if a transcript & consequence specified
#    $self->context->go('ReturnError', 'custom', [" No annotations available for this feature type: " . $req_feat->{featureType}->{name} ])
#      unless $req_feat->{featureType}->{name} eq "transcript";

    my $transcript = $tra->fetch_by_stable_id( $req_feat );  
    $c->go('ReturnError', 'custom', [" feature $req_feat not found"])
      if !$transcript;

    my $tvs;
    if(exists $data->{required_effects}){
      my @cons_terms = (keys %{$data->{required_effects} });

      my $constraint = " tv.consequence_types IN (";
      foreach my $cons(@cons_terms){ $constraint .= "\"$cons\",";}
      $constraint =~ s/\,$/\)/;
      $constraint .= " and somatic = 1 ";
#      print "Limiting effects for feature to : $constraint\n";
      $tvs = $tva->fetch_all_by_Transcripts_with_constraint([$transcript], $constraint );	
    }
    else{
      $tvs = $tva->fetch_all_somatic_by_Transcripts([$transcript]);
    }


    next unless scalar(@{$tvs}) > 0; 

    ## create an annotation record at the variation_feature level
    foreach my $tv (@{$tvs}){

      my $pos = $tv->variation_feature->seq_region_start();
      ## skip if already returned
      next if defined $next_pos && $pos < $next_pos;

      if($running_total == $data->{pageSize}){
        $nextPageToken = $req_feat . "_" . $pos;
        last;
      }

      my $var_ann;

      ## get consequences for each alt allele
      my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
      foreach my $tva(@{$tvas}) {
        my $ga_annotation  = $self->formatTVA( $tva, $tv ) ;
        push @{$var_ann->{transcriptEffects}}, $ga_annotation if defined $ga_annotation->{impact};
      }
      ## don't count or store until TV available
      next unless exists $var_ann->{transcriptEffects};

      $running_total++;
      $var_ann->{variantId} = $tv->variation_feature->variation_name();
      $var_ann->{variantAnnotationSetId} = $data->{current_set};
      $var_ann->{created} = 'FIXME_release_date';

      ## add co-located
      my $coLocated = $self->getColocated( $tv->variation_feature );
      $var_ann->{coLocatedVariants} = $coLocated if exists $coLocated->[0];

      if( defined $tv->variation_feature->minor_allele() ) {    
        $var_ann->{info}  = {  "1KG_minor_allele"           =>  $tv->variation_feature->minor_allele(),
                               "1KG_minor_allele_frequency" =>  $tv->variation_feature->minor_allele_frequency()
                            };
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

  my $sla = $c->model('Registry')->get_adaptor($species, 'Core',      'Slice');

  my $start = $data->{start};
  $start = $data->{pageToken} if defined $data->{pageToken};  
  
  my $location = "$data->{referenceName}\:$start\-$data->{end}";
  my $slice = $sla->fetch_by_toplevel_location($location);

  $c->go('ReturnError', 'custom', ["sequence  location: $location not available in this assembly"])
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

  my $vfa = $self->context->model('Registry')->get_adaptor($species, 'Variation', 'VariationFeature');
  $vfa->db->include_failed_variations(0); ## don't extract multi-mapping variants

  my $vfs;
  if( exists $data->{required_effects}){
    my @cons_terms = (keys %{$data->{required_effects} });

    my $constraint; 
    foreach my $cons(@cons_terms){ 
      $constraint .= "vf.consequence_types like \"%$cons\%\" and ";}
      $constraint =~ s/and $//;
#    print "limiting over region with $constraint\n";
    $vfs = $vfa->fetch_all_by_Slice_constraint($slice, $constraint);
  }
  else{
     $vfs = $vfa->fetch_all_by_Slice($slice);
  }

  my $next_pos; ## save this for paging
  foreach my $vf(@{$vfs}){
     next unless $vf->{variation_name} =~/rs|COS/; ## Exclude non dbSNP/ COSMIC for now
#    warn "seeking annot for " . $vf->variation_name() . " count is $count\n";
    ## use next variant location as page token
    if ($count == $data->{pageSize}){
      $next_pos = $vf->seq_region_start();
      last;
    }

    ## filter by variant name if required
    next if defined $data->{variantName} && $vf->{variation_name} ne  $data->{variantName};

    ## extract & format - may not get a result if filtering applied
    my $var_an = $self->fetchByVF($vf, $data);
    if (exists $var_an->{variantId}){
      push @response, $var_an ;
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
  $var_ann->{annotationSetId} = $data->{current_set};
  $var_ann->{created} = 'FIXME_release_date';


  my $tvs =  $vf->get_all_TranscriptVariations();
  return undef unless scalar(@{$tvs} > 0);

  foreach my $tv (@{$tvs}){   

    ## check if a feature list was specified
    next if scalar @{$data->{feature_ids}}>0 && !exists $data->{required_features}->{ $tv->transcript()->stable_id()} ; 

    my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
    foreach my $tva(@{$tvas}) {

      my $ga_annotation  = $self->formatTVA( $tva, $tv ) ;
      push @{$var_ann->{transcriptEffects}}, $ga_annotation;
    }
  }

  ## add 1KG global MAF as illustrative info
  if( defined $vf->minor_allele() ) {
    $var_ann->{info}  = {  "1KG_minor_allele"           =>  $vf->minor_allele(),
                           "1KG_minor_allele_frequency" =>  $vf->minor_allele_frequency()
                        };
  }

  ## add co-located
  my $coLocated = $self->getColocated( $vf );
  $var_ann->{coLocatedVariants} = $coLocated if exists $coLocated->[0];

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

  $ga_annotation->{feature_id} = $tv->transcript()->stable_id();
  $ga_annotation->{id} = "id";

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
  $sift_analysis->{analysisResult} = $tva->sift_prediction();
  $sift_analysis->{analysisScore}  = $tva->sift_score();
 
  if (defined $sift_analysis->{analysisResult}){
    ## move to anotationset?  
    $sift_analysis->{analysis}     = { id          => 'placeholder',
                                       description => 'SIFT',
                                       software    => 'SIFT.5.2.2',
                                     };
    push @analysisResults , $sift_analysis;
  }    
      
  my $polyphen_analysis;
  $polyphen_analysis->{analysisResult} = $tva->polyphen_prediction();
  $polyphen_analysis->{analysisScore}  = $tva->polyphen_score();

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
    $req_hash->{$required}         = 1 if $type eq 'features';
    $req_hash->{$required->{name}} = 1 if $type eq 'SO';
  }
  return $req_hash;
}

## create temp feature set name from curent db version
## replace with GA4GH id when format available
## needs to go in a utils
sub getSet{

  my $self = shift;

  my $var_ad   = $self->context->model('Registry')->get_DBAdaptor($species, 'variation');
  my $var_meta = $var_ad->get_MetaContainer();
  my $version  = $var_meta->schema_version();

  my $set = "Ensembl_" . $version; 

  return $set;
}

## get variants at same location
## add method to VFad for speed

sub getColocated{

  my $self   = shift;
  my $vf     = shift;

  my $vfa = $self->context->model('Registry')->get_adaptor($species, 'Variation', 'VariationFeature');

  my @colocatedNames;
  my $featureSlice = $vf->feature_Slice;

  my @colocated = (@{$vfa->fetch_all_by_Slice($featureSlice)}, @{$vfa->fetch_all_somatic_by_Slice($featureSlice)});

  return undef unless scalar(@colocated) >1;

  foreach my $colocated (@colocated){
    my $name = $colocated->variation_name();
    push @colocatedNames, $name unless $name eq  $vf->variation_name();
  }

  return \@colocatedNames;
}

sub searchVariantAnnotations_compliance{

  my $self = shift;
  my $data = shift;

  my $variantSetId = (split/\:/,  $data->{variantAnnotationSetId})[1];
 
  ##VCF collection object for the compliance annotation set
  $data->{vcf_collection} =  $self->context->model('ga4gh::ga4gh_utils')->fetch_VCFcollection_by_id( $variantSetId );
  $self->context()->go( 'ReturnError', 'custom', [ " Failed to find the specified variantSetId"])
    unless defined $data->{vcf_collection}; 

  ## look up filename from vcf collection object
  my $file  =  $data->{vcf_collection}->filename_template();

  $file =~ s/\#\#\#CHR\#\#\#/$data->{referenceName}/;

  # return these ordered by position for simple pagination
  my @var_ann;
  my $nextPageToken;
  $data->{pageToken} = $data->{start} unless defined $data->{pageToken}  && $data->{pageToken}  =~/\d+/;

  ## exits here if unsupported chromosome requested
  return (\@var_ann, $nextPageToken) unless -e $file;
 
  my $parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open( $file ) || die "Failed to get parser : $!\n";

  ## find column headers
  my $meta = $parser->get_metadata_by_pragma("INFO");
  my @headers;

  foreach my $m(@{$meta}){
    next unless $m->{Description} =~/Functional annotations/;
    my $headers = (split/\'/, $m->{Description})[1]; 
    @headers = split/\|/, $headers;
  }

  $parser->seek($data->{referenceName}, $data->{pageToken}, $data->{end});

  my $n = 0;

  while($n < ($data->{pageSize} + 1)){

    my $got_something = $parser->next();
    last if $got_something ==0;

    my $name = $parser->get_IDs->[0];

    if ($n == $data->{pageSize} ){
      ## batch complete 
      ## save next position for new page token
      $nextPageToken = $parser->get_raw_start -1;
      last;
    }

    $n++;
    my $variation_hash;
    $name = $parser->get_seqname ."_". $parser->get_raw_start if $name eq "." ; ## create placeholder name 
    $variation_hash->{id}                      = $data->{variantAnnotationSetId} .":". $variantSetId .":".$name; 
    $variation_hash->{variantAnnotationSetId}  = $data->{variantAnnotationSetId};
    $variation_hash->{variantId}               = $variantSetId .":".$name;

    $variation_hash->{created}                 = $data->{vcf_collection}->created || undef;
    $variation_hash->{updated}                 = $data->{vcf_collection}->updated || undef; 


    ## format array of transcriptEffects
    my $var_info = $parser->get_info();

    my @effects = split/\,/, $var_info->{ANN};
    foreach my $effect (@effects){

      my @annot = split/\|/, $effect;
      my %annotation;

      for (my $an =0; $an < scalar @annot; $an++ ){
        ## add filter for feature name if required
        $headers[$an] =~ s/^\s+|\s+$//g;
        $annotation{$headers[$an]} = $annot[$an] if defined $annot[$an];
      }

      ## check if a feature list was specified that this feature is required
      next if scalar @{$data->{feature_ids}}>0 && !exists $data->{required_features}->{ $annotation{Feature_ID} } ;

      my ($cdna_pos, $cdna_length ) = split /\//, $annotation{"cDNA.pos / cDNA.length"} 
        if defined $annotation{"cDNA.pos / cDNA.length"}; 

      my ($cds_pos, $cds_length )   = split /\//, $annotation{"CDS.pos / CDS.length"}
        if defined $annotation{"CDS.pos / CDS.length"} ;

      my ($aa_pos, $aa_length )     = split /\//, $annotation{"AA.pos / AA.length"}
        if defined $annotation{"AA.pos / AA.length"};

      my $transcriptEffect = { id              => "placeholder",
                               featureId       => $annotation{Feature_ID},
                               alternateBases  => $annotation{Allele},
                               impact          => $annotation{Annotation_Impact},
                               effects         => [{ id             => 'PLACEHOLDER',
                                                    name           => $annotation{Annotation},
                                                    ontologySource => 'Sequence Ontology'
                                                  }], 
                               HGVSg           => $annotation{"HGVS.g"} || undef,
                               HGVSc           => $annotation{"HGVS.c"} || undef,
                               HGVSp           => $annotation{"HGVS.p"} || undef,
                               cDNALocation    => {},
                               CDSLocation     => {},
                               proteinLocation => {},
                               analysisResults => []
                             };


    $transcriptEffect->{cDNALocation} = { overlapStart      => $cdna_pos -1,
                                          overlapEnd        => undef,
                                          referenceSequence => undef,
                                          alternateSequence => undef,
                                        } if defined $cdna_pos;

    $transcriptEffect->{CDSLocation}   = { overlapStart      => $cds_pos -1 ,
                                           overlapEnd        => undef,
                                           referenceSequence => undef,
                                           alternateSequence => undef,
                                         } if defined $cds_pos ;

    $transcriptEffect->{proteinLocation} =  { overlapStart      => $aa_pos -1,
                                              overlapEnd        => undef,
                                              referenceSequence => undef,
                                              alternateSequence => undef,
                                             } if defined $aa_pos;

      push @{$variation_hash->{transcriptEffects}}, $transcriptEffect;
    } 

    push @var_ann, $variation_hash if exists $variation_hash->{transcriptEffects}->[0];
  }

  $parser->close();

  return ( {"variantAnnotations"  => \@var_ann,
            "nextPageToken"       => $nextPageToken });
}

1;
