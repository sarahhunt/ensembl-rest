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

package EnsEMBL::REST::Model::GAvariant;

use Moose;
extends 'Catalyst::Model';
use Data::Dumper;


with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');



sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub fetch_gavariant {

  my ($self, $data ) = @_;

  print  localtime() . " Got request\n"; 

  ## only supporting human at the moment; species not specified in request
  $data->{species} = 'homo_sapiens';

  ## need these ordered for paging
  @{$data->{variantSetIds}} = sort @{$data->{variantSetIds}};

  ## load callSet to variantSet link 
  $data = $self->get_set_info($data);

  ## exclude samples outside specified variantSet and variantSet with no required samples
  $data = $self->check_sample_info($data) if defined $data->{callSetIds}->[0] ;


  ## get variation data 
  if ( $data->{variantName} ne ""){
    ## if variant name supplied look up by name   
    return $self->fetch_by_name($data);      
  }

  else{
    ## pull out set by region if variant name not supplied
    return $self->fetch_by_region($data);
  }

}


sub get_set_info{

  my ($self, $data ) = @_;

  ## link sample names to variantSet for later filtering
  if(defined $data->{variantSetIds}->[0] && $data->{variantSetIds}->[0] >0){ ##use 0 for everything?

    foreach my $set (@{$data->{variantSetIds}}){
      $data->{required_set}->{$set} = 1;
    }

    ## hack pending db update
    open my $ind_dump, "/tmp/ensrest/GAFiles/individual.dat" || die "failed to open individual - set file\n";
    while(<$ind_dump>){
      chomp;
      my @a = split;
      $data->{sample2set}->{$a[1]} = $a[2] if  defined $data->{required_set}->{$a[2]};
    }
  }
 return $data;
}


## format sample names if filtering by sample required
## variantSet limitation takes presidence
sub check_sample_info{

  my ($self, $data ) = @_;

  my %req_samples; ## store sample <-> variantSet link
  my %req_sets;    ## if samples are specified, don't look in sets not containing them

  foreach my $sample ( @{$data->{callSetIds}} ){
    if (defined  $data->{sample2set}->{$sample}){
      $req_samples{$sample} = 1;
      $req_sets{ $data->{sample2set}->{$sample} } = 1;
    }
  }
 
  ## exit if the samples & sets are incompatible
  $self->context()->go( 'ReturnError', 'custom', [ " The specified callSets are not available in the specified variantSets"])
    unless scalar(keys %req_samples) >0;

  $data->{samples} = \%req_samples;

  ## reset variantSets to only those with samples
  $data->{required_set} = \%req_sets;
  
  return $data;
}


sub fetch_by_region{

  my $self = shift;
  my $data = shift;

  my $c = $self->context();

  ## if this is a new request set first token to start of region and set 0  
  $data->{pageToken} = $data->{start} . "_0" if $data->{pageToken} eq "";

  ## get the next range of ids for the token
  my $ids = $self->get_next_by_token($data);
  ## exit if none found
  $c->go( 'ReturnError', 'custom', [ " No variants are available for this search token: $data->{pageToken}, batch: $data->{maxResults}" ] )  
     unless scalar(@{$ids}) >0;

  my $vfa = $c->model('Registry')->get_adaptor($data->{species}, 'Variation', 'VariationFeature');
  
  ## get details for variationfeature ids
  my @var_response;
  my $gavariant_count = 0;
  my $next_token;
  my ($last_pos, $last_set ) = split/\_/,  $data->{pageToken} ;

  foreach my $vfid( @{$ids}){

    last if $gavariant_count == $data->{maxResults};

    my $vf = $vfa->fetch_by_dbID($vfid);
     
    $c->go('ReturnError', 'custom', ["No data for next variation in list ($vf)"]) if !$vf;;
      
    ## extract & format - may return many GAvariants from different variantsets for a variant
    my $ga_vars = $self->to_hash($vf, $data );


    foreach my $ga_var (@{$ga_vars}){ 
      
      ## last batch may have ended mid-variant - skip any which have been returned
      next if $ga_var->{start} == $last_pos &&  $ga_var->{variantSetId} <= $last_set;

      $gavariant_count++;
      push @var_response, $ga_var;
      
      if($gavariant_count == $data->{maxResults}){
        ## set next token to position and set of last reported result  
        $next_token =  $ga_var->{start} ."_". $ga_var->{variantSetId} ;
        last;
      }
    }
  }
  push @var_response, ["nextPageToken" => $next_token];
  print  localtime() . " responding\n";
  return (\@var_response);
  
}

sub fetch_by_name{

  my $self = shift;
  my $data = shift;

  my $c = $self->context();

  ## look up variant info
  my $vfa = $c->model('Registry')->get_adaptor($data->{species}, 'Variation', 'VariationFeature');
  my $va  = $c->model('Registry')->get_adaptor($data->{species}, 'Variation', 'Variation');


  my @var_response;

  my $var  =  $va->fetch_by_name( $data->{variantName} );
  $c->go('ReturnError', 'custom', ["No data for a variation named $data->{variantName}" ]) unless defined $var;

  my $vfs  =  $vfa->fetch_all_by_Variation($var);
  $c->go('ReturnError', 'custom', ["No location for a variation named $data->{variantName}" ]) unless defined $vfs->[0];

  ## check position
  my $found_in_region = 0;

  foreach my $vf (@{$vfs}){
    ##print "Checking vf location: ".  $vf->slice->seq_region_name() .":".$vf->seq_region_start() . "  v $data->{referenceName} : $data->{start}  & $data->{end}\n"; 
    if( $vf->slice->seq_region_name() eq $data->{referenceName} &&
         $vf->seq_region_start() >= $data->{start}  &&
         $vf->seq_region_start() <= $data->{end}){

       ## extract & format
       my $ga_var = $self->to_hash($vf, $data);
       push @var_response, $ga_var;
       $found_in_region = 1 ;
    }
  }
  $c->go('ReturnError', 'custom', ["No data for variation named $data->{variantName} in required region"]) unless $found_in_region == 1 ;
  print  localtime() . " responding\n"; 
  return \@var_response;
}



sub to_hash {

  my ($self, $vf, $data) = @_;

  my @variants;

  ## get non-set specific data to share between GA variants

  my @alleles = (split/\//, $vf->allele_string); 

  my $genotype_calls = $self->Genotypes($vf, \@alleles, $data);

  ## ?? Exclude multi-mapping, non-sorted var [calc set specific here?] 
  my $referenceBases = shift @alleles;
  my $alternateBases = \@alleles;

  ## extract consequence SO terms
  my $consequences   = $self->Consequences($vf);

  ## loop over sets to divide up genotypes
  foreach my $set_required (@{$data->{variantSetIds}}){

    ## check there are genotypes to return
    next unless exists $genotype_calls->{$set_required}->[0];

    my $variation_hash;

    $variation_hash->{variantSetId}   = $set_required ;
    $variation_hash->{id}             = $vf->name();
    $variation_hash->{name}           = $vf->name();
  
    $variation_hash->{calls}          = $genotype_calls->{$set_required};
 
    $variation_hash->{referenceBases} = $referenceBases;
    $variation_hash->{alternateBases} = $alternateBases;
    $variation_hash->{consequences}   = $consequences;

    $variation_hash->{referenceName}  = $vf->seq_region_name();

    ## position is zero-based
    ## closed start of interval 
    $variation_hash->{start}          = $vf->seq_region_start() -1;                                           
    ## open end of interval 
    $variation_hash->{end}            = $vf->seq_region_end() ;

    push @variants, $variation_hash;
  }
  return \@variants;
}



## extract genotypes & apply filtering
sub Genotypes {
  my ($self, $vf, $alleles, $data) = @_;

  my %id;
  for (my $n=0; $n < scalar(@{$alleles}); $n++){
       $id{$alleles->[$n]} = $n;
  }  


  my %genotypes;
  my $genotypes = $vf->variation->get_all_IndividualGenotypes();
  foreach my $gen (@$genotypes) {

    ## filter by individual if required
    next if defined $data->{callSetIds}->[0] && 
      !defined $data->{samples}->{$gen->individual()->name()};

    ## store genotypes by variationSetId
    push @{$genotypes{ $data->{sample2set}->{ $gen->individual()->name() }} }, $self->gen_as_hash($gen, \%id)
      if defined $data->{sample2set}->{ $gen->individual()->name() };
  }
  return \%genotypes;
}

## format and code genotypes
sub gen_as_hash {
  my ($self, $gen, $id) = @_;

  my $gen_hash;
  $gen_hash->{callSetId}    = $gen->individual->name();
  $gen_hash->{callSetName}  = $gen->individual->name();

  my $genos =  $gen->genotype();
  foreach my $g( @{$genos}){
    push @{$gen_hash->{genotype}}, $id->{$g};
  }

  return $gen_hash;
}


## extract all SO terms for the variants consequences & return string
sub Consequences{ 

  my $self = shift;
  my $vf   = shift;

  my %cons_list;

  my $consequences = $vf->get_all_OverlapConsequences();
  
  return undef unless defined $consequences->[0];

  foreach my $cons(@{$consequences}){
     $cons_list{ $cons->SO_term() } = 1;
   }
   my $cons_string = join(",", keys %cons_list);  

  return $cons_string;

}


## extract a batch of ids for request and hold new start pos in token string
sub get_next_by_token{

  my ($self, $data) = @_;
 
  ## These are the last seq start and set reported
  my ($batch_start, $set_start) = (split/\_/,$data->{pageToken});

  ## filter by set (or return of everything?) 
  my $constraint = " and (find_in_set(\'$data->{variantSetIds}->[0]\', vf.variation_set_id)>0 ";
  foreach my $n (1 .. scalar(@{$data->{variantSetIds}})-1 ){
     $constraint .= " or find_in_set(\'$data->{variantSetIds}->[$n]\', vf.variation_set_id)>0 "
  }
  $constraint .= ") ";
  

  my $limit = $data->{maxResults} + 50; ## is there anything to come back for in another batch
  ## allowing extra here as variant in set may have no data for a chosen individual

  my $dbh = $self->context()->model('Registry')->get_DBAdaptor($data->{species}, 'Variation');

  $self->context()->go('ReturnError', 'custom', ["Failed to connect to database"]) unless defined $dbh;

  my $data_ext_sth = $dbh->prepare(qq[  select vf.variation_feature_id, vf.seq_region_start  
                                        from seq_region sr, variation_feature vf
                                        where vf.seq_region_id = sr.seq_region_id 
                                        and sr.name= ?
                                        and vf.seq_region_start between ? and ?
                                        $constraint
                                        limit $limit ]);


  $data_ext_sth->{'mysql_use_result'} = 1;
  print  localtime() . " seeking ids:" . $data->{referenceName} .":". $batch_start ."-". $data->{end} . " limit $limit\n";
  $data_ext_sth->execute($data->{referenceName}, $batch_start, $data->{end} ) ||die;

  my $n  = 0;
  my @object_ids;

  while(my $line = $data_ext_sth->fetchrow_arrayref()){	
    $n++;
    push @object_ids, $line->[0] ;
  }

  ## this should not happen
  $self->context()->go('ReturnError', 'custom', ["No data found in the required region"]) if $n ==0;
  

  return \@object_ids;

}



with 'EnsEMBL::REST::Role::Content';

__PACKAGE__->meta->make_immutable;

1;
