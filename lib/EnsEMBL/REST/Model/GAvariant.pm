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
use Data::UUID;

with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');



sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub fetch_gavariant {

  my ($self, $data ) = @_;

  print  localtime() . " Got request\n"; 
#  print Dumper $data;

  ## only supporting human at the moment; species not specified in request
  $data->{species} = 'homo_sapiens';


  ## format sample names if filtering by sample required
  if(defined $data->{callSetIds}->[0]){
    my %req_samples;
    foreach my $sample ( @{$data->{callSetIds}} ){
      $req_samples{$sample} = 1; 
    }
    $data->{samples} = \%req_samples;
  }

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

sub fetch_by_region{

  my $self = shift;
  my $data = shift;

  my $c = $self->context();

  ## if this is a new request look up the ids required & set first token   
  $data->{pageToken} = $self->get_first_token($c, $data)  if $data->{pageToken} eq "";

  $c->go('ReturnError', ' no data available') unless exists $data->{pageToken};


  ## get the next range of ids for the token
  my ($ids, $next_token) = $self->get_next_by_token($data);
  ## exit if none found
  $c->go( 'ReturnError', 'custom', [ " No variants are available for this search token: $data->{pageToken}, batch: $data->{maxResults}" ] )  unless exists $ids->[0];

  my $vfa = $c->model('Registry')->get_adaptor($data->{species}, 'Variation', 'VariationFeature');
  
  ## get details for variationfeature ids
  my @var_response;

  foreach my $vfid( @{$ids}){
    my $vf = $vfa->fetch_by_dbID($vfid);
     
    $c->go('ReturnError', 'custom', ["No data for next variation in list ($vf)"]) if !$vf;;
      
    ## extract & format
    my $ga_var = $self->to_hash($vf, $data);
    push @var_response, $ga_var;
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

  my $variation_hash;

  $variation_hash->{variantSetId}   = "Ensembl_version";
  $variation_hash->{id}             = $vf->name();
  $variation_hash->{name}           = $vf->name();
 
  my @alleles = (split/\//, $vf->allele_string); 

  $variation_hash->{calls}          = $self->Genotypes($vf, \@alleles, $data);


  ## ?? Exclude multi-mapping, non-sorted var 
  $variation_hash->{referenceBases} = shift @alleles;
  $variation_hash->{alternateBases} = \@alleles;


  $variation_hash->{referenceName}  = $vf->seq_region_name();
  ## position is zero-based
  ## closed start of interval 
  $variation_hash->{start}          = $vf->seq_region_start() -1;                                           
  ## open end of interval 
  $variation_hash->{end}            = $vf->seq_region_end() ; 

  ## extract consequence SO terms
  $variation_hash->{consequences}   = $self->Consequences($vf);

  return $variation_hash;
}



## extract genotypes & apply filtering
sub Genotypes {
  my ($self, $vf, $alleles, $data) = @_;

  my %id;
  for (my $n=0; $n < scalar(@{$alleles}); $n++){
       $id{$alleles->[$n]} = $n;
  }  


  my @genotypes;
  my $genotypes = $vf->variation->get_all_IndividualGenotypes();
  foreach my $gen (@$genotypes) {

    ## filter by individual if required
    next if defined $data->{callSetIds}->[0] && 
      !defined $data->{samples}->{$gen->individual()->name()};

    push @genotypes, $self->gen_as_hash($gen, \%id);
  }
  return \@genotypes;
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

sub get_first_token{

  my ($self, $c, $data) = @_;

  my $unique_id = get_unique_id();

  my $token = $unique_id ."_". $data->{start};

  return $token;
}


## extract a batch of ids for request and hold new start pos in token string
sub get_next_by_token{

  my ($self, $data) = @_;
 
  my ($request_id, $batch_start) = (split/\_/,$data->{pageToken});


  ## allow filtering by set or return of everything
  my $constraint = " ";
  print Dumper $data->{variantSetIds}; 
  $constraint = " and find_in_set(\'$data->{variantSetIds}->[0]\', vf.variation_set_id)>0 "
    if defined $data->{variantSetIds}->[0] && $data->{variantSetIds}->[0] > 0 ;
  print "Constraint is : $constraint\n";

  my $limit = $data->{maxResults} + 1; ## is there anything to come back for in another batch

  my $dbh = $self->context()->model('Registry')->get_DBAdaptor($data->{species}, 'Variation');
  my $data_ext_sth = $dbh->prepare(qq[  select vf.variation_feature_id, vf.seq_region_start  
                                        from seq_region sr, variation_feature vf
                                        where vf.seq_region_id = sr.seq_region_id 
                                        and sr.name= ?
                                        and vf.seq_region_start between ? and ?
                                        $constraint
                                        limit $limit ]);


  $data_ext_sth->{'mysql_use_result'} = 1;
  print  localtime() . " seeking " . $data->{referenceName} .":". $batch_start ."-". $data->{end} ."\n";
  $data_ext_sth->execute($data->{referenceName}, $batch_start, $data->{end} ) ||die;

  my $n  = 0;
  my @object_ids;
  my $last_pos;   ## next batch start

  while(my $line = $data_ext_sth->fetchrow_arrayref()){	
    $n++;
    push @object_ids, $line->[0] if $n <= $data->{maxResults} ;
    $last_pos = $line->[1];
  }
  ## this should not happen
  $self->context()->go('ReturnError', 'custom', ["No data found in the required region"]) if $n ==0;
    

  ## create next token to include the position of next the variant if an extra one is available
  my $token = $request_id ."_" . $last_pos  if $n>$data->{maxResults} ;


  return (\@object_ids, $token);

}


## is this OK - is md5 on request better?
sub get_unique_id{

  my $ug   = Data::UUID->new;
  my $uuid = $ug->create_str(); 
  return $uuid;
}




with 'EnsEMBL::REST::Role::Content';

__PACKAGE__->meta->make_immutable;

1;
