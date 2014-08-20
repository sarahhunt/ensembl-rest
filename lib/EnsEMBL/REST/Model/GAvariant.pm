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

our $tmp_directory = "/tmp/ensrest/GAFiles";

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
  $data->{pageToken} = $self->get_variantRequestIds($c, $data)  if $data->{pageToken} eq "";

  $c->go('ReturnError', ' no data available') unless exists $data->{pageToken};


  ## get the next range of ids for the token
  my ($ids, $next_token) = get_nextByToken($data->{pageToken}, $data->{maxResults});
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

## extract all ids for request and put in temp file to allow paging
sub get_variantRequestIds{

    my ($self, $c, $data) = @_;

    ## request id/ file name (is md5 on request better?)
    my $unique_id = get_unique_id();
 
    ## allow filtering by set or return of everything
    my $constraint = " ";
print Dumper $data->{variantSetIds}; 
    $constraint = " and find_in_set(\'$data->{variantSetIds}->[0]\', variation_feature.variation_set_id)>0 "
	if defined $data->{variantSetIds}->[0] && $data->{variantSetIds}->[0] > 0 ;
    print "Constraint is : $constraint\n";

    my $dbh = $c->model('Registry')->get_DBAdaptor($data->{species}, 'Variation');
    my $data_ext_sth = $dbh->prepare(qq[select variation_feature.variation_feature_id 
                                        from seq_region, variation_feature 
                                         where variation_feature.seq_region_id = seq_region.seq_region_id 
                                        and seq_region.name= ?
                                        and variation_feature.seq_region_start between ? and ?
                                        $constraint ]);


    $data_ext_sth->{'mysql_use_result'} = 1;
    print  localtime() . " seeking " . $data->{referenceName} .":". $data->{start} ."-". $data->{end} ."\n";
    $data_ext_sth->execute($data->{referenceName}, $data->{start}, $data->{end} ) ||die;


    ## create request specific temp file to handle paging
    my $filename =  $tmp_directory . "/GA_$unique_id\.bed";

    open my $out, ">$filename" ||die "Failed to open file ($filename) to write : $!\n";

    my $n   = 0;

    while(my $line = $data_ext_sth->fetchrow_arrayref()){	
	$n++;
	print $out "1\t$n\t$n\t$line->[0]\n";       ## faked BED

    }
    close $out || die "Failed to close tmp file : $!\n";

    $c->go('ReturnError', 'custom', ["No data found in the required region"]) if $n ==0;
    
    ## format file
    print  localtime() ." dumped - zipping\n"; 
    system("bgzip $filename ") == 0      || die "Failed to bgzip $filename: $!";
    print  localtime() ." zipped - indexing\n";
    system("tabix -p bed $filename\.gz ") == 0  || die "Failed to tabix $filename: $!";
    print  localtime() ." indexed\n";
    
    ## create first token
    my $token = "$unique_id\_1";

    return $token;

}

### Database it??
sub get_nextByToken{

    my $token      = shift;
    my $batch_size = shift;

    my @return_ids;

    my ($unique_id, $start_pos ) = split/\_/,$token;

    my $filename  = $tmp_directory . "/GA_$unique_id\.bed.gz";
    die  "$filename not found\n"   unless -e($filename) ;

    my $new_pos = $start_pos + $batch_size;

   
    my $c = 0;
    open my $tab, "tabix $filename 1:$start_pos\-$new_pos  |" 
	||die "Failed to extract vf from index file : $!\n";
    while(<$tab>){
	my @a = split;
	push @return_ids, $a[3];
	$c++;
    }

    ## increment to next in file for new token
    my $new_token;
    if( $c == $batch_size ){
	## how to handle end of file - need to not return at token if at the end (this doen't catch the last batch if it is full)

	$new_token = "$unique_id\_$new_pos";
    }
    return (\@return_ids, $new_token );
}

## is this OK - is md5 on request better?
sub get_unique_id{

    my $ug    = Data::UUID->new;
    my $uuid = $ug->create_str(); 
    return $uuid;
}




with 'EnsEMBL::REST::Role::Content';

__PACKAGE__->meta->make_immutable;

1;
