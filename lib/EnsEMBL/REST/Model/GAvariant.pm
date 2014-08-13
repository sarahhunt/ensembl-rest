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
  my ($self, $data) = @_;

  my $c = $self->context();
  my $species = 'homo_sapiens'; ## GA4GH only for human, species not specified in request

  ## is filtering by sample required
  if(defined $data->{callSetIds}->[0]){
    my %req_samples;
    foreach my $sample ( @{$data->{callSetIds}} ){
      $req_samples{$sample} = 1; 
    }
    $data->{samples} = \%req_samples;
  
print Dumper $data;

  }
  my $sla = $c->model('Registry')->get_adaptor($species, 'Core', 'Slice');
  my $vfa = $c->model('Registry')->get_adaptor($species, 'Variation', 'VariationFeature');

  my $location = "$data->{referenceName}\:$data->{start}\-$data->{end}";
  my $slice = $sla->fetch_by_toplevel_location($location);
  if (!$slice) {
    $c->go('ReturnError', 'custom', ["sequence $location not found for $species"]);
  }

  my @var_response;

  my $vfs = $vfa->fetch_all_by_Slice($slice);

  foreach my $vf(@{$vfs}){

    ## filter by variant name if required
    next if defined $data->{variantName} && $vf->{variation_name} ne  $data->{variantName};
   
    ## extract & format
    my $ga_var = $self->to_hash($vf, $data);
    push @var_response, $ga_var;
  }

  return (\@var_response);
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

with 'EnsEMBL::REST::Role::Content';

__PACKAGE__->meta->make_immutable;

1;
