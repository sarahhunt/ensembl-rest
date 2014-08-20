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

package EnsEMBL::REST::Model::GAcallSet;

use Moose;
extends 'Catalyst::Model';
use Data::Dumper;


with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');


sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}

sub fetch_ga_callSet {

  my ($self, $data ) = @_;

  print  localtime() . " Got request\n"; 
  print Dumper $data;

  ## only supporting human at the moment; species not specified in request
  $data->{species} = 'homo_sapiens';

  ## format set ids if filtering by set required
  if(defined $data->{variantSetIds}->[0]){

    my %req_variantset;
    foreach my $set ( @{$data->{variantSetIds}} ){
      $req_variantset{$set} = 1; 
    }
    $data->{req_variantsets} = \%req_variantset;
  } 

  ## extract required sets
  my $callsets = $self->fetch_callSets($data);

  return  $callsets;
}




sub fetch_callSets{

  my $self = shift;
  my $data = shift;

  my $c = $self->context();
  
  my $dbh = $c->model('Registry')->get_DBAdaptor($data->{species}, 'Variation');

=head
  ## limit by sample name if required
  my $constraint = " ";
  $constraint = " and ind.name like \'%$data->{name}%\'  " if defined $data->{name} && $data->{name} =~/\w+/;

  my $indset_ext_stmt  = (qq[ select ind.individual_id, ind.name, ind.variation_set_id from individual ind 
                             where ind.variation_set_id is not null $constraint ]);
  warn "preping $indset_ext_stmt \n";
  my $indset_ext_sth  = $dbh->prepare($indset_ext_stmt);

  $indset_ext_sth->execute()||die;
  my $set_data = $indset_ext_sth->fetchall_arrayref();
=cut

  ## hack pending db update
  my $set_data;
  open my $ind_dump, "/tmp/ensrest/GAFiles/individual.dat" || die "failed to open individual - set file\n";
  while(<$ind_dump>){
    chomp;
    my @a = split/\t/;
    next if defined $data->{name} && $a[1] !~ /$data->{name}/i ;
    push @{$set_data}, \@a;
 }
 ## 

  my @callsets;
  foreach my $l(@{$set_data}){
    my @sets = split/\,/,$l->[2];
    foreach my $set (@sets){
      
      ## limit by variant set if required
      next if defined $data->{req_variantsets} && ! defined $data->{req_variantsets}->{$set} ;
      my $callset;
      $callset->{sampleId} = $l->[1];
      $callset->{id}       = $l->[1]; ##$set . "_" . $l->[0];        ##concaternation of set id & sample id
      $callset->{name}     = $l->[1]; ##$name{$set} . "_". $l->[1];  ##concaternation of set name & sample name
      push @callsets, $callset;
    }
  }
  print  localtime() . " responding\n";
  return (\@callsets);
  
}


1;
