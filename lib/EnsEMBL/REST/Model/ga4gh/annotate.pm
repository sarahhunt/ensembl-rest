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

## fetch by ga var


package EnsEMBL::REST::Model::ga4gh::annotate;

use Moose;
extends 'Catalyst::Model';
use Time::HiRes qw(gettimeofday);
use Bio::EnsEMBL::Variation::Utils::VEP qw( get_version_data get_all_consequences vf_to_consequences parse_line read_cache_info );
use Bio::DB::Fasta;
use Data::Dumper;
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;
with 'Catalyst::Component::InstancePerContext';

has 'context' => (is => 'ro');

has 'fasta_db' => (
  isa => 'Bio::DB::Fasta',
  is => 'ro',
  lazy => 1,
  builder => '_find_fasta_cache',
);

has 'fasta' => (
  isa =>'Str',
  is =>'ro'
);

sub build_per_context_instance {
  my ($self, $c, @args) = @_;
  return $self->new({ context => $c, %$self, @args });
}



=head2 annotateVariants

 handles POST /annotate/variant
 requires array of variants in GA4GH format

=cut
sub annotateVariants{

  my ($self, $data ) = @_; 
  print "Starting at ". localtime() ."\n";
  my @gavar_an;

  ## create VEP config 
  my $config = $self->create_config();

## FIX get all version info 
#  my $vep_version = get_version_data($config);

  my $creationtime =  int (gettimeofday * 1000);


  foreach my $gavar (@{$data->{variants}}){

    my @ga_annotation;
    my $vf = $self->convert_to_vf( $config, $gavar);

    no warnings 'redefine';
    local *Bio::EnsEMBL::Slice::seq = $self->_new_slice_seq();

    my $varfeat =  $vf->[0];
    my $consequences = get_all_consequences( $config, $vf);

    my $var_annotation;
    ## copy id from input variant
    $var_annotation->{id}      = $gavar->{id};
    $var_annotation->{created} = $creationtime ; 

    ##save colocated 
    my $existing_var;

    foreach my $ann (@{$consequences}){

      $existing_var =  $ann->{Existing_variation}  
       if defined $ann->{Existing_variation} && $ann->{Existing_variation} =~/^rs|COSMIC|HGMD/;

      my $ga_annotation; 

      $ga_annotation->{id}       = "ph_" .$gavar->{id} . "_" . $ann->{Feature};

      $ga_annotation->{annotationSetId}  = 'Ensembl_80_GRCh38'; ## derive from version
      $ga_annotation->{feature}  = $ann->{Feature};

      push @{$ga_annotation->{effects}}, $ann->{Consequence};

      $ga_annotation->{impact}   = $ann->{IMPACT} if defined $ann->{IMPACT} ;

      $ga_annotation->{feature}  = $ann->{Feature};

      eval{ 
        ## current problem getting HGVSg for intergenic/ intronic
        my $hgvsg = $varfeat->get_all_hgvs_notations('', 'g') ;
        $ga_annotation->{HGVSg}    = $hgvsg->{$ann->{Allele}};
      };
      warn "Problem getting HGVSg : $@\n" if defined $@; 

      $ga_annotation->{HGVSc}    = $ann->{Extra}->{HGVSc} if defined $ann->{Extra}->{HGVSc};
      $ga_annotation->{HGVSp}    = $ann->{Extra}->{HGVSp} if defined $ann->{Extra}->{HGVSp};
      $ga_annotation->{impact}   = $ann->{Extra}->{IMPACT} if defined $ann->{Extra}->{IMPACT};
    
      $ga_annotation->{alternateSequence} =  $ann->{Allele};

      if(defined $ann->{Amino_acids}){
        my @aa   = split /\//, $ann->{Amino_acids};
        my @ppos = split(/\-/, $ann->{Protein_position} ) if $ann->{Protein_position} =~/\-/;
        $ga_annotation->{proteinLocation}->{referenceSequence} = $aa[0];
        $ga_annotation->{proteinLocation}->{alternateSequence} = $aa[1];
        $ga_annotation->{proteinLocation}->{overlapStart}      = $ppos[0] || $ann->{Protein_position};
        $ga_annotation->{proteinLocation}->{overlapEnd}        = $ppos[1] || $ann->{Protein_position};
        $ga_annotation->{proteinLocation}->{overlapStart}      = $ga_annotation->{proteinLocation}->{overlapStart}  - 2;
        $ga_annotation->{proteinLocation}->{overlapEnd}        = $ga_annotation->{proteinLocation}->{overlapEnd}  -1;
      }

      if( defined $ann->{Codons}){
        my @codons = split /\//, $ann->{Codons};
        my @cpos   = split/\-/,$ann->{CDS_position} if $ann->{CDS_position} =~/\-/;
        $ga_annotation->{cdsLocation}->{referenceSequence} = $codons[0];
        $ga_annotation->{cdsLocation}->{alternateSequence} = $codons[1];
        $ga_annotation->{cdsLocation}->{overlapStart}      = $cpos[0]  || $ann->{CDS_position};
        $ga_annotation->{cdsLocation}->{overlapEnd}        = $cpos[1] || $ann->{CDS_position};
        $ga_annotation->{cdsLocation}->{overlapStart}      = $ga_annotation->{cdsLocation}->{overlapStart} -2;
        $ga_annotation->{cdsLocation}->{overlapEnd}        = $ga_annotation->{cdsLocation}->{overlapEnd} -1;
      }

      if(defined $ann->{cDNA_position} && $ann->{cDNA_position}=~/\d+/){

        my @cdnapos = split/\-/,$ann->{cDNA_position} if $ann->{cDNA_position} =~/\-/;

	$ga_annotation->{cDNALocation}->{overlapStart}      =  $cdnapos[0] || $ann->{cDNA_position};
	$ga_annotation->{cDNALocation}->{overlapEnd}        =  $cdnapos[1] || $ann->{cDNA_position};
        $ga_annotation->{cDNALocation}->{overlapStart}      =  $ga_annotation->{cDNALocation}->{overlapStart} - 2;
        $ga_annotation->{cDNALocation}->{overlapEnd}        =  $ga_annotation->{cDNALocation}->{overlapEnd} -1;

      }

      if(defined $ann->{Extra}->{SIFT}){
        my $sift_analysis;
        $sift_analysis->{analysisResult} = (split/\(|\)/, $ann->{Extra}->{SIFT})[0];
        $sift_analysis->{analysisScore}  = (split/\(|\)/, $ann->{Extra}->{SIFT})[1];
        $sift_analysis->{analysis}       = { id          => 'ph_sift_date',
                                             description => 'SIFT',
                                             software    => 'SIFT.5.2.2',
                                             info        => { 'protein database' => 'UniRef90 2014_11'}
                                           };
        push @{$ga_annotation->{analysisResults}} , $sift_analysis;
      }

      if(defined $ann->{Extra}->{PolyPhen}){
        my $polyphen_analysis;
        $polyphen_analysis->{analysisResult} = (split/\(|\)/, $ann->{Extra}->{PolyPhen})[0];
        $polyphen_analysis->{analysisScore}  = (split/\(|\)/, $ann->{Extra}->{PolyPhen})[1];
        $polyphen_analysis->{analysis}       = { id          => 'ph_polyphen_date',
                                                 description => 'Polyphen',
                                                 software    => 'Polyphen.2.2.2_r405',
                                                 info        => { 'protein database' => 'UniRef100 Release 2011_12'}
                                               };

        push @{$ga_annotation->{analysisResults}} , $polyphen_analysis;
      }
  
     push @ga_annotation, $ga_annotation;

    }
    @{$var_annotation->{featureAnnotation}} = @ga_annotation;

    push @{$var_annotation->{coLocatedVariants}}, $existing_var if defined $existing_var;

    push  @gavar_an,  $var_annotation;   
  }

#local $Data::Dumper::Terse = 1;
#local $Data::Dumper::Indent = 1;
# $self->context->log->debug(Dumper {"variantAnnotations" => \@gavar_an});
  print "\nEnding at ". localtime() ."\n";
  return ({ "variantAnnotations"   => \@gavar_an});

}

sub convert_to_vf{

  my ($self, $config, $gavar ) = @_;

  ## convert to VCF line & put through ususal VEP parser
  my $vcf_start = $gavar->{start} + 1;
  my $vcf_alt   = join(",",@{$gavar->{alternateBases}});
  my $vcf_line  = "$gavar->{referenceName}\t$vcf_start\t$gavar->{name}\t$gavar->{referenceBases}\t$vcf_alt";

  no warnings 'redefine';
  local *Bio::EnsEMBL::Slice::seq = $self->_new_slice_seq();
 
  return parse_line($config,$vcf_line) ;
}

=head
## Quicker to handle in bilk 
sub convert_list_to_vf{

  my ($self, $config, $data ) = @_;

  warn "In annotateVariants\n";

  my @vfs;

  foreach my $gavar (@{$data->{variants}}){

    ## convert to VCF line & put through ususal VEP parser
    my $vcf_start = $gavar->{start} + 1;
    my $vcf_alt   = join(",",@{$gavar->{alternateBases}});
    my $vcf_line  = "$gavar->{referenceName}\t$vcf_start\t$gavar->{name}\t$gavar->{referenceBases}\t$vcf_alt";

    push @vfs, @{ parse_line($config,$vcf_line) };
  }
  return  \@vfs;
}
=cut

sub create_config {
  my ($self) = @_;
  
  my $c = $self->context();

  my %vep_params = %{ $c->config->{'Controller::VEP'} };
  read_cache_info(\%vep_params);

  $vep_params{hgvs}   = 1;
  $vep_params{format} = 'vcf';
  undef $vep_params{rest} ;

  return \%vep_params;
}



#### copied from VEP - to be remerged
sub _new_slice_seq {
  # replacement seq method to read from FASTA DB
  my $self = shift;

  my $fasta_db = $self->fasta_db;
  return sub {
    my $self = shift;
    my $seq = $fasta_db->seq( $self->seq_region_name, $self->start => $self->end );
    $seq ||= 'N' x $self->length();
    #reverse_comp( \$seq ) if $self->strand < 0;
    # default to a string of Ns if we couldn't get sequence
    return $seq;
  };
};

sub _find_fasta_cache {
  my $self = shift;
  my $fasta_db = Bio::DB::Fasta->new($self->fasta);
  return $fasta_db;
}



with 'EnsEMBL::REST::Role::Content';

__PACKAGE__->meta->make_immutable;

1;




1;
