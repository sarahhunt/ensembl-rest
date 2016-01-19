# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http=>//www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

BEGIN {
  use FindBin qw/$Bin/;
  use lib "$Bin/lib";
  use RestHelper;
  $ENV{CATALYST_CONFIG} = "$Bin/../ensembl_rest_testing.conf";
  $ENV{ENS_REST_LOG4PERL} = "$Bin/../log4perl_testing.conf";
}

use Test::More;
use Test::Differences;
use Catalyst::Test ();
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Data::Dumper;

Catalyst::Test->import('EnsEMBL::REST');

my $dba = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $base = '/ga4gh/references/search';

## search by reference set id
my $post_data1  = '{ "referenceSetId": "GRCh37.p13", "pageSize": 1  }';

my $expected_post_data =  {
  nextPageToken => 1, 
   references => [       
     { isPrimary => "true",
      sourceURI => 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.11.fa.gz',
      name => "11",
      sourceAccessions => [
         "NC_000011.9"
      ],
      md5checksum =>  "98c59049a2df285c76ffb1c6db8f8b96",
      ncbiTaxonId => 9609,
      length => 135006516,
      sourceDivergence => undef,
      isDerived => 'true',
      id => "98c59049a2df285c76ffb1c6db8f8b96"}
    ]
};


my $json = json_POST($base, $post_data1, 'references by referenceSetId');
eq_or_diff($json, $expected_post_data, "Checking the result from the GA4GH references endpoint by referenceSetId");



## GET

$base =~ s/\/search//;
my $id = '98c59049a2df285c76ffb1c6db8f8b96';
my $json_get = json_GET("$base/$id", 'get references');

my $expected_get_data =  { 
      isPrimary => "true",
      sourceURI => 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.11.fa.gz',
      name => "11",
      sourceAccessions => [
         "NC_000011.9"
      ],
      md5checksum =>  "98c59049a2df285c76ffb1c6db8f8b96",
      ncbiTaxonId => 9609,
      length => 135006516,
      sourceDivergence => undef,
      isDerived => 'true',
      id => "98c59049a2df285c76ffb1c6db8f8b96"
} ; 

eq_or_diff($json_get, $expected_get_data, "Checking the get result from the references endpoint");

my $seq_id = '1d3a93a248d92a729ee764823acbbc6b';
my $query = "$base/$seq_id/bases?start=1080164&end=1080194";
my $json_seq_get = json_GET( $query, 'get references');

my $expected_get_seq_data =  { sequence => 'CTCAAATAAGAGCCACAAACGTGGAAGATA',
                               offset   => 1080164,
                               nextPageToken  => undef
                              };

eq_or_diff($json_seq_get, $expected_get_seq_data, "Checking the get bases result from the references endpoint");





done_testing();

