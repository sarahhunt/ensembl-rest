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
use Data::Dumper;

Catalyst::Test->import('EnsEMBL::REST');

my $base = '/ga4gh/variants/search';

my $post_data1 = '{ "referenceName": 22,"start": 16050150 ,"end": 16060170 ,"pageSize": 1, "callSetIds": ["NA12878"], "variantSetIds":[65] }';
my $post_data2 = '{ "referenceName": 22,"start": 16132100 ,"end": 16132110 ,"pageSize": 1,  "variantSetIds":[22], "callSetIds": ["NA19060", "NA18990"] ,"variantName": "rs150753069" }';
my $post_data3 = '{ "referenceName": 22,"start": 16132100 ,"end": 16132110 ,"pageSize": 1,  "variantSetIds":[22, 65], "callSetIds": ["NA19060","NA12878"] ,"variantName": "rs150753069" }';
my $post_data4 = '{ "referenceName": 22,"start": 16050150 ,"end": 16060170 ,"pageSize": 1, "callSetIds": ["NA12878"], "variantSetIds":[65], "pageToken":"16050158_65_2" }';

my $expected_data1 = {  nextPageToken => "16050158_65_2",
  variants => [                 
    {                          
     alternateBases => [     
            'T'                  
          ],                    
          calls => [           
            {                         
              callSetId => 'NA12878', 
             callSetName => 'NA12878',
             genotype => [    
               0,          
               1          
             ]             
           }              
         ],              
         end => 16050159,      
         id => '22_16050159', 
         info => { AC => ['1'],
                   AF => ['0.50'], 
                   AN => ['2']}, 
         names => ['22_16050159'],  
         referenceBases => 'C', 
         referenceName => '22',
         start => '16050158', 
         variantSetId => '65',
         created => '1419292800000',
         updated => '1419292800000', 
       }                    
     ]                     
   };    
            

my $expected_data2 = {                                   
  nextPageToken => '16132100_22_1',   
  variants => [                     
    {                               
      alternateBases => [           
        'A'                         
      ],                            
      calls => [                    
        {                           
          callSetId => 'NA18990',   
          callSetName => 'NA18990', 
          genotype => [             
            1,                    
            1                     
          ],
          genotypeLikelihood => ['-0.48', '-0.48', '-0.48'],
          info => { 
            DS => [ 
              '2.000' 
            ]         
          },   
          phaseset => ''                           
        },                          
        {                           
          callSetId => 'NA19060',   
          callSetName => 'NA19060', 
          genotype => [             
            1,                    
            1                     
          ]                         
        }                           
      ],                            
      end => 16132101,              
      id => 'rs150753069',          
      name => 'rs150753069',        
      referenceBases => 'G',        
      referenceName => '22',        
      start => '16132100',          
      variantSetId => '22'          
          ],
          genotypeLikelihood => ['-0.48', '-0.48', '-0.48'],
          info => {             
            DS => [
              '1.950' 
            ]         
          },
          phaseset => ''                            
        }                           
      ],                            
      end => 16132101,              
      id => 'rs150753069',
      info => {AA => ['N'],    
        AC => ['2143'],
        AF => ['0.98'], 
        AFR_AF => ['0.99'], 
        AMR_AF => ['0.99'], 
        AN => ['2184'],     
        ASN_AF => ['0.99'], 
        AVGPOST => ['0.8755'], 
        ERATE => ['0.0050'],   
        EUR_AF => ['0.96'],   
        LDAF => ['0.9238'],  
        RSQ => ['0.2636'], 
        SNPSOURCE => ['LOWCOV'],
        THETA => ['0.0202'],  
        VT => ['SNP']},          
      names => ['rs150753069'],        
      referenceBases => 'G',        
      referenceName => '22',        
      start => '16132100',          
      variantSetId => '22',
      created => '1432745640000',
      updated => '1432745640000',  
    }                               
  ]                                
};

my $expected_data3 = {
  nextPageToken => '16132100_22_1',
  variants => [
    {
      alternateBases => [
        'A'
      ],
      calls => [
        {
          callSetId => 'NA19060',
          callSetName => 'NA19060',
          genotype => [
            1,
            1
          ],
          genotypeLikelihood => ['-0.48', '-0.48', '-0.48'],
          info => {
            DS => [
              '1.950'
            ]
          },
          phaseset => ''
        }
      ],
      end => 16132101,
      id => 'rs150753069',
      info => {AA => ['N'],
        AC => ['2143'],
        AF => ['0.98'],
        AFR_AF => ['0.99'],
        AMR_AF => ['0.99'],
        AN => ['2184'],
        ASN_AF => ['0.99'],
        AVGPOST => ['0.8755'],
        ERATE => ['0.0050'],
        EUR_AF => ['0.96'],
        LDAF => ['0.9238'],
        RSQ => ['0.2636'],
        SNPSOURCE => ['LOWCOV'],
        THETA => ['0.0202'],
        VT => ['SNP']},
      names => ['rs150753069'],
      referenceBases => 'G',
      referenceName => '22',
      start => '16132100',
      variantSetId => '22',
      created => '1432745640000',
      updated => '1432745640000',
    }
  ]
};

my $expected_data4 = {  nextPageToken => "16050251_65_2",
  variants => [
    {
     alternateBases => [
            'T'
          ],
          calls => [
            {
              callSetId => 'NA12878',
             callSetName => 'NA12878',
             genotype => [
               0,
               1
             ]
           }
         ],
         end => 16050252,
         id => '22_16050252',
         info => { AC => ['1'],                  
                   AF => ['0.50'],
                   AN => ['2']},
         names => ['22_16050252'],
         referenceBases => 'A',
         referenceName => '22',
         start => '16050251',
         variantSetId => '65',
         created => '1419292800000',
         updated => '1419292800000',
       }
     ]
   };

my $json1 = json_POST( $base, $post_data1, 'variants by callset & varset' );
eq_or_diff($json1, $expected_data1, "Checking the result from the gavariant endpoint - varset & callset");

my $json2 = json_POST($base, $post_data2, 'variants by callset & varset & var name');
eq_or_diff($json2, $expected_data2, "Checking the result from the gavariant endpoint - varset & callset & vaname");

my $json3 = json_POST($base, $post_data3, 'variants by callset & 2 varsets & var name');
eq_or_diff($json3, $expected_data3, "Checking the result from the gavariant endpoint - 2 datasets (one empty) & callset & vaname");

my $json4 = json_POST($base, $post_data4, 'variants by callset & token');
eq_or_diff($json4, $expected_data4, "Checking the result from the gavariant endpoint - with token");


my $bad_post = q/{ "referenceName": 22,"start": 16050150 ,"end": 16050150 ,"pageSize": 1, "callSetIds": ["NA12878"], "variantSetIds":[65], "pageToken":"16050158_65_2" }/;

action_bad_post($base, $bad_post, qr/must not equal/, 'Throw nasty data at endpoint' );

my $bad_post2 = q/{ "referenceName": 22,"start": 1 ,"end": 10 ,"pageSize": 1, "callSetIds": ["NA12878"], "variantSetIds":[65] }/;

action_bad_post($base, $bad_post2, qr/No variants are available for this region/, 'Throw if no data' );



done_testing();


