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


Catalyst::Test->import('EnsEMBL::REST');

my $base = '/variants/search';

my $post_data1 = '{ "referenceName": 22,"start": 16050150 ,"end": 16050170 ,"pageSize": 1, "callSetIds": ["NA12878"], "variantSetIds":[65] }';
my $post_data2 = '{ "referenceName": 22,"start": 16132100 ,"end": 16132110 ,"pageSize": 1,  "variantSetIds":[22], "callSetIds": ["NA19060", "NA18990"] ,"variantName": "rs150753069" }';


my $expected_data1 = {  nextPageToken => "16050158_65",
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
               '0',          
               '1'          
             ]             
           }              
         ],              
         end => 16050159,      
         id => '22_16050159', 
         name => '22_16050159',  
         referenceBases => 'C', 
         referenceName => '22',
         start => '16050158', 
         variantSetId => '65'
       }                    
     ]                     
   };    
            

my $expected_data2 = {                                   
  nextPageToken => '16132100_22',   
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
            '1',                    
            '1'                     
          ]                         
        },                          
        {                           
          callSetId => 'NA19060',   
          callSetName => 'NA19060', 
          genotype => [             
            '1',                    
            '1'                     
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
    }                               
  ]                                
};



my $json1 = json_POST( $base, $post_data1, 'variants by callset & varset' );
eq_or_diff($json1, $expected_data1, "Checking the result from the gavariant endpoint - varset & callset");

my $json2 = json_POST($base, $post_data2, 'variants by callset & varset & var name');
eq_or_diff($json2, $expected_data2, "Checking the result from the gavariant endpoint - varset & callset & vaname");


done_testing();


