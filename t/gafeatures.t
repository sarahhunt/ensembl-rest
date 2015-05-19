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
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Data::Dumper;

my $dba = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('multi');
Catalyst::Test->import('EnsEMBL::REST');


my $base = '/ga4gh/features';

my $post_data1 = '{"pageSize": 2, "featureSetIds": [1],  "features":[{"source":"SO","name":"transcript","id":"SO:0000673"}]}';

my $post_data2 = '{"pageSize": 2, "featureSetIds": [1],  "features":[{"source":"SO","name":"transcript","id":"SO:0000673"}], "pageToken": 1384970}';

my $post_data3 = '{"pageSize": 2, "featureSetIds": [1],  "range":[{"length":10000,"start":{"base":{"position":1080164,"referenceName":"6"}}}] }';


my $expected_data1 = {                                               
  features => [                                 
    {                                           
      attributes => {                           
        biotype => 'protein_coding',            
        created => '1209467861',                
        gene => 'ENSG00000176515',              
        updated => '1209467861',                
        version => '1'                          
      },                                        
      featureSetId => 'placeholder_Ensembl79',  
      featureType => {                          
        id => 'SO:0000673',                     
        name => 'transcript',                   
        source => 'SO'                          
      },                                        
      id => 'ENST00000314040',                  
      path => [                                 
        {                                       
          length => 25017,                      
          start => {                            
            base => {                           
              position => 1080164,              
              referenceName => '6'              
            }                                   
          }                                     
        }                                       
      ]                                         
    },                                          
    {                                           
      attributes => {                           
        biotype => 'snoRNA',                    
        created => '1268996515',                
        gene => 'ENSG00000238438',              
        updated => '1268996515',                
        version => '1'                          
      },                                        
      featureSetId => 'placeholder_Ensembl79',  
      featureType => {                          
        id => 'SO:0000673',                     
        name => 'transcript',                   
        source => 'SO'                          
      },                                        
      id => 'ENST00000459140',                  
      path => [                                 
        {                                       
          length => 102,                        
          start => {                            
            base => {                           
              position => 1186753,              
              referenceName => '6'              
            }                                   
          }                                     
        }                                       
      ]                                         
    }                                           
  ],                                            
  nextPageToken => '1384970'   
}; 
  
            

my $expected_data2 = {
  features => [                                 
    {                                           
      attributes => {                           
        biotype => 'protein_coding',            
        created => '1209467861',                
        gene => 'ENSG00000164379',              
        updated => '1242726437',                
        version => '2'                          
      },                                        
      featureSetId => 'placeholder_Ensembl79',  
      featureType => {                          
        id => 'SO:0000673',                     
        name => 'transcript',                   
        source => 'SO'                          
      },                                        
      id => 'ENST00000296839',                  
      path => [                                 
        {                                       
          length => 2317,                       
          start => {                            
            base => {                           
              position => 1312675,              
              referenceName => '6'              
            }                                   
          }                                     
        }                                       
      ]                                         
    },                                                            
    {                                           
      attributes => {                           
        biotype => 'antisense',                 
        created => '1321005463',                
        gene => 'ENSG00000261730',              
        updated => '1321005463',                
        version => '1'                          
      },                                        
      featureSetId => 'placeholder_Ensembl79',  
      featureType => {                          
        id => 'SO:0000673',                     
        name => 'transcript',                   
        source => 'SO'                          
      },                                        
      id => 'ENST00000568244',                  
      path => [                                 
        {                                       
          length => 1276,                       
          start => {                            
            base => {                           
              position => 1384025,              
              referenceName => '6'              
            }                                   
          }                                     
        }                                       
      ]                                         
    }                                           
  ],                                            
  nextPageToken => '1384972'
};


my $expected_data3 = { 
features => [             
    {                                           
      attributes => {                           
        biotype => 'protein_coding',            
        created => '1209467861',                
        gene => 'ENSG00000176515',              
        updated => '1209467861',                
        version => '1'                          
      },                                        
      featureSetId => 'placeholder_Ensembl79',  
      featureType => {                          
        id => 'SO:0000673',                     
        name => 'transcript',                   
        source => 'SO'                          
      },                                        
      id => 'ENST00000314040',                  
      path => [                                 
        {                                       
          length => 25017,                      
          start => {                            
            base => {                           
              position => 1080164,              
              referenceName => '6'              
            }                                   
          }                                     
        }                                       
      ]                                         
    }                                           
  ],                                            
  nextPageToken => undef       
};




my $postbase = $base ."/search";
my $json1 = json_POST( $postbase, $post_data1, 'sequence annotations ' );
eq_or_diff($json1, $expected_data1, "Checking the result from the sequence annotations endpoint");

my $json2 = json_POST($postbase, $post_data2, 'sequence annotations by page token');
eq_or_diff($json2, $expected_data2, "Checking the result from the sequence annotations endpoint - by page token");

my $json3 = json_POST($postbase, $post_data3, 'sequence annotations by region');
eq_or_diff($json3, $expected_data3, "Checking the result from the sequence annotations endpoint - by region");




## check rubbish handled
my $bad_post1 = '{ "pageSize": 2, "SetIds": [1],"features":[{"source":"SO","name":"transcript","id":"SO:0000673"}]}';

action_bad_post($postbase, $bad_post1, qr/Cannot find/, 'Throw bad request at endpoint' );



### check get

my $id = 'ENST00000568244';
my $json4 = json_GET("$base/$id", 'get transcript');

my $expected_data4 =  {                                           
  features => [  
    {
      attributes => {                           
        biotype => 'antisense',                 
        created => '1321005463',                
        gene => 'ENSG00000261730',              
        updated => '1321005463',                
        version => '1'                          
      },                                        
      featureSetId => 'placeholder_Ensembl79',  
      featureType => {                          
        id => 'SO:0000673',                     
        name => 'transcript',                   
        source => 'SO'                          
      },                                        
      id => 'ENST00000568244',                  
      path => [                                 
        {                                       
          length => 1276,                       
          start => {                            
            base => {                           
              position => 1384025,              
              referenceName => '6'              
            }                                   
          }                                     
        } 
      ]
    }
 ]
};

eq_or_diff($json4, $expected_data4, "Checking the get result from the sequence annotation endpoint");



done_testing();


