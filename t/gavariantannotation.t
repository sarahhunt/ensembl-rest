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


my $base = '/ga4gh/variantannotations/search';

my $post_data1 = '{"pageSize": 1, "annotationSetIds": [1], "features":[ {"id": "ENST00000381657","featureType": {"source":"SO","name":"transcript","id":"SO:0000673"}} ] }';

my $post_data2 = '{"pageSize": 1, "annotationSetIds": [1], "features":[ {"id": "ENST00000381657","featureType": {"source":"SO","name":"transcript","id":"SO:0000673"}} ],"effects" : [ {"source":"SO","name":"missense_variant","id":"SO:0001583"}]  }';
 
my $post_data3 = '{ "referenceName": "X","start": 215790 ,"end": 215978, "annotationSetIds": [1] ,"pageSize": 1}';

my $post_data4 = '{ "referenceName": "X","start": 215000 ,"end": 217000 ,"pageSize": 1, "annotationSetIds": [1],  "pageToken": 215811  }';

my $post_data5 = '{ "referenceName": "X","start": 215000 ,"end": 217000 ,"pageSize": 1, "annotationSetIds": [1],  "effects" : [ {"source":"SO","name":"missense_variant","id":"SO:0001583"}]  }';



my $expected_data1 = {                                             
  nextPageToken => 'ENST00000381657_208197',                      
  variant_annotation => [                                        
    {                                                           
      transcriptEffects => [                                        
        {                                                           
          HGVSc => 'ENST00000381657.1:c.-205+4695_-205+4696insGCT', 
          HGVSg => '6:g.1084924_1084925insGCT',                    
          HGVSp => undef,                                         
          alternateBase => 'GCT',                                
          effects => [                                          
            {   
              id => 'SO:0001627',                              
              name => 'intron_variant',                       
              ontologySource => 'Sequence Ontology'          
            },                                              
            {                                              
              id => 'SO:0001906',                         
              name => 'feature_truncation',            
              ontologySource => 'Sequence Ontology'     
            }                                          
          ],                                          
          impact => 'MODIFIER'                       
        }                                           
    ],                                           
    variantId => 'tmp__',                                                                                                           
    annotationSetId => 'Ensembl_79',         
    created => 'FIXME_release_date',                                            
   }
 ]
}; 
  
            

my $expected_data2 = { nextPageToken => 'ENST00000381657_208197',
  variant_annotation => [                          
    {                                                       
      transcriptEffects => [                   
        {                                         
          HGVSc => 'ENST00000381657.2:c.943G>A', 
          HGVSg => 'X:g.215973G>A',             
          HGVSp => 'ENSP00000371073.2:p.Ala315Thr',  
          alternateBase => 'A',                     
          analysisResults => [                     
            {                                     
              analysis => {                      
                description => 'SIFT',          
                id => 'placeholder',           
                software => 'SIFT.5.2.2'      
              },                             
              analysisResult => 'tolerated',
              analysisScore => '0.51'      
            }                             
          ],                             
          cDNALocation => {             
            overlapEnd => 1456,        
            overlapStart => 1455      
          },                         
          cdsLocation => {          
            alternateSequence => 'ACG',            
            overlapEnd => 942,                    
            overlapStart => 941,                 
            referenceSequence => 'GCG'          
          },                                   
          effects => [                       
            {                              
              id => 'SO:0001583',           
              name => 'missense_variant',           
              ontologySource => 'Sequence Ontology'
            }                                    
          ],                                
          impact => 'MODERATE',                 
          proteinLocation => {                 
            alternateSequence => 'T',         
            overlapEnd => 314,               
            overlapStart => 313,           
            referenceSequence => 'A'      
         }                              
        }                               
      ],                               
      variantId => 'COSM1119154',                                            
      annotationSetId => 'Ensembl_79',       
      created => 'FIXME_release_date',  
    }
  ]
};


my $expected_data3 = { 
  nextPageToken => 215811,                             
  variant_annotation => [                              
    {                                                  
      annotationSetId => 'Ensembl_79',                 
      created => 'FIXME_release_date',                 
      transcriptEffects => [                           
        {                                              
          HGVSc => 'ENST00000381657.2:c.765G>A',       
          HGVSg => undef,                              
          HGVSp => 'ENST00000381657.2:c.765G>A(p.=)',  
          alternateBase => 'A',                        
          cDNALocation => {                            
            overlapEnd => 1278,                        
            overlapStart => 1277                       
          },                                           
          cdsLocation => {                             
            alternateSequence => 'ACA',                
            overlapEnd => 764,                         
            overlapStart => 763,                       
            referenceSequence => 'ACG'                 
          },                                           
          effects => [                                 
            {                                          
              id => 'SO:0001819',                      
              name => 'synonymous_variant',            
              ontologySource => 'Sequence Ontology'    
            }                                          
          ],                                           
          impact => 'LOW',                             
          proteinLocation => {                         
            alternateSequence => 'T',                  
            overlapEnd => 254,                         
            overlapStart => 253,                       
            referenceSequence => 'T'                   
          }                                            
        }                                              
      ],                                               
      variantId => 'rs370879507'                       
    }
  ]
};

my $expected_data4 = {
  nextPageToken => 215818,                           
  variant_annotation => [                            
    {                                                
      annotationSetId => 'Ensembl_79',               
      created => 'FIXME_release_date',               
      transcriptEffects => [                         
        {                                            
          HGVSc => 'ENST00000381657.2:c.781G>A',     
          HGVSg => undef,                            
          HGVSp => 'ENSP00000371073.2:p.Val261Ile',  
          alternateBase => 'A',                      
          analysisResults => [                       
            {                                        
              analysis => {                          
                description => 'SIFT',               
                id => 'placeholder',                 
                software => 'SIFT.5.2.2'             
              },                                     
              analysisResult => 'tolerated',         
              analysisScore => '1'                   
            }                                        
          ],                                         
          cDNALocation => {                          
            overlapEnd => 1294,                      
            overlapStart => 1293                     
          },                                         
          cdsLocation => {                           
            alternateSequence => 'ATT',              
            overlapEnd => 780,                       
            overlapStart => 779,                     
            referenceSequence => 'GTT'               
          },                                         
          effects => [                               
            {                                        
              id => 'SO:0001583',                    
              name => 'missense_variant',            
              ontologySource => 'Sequence Ontology'  
            }                                        
          ],                                         
          impact => 'MODERATE',                      
          proteinLocation => {                       
            alternateSequence => 'I',                
            overlapEnd => 260,                       
            overlapStart => 259,                     
            referenceSequence => 'V'                 
          }                                          
        }                                            
      ],                                             
      variantId => 'rs368117410'                     
    }                                                
  ]
};


my $json1 = json_POST( $base, $post_data1, 'variantannotations by transcript' );
eq_or_diff($json1, $expected_data1, "Checking the result from the variantannotation endpoint - by transcript");

my $json2 = json_POST($base, $post_data2, 'variantannotations by transcript & effect');
eq_or_diff($json2, $expected_data2, "Checking the result from the variant annotation endpoint - by transcript & effect");

my $json3 = json_POST($base, $post_data3, 'variantannotations by region & effect');
eq_or_diff($json3, $expected_data3, "Checking the result from the variant annotation endpoint - by region");

my $json4 = json_POST($base, $post_data4, 'variantannotations by region');
eq_or_diff($json4, $expected_data4, "Checking the result from the variant annotation endpoint - by region & token");

## re-using result above
my $json5 = json_POST($base, $post_data5, 'variantannotations by region');
eq_or_diff($json5, $expected_data4, "Checking the result from the variant annotation endpoint - by region & effect");


my $bad_post1 = '{ "pageSize": 2, "annotationSetIds": [1], "features":[ {"id": "ENST00000381657","featureType": {"source":"SO","name":"fish","id":"SO:0000000"}}]}';

action_bad_post($base, $bad_post1, qr/No annotations available for this feature type/, 'Throw bad feature request at endpoint' );

my $bad_post2 = '{ "pageSize": 1, "annotationSetIds": [1], "features":[ {"id": "ENST00000000000","featureType": {"source":"SO","name":"transcript","id":"SO:0000673"}}]}';

action_bad_post($base, $bad_post2, qr/feature ENST00000000000 not found/, 'Throw bad transcript request at endpoint' );




done_testing();


