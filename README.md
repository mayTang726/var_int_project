# var_int_project
this is a internship project, using R language and Rstudio as development tool

# project structure

1. read_bed.R
  - read .bed file
  - create the query parameter for request link : " chr:position:ref_allele:alt_allele"
  - put all of the query into different .txt file based on different chromosome for varsome   request using 
  
2. liftover_hg19_hg38.R --- limited request usage
  - change all of the variant position in hg19 to hg38 (this code based on 4-part genomic variant specification, it also accept HGVS Protein-level variant, HGVS DNA-level variant, rs_id, variant_id)
  - set changed variant position information in the local file -> "chr17_position_varriant_array_hg38.txt"
  
3. hg19_request.R( based on chromosome17)
  
4. hg38_request.R( based on chromosome17)
 
5. resolve_varsome_response.R
  - unset the data from varsome based on hg19 and hg38
  - save them in the dataframe
  
6. resolve_varsome_response_all.R
  - connect all information from varsome
  - connect cosmic database get sample type for each variant
  - connect sift / mutationtaster database get prediction
  - write data to local file (one chromosome one csv file)
  
### separate token to config file



