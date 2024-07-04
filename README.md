# var_int_project
this is a internship project, using R language and Rstudio as development tool

# project structure

1. read_bed.R
  - read .bed file
  - create the query parameter for request link : " chr:position:ref_allele:alt_allele"
  - put all of the query into a file and set it to the local for request using 
  
2. liftover_hg19_hg38.R
  - change all of the variant position in hg19 to hg38 (this code based on 4-part genomic variant specification, it also accept HGVS Protein-level variant, HGVS DNA-level variant, rs_id, variant_id)
  - set changed variant position information in the local file -> "chr17_position_varriant_array_hg38.txt"
  
3. hg19_request.R
  - get all of the variant information in the different database based on hg19
  - write the information to local file
  
4. hg38_request.R
  - get all of the variant information in the different database baed on hg38
  - write the information to local file
 
5. resolve_varsome_response.R
  - unset the data from varsome based on hg19 and hg38
  - save them in the dataframe

6. scoring_system.R - have not complete
  - connect different database to get the score
  - calculate the score to a total score
  - give a final result for each variant

### separate token to config file