# liftover page
getwd()
# import request files
source("Rcode/db_connect_common.R") # request file
source("utils/request_body_resolve.R") # request body resolve function

# import variant position file
varaint_array_hg19 <- readLines('data/chr17_position_varriant_array_hg19.txt')

varsome_19_obj <- list(
  API = TRUE,
  param_type = "varsome",
  url = "https://api.varsome.com/lookup/batch/hg19?add-source-databases=all&add-ACMG-annotation=1&add-AMP-annotation=1",
  req_type = "POST",
  token = "Token W0beUnL?kAjoieM4ueklHW%p5AzMVaG@ju@8dD@y",
  body = request_body(varaint_array_hg19)
)

# result_resolve
result_resolve <- function(search_result,db_type) {
  write(search_result$response_data,'data/variant_response_hg19.json')
}

# request function
search_result <- list()
db_check_fun <- function(db_type) {
  search_result <- connect_function(varsome_19_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result,db_type)
}

db_check_fun("varsome_19")

