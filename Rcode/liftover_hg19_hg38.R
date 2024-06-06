# liftover page
getwd()
# import request function
source("Rcode/db_connect_common.R") # request file
source("utils/request_body_resolve.R") # request body resolve function

# import variant position file
varaint_array_hg19 <- readLines('data/chr17_position_varriant_array_hg19.txt')

# entry object
varsome_url <- "https://api.varsome.com/lookup/lifted-over-variant/"

varsome_liftover_obj <- list(
  API = TRUE,
  param_type = "varsome",
  url = varsome_url,
  req_type = "GET",
  token = "Token W0beUnL?kAjoieM4ueklHW%p5AzMVaG@ju@8dD@y",
  body = liftover_body_fun(varaint_array_hg19)
)

# result_resolve
result_resolve <- function(search_result,db_type) {
  library(jsonlite)
  json_data <- toJSON(search_result,pretty = TRUE)
  write(json_data,file="data/chr17_position_varriant_array_hg38.txt")
  
}

# request function
search_result <- list()
db_check_fun <- function(db_type) {
  search_result <- connect_function(varsome_liftover_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result,db_type)
}

db_check_fun("varsome_liftover")

