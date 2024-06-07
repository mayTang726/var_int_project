# this version just based on grch37 search result
# grch37 request link : https://grch37.rest.ensembl.org/vep/human/hgvs
# grch38 request link: https://rest.ensembl.org/vep/human/hgvs

getwd()
# source("utils/request_body_resolve.R")
source("Rcode/resolve_varsome_response.R")
source("Rcode/db_connect_common.R")

## connect cosmic
db_type <- "SIFT"  # 到时候要修改为动态的
sift_url <- "https://grch37.rest.ensembl.org/vep/human/hgvs"

#request body prepare
## combine ensembl_name and hgvs as a string and put it to the list
sift_request_body_list <- hg19_df %>%
  filter(!is.na(ensembl_name) & !is.na(ensembl_hgvs)) %>%
  mutate(combined = paste(ensembl_name, ensembl_hgvs, sep = ":")) %>%
  pull(combined)

#dbname, host, port, user, password, API, parameter_num
sift_obj <- list(
  API = TRUE,  
  param_type = "SIFT",
  url = sift_url,
  req_type = "POST",
  body = sift_request_body_list
)
# calculate score
response_result <- list()
result_resolve <- function(search_result,db_type) {
  # print('进入回调中')
  # print(search_result$response_data)
  response_result <- search_result$response_data
  response_result <- as.data.frame(response_result)
  print(response_result)
}


db_check_fun <- function(db_type) {
  search_result <- connect_function(sift_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result,db_type)
}

db_check_fun(db_type)

