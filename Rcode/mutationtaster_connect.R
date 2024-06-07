getwd()
source("Rcode/resolve_varsome_response.R")
source("Rcode/db_connect_common.R")

## connect cosmic
db_type <- "mutation"  # 到时候要修改为动态的
mutation_url <- "https://www.genecascade.org/MT2021/MT_API102.cgi?variants="

#request body prepare
## the param need format - "<chromosome_number>:<position><pre>%3E<alt>"
## convert original tring to param format
convert_variant <- function(variant) {
  # 移除 "chr" 前缀
  variant <- sub("^chr", "", variant)
  # 保留第一个冒号
  first_char <- regexpr(":", variant)
  first_part <- substr(variant, 1, first_char)
  remaining_part <- substr(variant, first_char + 1, nchar(variant))
  remaining_part <- sub(":(?=[^:]*$)", "%3E", remaining_part, perl = TRUE)
  remaining_part <- gsub(":", "", remaining_part)
  variant <- paste0(first_part, remaining_part)
  return(variant)
}
# delete duplicated elements
mutation_request_body_list <- unique(convert_variant(hg19_df$original_variant))
print(mutation_request_body_list)
#dbname, host, port, user, password, API, parameter_num
mutation_obj <- list(
  API = TRUE,  
  param_type = "mutation",
  url = mutation_url,
  req_type = "GET",
  body = mutation_request_body_list
)
# calculate score
response_result <- list()
result_resolve <- function(search_result,db_type) {
  # print('进入回调中')
  # print(search_result$response_data)
  response_result <- search_result$response_data
  print(response_result)
  # parse the response 
    # prediction_list <- list()
    # obj <- list()
    # for (i in seq_along(search_result$response_data)) {
    #   obj["transcript_stable"] <- data_list[[i]]$transcript_stable
    #   obj["prediction"] <- data_list[[i]]$prediction
    #   prediction_list <- c(prediction_list, list(obj))
    # }
}


db_check_fun <- function(db_type) {
  search_result <- connect_function(mutation_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result,db_type)
}

db_check_fun(db_type)

