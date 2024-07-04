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
response_result_mutation <- data.frame()
result_resolve <- function(search_result,db_type) {
  # print('进入回调中')
  # print(search_result$response_data)
  # output: polyphen_prediction(有害), Polymorphism(可能是无害的)
  # resolve response format to df
  response_result_mutation <<- as.data.frame(search_result$response_data)
  response_result_mutation_new <- response_result_mutation %>%
    pivot_longer(
      cols = everything(),
      names_to = c(".value", "set"),
      names_sep = "\\."
    ) %>%
    select(-set)
  
  # add chromosome variant column
  response_result_mutation_new <- response_result_mutation_new %>%
    mutate(chr_variant = paste0('chr',chr,':',pos,':',ref,':', alt))
  
  # left useful column for matching gh19_df
  response_result_mutation_new <- response_result_mutation_new %>%
    select(chr_variant, prediction,transcript_stable)
  # delete rows which are full NA
  response_result_mutation_new <- response_result_mutation_new %>%
    filter(rowSums(is.na(.)) < ncol(.))
  # match with hg19_df 
  hg19_df <- left_join(hg19_df, response_result_mutation_new,
                         by = c("original_variant" = "chr_variant", "ensembl_name" = "transcript_stable"),
                         relationship = "many-to-many"
  )
  hg19_df <- hg19_df %>% rename(mutation_taster_prediction = prediction)
}


db_check_fun <- function(db_type) {
  search_result <- connect_function(mutation_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result,db_type)
}

db_check_fun(db_type)

