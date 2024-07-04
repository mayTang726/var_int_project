getwd()
source("Rcode/resolve_varsome_response.R")
source("Rcode/db_connect_common.R")

## connect cosmic
db_type <- "oncokb"  # 到时候要修改为动态的

#request body prepare
oncokb_queries_list <- list(
  referenceGenome = "",  
  hugoSymbol = "BRAF",
  entrezGeneId = "",
  alteration = "V600E",
  consequence = "",
  proteinStart = "",
  proteinEnd = "",
  tumorType = ""
)
oncokb_queries_params <- paste0("referenceGenome=", oncokb_queries_list$referenceGenome,"&",   
                              "hugoSymbol=", oncokb_queries_list$hugoSymbol, "&", 
                              "entrezGeneId=", oncokb_queries_list$entrezGeneId, "&",  
                              "alteration=", oncokb_queries_list$alteration,"&",
                              "consequence=", oncokb_queries_list$consequence,"&",
                              "proteinStart=", oncokb_queries_list$proteinStart,"&",  
                              "proteinEnd=", oncokb_queries_list$proteinEnd,"&",
                              "tumorType=", oncokb_queries_list$tumorType)
print(oncokb_queries_params)
oncokb_url <- "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange" 

oncokb_obj <- list(
  API = TRUE,  
  param_type = "oncokb",
  url = oncokb_url,
  req_type = "POST",
  body = oncokb_queries_params,
  token = "Bearer 6481744b-096a-493e-a8ef-5617b848c1cf"   ## token expired in 344 days
)
# calculate score
response_result_oncokb <- data.frame()
result_resolve <- function(search_result,db_type) {
  # print('进入回调中')
  # print(search_result$response_data)
  # output: Oncogenic / Likely Oncogenic
  # resolve response format to df
  response_result_oncokb <<- search_result$response_data
  
}


db_check_fun <- function(db_type) {
  search_result <- connect_function(oncokb_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result,db_type)
}

db_check_fun(db_type)

