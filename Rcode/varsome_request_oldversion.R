# this is the process of request the VarSome API, from request hg19, listover hg19 -> hg38, request hg38

getwd()
source("Rcode/db_connect_common.R") # request file
source("utils/request_body_resolve.R") # request body resolve function
# whole chromosome
varaint_array_hg19_whole <- readLines('data/chr_position_varriant_array_hg19.txt')

# chromosome 17
varaint_array_hg19 <- readLines('data/chr17_position_varriant_array_hg19.txt')
varaint_array_hg38 <- readLines('data/chr17_position_varriant_array_hg38.txt')


# entry object
db_type <- "varsome_19"  # 到时候要修改为动态的
varsome_url <- ""
if(db_type == 'varsome_19') {
  varsome_url <- "https://api.varsome.com/lookup/batch/hg19?add-source-databases=all&add-ACMG-annotation=1&add-AMP-annotation=1"
}else if(db_type == 'varsome_38') {
  varsome_url <- "https://api.varsome.com/lookup/batch/hg38?add-source-databases=all&add-ACMG-annotation=1&add-AMP-annotation=1"
}

varsome_19_obj <- list(
  API = TRUE,
  param_type = "varsome",
  url = varsome_url,
  req_type = "POST",
  token = "Token W0beUnL?kAjoieM4ueklHW%p5AzMVaG@ju@8dD@y",
  body = ''
)

varsome_38_obj <- list(
  API = TRUE,
  param_type = "varsome",
  url = varsome_url,
  req_type = "POST",
  token = "Token W0beUnL?kAjoieM4ueklHW%p5AzMVaG@ju@8dD@y",
  body = ''
)

# calculate score
search_result <- list()
calculate_score <- function(search_result,db_type) {
  print('callback')
  if(db_type == 'varsome_liftover'){
    library(jsonlite)
    json_data <- toJSON(search_result,pretty = TRUE)
    write(json_data,file="data/chr17_position_varriant_array_hg38.txt")
  }else if (db_type == 'varsome_19') {
    # write(search_result$response_data,'data/variant_response_hg19.json')
    write(search_result$response_data,'data/whole_variant_response_hg19.json')
  }else {
    print(search_result$response_data)
    write(search_result$response_data,'data/variant_response_hg38.json')
  }
  # search_result <- search_result
}


param <- ""
search_result <- ""

db_check_fun <- function(db_type) {
  if(db_type == 'varsome_19') {
    # varsome_19_obj$body <- request_body(varaint_array_hg19)
    varsome_19_obj$body <- request_body(varaint_array_hg19_whole)
    param <- varsome_19_obj
  }else if(db_type == 'varsome_38') {
    varsome_38_obj$body <- request_body(varaint_array_hg38)
    param <- varsome_38_obj
  }
  
  search_result <- connect_function(param) # list[[1]]: param_type, list[[2]]: list ---> searched information
  calculate_score(search_result,db_type)
}

db_check_fun(db_type)

