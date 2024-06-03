# this is the process of request the VarSome API, from request hg19, listover hg19 -> hg38, request hg38

getwd()
source("Rcode/db_connect_common.R") # request file
varaint_array_hg19 <- readLines('data/chr17_position_varriant_array_hg19.txt')
varaint_array_hg38 <- readLines('data/chr17_position_varriant_array_hg38.txt')

# function for varsome variant post request body
request_body <- function (array) {
  # change list to json charactor for post request use
  variants_char <- character()
  for (i in 2:length(array)-1) {
    if( i != 1){
      character_ <- array[[i]]
      character_ <- gsub("\\["," ",character_)
      character_ <- gsub("\\]"," ",character_)
      character_ <- gsub(","," ",character_)
      variants_char <- append(variants_char, character_)
    }
  }
  print(variants_char)
  variants_char <- paste(variants_char, collapse = ',')
  variants_char <- paste("{'variants':","[",variants_char,"]}")
  variants_char <- gsub(" ","",variants_char)
  variants_char <- gsub("'",'"',variants_char)
  return (variants_char)
}

# function for liftover get request body
liftover_body_fun <- function(array){
  liftover_char <- list()
  for (i in 2:47) {
    item_ <- array[[i]]
    item_ <- gsub("\\["," ",item_)
    item_ <- gsub("\\]"," ",item_)
    item_ <- gsub(","," ",item_)
    item_ <- gsub("\"","",item_)
    item_ <- gsub(" ","",item_)
    liftover_char <- append(liftover_char, item_)
  }
  return(liftover_char)
}

# entry object
db_type <- "varsome_liftover"  # 到时候要修改为动态的
varsome_url <- ""
if(db_type == "varsome_liftover"){
  varsome_url <- "https://api.varsome.com/lookup/lifted-over-variant/"
}else if(db_type == 'varsome_19') {
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

varsome_liftover_obj <- list(
  API = TRUE,
  param_type = "varsome",
  url = varsome_url,
  req_type = "GET",
  token = "Token W0beUnL?kAjoieM4ueklHW%p5AzMVaG@ju@8dD@y",
  body = variant_list_liftover
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
    write(search_result$response_data,'data/variant_response_hg19.json')
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
    varsome_19_obj$body <- request_body(varaint_array_hg19)
    param <- varsome_19_obj
  }else if(db_type == 'varsome_liftover') {
    varsome_liftover_obj$body <- liftover_body_fun(varaint_array_hg19)
    param <- varsome_liftover_obj
  }else if(db_type == 'varsome_38') {
    varsome_38_obj$body <- request_body(varaint_array_hg38)
    param <- varsome_38_obj
  }
  
  search_result <- connect_function(param) # list[[1]]: param_type, list[[2]]: list ---> searched information
  calculate_score(search_result,db_type)
}

db_check_fun(db_type)

