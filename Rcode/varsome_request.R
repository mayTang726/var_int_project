# this is the process of request the VarSome API, from request hg19, listover hg19 -> hg38, request hg38

getwd()
source("Rcode/db_connect_common.R")
source("Rcode/read_bed.R")

variant_list_liftover <- varraint_array # for liftover
# set post request body
request_body <- function (array) {
  # change list to json charactor for post request use
  for (i in 1:length(array)) {
    character_ <- paste("'",array[[i]],"'")
    variants_char <- append(variants_char, character_)
  }
  variants_char <- paste(variants_char, collapse = ',')
  variants_char <- paste("{'variants':","[",variants_char,"]}")
  variants_char <- gsub(" ","",variants_char)
  variants_char <- gsub("'",'"',variants_char)
  return (variants_char)
}

db_type <- "varsome_liftover"  # 到时候要修改为动态的
varsome_url <- ""
if(db_type == "varsome_liftover"){
  varsome_url <- "https://api.varsome.com/lookup/lifted-over-variant/"
}else if(db_type == 'varsome_19') {
  varsome_url <- "https://api.varsome.com/lookup/batch/hg19?add-source-databases=all"
}else if(db_type == 'varsome_38') {
  varsome_url <- "https://api.varsome.com/lookup/batch/hg38?add-source-databases=all"
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
  print('进入回调中')
  if (db_type == 'varsome_19') {
    db_type = 'varsome_liftover'
  } else if(db_type == 'varsome_liftover'){
    write(search_result,file="data/chr_list_liftover.txt")
    # db_type = 'varsome_38'
    # request_body(search_result)
    # varsome_38_obj$body <- variants_char
    # param <- varsome_19_obj
    # db_check_fun(db_type)
  }else {
    # hg38 chromosome format one to one corresponding to hg19
    # 当前的search result已经是在liftover之后的内容了，现在调用hg38的内容
    # request_body(search_result)
    db_check_fun(db_type)
  }
  # search_result <- search_result
}


param <- ""
search_result <- ""

db_check_fun <- function(db_type) {
  if(db_type == 'varsome_19') {
    varsome_19_obj$body <- request_body(varraint_array)
    param <- varsome_19_obj
  }else if(db_type == 'varsome_liftover') {
    param <- varsome_liftover_obj
  }else if(db_type == 'varsome_38') {
    varsome_38_obj$body <- request_body(search_result)
    param <- varsome_38_obj
  }
  
  search_result <- connect_function(param) # list[[1]]: param_type, list[[2]]: list ---> searched information
  calculate_score(search_result,db_type)
}

db_check_fun(db_type)

