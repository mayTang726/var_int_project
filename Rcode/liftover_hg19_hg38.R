# liftover page
getwd()
# import request function
source("Rcode/db_connect_common.R") # request file
source("utils/request_body_resolve.R") # request body resolve function
library(jsonlite)
# # import variant file
# varaint_array_hg19 <- readLines('data/chr17_position_varriant_array_hg19.txt')


# entry object
varsome_url <- "https://api.varsome.com/lookup/lifted-over-variant/"

# result_resolve
result_resolve <- function(search_result,db_type,file_name) {
  json_data <- toJSON(search_result,pretty = TRUE)
  name <- gsub("hg19", "hg38", file_name)
  path <- paste('data/bed_parsed_hg38/',name)
  path <- gsub(" ","", path)
  write(json_data,file = path)
}

# request function
search_result <- list()
db_check_fun <- function(db_type, file_name) {
  path <- paste('data/bed_parsed_hg19/',file_name)
  path <- gsub(" ","", path)
  varaint_array <- readLines(path)
  varsome_liftover_obj <- list(
    API = TRUE,
    param_type = "varsome",
    url = varsome_url,
    req_type = "GET",
    token = "Token W0beUnL?kAjoieM4ueklHW%p5AzMVaG@ju@8dD@y",
    body = liftover_body_fun(varaint_array)
  )
  search_result <- connect_function(varsome_liftover_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result,db_type, file_name)
}

# 循环处理当前文件夹下面的所有文件
file_list <-as.list(list.files(path = 'data/bed_parsed_hg19'))
for (i in 1:length(file_list)) {
  print(file_list[i])
  db_check_fun("varsome_liftover",file_list[i])
}



