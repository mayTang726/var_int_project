# this version just based on grch37 search result
# grch37 request link : https://grch37.rest.ensembl.org/vep/human/hgvs
# grch38 request link: https://rest.ensembl.org/vep/human/hgvs

#  sift prediction based on amino acids variation
getwd()
# source("Rcode/resolve_varsome_response.R")

#import request function 
source("Rcode/db_connect_common.R") 

# connect sift, get sift prediction
sift_url <- "https://grch37.rest.ensembl.org/vep/human/hgvs"
# combine ensembl_name and hgvs as a string and put it to the list

sift_request_body_fun <- function(df_name) {
  sift_request_body_list <- df_name %>%
    filter(!is.na(ensembl_name) & !is.na(ensembl_hgvs)) %>%
    mutate(combined = paste(ensembl_name, ensembl_hgvs, sep = ":")) %>%
    distinct(combined) %>%
    pull(combined)
  return(sift_request_body_list)
}

sift_obj_fun <- function(df_name){
  sift_obj <- list(
    API = TRUE,
    param_type = "SIFT",
    url = sift_url,
    req_type = "POST",
    body = sift_request_body_fun(df_name)
  )
  return(sift_obj)
}
response_result <- data.frame()
result_resolve <- function(search_result, db_type, df_name) {
  print("进入回调中")
  # parse response to dataframe
  response_result <- search_result$response_data[[1]]
  response_result <- fromJSON(response_result)
  response_result <<- as.data.frame(response_result)
}

# request function
db_check_fun <- function(db_type,df_name) {
  search_result <- connect_function(sift_obj_fun(df_name)) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result, db_type, df_name)
}
sift_pre_function <- function(df_name){
  # print(df_name)
  # abstract transcription column for matching varsome response
  response_result_new <- data.frame()
  response_result_new <- response_result %>%
    select(id, transcript_consequences)
  response_result_new <- response_result_new %>%
    unnest(transcript_consequences) %>% filter(!is.na(amino_acids)) 
  response_result_new <- unique(response_result_new)
  # because the presiction based on amino acids variation, just need to match transcriptid 和 amino acids
  # connect transcriptid, protein_start, amino_acids value to hgvsp str for matching hg19_df
  create_hgvsp <- function(amino_acids, start_pos) {
    # 1. 先判断是否包含 / -> 不包含的话，直接取第一个字母
    # 2. 再判断是否第一个和第三个字母是否相同，相同的话，直接取第一个字母，不相同的话，第一位和第三位都取
    print(amino_acids)
    if(is.na(amino_acids)){
      paste0("NA")
    }else{
      if(nchar(amino_acids) >= 3){
        if(nchar(amino_acids) == 3){
          if(substr(amino_acids, 1, 1) == substr(amino_acids, 3, 3)){
            paste0("p.",substr(amino_acids, 1, 1), start_pos, "=")
          }else{
            paste0("p.",substr(amino_acids, 1, 1), start_pos, substr(amino_acids, 3, 3))
          }
        }else{
          paste0("p.",substr(amino_acids, 1, 1), start_pos, substr(amino_acids, 3, nchar(amino_acids)))
        }
      }else{
        paste0("p.",amino_acids, start_pos, "=")
      }
    }
  }
  # check if protein_start equal protein_end，create hgvsp column
  # mutate() from dplyr package, function: add new column / modify existed column
  response_result_new <- response_result_new %>%
    mutate(hgvsp = mapply(create_hgvsp, amino_acids, protein_start))
  #### match varsome resposne based on transcriptid, amino acids variatoin
  # abstract sift information and transcript_id, hgvsp(one letter)
  response_result_new_select <- response_result_new %>%
    select(transcript_id, hgvsp, sift_prediction)
  df_name <- df_name %>%
    mutate(ensembl_name_no_version = sub("\\..*", "", ensembl_name))
  print(typeof(df_name$ensembl_name_no_version))
  df_name <- merge(df_name, response_result_new_select,
                   by.x = c("ensembl_name_no_version","ensembl_hgvs_p1"),
                   by.y = c("transcript_id","hgvsp"),
                   all.x = TRUE)
  df_name$ensembl_name_no_version <- NULL
  df_name <- unique(df_name)
  print(colnames(df_name))
  return(df_name)
  
}


# source("Rcode/db_connect_common.R")
# 
# sift_url <- "https://grch37.rest.ensembl.org/vep/human/hgvs"
# 
# # request body prepare
# ## combine ensembl_name and hgvs as a string and put it to the list
# sift_request_body_list <- hg19_df %>%
#   filter(!is.na(ensembl_name) & !is.na(ensembl_hgvs)) %>%
#   mutate(combined = paste(ensembl_name, ensembl_hgvs, sep = ":")) %>%
#   distinct(combined) %>%
#   pull(combined)
# 
# # dbname, host, port, user, password, API, parameter_num
# sift_obj <- list(
#   API = TRUE,
#   param_type = "SIFT",
#   url = sift_url,
#   req_type = "POST",
#   body = sift_request_body_list
# )
# # calculate score
# response_result <- data.frame()
# result_resolve <- function(search_result, db_type) {
#   print("进入回调中")
#   # parse response to dataframe
#   response_result <- search_result$response_data[[1]]
#   response_result <- fromJSON(response_result)
#   response_result <- as.data.frame(response_result)
#   # abstract transcription column for matching varsome response
#   response_result_new <- data.frame()
#   response_result_new <- response_result %>%
#     select(id, transcript_consequences)
#   response_result_new <- response_result_new %>%
#     unnest(transcript_consequences)
#   # because the presiction based on amino acids variation, just need to match transcriptid 和 amino acids
#   # connect transcriptid, protein_start, amino_acids value to hgvsp str for matching hg19_df
#   create_hgvsp <- function(amino_acids, start_pos) {
#     codon_split <- strsplit(amino_acids, "/")[[1]]
#     paste0(codon_split[1], start_pos, codon_split[2])
#   }
#   # check if protein_start equal protein_end，create hgvsp column
#   # mutate() from dplyr package, function: add new column / modify existed column
#   
#   response_result_new <- response_result_new %>%
#     mutate(hgvsp = ifelse(protein_start == protein_end,
#                           create_hgvsp(amino_acids, protein_start),
#                           NA
#     ))
#   
#   #### match varsome resposne based on transcriptid, amino acids variatoin
#   # abstract sift information and transcript_id, hgvsp(one letter)
#   response_result_new_select <- response_result_new %>%
#     select(transcript_id, hgvsp, sift_prediction) 
#   hg19_df <- left_join(hg19_df, response_result_new_select,
#                        by = c("ensembl_name" = "transcript_id", "ensembl_hgvs_p1" = "hgvsp"),
#                        relationship = "many-to-many"
#   )
# }
# 
# # request function
# db_check_fun <- function(db_type) {
#   search_result <- connect_function(sift_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
#   result_resolve(search_result, db_type)
# }
# 
# 
# db_check_fun('SIFT')
