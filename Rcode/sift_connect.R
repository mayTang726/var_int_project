# this version just based on grch37 search result
# grch37 request link : https://grch37.rest.ensembl.org/vep/human/hgvs
# grch38 request link: https://rest.ensembl.org/vep/human/hgvs

#  sift prediction based on amino acids variation

getwd()
# source("utils/request_body_resolve.R")
source("Rcode/resolve_varsome_response.R")
source("Rcode/db_connect_common.R")

## connect cosmic
db_type <- "SIFT" # 到时候要修改为动态的
# API address ： https://rest.ensembl.org/documentation/info/vep_hgvs_post
sift_url <- "https://grch37.rest.ensembl.org/vep/human/hgvs"

# request body prepare
## combine ensembl_name and hgvs as a string and put it to the list
sift_request_body_list <- hg19_df %>%
  filter(!is.na(ensembl_name) & !is.na(ensembl_hgvs)) %>%
  mutate(combined = paste(ensembl_name, ensembl_hgvs, sep = ":")) %>%
  distinct(combined) %>%
  pull(combined)
print(sift_request_body_list)

# dbname, host, port, user, password, API, parameter_num
sift_obj <- list(
  API = TRUE,
  param_type = "SIFT",
  url = sift_url,
  req_type = "POST",
  body = sift_request_body_list
)
# calculate score
response_result <- data.frame()
result_resolve <- function(search_result, db_type) {
  print("进入回调中")
  # parse response to dataframe
  response_result <- search_result$response_data[[1]]
  response_result <- fromJSON(response_result)
  response_result <<- as.data.frame(response_result)
  # abstract transcription column for matching varsome response
  response_result_new <- data.frame()
  response_result_new <- response_result %>%
    select(id, transcript_consequences)
  response_result_new <- response_result_new %>%
    unnest(transcript_consequences)
  # because the presiction based on amino acids variation, so just need to match transcriptid 和 amino acids
  # connect transcriptid, protein_start, amino_acids value to hgvsp str for matching hg19_df
  create_hgvsp <- function(amino_acids, start_pos) {
    codon_split <- strsplit(amino_acids, "/")[[1]]
    paste0(codon_split[1], start_pos, codon_split[2])
  }
  # check if protein_start equal protein_end，create hgvsp column
  # mutate() from dplyr package, function: add new column / modify existed column

  response_result_new <- response_result_new %>%
    mutate(hgvsp = ifelse(protein_start == protein_end,
      create_hgvsp(amino_acids, protein_start),
      NA
    ))

  #### match varsome resposne based on transcriptid, amino acids variatoin
  # abstract sift information and transcript_id, hgvsp(one letter)
  response_result_new_select <- response_result_new %>%
    select(transcript_id, hgvsp, sift_prediction) 
  # %>% rename(
  #       score_original = sift_score,
  #       pre_original = sift_prediction
  #     ) %>%
  #     distinct()
  
  # left_join combine matched rows
  # hg19_df <- hg19_df %>%
  #   mutate(
  #     sift_score_original = NA,
  #     sift_pre_original = NA
  #   )
  
  hg19_df <- left_join(hg19_df, response_result_new_select,
    by = c("ensembl_name" = "transcript_id", "ensembl_hgvs_p1" = "hgvsp"),
    relationship = "many-to-many"
  )
  # hg19_df <- merged_df %>%
  #   mutate(
  #     sift_score_original = if_else(is.na(sift_score_original), score_original, sift_score_original),
  #     sift_pre_original = if_else(is.na(sift_pre_original), pre_original, sift_pre_original)
  #   ) %>%
  #   select(-pre_original, -score_original)
  
  #### 有个问题： hg19_df中有多列，但是目前每个相同的值只在hg19_df中只match了一次，其他相同的都是NA，
  #### 这个需要处理一下
}

# request function
db_check_fun <- function(db_type) {
  search_result <- connect_function(sift_obj) # list[[1]]: param_type, list[[2]]: list ---> searched information
  result_resolve(search_result, db_type)
}


db_check_fun(db_type)
