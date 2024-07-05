getwd()

library(jsonlite)
library(tidyr)
library(dplyr)
library(purrr)

# 1. read json file
# the direction path of result json files
folder_path <- "data/varsome_result_all"
file_list <- list.files(path = folder_path, pattern = "\\.json$", full.names = TRUE)
data_list <- list()
# put all json file to list for changing to df
for (file in file_list) {
  json_data <- fromJSON(file)
  data_df <- as.data.frame(json_data)
  data_list <- append(data_list, list(data_df))
}
# find all columns
all_colnames <- unique(unlist(lapply(data_list, names)))
print(all_colnames)
# fill in missing columns
data_list <- lapply(data_list, function(df) {
  missing_cols <- setdiff(all_colnames, names(df))
  for (col in missing_cols) {
    df[[col]] <- NA
  }
  return(df)
})
# combine all data to df
varsome_all_result_df <- bind_rows(data_list)








