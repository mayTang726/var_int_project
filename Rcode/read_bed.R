
# source("Rcode/db_connect_common.R")

#parse .bed file
bed_lines <- readLines('./data/hotspot_region_Oncomine_Focus_DNA_Hotspots_v1.4_hg19.bed', encoding = "UTF-8")[-1]
bed_data <- read.table(text = bed_lines, header = FALSE, sep = "\t")
# change list to dataframe
bed_df <- data.frame(
  chromosome = character(),
  start_position = numeric(),
  end_position = numeric(),
  cosmic_id = character(),
  score = numeric(),
  strand = character(),
  variant_info = character(),
  stringsAsFactors = FALSE  # 确保字符型数据不会转换为因子
)
chromosome <- c(bed_data[1])
start_position <- c(bed_data[2])
end_position <- c(bed_data[3])
cosmic_id <- c(bed_data[4])
score <- c(bed_data[5])
strand <- c(bed_data[6])
variant_info <- c(bed_data[7])

# whole df
bed_df = data.frame(chromosome,start_position,end_position,cosmic_id,score,strand,variant_info)
colnames(bed_df) <- c("chromosome", "start_position", "end_position", "cosmic_id", "score", "strand", "variant_info")

#get chr17 df
# chr17_df <- subset(bed_df, chromosome == 'chr17')


# 使用strsplit()函数将variant_info列的值分割成一个列表
split_info <- strsplit(bed_df$variant_info, ";")

# 初始化用于存储拆分值的向量
REF <- character(length = nrow(bed_df))
OBS <- character(length = nrow(bed_df))
ANCHOR <- character(length = nrow(bed_df))

# 从列表中提取REF，OBS和ANCHOR的值
for (i in 1:length(split_info)) {
  for (j in 1:length(split_info[[i]])) {
    if (grepl("^REF=", split_info[[i]][j])) {
      REF[i] <- gsub("^REF=(.*)$", "\\1", split_info[[i]][j])
    } else if (grepl("^OBS=", split_info[[i]][j])) {
      OBS[i] <- gsub("^OBS=(.*)$", "\\1", split_info[[i]][j])
    } else if (grepl("^ANCHOR=", split_info[[i]][j])) {
      ANCHOR[i] <- gsub("^ANCHOR=(.*)$", "\\1", split_info[[i]][j])
    }
  }
}

# 创建REF，OBS和ANCHOR列，用于合成搜索内容
new_df <- data.frame(REF = REF, OBS = OBS, ANCHOR = ANCHOR)


# 将新的数据框添加到原始数据框中
bed_df <- cbind(bed_df, new_df)
bed_df$variant_info <- paste0(bed_df$ANCHOR, ":", bed_df$OBS)

bed_df$search_column <- apply(bed_df, 1, function(row){
  paste(row['chromosome'],row['start_position'],row['variant_info'], sep = ':')
})

# chromosome variant position - based on hg19
varraint_array <- as.list(bed_df$search_column)

library(jsonlite)
json_data <- toJSON(varraint_array,pretty = TRUE)
# store all chromosome variant position to the .txt file for varsome using
write(json_data,file="data/chr_position_varriant_array_hg19.txt")

#get chr17 df
chr17_df <- subset(bed_df, chromosome == 'chr17')
varraint_array_chr17 <- as.list(chr17_df$search_column)
json_data_chr17 <- toJSON(varraint_array_chr17,pretty = TRUE)

# store chromosome 17 variant position to the .txt file for varsome using
write(json_data_chr17,file="data/chr17_position_varriant_array_hg19.txt")
