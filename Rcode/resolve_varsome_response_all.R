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
# unique(varsome_all_result_df$chromosome)
varsome_all_result_df <- varsome_all_result_df[!is.na(varsome_all_result_df$chromosome),]

varsome_all_result_df['regions'] <- NULL 
varsome_all_result_df['gerp'] <- NULL 
varsome_all_result_df['phastcons100way'] <- NULL 
varsome_all_result_df['phylop100way'] <- NULL 
varsome_all_result_df['maxentscan'] <- NULL 
varsome_all_result_df['publications'] <- NULL
varsome_all_result_df['cbio_portal'] <- NULL
varsome_all_result_df['detail'] <- NULL
varsome_all_result_df['weill_cornell_medicine_pmkb'] <- NULL
varsome_all_result_df['icgc_somatic'] <- NULL
varsome_all_result_df['nih_gdc'] <- NULL
varsome_all_result_df['wustl_docm'] <- NULL
varsome_all_result_df['variant_pubmed_automap'] <- NULL
varsome_all_result_df['cadd'] <- NULL
varsome_all_result_df['dann_snvs'] <- NULL
varsome_all_result_df['error'] <- NULL
varsome_all_result_df['cancer_hotspots'] <- NULL
varsome_all_result_df['uniprot_variants'] <- NULL
varsome_all_result_df['lumc_lovd'] <- NULL
varsome_all_result_df['omim'] <- NULL
varsome_all_result_df['gnomad_exomes'] <- NULL
varsome_all_result_df['gnomad_genomes'] <- NULL
varsome_all_result_df['bravo'] <- NULL
varsome_all_result_df['dbnsfp_dbscsnv'] <- NULL
varsome_all_result_df['isb_kaviar3'] <- NULL
varsome_all_result_df['nih_clingen_variants'] <- NULL
varsome_all_result_df['uoi_dvd'] <- NULL
varsome_all_result_df['pharmgkb'] <- NULL

# 4. add "amp_annotation_verdict_tier" as a new column
varsome_all_result_df['amp_annotation_verdict_tier'] <- varsome_all_result_df$amp_annotation$verdict$tier
varsome_all_result_df['amp_annotation'] <- NULL

# add "acmg_annotation_verdict" as a new column
varsome_all_result_df['acmg_annotation_verdict'] <- varsome_all_result_df$acmg_annotation$verdict$ACMG_rules$verdict
varsome_all_result_df['acmg_annotation'] <- NULL
varsome_all_result_df['dbnsfp'] <- NULL

varsome_all_result_df <- varsome_all_result_df %>%
  mutate(ncbi_dbsnp = if_else(map_lgl(ncbi_dbsnp, is.null), 
                              list(list(version = NA, rsid = list(NA))), 
                              ncbi_dbsnp))
varsome_all_result_df <- varsome_all_result_df %>%
  unnest_wider(ncbi_dbsnp) %>%
  unnest(rsid, names_repair = "unique")
varsome_all_result_df <- varsome_all_result_df %>%
  select(-version)

varsome_all_result_df <- varsome_all_result_df %>%
  unnest(gnomad_exomes_coverage) %>% 
  unnest(coverage_mean) %>%
  unnest(coverage_median) %>%
  unnest(coverage_20_frequency) %>%
  rename_with(~ paste0("gnomad_exomes_", .), .cols = c("coverage_mean", "coverage_median","coverage_20_frequency"))  # 为所有列添加前缀
varsome_all_result_df['version'] <- NULL

varsome_all_result_df <- varsome_all_result_df %>%
  unnest(gnomad_genomes_coverage) %>% 
  unnest(coverage_mean) %>%
  unnest(coverage_median) %>%
  unnest(coverage_20_frequency) %>%
  rename_with(~ paste0("gnomad_genomes_", .), .cols = c("coverage_mean", "coverage_median","coverage_20_frequency"))
varsome_all_result_df['version'] <- NULL

# ncbi_clinvar2 database
varsome_all_result_df$ncbi_clinvar2_clinical_significance <- lapply(varsome_all_result_df$ncbi_clinvar2, function(x) {
  if (!is.null(x$clinical_significance) && length(x$clinical_significance) > 0) {
    # get clinical_significance
    clinical_significance_list <- x$clinical_significance[[1]]
    if (length(clinical_significance_list) > 1) {
      return(paste(clinical_significance_list, collapse = ","))
    } else {
      return(clinical_significance_list)
    }
  } else {
    return(NA)  # 如果 clinical_significance 为空或不存在，返回 NA
  }
})
varsome_all_result_df$ncbi_clinvar2_clinical_significance <- unlist(varsome_all_result_df$ncbi_clinvar2_clinical_significance)
varsome_all_result_df$ncbi_clinvar2 <- NULL
# varsome database
varsome_all_result_df$saphetor_known_pathogenicity <- lapply(varsome_all_result_df$saphetor_known_pathogenicity, function(x) {
  if (!is.null(x$items) && length(x$items) > 0 && !is.null(x$items[[1]]$annotations$`Saphetor PubMedUserEntry`)) {
    # get pathogenicity value
    pathogenicity_values <- sapply(x$items[[1]]$annotations$`Saphetor PubMedUserEntry`, function(entry) {
      entry$pathogenicity
    })
    # use comma connect each value
    return(paste(pathogenicity_values, collapse = ", "))
  } else {
    return(NA)  
  }
})
varsome_all_result_df$saphetor_known_pathogenicity <- unlist(varsome_all_result_df$saphetor_known_pathogenicity)

# civic -  返回中没有表示 pathogenetic的字段
clinical_significance_list <- lapply(varsome_all_result_df$wustl_civic, function(row) {
  if (!is.null(row$items) && length(row$items) > 0) {
    # 提取第一个 clinical_significance 或按需调整
    if (!is.null(row$items[[1]]$clinical_significance)) {
      return(row$items[[1]]$clinical_significance[1])  # 提取第一个值
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
})
varsome_all_result_df$civic_clinical_significance <- unlist(clinical_significance_list)
varsome_all_result_df$wustl_civic <- NULL

# read cosmic file - change file address to your local file address
cosmic_mutation_census_df <- read.delim('/Users/stan/Desktop/internship_project/database/params 1/Cosmic_MutantCensus_Tsv_v99_GRCh37/Cosmic_MutantCensus_v99_GRCh37.tsv', sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cosmic_sample_df <- read.delim('/Users/stan/Desktop/internship_project/database/params 1/Cosmic_Sample_Tsv_v99_GRCh37/Cosmic_Sample_v99_GRCh37.tsv', sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# cosmic_mutation_census_df: CHROMOSOME, GENOME_START, MUTATION_CDS, MUTATION_AA
# CHROMOSOME -> chromosome, GENOME_START -> pos, MUTATION_CDS -> hgvs, MUTATION_AA -> hgvs_p1
cosmic_mutation_census_df <- cosmic_mutation_census_df %>% select('CHROMOSOME', 'GENOME_START', 'MUTATION_CDS', 'MUTATION_AA', 'COSMIC_SAMPLE_ID')
cosmic_mutation_census_df$CHROMOSOME <- paste0("chr", cosmic_mutation_census_df$CHROMOSOME)
cosmic_sample_df <- cosmic_sample_df %>% select("COSMIC_SAMPLE_ID","SAMPLE_TYPE")

############### match  ensembl_transcripts and refseq_transcripts #######################################
# resolve ensembl_transcripts and refseq_transcripts matching problem

ensemb_refseq_resolve <- function(chr_df){
  ensembl_df <- chr_df %>% select("original_variant","chromosome","pos","ensembl_transcripts")
  refseq_df <- chr_df %>% select("original_variant","chromosome","pos","refseq_transcripts")
  # parse variables
  ensembl_df <- ensembl_df %>%
    unnest(ensembl_transcripts) %>%
    unnest(items) %>%
    rename_with(~ paste0("ensembl_", .), .cols = c("name", "strand", "coding_impact", "function", "hgvs", "hgvs_p1", "hgvs_p3", "location", "coding_location", "canonical", "gene_symbol", "splice_distance", "ensembl_support_level", "ensembl_appris", "mane_select", "mane_plus", "uniprot_id"))  # 为所有列添加前缀
  ensembl_df['version'] <- NULL
  ensembl_df <- subset(ensembl_df, select = -c(ensembl_splice_distance, ensembl_ensembl_appris, ensembl_mane_plus, ensembl_uniprot_id, ensembl_location, ensembl_coding_location)) 
  ensembl_df <- ensembl_df %>% filter(!is.na(ensembl_hgvs) & !is.na(ensembl_hgvs_p1) & !ensembl_hgvs_p1 == 'p.?')
  # parse variables
  refseq_df <- refseq_df %>%
    unnest(refseq_transcripts) %>%
    unnest(items) %>%
    rename_with(~ paste0("refseq_", .), .cols = c("name", "strand", "coding_impact", "function", "hgvs", "hgvs_p1", "hgvs_p3", "location", "coding_location", "canonical", "gene_symbol", "splice_distance", "ensembl_support_level", "ensembl_appris", "mane_select", "mane_plus", "uniprot_id"))  # 为所有列添加前缀
  refseq_df['version'] <- NULL
  refseq_df <- subset(refseq_df, select = -c(refseq_splice_distance, refseq_ensembl_appris, refseq_mane_plus, refseq_uniprot_id, refseq_location, refseq_coding_location)) 
  # filter rows which hgvs information is not empty 
  refseq_df <- refseq_df %>% filter(!is.na(refseq_hgvs) & !is.na(refseq_hgvs_p1) & !refseq_hgvs_p1 == 'p.?')
  # delete repeat rows
  ensembl_df <- unique(ensembl_df)
  refseq_df <- unique(refseq_df)
  # merge refseq_df and ensembl_df together by original_variant, hgvs, hgvs_p1
  matched_rows <- ensembl_df %>%
    inner_join(refseq_df,
               by = c("original_variant" = "original_variant",
                      "ensembl_hgvs" = "refseq_hgvs",
                      "ensembl_hgvs_p1" = "refseq_hgvs_p1"),
               relationship = "many-to-many")
  # save refseq information for cosmic connect using
  matched_rows <- matched_rows %>% mutate(refseq_hgvs = ensembl_hgvs,
                                          refseq_hgvs_p1 = ensembl_hgvs_p1)
  matched_rows <- matched_rows %>% select(-chromosome.y,-pos.y) %>% rename('chromosome' = 'chromosome.x','pos' = 'pos.x')
  # delete repeat rows
  matched_rows <- matched_rows %>% distinct(original_variant,ensembl_name, ensembl_hgvs, ensembl_hgvs_p1, .keep_all = TRUE)
  # select rows which didnt match
  ensembl_not_in_refseq <- anti_join(
    ensembl_df, refseq_df,
    by = c("original_variant" = "original_variant", 
           "ensembl_hgvs" = "refseq_hgvs", 
           "ensembl_hgvs_p1" = "refseq_hgvs_p1")
  )
  refseq_not_in_ensembl <- anti_join(
    refseq_df, ensembl_df,
    by = c("original_variant" = "original_variant", 
           "refseq_hgvs" = "ensembl_hgvs", 
           "refseq_hgvs_p1" = "ensembl_hgvs_p1")
  )
  # 针对ensemble 和 refsqe 的共有的和 特有的数据都进行cosmic 匹配，然后再将所有的数据合并到hg19中去
  # search function defining
  sample_match <- function(type,df){
    x_colum <- list()
    merged_df <- NULL
    if(type == 'matched' | type == 'ensembl'){
      df$ensembl_hgvs_p1 <- ifelse(is.na(df$ensembl_hgvs_p1), NA, paste0("p.", df$ensembl_hgvs_p1)) # set same hgvsp. format with cosmic
      x_colum <- c("chromosome", "pos", "ensembl_hgvs", "ensembl_hgvs_p1")
    }else{
      x_colum <- c("chromosome", "pos", "refseq_hgvs", "refseq_hgvs_p1")
      df$refseq_hgvs_p1 <- ifelse(is.na(df$refseq_hgvs_p1), NA, paste0("p.", df$refseq_hgvs_p1))
    }
    merged_df <- merge(df, cosmic_mutation_census_df,
                       by.x = x_colum,
                       by.y = c("CHROMOSOME", "GENOME_START", "MUTATION_CDS", "MUTATION_AA"),
                       all.x = TRUE)
    merged_df <- unique(merged_df) #delete repeat rows
    # merge sample file and match sample_type to each variant
    merged_df <- merge(merged_df, cosmic_sample_df,
                       by.x = c("COSMIC_SAMPLE_ID"),
                       by.y = c("COSMIC_SAMPLE_ID"),
                       all.x = TRUE)
    print(7)
    return(merged_df)
  }
  #  1. 对matched 数据进行smaple type的查询
  # matched_rows_add_sample <- NULL
  # ensembl_not_in_refseq_add_sample <- NULL
  # refseq_not_in_ensembl_add_sample <- NULL
  matched_rows_add_sample <- sample_match('matched', matched_rows)
  #  2. 对ensemble特有的数据进行sample type的查询
  ensembl_not_in_refseq_add_sample <- sample_match('ensembl', ensembl_not_in_refseq)
  #  3. 对refuse特有的数据进行sample type的查询
  refseq_not_in_ensembl_add_sample <- sample_match('refseq', refseq_not_in_ensembl)
  
  transfer_df <- bind_rows(matched_rows_add_sample, ensembl_not_in_refseq_add_sample, refseq_not_in_ensembl_add_sample)
  
  chr_df <- chr_df %>% inner_join(transfer_df, by = c("original_variant" = "original_variant"),
                                  relationship = "many-to-many") %>% select(-refseq_transcripts, -ensembl_transcripts)
  chr_df <- unique(chr_df) #delete repeat rows
  return(chr_df)
}

# separate all variant to different group based on different chromosome position, 
# this step for parse ensemble and refseq data, because this two database will make database bigger more
split_dfs <- split(varsome_all_result_df, varsome_all_result_df$chromosome)
for (i in names(split_dfs)) {
  df_name <- paste0(i, "_df")
  as.data.frame(assign(df_name, split_dfs[[i]])) 
}

chr1_df <- ensemb_refseq_resolve(chr1_df) 
chr1_df <- chr1_df %>% select(-chromosome.y,COSMIC_SAMPLE_ID) %>% rename(chromosome = chromosome.x)
chr2_df <- ensemb_refseq_resolve(chr2_df)
chr2_df <- chr2_df %>% select(-chromosome.y,COSMIC_SAMPLE_ID) %>% rename(chromosome = chromosome.x)
chr3_df <- ensemb_refseq_resolve(chr3_df)
chr3_df <- chr3_df %>% select(-chromosome.y,COSMIC_SAMPLE_ID) %>% rename(chromosome = chromosome.x)
chr4_df <- ensemb_refseq_resolve(chr4_df)
chr4_df <- chr4_df %>% select(-chromosome.y,COSMIC_SAMPLE_ID) %>% rename(chromosome = chromosome.x)
# chr5_df <- ensemb_refseq_resolve(chr5_df) #----- 有问题，需要检查内容
chr6_df <- ensemb_refseq_resolve(chr6_df)
chr6_df <- chr6_df %>% select(-chromosome.y,COSMIC_SAMPLE_ID) %>% rename(chromosome = chromosome.x)
chr7_df <- ensemb_refseq_resolve(chr7_df)
chr7_df <- chr7_df %>% select(-chromosome.y,COSMIC_SAMPLE_ID) %>% rename(chromosome = chromosome.x)
# chr8_df <- ensemb_refseq_resolve(chr8_df) #----- 有问题，需要检查内容

# resolve missed list column ---- refseq_function
resolve_list_problem <- function(df){
  df[] <- lapply(df, function(x) {
    if (is.list(x)) {
      sapply(x, function(y) paste(y, collapse = ","))
    } else {
      x
    }
  })
  return(df)
}
chr_df_all <- list(chr1_df, chr2_df, chr3_df,chr4_df, chr6_df, chr7_df)  
resolve_every_df <- lapply(chr_df_all, resolve_list_problem)
write.csv(chr1_df, file = "/Users/stan/Desktop/internship_project/project/data_collection_table/chr1_df.csv", row.names = FALSE)
for (i in seq_along(resolve_every_df)) {
  write.csv(resolve_every_df[[i]], file = paste0("/Users/stan/Desktop/internship_project/project/data_collection_table/chr", i, "_df.csv"), row.names = FALSE)
}

source('Rcode/sift_connect.R')
db_check_fun('SIFT',chr1_df)
db_check_fun('SIFT',chr2_df)
db_check_fun('SIFT',chr3_df)
db_check_fun('SIFT',chr4_df)
db_check_fun('SIFT',chr5_df)
db_check_fun('SIFT',chr6_df)
db_check_fun('SIFT',chr7_df)
db_check_fun('SIFT',chr8_df)
# chr2_df$sift_prediction <- NULL
chr1_df <- sift_pre_function(chr1_df)
chr2_df <- sift_pre_function(chr2_df)
chr2_df <- sift_pre_function(chr3_df)
chr2_df <- sift_pre_function(chr4_df)
chr2_df <- sift_pre_function(chr5_df)
chr2_df <- sift_pre_function(chr6_df)
chr2_df <- sift_pre_function(chr7_df)
chr2_df <- sift_pre_function(chr8_df)
print(colnames(chr1_df))
# connect mutation master, get mutation master prediction


