# Read response from varsome variant request based on hg19 and hg38  
# in the last step , i stored the response as json file, so that we dont need to request every time, except the detail of .bed file changed
# change json file to dataframe, and unset the attribute to get the useful information

getwd()

library(jsonlite)
library(tidyr)
library(dplyr)
library(purrr)
# 当前是基于chr17 做的，需要将此部分改为循环

# 1. read json file
hg19_info <- 'data/variant_response_hg19.json'
# hg38_info <- 'data/variant_response_hg38.json'

# 2. parse json file -> dataframe
hg19_df <- fromJSON(hg19_info)
# hg38_df <- fromJSON(hg38_info)
print(colnames(hg19_df))

# resolve hg19 response
# 3. delete unuseful columns
## if you want use the information below, you can add comment for the line which you wanna use
hg19_df <- hg19_df[!is.na(hg19_df$chromosome),]

hg19_df['regions'] <- NULL 
hg19_df['gerp'] <- NULL 
hg19_df['phastcons100way'] <- NULL 
hg19_df['phylop100way'] <- NULL 
hg19_df['maxentscan'] <- NULL 
hg19_df['publications'] <- NULL
hg19_df['cbio_portal'] <- NULL
hg19_df['detail'] <- NULL
hg19_df['weill_cornell_medicine_pmkb'] <- NULL
hg19_df['icgc_somatic'] <- NULL
hg19_df['nih_gdc'] <- NULL
hg19_df['wustl_docm'] <- NULL
hg19_df['variant_pubmed_automap'] <- NULL
hg19_df['cadd'] <- NULL
hg19_df['dann_snvs'] <- NULL
hg19_df['error'] <- NULL

# 4. add "amp_annotation_verdict_tier" as a new column
hg19_df['amp_annotation_verdict_tier'] <- hg19_df$amp_annotation$verdict$tier
hg19_df['amp_annotation'] <- NULL

# add "acmg_annotation_verdict" as a new column
hg19_df['acmg_annotation_verdict'] <- hg19_df$acmg_annotation$verdict$ACMG_rules$verdict
hg19_df['acmg_annotation'] <- NULL
hg19_df['dbnsfp'] <- NULL

hg19_df <- hg19_df %>%
  mutate(ncbi_dbsnp = if_else(map_lgl(ncbi_dbsnp, is.null), 
                              list(list(version = NA, rsid = list(NA))), 
                              ncbi_dbsnp))
hg19_df <- hg19_df %>%
  unnest_wider(ncbi_dbsnp) %>%
  unnest(rsid, names_repair = "unique")
hg19_df <- hg19_df %>%
  select(-version)

hg19_df <- hg19_df %>%
  unnest(gnomad_exomes_coverage) %>% 
  unnest(coverage_mean) %>%
  unnest(coverage_median) %>%
  unnest(coverage_20_frequency) %>%
  rename_with(~ paste0("gnomad_exomes_", .), .cols = c("coverage_mean", "coverage_median","coverage_20_frequency"))  # 为所有列添加前缀
hg19_df['version'] <- NULL

hg19_df <- hg19_df %>%
  unnest(gnomad_genomes_coverage) %>% 
  unnest(coverage_mean) %>%
  unnest(coverage_median) %>%
  unnest(coverage_20_frequency) %>%
  rename_with(~ paste0("gnomad_genomes_", .), .cols = c("coverage_mean", "coverage_median","coverage_20_frequency"))
hg19_df['version'] <- NULL

# 5. unset the columns which have multiple layers
############### match  ensembl_transcripts and refseq_transcripts #######################################
# resolve ensembl_transcripts and refseq_transcripts matching problem

ensembl_df <- hg19_df %>% select("original_variant","chromosome","pos","ensembl_transcripts")
refseq_df <- hg19_df %>% select("original_variant","chromosome","pos","refseq_transcripts")
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


###############  connect cosmic #######################################
# read cosmic file - change file address to your local file address
cosmic_mutation_census_df <- read.delim('/Users/stan/Desktop/internship_project/database/params 1/Cosmic_MutantCensus_Tsv_v99_GRCh37/Cosmic_MutantCensus_v99_GRCh37.tsv', sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cosmic_sample_df <- read.delim('/Users/stan/Desktop/internship_project/database/params 1/Cosmic_Sample_Tsv_v99_GRCh37/Cosmic_Sample_v99_GRCh37.tsv', sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# cosmic_mutation_census_df: CHROMOSOME, GENOME_START, MUTATION_CDS, MUTATION_AA
# CHROMOSOME -> chromosome, GENOME_START -> pos, MUTATION_CDS -> hgvs, MUTATION_AA -> hgvs_p1
cosmic_mutation_census_df <- cosmic_mutation_census_df %>% select('CHROMOSOME', 'GENOME_START', 'MUTATION_CDS', 'MUTATION_AA', 'COSMIC_SAMPLE_ID')
cosmic_mutation_census_df$CHROMOSOME <- paste0("chr", cosmic_mutation_census_df$CHROMOSOME)
cosmic_sample_df <- cosmic_sample_df %>% select("COSMIC_SAMPLE_ID","SAMPLE_TYPE")

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

hg19_df <- hg19_df %>% inner_join(transfer_df, by = c("original_variant" = "original_variant"),
                                  relationship = "many-to-many") %>% 
  select(-refseq_transcripts, -ensembl_transcripts)
hg19_df <- unique(hg19_df) #delete repeat rows








###############  resolve hg38 #######################################
# resolve hg38 response (same step as hg19)
hg38_df['regions'] <- NULL 
hg38_df['gerp'] <- NULL 
hg38_df['phastcons100way'] <- NULL 
hg38_df['phylop100way'] <- NULL 
hg38_df['maxentscan'] <- NULL 
hg38_df['publications'] <- NULL
hg38_df['cbio_portal'] <- NULL
hg38_df['detail'] <- NULL
hg38_df['weill_cornell_medicine_pmkb'] <- NULL
hg38_df['icgc_somatic'] <- NULL
hg38_df['nih_gdc'] <- NULL
hg38_df['wustl_docm'] <- NULL
hg38_df['variant_pubmed_automap'] <- NULL
hg19_df['dbnsfp'] <- NULL

hg38_df['amp_annotation_verdict_tier'] <- hg38_df$amp_annotation$verdict$tier
hg38_df['amp_annotation'] <- NULL
hg38_df['acmg_annotation_verdict'] <- hg38_df$acmg_annotation$verdict$ACMG_rules$verdict
hg38_df['acmg_annotation'] <- NULL

hg38_df <- hg38_df %>%
  unnest(refseq_transcripts) %>%
  unnest(items) %>%
  rename_with(~ paste0("refseq_", .), .cols = c("name", "strand", "coding_impact", "function", "hgvs", "hgvs_p1", "hgvs_p3", "location", "coding_location", "canonical", "gene_symbol", "splice_distance", "ensembl_support_level", "ensembl_appris", "mane_select", "mane_plus", "uniprot_id"))  # 为所有列添加前缀
hg38_df['version'] <- NULL


hg38_df <- hg38_df %>%
  unnest(ensembl_transcripts) %>% 
  unnest(items) %>%
  rename_with(~ paste0("ensembl_", .), .cols = c("name", "strand", "coding_impact", "function", "hgvs", "hgvs_p1", "hgvs_p3", "location", "coding_location", "canonical", "gene_symbol", "splice_distance", "ensembl_support_level", "ensembl_appris", "mane_select", "mane_plus", "uniprot_id"))  # 为所有列添加前缀
hg38_df['version'] <- NULL

hg38_df <- hg38_df %>%
  unnest(gnomad_exomes_coverage) %>% 
  unnest(coverage_mean) %>%
  unnest(coverage_median) %>%
  unnest(coverage_20_frequency) %>%
  rename_with(~ paste0("gnomad_exomes_", .), .cols = c("coverage_mean", "coverage_median","coverage_20_frequency"))  # 为所有列添加前缀
hg38_df['version'] <- NULL

hg38_df <- hg38_df %>%
  unnest(gnomad_genomes_coverage) %>% 
  unnest(coverage_mean) %>%
  unnest(coverage_median) %>%
  unnest(coverage_20_frequency) %>%
  rename_with(~ paste0("gnomad_genomes_", .), .cols = c("coverage_mean", "coverage_median","coverage_20_frequency"))
hg38_df['version'] <- NULL

hg38_df <- hg38_df %>%
  mutate(ncbi_dbsnp = if_else(map_lgl(ncbi_dbsnp, is.null), 
                              list(list(version = NA, rsid = list(NA))), 
                              ncbi_dbsnp))
hg38_df <- hg38_df %>%
  unnest_wider(ncbi_dbsnp) %>%
  unnest(rsid, names_repair = "unique")
hg38_df <- hg38_df %>%
  select(-version)
print(colnames(hg38_df))

hg38_df$original_variant <- gsub("\\-", ":", hg38_df$original_variant)
