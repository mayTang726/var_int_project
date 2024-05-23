# this page just for separately reading hg19 .bed file and hg38 .bed file
# because the VarSome API already expired, so I stored them as a file, and 
# now just need to collect useful information 

getwd()
# 1. read json file
# 2. parse json file 
# 3. connect hg19 and hg38 details
# 4. select useful information
library(jsonlite)
library(tidyr)
library(dplyr)
library(purrr)

# resolve hg19
hg19_info <- 'data/variant_response_hg19.json'
hg38_info <- 'data/variant_response_hg38.json'

hg19_df <- fromJSON(hg19_info)
hg38_df <- fromJSON(hg38_info)

print(colnames(hg19_df))

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

hg19_df['amp_annotation_verdict_tier'] <- hg19_df$amp_annotation$verdict$tier
hg19_df['amp_annotation'] <- NULL
hg19_df['acmg_annotation_verdict'] <- hg19_df$acmg_annotation$verdict$ACMG_rules$verdict
hg19_df['acmg_annotation'] <- NULL

hg19_df <- hg19_df %>%
  unnest(refseq_transcripts) %>%
  unnest(items) %>%
  rename_with(~ paste0("refseq_", .), .cols = c("name", "strand", "coding_impact", "function", "hgvs", "hgvs_p1", "hgvs_p3", "location", "coding_location", "canonical", "gene_symbol", "splice_distance", "ensembl_support_level", "ensembl_appris", "mane_select", "mane_plus", "uniprot_id"))  # 为所有列添加前缀
hg19_df['version'] <- NULL


hg19_df <- hg19_df %>%
  unnest(ensembl_transcripts) %>% 
  unnest(items) %>%
  rename_with(~ paste0("ensembl_", .), .cols = c("name", "strand", "coding_impact", "function", "hgvs", "hgvs_p1", "hgvs_p3", "location", "coding_location", "canonical", "gene_symbol", "splice_distance", "ensembl_support_level", "ensembl_appris", "mane_select", "mane_plus", "uniprot_id"))  # 为所有列添加前缀
hg19_df['version'] <- NULL

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

hg19_df <- hg19_df %>%
  mutate(dbnsfp = map_if(dbnsfp, ~ !is.null(.x), ~ transmute(.x, mutationtaster_pred, mutationtaster_score, sift_score, sift_pred))) %>%
  unnest(dbnsfp)
hg19_df['version'] <- NULL
print(colnames(hg19_df))


# resolve hg38
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
  mutate(dbnsfp = map_if(dbnsfp, ~ !is.null(.x), ~ transmute(.x, mutationtaster_pred, mutationtaster_score, sift_score, sift_pred))) %>%
  unnest(dbnsfp)
hg38_df['version'] <- NULL
print(colnames(hg19_df))

