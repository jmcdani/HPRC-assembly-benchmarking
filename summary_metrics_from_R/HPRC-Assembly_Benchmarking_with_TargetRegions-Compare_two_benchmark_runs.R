##########################################
## This script is to summarize the output metrics from hap.py results_extended.csv using two different dipcall contig length for comparison.  
## This script assumes hap.py was run with the following hap.py parameters: 
## with and without target-region and stratification with GIAB stratification BEDs.
############################################

library(tidyverse)
library(here)
library(fs)
library(ggplot2)

# directory paths for results_extended.csv files
data_dir10 <- "/Users/jmcdani/Documents/GiaB/Benchmarking/assembly-benchmarking/hpp/summary_metrics_from_R/ext_results_for_R/10kb/"
data_dir50 <- "/Users/jmcdani/Documents/GiaB/Benchmarking/assembly-benchmarking/hpp/summary_metrics_from_R/ext_results_for_R/50kb/"
data_dir10_nontargeted <- "/Users/jmcdani/Documents/GiaB/Benchmarking/assembly-benchmarking/hpp/summary_metrics_from_R/ext_results_for_R/10kb_nontargeted/"

############################################
## Read in non-targeted summary metrics for use with calculating targeted summary metrics. The targeted calculated recall metrics use the total number of 
## benchmark variants as the denominator (non-targeted) before targeting to the dip.bed, so it gives an idea of how many variants are excluded by the dip.bed
############################################
## funtion to read in only the metrics needed to calculate new targeted metrics
## nontargeted_ext_results(dir = path to directory w/ non-targeted extended results, contig = length of contig used for threshold in dipcall)
############################################

nontargeted_ext_results <- function(dir, contig) {
  extended_files <- fs::dir_ls(dir, regexp = "extended.csv$")
  ext_results <- extended_files %>%  set_names() %>% map_dfr(read_csv, .id="source") %>%
    filter(Filter == "PASS", Subtype == "*")  %>%
    select(c(source,
             Filter,
             Type,
             Subtype,
             Subset, 
             TRUTH.TOTAL.het,
             TRUTH.TOTAL.homalt,
             TRUTH.TOTAL)) %>% 
    rename(assembly = source,
           ntTRUTH.TOTAL.homalt = TRUTH.TOTAL.homalt,
           ntTRUTH.TOTAL.het = TRUTH.TOTAL.het,
           ntTRUTH.TOTAL = TRUTH.TOTAL)
  ext_results$assembly <- str_split(ext_results$assembly, "/") %>% sapply(tail,1)
  ext_results$assembly <- str_split(ext_results$assembly, "_") %>% sapply (head, 1)
  return(ext_results)
}
#mutate(dipcall_contig = contig) %>% 
#ext_results <- ext_results[,c(1,9,2,3,4,5,6,7,8)] 

nt10kb <- nontargeted_ext_results(data_dir10_nontargeted, 10)

############################################
## Read in targeted summary metrics for use with calculating targeted summary metrics. Only Filter == "PASS" metrics are needed since
## filtered sites should be outside dip.bed used for targeting. 
###########################################
## function to read in only the metrics needed to calculate new targeted metrics
## targeted_ext_results(dir = path to directory w/ targeted extended results, contig = length of contig used for threshold in dipcall)
############################################

targeted_ext_results <- function(dir, contig) {
  extended_files <- fs::dir_ls(dir, regexp = "extended.csv$")
  ext_results <- extended_files %>%  set_names() %>% map_dfr(read_csv, .id="source") %>%
    filter(Filter == "PASS", Subtype == "*")  %>%
    select(c(source,
             Filter,
             Type,
             Subtype,
             Subset, 
             QUERY.FP, 
             Subset.IS_CONF.Size, 
             FP.gt, 
             FP.al,
             METRIC.Recall,
             TRUTH.TOTAL,
             TRUTH.TP,
             TRUTH.TP.het,
             TRUTH.TP.homalt,
             TRUTH.TOTAL.het,
             TRUTH.TOTAL.homalt)) %>% 
    mutate(dipcall_contig = contig) %>% 
    rename(assembly = source)
  ext_results <- ext_results[,c(1,17,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)] 
  ext_results$assembly <- str_split(ext_results$assembly, "/") %>% sapply(tail,1)
  ext_results$assembly <- str_split(ext_results$assembly, "_") %>% sapply (head, 1)
  return(ext_results)
}

## run targeted_ext_results() for each contig. Join separately to non-targeted, same non-targeted values used for both contig lengths.
## Once targeted/non-targeted joined, combine the two dataframes
targeted_ext_results10 <- targeted_ext_results(data_dir10, 10)
targeted_ext_results50 <- targeted_ext_results(data_dir50, 50)

ext_results10 <- left_join(targeted_ext_results10, nt10kb, 
                           by = c("assembly", "Filter", "Type", "Subtype", "Subset"))
ext_results50 <- left_join(targeted_ext_results50, nt10kb, 
                           by = c("assembly", "Filter", "Type", "Subtype", "Subset"))
ext_results_combined <- bind_rows(ext_results10, ext_results50)

## Create list of all "Subset" stratifications for which metrics will be calculated
subsets <- c("*" ,
             "GRCh38_HG002_GIABv4.1_notin_complexandSVs_alldifficultregions.bed.gz" ,
             "GRCh38_notinalldifficultregions.bed.gz",
             "GRCh38_notinalllowmapandsegdupregions.bed.gz",
             "GRCh38_notinsegdups.bed.gz",
             "GRCh38_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz",
             "GRCh38_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
             "GRCh38_segdups.bed.gz",
             "GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz",
             "GRCh38_MHC.bed.gz")

############################################
## Filter the orginial data frame for Type (SNP/INDEL) and rename colnames to note metric Type in larger table 
## once combined. Tables are then joined by colnames in common.
############################################
## fuction to rename colnames to specific metric and combine metric types into single table where all type.metric
## cols are present. 
## type_rename(combined_results = table where ext result from runs where different dipcall contig lengths were used)
############################################

type_rename <- function(combined_results) {
SNP_ext_results<- combined_results %>% 
  filter(Type == "SNP" & Subset %in% subsets) %>% 
  rename(SNP.QUERY.FP = QUERY.FP, 
         SNP.Subset.IS_CONF.Size = Subset.IS_CONF.Size, 
         SNP.FP.gt = FP.gt, 
         SNP.FP.al = FP.al, 
         SNP.Recall = METRIC.Recall, 
         SNP.TRUTH.TOTAL= TRUTH.TOTAL, 
         SNP.TRUTH.TP = TRUTH.TP,
         SNP.TRUTH.TP.het = TRUTH.TP.het,
         SNP.TRUTH.TP.homalt = TRUTH.TP.homalt,
         SNP.TRUTH.TOTAL.het = TRUTH.TOTAL.het,
         SNP.TRUTH.TOTAL.homalt = TRUTH.TOTAL.homalt,
         SNP.ntTRUTH.TOTAL.het = ntTRUTH.TOTAL.het,
         SNP.ntTRUTH.TOTAL.homalt = ntTRUTH.TOTAL.homalt,
         SNP.ntTRUTH.TOTAL = ntTRUTH.TOTAL) %>%
  select(-Type)

INDEL_ext_results <- combined_results %>% 
  filter(Type == "INDEL" & Subset %in% subsets ) %>%
  rename(INDEL.QUERY.FP = QUERY.FP,
         INDEL.Subset.IS_CONF.Size = Subset.IS_CONF.Size,
         INDEL.FP.gt = FP.gt,
         INDEL.FP.al = FP.al,
         INDEL.Recall = METRIC.Recall,
         INDEL.TRUTH.TOTAL = TRUTH.TOTAL,
         INDEL.TRUTH.TP = TRUTH.TP,
         INDEL.TRUTH.TP.het = TRUTH.TP.het,
         INDEL.TRUTH.TP.homalt = TRUTH.TP.homalt,
         INDEL.TRUTH.TOTAL.het = TRUTH.TOTAL.het,
         INDEL.TRUTH.TOTAL.homalt = TRUTH.TOTAL.homalt,
         INDEL.ntTRUTH.TOTAL.het = ntTRUTH.TOTAL.het,
         INDEL.ntTRUTH.TOTAL.homalt = ntTRUTH.TOTAL.homalt,
         INDEL.ntTRUTH.TOTAL = ntTRUTH.TOTAL) %>%
  select(-Type)

ext_results_combined_renamed <- left_join(INDEL_ext_results,
                                          SNP_ext_results,
                                          by = c("assembly", "dipcall_contig", "Filter", "Subtype", "Subset"))
return(ext_results_combined_renamed)
}

ext_results_combined_renamed <- type_rename(ext_results_combined)

############################################
## Calcluate new summary metrics provided by JZ
############################################
## Function calculates new summary metrics using values form exteneded results
## summary_metrics(results = DF that provides the appropriate fields for which summary metics should be calculated)
############################################

summary_metrics <- function(results) {                                                                    
  metrics_table <- results %>%
    ## all QV metric denominators multiplied by 2 to be more representative of diploid value
    mutate(QV_dip_snp = -10*log10(SNP.QUERY.FP/(2 * SNP.Subset.IS_CONF.Size))) %>%
    mutate(QV_dip_indel = -10*log10(INDEL.QUERY.FP/(2 * INDEL.Subset.IS_CONF.Size))) %>%
    mutate(QV_dip_snp_indel = -10*log10((SNP.QUERY.FP + INDEL.QUERY.FP)/(2 * SNP.Subset.IS_CONF.Size))) %>%
    mutate(QV_hap_snp = -10*log10((SNP.QUERY.FP - SNP.FP.gt -SNP.FP.al)/(2 * SNP.Subset.IS_CONF.Size))) %>%
    mutate(QV_ignoreGT_snp = -10*log10((SNP.QUERY.FP - SNP.FP.gt)/(2 * SNP.Subset.IS_CONF.Size))) %>%
    mutate(QV_hap_indel = -10*log10((INDEL.QUERY.FP - INDEL.FP.gt - INDEL.FP.al)/(2 * INDEL.Subset.IS_CONF.Size))) %>%
    mutate(QV_ignoreGT_indel = -10*log10((INDEL.QUERY.FP - INDEL.FP.gt)/(2 * INDEL.Subset.IS_CONF.Size))) %>%
    mutate(QV_hap_snp_indel = -10*log10(((SNP.QUERY.FP - SNP.FP.gt - SNP.FP.al) + (INDEL.QUERY.FP - INDEL.FP.gt - INDEL.FP.al))/(2 * SNP.Subset.IS_CONF.Size))) %>%
    mutate(QV_ignoreGT_snp_indel = -10*log10(((SNP.QUERY.FP - SNP.FP.gt) + (INDEL.QUERY.FP - INDEL.FP.gt))/(2 * SNP.Subset.IS_CONF.Size))) %>%
    mutate(SNP.Recall_ignoreGT = (SNP.TRUTH.TP + SNP.FP.gt) / SNP.TRUTH.TOTAL) %>%
    mutate(INDEL.Recall_ignoreGT = (INDEL.TRUTH.TP + INDEL.FP.gt) / INDEL.TRUTH.TOTAL) %>%
    mutate(SNP.Recall.het = SNP.TRUTH.TP.het / SNP.TRUTH.TOTAL.het) %>%
    mutate(SNP.Recall.hom = SNP.TRUTH.TP.homalt / SNP.TRUTH.TOTAL.homalt) %>%
    mutate(INDEL.Recall.het = SNP.TRUTH.TP.het / SNP.TRUTH.TOTAL.het) %>%
    mutate(INDEL.Recall.hom = SNP.TRUTH.TP.homalt / SNP.TRUTH.TOTAL.homalt) %>%
    ## for all non-targeted metrics (SNP.nt) the dipcall 10kb output will be used for both 10kb and 50kb calculated metrics. JZ notes these should be the same
    mutate(SNP.Recall.het.TOTALnontarget = SNP.TRUTH.TP.het / SNP.ntTRUTH.TOTAL.het) %>%
    mutate(SNP.Recall.hom.TOTALnontarget = SNP.TRUTH.TP.homalt / SNP.ntTRUTH.TOTAL.homalt) %>%
    mutate(INDEL.Recall.het.TOTALnontarget = SNP.TRUTH.TP.het / SNP.ntTRUTH.TOTAL.het) %>%
    mutate(INDEL.Recall.hom.TOTALnontarget = SNP.TRUTH.TP.homalt / SNP.ntTRUTH.TOTAL.homalt) %>%
    mutate(SNP.Recall_ignoreGT.TOTALnontarget = (SNP.TRUTH.TP + SNP.FP.gt) / SNP.ntTRUTH.TOTAL) %>%
    mutate(SNP.Recall.TOTALnontarget = (SNP.TRUTH.TP) / SNP.ntTRUTH.TOTAL)
  return(metrics_table)
}

benchmarking_summary <- summary_metrics(ext_results_combined_renamed)

write.csv(benchmarking_summary, "dipcall_10kbAND50kb_happy312_with_gap2homvarbutfiltered_target_regions.csv")
