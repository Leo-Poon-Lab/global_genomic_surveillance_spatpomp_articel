##
library(tidyverse)
library(parallel)
library(here)
n_cores = 16
setwd("global_genomic_surveillance_spatpomp/")
##

BA1<-read_tsv("results/model_data/origin_infer/vn_clean.csv")
colnames(BA1)<-"Virus name"
gisaid_BA1<-left_join(BA1,df_meta,by="Virus name")
colnames(gisaid_BA1)<-c("strain","virus","date","host","clade","lineage","continent","country","region","date_decimal","code")
write_tsv(gisaid_BA1, "results/model_data/origin_infer/vn_clean_gisaid.tsv")

BA2<-read_tsv("results/model_data/origin_infer/vn_clean.csv")
colnames(BA2)<-"Virus name"
gisaid_BA2<-left_join(BA2,df_meta,by="Virus name")
colnames(gisaid_BA2)<-c("strain","virus","date","host","clade","lineage","continent","country","region","date_decimal","code")
write_tsv(gisaid_BA2, "results/model_data/origin_infer/vn_clean_gisaid.tsv")

others<-read_tsv("results/model_data/origin_infer/vn_clean.csv")
colnames(others)<-"Virus name"
gisaid_others<-left_join(others,df_meta,by="Virus name")
colnames(gisaid_others)<-c("strain","virus","date","host","clade","lineage","continent","country","region","date_decimal","code")
write_tsv(gisaid_others, "results/model_data/origin_infer/vn_clean_gisaid.tsv")


Delta<-read_tsv("results/model_data/origin_infer/vn_clean.csv")
colnames(Delta)<-"Virus name"
gisaid_Delta<-left_join(Delta,df_meta,by="Virus name")
colnames(gisaid_Delta)<-c("strain","virus","date","host","clade","lineage","continent","country","region","date_decimal","code")
write_tsv(gisaid_Delta, "results/model_data/origin_infer/vn_clean_gisaid.tsv")


