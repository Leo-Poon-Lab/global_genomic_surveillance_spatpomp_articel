#r
library(tidyverse)

library(ape)
library(treeio)
library(tidytree)

tree<-read.nhx("result/infer_origin/BA1/data/usher/imported_align/imported/mugration/BA1_final_annotated_tree.nexus")
tree<-as_tibble(tree)
id<-tree[grep("hCoV-19/", tree$label),"label"]
results<-as.data.frame(id)
metadata<-read.table("~/Project/04.SARS-COV-2/01.global_epi/00.data/BA2/data/metadata.tsv",head=T)
for (i in 1:nrow(results)) {
    label <- results$label[i]
    parent_node <- parent(tree, label)
    results$parent[i] <- parent_node$code
    results$node[i] <- parent_node$node
    matching_row <- which(tree$label == label)
   # results$country[i] <- tree$code[matching_row]
    results$country[i] <- metadata[metadata$strain == label, "code"]
    results$date[i] <- metadata[metadata$strain == label, "date_decimal"]
}
results[results$parent==results$country,5]<-"Community"
results[results$parent!=results$country,5]<-"Imported"
write.table(results,"result/infer_origin/BA1/data/usher/imported_align/imported/imported_result.txt",quote=FALSE,row.names = FALSE)
