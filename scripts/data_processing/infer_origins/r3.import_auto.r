library(tidyverse)
# Define your base directory and file name
base_dir <- "results/BA1/data/usher/newsample"  # replace with your actual base directory
file_name <- "import_result1"  # replace with your actual file name

# Find all files with the specified name in the base directory and its subdirectories
file_paths <- list.files(base_dir, pattern = paste0("^", "parsimony-scores2.tsv", "$"), full.names = TRUE, recursive = TRUE)
country<-read.csv("results/BA1/data/usher/mugration/colname.csv",head=F)
nodes_country_all<-read.csv("results/BA1/data/usher/mugration/confidence.csv",head=T)
rownames(nodes_country_all)<-nodes_country_all[,1]
nodes_country_all<-nodes_country_all[,-1]
colnames(nodes_country_all)<-t(country)
ID <- nodes_country_all[grep("hCoV-19", rownames(nodes_country_all)), ]
metadata <- read.table("results/BA1/data/metadata.tsv",head=TRUE)
FUN <- function(x) { country[which.max(x), ] }

# Iterate over each file in the list
for (file_path in file_paths) {
  usher_result <- read.table(file_path,head=F)
  usher_result$country_infer<-apply(nodes_country_all[usher_result$V2,],1,FUN)
colnames(usher_result)<-c("strain","node","country_infer")
usher_result2<-left_join(usher_result,metadata,by="strain")
usher_result2[usher_result2$country_infer==usher_result2$code,14]<-"Community"
usher_result2[usher_result2$country_infer!=usher_result2$code,14]<-"Imported"
usher_result2[which(duplicated(usher_result2[,c("node","code")])),14]<-"Community"
  # Write the output to a file in the same directory as the input file, with a different name
output_file <- sub("\\.tsv$", "_import_result.csv", file_path)
write.table(usher_result2, output_file, quote = FALSE, row.names = FALSE)
}
