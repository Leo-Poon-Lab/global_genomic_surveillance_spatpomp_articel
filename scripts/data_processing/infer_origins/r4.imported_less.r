library(dplyr)
library(lubridate)
args<-commandArgs(TRUE)
usher_result2<-read.table(args[1],head=T)
# Ensure your date column is in the correct date format
usher_result2$date <- as.Date(usher_result2$date)
usher_result2$Imported_modify<-usher_result2$V14
usher_result2 <- usher_result2 %>%
  select(strain,date, date_decimal,node, code,country_infer,country,virus,V14) %>%
  arrange(date, node, code) %>%
  group_by(node,code) %>%
  mutate(date_diff = date - lag(date),
         Imported_modify = ifelse(is.na(date_diff) | date_diff >= 30 | V14 != "Imported", V14, "Community")) %>%
  ungroup()
write.table(usher_result2,args[2],quote=FALSE,row.names = FALSE)
