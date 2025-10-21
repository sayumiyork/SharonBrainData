#Load libraries
library("tidyverse")
library("DESeq2")

#processed files
brain<-readRDS("data-raw/brain.RDS") #DESeq2 object
brain_metadata<-read.csv("data-raw/brain_metadata.csv", row.names = 1) #Metadata for DESeq2 objects

#usethis processed files
usethis::use_data(brain, overwrite = TRUE)
usethis::use_data(brain_metadata, overwrite = TRUE)
