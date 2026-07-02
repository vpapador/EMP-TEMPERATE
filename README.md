# FINAL-FINAL-FINAL

### ΌΛΑ ΤΑ ΔΕΔΟΜΈΝΑ ΜΟΥ
```
#COMMUNITY MATRIX ALL

#SUBSET TEMPERATE

#SUBSET 1500
FINAL_temp<-read.table(
  file.choose(),
  header = TRUE,
  sep = ";",
  check.names = FALSE
)
rownames(FINAL_temp) <- FINAL_temp[,1]
FINAL_temp<- FINAL_temp[, -1]
dataframe_FINAL_temp<-FINAL_temp
dataframe_FINAL_temp <- rownames_to_column(dataframe_FINAL_temp, var = "ASV")
dim(dataframe_FINAL_temp)


#SUBSET 500
library(tibble)
library(tidyverse)
FINAL_500_all<-read.table(
  file.choose(),
  header = TRUE,
  sep = ";",
  check.names = FALSE
)
rownames(FINAL_500_all) <- FINAL_500_all[,1]
FINAL_500_all<- FINAL_500_all[, -1]
dataframe_FINAL_500_all<-FINAL_500_all
dataframe_FINAL_500_all <- rownames_to_column(dataframe_FINAL_500_all, var = "ASV")
dim(dataframe_FINAL_500_all)

#SUBSET 1000
FINAL_1000<-read.table(
  file.choose(),
  header = TRUE,
  sep = ";",
  check.names = FALSE
)
rownames(FINAL_1000) <- FINAL_1000[,1]
FINAL_1000<- FINAL_1000[, -1]
dataframe_FINAL_1000<-FINAL_1000
dataframe_FINAL_1000 <- rownames_to_column(dataframe_FINAL_1000, var = "ASV")
dim(dataframe_FINAL_1000)

#SUBSET 1500
FINAL_1500<-read.table(
  file.choose(),
  header = TRUE,
  sep = ";",
  check.names = FALSE
)
rownames(FINAL_1500) <- FINAL_1500[,1]
FINAL_1500<- FINAL_1500[, -1]
dataframe_FINAL_1500<-FINAL_1500
dataframe_FINAL_1500 <- rownames_to_column(dataframe_FINAL_1500, var = "ASV")
dim(dataframe_FINAL_1500)

################################################################

#METADATA

#περνάω μέσα το σωστό temperate
temperate_metadata_final<-read.csv(file.choose(), sep=";" , header=TRUE, dec=".")
View(temperate_metadata_final)


#περνάω τα 500
metadata_500_final<-read.csv(file.choose(), sep="," , header=TRUE, dec=".")
View(metadata_500_final)

#περνάω τα 1000
metadata_1000_final<-read.csv(file.choose(), sep="," , header=TRUE, dec=".")
View(metadata_1000_final)


#περνάω τα 1500
metadata_1500_final<-read.csv(file.choose(), sep="," , header=TRUE, dec=".")
View(metadata_1500_final)

#############################################################

#Phylogeny

library(ape)
tree <- read.tree(file.choose())
tree_temperate <- read.tree(file.choose())
tree_500 <- read.tree(file.choose())
tree_1000 <- read.tree(file.choose())
tree_1500 <- read.tree(file.choose())
length(tree_temperate$tip.label)


########################################################

#Taxonomy


asv_taxa_temperate_subset_FINAL<- read.table(
  file.choose(),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
asv_taxa_temperate_subset_FINAL<-asv_taxa_temperate_subset_FINAL[,-1]
View(asv_taxa_temperate_subset_FINAL)



asv_taxa_subset_500_all<- read.table(
  file.choose(),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
asv_taxa_subset_500_all<-asv_taxa_subset_500_all[,-1]
View(asv_taxa_subset_500_all)


asv_taxa_subset_1000_all<- read.table(
  file.choose(),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
asv_taxa_subset_1000_all<-asv_taxa_subset_1000_all[,-1]
View(asv_taxa_subset_1000_all)



asv_taxa_subset_1500_all<- read.table(
  file.choose(),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
asv_taxa_subset_1500_all<-asv_taxa_subset_1500_all[,-1]
View(asv_taxa_subset_1500_all)


#####################################################

```
