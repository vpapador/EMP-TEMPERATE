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

```
### ΔΗΜΙΟΥΡΓΊΑ ΔΟΜΉΣ ΓΙΑ MICROECO OBJECT ΓΙΑ ΤΟ TEMPERATE

```
library(microeco)
#community matrix
otu_table_temperate<-dataframe_FINAL_temp
View(otu_table_temperate)
rownames(otu_table_temperate) <- otu_table_temperate[,1]
otu_table_temperate <- otu_table_temperate[,-1]


#metadata
sample_table_temperate<-temperate_metadata_final
rownames(sample_table_temperate) <- sample_table_temperate[,1]
View(sample_table_temperate)

#taxonomy
tax_table_temperate<-asv_taxa_temperate_subset_FINAL
rownames(tax_table_temperate) <- tax_table_temperate[,1]
tax_table_temperate <- tax_table_temperate[,-1]
View(tax_table_temperate)

#phylogeny
tree_temperate

```
### ΔΗΜΙΟΥΡΓΊΑ MICROECO OBJECT ΓΙΑ ΤΟ TEMPERATE ΚΑΙ RAREFACTION
```
#Δημιουργώ ένα microeco object 
mt_temperate <- microtable$new(otu_table = otu_table_temperate, sample_table = sample_table_temperate, tax_table = tax_table_temperate, phylo_tree = tree_temperate)
mt_temperate1<-mt_temperate
mt_temperate1$tidy_dataset()
mt_temperate1
head(mt_temperate1$tax_table)


# default parameter (rel = TRUE) denotes relative abundance
mt_temperate1$cal_abund()

# first clone the data
mt_rarefied_temp <- clone(mt_temperate1)
# use sample_sums to check the sequence numbers in each sample
mt_rarefied_temp$sample_sums() %>% range

# As an example, use 10000 sequences in each sample
mt_rarefied_temp$rarefy_samples(sample.size = 6500)

```
### ΔΗΜΙΟΥΡΓΙΑ SUBSETS ΓΙΑ ΤΑ TEMPERATE ANTHROPOGENIC & NATURAL
```
#antropo
group_anthropo_temp <- clone(mt_temperate)
# select 'anthropogenic biome'
group_anthropo_temp$sample_table <- subset(group_anthropo_temp$sample_table, temperate == "anthropogenic biome")
# trim all the data
group_anthropo_temp$tidy_dataset()

group_anthropo_temp$cal_abund()

# first clone the data
group_anthropo_temp_rarefied <- clone(group_anthropo_temp)
# use sample_sums to check the sequence numbers in each sample
group_anthropo_temp_rarefied$sample_sums() %>% range

# As an example, use 10000 sequences in each sample
group_anthropo_temp_rarefied$rarefy_samples(sample.size = 6500)

#natural 
group_natural_temp <- clone(mt_temperate)
# select 'anthropogenic biome'
group_natural_temp$sample_table <- subset(group_natural_temp$sample_table, temperate == "natural habitat")
# trim all the data
group_natural_temp$tidy_dataset()

group_natural_temp$cal_abund()

# first clone the data
group_natural_temp_rarefied <- clone(group_natural_temp)
# use sample_sums to check the sequence numbers in each sample
group_natural_temp_rarefied$sample_sums() %>% range

# As an example, use 10000 sequences in each sample
group_natural_temp_rarefied$rarefy_samples(sample.size = 6500)

```


### ALPHA DIVERSITY ΓΙΑ ΤΟ TEMPERATE
```
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
mt_rarefied_temp$cal_alphadiv(PD = TRUE)

# return alpha_diversity in the object
class(mt_rarefied_temp$alpha_diversity)

```
### BETA DIVERSITY ΓΙΑ ΤΟ TEMPERATE
```
# require GUniFrac package installed
mt_rarefied_temp$cal_betadiv(unifrac = TRUE)
# return beta_diversity list in the object
class(mt_rarefied_temp$beta_diversity)

```

### RELATIVE ABUNDANCE
```
phyla_temp <- trans_abund$new(dataset = mt_temperate, taxrank = "Phylum", ntaxa = 8)
phyla_temp$plot_bar(others_color = "grey70", facet = "temperate", xtext_keep = FALSE, legend_text_italic = FALSE)
# return a ggplot2 object


# require package ggh4x, first run install.packages("ggh4x") if not installed
phyla_temp$plot_bar(others_color = "grey70", facet = c("temperate", "envo_biome_2"), xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)
```
<img width="1746" height="868" alt="image" src="https://github.com/user-attachments/assets/47de26ac-a259-4b32-942c-2d383dffe610" />

```
pie_phyla_temp <- trans_abund$new(dataset = mt_temperate, taxrank = "Phylum", ntaxa = 9, groupmean = "temperate")
pie_phyla_temp$plot_donut(label = FALSE)
pie_phyla_temp$plot_donut(label = TRUE)

```
<img width="1006" height="688" alt="image" src="https://github.com/user-attachments/assets/3a4512c6-0684-43f0-95d4-ee997d2d7160" />


