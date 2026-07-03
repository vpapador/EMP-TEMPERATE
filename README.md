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




```
pie_phyla_temp <- trans_abund$new(dataset = mt_temperate, taxrank = "Phylum", ntaxa = 9, groupmean = "envo_biome_2")
pie_phyla_temp$plot_donut(label = FALSE)
pie_phyla_temp$plot_donut(label = TRUE)

```



<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/480e7c42-ea2f-430a-960a-fb359945e134" />




```

# merge samples as one community for each group
tmp <- mt_rarefied_temp$merge_samples("envo_biome_2")
# tmp is a new microtable object
# create trans_venn object
tt_ven <- trans_venn$new(tmp, ratio = NULL)
tt_ven$plot_venn()

```

<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/dc5f5a53-2c0f-41eb-a485-48efb9da36f8" />


### ALPHA DIVERSITY GRAPHS

```
library(ggplot2)
library(patchwork)


trans_alpha_temp <- trans_alpha$new(dataset = mt_rarefied_temp, group = "temperate")
# return t1$data_stat
head(trans_alpha_temp$data_stat)

trans_alpha_temp$cal_diff(method = "wilcox")
head(trans_alpha_temp$res_diff)



p1 <- trans_alpha_temp$plot_alpha(
  measure = "Chao1",
  shape = "temperate",
  add = "jitter",
  y_start = 0.1,
  y_increase = 0.1
) 
p2 <- trans_alpha_temp$plot_alpha(
  measure = "Shannon",
  shape = "temperate",
  add = "jitter",
  y_start = 0.1,
  y_increase = 0.1
) 


p3 <- trans_alpha_temp$plot_alpha(
  measure = "Simpson",
  shape = "temperate",
  add = "jitter",
  y_start = 0.1,
  y_increase = 0.1
)


p4 <- trans_alpha_temp$plot_alpha(
  measure = "PD",
  shape = "temperate",
  add = "jitter",
  y_start = 0.1,
  y_increase = 0.1,
  group_order = c("anthropogenic biome", "natural habitat")
)

(p1 + p2) /
(p3 + p4) 

```



<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/8d0e18fc-5311-4f6e-b824-a852e410bb76" />







```


trans_alpha_temp2 <- trans_alpha$new(dataset = mt_rarefied_temp, group = "envo_biome_2", by_group = "temperate")

trans_alpha_temp2$cal_diff(method = "wilcox")
library(magrittr)
trans_alpha_temp2$res_diff %<>% base::subset(Significance != "ns")
trans_alpha_temp2$plot_alpha(measure = "Chao1")
trans_alpha_temp2$res_diff %<>% base::subset(Significance != "ns")
trans_alpha_temp2$plot_alpha(measure = "Shannon")
trans_alpha_temp2$res_diff %<>% base::subset(Significance != "ns")
trans_alpha_temp2$plot_alpha(measure = "Simpson")
trans_alpha_temp2$res_diff %<>% base::subset(Significance != "ns")
trans_alpha_temp2$plot_alpha(measure = "PD")



library(patchwork)

# Κράτα μόνο τα σημαντικά αποτελέσματα
trans_alpha_temp2$res_diff <- subset(
  trans_alpha_temp2$res_diff,
  Significance != "ns"
)

p5 <- trans_alpha_temp2$plot_alpha(measure = "Chao1")
p6 <- trans_alpha_temp2$plot_alpha(measure = "Shannon")
p7 <- trans_alpha_temp2$plot_alpha(measure = "Simpson")
p8 <- trans_alpha_temp2$plot_alpha(measure = "PD")

(p5 + p6) /
(p7 + p8)



```



<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/a7298659-a5e2-4d89-a314-f6765ecea020" />




### BETA DIVERSITY GRAPHS
```

#TRANSBETA OBJECT

trans_beta_temp <- trans_beta$new(dataset = mt_rarefied_temp, group = "temperate", measure = "bray")
trans_beta_temp2 <- trans_beta$new(dataset = mt_rarefied_temp, group = "envo_biome_2", measure = "bray")


trans_beta_temp$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(trans_beta_temp$res_ordination)
# plot the PCoA result with confidence ellipse
trans_beta_temp$plot_ordination(plot_color = "temperate", plot_shape = "temperate", plot_type = c("point", "ellipse"))


```


<img width="882" height="783" alt="image" src="https://github.com/user-attachments/assets/bca04720-9bf5-4ddd-b630-9804931220a4" />


```

trans_beta_temp2$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(trans_beta_temp2$res_ordination)
# plot the PCoA result with confidence ellipse
trans_beta_temp2$plot_ordination(plot_color = "envo_biome_2", plot_shape = "envo_biome_2", plot_type = c("point", "ellipse"))


```


<img width="882" height="783" alt="image" src="https://github.com/user-attachments/assets/4e5a6047-5db7-40d2-b75d-c5ab34188b49" />



### PERMANOVA


```
#PERMANOVA
trans_beta_temp$cal_manova(manova_all = TRUE)
trans_beta_temp$res_manova


```




| Source | Df | Sum of Squares | R² | F | P-value | Significance |
|--------|---:|---------------:|---:|---:|:--------:|:------------:|
| Temperate | 1 | 13.459 | 0.131 | 40.979 | 0.001 | *** |
| Residual | 273 | 89.666 | 0.869 | – | – | – |
| Total | 274 | 103.126 | 1.000 | – | – | – |



### LEfSe


```

lda_temp <- trans_diff$new(dataset = mt_temperate, method = "lefse", group = "temperate", alpha = 0.01, lefse_subgroup = NULL)
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
lda_temp$plot_diff_bar(threshold = 4)
# we show 20 taxa with the highest LDA (log10)
lda_temp$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("anthropogenic biome", "natural habitat"))


# clade_label_level 5 represent phylum level in this analysis
# require ggtree package



```



<img width="672" height="671" alt="image" src="https://github.com/user-attachments/assets/df8c6918-494c-4cee-9c3b-d0623db3a6b0" />





```

mt_2 <- mt_temperate
tax <- mt_2$tax_table

tax[tax == "None"] <- NA

prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

for (i in seq_along(prefixes)) {
  tax[, i] <- ifelse(
    is.na(tax[, i]) | tax[, i] == "",
    paste0(prefixes[i], "unclassified"),
    tax[, i]
  )
}

mt_2$tax_table <- tax

mt_2$cal_abund()
mt_2 <- tidy_taxonomy(mt_2, na_fill = TRUE)
lda_temp3 <- trans_diff$new(dataset = mt_2, method = "lefse", group = "temperate", alpha = 0.01, lefse_subgroup = NULL)

lda_temp3$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 30, clade_label_level = 5, group_order = c("anthropogenic biome", "natural habitat"))




important <- c(
  # πάνω (Bacteroidetes / Verrucomicrobia περιοχή)
  "p__Verrucomicrobia", 
#"c__[Spartobacteria]",
  "p__Bacteroidetes", 
#"c__[Saprospirae]", "c__Cytophagia",
  #"f__Chitinophagaceae",

  # δεξιά πάνω (Nitrospirae / Archaea)
  "p__Nitrospirae", 
#"c__Nitrospira",
  #"k__Archaea",
 "p__Crenarchaeota",

  # δεξιά κάτω (Proteobacteria)
  "p__Proteobacteria", 
#"c__Alphaproteobacteria", "o__Rhizobiales",
 # "c__Betaproteobacteria", "c__Gammaproteobacteria",
  #"o__Xanthomonadales", "c__Deltaproteobacteria", "o__Syntrophobacterales",

  # κάτω αριστερά (Actinobacteria)
  "p__Actinobacteria",
# "c__Actinobacteria", "o__Actinomycetales",

  # αριστερά (Acidobacteria / Gemmatimonadetes)
  #"k__Bacteria",
 "p__Acidobacteria",
  "p__Gemmatimonadetes"
#, "o__Ellin6513"
)

lda_temp3$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 30, select_show_labels = important)

```



<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/f022a57f-b921-4b8b-9371-33a285228ad1" />


### NETWORK ANALYSIS

```

#network analysis
library(igraph)
anthropo_network <- trans_network$new(dataset = group_anthropo_temp, cor_method = "spearman", filter_thres = 0.001)
natural_network <- trans_network$new(dataset = group_natural_temp, cor_method = "spearman", filter_thres = 0.001)


# construct network; require igraph package
anthropo_network $cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
#anthropo_network $cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return anthropo_network $res_network

# construct network; require igraph package
natural_network $cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)
# use arbitrary coefficient threshold to contruct network
#natural_network$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
# return anthropo_network $res_network

# invoke igraph cluster_fast_greedy function for this undirected network 
anthropo_network$cal_module(method = "cluster_fast_greedy")

natural_network$cal_module(method = "cluster_fast_greedy")



V(anthropo_network$res_network)$module  # έλεγξε ότι δεν είναι NULL

# Εξαγωγή
write_graph(anthropo_network$res_network,
            file = "anthropo_network.graphml",
            format = "graphml")

vcount(anthropo_network$res_network)  # αριθμός κόμβων (OTUs)
ecount(anthropo_network$res_network)  



V(natural_network$res_network)$module  # έλεγξε ότι δεν είναι NULL

# Εξαγωγή
write_graph(natural_network$res_network,
            file = "natural_network.graphml",
            format = "graphml")

```





















<img width="1024" height="1024" alt="anthropo" src="https://github.com/user-attachments/assets/1022ceec-6c2d-4df9-b397-e1b0ef6f9719" />  <img width="1024" height="1024" alt="natural" src="https://github.com/user-attachments/assets/a1d6bb96-3b4e-4272-87dc-a98828b76e5a" />

















```
anthropo_network$plot_taxa_roles(use_type = 2)


```




<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/80c5fc67-01b4-4cc4-b745-ccc190be63be" />






```
natural_network$plot_taxa_roles(use_type = 2)


```



<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/9f14c7cd-5c37-48f2-8f6f-052a2c63f1b9" />






















### ΑΛΛΗΛΕΠΙΔΡΑΣΕΙΣ ΜΙΚΡΟΒΙΩΝ ΣΕ ΕΠΙΠΕΔΟ ΚΛΑΣΗΣ ΜΕΣΑ ΣΤΟ ΔΙΚΤΥΟ

Επιλέγω να κάνω τις θετικές, αντίστοιχα μπορώ να κάνω και τις αρνητικές αλληλεπιδράσεις

```
anthropo_network$cal_sum_links(taxa_level = "Class")
anthropo_network$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
anthropo_network$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))
```









<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/7c645b4b-640a-4a74-9c6e-6aac1e978609" />









```
natural_network$cal_sum_links(taxa_level = "Class")
natural_network$plot_sum_links(plot_pos = TRUE, plot_num = 10, color_values = RColorBrewer::brewer.pal(10, "Paired"))
natural_network$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))


```



<img width="1920" height="954" alt="image" src="https://github.com/user-attachments/assets/8fe2a10b-f6cd-43a7-9e2d-684856158dd1" />










### ΔΗΜΙΟΥΡΓΙΑ NULL MODELS ΜΕ ΤΟ DATASET ΤΩΝ TEMPERATE
```
# generate trans_nullmodel object
# as an example, we only use high abundance OTU with mean relative abundance > 0.0005
t1 <- trans_nullmodel$new(mt_temperate1, filter_thres = 0.0005)

# see null.model parameter for other null models
# null model run 500 times for the example
t1$cal_ses_betampd(runs = 500, abundance.weighted = TRUE)
# return t1$res_ses_betampd
```

```
mt_temperate1$beta_diversity[["betaNRI"]] <- t1$res_ses_betampd
# create trans_beta class, use measure "betaNRI"
t2 <- trans_beta$new(dataset = mt_temperate1, group = "temperate", measure = "betaNRI")
# transform the distance for each group
t2$cal_group_distance()
# see the help document for more methods, e.g. "anova" and "KW_dunn"
t2$cal_group_distance_diff(method = "wilcox")
# plot the results
g1 <- t2$plot_group_distance(add = "mean")
g1 + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)

```




<img width="672" height="671" alt="image" src="https://github.com/user-attachments/assets/4bc84134-9c1f-46d2-a75e-5677e7aa83e3" />









```

tmp2 <- "./test1"; dir.create(tmp2)
t1$cal_ses_betamntd(runs = 1000, abundance.weighted = TRUE, use_iCAMP = TRUE, iCAMP_tempdir = tmp2)


# result stored in t1$res_rcbray
t1$cal_rcbray(runs = 1000)
# return t1$res_rcbray
object$res_rcbray


# use betaNTI and rcbray to evaluate processes
t1$cal_process(use_betamntd = TRUE, group = "temperate")
t1$cal_process(use_betamntd = TRUE)
t1$res_process



```




| Community Type          | Process                 | Percentage (%) |
|------------------------|------------------------|----------------|
| Natural habitat        | Variable selection     | 34.75          |
| Natural habitat        | Homogeneous selection  | 1.59           |
| Natural habitat        | Dispersal limitation   | 8.90           |
| Natural habitat        | Homogeneous dispersal  | 1.24           |
| Natural habitat        | Drift                  | 53.52          |
| Anthropogenic biome    | Variable selection     | 7.63           |
| Anthropogenic biome    | Homogeneous selection  | 0.28           |
| Anthropogenic biome    | Dispersal limitation   | 70.94          |
| Anthropogenic biome    | Homogeneous dispersal  | 0.00           |
| Anthropogenic biome    | Drift                  | 21.15          |



```
# require NST package to be installed
t1$cal_NST(method = "tNST", group = "temperate", dist.method = "bray", abundance.weighted = TRUE, output.rand = TRUE, SES = TRUE)
t1$res_NST$index.grp

```





| Group                | Size | ST.i.bray | NST.i.bray | MST.i.bray | SES.i.bray |
|---------------------|------|-----------|------------|------------|------------|
| Natural habitat     | 8385 | 0.776     | 0.517      | 0.449      | 3.756      |
| Anthropogenic biome | 10440 | 0.724     | 0.603      | 0.496      | -0.809     |


```
#αυτο δε το ετρεξα 
t1$cal_NST_test(method = "nst.boot")


# for pNST method, phylogenetic tree is needed
t1$cal_NST(method = "pNST", group = "temperate", output.rand = TRUE, SES = TRUE)
t1$cal_NST_test(method = "nst.boot")

#αυτο δε το ετρεξα
t1$cal_NRI(null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
t1$cal_NTI(null.model = "taxa.labels", abundance.weighted = TRUE, runs = 999)

```

