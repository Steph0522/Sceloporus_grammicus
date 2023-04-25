setwd("~/Documents/MICROBIOTA REPTILES/PAPERS 2022/Second Article/Integrative_Zoology")
library(phyloseq)
library(ggplot2)
library(vegan)
library(picante)
library(ape)
library(devtools)
library("DESeq2")
library(vegan)
library("edgeR")
library("RColorBrewer")
library(scales)
library(grid)
library(reshape2)
library(dplyr)
library(scales)
library(viridis)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(dplyr)

# Nota: merge_phyloseq combina la información de todas las tablas.

metadata <- read.csv(file = "metadata.csv", header = TRUE, row.names = 1)
otu_table <- read.csv("feature_table.csv", header = TRUE, row.names = 1)
taxonomy <- read.csv("taxonomy.csv", header = TRUE, row.names = 1) %>% 
  mutate_at(c("Phylum"), str_replace,"p__", "")
#phylo<-read.tree(file = "tree.nwk") 

# Crear objeto de categoría phyloseq
SAM <- sample_data(metadata)
TAX <- tax_table(as.matrix(taxonomy))
OTU <- otu_table(otu_table, taxa_are_rows = TRUE)  
#PHY<-phy_tree(phylo)
physeq <- merge_phyloseq(OTU, TAX, SAM)

# Visualización de los datos
sample_names(physeq)
rank_names(physeq)
sample_variables(physeq)

# Convertir los recuentos de OTUs en abundancias relativas y normalizar  
# la abundacia de cada OTU
relative  = transform_sample_counts(physeq = physeq, function(OTU) OTU / sum(OTU))

# Remover las bacterias sin identificacion a nivel Kingdom
physeq_sub1 <- subset_taxa(physeq, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned"))
physeq_sub2 <- subset_taxa(physeq, !is.na(Genus) & !Genus %in% c("", "Unassigned"))

################################################################################
# Abundancia relativa por muestra a nivel de phylum (Original paper)
# Genera paleta de colores
library(rcartocolor)
my_colors = carto_pal(12, "Safe")
my_colors
blind_pal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
               "#44AA99", "#999933", "#888888", "#661100", "#6699CC", "#E69F00")

Phylum_RelAbund <- plot_bar(physeq = relative, "Sample", fill = "Phylum") + 
  facet_grid(~factor(SampleType, levels = c("Stomach", "Small intestine", "Rectum", "Cloaca","Feces"),
                                            labels= c("Stomach", "Small intestine", "Rectum", "Cloaca", "Feces")), 
             scales = "free", space = "free") +
  labs(y="Relative abundance") +
  geom_bar(stat = "identity", position="stack", res=300) +
  scale_fill_manual(values = blind_pal)+ theme(strip.text.x = element_text(face = "bold"),
                                           axis.title.y = element_text(face = "bold")) +
  theme(text = element_text(size = 10))  

print(Phylum_RelAbund)
ggsave("Phylum_RelAbund.jpeg", width=7.2, height=4.5, dpi=300)

################################################################################
# Calcular los porcentajes de abundacia relativa por phylum y por muestra

# Load files
otutable <- read.csv("feature_table.csv", row.names = 1)
metadata <- read.csv("metadata.csv", check.names = F)
metadata$Ind <- as.factor(metadata$Ind)
metadata$Library <- as.factor(metadata$Library)
metadata$SampleType <- as.factor(metadata$SampleType)
taxonomy <- read.csv("taxonomy.csv", check.names = F) %>% unite(
  taxa, Kingdom:Species, remove = F, sep = ";")

otutable_metadata <- otu_table %>% rownames_to_column(var="OTUID") %>% 
  inner_join(taxonomy)

phyl_01 <- otutable_metadata %>% group_by(Phylum) %>% 
  summarise_if(is.numeric, sum)

print(phyl_01)
na.omit(phyl_01)

phyl_01 <- phyl_01 %>% column_to_rownames(var = "Phylum")
phy.ra <- t(t(phyl_01)/colSums(phyl_01)*100)

rowMeans(phy.ra) %>% as.data.frame() %>% arrange(desc(.))
apply(phy.ra,1,sd)

phylum <- phy.ra  %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata)

prom <- phylum %>% group_by(SampleType) %>% summarise_if(is.numeric, mean)
sd <- phylum %>% group_by(SampleType) %>% summarise_if(is.numeric, sd)

aggregate(phylum[ ,2:13], list(phylum$SampleType), mean)
aggregate(phylum[ ,2:13], list(phylum$SampleType), sd)
