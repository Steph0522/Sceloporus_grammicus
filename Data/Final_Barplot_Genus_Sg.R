# Get work directory 
setwd("~/Documents/MICROBIOTA REPTILES/PAPERS 2022/Second Article/Integrative_Zoology")
library(tidyverse)

# Load files
otutable <- read.csv("feature_table.csv", row.names = 1)
metadata <- read.csv("metadata.csv", check.names = F)
metadata$Ind <- as.factor(metadata$Ind)
metadata$Library <- as.factor(metadata$Library)
metadata$SampleType <- as.factor(metadata$SampleType)
taxonomy <- read.csv("taxonomy.csv", check.names = F) %>% unite(
  taxa, Kingdom:Species, remove = F, sep = ";")

otutable_metadata <- otutable %>% rownames_to_column(var="OTUID") %>% 
  inner_join(taxonomy)

################################################################################
# Calcular la abundancia relativa por género y SampleID
Genus_01 <- otutable_metadata %>% group_by(Genus) %>% summarise_if(is.numeric, sum)
Genus_01 <- Genus_01[c(-1:-2),]

print(Genus_01)
na.omit(Genus_01)

Genus_01 <- Genus_01 %>% column_to_rownames(var = "Genus")
Genus.ra <- t(t(Genus_01)/colSums(Genus_01)*100)

ver <- rowMeans(Genus.ra) %>% as.data.frame() %>% arrange(desc(.))
apply(Genus.ra, 1, sd)

Genus <- Genus.ra %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata)

Prom_Genus <- Genus %>% group_by(SampleType) %>% summarise_if(is.numeric, mean)
SD_Genus <- Genus %>% group_by(SampleType) %>% summarise_if(is.numeric, sd)

aggregate(Genus[ ,2:141], list(Genus$SampleType), mean)
aggregate(Genus[ ,2:141], list(Genus$SampleType), sd)
#write.table(Prom_Genus, file="./Genus_Pro_ST.txt", sep = "\t")

################################################################################
# Lo primero que vamos a realizar es cargar los datos que corresponden a:
# 1. Metadatos
# 2. OTU_table
# 3. Taxonomía
# 4. Phylo tree ()

library(phyloseq)
library(ggplot2)
library(vegan)
library(picante)
library(ape)
library(devtools)
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

# Importar los datos de las tablas y generar objetos o categorías phyloseq
# Las tablas contienten información de la abundancia de taxones bacterianos (OTU), asignación 
# de taxones (TAX) y datos de muestra (SAM). 

# Nota: merge_phyloseq combina la información de todas las tablas.

metadata <- read.csv(file = "metadata.csv", header = TRUE, row.names = 1) 
otu_table <- read.csv(file = "feature_table.csv", check.names = F) 
#taxonomy_raw<- read.csv(file = "Genus_Abun_Rel_Sg.csv", check.names = F)
taxonomy <-  read.csv("taxonomy.csv", check.names = F) %>% mutate_at(
  c("Genus"), str_replace,"g__", "")

lista <- rowMeans(Genus.ra) %>% as.data.frame() %>% arrange(desc(.)) %>% 
  slice_head(n=13) %>% rownames_to_column(var = "Genus") %>% 
  filter(!Genus =="g__") %>% filter(!Genus == "Unassigned") %>% 
  mutate_at(c("Genus"), str_replace,"g__", "")
list <- lista$Genus
lista01 <- read.csv(file = "lista.csv", check.names = F)
list02 <- lista01$Genus
#write.table(lista, file="./lista.txt", sep = "\t")

taxonomy_filter <- taxonomy %>% filter(Genus %in% list02)
taxonomy %>% inner_join(lista)
taxonomy_1 <- taxonomy_filter %>% inner_join(otu_table, by =c(
  "OTUID"="OTUID")) %>% dplyr::select(1:8)

otu_table_1 <- read.csv(file = "feature_table.csv", header = TRUE,
                        row.names = 1) %>% rownames_to_column(var = "OTUID") %>% 
  inner_join(taxonomy_1, by = "OTUID") %>% dplyr::select(-52:-58) %>% 
  column_to_rownames(var = "OTUID")

taxo <- taxonomy_1 %>% column_to_rownames(var = "OTUID")

# Crear objeto de categoría phyloseq
SAM <- sample_data(metadata)
TAX <- tax_table(as.matrix(taxo)) 
OTU <- otu_table(otu_table_1, taxa_are_rows=TRUE)  
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
physeq_sub <- subset_taxa(physeq, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned"))
physeq_sub <- subset_taxa(physeq, !is.na(Genus) & !Genus %in% c("", "Unassigned"))

################################################################################
paleta <- c("#6699CC", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
            "#44AA99", "#999933", "#888888", "#661100", "#88CCEE", "#E69F00",
            "#004949")
print(paleta)

Final_Genus_Sg <- plot_bar(physeq = relative, "Sample", fill = "Genus") +
  facet_grid(~factor(SampleType, levels = c("Stomach", "Small intestine", "Rectum", 
                                            "Cloaca", "Feces"),
                     labels = c("Stomach", "Small intestine", "Rectum", "Cloaca", 
                                "Feces")), scales = "free", space = "free") +
  labs(y="Relative abundance") +
  geom_bar(stat = "identity", position = "stack", res=300) +
  scale_fill_manual(values = paleta) + 
  theme(legend.text = element_text(face = "italic")) +
  scale_fill_manual(values = paleta) + 
  theme(strip.text.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  theme(text = element_text(size = 10))

print(Final_Genus_Sg)
ggsave("Final_Genus_Sg.jpeg", width=7.2, height=4.5, dpi=300)