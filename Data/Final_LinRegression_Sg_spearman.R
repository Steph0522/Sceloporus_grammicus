# Get work directory 
setwd("~/Documents/MICROBIOTA REPTILES/PAPERS 2022/Second Article/Integrative_Zoology")

# Linear regression
# Laod packages
library(tidyverse)
library(CoDaSeq)
library(zCompositions)
library(compositions)
library(propr)
library(CoDaSeq)

#Load files
phyl <- read_csv("level-2.csv")
phyl2 <- phyl %>% dplyr::select(index, contains("d__")) %>% 
  column_to_rownames(var = "index")

d.pro <- cmultRepl(t(phyl2), method = "CZM", output = "p-counts")
d.clr.abund.codaseq <- codaSeq.clr(x= d.pro, samples.by.row = F)

############################
### Cloaca versus Rectum ###
############################
phyl_S_R <- data.frame(t(d.clr.abund.codaseq)) %>% 
  rownames_to_column(var = "index") %>% inner_join(phyl) %>% 
  dplyr::select(c(1:13), SampleType) %>% filter(
      SampleType=="Rectum"|SampleType=="Swab") %>% 
  #dplyr::select(contains(c("HC", "R"))) %>% 
  pivot_longer(cols = starts_with("d__"),
               names_to = "names", values_to = "values") %>% pivot_wider(
                 names_from = SampleType, values_from = values) %>% 
  replace(is.na(.), 0)

otu_S_R <- phyl_S_R %>% dplyr::select(-index)
namesotu <- otu_S_R$names
write_tsv(phyl_S_R, "phyl_S_R.tsv")

# Create pallete color 
blind_pal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
               "#44AA99", "#999933", "#888888", "#661100", "#6699CC", "#E69F00")
#Load file
SR <- read.csv("Swab_Rectum.csv")
Cloaca_Rectum <- SR %>% ggplot(aes(x=Rectum, y=Cloaca, color=Phylum)) + 
  geom_point() +
  scale_color_manual(values = blind_pal) +
  #stat_summary(fun.data= mean_cl_normal) + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x=7, y=-5, size=5,label=bquote(paste('r'['s']*'=',.(round(
    cor(SR$Cloaca, SR$Rectum, method = "spearman"),digits = 2)))))
#  labs(title = paste("Adj R2 = ",signif(summary(data.lm_SR)$adj.r.squared, 5),
 #                   "Intercept =",signif(data.lm_SR$coef[[1]],5 ),
  #                   " Slope =",signif(data.lm_SR$coef[[2]], 5),
   #                  " P =",signif(summary(data.lm_SR)$coef[2,4], 5)))
print(Cloaca_Rectum)
ggsave("Cloaca_Rectum.jpeg", width=7, height=4.5, dpi=300)
cor1 <- cor.test(SR$Cloaca, SR$Rectum, method = "spearman")
cor1

#####################################
### Cloaca versus Small intestine ###
#####################################
phyl_S_SI <- data.frame(t(d.clr.abund.codaseq)) %>% 
  rownames_to_column(var = "index") %>% 
  inner_join(phyl) %>% dplyr::select(c(1:13), SampleType) %>% filter(
      SampleType=="Smallintestine"|SampleType=="Swab") %>% 
  #dplyr::select(contains(c("HC", "R"))) %>% 
  pivot_longer(cols = starts_with("d__"),
               names_to = "names", values_to = "values") %>% pivot_wider(
                 names_from = SampleType, values_from = values) %>% 
  replace(is.na(.), 0)

otu_S_SI <- phyl_S_SI %>% dplyr::select(-index)
namesotu <- otu_S_SI$names
write_tsv(phyl_S_SI, "phyl_S_SI.tsv")

#Load file
SI <- read.csv("Swab_Intestine.csv")
Cloaca_Intest <- SI %>% ggplot(aes(x=Small.intestine, y=Cloaca, color=Phylum)) + 
  geom_point() +
  scale_color_manual(values = blind_pal) +
  #stat_summary(fun.data= mean_cl_normal) + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x=7, y=-5, size=5,label=bquote(paste('r'['s']*'=',.(round(
    cor(SI$Cloaca, SI$Small.intestine, method = "spearman"),digits = 2)))))
# labs(title = paste("Adj R2 = ",signif(summary(data.lm_SR)$adj.r.squared, 5),
#                   "Intercept =",signif(data.lm_SR$coef[[1]],5 ),
#                  " Slope =",signif(data.lm_SR$coef[[2]], 5),
#                 " P =",signif(summary(data.lm_SR)$coef[2,4], 5)))
print(Cloaca_Intest)
ggsave("Cloaca_Intest.jpeg", width=7, height=4.5, dpi=300)
cor2 <- cor.test(SI$Cloaca, SI$Small.intestine, method = "spearman")
cor2

#############################
### Cloaca versus Stomach ###
#############################
phyl_S_Sto <- data.frame(t(d.clr.abund.codaseq)) %>% 
  rownames_to_column(var = "index") %>% 
  inner_join(phyl) %>% dplyr::select(c(1:13), SampleType) %>% filter(
    SampleType=="Stomach"|SampleType=="Swab") %>% 
  #dplyr::select(contains(c("HC", "R"))) %>% 
  pivot_longer(cols = starts_with("d__"),
               names_to = "names", values_to = "values") %>% pivot_wider(
                 names_from = SampleType, values_from = values) %>% 
  replace(is.na(.), 0)

otu_S_Sto <- phyl_S_Sto %>% dplyr::select(-index)
namesotu <- otu_S_Sto$names
write_tsv(phyl_S_Sto, "phyl_S_Sto.tsv")

#Load file
SStom <- read.csv("Swab_Stomach.csv")
Cloaca_Stomach <- SStom %>% ggplot(aes(x=Stomach, y=Cloaca, color=Phylum)) + 
  geom_point() +
  scale_color_manual(values = blind_pal) +
  #stat_summary(fun.data= mean_cl_normal) + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x=7, y=-5, size=5,label=bquote(paste('r'['s']*'=',.(round(
    cor(SStom$Cloaca, SStom$Stomach, method = "spearman"),digits = 2)))))
# labs(title = paste("Adj R2 = ",signif(summary(data.lm_SR)$adj.r.squared, 5),
#                   "Intercept =",signif(data.lm_SR$coef[[1]],5 ),
#                  " Slope =",signif(data.lm_SR$coef[[2]], 5),
#                 " P =",signif(summary(data.lm_SR)$coef[2,4], 5)))
print(Cloaca_Stomach)
ggsave("Cloaca_Stomach.jpeg", width=7, height=4.5, dpi=300)
cor3 <- cor.test(SStom$Cloaca, SStom$Stomach, method = "spearman")
cor3

###########################
### Feces versus Rectum ###
###########################
phyl_F_R <- data.frame(t(d.clr.abund.codaseq)) %>% 
  rownames_to_column(var = "index") %>% 
  inner_join(phyl) %>% dplyr::select(c(1:13), SampleType) %>% filter(
    SampleType=="Rectum"|SampleType=="Feces") %>% 
  #dplyr::select(contains(c("HC", "R"))) %>% 
  pivot_longer(cols = starts_with("d__"),
               names_to = "names", values_to = "values") %>% pivot_wider(
                 names_from = SampleType, values_from = values) %>% 
  replace(is.na(.), 0)

otu_F_R <- phyl_F_R %>% dplyr::select(-index)
namesotu <- otu_F_R$names
write_tsv(phyl_F_R, "phyl_F_R.tsv")

# Load file
FR <- read.csv("Feces_Rectum.csv")
Feces_Rectum <- FR %>% ggplot(aes(x=Rectum, y=Feces, color=Phylum)) +
  geom_point() +
  scale_color_manual(values = blind_pal) +
  #stat_summary(fun.data= mean_cl_normal) + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x=7, y=-5, size=5,label=bquote(paste('r'['s']*'=',.(round(
    cor(FR$Feces, FR$Rectum, method = "spearman"),digits = 2)))))
# labs(title = paste("Adj R2 = ",signif(summary(data.lm_SR)$adj.r.squared, 5),
#                   "Intercept =",signif(data.lm_SR$coef[[1]],5 ),
#                  " Slope =",signif(data.lm_SR$coef[[2]], 5),
#                 " P =",signif(summary(data.lm_SR)$coef[2,4], 5)))
print(Feces_Rectum)
ggsave("Feces_Rectum.jpeg", width=7, height=4.5, dpi=300)
cor4 <- cor.test(FR$Feces, FR$Rectum, method = "spearman")
cor4

####################################
### Feces versus Small intestine ###
####################################
phyl_F_SI <- data.frame(t(d.clr.abund.codaseq)) %>% 
  rownames_to_column(var = "index") %>% 
  inner_join(phyl) %>% dplyr::select(c(1:13), SampleType) %>% filter(
    SampleType=="Smallintestine"|SampleType=="Feces") %>% 
  #dplyr::select(contains(c("HC", "R"))) %>% 
  pivot_longer(cols = starts_with("d__"),
               names_to = "names", values_to = "values") %>% pivot_wider(
                 names_from = SampleType, values_from = values) %>% 
  replace(is.na(.), 0)

otu_F_SI <- phyl_F_SI %>% dplyr::select(-index)
namesotu <- otu_F_SI$names
write_tsv(phyl_F_SI, "phyl_F_SI.tsv")

# Load file
FI <- read.csv("Feces_Intestine.csv")
Feces_Intest <- FI %>% ggplot(aes(x=Small.intestine, y=Feces, color=Phylum))+
  geom_point() +
  scale_color_manual(values = blind_pal) +
  #stat_summary(fun.data= mean_cl_normal) + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x=7, y=-5, size=5,label=bquote(paste('r'['s']*'=',.(round(
    cor(FI$Feces, FI$Small.intestine, method = "spearman"),digits = 2)))))
# labs(title = paste("Adj R2 = ",signif(summary(data.lm_SR)$adj.r.squared, 5),
#                   "Intercept =",signif(data.lm_SR$coef[[1]],5 ),
#                  " Slope =",signif(data.lm_SR$coef[[2]], 5),
#                 " P =",signif(summary(data.lm_SR)$coef[2,4], 5)))
print(Feces_Intest)
ggsave("Feces_Intest.jpeg", width=7, height=4.5, dpi=300)
cor5 <- cor.test(FI$Feces, FI$Small.intestine, method = "spearman")
cor5

############################
### Feces versus Stomach ###
############################
phyl_F_Sto <- data.frame(t(d.clr.abund.codaseq)) %>% 
  rownames_to_column(var = "index") %>% 
  inner_join(phyl) %>% dplyr::select(c(1:11), SampleType) %>% filter(
    SampleType=="Stomach"|SampleType=="Feces") %>% 
  #dplyr::select(contains(c("HC", "R"))) %>% 
  pivot_longer(cols = starts_with("d__"),
               names_to = "names", values_to = "values") %>% pivot_wider(
                 names_from = SampleType, values_from = values) %>% 
  replace(is.na(.), 0)

otu_F_Sto <- phyl_F_Sto %>% dplyr::select(-index)
namesotu <- otu_F_Sto$names
write_tsv(phyl_F_Sto, "phyl_F_Sto.tsv")

# Load file
FStom <- read.csv("Feces_Stomach.csv")
Feces_Stomach <- FStom %>% ggplot(aes(x=Stomach, y=Feces, color=Phylum)) + 
  geom_point() +
  scale_color_manual(values = blind_pal) +
  #stat_summary(fun.data= mean_cl_normal) + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x=7, y=-5, size=5,label=bquote(paste('r'['s']*'=',.(round(
    cor(FStom$Feces, FStom$Stomach, method = "spearman"),digits = 2)))))
# labs(title = paste("Adj R2 = ",signif(summary(data.lm_SR)$adj.r.squared, 5),
#                   "Intercept =",signif(data.lm_SR$coef[[1]],5 ),
#                  " Slope =",signif(data.lm_SR$coef[[2]], 5),
#                 " P =",signif(summary(data.lm_SR)$coef[2,4], 5)))
print(Feces_Stomach)
ggsave("Feces_Stomach.jpeg", width=7, height=4.5, dpi=300)
cor6 <- cor.test(FStom$Feces, FStom$Stomach, method = "spearman")
cor6
