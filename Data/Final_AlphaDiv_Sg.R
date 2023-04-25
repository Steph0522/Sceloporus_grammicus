setwd("~/Documents/MICROBIOTA REPTILES/PAPERS 2022/Second Article/Integrative_Zoology")

# Cargar la base de datos
Hill <- read.csv(file = "feature_table.csv", header = TRUE, row.names = 1)

# Calcular la diversidad a los diferentes órdenes
library(hillR) 
q0 <- hill_taxa(comm = t(Hill), q = 0)
q1 <- hill_taxa(comm = t(Hill), q = 1)
q2 <- hill_taxa(comm = t(Hill), q = 2)
Hill_div <- cbind(q0, q1, q2)
write.table(Hill_div, file="./HillTaxa_Div.txt", sep = "\t")

## ALPHA DIVERSITY
library(tidyverse)
library(ggpubr)

#loading files
alpha <- read.csv("Hill_Numbers_q012.csv") %>% dplyr::select(SampleID, q0, q1, q2)
metadata <- read.csv("metadata.csv", check.names = F)
alpha <- alpha %>% inner_join(metadata, by = c("SampleID"="SampleID"))

# Normality test
shapiro.test(x =alpha$q0)
shapiro.test(x =alpha$q1)
shapiro.test(x =alpha$q2)

hist(alpha$q0)
hist(alpha$q1)
hist(alpha$q2)
# Data are not normal

# Paired test (Wilcoxon) (q0)
feces <- subset(alpha, SampleType== "Feces", q0, drop = TRUE)
cloaca <- subset(alpha, SampleType== "Cloaca", q0, drop = TRUE)
stomach <-  subset(alpha, SampleType== "Stomach", q0, drop = TRUE)
intestine <- subset(alpha, SampleType== "Small intestine", q0, drop = TRUE)
rectum <-  subset(alpha, SampleType== "Rectum", q0, drop = TRUE)

FvsS_q0 <- wilcox.test(x= feces, y= stomach, paired = TRUE)
FvsI_q0 <- wilcox.test(x= feces, y= intestine, paired = TRUE)
FvsR_q0 <- wilcox.test(x= feces, y= rectum, paired = TRUE)
CvsS_q0 <- wilcox.test(x= cloaca, y= stomach, paired = TRUE)
CvsI_q0 <- wilcox.test(x= cloaca, y= intestine, paired = TRUE)
CvsR_q0 <- wilcox.test(x= cloaca, y= rectum, paired = TRUE)

# Paired test (Wilcoxon) (q1)
feces1 <- subset(alpha, SampleType== "Feces", q1, drop = TRUE)
cloaca1 <- subset(alpha, SampleType== "Cloaca", q1, drop = TRUE)
stomach1 <-  subset(alpha, SampleType== "Stomach", q1, drop = TRUE)
intestine1 <- subset(alpha, SampleType== "Small intestine", q1, drop = TRUE)
rectum1 <-  subset(alpha, SampleType== "Rectum", q1, drop = TRUE)

FvsS_q1 <- wilcox.test(x= feces1, y= stomach1, paired = TRUE)
FvsI_q1 <- wilcox.test(x= feces1, y= intestine1, paired = TRUE)
FvsR_q1 <- wilcox.test(x= feces1, y= rectum1, paired = TRUE)
CvsS_q1 <- wilcox.test(x= cloaca1, y= stomach1, paired = TRUE)
CvsI_q1 <- wilcox.test(x= cloaca1, y= intestine1, paired = TRUE)
CvsR_q1 <- wilcox.test(x= cloaca1, y= rectum1, paired = TRUE)

# Paired test (Wilcoxon) (q2)
feces2 <- subset(alpha, SampleType== "Feces", q2, drop = TRUE)
cloaca2 <- subset(alpha, SampleType== "Cloaca", q2, drop = TRUE)
stomach2 <-  subset(alpha, SampleType== "Stomach", q2, drop = TRUE)
intestine2 <- subset(alpha, SampleType== "Small intestine", q2, drop = TRUE)
rectum2 <-  subset(alpha, SampleType== "Rectum", q2, drop = TRUE)

FvsS_q2 <- wilcox.test(x= feces2, y= stomach2, paired = TRUE)
FvsI_q2 <- wilcox.test(x= feces2, y= intestine2, paired = TRUE)
FvsR_q2 <- wilcox.test(x= feces2, y= rectum2, paired = TRUE)
CvsS_q2 <- wilcox.test(x= cloaca2, y= stomach2, paired = TRUE)
CvsI_q2 <- wilcox.test(x= cloaca2, y= intestine2, paired = TRUE)
CvsR_q2 <- wilcox.test(x= cloaca2, y= rectum2, paired = TRUE)

## ORIGINAL PAPER ##
#Visualize: Specify the comparisons you want with paired tests

# Diversity order q=0
my_comparisons_q0 <- list(c("Feces", "Stomach"),
                          c("Cloaca", "Small intestine"),
                          c("Cloaca", "Rectum"))

# Crear un argumento para dejar la (q) en italica.
titulo0 <- expression(paste("Effective number of ASVs (", italic("q"), "=0)"))
HillNumb_q0 <- ggboxplot(alpha, x= "SampleType", y= "q0",
                         color = "black",  width = 0.6, lwd=0.3,
                         order = c("Stomach", "Small intestine", "Rectum", "Feces", "Cloaca"),
                         fill = c("#CC6677","#DDCC77","#44AA99", "#D55E00", "#888888")) +
  labs(x = element_blank(), y = titulo0) +
  theme_gray() + theme(text = element_text (size = 14)) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons_q0, label = "p.signif")
print(HillNumb_q0)

# Diversity order q=1
my_comparisons_q1 <- list(c("Feces", "Stomach"),
                          c("Feces", "Rectum"),
                          c("Cloaca", "Rectum"))

titulo1 <- expression(paste("Effective number of ASVs (", italic("q"), "=1)"))
HillNumb_q1 <- ggboxplot(alpha, x= "SampleType", y= "q1",
                         color = "black",  width = 0.6, lwd=0.3,
                         order = c("Stomach", "Small intestine", "Rectum", "Feces", "Cloaca"),
                         fill = c("#CC6677","#DDCC77","#44AA99", "#D55E00", "#888888")) +
  labs(x = element_blank(), y = titulo1) +
  theme_gray() + theme(text = element_text (size = 14)) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons_q1, label = "p.signif")
print(HillNumb_q1)

# Diversity order q=2
my_comparisons_q2 <- list(c("Feces", "Stomach"),
                          c("Feces", "Rectum"),
                          c("Cloaca", "Rectum"))

titulo2 <- expression(paste("Effective number of ASVs (", italic("q"), "=2)"))
HillNumb_q2 <- ggboxplot(alpha, x= "SampleType", y= "q2",
                         color = "black",  width = 0.6, lwd=0.3,
                         order = c("Stomach", "Small intestine", "Rectum", "Feces", "Cloaca"),
                         fill = c("#CC6677","#DDCC77","#44AA99", "#D55E00", "#888888")) +
  labs(x = element_blank(), y = titulo2) +
  theme_gray() + theme(text = element_text (size = 14)) +
  #theme(legend.position = "none",
   #     axis.ticks.x = element_blank(),
    #    axis.text.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons_q2, label = "p.signif")
print(HillNumb_q2)

library(cowplot)
Graphics_boxplot_final <- plot_grid(HillNumb_q0, HillNumb_q1, HillNumb_q2, 
                              nrow = 3, ncol = 1,
                              label_size = 13, rel_heights = c(1, 1, 1))
print(Graphics_boxplot_final)
ggsave("Graphics_boxplot_final.jpeg", width=6.0, height=11.0, dpi=300)
