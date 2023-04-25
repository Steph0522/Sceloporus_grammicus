# Get work directory
setwd("~/Documents/MICROBIOTA REPTILES/PAPERS 2022/Second Article/Integrative_Zoology")

# Load libraries and files
library(qiime2R)
library(hillR)
library(hilldiv)
library(tidyverse)
library(ggpubr)

# FUNCTIONAL DIVERSITY
funciones <- read_tsv("EC_predicted.tsv") %>% 
  arrange(desc(sequence)) %>% column_to_rownames(var = "sequence")
tabla <- data.frame(read_q2biom("feature-table.biom")) %>%
  rownames_to_column(var = "sequence") %>% arrange(desc(sequence)) %>% 
  column_to_rownames(  var = "sequence") %>%t() %>% as.data.frame()

# Run functional diversity by Hill numbers (at different q orders)
func_q0 <- hill_func(comm = tabla, traits = funciones, q = 0)
func_q1 <- hill_func(comm = tabla, traits = funciones, q = 1)
func_q2 <- hill_func(comm = tabla, traits = funciones, q = 2)

# Saved as table
Hill_Funct <- cbind(func_q0, func_q1, func_q2)
write.table(func_q0, file="./hill_funct_q0.txt", sep = "\t")
write.table(func_q1, file="./hill_funct_q1.txt", sep = "\t")
write.table(func_q2, file="./hill_funct_q2.txt", sep = "\t")

## ALPHA FUNCTIONAL DIVERSITY
alpha <- read.csv("Functional_div.csv", header = TRUE, check.names = F)
metadata <- read.csv("metadata.csv", check.names = F)
Div_Funct <- alpha %>% inner_join(metadata, by = c("SampleID"="SampleID"))

# Normality test
shapiro.test(x = Div_Funct$MD_q0)
shapiro.test(x = Div_Funct$MD_q1)
shapiro.test(x = Div_Funct$MD_q2)

hist(alpha$MD_q0)
hist(alpha$MD_q1)
hist(alpha$MD_q2)
# Data are not normal

# Paired test (Wilcoxon) (q0)
feces <- subset(Div_Funct, SampleType== "Feces", MD_q0, drop = TRUE)
cloaca <- subset(Div_Funct, SampleType== "Cloaca", MD_q0, drop = TRUE)
stomach <-  subset(Div_Funct, SampleType== "Stomach", MD_q0, drop = TRUE)
intestine <- subset(Div_Funct, SampleType== "Small intestine", MD_q0, drop = TRUE)
rectum <-  subset(Div_Funct, SampleType== "Rectum", MD_q0, drop = TRUE)


FvsS_q0 <- wilcox.test(x= feces, y= stomach, paired = TRUE)
FvsI_q0 <- wilcox.test(x= feces, y= intestine, paired = TRUE)
FvsR_q0 <- wilcox.test(x= feces, y= rectum, paired = TRUE)
CvsS_q0 <- wilcox.test(x= cloaca, y= stomach, paired = TRUE)
CvsI_q0 <- wilcox.test(x= cloaca, y= intestine, paired = TRUE)
CvsR_q0 <- wilcox.test(x= cloaca, y= rectum, paired = TRUE)

# Paired test (Wilcoxon) (q1)
feces1 <- subset(Div_Funct, SampleType== "Feces", MD_q1, drop = TRUE)
cloaca1 <- subset(Div_Funct, SampleType== "Cloaca", MD_q1, drop = TRUE)
stomach1 <-  subset(Div_Funct, SampleType== "Stomach", MD_q1, drop = TRUE)
intestine1 <- subset(Div_Funct, SampleType== "Small intestine", MD_q1, drop = TRUE)
rectum1 <-  subset(Div_Funct, SampleType== "Rectum", MD_q1, drop = TRUE)

FvsS_q1 <- wilcox.test(x= feces1, y= stomach1, paired = TRUE)
FvsI_q1 <- wilcox.test(x= feces1, y= intestine1, paired = TRUE)
FvsR_q1 <- wilcox.test(x= feces1, y= rectum1, paired = TRUE)
CvsS_q1 <- wilcox.test(x= cloaca1, y= stomach1, paired = TRUE)
CvsI_q1 <- wilcox.test(x= cloaca1, y= intestine1, paired = TRUE)
CvsR_q1 <- wilcox.test(x= cloaca1, y= rectum1, paired = TRUE)

# Paired test (Wilcoxon) (q2)
feces2 <- subset(Div_Funct, SampleType== "Feces", MD_q2, drop = TRUE)
cloaca2 <- subset(Div_Funct, SampleType== "Cloaca", MD_q2, drop = TRUE)
stomach2 <-  subset(Div_Funct, SampleType== "Stomach", MD_q2, drop = TRUE)
intestine2 <- subset(Div_Funct, SampleType== "Small intestine", MD_q2, drop = TRUE)
rectum2 <-  subset(Div_Funct, SampleType== "Rectum", MD_q2, drop = TRUE)

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
                          c("Cloaca", "Rectum"))

# Crear un argumento para dejar la (q) en italica.
tituloA <- expression(paste("Mean functional diversity (", italic("q"), "=0)"))
HillNumb_q0 <- ggboxplot(Div_Funct, x= "SampleType", y= "MD_q0",
                         color = "black",  width = 0.6, lwd=0.3,
                         order = c("Stomach", "Small intestine", "Rectum", "Feces", "Cloaca"),
                         fill = c("#CC6677","#DDCC77","#44AA99", "#D55E00", "#888888")) +
  labs(x = element_blank(), y = tituloA) +
  theme_gray() + theme(text = element_text (size = 14)) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_q0, label = "p.signif")
print(HillNumb_q0)

# Diversity order q=1
my_comparisons_q1 <- list(c("Feces", "Stomach"),
                          c("Cloaca", "Small intestine"),
                          c("Cloaca", "Rectum"))

tituloB <- expression(paste("Mean functional diversity (", italic("q"), "=1)"))
HillNumb_q1 <- ggboxplot(Div_Funct, x= "SampleType", y= "MD_q1",
                         color = "black",  width = 0.6, lwd=0.3,
                         order = c("Stomach", "Small intestine", "Rectum", "Feces", "Cloaca"),
                         fill = c("#CC6677","#DDCC77","#44AA99", "#D55E00", "#888888")) +
  labs(x = element_blank(), y = tituloB) +
  theme_gray() + theme(text = element_text (size = 14)) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons_q1, label = "p.signif")
print(HillNumb_q1)

# Diversity order q=2
my_comparisons_q2 <- list(c("Feces", "Stomach"),
                          c("Cloaca", "Small intestine"),
                          c("Cloaca", "Rectum"))

tituloC <- expression(paste("Mean functional diversity (", italic("q"), "=2)"))
HillNumb_q2 <- ggboxplot(Div_Funct, x= "SampleType", y= "MD_q2",
                         color = "black",  width = 0.6, lwd=0.3,
                         order = c("Stomach", "Small intestine", "Rectum", "Feces", "Cloaca"),
                         fill = c("#CC6677","#DDCC77","#44AA99", "#D55E00", "#888888")) +
  labs(x = element_blank(), y = tituloC) +
  theme_gray() + theme(text = element_text (size = 14)) +
  #theme(legend.position = "none",
  #     axis.ticks.x = element_blank(),
  #    axis.text.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons_q2, label = "p.signif")
print(HillNumb_q2)

library(cowplot)
Boxplot_Final_FunctDiv <- plot_grid(HillNumb_q0, HillNumb_q1, HillNumb_q2, 
                                    nrow = 3, ncol = 1,
                                    label_size = 13, rel_heights = c(1, 1, 1))
print(Boxplot_Final_FunctDiv)
ggsave("Boxplot_Final_FunctDiv.jpeg", width=6.0, height=11.0, dpi=300)
