setwd("~/Documents/MICROBIOTA REPTILES/PAPERS 2022/Second Article/Integrative_Zoology")

alpha <- read.csv("Hill_Numbers_q012.csv") %>% dplyr::select(SampleID, q0, q1, q2)
metadata <- read.csv("metadata.csv", check.names = F)
alpha <- alpha %>% inner_join(metadata, by = c("SampleID"="SampleID"))
#write.table(alpha, file="./Correlat_AlphaDiv.txt", sep = "\t")
alpha_div <- read.csv("Cor_AlphaDiv.csv", header = TRUE, check.names = F)

# Graficar la regresion (Stomach vs Feces) (q0)
Stom <- alpha_div %>% filter(SampleType== "Stomach") %>% mutate(stomach=q0)
Fec <- alpha_div %>% filter(SampleType== "Feces") %>% mutate(feces=q0)
join <- cbind(Stom, Fec)
join$ids <- c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10")

#Plot
data.Stom.Fec <- lm(stomach ~ feces, join)
Stom_Feces_q0 <- join %>% dplyr::select(stomach, feces, ids) %>% 
  ggplot(aes(x = stomach, y = feces, color = ids)) + 
  geom_point() +
  geom_abline(slope = coef(data.Stom.Fec)[[2]], 
              intercept = coef(data.Stom.Fec)[[1]]) +
  labs(title = paste("Adj R2 = ", signif(summary(data.Stom.Fec)$adj.r.squared, 5),
                     "Intercept =",signif(data.Stom.Fec$coef[[1]],5),
                     " Slope =",signif(data.Stom.Fec$coef[[2]], 5),
                     " P =",signif(summary(data.Stom.Fec)$coef[2,4], 5)))
print(Stom_Feces_q0)
ggsave("Stom_Feces_q0.jpeg", width=7, height=4.5, dpi=300)

# Graficar la regresion (Feces vs Small intestine) (q0)
Small <- alpha_div %>% filter(SampleType== "Small intestine") %>% 
  mutate(small_intestine=q0)
Fec <- alpha_div %>% filter(SampleType== "Feces") %>% mutate(feces=q0)
join_1 <- cbind(Fec, Small)
join_1$ids <- c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10")

#Plot
data.Small.Fec <- lm(small_intestine ~ feces, join_1)
Small_Feces_q0 <- join_1 %>% dplyr::select(small_intestine, feces, ids) %>% 
  ggplot(aes(x = small_intestine, y = feces, color = ids)) +
  geom_point() +
  geom_abline(slope = coef(data.Small.Fec)[[2]], 
              intercept = coef(data.Small.Fec)[[1]]) +
  labs(title = paste("Adj R2 = ",signif(summary(data.Small.Fec)$adj.r.squared, 5),
                     "Intercept =",signif(data.Small.Fec$coef[[1]],5),
                     " Slope =",signif(data.Small.Fec$coef[[2]], 5),
                     " P =",signif(summary(data.Small.Fec)$coef[2,4], 5)))
print(Small_Feces_q0)
ggsave("Small_Feces_q0.jpeg", width=7, height=4.5, dpi=300)

# Graficar la regresion (Feces vs Rectum) (q0)
Rectum <- alpha_div %>% filter(SampleType== "Rectum") %>% mutate(rectum=q0)
Fec <- alpha_div %>% filter(SampleType== "Feces") %>% mutate(feces=q0)
join_2 <- cbind(Fec, Rectum)
join_2$ids <- c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10")

#Plot
data.Rectum.Fec <- lm(rectum ~ feces, join_2)
Rectum_Feces_q0 <- join_2 %>% dplyr::select(rectum, feces, ids) %>% 
  ggplot(aes(x = rectum, y = feces, color = ids)) +
  geom_point() +
  geom_abline(slope = coef(data.Rectum.Fec)[[2]], 
              intercept = coef(data.Rectum.Fec)[[1]]) +
  labs(title = paste("Adj R2 = ",signif(summary(data.Rectum.Fec)$adj.r.squared, 5),
                     "Intercept =",signif(data.Rectum.Fec$coef[[1]],5),
                     " Slope =",signif(data.Rectum.Fec$coef[[2]], 5),
                     " P =",signif(summary(data.Rectum.Fec)$coef[2,4], 5)))
print(Rectum_Feces_q0)
ggsave("Rectum_Feces_q0.jpeg", width=7, height=4.5, dpi=300)

###############################################################################
# Graficar la regresion (Stomach vs Cloaca) (q0)
Stom <- alpha_div %>% filter(SampleType== "Stomach") %>% mutate(stomach=q0)
Cloa <- alpha_div %>% filter(SampleType== "Cloaca") %>% mutate(cloaca=q0)
join_3 <- cbind(Stom, Cloa)
join_3$ids <- c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10")

#Plot
data.Stom.Cloa <- lm(stomach ~ cloaca, join_3)
Stom_Cloaca_q0 <- join_3 %>% dplyr::select(stomach, cloaca, ids) %>% 
  ggplot(aes(x = stomach, y = cloaca, color = ids)) +
  geom_point() +
  geom_abline(slope = coef(data.Stom.Cloa)[[2]], 
              intercept = coef(data.Stom.Cloa)[[1]])+
  labs(title = paste("Adj R2 = ",signif(summary(data.Stom.Cloa)$adj.r.squared, 5),
                     "Intercept =",signif(data.Stom.Cloa$coef[[1]], 5),
                     " Slope =",signif(data.Stom.Cloa$coef[[2]], 5),
                     " P =",signif(summary(data.Stom.Cloa)$coef[2,4], 5)))
print(Stom_Cloaca_q0)
ggsave("Stom_Cloaca_q0.jpeg", width=7, height=4.5, dpi=300)

# Graficar la regresion (Samll intestine vs Cloaca) (q0)
Small <- alpha_div %>% filter(SampleType== "Small intestine") %>% mutate(small_intestine=q0)
Cloa <- alpha_div %>% filter(SampleType== "Cloaca") %>% mutate(cloaca=q0)
join_4 <- cbind(Small, Cloa)
join_4$ids <- c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10")

#Plot
data.Small.Cloa <- lm(small_intestine ~ cloaca, join_4)
Small_Cloaca_q0 <- join_4 %>% dplyr::select(small_intestine, cloaca, ids) %>% 
  ggplot(aes(x = small_intestine, y = cloaca, color = ids)) +
  geom_point() +
  geom_abline(slope = coef(data.Small.Cloa)[[2]], 
              intercept = coef(data.Small.Cloa)[[1]])+
  labs(title = paste("Adj R2 = ",signif(summary(data.Small.Cloa)$adj.r.squared, 5),
                     "Intercept =",signif(data.Small.Cloa$coef[[1]], 5),
                     " Slope =",signif(data.Small.Cloa$coef[[2]], 5),
                     " P =",signif(summary(data.Small.Cloa)$coef[2,4], 5)))
print(Small_Cloaca_q0)
ggsave("Small_Cloaca_q0.jpeg", width=7, height=4.5, dpi=300)

# Graficar la regresion (Rectum vs Cloaca) (q0)
Rectum <- alpha_div %>% filter(SampleType== "Rectum") %>% mutate(rectum=q0)
Cloa <- alpha_div %>% filter(SampleType== "Cloaca") %>% mutate(cloaca=q0)
join_5 <- cbind(Rectum, Cloa)
join_5$ids <- c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10")

#Plot
data.Rectum.Cloa <- lm(rectum ~ cloaca, join_5)
Rectum_Cloaca_q0 <- join_5 %>% dplyr::select(rectum, cloaca, ids) %>% 
  ggplot(aes(x = rectum, y = cloaca, color = ids)) +
  geom_point() +
  geom_abline(slope = coef(data.Rectum.Cloa)[[2]], 
              intercept = coef(data.Rectum.Cloa)[[1]]) +
  labs(title = paste("Adj R2 = ",signif(summary(data.Rectum.Cloa)$adj.r.squared, 5),
                     "Intercept =",signif(data.Rectum.Cloa$coef[[1]], 5),
                     " Slope =",signif(data.Rectum.Cloa$coef[[2]], 5),
                     " P =",signif(summary(data.Rectum.Cloa)$coef[2,4], 5)))
print(Rectum_Cloaca_q0)
ggsave("Rectum_Cloaca_q0.jpeg", width=7, height=4.5, dpi=300)
