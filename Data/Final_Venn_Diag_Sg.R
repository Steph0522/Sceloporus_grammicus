library(tidyverse)
# Loading files
setwd("~/Documents/MICROBIOTA REPTILES/PAPERS 2022/Second Article/Integrative_Zoology")

# Core microbiota (50%)
cloaca_50 <- read.csv("Core_50_Cloaca.csv",
                      check.names = F)
feces_50 <- read.csv("Core_50_Feces.csv",
                     check.names = F)
rectum_50 <- read.csv("Core_50_Rectum.csv",
                      check.names = F)
intestine_50 <- read.csv("Core_50_SmallIntestine.csv",
                         check.names = F) 
stomach_50 <- read.delim("Core_50_Stomach.csv",
                         check.names = F)

# Create Venn Diagramm
library(VennDiagram)

venn.plot_50 <- venn.diagram(
  x = list(Cloaca = cloaca_50$OTUID,
           Feces = feces_50$OTUID,
           Rectum = rectum_50$OTUID,
           Intestine = intestine_50$OTUID,
           Stomach = stomach_50$OTUID),
  category.names = c(
    expression(bold("Cloaca")),
    expression(bold("Feces")),
    expression(bold("Rectum")),
    expression(bold("Small intestine")),
    expression(bold("Stomach"))),
  filename = "viendo_50.tiff",
  output = TRUE,
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  units = "px",
  lwd = 6,
  lty = "blank",
  fill = c("#888888", "#D55E00", "#44AA99", "#DDCC77", "#CC6677"),
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 115, -125, -155),
  cat.dist = c(0.055, 0.055, 0.075, 0.060, 0.04),
  cat.fontfamily = "sans")
