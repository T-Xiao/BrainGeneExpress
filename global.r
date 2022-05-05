#options(repos = BiocManager::repositories())
#rsconnect::appDependencies()
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.10/bioc"))
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.10/bioc"))
library(BiocManager)
options(repos = BiocManager::repositories())
#getOption("repos")
library("shiny")
library(gridExtra)
library(grid)
library(ggplot2)
library(candisc)
library(shiny)
library(shinydashboard)
library(openxlsx)
library(rms)
library(zoo)
library("rpart")
library("partykit")
library(dplyr)
library(magrittr)
library(formattable)
library(yacca)
library(DT)
library(DEGreport)
library(diceR)
library(pheatmap)
library(tidyr)
library(Hmisc)
#library(readxl)
#library(xlsx)
library(CCA)
#library("NMF")
source("http://www.statpower.net/R312/CanCorr.r")
library("rgl")
library("scatterplot3d") # load
library("xtable")
library("rpart")
library("partykit")
library("survminer")
library("rlist")
library("MASS")
library("BBmisc") # is.error

TCGA_GBM_u133a_overall= readRDS(
  "data/TCGA_GBM_u133a_overall.rds"
)
# read in GBM rna seq data: 
TCGA_GBM_rna_seq_overall = readRDS(
  "data/TCGA_GBM_rna_seq_overall.rds"
)
#RNAoverallgenenames <- colnames(TCGA_GBM_rna_seq_overall)[1:20028]
#check how many columns are numeric:
#dim(dplyr::select_if(TCGA_GBM_rna_seq_overall[,1:20028], is.numeric))

# read in GBM agilent data:
TCGA_GBM_agilent_overall=readRDS("data/TCGA_GBM_agilent_overall.rds")
#agilentoverallgenenames <- colnames(TCGA_GBM_agilent_overall)[1:17814]
# read in LGG data:(" https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.6.100.0.0.tar.gz  repos=NULL, type="source")
TCGA_LGG_overall=readRDS("data/TCGA_LGG_overall.rds")
dim(dplyr::select_if(TCGA_LGG_overall[,1:20225], is.numeric))

Ivygap=readRDS("data/all_numeric.rds")

U133aoverallgenenames <- read.csv("data/U133aoverallgenenames.csv")
RNAoverallgenenames <- read.csv("data/RNAoverallgenenames.csv")
agilentoverallgenenames <- read.csv("data/agilentoverallgenenames.csv")
lggoverallgenenames <- read.csv("data/lggoverallgenenames.csv")
ivygapgenenames<-read.csv("data/ivygapgenenames.csv")

