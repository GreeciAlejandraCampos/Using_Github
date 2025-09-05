# Using_Github
tarea jeje 

title: "tarea 2"
author: "greecia"
date: "2025-09-05"
output: html_document
---

```{r}
# UBICACIÓN EN CARPETA
getwd()
```

# INSTALACIÓN DE PAQUETES

# BIOCOUNDUTOR
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "DOSE", "enrichplot"), ask = FALSE)
```

# CRAN
```{r}
install.packages(c("tidyverse", "readr", "dplyr"))
install.packages("tidyverse")
```

# CLUSTERPROFILE
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("DOSE")  

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("kohonen")
BiocManager::install("gplots")
BiocManager::install("VennDiagram")
BiocManager::install("dendsort")
```


# APLICANDO PAQUETERÍAS
```{r}
library(DESeq2);	#Normalization and everything related to that
library(tidyverse)  # For data manipulation (optional, useful for merging and plotting)
library(readr)      # For reading TSV/CSV files
library(dplyr)      # For filtering and data manipulation
library(clusterProfiler)  # For enrichment analysis
library(DOSE)         # Helper functions for enrichment
library(enrichplot)   # Visualization of enrichment results

library(limma);
library(edgeR);
library(pheatmap); # Pretty heatmaps
library(ggplot2); # For transparency in the colors
library(ggrepel);
library(matrixStats);
library(kohonen);	# This is the library for the SOM (Self-Organizing Maps)
library(gplots);	# Easy heatmaps también pheat más facile
library(VennDiagram);	
library(dendsort);	# Sorting (clasificación) dendrograms → Rotate dendograms
```

# CARGA DE TABLAS ----
```{r}
diff_genes <- read.csv("https://raw.githubusercontent.com/lucinoboa/Rust_Transcriptomics_Practice/refs/heads/main/Up_DEG_G2_vs_G1_strict.csv")

colnames(diff_genes)[1] <- "ID"       # Rename first column
diff_genes <- diff_genes[, c("ID", "log2FoldChange")]

#annotation <- readr::read_delim("C:/Users/edgar/OneDrive - Instituto Tecnologico y de Estudios Superiores de #Monterrey/Documents/R_learning/functional_anotation/fullAnnotation2.txt",
#                                delim = "\t",
#                                col_types = readr::cols()  # o simplemente omítelo para que lo infiera
#)
```
## ANÁLISIS DE ENRIQUESIMIENTO
```{r}
head(annotation, 100)
```
# Assume at least these columns exist: "ID", "GO", "description"
```{r}
head(annotation)

deg_annot <- annotation %>%
  filter(ID %in% diff_genes$ID)

deg_genes <- deg_annot$ID  

all_genes <- annotation$ID  

deg_ID <- deg_annot$ID
term2gene <- annotation[, c("GOs", "ID")]  # ajusta "GOs" al nombre exacto de tu columna GO
term2name <- annotation[, c("GOs", "Description")]  # opcional

ora_GOs <- enricher(
  gene = deg_ID,         # vector de genes, no GO
  universe = all_genes,        # vector de todos los genes detectados
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  TERM2GENE = term2gene,
  TERM2NAME = term2name
)
```

## GRAFICAR

sum(ora_GOs@result$p.adjust < 0.05)

# gráfico 1
```{r}
dotplot(ora_GOs, showCategory = 12)
#Qué tan significativo es el dato, síntesis de metabolitos secundarios
```


# gráfico 2
```{r}
barplot(ora_GOs, showCategory = 12)



ora_GOs <- pairwise_termsim(ora_GOs, method = "JC")

cog_dict <- c(
  "C"="Energy production and conversion", "D"="Cell cycle control, cell division, chromosome partitioning",
  "E"="Amino acid transport and metabolism", "F"="Nucleotide transport and metabolism",
  "G"="Carbohydrate transport and metabolism", "H"="Coenzyme transport and metabolism",
  "I"="Lipid transport and metabolism", "J"="Translation, ribosomal structure and biogenesis",
  "K"="Transcription", "L"="Replication, recombination and repair",
  "M"="Cell wall/membrane/envelope biogenesis", "N"="Cell motility",
  "O"="Posttranslational modification, protein turnover, chaperones",
  "P"="Inorganic ion transport and metabolism", "Q"="Secondary metabolites biosynthesis, transport and catabolism",
  "R"="General function prediction only", "S"="Function unknown",
  "T"="Signal transduction mechanisms", "U"="Intracellular trafficking, secretion, vesicular transport",
  "V"="Defense mechanisms", "W"="Extracellular structures", "Y"="Nuclear structure", "Z"="Cytoskeleton"
)

annotation <- annotation %>%
  mutate(COG_name = cog_dict[COG_category])
```


## Exportando datos

```{r}
deg_G2vsG1 <- read.csv("https://raw.githubusercontent.com/lucinoboa/Rust_Transcriptomics_Practice/refs/heads/main/DEG_G2_vs_G1_strict.csv")

deg_G3vsG1 <- read.csv("https://raw.githubusercontent.com/lucinoboa/Rust_Transcriptomics_Practice/refs/heads/main/DEG_G3_vs_G1_strict.csv")

deg_G2vsG3 <- read.csv("https://raw.githubusercontent.com/lucinoboa/Rust_Transcriptomics_Practice/refs/heads/main/DEG_G2_vs_G3_strict.csv")


deg_up <- subset(deg_G2vsG1, log2FoldChange > 1)
deg_down <- subset(deg_G2vsG1, log2FoldChange < -1)

deg_up_annot <- deg_up %>%
  as_tibble(rownames = "ID") %>%  # convierte rownames a columna "ID"
  left_join(dplyr::select(annotation, ID, COG_category, COG_name), by = "ID")

deg_down_annot <- deg_down %>%
  as_tibble(rownames = "ID") %>%
  left_join(dplyr::select(annotation, ID, COG_category, COG_name), by = "ID")


universe_summary <- annotation %>%
  filter(!is.na(COG_name)) %>%
  group_by(COG_name) %>%
  summarise(Universe = n_distinct(ID))

up_summary <- deg_up_annot %>%
  filter(!is.na(COG_name)) %>%
  group_by(COG_name) %>%
  summarise(Up = n_distinct(ID))

down_summary <- deg_down_annot %>%
  filter(!is.na(COG_name)) %>%
  group_by(COG_name) %>%
  summarise(Down = n_distinct(ID))
```

# Merge the summaries.
```{r}
summary_table <- universe_summary %>%
  full_join(up_summary, by="COG_name") %>%
  full_join(down_summary, by="COG_name") %>%
  replace(is.na(.), 0) %>%
  arrange(desc(Universe))



print(summary_table, n=21)
```

