Cnidarian Marker Gene Analysis - Updated
================
Zoe Dellaert
2026-11-07

- [Overview](#overview)
  - [Inputs](#inputs)
- [paths relative to
  `HI_genome_annotations`](#paths-relative-to-hi_genome_annotations)
  - [Outputs](#outputs)
- [1. Setup](#1-setup)
- [2. Marker Gene Sources](#2-marker-gene-sources)
  - [Sources](#sources)
  - [File locations](#file-locations)
  - [Load file metadata](#load-file-metadata)
- [3. Load and standardize marker gene
  lists](#3-load-and-standardize-marker-gene-lists)
  - [3.1. Read in all marker files](#31-read-in-all-marker-files)
  - [3.2: Species-specific dataframe
    cleaning](#32-species-specific-dataframe-cleaning)
  - [3.3 First-pass a combined marker
    list](#33-first-pass-a-combined-marker-list)
- [4. Match marker list gene names to protein IDs used in orthology
  analyses](#4-match-marker-list-gene-names-to-protein-ids-used-in-orthology-analyses)
  - [4.1 Load orthologs from HI_genome_annotations (see
    “Overview”)](#41-load-orthologs-from-hi_genome_annotations-see-overview)
  - [4.2 Load source species annotation
    files](#42-load-source-species-annotation-files)
  - [4.3 Match marker list gene names to protein IDs used in orthology
    analyses](#43-match-marker-list-gene-names-to-protein-ids-used-in-orthology-analyses)
  - [4.4 Join to standardized cluster names (manually standardized by
    me, csv key
    provided)](#44-join-to-standardized-cluster-names-manually-standardized-by-me-csv-key-provided)
  - [4.5 Achieved a *better* combined marker
    list!](#45-achieved-a-better-combined-marker-list)
- [5. Cross-species marker analysis via
  orthology](#5-cross-species-marker-analysis-via-orthology)
  - [5.1 Per-species alluvial plots: identify shared marker orthologs
    with other scRNA
    datasets](#51-per-species-alluvial-plots-identify-shared-marker-orthologs-with-other-scrna-datasets)
  - [5.2 Generate ortholog marker gene tables for my focal
    species](#52-generate-ortholog-marker-gene-tables-for-my-focal-species)

# Overview

I am compiling marker gene lists from published and preprinted cnidarian
single-cell atlases and using orthology analysis to try to find
consensus across studies and use as predictions for my three focus
species (*Pocillopora acuta* (Pacuta), *Montipora capitata* (Mcap), and
*Porites compressa* (Pcomp)).

## Inputs

# paths relative to `HI_genome_annotations`

| Input | Location |
|----|----|
| Published marker tables | `references/cnidarian_marker_genes/cnidarian_sc_marker_lists/` |
| Source-species annotation gff/gtf files (for ID mapping) | `references/cnidarian_marker_genes/cnidarian_sc_marker_lists/annotation_files/` |
| File metadata (paper, species, file names) | `references/cnidarian_marker_genes/sequence_download.csv` |
| Cluster name standardization key | `references/cnidarian_marker_genes/cnidarian_sc_marker_lists/cluster_key.csv` |
| Ortholog pairs | `references/cnidarian_orthologs/broccoli/dir_step4/orthologous_pairs.txt` |

## Outputs

| Output | Location |
|----|----|
| Combined cnidarian marker table | `references/cnidarian_marker_genes/combined_cnidarian_markers.csv` |
| ID-cleaned ortholog pairs | `references/cnidarian_marker_genes/output/orthologous_pairs_modified.txt` |
| Our species marker pairs | `references/cnidarian_marker_genes/output/{Mcap,Pcomp,Pacuta}/` |
| Cross-species cluster alluvial plots | `references/cnidarian_marker_genes/output/combined_markers/` |
| Marker-orthogroup UpSet plot | `references/cnidarian_marker_genes/output/combined_markers/upset_markers.png` |

# 1. Setup

``` r
library(here)
library(stringr)
library(readxl)
library(janitor)
library(rtracklayer)
library(ggalluvial)
library(circlize)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("krassowski/complex-upset")
#library(ComplexUpset)
library(data.table)
library(tidyverse)

getwd()
```

    ## [1] "/Users/zoedellaert/Documents/URI/HI_genome_annotations/code"

``` r
here::i_am("code/scRNA_markers.Rmd")

markers_dir <- here("references/cnidarian_marker_genes")
```

# 2. Marker Gene Sources

I gathered marker genes from the supplementary tables and provided
datasets (Github, Dryad) of nine published single-cell RNA-seq atlases
from cnidarians as of March 2025. **I am working on adding more recent
preprints + papers.** Where pre-computed marker lists were not provided
(Hydra, Xenia), markers were derived from deposited Seurat objects and
stored in this repository.

## Sources

| \# | Species | Citation | Source of markers |
|----|----|----|----|
| 1 | *Exaiptasia diaphana* (Edia) | [Cui et al. 2023](https://doi.org/10.1126/sciadv.adf7108) | Supplementary Data S3 |
| 2 | *Hydra vulgaris* (Hvul) | [Siebert et al. 2019](https://doi.org/10.1126/science.aav9314) | Re-derived markers from genome-mapped [Seurat RDS file (Hydra_Seurat_Whole_Genome.rds)](https://datadryad.org/dataset/doi:10.5061/dryad.v5r6077)) |
| 3 | *Hydractinia symbiolongicarpus* (Hsym) | [Salamanca-Díaz et al. 2025](https://doi.org/10.1038/s41467-025-57168-z) | Supplementary Data 3 |
| 4 | *Nematostella vectensis* (Nvec) | [Cole et al. 2024](https://doi.org/10.1186/s12983-024-00529-z) | Supplementary Information: Additional files 2-3, Supplementary Data 1-2 |
| 5 | *Astrangia poculata* (Apoc) | [Valadez-Ingersoll et al. 2024](https://doi.org/10.1101/2024.06.20.599951) | GitHub repository; Oculina mapped to *Astrangia* genome |
| 6 | *Orbicella faveolata* (Ofav) | [del Campo et al. 2024](https://doi.org/10.21203/rs.3.rs-5397639/v1) | Extended Data Tables 1-2 + re-derived in R |
| 7 | *Pocillopora damicornis* (Pdam) | [Diao et al. 2025](https://doi.org/10.3389/fmars.2025.1551763) | Supplementary Table S1 |
| 8 | *Stylophora pistillata* (Spis) | [Levy et al. 2021](https://doi.org/10.1016/j.cell.2021.04.005) | GitHub repository |
| 9 | *Xenia sp.* (Xesp) | [Hu et al. 2020](https://doi.org/10.1038/s41586-020-2385-7) | Converted [RDS file](https://cmo.carnegiescience.edu/endosymbiosis/data/cluster_markers.RDS) into csv |
| 10 | *Acropora millepora* (Amil) | [Levy et al. 2025](10.1038/s41586-025-09623-6) | GitHub repository |
| 11 | *Acropora muricata* (Amur) | [Han et al. 2025](10.1038/s42003-025-08089-6) | Supplementary Data 2 |

### Additional processing:

#### Siebert et al. Hydra markers

Siebert et al. published markers based on a *de novo* transcriptome but
also deposited a Seurat object mapped to the genome. Since downstream
ortholog analysis requires genome-based gene IDs, I re-derived markers
from the genome-mapped object. However, I did have to run
“UpdateSeuratObject”, “FindAllMarkers”, and filter by p-value to match
their published analysis. See code below, which is based on their
[github file
SA06_ClustGenome.Rmd](https://github.com/cejuliano/hydra_single_cell/blob/master/SA06_ClustGenome.Rmd).

``` r
#wget https://datadryad.org/downloads/file_stream/110292
Hvul_Genome <- readRDS("Hydra_Seurat_Whole_Genome.rds")
Hvul_Genome <- UpdateSeuratObject(object = Hvul_Genome)
Hvul_Genome

DimPlot(Hvul_Genome, reduction = "tsne",label=TRUE)

Hvul_genome_markers <- FindAllMarkers(object = Hvul_Genome, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
Hvul_genome_markers_p05 <- Hvul_genome_markers[Hvul_genome_markers$p_val_adj<0.05,]

write_csv(Hvul_genome_markers_p05,
file.path(markers_dir,"/cnidarian_sc_marker_lists/02-Hydra_Siebert2019_genome_markers.csv"))
```

## File locations

All of these files have been downloaded and stored in
`cnidarian_marker_genes/cnidarian_sc_marker_lists/`. All csv files were
converted to excel format for consistency to account for multiple sheets
in published supplementary tables. Extensive file metadata (species,
abbreviations, paper DOIs, file download links, original file names,
renamed files) for the marker genes, genomes, and genome annotations for
all of these species are recorded in
`cnidarian_marker_genes/sequence_download.csv`.

## Load file metadata

``` r
file_info <- read.csv(file.path(markers_dir, "/sequence_download.csv"))
file_info %>% dplyr::select(species_name,species_abbreviation,marker_file_name)
```

    ##                     species_name species_abbreviation
    ## 1            Exaiptasia_diaphana                 Edia
    ## 2                 Hydra_vulgaris                 Hvul
    ## 3  Hydractinia_symbiolongicarpus                 Hsym
    ## 4         Nematostella_vectensis                 Nvec
    ## 5             Astrangia_poculata                 Apoc
    ## 6            Orbicella_faveolata                 Ofav
    ## 7         Pocillopora_damicornis                 Pdam
    ## 8      Pocillopora_damicornis_RG               PdamRG
    ## 9          Stylophora_pistillata                 Spis
    ## 10                      Xenia_sp                 Xesp
    ## 11            Acropora_millepora                 Amil
    ## 12             Acropora_muricata                 Amur
    ## 13             Pocillopora_acuta                 Pacu
    ## 14             Porites_compressa                 Pcom
    ## 15            Montipora_capitata                 Mcap
    ##                                marker_file_name
    ## 1      01-Aiptasia_Cui2023_cluster_markers.xlsx
    ## 2      02-Hydra_Siebert2019_genome_markers.xlsx
    ## 3      03-Hydrac_SalamancaDiaz2025_markers.xlsx
    ## 4             04-Nvec_Cole2024_markers_all.xlsx
    ## 5  05-Oculina_Ingersoll2024_markerswithsym.xlsx
    ## 6  06-Ofav_Bonacolta2024_Ranalysis_markers.xlsx
    ## 7          07-Pdam_Diao2025_clustermarkers.xlsx
    ## 8                 07-PdamRG_clustermarkers.xlsx
    ## 9                 08-Spis_Levy2021_markers.xlsx
    ## 10                09-Xenia_cluster_markers.xlsx
    ## 11                10-Amil_celltype_markers.xlsx
    ## 12                11-Amur_celltype_markers.xlsx
    ## 13                                         <NA>
    ## 14                                         <NA>
    ## 15                                         <NA>

# 3. Load and standardize marker gene lists

All of these marker gene lists use different column names, cut-offs for
a gene to be called a “marker gene”, gene ID formats, etc. The code I
have written does the following:

1.  Read in each file based on sequence_download.csv, selecting only
    necessary sheets for some species
2.  For each file, systematically go through and rename columns to:
    `gene_id` \| `cluster` \| `source_sp` \| `protein_id`
3.  When log2FC and p-values are given, filter to high-confidence
    markers (log2FC \> 1, padj \< 0.05)
4.  For marker genes appearing in multiple clusters, keep only the
    cluster with the highest log2FC

## 3.1. Read in all marker files

Semi-standardized file import method that works for excel files with one
sheet or clusters spread across many sheets:

``` r
marker_files <- unlist(file_info$marker_file_name[!is.na(file_info$marker_file_name)])
file_names <- unlist(file_info$species_abbreviation[!is.na(file_info$marker_file_name)])
file_names <- paste0(file_names,"_markers")

for (i in 1:length(marker_files)){
  # find name of first marker gene file
  file_name <- paste0(markers_dir, "/cnidarian_sc_marker_lists/",file_info$marker_file_name[i])

  # find name of corresponding species (species which Genome was used, if different)
  species <- file_info$species_abbreviation[i]
  name <- paste0(species,"_markers")
  
    # fix the path to new Pdam file
  if (file_name==file.path(markers_dir, "cnidarian_sc_marker_lists/07-PdamRG_clustermarkers.xlsx")){
         file_name <- paste0("../../07-PdamRG_clustermarkers.xlsx")
  }

  # identify sheets in file; for any of these files with multiple sheets, each sheet is a different cluster/marker category
  sheets <- excel_sheets(file_name) 
  
  # remove sheets "sample.metadata", "ao.DEGs", and "ao.DETFs" from one problematic file that has three sheets formatted differently from its others
  if (file_name==file.path(markers_dir, "cnidarian_sc_marker_lists/04-Nvec_Cole2024_markers_all.xlsx")){sheets <- sheets[grep("ao",sheets, invert=TRUE)]}
  if (file_name==file.path(markers_dir, "cnidarian_sc_marker_lists/04-Nvec_Cole2024_markers_all.xlsx")){sheets <- sheets[grep("sample.metadata",sheets, invert=TRUE)]}

  assign(name, data.frame())

  # read in file
  for (sheet in sheets){
    df <- read_xlsx(file_name,sheet = sheet)
    
    df$sheet <- sheet
    assign(name, rbind(get(name),df))
  }
}
```

## 3.2: Species-specific dataframe cleaning

### Edia (Aiptasia)

``` r
name <- file_names[1] #"Processing file: Edia_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Edia_markers"

``` r
clean <- get(name)

# current header is junk, set first row to column names and clean column names with janitor package
clean <- clean %>%
  row_to_names(row_number = 1) %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster
clean <- clean %>% mutate_at(c("average_log_fc","adjusted_p_value"),as.double) %>%
                   filter(average_log_fc>1)

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(gene_id) %>%
                   filter(average_log_fc == max(average_log_fc)) %>%
                   ungroup()

# select columns with gene and cluster names
clean <- clean %>% select(gene_id,cluster)

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id      cluster source_sp protein_id  
    ##   <chr>        <chr>   <chr>     <chr>       
    ## 1 AIPGENE8252  2       Edia      AIPGENE8252 
    ## 2 AIPGENE25192 2       Edia      AIPGENE25192
    ## 3 AIPGENE5232  2       Edia      AIPGENE5232 
    ## 4 AIPGENE11887 2       Edia      AIPGENE11887
    ## 5 AIPGENE8505  3       Edia      AIPGENE8505

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Edia_markers"

### Hvul (Hydra)

``` r
name <- file_names[2] #"Processing file: Hvul_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Hvul_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster
clean <- clean %>% filter(avg_log2fc>1)

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(gene) %>%
                   filter(avg_log2fc == max(avg_log2fc)) %>%
                   ungroup()

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(gene,cluster) %>%
                   dplyr::rename(gene_id=gene)

#clean gene names to only have gene name from genome (delete everything after the "-")
clean <- clean %>% mutate(gene_id = sub("\\-.*", "", gene_id))

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id   cluster        source_sp protein_id
    ##   <chr>     <chr>          <chr>     <chr>     
    ## 1 g11134.t1 ecEp_stem_cell Hvul      g11134.t1 
    ## 2 g11543.t1 ecEp_stem_cell Hvul      g11543.t1 
    ## 3 g20377.t1 ecEp_stem_cell Hvul      g20377.t1 
    ## 4 g31502.t1 ecEp_stem_cell Hvul      g31502.t1 
    ## 5 g11131.t1 ecEp_stem_cell Hvul      g11131.t1

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Hvul_markers"

### Hsym (Hydractinia)

``` r
name <- file_names[3] #"Processing file: Hsym_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Hsym_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster and adjusted p-val < 0.05
clean <- clean %>% filter(logfoldchanges>1)  %>% filter(pvals_adj < 0.05) 

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(names) %>%
                   filter(logfoldchanges == max(logfoldchanges)) %>%
                   ungroup()

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(names,sheet) %>%
                   dplyr::rename(gene_id=names) %>%
                   dplyr::rename(cluster=sheet)

#clean cluster names to replace spaces with underscores
clean <- clean %>% mutate(cluster = sub(" ", "_", cluster))

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id      cluster   source_sp protein_id  
    ##   <chr>        <chr>     <chr>     <chr>       
    ## 1 LOC130630194 Cluster_1 Hsym      LOC130630194
    ## 2 LOC130641826 Cluster_1 Hsym      LOC130641826
    ## 3 LOC130613966 Cluster_1 Hsym      LOC130613966
    ## 4 LOC130629617 Cluster_1 Hsym      LOC130629617
    ## 5 LOC130635707 Cluster_1 Hsym      LOC130635707

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Hsym_markers"

### Nvec (Nematostella)

``` r
name <- file_names[4] #"Processing file: Nvec_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Nvec_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

# select only "partitionMarkers" sheet for now- there are other sheets for each cluster but they don't seem to be filtered for pval or log-foldchange.
#  select only markers with avg_log2FC > 1 in the cluster and adjusted p-val < 0.05
clean <- clean %>% filter(sheet=="PartitionMarkers") %>% filter(avg_log2fc>1)  %>% filter(p_val_adj< 0.05) 

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(gene_id) %>%
                   filter(avg_log2fc == max(avg_log2fc)) %>%
                   ungroup()
       
# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(gene_id,cluster)

#remove space in cluster name
clean <- clean %>% mutate(cluster = gsub(" ",".",cluster))

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id   cluster source_sp protein_id
    ##   <chr>     <chr>   <chr>     <chr>     
    ## 1 NV2.13863 pSC     Nvec      NV2.13863 
    ## 2 NV2.18050 pSC     Nvec      NV2.18050 
    ## 3 NV2.6574  pSC     Nvec      NV2.6574  
    ## 4 NV2.16402 pSC     Nvec      NV2.16402 
    ## 5 NV2.6575  pSC     Nvec      NV2.6575

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Nvec_markers"

### Apoc (Astrangia/Oculina)

``` r
name <- file_names[5] #"Processing file: Apoc_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Apoc_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster
clean <- clean %>% filter(avg_log2fc>1)

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(protein_id) %>%
                   filter(avg_log2fc == max(avg_log2fc)) %>%
                   ungroup()

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% filter(!(protein_id=="NA")) %>% #there are some weird NAs I am removing
                   select(protein_id,cluster) %>%
                   dplyr::rename(gene_id=protein_id) 

# set cluster column to be a character column
clean <- clean %>% mutate(cluster = as.character(cluster))

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id             cluster source_sp protein_id         
    ##   <chr>               <chr>   <chr>     <chr>              
    ## 1 evm.model.Ap4.5489  0       Apoc      evm.model.Ap4.5489 
    ## 2 evm.model.Ap13.2163 0       Apoc      evm.model.Ap13.2163
    ## 3 evm.model.Ap3.1613  0       Apoc      evm.model.Ap3.1613 
    ## 4 evm.model.Ap1.3218  0       Apoc      evm.model.Ap1.3218 
    ## 5 evm.model.Ap14.1598 0       Apoc      evm.model.Ap14.1598

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Apoc_markers"

### Ofav (Orbicella)

``` r
name <- file_names[6] #"Processing file: Ofav_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Ofav_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster
clean <- clean %>% filter(avg_log2fc>1)

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(gene) %>%
                   filter(avg_log2fc == max(avg_log2fc)) %>%
                   ungroup()

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(gene,cluster) %>%
                   dplyr::rename(gene_id=gene)

#clean cluster names to replace parentheses with underscores
clean <- clean %>% mutate(cluster = gsub(" ", "", cluster))
clean <- clean %>% mutate(cluster = gsub("\\(", "_", cluster))
clean <- clean %>% mutate(cluster = gsub("\\)", "", cluster))

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id      cluster      source_sp protein_id  
    ##   <chr>        <chr>        <chr>     <chr>       
    ## 1 QW917-006349 Unidentified Ofav      QW917-006349
    ## 2 QW917-002686 Unidentified Ofav      QW917-002686
    ## 3 QW917-002688 Unidentified Ofav      QW917-002688
    ## 4 QW917-025087 Unidentified Ofav      QW917-025087
    ## 5 QW917-010184 Unidentified Ofav      QW917-010184

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Ofav_markers"

### Pdam (Pocillopora)

``` r
name <- file_names[7] #"Processing file: Pdam_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Pdam_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster
clean <- clean %>% filter(avg_log2fc>1)

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(x1) %>%
                   filter(avg_log2fc == max(avg_log2fc)) %>%
                   ungroup()

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(x1,sheet) %>%
                   dplyr::rename(gene_id=x1) %>%
                   dplyr::rename(cluster=sheet)

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id      cluster  source_sp protein_id  
    ##   <chr>        <chr>    <chr>     <chr>       
    ## 1 LOC113666310 Cluster0 Pdam      LOC113666310
    ## 2 LOC131792517 Cluster0 Pdam      LOC131792517
    ## 3 LOC113673831 Cluster0 Pdam      LOC113673831
    ## 4 LOC136280560 Cluster0 Pdam      LOC136280560
    ## 5 LOC113679584 Cluster0 Pdam      LOC113679584

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Pdam_markers"

### Pdam (Pocillopora) - ReefGenomics Version

``` r
name <- file_names[8] #"Processing file: PdamRG_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: PdamRG_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster
clean <- clean %>% filter(marker_log2fc>1)

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(gene) %>%
                   filter(marker_log2fc == max(marker_log2fc)) %>%
                   ungroup()

Pdam_RG_filtered <- clean 
# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(gene,cluster) %>%
                   dplyr::rename(gene_id=gene)  %>%
                   mutate(cluster=as.character(cluster))

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- str_replace_all(clean$gene_id,"-", "_")
Pdam_RG_filtered$protein_id <- str_replace_all(Pdam_RG_filtered$gene,"-", "_") 
Pdam_RG_filtered <- Pdam_RG_filtered %>% mutate(cluster=as.character(cluster))

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id       cluster source_sp protein_id   
    ##   <chr>         <chr>   <chr>     <chr>        
    ## 1 pdam-00007194 0       PdamRG    pdam_00007194
    ## 2 pdam-00007071 0       PdamRG    pdam_00007071
    ## 3 pdam-00012072 0       PdamRG    pdam_00012072
    ## 4 pdam-00013116 0       PdamRG    pdam_00013116
    ## 5 pdam-00017857 0       PdamRG    pdam_00017857

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: PdamRG_markers"

### Spis (Stylophora)

``` r
name <- file_names[9] #"Processing file: Spis_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Spis_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(gene,cell_type) %>%
                   dplyr::rename(gene_id=gene) %>%
                   dplyr::rename(cluster=cell_type)

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id             cluster   source_sp protein_id         
    ##   <chr>               <chr>     <chr>     <chr>              
    ## 1 Spis_XP_022797422_1 cnidocyte Spis      Spis_XP_022797422_1
    ## 2 Spis_XP_022789499_1 cnidocyte Spis      Spis_XP_022789499_1
    ## 3 Spis_XP_022795719_1 cnidocyte Spis      Spis_XP_022795719_1
    ## 4 Spis_XP_022787371_1 cnidocyte Spis      Spis_XP_022787371_1
    ## 5 Spis18343_1         cnidocyte Spis      Spis18343_1

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Spis_markers"

### Xesp (Xenia)

``` r
name <- file_names[10] #"Processing file: Xesp_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Xesp_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster
clean <- clean %>% filter(avg_log_fc>1)

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(gene) %>%
                   filter(avg_log_fc == max(avg_log_fc)) %>%
                   ungroup()

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(gene,cluster) %>%
                   dplyr::rename(gene_id=gene) 

# set cluster column to be a character column
clean <- clean %>% mutate(cluster = as.character(cluster))
  
# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id                               cluster source_sp protein_id            
    ##   <chr>                                 <chr>   <chr>     <chr>                 
    ## 1 Xe-006061                             1       Xesp      Xe-006061             
    ## 2 Xe-006047                             1       Xesp      Xe-006047             
    ## 3 Common central domain of tyrosinase-1 1       Xesp      Common central domain…
    ## 4 Xe-014936                             1       Xesp      Xe-014936             
    ## 5 SSPO-8                                1       Xesp      SSPO-8

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Xesp_markers"

### Amil (Acropora)

``` r
name <- file_names[11] #"Processing file: Amil_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Amil_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with avg_log2FC > 1 in the cluster
clean <- clean %>% filter(avg_log2fc>1) %>% filter(p_val_adj<0.05) 

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(gene,focus_node) %>%
                   dplyr::rename(gene_id=gene) %>%
                   dplyr::rename(cluster=focus_node)

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id              cluster        source_sp protein_id          
    ##   <chr>                <chr>          <chr>     <chr>               
    ## 1 Amil_Amillepora26052 gastrodermis_1 Amil      Amil_Amillepora26052
    ## 2 Amil_Amillepora21884 gastrodermis_1 Amil      Amil_Amillepora21884
    ## 3 Amil_Amillepora10946 gastrodermis_1 Amil      Amil_Amillepora10946
    ## 4 Amil_Amillepora14445 gastrodermis_1 Amil      Amil_Amillepora14445
    ## 5 Amil_Amillepora12321 gastrodermis_1 Amil      Amil_Amillepora12321

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Amil_markers"

### Amur (Acorpora muricata)

``` r
name <- file_names[12] #"Processing file: Amur_markers"
print(paste0("Processing file: ", name))
```

    ## [1] "Processing file: Amur_markers"

``` r
clean <- get(name)

#  clean column names with janitor package
clean <- clean %>%
  clean_names() #clean column names with janitor package

#  select only markers with *Unwounded* avg_log2FC > 1 in the cluster and *Unwounded* adjusted p-val < 0.05
clean <- clean %>% filter(unwounded_avg_log2fc>1)  %>% filter(unwounded_p_val_adj < 0.05) 

# identify duplicate genes (pointing to > 1 cluster), and keep only the cluster to which that gene has the highest log2foldchange
clean <- clean %>% group_by(gene_id) %>%
                   filter(unwounded_avg_log2fc == max(unwounded_avg_log2fc)) %>%
                   ungroup()

# select columns with gene and cluster names and rename to standard format
clean <- clean %>% select(gene_id,sheet) %>%
                   dplyr::rename(cluster=sheet)

#clean cluster names to replace spaces with underscores
clean <- clean %>% mutate(cluster = sub(" ", "_", cluster))

# add species source column
clean$source_sp <- sub("_markers","",name)
clean$protein_id <- clean$gene_id

# assign the species name to the dataframe
assign(name,clean)
head(get(name),5)
```

    ## # A tibble: 5 × 4
    ##   gene_id             cluster     source_sp protein_id         
    ##   <chr>               <chr>       <chr>     <chr>              
    ## 1 evm.model.Chr12.6   Symbiocytes Amur      evm.model.Chr12.6  
    ## 2 evm.model.Chr4.674  Symbiocytes Amur      evm.model.Chr4.674 
    ## 3 evm.model.Chr1.1138 Symbiocytes Amur      evm.model.Chr1.1138
    ## 4 evm.model.Chr3.417  Symbiocytes Amur      evm.model.Chr3.417 
    ## 5 evm.model.Chr8.332  Symbiocytes Amur      evm.model.Chr8.332

``` r
print(paste0("Processed file: ", name))
```

    ## [1] "Processed file: Amur_markers"

## 3.3 First-pass a combined marker list

``` r
glimpse(lapply(file_names,get))
```

    ## List of 12
    ##  $ : tibble [649 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:649] "AIPGENE8252" "AIPGENE25192" "AIPGENE5232" "AIPGENE11887" ...
    ##   ..$ cluster   : chr [1:649] "2" "2" "2" "2" ...
    ##   ..$ source_sp : chr [1:649] "Edia" "Edia" "Edia" "Edia" ...
    ##   ..$ protein_id: chr [1:649] "AIPGENE8252" "AIPGENE25192" "AIPGENE5232" "AIPGENE11887" ...
    ##  $ : tibble [5,589 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:5589] "g11134.t1" "g11543.t1" "g20377.t1" "g31502.t1" ...
    ##   ..$ cluster   : chr [1:5589] "ecEp_stem_cell" "ecEp_stem_cell" "ecEp_stem_cell" "ecEp_stem_cell" ...
    ##   ..$ source_sp : chr [1:5589] "Hvul" "Hvul" "Hvul" "Hvul" ...
    ##   ..$ protein_id: chr [1:5589] "g11134.t1" "g11543.t1" "g20377.t1" "g31502.t1" ...
    ##  $ : tibble [5,031 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:5031] "LOC130630194" "LOC130641826" "LOC130613966" "LOC130629617" ...
    ##   ..$ cluster   : chr [1:5031] "Cluster_1" "Cluster_1" "Cluster_1" "Cluster_1" ...
    ##   ..$ source_sp : chr [1:5031] "Hsym" "Hsym" "Hsym" "Hsym" ...
    ##   ..$ protein_id: chr [1:5031] "LOC130630194" "LOC130641826" "LOC130613966" "LOC130629617" ...
    ##  $ : tibble [955 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:955] "NV2.13863" "NV2.18050" "NV2.6574" "NV2.16402" ...
    ##   ..$ cluster   : chr [1:955] "pSC" "pSC" "pSC" "pSC" ...
    ##   ..$ source_sp : chr [1:955] "Nvec" "Nvec" "Nvec" "Nvec" ...
    ##   ..$ protein_id: chr [1:955] "NV2.13863" "NV2.18050" "NV2.6574" "NV2.16402" ...
    ##  $ : tibble [4,210 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:4210] "evm.model.Ap4.5489" "evm.model.Ap13.2163" "evm.model.Ap3.1613" "evm.model.Ap1.3218" ...
    ##   ..$ cluster   : chr [1:4210] "0" "0" "0" "0" ...
    ##   ..$ source_sp : chr [1:4210] "Apoc" "Apoc" "Apoc" "Apoc" ...
    ##   ..$ protein_id: chr [1:4210] "evm.model.Ap4.5489" "evm.model.Ap13.2163" "evm.model.Ap3.1613" "evm.model.Ap1.3218" ...
    ##  $ : tibble [19,773 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:19773] "QW917-006349" "QW917-002686" "QW917-002688" "QW917-025087" ...
    ##   ..$ cluster   : chr [1:19773] "Unidentified" "Unidentified" "Unidentified" "Unidentified" ...
    ##   ..$ source_sp : chr [1:19773] "Ofav" "Ofav" "Ofav" "Ofav" ...
    ##   ..$ protein_id: chr [1:19773] "QW917-006349" "QW917-002686" "QW917-002688" "QW917-025087" ...
    ##  $ : tibble [8,433 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:8433] "LOC113666310" "LOC131792517" "LOC113673831" "LOC136280560" ...
    ##   ..$ cluster   : chr [1:8433] "Cluster0" "Cluster0" "Cluster0" "Cluster0" ...
    ##   ..$ source_sp : chr [1:8433] "Pdam" "Pdam" "Pdam" "Pdam" ...
    ##   ..$ protein_id: chr [1:8433] "LOC113666310" "LOC131792517" "LOC113673831" "LOC136280560" ...
    ##  $ : tibble [2,331 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:2331] "pdam-00007194" "pdam-00007071" "pdam-00012072" "pdam-00013116" ...
    ##   ..$ cluster   : chr [1:2331] "0" "0" "0" "0" ...
    ##   ..$ source_sp : chr [1:2331] "PdamRG" "PdamRG" "PdamRG" "PdamRG" ...
    ##   ..$ protein_id: chr [1:2331] "pdam_00007194" "pdam_00007071" "pdam_00012072" "pdam_00013116" ...
    ##  $ : tibble [491 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:491] "Spis_XP_022797422_1" "Spis_XP_022789499_1" "Spis_XP_022795719_1" "Spis_XP_022787371_1" ...
    ##   ..$ cluster   : chr [1:491] "cnidocyte" "cnidocyte" "cnidocyte" "cnidocyte" ...
    ##   ..$ source_sp : chr [1:491] "Spis" "Spis" "Spis" "Spis" ...
    ##   ..$ protein_id: chr [1:491] "Spis_XP_022797422_1" "Spis_XP_022789499_1" "Spis_XP_022795719_1" "Spis_XP_022787371_1" ...
    ##  $ : tibble [965 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:965] "Xe-006061" "Xe-006047" "Common central domain of tyrosinase-1" "Xe-014936" ...
    ##   ..$ cluster   : chr [1:965] "1" "1" "1" "1" ...
    ##   ..$ source_sp : chr [1:965] "Xesp" "Xesp" "Xesp" "Xesp" ...
    ##   ..$ protein_id: chr [1:965] "Xe-006061" "Xe-006047" "Common central domain of tyrosinase-1" "Xe-014936" ...
    ##  $ : tibble [15,371 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:15371] "Amil_Amillepora26052" "Amil_Amillepora21884" "Amil_Amillepora10946" "Amil_Amillepora14445" ...
    ##   ..$ cluster   : chr [1:15371] "gastrodermis_1" "gastrodermis_1" "gastrodermis_1" "gastrodermis_1" ...
    ##   ..$ source_sp : chr [1:15371] "Amil" "Amil" "Amil" "Amil" ...
    ##   ..$ protein_id: chr [1:15371] "Amil_Amillepora26052" "Amil_Amillepora21884" "Amil_Amillepora10946" "Amil_Amillepora14445" ...
    ##  $ : tibble [1,017 × 4] (S3: tbl_df/tbl/data.frame)
    ##   ..$ gene_id   : chr [1:1017] "evm.model.Chr12.6" "evm.model.Chr4.674" "evm.model.Chr1.1138" "evm.model.Chr3.417" ...
    ##   ..$ cluster   : chr [1:1017] "Symbiocytes" "Symbiocytes" "Symbiocytes" "Symbiocytes" ...
    ##   ..$ source_sp : chr [1:1017] "Amur" "Amur" "Amur" "Amur" ...
    ##   ..$ protein_id: chr [1:1017] "evm.model.Chr12.6" "evm.model.Chr4.674" "evm.model.Chr1.1138" "evm.model.Chr3.417" ...

``` r
markers_all <- bind_rows(lapply(file_names,get))
```

**Note:** Diao et al Pdam markers report only unadjusted p-values; Spis
markers have no fold-change/p-value info. All other lists are filtered
to LFC \>1 and padj \<0.05

# 4. Match marker list gene names to protein IDs used in orthology analyses

Because things can never be simple, many of the gene IDs used in the
marker gene lists are not the same as the protein IDs in the protein
FASTA files used to run the ortholog analysis. Therefore, here I gather
all of the genome annotation gtf/gff files for these genomes and create
a reference so marker gene lists can be matched back to the genome
annotations and ortholog analysis output files.

## 4.1 Load orthologs from HI_genome_annotations (see “Overview”)

``` r
broc_groups <- data.table::fread(here("references/cnidarian_orthologs/broccoli/dir_step3/table_OGs_protein_names.txt"), header=T)
broc_unclassified <- read.delim(here("references/cnidarian_orthologs/broccoli/dir_step3/unclassified_proteins.txt"),header=F) 
broc_pairs <- data.table::fread(here("references/cnidarian_orthologs/broccoli/dir_step4/orthologous_pairs.txt"), header=F)
colnames(broc_pairs) <- c("Gene_spA","Gene_spB")
```

## 4.2 Load source species annotation files

### Code used to download all files

``` bash
# make directory
cd cnidarian_marker_genes/cnidarian_sc_marker_lists/
mkdir annotation_files
cd annotation_files

# 1 - Aiptasia, assembled by Baumgarten et al., 2015 (10.1073/pnas.1513318112)
wget http://aiptasia.reefgenomics.org/download/aiptasia_genome.gff3.gz

# 2 - Hydra, assembled by Steele et al 2015 (https://research.nhgri.nih.gov/hydra/)
wget https://research.nhgri.nih.gov/hydra/download/genemodels_gff3/hydra2.0_genemodels.gff3.gz

# 3 - Hydractinia, assembled by Kon-Nanjo et al 2023 (10.1093/g3journal/jkad107)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/227/915/GCF_029227915.1_HSymV2.1/GCF_029227915.1_HSymV2.1_genomic.gtf.gz

# 4 - Nematostella, assembled by Zimmermann et al., 2023 (10.1038/s41467-023-44080-7)
wget https://simrbase.stowers.org/files/pub/nematostella/Nvec/genomes/Nvec200/aligned/tcs_v2/20240221/NV2g.20240221.gff
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs12983-024-00529-z/MediaObjects/12983_2024_529_MOESM3_ESM.xlsx

# 5 - Astrangia, assembled by Stankiewicz et al., 2025 (10.1093/g3journal/jkaf033)
wget -O apoculata.gtf.gz https://zenodo.org/records/14110456/files/apoculata.gtf.gz?download=1

# 6 - Orbicella, assembled by Young et al., 2024 (10.1186/s12864-024-10092-w)
wget -O Orbicella_faveolata_gen_17.gff3 https://zenodo.org/records/10151798/files/Orbicella_faveolata_gen_17.gff3?download=1

# 7 - Pocillopora damicornis, assembled by Cunning et al., 2018 (10.1038/s41598-018-34459-8)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_genomic.gtf.gz
#reefgenomics version
wget http://pdam.reefgenomics.org/download/pdam_annotation.gff3.gz

# 8 - Stylophora, assembled by Voolstra et al., 2017 (10.1038/s41598-017-17484-x)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/571/385/GCF_002571385.2_Stylophora_pistillata_v1.1/GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gtf.gz

# 9 - Xenia, assembled by Hu et al., 2020 (10.1038/s41586-020-2385-7)
wget https://cmo.carnegiescience.edu/endosymbiosis/genome/xenSp1.gtf
wget https://cmo.carnegiescience.edu/endosymbiosis/data/id2symbol.txt

# 10 - A. millepora - download key used by Grau-Bove et al 2025 to compare to Refseq IDs
wget https://raw.githubusercontent.com/sebepedroslab/oculina-coral-sc-atlas/refs/heads/master/data/reference/Amil_Genes_ncbi_dataset.tsv

# 11 - A. muricata - download genome assembly preprint key for original "evm" gene IDs to match Refseq IDs
wget https://assets-eu.researchsquare.com/files/rs-9352493/v1/3f7379b55e3fce8c4e24fd83.xlsx # sheet 2 of this (S2.1_Am_gene annotation) saved as csv -- called 11-Amur_annot.csv

gunzip *

mv aiptasia_genome.gff3 01-Edia_annot.gff3
mv hydra2.0_genemodels.gff3 02-Hvul_annot.gff3
mv GCF_029227915.1_HSymV2.1_genomic.gtf 03-Hsym_annot.gtf
mv NV2g.20240221.gff 04-Nvec_annot.gff
mv apoculata.gtf 05-Apoc_annot.gtf
mv Orbicella_faveolata_gen_17.gff3 06-Ofav_annot.gff3
mv GCF_003704095.1_ASM370409v1_genomic.gtf 07-Pdam_annot.gtf
mv pdam_annotation.gff3 07-PdamRG_annot.gff3
mv GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gtf 08-Spis_annot.gtf
mv xenSp1.gtf 09-Xesp_annot.gtf

mv 12983_2024_529_MOESM3_ESM.xlsx 04-Nvec_Cole2024_markers_top.xlsx
mv id2symbol.txt 09-Xesp_annot.txt

mv Amil_Genes_ncbi_dataset.tsv 10-Amil_annot.tsv
```

### Species-specific processing

``` r
# 2 - Nematostella - read in gene ID annotation file
Nvec_names <- read_xlsx(file.path(markers_dir, "cnidarian_sc_marker_lists/annotation_files/04-Nvec_Cole2024_markers_top.xlsx"),sheet = 1)

# 3 - Hydractinia, assembled by Kon-Nanjo et al 2023 (10.1093/g3journal/jkad107)
Hsym_annot <- rtracklayer::import(paste0(markers_dir, "/cnidarian_sc_marker_lists/annotation_files/",file_info$gtf_file_name[3]))
Hsym_prot <- as.data.frame(Hsym_annot) %>%
                        filter(type == "CDS") %>% 
                        select(protein_id, gene_id) %>% 
                        distinct()

# 7 - Pocillopora damicornis, assembled by Cunning et al., 2018 (10.1038/s41598-018-34459-8)
Pdam_annot <- rtracklayer::import(paste0(markers_dir, "/cnidarian_sc_marker_lists/annotation_files/",file_info$gtf_file_name[7]))
Pdam_prot <- as.data.frame(Pdam_annot) %>%
                        filter(type == "CDS") %>% 
                        select(protein_id, gene_id) %>% 
                        distinct()

Pdam_transcripts <- as.data.frame(Pdam_annot) %>%
                        filter(type == "CDS") %>% 
                        select(protein_id, transcript_id) %>% 
                        distinct()

# 9 - Xenia, assembled by Hu et al., 2020 (10.1038/s41586-020-2385-7)
Xesp_annot <- rtracklayer::import(paste0(markers_dir, "/cnidarian_sc_marker_lists/annotation_files/",file_info$gtf_file_name[10]))
Xesp_prot <- as.data.frame(Xesp_annot) %>%
                        filter(type == "CDS") %>% 
                        select(transcript_id, gene_id,gene_name) %>% 
                        distinct()
    ## join gtf info to annotation and fix symbol names for the few names that have dashes where there should be underscores
    Xesp_annot_names <- read_delim(paste0(markers_dir, "/cnidarian_sc_marker_lists/annotation_files/",file_info$gene_annotation_file_name[10]), 
        delim = "\t", escape_double = FALSE, 
        trim_ws = TRUE)
    
    Xesp_prot <- Xesp_prot %>% left_join(Xesp_annot_names, by=c("gene_id"="ID")) %>% 
                             mutate(symbol=str_replace_all(symbol,"-","_")) %>% select(transcript_id,symbol)

    
## 10 - Amil
Amil_annot_names <- read_delim(paste0(markers_dir, "/cnidarian_sc_marker_lists/annotation_files/",file_info$gene_annotation_file_name[11]), 
        delim = "\t", escape_double = FALSE, skip = 5, col_names=FALSE,
        trim_ws = TRUE) 
   
colnames(Amil_annot_names) = c("old_ID", "protein_id", "mRNA_id", "annotation")
Amil_annot_names$old_ID <- str_replace_all(Amil_annot_names$old_ID, "-RA", "")

## 11 - Amur
Amur_annot_names <- read_delim(paste0(markers_dir, "/cnidarian_sc_marker_lists/annotation_files/",file_info$gene_annotation_file_name[12]),trim_ws = TRUE) %>% select(gene_id,NCBI_protein_ID)

# remove prefixes such as .1 and .2 from the end of certain gene ids to match the evm ids
Amur_annot_names <- Amur_annot_names %>% mutate(gene_id=
                                                  paste(str_split_i(gene_id,"\\.",1),
                                                         str_split_i(gene_id,"\\.",2),
                                                         str_split_i(gene_id,"\\.",3),
                                                        str_split_i(gene_id,"\\.",4),
                                                        sep="."))


# remove all the large gtf objects 
rm(list = ls(pattern = "_annot$"))
```

## 4.3 Match marker list gene names to protein IDs used in orthology analyses

### What do we have before cleaning

``` r
for (i in 1:length(marker_files)){
  # find name of corresponding species (species which Genome was used, if different)
  species <- file_info$species_abbreviation[i]
  name <- paste0(species,"_markers")
  protein_file <- file_info$protein_file_shortname[i]
  
  print(paste0(name," marker gene format:"))
  
  #get an idea of what the marker gene names look like
  print(head(get(name)$protein_id))
  
  print(paste0(name," fasta format:"))
  #get an idea of what the protein names look like from the protein fasta (read into Broccoli)
  col <- which(colnames(broc_groups) == protein_file)
  print(unlist(head(broc_groups[,..col],3),use.names = F))
  
  print(paste0(name," matches of gene id to broccoli protein names:"))
  print(length(unique(intersect(intersect(get(name)$protein_id,broc_pairs$Gene_spA),intersect(get(name)$protein_id,broc_pairs$Gene_spB)))))
}
```

    ## [1] "Edia_markers marker gene format:"
    ## [1] "AIPGENE8252"  "AIPGENE25192" "AIPGENE5232"  "AIPGENE11887" "AIPGENE8505" 
    ## [6] "AIPGENE4933" 
    ## [1] "Edia_markers fasta format:"
    ## [1] ""                                      
    ## [2] "AIPGENE14349 AIPGENE14348 AIPGENE14327"
    ## [3] ""                                      
    ## [1] "Edia_markers matches of gene id to broccoli protein names:"
    ## [1] 420
    ## [1] "Hvul_markers marker gene format:"
    ## [1] "g11134.t1" "g11543.t1" "g20377.t1" "g31502.t1" "g11131.t1" "g28273.t1"
    ## [1] "Hvul_markers fasta format:"
    ## [1] "Sc4wPfr_749.1.g22025.t1" ""                       
    ## [3] ""                       
    ## [1] "Hvul_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Hsym_markers marker gene format:"
    ## [1] "LOC130630194" "LOC130641826" "LOC130613966" "LOC130629617" "LOC130635707"
    ## [6] "LOC130644998"
    ## [1] "Hsym_markers fasta format:"
    ## [1] ""               ""               "XP_057301293.1"
    ## [1] "Hsym_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Nvec_markers marker gene format:"
    ## [1] "NV2.13863"  "NV2.18050"  "NV2.6574"   "NV2.16402"  "NV2.6575"  
    ## [6] "PCNA_XENLA"
    ## [1] "Nvec_markers fasta format:"
    ## [1] ""                "NV2t017494001.1" ""               
    ## [1] "Nvec_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Apoc_markers marker gene format:"
    ## [1] "evm.model.Ap4.5489"  "evm.model.Ap13.2163" "evm.model.Ap3.1613" 
    ## [4] "evm.model.Ap1.3218"  "evm.model.Ap14.1598" "evm.model.Ap9.2740" 
    ## [1] "Apoc_markers fasta format:"
    ## [1] "protein|evm.model.Ap1.1821"                           
    ## [2] "protein|evm.model.Ap1.3252 protein|evm.model.Ap1.3251"
    ## [3] "protein|evm.model.Ap1.3733"                           
    ## [1] "Apoc_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Ofav_markers marker gene format:"
    ## [1] "QW917-006349" "QW917-002686" "QW917-002688" "QW917-025087" "QW917-010184"
    ## [6] "QW917-018150"
    ## [1] "Ofav_markers fasta format:"
    ## [1] "QW917_029422-T1"                 "QW917_029954-T1 QW917_029955-T1"
    ## [3] "QW917_029580-T1"                
    ## [1] "Ofav_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Pdam_markers marker gene format:"
    ## [1] "LOC113666310" "LOC131792517" "LOC113673831" "LOC136280560" "LOC113679584"
    ## [6] "LOC131780936"
    ## [1] "Pdam_markers fasta format:"
    ## [1] "XP_027053049.1"                              
    ## [2] "XP_027055855.1 XP_027055862.1 XP_027055571.1"
    ## [3] "XP_027054032.1"                              
    ## [1] "Pdam_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "PdamRG_markers marker gene format:"
    ## [1] "pdam_00007194" "pdam_00007071" "pdam_00012072" "pdam_00013116"
    ## [5] "pdam_00017857" "pdam_00020216"
    ## [1] "PdamRG_markers fasta format:"
    ## [1] "pdam_00020748-RA" "pdam_00010709-RA" "pdam_00007905-RA"
    ## [1] "PdamRG_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Spis_markers marker gene format:"
    ## [1] "Spis_XP_022797422_1" "Spis_XP_022789499_1" "Spis_XP_022795719_1"
    ## [4] "Spis_XP_022787371_1" "Spis18343_1"         "Spis_XP_022789491_1"
    ## [1] "Spis_markers fasta format:"
    ## [1] "XP_022805022.1"                                                                                          
    ## [2] "XP_022805551.1 XP_022805557.1 XP_022805558.1 XP_022805559.1 XP_022805560.1 XP_022805563.1 XP_022805561.1"
    ## [3] "XP_022794131.1"                                                                                          
    ## [1] "Spis_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Xesp_markers marker gene format:"
    ## [1] "Xe-006061"                            
    ## [2] "Xe-006047"                            
    ## [3] "Common central domain of tyrosinase-1"
    ## [4] "Xe-014936"                            
    ## [5] "SSPO-8"                               
    ## [6] "tyrosinase-2"                         
    ## [1] "Xesp_markers fasta format:"
    ## [1] "Xe_002198-T1" ""             "Xe_000860-T1"
    ## [1] "Xesp_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Amil_markers marker gene format:"
    ## [1] "Amil_Amillepora26052" "Amil_Amillepora21884" "Amil_Amillepora10946"
    ## [4] "Amil_Amillepora14445" "Amil_Amillepora12321" "Amil_Amillepora23199"
    ## [1] "Amil_markers fasta format:"
    ## [1] "XP_044164386.1"                                                            
    ## [2] "XP_029192641.2 XP_029192647.2 XP_029192648.2 XP_029192568.2 XP_029192569.2"
    ## [3] "XP_029179328.2 XP_029179335.2"                                             
    ## [1] "Amil_markers matches of gene id to broccoli protein names:"
    ## [1] 0
    ## [1] "Amur_markers marker gene format:"
    ## [1] "evm.model.Chr12.6"    "evm.model.Chr4.674"   "evm.model.Chr1.1138" 
    ## [4] "evm.model.Chr3.417"   "evm.model.Chr8.332"   "evm.model.Chr10.1157"
    ## [1] "Amur_markers fasta format:"
    ## [1] "XP_067056449.1 XP_067056450.1"               
    ## [2] "XP_067058083.1 XP_067058084.1 XP_067055809.1"
    ## [3] "XP_067016878.1 XP_067016879.1"               
    ## [1] "Amur_markers matches of gene id to broccoli protein names:"
    ## [1] 0

### Two-fold cleaning

This chunk combines two types of cleaning I should probably clarify:

1.  Make “broc_pairs_clean” to update the ortholog gene IDs to match the
    marker gene names
2.  Change marker gene names to match protein IDs using gtf-based tables

``` r
# Correct protein IDs used in broccoli and orthofinder to match gene ids 
broc_pairs_clean <- broc_pairs

# 1 - Edia: all gene ids match 

# 2 - Hvul: remove chromosome name "Sc4wPfr_" from the beginning of each protein name in the fasta
broc_pairs_clean <- broc_pairs_clean %>% mutate(Gene_spA = case_when(str_detect(Gene_spA,"Sc4wPfr_") ~ 
                                                                str_remove(Gene_spA, ".*(?=g)"),
                                                               .default = Gene_spA)) %>% 
                                         mutate(Gene_spB = case_when(str_detect(Gene_spB,"Sc4wPfr_") ~ 
                                                                str_remove(Gene_spB, ".*(?=g)"),
                                                               .default = Gene_spB))

# 3 - Hsym: used NCBI GTF file to match gene id to protein id
Hsym_markers <- Hsym_markers %>% select(-protein_id) %>% left_join(Hsym_prot)
  
# 4 - Nvec: gene symbol interpretation info found at https://simrbase.stowers.org/starletseaanemone
broc_pairs_clean <- broc_pairs_clean %>%
                                 #change gene name prefix and remove any leading zeros before gene number
                                   mutate(Gene_spA = str_replace(Gene_spA,"NV2t0*","NV2.")) %>% 
                                   mutate(Gene_spB = str_replace(Gene_spB,"NV2t0*","NV2.")) %>% 
                                #remove 00X.1 transcript ID not present in gene ids (but only for Nvec genes)
                                   mutate(Gene_spA = case_when(str_detect(Gene_spA,"NV2") ~ 
                                                                str_remove(Gene_spA, "\\d\\d\\d\\.\\d"),
                                                               .default = Gene_spA)) %>%
                                   mutate(Gene_spB = case_when(str_detect(Gene_spB,"NV2") ~ 
                                                                str_remove(Gene_spB, "\\d\\d\\d\\.\\d"),
                                                               .default = Gene_spB))

# 5 - Apoc: remove "protein|" prefix from protein IDs and then they match the marker gene ids
broc_pairs_clean <- broc_pairs_clean %>% mutate(Gene_spA = str_remove(Gene_spA,"protein\\|")) %>% 
                                   mutate(Gene_spB = str_remove(Gene_spB,"protein\\|")) 

# 6 - Ofav: remove "-T#" suffix from protein IDs and then they match the marker gene ids (but only for Ofav ids beginning with QW)
broc_pairs_clean <- broc_pairs_clean %>% mutate(Gene_spA = case_when(str_detect(Gene_spA,"QW917") ~ 
                                                                str_remove(Gene_spA, "\\-T\\d"),
                                                               .default = Gene_spA)) %>% 
                                         mutate(Gene_spB = case_when(str_detect(Gene_spB,"QW917") ~ 
                                                                str_remove(Gene_spB, "\\-T\\d"),
                                                               .default = Gene_spB))
    #replace all the dashes with underscores
    Ofav_markers <- Ofav_markers %>% mutate(protein_id = str_replace_all(gene_id, "-", "_"))

# 7 - Pdam - used NCBI GTF file to match gene id to protein id

      # NOTE, come back to this, because many genes were not assigned protein IDs. this is likely because the single cell paper mapped to transcriptome data in addition to the genome
      Pdam_markers <- Pdam_markers %>% select(-protein_id) %>% left_join(Pdam_prot)

# 8 - PdamRG: remove "-RA#" suffix from protein IDs and then they match the marker gene ids
broc_pairs_clean <- broc_pairs_clean %>% mutate(Gene_spA = case_when(str_detect(Gene_spA,"pdam_") ~ 
                                                                str_remove(Gene_spA, "\\-RA"),
                                                               .default = Gene_spA)) %>% 
                                         mutate(Gene_spB = case_when(str_detect(Gene_spB,"pdam_") ~ 
                                                                str_remove(Gene_spB, "\\-RA"),
                                                               .default = Gene_spB))

# 9 - Spis
      # NOTE, come back to this, because some genes were not assigned protein IDs. these inclue orphan peaks and gene ids with "LOC____"
      
Spis_markers <- Spis_markers %>% mutate(protein_id = str_remove(gene_id,"Spis_"))
Spis_markers <- Spis_markers %>% mutate(protein_id = str_replace(protein_id,"_1$",".1"))

# 10 - Xesp: replace all the dashes with underscores and join to annotation dataframe
Xesp_markers <- Xesp_markers %>% mutate(protein_id = str_replace_all(gene_id, "-", "_")) %>%
                                 left_join(Xesp_prot, by=c("protein_id"="symbol")) %>%                           
                                 select(gene_id,cluster,source_sp,transcript_id)  %>% dplyr::rename(protein_id=transcript_id)

# 11 - Amil
Amil_markers <- Amil_markers %>% select(-protein_id) %>% left_join(Amil_annot_names, by=c("gene_id"="old_ID")) %>%                 
                                 select(gene_id,cluster,source_sp,protein_id) %>% filter(!is.na(protein_id))

# 12 - Amur

# 15 rows have weird IDs where there are two ids concatenated -- evm.model.Chr6.198_evm.model.Chr6.199 is evm.model.Chr6.198_evm in the key
# so we just need to replace any _evm*** with _evm

Amur_markers <- Amur_markers %>% mutate(gene_id = str_replace(gene_id, "_evm\\.model.*$","_evm"))

Amur_markers <- Amur_markers %>% select(-protein_id) %>% left_join(Amur_annot_names) %>%  dplyr::rename(protein_id=NCBI_protein_ID) %>% 
                                 select(gene_id,cluster,source_sp,protein_id) %>% filter(!is.na(protein_id))
```

### What do we have after cleaning

``` r
for (i in 1:length(marker_files)){
  # find name of corresponding species (species which Genome was used, if different)
  species <- file_info$species_abbreviation[i]
  name <- paste0(species,"_markers")
  protein_file <- file_info$protein_file_shortname[i]

  print(paste0(name," number of unique gene ids in marker list:"))
  print(length(unique(get(name)$protein_id)))
  print(paste0(name," matches of gene id to broccoli protein names:"))
  print(length(unique(union(intersect(get(name)$protein_id,broc_pairs_clean$Gene_spA),intersect(get(name)$protein_id,broc_pairs_clean$Gene_spB)))))
  
  print(paste0(name," matches of gene id to *unclassified* broccoli protein names:"))
  print(length(unique(intersect(get(name)$protein_id,broc_unclassified$V1))))
}
```

    ## [1] "Edia_markers number of unique gene ids in marker list:"
    ## [1] 649
    ## [1] "Edia_markers matches of gene id to broccoli protein names:"
    ## [1] 494
    ## [1] "Edia_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 142
    ## [1] "Hvul_markers number of unique gene ids in marker list:"
    ## [1] 5589
    ## [1] "Hvul_markers matches of gene id to broccoli protein names:"
    ## [1] 4582
    ## [1] "Hvul_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 0
    ## [1] "Hsym_markers number of unique gene ids in marker list:"
    ## [1] 8147
    ## [1] "Hsym_markers matches of gene id to broccoli protein names:"
    ## [1] 6942
    ## [1] "Hsym_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 1112
    ## [1] "Nvec_markers number of unique gene ids in marker list:"
    ## [1] 955
    ## [1] "Nvec_markers matches of gene id to broccoli protein names:"
    ## [1] 764
    ## [1] "Nvec_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 0
    ## [1] "Apoc_markers number of unique gene ids in marker list:"
    ## [1] 4210
    ## [1] "Apoc_markers matches of gene id to broccoli protein names:"
    ## [1] 4010
    ## [1] "Apoc_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 0
    ## [1] "Ofav_markers number of unique gene ids in marker list:"
    ## [1] 19773
    ## [1] "Ofav_markers matches of gene id to broccoli protein names:"
    ## [1] 4761
    ## [1] "Ofav_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 0
    ## [1] "Pdam_markers number of unique gene ids in marker list:"
    ## [1] 2801
    ## [1] "Pdam_markers matches of gene id to broccoli protein names:"
    ## [1] 2790
    ## [1] "Pdam_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 5
    ## [1] "PdamRG_markers number of unique gene ids in marker list:"
    ## [1] 2331
    ## [1] "PdamRG_markers matches of gene id to broccoli protein names:"
    ## [1] 2288
    ## [1] "PdamRG_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 0
    ## [1] "Spis_markers number of unique gene ids in marker list:"
    ## [1] 491
    ## [1] "Spis_markers matches of gene id to broccoli protein names:"
    ## [1] 356
    ## [1] "Spis_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 4
    ## [1] "Xesp_markers number of unique gene ids in marker list:"
    ## [1] 1030
    ## [1] "Xesp_markers matches of gene id to broccoli protein names:"
    ## [1] 741
    ## [1] "Xesp_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 278
    ## [1] "Amil_markers number of unique gene ids in marker list:"
    ## [1] 5019
    ## [1] "Amil_markers matches of gene id to broccoli protein names:"
    ## [1] 5004
    ## [1] "Amil_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 10
    ## [1] "Amur_markers number of unique gene ids in marker list:"
    ## [1] 1003
    ## [1] "Amur_markers matches of gene id to broccoli protein names:"
    ## [1] 1002
    ## [1] "Amur_markers matches of gene id to *unclassified* broccoli protein names:"
    ## [1] 0

## 4.4 Join to standardized cluster names (manually standardized by me, csv key provided)

``` r
cluster_key <- read_csv(file.path(markers_dir, "cnidarian_sc_marker_lists/cluster_key.csv"))

for (name in file_names) {
  assign(name, 
         get(name) %>% 
           left_join(cluster_key, join_by(source_sp == Species, cluster == Cluster_ID)), 
         envir = .GlobalEnv)
}

markers_all_annot <- bind_rows(lapply(file_names,get))
```

## 4.5 Achieved a *better* combined marker list!

``` r
glimpse(markers_all_annot)
```

    ## Rows: 62,840
    ## Columns: 5
    ## $ gene_id           <chr> "AIPGENE8252", "AIPGENE25192", "AIPGENE5232", "AIPGE…
    ## $ cluster           <chr> "2", "2", "2", "2", "3", "3", "3", "4", "4", "4", "4…
    ## $ source_sp         <chr> "Edia", "Edia", "Edia", "Edia", "Edia", "Edia", "Edi…
    ## $ protein_id        <chr> "AIPGENE8252", "AIPGENE25192", "AIPGENE5232", "AIPGE…
    ## $ Standardized_Name <chr> "gastrodermis", "gastrodermis", "gastrodermis", "gas…

``` r
# Save the files generated above to be able to quickly read them in from here to run chunks below
write.csv(markers_all_annot, file.path(markers_dir, "/combined_cnidarian_markers_full.csv"), row.names =FALSE)

write.csv(markers_all_annot %>% filter(source_sp != "PdamRG"), file.path(markers_dir, "/combined_cnidarian_markers.csv"), row.names =FALSE)

write.csv(broc_pairs_clean, file.path(markers_dir, "/output/orthologous_pairs_modified.txt"), row.names =FALSE)
```

# 5. Cross-species marker analysis via orthology

## 5.1 Per-species alluvial plots: identify shared marker orthologs with other scRNA datasets

``` r
broc_pairs_clean <-  data.table::fread(file.path(markers_dir, "/output/orthologous_pairs_modified.txt"), header=T)
markers_all_annot <- read.csv(file.path(markers_dir, "/combined_cnidarian_markers_full.csv"))

for (i in 1:length(file_names)){
  species <- sub("_markers","",file_names[i])
  
  combined <- get(file_names[i]) %>% 
    left_join(broc_pairs_clean, by = c("protein_id"="Gene_spA")) %>%
    bind_rows(
      get(file_names[i]) %>% left_join(broc_pairs_clean, by = c("protein_id"="Gene_spB")) %>%
        dplyr::rename(Gene_spB = Gene_spA)) %>% 
    filter(!is.na(Gene_spB)) %>%  
    dplyr::rename(protein_id_spB = Gene_spB) %>%
    rename_with(~paste0(.x,"_spA"),c(gene_id,Standardized_Name,protein_id))
  
  assign(paste0("combined_",file_names[i]),combined)
  write.csv(get(paste0("combined_",file_names[i])),paste0(markers_dir, "/output/combined_markers/",paste0("combined_",file_names[i]),".csv"))

  combined_clust <- combined %>%
    left_join(markers_all_annot, by = c("protein_id_spB"="protein_id")) %>% 
    filter(!is.na(source_sp.y)) %>% #keep only genes in other species that are marker genes
    dplyr::rename(source_sp = source_sp.x) %>%
    dplyr::rename(source_spB = source_sp.y) %>%
    rename_with(~paste0(.x,"_spB"),c(gene_id,Standardized_Name)) 
  
  wide <- combined_clust %>% 
    group_by(across(contains("spA")),source_spB) %>%
    summarise(Standardized_Name_spB = ifelse(length(unique(Standardized_Name_spB)) > 1,
                            "multiple_clusters",
                            unique(Standardized_Name_spB)[1]), .groups = "drop") %>% 
    pivot_wider(names_from = source_spB, values_from = Standardized_Name_spB, values_fill = NA) 

  # figure out which species this shares the most marker genes with and make plots for the top 3

close_sp <- wide %>%
    summarize(across(-contains("spA"), ~sum(. != "unidentified" & !is.na(.)))) %>%
    pivot_longer(everything(), names_to = "species", values_to = "non_NA_count") %>%
    arrange(desc(non_NA_count)) %>% slice_head(n = 3) %>% pull(species)

g <- ggplot(data = wide %>% filter(if_all(close_sp, ~!is.na(.))),  
       aes(axis1 = Standardized_Name_spA, axis2 = get(close_sp[1]), axis3 = get(close_sp[2]), axis4 = get(close_sp[3]),
           y = 1)) +  
  scale_x_discrete(limits = c(species,close_sp)) +
  geom_alluvium(aes(fill = Standardized_Name_spA)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  theme_minimal() +
  ggtitle(paste0("Cluster Markers of ", species, " shared by top 3 most similar species")) +
  scale_fill_discrete(name = species,guide = guide_legend(nrow = 3)) +
  theme(legend.position = "NA", aspect.ratio = 0.5, 
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

print(g)

ggsave(paste0(markers_dir, "/output/combined_markers/","combined_",file_names[i],".png"), g, width = 8, height = 5, units = c("in"))
}
```

![](scRNA_markers_files/figure-gfm/alluvial-1.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-2.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-3.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-4.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-5.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-6.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-7.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-8.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-9.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-10.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-11.png)<!-- -->![](scRNA_markers_files/figure-gfm/alluvial-12.png)<!-- -->

## 5.2 Generate ortholog marker gene tables for my focal species

``` r
for (i in 1:length(file_names)){
  species <- sub("_markers","",file_names[i])
  
  markers_pairs <- get(paste0("combined_",file_names[i]))
  
  markers_pairs %>% 
    filter(grepl("Montipora",protein_id_spB)) %>% 
    write.csv(paste0(markers_dir, "/output/Mcap/","Mcap_",species,"_Markers_pairs.csv"))
  
   markers_pairs %>% 
    filter(grepl("Porites",protein_id_spB)) %>% 
    write.csv(paste0(markers_dir, "/output/Pcomp/","Pcomp_",species,"_Markers_pairs.csv"))
  
  markers_pairs %>% 
    filter(grepl("Pocillopora_acuta",protein_id_spB)) %>% 
    write.csv(paste0(markers_dir, "/output/Pacuta/","Pacuta_",species,"_Markers_pairs.csv"))
}

Kevin_all <- combined_PdamRG_markers %>%
    left_join(Pdam_RG_filtered, by = join_by(cluster,gene_id_spA==gene)) %>% 
    filter(grepl("Montipora",protein_id_spB)|
             grepl("Porites",protein_id_spB) | 
           grepl("Pocillopora",protein_id_spB)) 
write.csv(Kevin_all,paste0(markers_dir, "/output/","Kevin_Pdam_combined_markers.csv"))
```
