


rm(list = ls())

library(magrittr)
library(ggplot2)

my_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/NCBI_Blast/ABYAL1514-NCBI-RefProt-Alignment_hit.txt"

my_data <- data.table::fread(
    input = my_f, sep = "\t", header = T, fill = T) %>%
    dplyr::mutate(., `Query cover` = as.numeric(sub("%", "", `Query cover`)))

pdf(sub(".txt", ".pdf", my_f), 10, 10)

ggplot(my_data, aes(x = `Perc. Ident`)) +
    geom_histogram(binwidth = 3, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 3, vjust = -0.5) +
    ylab("Number of sequence")

ggplot(my_data, aes(x = `Query cover`)) +
    geom_histogram(binwidth = 5, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 5, vjust = -0.5) +
    ylab("Number of sequence")

ggplot(my_data, aes(x = `E-Value`)) +
    geom_histogram(binwidth = 0.001, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 0.001, vjust = -0.5) +
    ylab("Number of sequence")

ggplot(my_data, aes(x = `Total Score`)) +
    geom_histogram(binwidth = 20, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 20, vjust = -0.5) +
    ylab("Number of isolates")

my_data_filt <- my_data %>%
    dplyr::filter(., `E-Value` < 0.0001 & `Perc. Ident` > 35 & `Query cover` > 25 & !(`Perc. Ident` > 95 & `Query cover` > 90))

ggplot(my_data_filt, aes(x = `Perc. Ident`)) +
    geom_histogram(binwidth = 3, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 3, vjust = -0.5) +
    ylab("Number of sequence")

ggplot(my_data_filt, aes(x = `Query cover`)) +
    geom_histogram(binwidth = 5, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 5, vjust = -0.5) +
    ylab("Number of sequence")

my_lineage <- taxize::classification(unique(my_data_filt$Taxid), db = "ncbi") %>%
    plyr::ldply(., data.table::data.table, .id = "Taxid")

unique(my_lineage[my_lineage$rank == "genus", ][["name"]])
unique(my_lineage[my_lineage$rank == "family", ][["name"]])
unique(my_lineage[my_lineage$rank == "order", ][["name"]])
unique(my_lineage[my_lineage$rank == "class", ][["name"]])
unique(my_lineage[my_lineage$rank == "phylum", ][["name"]])
unique(my_lineage[my_lineage$rank == "superkingdom", ][["name"]])




dev.off()


