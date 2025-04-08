


rm(list = ls())

library(magrittr)
library(ggplot2)
library(phangorn)

my_big_cols <- RColorBrewer::brewer.pal(n = 12, name = "Paired") %>%
    sample(x = ., size = length(.), replace = F) %>%
    c(., c("#29FF3B", "#FF00BF", "#0000FA"))

my_plots <- list()

my_files <- c(
    `SIR2` = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/SIRim/Rblast_Sirtuins_Bonhomme_SIR2/Best_Reciproc_Blast_cross-map_Sirtuins_Bonhomme_SIR2_recip",
    `SIRim` = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/SIRim/Rblast_Sirtuins_Bonhomme_SIRim/Best_Reciproc_Blast_cross-map_Sirtuins_Bonhomme_SIRim_recip")

my_data <- lapply(my_files, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "",
        header = T, colClasses = "character")
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "Bonhomme_2025") %>%
    dplyr::group_by(., qseqid) %>%
    dplyr::summarise(., Bonhomme_2025 = paste0(
        sort(unique(Bonhomme_2025)), collapse = ";")) %>%
    dplyr::ungroup(.)

my_gene_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25_selected2.txt"

my_gene <- data.table::fread(
    input = my_gene_f, sep = "\t", header = T) %>%
    dplyr::left_join(x = ., y = my_data, by = c("Cluster ID" = "qseqid"))

data.table::fwrite(
    x = my_gene, my_gene_f, append = F, quote = F,
    sep = "\t", row.names = F, col.names = T)


