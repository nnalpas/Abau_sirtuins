


rm(list = ls())

library(magrittr)
library(ggplot2)
library(phangorn)

my_big_cols <- RColorBrewer::brewer.pal(n = 4, name = "Dark2") %>%
    sample(x = ., size = length(.), replace = F)

my_plots <- list()

my_tree_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf/All_consurf_sequences_msa.tree"
my_tree_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf/PhyloSuite/FastTree_results/2024_05_28-12_08_20/all_gene_trees.nwk"

my_gene_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf/All_consurf_sequences_header.txt"



my_tree <- ape::read.tree(file = my_tree_f)
plot(my_tree, "unrooted", show.tip.label = F)

if (all(grepl("_", my_tree$tip.label))) {
    my_ids <- my_tree$tip.label %>%
        data.table::data.table(ID = .) %>%
        tidyr::separate(
            data = ., col = ID, into = c("Number", "UniProtID"),
            sep = "_", remove = FALSE, extra = "merge")
} else {
    my_ids <- data.table::data.table(UniProtID = my_tree$tip.label)
}



my_gene <- data.table::fread(
    input = my_gene_f, sep = "\t", header = F, col.names = c("Header")) %>%
    tidyr::separate(data = ., col = "Header", into = c("Type", "UniProtID", "Match", "E-value", "Organism"), sep = "\\|") %>%
    dplyr::mutate(., Type = sub("\\..*", "", Type)) %>%
    tidyr::separate(data = ., col = "Match", into = c("Start", "End"), sep = "_")

#data.table::fwrite(
#    x = my_gene, file = sub("\\.txt", "_formatted.txt", my_gene_f),
#    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

my_name <- paste(
    dirname(my_gene_f),
    "All_consurf_sequences_header_genename.txt", sep = "/") %>%
    data.table::fread(
        input = ., sep = "\t", header = T,
        col.names = c("UniProtID", "GeneName"))

my_gene %<>%
    dplyr::left_join(x = ., y = my_name, by = "UniProtID") %>%
    dplyr::mutate(., Label = dplyr::case_when(
        !is.na(GeneName) & nchar(GeneName) < 7 ~ GeneName,
        is.na(GeneName) & is.na(UniProtID) ~ Type,
        TRUE ~ UniProtID)) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(., UniProtID = dplyr::case_when(
            is.na(UniProtID) ~ paste0(c(Type, "_Input_seq"), collapse = ""),
            TRUE ~ UniProtID)
    ) %>%
    dplyr::ungroup(.)

table(my_gene$Label)

tocolour <- my_gene %>%
    dplyr::filter(., !is.na(Type)) %>%
    dplyr::select(., Type) %>%
    unique(.) %>%
    dplyr::bind_cols(., Colour = my_big_cols[1:nrow(.)])

colors <- my_gene %>%
    dplyr::filter(., UniProtID %in% my_ids$UniProtID) %>%
    dplyr::left_join(
        x = ., y = tocolour,
        by = "Type") %>%
    dplyr::mutate(
        ., Colour = tidyr::replace_na(
            data = Colour, replace = "#404040")) %>%
    dplyr::arrange(., match(UniProtID, my_ids$UniProtID))

tipcol <- colors[["Colour"]]

legend <- colors %>%
    dplyr::select(., Type, Colour)
legend %<>%
    unique(.)

plot(x = my_tree, main = "UPGMA", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
legend("bottomright", legend = legend[["Type"]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)

plot(x = my_tree, "unrooted", main = "NJ", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.4)
legend("bottomright", legend = legend[["Type"]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)





my_taxid <- my_gene$Organism %>%
    unique(.) %>%
    taxize::get_ids(sci_com = ., db = "ncbi")

my_taxid_format <- data.table::data.table(
    ID = as.integer(my_taxid[[1]]), Organism = names(my_taxid[[1]]))

my_lineage <- unique(my_taxid_format$ID) %>%
    taxize::classification(sci_id = ., db = "ncbi")

my_lineage_format <- my_lineage %>%
    plyr::ldply(., data.table::data.table, .id = "ID") %>%
    dplyr::filter(., !is.na(ID)) %>%
    dplyr::select(., -`V1`) %>%
    dplyr::mutate(., ID = as.integer(as.character(ID)))

my_gene_lineage <- my_taxid_format %>%
    dplyr::left_join(x = ., y = my_lineage_format, by = c()) %>%
    dplyr::filter(., rank == "class") %>%
    dplyr::select(., Organism, name, id) %>%
    unique(.) %>%
    dplyr::left_join(x = my_gene, y = .)

my_lineage_stats <- my_gene_lineage %>%
    dplyr::filter(., !is.na(name)) %>%
    dplyr::group_by(., Type, name) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::group_by(., Type) %>%
    dplyr::summarise(
        ., name = name, Count = Count, Percent = Count/sum(Count)*100) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., Type, dplyr::desc(Percent))

my_lineage_stats_format <- my_lineage_stats %>%
    dplyr::


for ()

ggplot(
    my_lineage_stats %>% dplyr::filter(., Type == "Sir2" & Count > 1),
    aes(x = "", y = Percent, fill = name)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_text(
        aes(label = paste0(format(Percent, digits = 2), "%")),
        position = position_stack(vjust = 0.5)) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_classic() +
    theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
