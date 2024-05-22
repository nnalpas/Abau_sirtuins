


rm(list = ls())

library(magrittr)
library(ggplot2)
library(phangorn)

my_big_cols <- RColorBrewer::brewer.pal(n = 12, name = "Set3") %>%
    sample(x = ., size = length(.), replace = F)

my_plots <- list()

my_tree_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/EggNOG/COG0846_tree.nwk"

my_gene_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/EggNOG/COG0846_nog_orthologs.txt"

my_taxid_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/EggNOG/e5.taxid_info.tsv"



my_tree <- ape::read.tree(file = my_tree_f)

plot(my_tree, "unrooted")


my_gene <- data.table::fread(
    input = my_gene_f, sep = "\t", header = F,
    col.names = c("EggNOG_ID", "Gene_name", "Taxon_name", "Taxon_ID", "Gene_synonym"))

my_taxid <- data.table::fread(
    input = my_taxid_f, sep = "\t", header = T)

table(my_gene$Gene_name[nchar(my_gene$Gene_name) < 5])

tmp <- grep("(sir|srt)", my_gene$Gene_name, ignore.case = T, value = T)
tmp[nchar(tmp)>=5]
unique(tmp)

my_gene %<>%
    dplyr::mutate(., Label = dplyr::case_when(
        grepl("cobb", Gene_name, ignore.case = T) ~ "CobB",
        grepl("npd", Gene_name, ignore.case = T) ~ "NpdA",
        grepl("sirt", Gene_name, ignore.case = T) ~ sub("sirt([0-9]).*", "SIRT\\1", Gene_name, ignore.case = T),
        grepl("(sir2|sir-2)", Gene_name, ignore.case = T) ~ "Sir2",
        TRUE ~ NA_character_
    ))

my_gene[my_gene$EggNOG_ID == "470.IX87_02480", "Label"] <- "Abau CobB"
my_gene[my_gene$EggNOG_ID == "1217715.F994_02012", "Label"] <- "Abau Sir2"

table(my_gene$Label)

tocolour <- my_gene %>%
    dplyr::filter(., !is.na(Label)) %>%
    dplyr::select(., Label) %>%
    unique(.) %>%
    dplyr::bind_cols(., Colour = my_big_cols[1:nrow(.)])

my_ids <- my_tree$tip.label %>%
    data.table::data.table(EggNOG_ID = .) %>%
    tidyr::separate(
        data = ., col = EggNOG_ID, into = c("taxid", "gene"),
        sep = "\\.", remove = FALSE)

colors <- my_gene %>%
    dplyr::filter(., EggNOG_ID %in% my_ids$EggNOG_ID) %>%
    dplyr::left_join(
        x = ., y = tocolour,
        by = "Label") %>%
    dplyr::mutate(
        ., Colour = tidyr::replace_na(
            data = Colour, replace = "#404040")) %>%
    dplyr::arrange(., match(EggNOG_ID, my_ids$EggNOG_ID))

tipcol <- colors[["Colour"]]

legend <- colors %>%
    dplyr::select(., Label, Colour)
legend %<>%
    unique(.)


plot(x = my_tree, main = "UPGMA", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
legend("bottomright", legend = legend[["Label"]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)

plot(x = my_tree, "unrooted", main = "NJ", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.4)
legend("bottomright", legend = legend[["Label"]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)



