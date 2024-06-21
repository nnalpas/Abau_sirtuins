


rm(list = ls())

library(magrittr)
library(ggplot2)
library(phangorn)

my_big_cols <- RColorBrewer::brewer.pal(n = 12, name = "Paired") %>%
    sample(x = ., size = length(.), replace = F)

my_plots <- list()

#my_tree_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf/All_consurf_sequences_msa.tree"
#my_tree_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf/PhyloSuite/FastTree_results/2024_05_28-12_08_20/all_gene_trees.nwk"
my_tree_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf_optimised/TrimAl_FastTree/All_consurf_sequences_downloaded_2024-06-12_14_msa_fasttree_noTrim.tree"
#my_tree_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf_optimised/TrimAl_FastTree/All_consurf_sequences_downloaded_2024-06-12_14_msa_fasttree.tree"
my_tree_f <- "C:/Users/nalpanic/Downloads/All_consurf_sequences_downloaded_2024-06-12_14_msa.tree"

#my_gene_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf_optimised/All_consurf_sequences_header.txt"
my_gene_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf_optimised/All_consurf_sequences_header_2024-06-12_14.txt"



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

my_name <- sub("\\.txt", "_genename.txt", my_gene_f) %>%
    data.table::fread(
        input = ., sep = "\t", header = T, select = c("From", "Gene Names")) %>%
    set_colnames(c("UniProtID", "GeneName"))

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
    dplyr::filter(., UniProtID %in% sub("UniRef90_", "", my_ids$UniProtID)) %>%
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
    dplyr::arrange(., dplyr::desc(Count)) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        Label = dplyr::case_when(
            Rank < 5 ~ name,
            TRUE ~ "Others"
        )) %>%
    dplyr::group_by(., Type, Label) %>%
    dplyr::summarise(
        ., name = paste0(name, collapse = "; "),
        Count = sum(Count)) %>%
    dplyr::group_by(., Type) %>%
    dplyr::summarise(
        ., name = name, Count = Count, Percent = Count/sum(Count)*100) %>%
    dplyr::ungroup(.)




pdf(sub("\\.(nwk|tree)", "_tree.pdf", my_tree_f), 10, 10)

plot(x = my_tree, main = "UPGMA", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
legend("bottomright", legend = legend[["Type"]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)
add.scale.bar(cex = 1, font = 2, col = "black")

plot(x = my_tree, "unrooted", main = "NJ", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
legend("bottomright", legend = legend[["Type"]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)
add.scale.bar(cex = 1, font = 2, col = "black")

my_distances <- ape::dist.nodes(x = my_tree)

dev.off()


pdf(sub("\\.(nwk|tree)", "_taxon.pdf", my_tree_f), 10, 10)

color_shades <- list(
    CobB = c("#052018", "#105f47", "#1b9e77", "#76c5ad", "#d1ece4"),
    CobB_long = c("#2e081c", "#8b1953", "#e7298a", "#f17fb9", "#fad4e8"),
    SrtN = c("#171624", "#46436b", "#7570b3", "#aca9d1", "#e3e2f0"),
    Sir2 = c("#2b1300", "#823901", "#d95f02", "#e89f67", "#f7dfcc"))

my_lineage_stats$Type <- sub("_optimised", "", my_lineage_stats$Type)

for (x in unique(my_lineage_stats$Type)) {
    
    toplot <- my_lineage_stats %>%
        dplyr::filter(., Type == x) %>%
        dplyr::arrange(., dplyr::desc(Count))
    
    toplot$name <- factor(x = toplot$name, levels = toplot$name, ordered = T)
    
    pl <- ggplot(
        toplot,
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
            axis.ticks = element_blank()) +
        scale_fill_manual(values = color_shades[[x]])
    
    print(pl)
    
}

dev.off()



my_gene_lineage_complete <- my_taxid_format %>%
    dplyr::left_join(x = ., y = my_lineage_format, by = c()) %>%
    #dplyr::filter(., rank == "species") %>%
    dplyr::select(., Organism, name, id, rank) %>%
    unique(.) %>%
    dplyr::left_join(x = my_gene, y = .)

priority_species <- c(
    "Acinetobacter baumannii", "Enterobacterales", "Salmonella Typhi",
    "Shigella spp.", "Mycobacterium tuberculosis", "Enterococcus faecium",
    "Pseudomonas aeruginosa", "Salmonella", "Neisseria gonorrhoeae",
    "Staphylococcus aureus", "Streptococci", "Streptococcus pneumoniae",
    "Haemophilus influenzae")
priority_genus <- c(
    "Acinetobacter", "Enterobacterales", "Salmonella",
    "Shigella", "Mycobacterium", "Enterococcus",
    "Pseudomonas", "Salmonella", "Neisseria",
    "Staphylococcus", "Streptococci", "Streptococcus",
    "Haemophilus")

my_gene_lineage_complete_filt <- my_gene_lineage_complete %>%
    dplyr::filter(., grepl(
        paste("^(", paste0(priority_genus, collapse = "|"), ")", sep = ""), name, ignore.case = T))

View(unique(my_gene_lineage_complete_filt[, c("Type", "name")]))

data.table::fwrite(
    x = my_gene_lineage_complete_filt,
    file = sub("\\.(nwk|tree)", "_pathogens.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

my_kingdom <- my_gene_lineage_complete %>%
    dplyr::filter(., rank == "superkingdom") %>%
    dplyr::select(., Type, name) %>%
    unique(.)

save.image(file = sub("\\.(nwk|tree)", "_session.RData", my_tree_f))


