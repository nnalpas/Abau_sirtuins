


rm(list = ls())

library(magrittr)
library(ggplot2)
library(phangorn)

my_big_cols <- RColorBrewer::brewer.pal(n = 12, name = "Paired") %>%
    sample(x = ., size = length(.), replace = F) %>%
    c(., c("#29FF3B", "#FF00BF", "#0000FA"))

my_plots <- list()

my_tree_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/Tree_2024-06-27/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25_selected2_msa_fasttree.tree"

my_gene_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25_selected2.txt"

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
    input = my_gene_f, sep = "\t", header = T) %>%
    dplyr::rename(., UniProtID = `Cluster ID`) %>%
    dplyr::mutate(., UniProtID = sub(":.*", "", UniProtID))

other_type <- my_gene %>%
    dplyr::filter(., Type == "Sirtuin family") %>%
    .[["Cluster Name"]] %>% 
    table(.) %>%
    data.table::as.data.table(.) %>%
    dplyr::arrange(., dplyr::desc(N))

my_gene$Other_type <- ifelse(
    my_gene$`Cluster Name` %in% other_type[1:10][[1]],
    my_gene$`Cluster Name`, NA)

my_gene$Target_type <- ifelse(
    my_gene$`Cluster Name` %in% unique(my_gene$`Cluster Name`)[1:13],
    my_gene$`Cluster Name`, NA)

label_col <- "Target_type"

table(my_gene[[label_col]])

tocolour <- my_gene %>%
    dplyr::filter(., !is.na(!!as.name(label_col))) %>%
    dplyr::select(., !!as.name(label_col)) %>%
    unique(.) %>%
    dplyr::bind_cols(., Colour = my_big_cols[1:nrow(.)])

colors <- my_gene %>%
    dplyr::filter(., UniProtID %in% my_ids$UniProtID) %>%
    dplyr::left_join(
        x = ., y = tocolour,
        by = label_col) %>%
    #dplyr::mutate(
    #    ., Colour = tidyr::replace_na(
    #        data = Colour, replace = "#404040")) %>%
    dplyr::arrange(., match(UniProtID, my_ids$UniProtID))

tipcol <- colors[["Colour"]]

legend <- colors %>%
    dplyr::select(., !!as.name(label_col), Colour)
legend %<>%
    unique(.)

my_lineage_stats <- my_gene %>%
    dplyr::filter(., !is.na(phylum)) %>%
    dplyr::group_by(., !!as.name(label_col), phylum) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::group_by(., !!as.name(label_col)) %>%
    dplyr::arrange(., dplyr::desc(Count), .by_group = T) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        Label = dplyr::case_when(
            Rank < 5 ~ phylum,
            TRUE ~ "Others"
        )) %>%
    dplyr::group_by(., !!as.name(label_col), Label) %>%
    dplyr::summarise(
        ., name = paste0(phylum, collapse = "; "),
        Count = sum(Count)) %>%
    dplyr::group_by(., !!as.name(label_col)) %>%
    dplyr::summarise(
        ., name = name, Count = Count, Percent = Count/sum(Count)*100) %>%
    dplyr::ungroup(.)

pdf(sub("\\.(nwk|tree)", "_tree.pdf", my_tree_f), 10, 10)

plot(x = my_tree, main = "UPGMA", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
legend("bottomright", legend = legend[[label_col]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)
add.scale.bar(cex = 1, font = 2, col = "black")

plot(x = my_tree, "unrooted", main = "NJ", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
i <- which(my_tree$tip.label %in% my_gene$UniProtID[!is.na(my_gene$Target_type)])
tiplabels(my_tree$tip.label[i], i, col = tipcol, adj = 0.5, frame = "none")
legend("bottomright", legend = legend[[label_col]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)
add.scale.bar(cex = 1, font = 2, col = "black")

my_distances <- ape::dist.nodes(x = my_tree)

my_dist_to_root <- adephylo::distRoot(x = my_tree)

data.table::fwrite(
    x = data.table::data.table(UniProtID = names(my_dist_to_root), Distance = my_dist_to_root),
    file = sub(".tree$", ".outliers", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

plot(x = my_tree, "unrooted", main = "NJ", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
i <- which(my_tree$tip.label %in% names(my_dist_to_root)[my_dist_to_root > 6.5])
tiplabels(my_tree$tip.label[i], i, col = tipcol, adj = 0.5, frame = "none")
legend("bottomright", legend = legend[[label_col]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)
add.scale.bar(cex = 1, font = 2, col = "black")


dev.off()


pdf(sub("\\.(nwk|tree)", "_taxon.pdf", my_tree_f), 10, 10)

color_shades <- list(
    CobB = c("#052018", "#105f47", "#1b9e77", "#76c5ad", "#d1ece4"),
    CobB_long = c("#2e081c", "#8b1953", "#e7298a", "#f17fb9", "#fad4e8"),
    SrtN = c("#171624", "#46436b", "#7570b3", "#aca9d1", "#e3e2f0"),
    Sir2 = c("#2b1300", "#823901", "#d95f02", "#e89f67", "#f7dfcc"))

my_lineage_stats[[label_col]] <- sub("_optimised", "", my_lineage_stats[[label_col]])

for (x in unique(my_lineage_stats[[label_col]])) {
    
    toplot <- my_lineage_stats %>%
        dplyr::filter(., !!as.name(label_col) == x) %>%
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

View(unique(my_gene_lineage_complete_filt[, c(label_col, "name")]))

data.table::fwrite(
    x = my_gene_lineage_complete_filt,
    file = sub("\\.(nwk|tree)", "_pathogens.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

my_kingdom <- my_gene_lineage_complete %>%
    dplyr::filter(., rank == "superkingdom") %>%
    dplyr::select(., !!as.name(label_col), name) %>%
    unique(.)

save.image(file = sub("\\.(nwk|tree)", "_session.RData", my_tree_f))


