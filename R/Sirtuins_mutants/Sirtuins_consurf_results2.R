


rm(list = ls())

library(magrittr)
library(ggplot2)
library(phangorn)

my_big_cols <- RColorBrewer::brewer.pal(n = 12, name = "Paired") %>%
    sample(x = ., size = length(.), replace = F) %>%
    c(., c("#29FF3B", "#FF00BF", "#0000FA"))

my_plots <- list()



### Customise the phylogenetic tree visualisation ------------------------

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

label_col <- "Type"

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

my_dist_to_root <- adephylo::distRoot(x = my_tree) %>%
    data.table::data.table(UniProtID = names(.), Distance = .)

my_dist_to_all <- cophenetic(x = my_tree) %>%
    data.table::as.data.table(.) %>%
    set_colnames(my_tree$tip.label) %>%
    dplyr::mutate(., Tip1 = my_tree$tip.label) %>%
    tidyr::pivot_longer(
        data = ., cols = -`Tip1`, names_to = "Tip2",
        values_to = "Cophenetic_Distances") %>%
    dplyr::filter(., `Tip1` != `Tip2`)

data.table::fwrite(
    x = my_dist_to_root,
    file = sub(".tree$", "_distance_to_root.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

data.table::fwrite(
    x = my_dist_to_all,
    file = sub(".tree$", "_distance_to_all.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

pdf(sub("\\.(nwk|tree)", "_tree.pdf", my_tree_f), 10, 10)

plot(x = my_tree, main = "UPGMA", tip.color = tipcol, show.tip.label = F, show.node.label = F)
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

plot(x = my_tree, "unrooted", main = "NJ", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
i <- which(my_tree$tip.label %in% my_dist_to_root[my_dist_to_root$Distance > 6.5][["UniProtID"]])
tiplabels(my_tree$tip.label[i], i, col = tipcol, adj = 0.5, frame = "none")
legend("bottomright", legend = legend[[label_col]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)
add.scale.bar(cex = 1, font = 2, col = "black")

dev.off()



### Get some statistics on the different cluster -------------------------

my_cluster_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/Tree_2024-06-27/Clusters.txt"

my_cluster <- data.table::fread(input = my_cluster_f, sep = "\t", header = T) %>%
    dplyr::mutate(., UniProtID = sub(" ", "_", UniProtID))

my_dist_clust <- my_dist_to_all %>%
    dplyr::left_join(x = ., y = my_cluster, by = c("Tip1" = "UniProtID")) %>%
    dplyr::rename(., Cluster1 = Cluster) %>%
    dplyr::left_join(x = ., y = my_cluster, by = c("Tip2" = "UniProtID")) %>%
    dplyr::rename(., Cluster2 = Cluster)

my_dist_per_clust <- my_dist_clust %>%
    dplyr::filter(., Cluster1 == Cluster2) %>%
    dplyr::group_by(., Cluster1, Cluster2) %>%
    dplyr::summarise(
        ., Mean_distance = mean(Cophenetic_Distances, na.rm = T),
        Max_distance = max(Cophenetic_Distances, na.rm = T))

data.table::fwrite(
    x = my_dist_per_clust, file = sub(".tree$", "_distance_per_cluster.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

my_gene_clust <- my_gene %>%
    dplyr::left_join(x = ., y = my_cluster)

tax_col <- "phylum"

cluster_col <- "Cluster"

my_phylum_stats <- my_gene_clust %>%
    dplyr::filter(., !is.na(!!as.name(tax_col))) %>%
    dplyr::group_by(., !!as.name(cluster_col), !!as.name(tax_col)) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::group_by(., !!as.name(cluster_col)) %>%
    dplyr::arrange(., dplyr::desc(Count), .by_group = T) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        Label = dplyr::case_when(
            Rank < 5 ~ !!as.name(tax_col),
            TRUE ~ "Others"
        )) %>%
    dplyr::group_by(., !!as.name(cluster_col), Label) %>%
    dplyr::summarise(
        ., name = paste0(!!as.name(tax_col), collapse = "; "),
        Count = sum(Count)) %>%
    dplyr::group_by(., !!as.name(cluster_col)) %>%
    dplyr::summarise(
        ., name = name, Count = Count, Percent = Count/sum(Count)*100) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Level = tax_col)

tax_col <- "superkingdom"

my_kingdom_stats <- my_gene_clust %>%
    dplyr::filter(., !is.na(!!as.name(tax_col)) & !!as.name(tax_col) != "") %>%
    dplyr::group_by(., !!as.name(cluster_col), !!as.name(tax_col)) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::group_by(., !!as.name(cluster_col)) %>%
    dplyr::arrange(., dplyr::desc(Count), .by_group = T) %>%
    dplyr::mutate(
        ., Rank = 1:dplyr::n(),
        Label = dplyr::case_when(
            Rank < 5 ~ !!as.name(tax_col),
            TRUE ~ "Others"
        )) %>%
    dplyr::group_by(., !!as.name(cluster_col), Label) %>%
    dplyr::summarise(
        ., name = paste0(!!as.name(tax_col), collapse = "; "),
        Count = sum(Count)) %>%
    dplyr::group_by(., !!as.name(cluster_col)) %>%
    dplyr::summarise(
        ., name = name, Count = Count, Percent = Count/sum(Count)*100) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Level = tax_col)

my_lineage_stats <- dplyr::bind_rows(my_kingdom_stats, my_phylum_stats)

data.table::fwrite(
    x = my_lineage_stats, file = sub(".tree$", "_lineage_stats.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

pdf(sub("\\.(nwk|tree)", "_taxon.pdf", my_tree_f), 10, 10)

color_shades <- list(
    `Class I` = c("#231208", "#6a3518", "#b15928", "#c17a53", "#e0bda9"),
    `Class II` = c("#231f08", "#6a5c18", "#b19928", "#c1ad53", "#e0d6a9"),
    `Class III` = c("#00231f", "#00695e", "#00af9c", "#33bfb0", "#99dfd7"),
    `Class IV` = c("#0a1c05", "#1e5510", "#328d1a", "#5ba448", "#add1a3"),
    `Class U` = c("#000b1c", "#002255", "#00388d", "#3360a4", "#99afd1"),
    `Class TM` = c("#210819", "#62194a", "#a3297b", "#b55495", "#daa9ca"),
    `Unknown` = c("#1a1a1a", "#333333", "#4d4d4d", "#808080", "#b3b3b3"))

for (x in unique(my_phylum_stats[[cluster_col]])) {
    
    toplot <- my_phylum_stats %>%
        dplyr::filter(., !!as.name(cluster_col) == x) %>%
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

#priority_species <- c(
#    "Acinetobacter baumannii", "Enterobacterales", "Salmonella Typhi",
#    "Shigella spp.", "Mycobacterium tuberculosis", "Enterococcus faecium",
#    "Pseudomonas aeruginosa", "Salmonella", "Neisseria gonorrhoeae",
#    "Staphylococcus aureus", "Streptococci", "Streptococcus pneumoniae",
#    "Haemophilus influenzae")
#priority_genus <- c(
#    "Acinetobacter", "Enterobacterales", "Salmonella",
#    "Shigella", "Mycobacterium", "Enterococcus",
#    "Pseudomonas", "Salmonella", "Neisseria",
#    "Staphylococcus", "Streptococci", "Streptococcus",
#    "Haemophilus")

my_domain_stats <- my_gene_clust %>%
    dplyr::select(., UniProtID, Cluster, Domain) %>%
    tidyr::separate_rows(data = ., Domain, sep = ";") %>%
    dplyr::filter(., Domain != "" & Domain != "-") %>%
    dplyr::group_by(., Cluster) %>%
    dplyr::mutate(., Total = dplyr::n_distinct(UniProtID)) %>%
    dplyr::group_by(., Cluster, Domain) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(UniProtID), Total = unique(Total)) %>%
    dplyr::mutate(., Percentage = Count * 100 / Total) %>%
    dplyr::arrange(., dplyr::desc(Count), .by_group = T) %>%
    dplyr::slice(., 1:20)

my_domain_cluster <- my_domain_stats %>%
    dplyr::group_by(., Domain) %>%
    dplyr::summarise_all(~paste0(., collapse = ";")) %>%
    dplyr::ungroup(.)

data.table::fwrite(
    x = my_domain_stats,
    file = sub("\\.(nwk|tree)", "_domain_stats.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

data.table::fwrite(
    x = my_domain_cluster,
    file = sub("\\.(nwk|tree)", "_domain_cluster.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

data.table::fwrite(
    x = my_gene_clust,
    file = sub("\\.(nwk|tree)", "_annotation.txt", my_tree_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

save.image(file = sub("\\.(nwk|tree)", "_session.RData", my_tree_f))


