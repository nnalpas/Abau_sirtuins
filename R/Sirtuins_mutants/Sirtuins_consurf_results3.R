


rm(list = ls())

library(magrittr)
library(ggplot2)

my_plots <- list()

my_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf_optimised/Sir2_optimised_Consurf_Outputs_1716881395/msa_fasta.aln"

my_fasta <- seqinr::read.fasta(file = my_f, seqtype = "AA", as.string = T)

my_df <- names(my_fasta) %>%
    data.table::data.table(Header = .) %>%
    dplyr::filter(., Header != "Input_seq") %>%
    dplyr::mutate(
        ., Evalue = sub("(.+)\\|.+?$", "\\1", Header) %>%
            sub(".+\\|", "", .) %>% as.numeric(.),
        Taxon = sub(".+\\|", "", Header))

my_taxon_hierar <- taxize::classification(my_df$Taxon, db = "ncbi") %>%
    plyr::ldply(., data.table::data.table, .id = "Taxon") %>%
    dplyr::select(., -`V1`)

my_taxon_hierar_wide <- my_taxon_hierar %>%
    dplyr::filter(
        ., !is.na(Taxon) & !is.na(name) &
            rank %in% c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
    dplyr::select(., Taxon, name, rank) %>%
    unique(.) %>%
    tidyr::pivot_wider(data = ., names_from = "rank", values_from = "name") %>%
    dplyr::left_join(x = my_df, y = .)

for (x in c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) {
    toplot <- my_taxon_hierar_wide %>%
        dplyr::group_by(., !!as.name(x)) %>%
        dplyr::summarise(., Count = dplyr::n()) %>%
        dplyr::ungroup(.)
    my_plots[[paste0("hist_", x)]] <- ggplot(toplot, aes(x = !!as.name(x), y = Count, fill = !!as.name(x), label = Count)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(stat = "identity", position = position_dodge(width = 0.9), hjust = -0.3) +
        ggpubr::theme_pubr() +
        coord_flip() +
        theme(legend.position = "none")
    my_plots[[paste0("box_", x)]] <- ggplot(my_taxon_hierar_wide, aes(x = !!as.name(x), y = -log(Evalue, 10), colour = !!as.name(x))) +
        geom_boxplot(linewidth = 1)+
        geom_jitter(position = position_jitter(width = 0.3), size = 1.5, alpha = 0.5) +
        ggpubr::theme_pubr() +
        coord_flip() +
        theme(legend.position = "none")
}

toplot <- my_taxon_hierar_wide %>%
    dplyr::group_by(., superkingdom, phylum, class, order) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(phylum)

my_plots[["treemap"]] <- ggplot(
    toplot, aes(area = Count, fill = phylum, label = order)) +#, subgroup = phylum)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white", reflow = T)# +
    #treemapify::geom_treemap_subgroup_border(colour = "red", size = 5) +
    #treemapify::geom_treemap_subgroup_text(
    #   place = "middle",
    #    colour = "red",
    #    alpha = 0.5,
    #    grow = T
    #)

my_plots[["box_phylum_class"]] <- ggplot(
    my_taxon_hierar_wide %>%
        dplyr::filter(., phylum %in% c("Pseudomonadota", "Bacteroidota")) %>%
        dplyr::filter(., !is.na(class)),
    aes(x = class, y = -log(Evalue, 10), colour = phylum)) +
    geom_boxplot(linewidth = 1)+
    geom_jitter(position = position_jitter(width = 0.3), size = 1.5, alpha = 0.5) +
    ggpubr::theme_pubr() +
    coord_flip() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(phylum), scales = "free_y", space = "free_y")

my_plots[["box_class_order"]] <- ggplot(
    my_taxon_hierar_wide %>%
        dplyr::filter(., class %in% c("Gammaproteobacteria", "Betaproteobacteria", "Bacteroidia")) %>%
        dplyr::filter(., !is.na(order)),
    aes(x = order, y = -log(Evalue, 10), colour = class)) +
    geom_boxplot(linewidth = 1)+
    geom_jitter(position = position_jitter(width = 0.3), size = 1.5, alpha = 0.5) +
    ggpubr::theme_pubr() +
    coord_flip() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(class), scales = "free_y", space = "free_y")

my_plots[["box_order_family"]] <- ggplot(
    my_taxon_hierar_wide %>%
        dplyr::filter(., order %in% c("Moraxellales", "Vibrionales", "Enterobacterales", "Burkohlderiales", "Bacteroidales")) %>%
        dplyr::filter(., !is.na(family)),
    aes(x = family, y = -log(Evalue, 10), colour = order)) +
    geom_boxplot(linewidth = 1)+
    geom_jitter(position = position_jitter(width = 0.3), size = 1.5, alpha = 0.5) +
    ggpubr::theme_pubr() +
    coord_flip() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(order), scales = "free_y", space = "free_y")

my_plots[["box_family_genus"]] <- ggplot(
    my_taxon_hierar_wide %>%
        dplyr::filter(., family %in% c("Moraxellaceae", "Vibrionaceae", "Erwiniaceae", "Enterobacteriaceae")) %>%
        dplyr::filter(., !is.na(genus)),
    aes(x = genus, y = -log(Evalue, 10), colour = family)) +
    geom_boxplot(linewidth = 1)+
    geom_jitter(position = position_jitter(width = 0.3), size = 1.5, alpha = 0.5) +
    ggpubr::theme_pubr() +
    coord_flip() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(family), scales = "free_y", space = "free_y")

pl <- cowplot::plot_grid(
    my_plots$box_phylum, my_plots$box_phylum_class,
    my_plots$box_class_order, my_plots$box_order_family,
    nrow = 1, ncol = 4)

pdf(sub("\\.aln", "_analysis.pdf", my_f), 12, 12)
my_plots
pl
dev.off()

pdf(sub("\\.aln", "_analysis_combined.pdf", my_f), width = 20, height = 10)
pl
dev.off()


