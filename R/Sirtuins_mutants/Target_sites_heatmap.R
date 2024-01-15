


rm(list = ls())

library(magrittr)
library(ggplot2)

my_plot <- list()

my_colo <- c(
    "#5757f9ff", "#f94040ff", "#00c000ff", "#fdd61aff",
    "#6ecdcdff", "#a049beff") %>%
    set_names(c(
        "WT", "\u0394Sir2-Ab17",
        "\u0394CobB", "\u0394Sir2-Ab17\u0394CobB", 
        "Characteristic", "Function"))

#my_colo <- c(
#    "#5757f9ff", "#5757f9ff", "#f94040ff", "#f94040ff",
#    "#00c000ff", "#00c000ff", "#fdd61aff", "#fdd61aff",
#    "#6ecdcdff", "#6ecdcdff", "#6ecdcdff",
#    "#a049beff", "#a049beff", "#a049beff", "#a049beff") %>%
#    set_names(c(
#        "A_WT", "S_WT", "A_\u0394Sir2-Ab17", "S_\u0394Sir2-Ab17",
#        "A_\u0394CobB", "S_\u0394CobB",
#        "A_\u0394Sir2-Ab17\u0394CobB", "S_\u0394Sir2-Ab17\u0394CobB",
#        "In_actsite", "In_Operon", "In_Essential",
#        "In_biofilm", "In_motility", "In_adhesion", "In_resistance"))

my_data_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Acetyl_formatted_per_sites.txt", "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Succinyl_formatted_per_sites.txt")

my_annot_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Annotation/2023-05-16/Acinetobacter_baumannii_ATCC_17978_full_annotation_2023-05-16_manual_review_withOperon.txt"

my_domain_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/InterPro/PTM_on_activesite_2024-01-10.txt"

my_data <- lapply(my_data_f, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = T)
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = NULL) %>%
    tidyr::separate_rows(data = ., Condition, sep = ";")

if (any(grepl("^Pel", my_data$Condition))) {
    my_data %<>%
        dplyr::filter(., grepl("^Pel", Condition))
}

my_data %<>%
    dplyr::mutate(
        ., Condition = sub("_(Acet|Succi)_Rep.", "", Condition)) %>%
    dplyr::mutate(., Condition = dplyr::case_when(
        Condition %in% c("Pel_WT", "WT", "BF_WT") ~ names(my_colo)[1],
        Condition %in% c("Pel_dKDAC", "dKDAC (=Ab17Sir2)", "BF_dKDAC") ~ names(my_colo)[2],
        Condition %in% c("Pel_dNpdA", "dNpdA (=CobB)", "BF_dNpdA") ~ names(my_colo)[3],
        Condition %in% c("Pel_dNpdA_dKDAC", "dKDAC_NpdA (=double mutant)", "BF_dNpdA_dKDAC") ~ names(my_colo)[4])) %>%
    dplyr::select(., -Modifications, -`Modifications protein`, -Conserved, -ID, -`Accessions A1S`) %>%
    unique(.) %>%
    dplyr::mutate(., PTM = sub("^(.).+", "\\1", PTM))

my_data$Condition <- factor(
    x = my_data$Condition, levels = names(my_colo), ordered = TRUE)

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "", header = T)

my_annot_targets <- my_annot
for (w in c("biofilm", "motility", "adhesion", "resistance")) {
    my_annot_targets[[paste0("In_", w)]] <- apply(my_annot, MARGIN = 1, function(x) {
        any(grepl(w, x, ignore.case = T))
    })
    
}

my_annot_targets$In_Operon <- !is.na(my_annot_targets$OperonID)

my_annot_targets$In_Essential <- !is.na(my_annot_targets$essentiality) & my_annot_targets$essentiality == "E"

my_domain <- data.table::fread(
    input = my_domain_f, sep = "\t", quote = "", header = T) %>%
    dplyr::filter(., grepl("Robin", Study))

my_domain_sum <- my_domain %>%
    dplyr::group_by(., `Accessions ABYAL`, Position) %>%
    dplyr::summarise(., `locs.description` = paste0(unique(`locs.description`), collapse = ";")) %>%
    dplyr::ungroup(.)

interesting_sites <- my_data %>%
    dplyr::left_join(x = ., y = my_domain_sum)

interesting_sites$In_actsite <- !is.na(interesting_sites$`locs.description`) & interesting_sites$`locs.description` != ""

interesting_sites <- my_annot_targets %>%
    dplyr::select(
        ., `Locus Tag`, `Gene Name`, `GenBank_gene_accession`,
        tidyselect::starts_with("In_")) %>%
    dplyr::left_join(
        x = interesting_sites,
        y = .,
        by = c("Accessions ABYAL" = "Locus Tag")) %>%
    dplyr::select(., -`locs.description`)

interesting_sites %<>%
    dplyr::mutate(
        ., Name = dplyr::case_when(
            !is.na(`Gene Name`) & `Gene Name` != "" ~ stringr::str_to_title(`Gene Name`),
            !is.na(`GenBank_gene_accession`) & `GenBank_gene_accession` != "" ~ `GenBank_gene_accession`,
            TRUE ~ `Accessions ABYAL`),
        Modified = paste0(AA, Position)) %>%
    tidyr::unite(data = ., col = "Site", Name, Modified, sep = "_", remove = F) %>%
    tidyr::unite(data = ., col = "In_condition", PTM, Condition, sep = "_")

interesting_sites_long <- interesting_sites %>%
    dplyr::mutate(., value = TRUE) %>%
    tidyr::pivot_wider(data = ., names_from = In_condition, values_from = value, values_fill = FALSE) %>%
    tidyr::pivot_longer(data = ., cols = tidyselect::starts_with(c("In_", "A_", "S_"))) %>%
    unique(.) %>%
    dplyr::mutate(., Type = sub("^(A|S)_", "", name)) %>%
    dplyr::mutate(
        ., Type = dplyr::case_when(
            !grepl("^In_", Type) ~ Type,
            Type %in% c("In_biofilm", "In_motility", "In_adhesion", "In_resistance") ~ "Function",
            TRUE ~ "Characteristic")) %>%
    dplyr::mutate(
        ., Type = ifelse(value == TRUE, Type, NA)) %>%
#    dplyr::filter(., value == TRUE) %>%
    dplyr::arrange(., Name, Position)

interesting_sites_long$Site <- factor(
    x = interesting_sites_long$Site,
    levels = rev(unique(interesting_sites_long$Site)),
    ordered = TRUE)

interesting_sites_long$name <- factor(
    x = interesting_sites_long$name,
    levels = c(
        "A_WT", "A_ΔSir2-Ab17", "A_ΔCobB", "A_ΔSir2-Ab17ΔCobB", "S_WT", "S_ΔSir2-Ab17", "S_ΔCobB", "S_ΔSir2-Ab17ΔCobB",
        "In_Operon", "In_Essential", "In_actsite", "In_biofilm", "In_motility", "In_adhesion", "In_resistance"),
    ordered = TRUE)

interesting_sites_long$Type <- factor(
    x = interesting_sites_long$Type,
    levels = c(
        "WT", "ΔSir2-Ab17", "ΔCobB", "ΔSir2-Ab17ΔCobB", "Characteristic", "Function", NA),
    ordered = TRUE)

my_targets <- my_annot_targets %>%
    dplyr::filter(., In_biofilm == TRUE | In_motility == TRUE | In_adhesion == TRUE | In_resistance == TRUE) %>%
    .[["Locus Tag"]] %>%
    unique(.)

toplot <- interesting_sites_long %>%
    dplyr::filter(., `Accessions ABYAL` %in% my_targets)

toplot_panel <- split(x = unique(toplot$Name), ceiling(seq_along(unique(toplot$Name))/14)) %>%
    plyr::ldply(., data.table::data.table, .id = "Panel") %>%
    set_colnames(c("Panel", "Name"))

toplot %<>%
    dplyr::left_join(x = ., y = toplot_panel, by = "Name")

my_plot[["heatmap_biofilm"]] <- ggplot(
    toplot,
    aes(x = name, y = Site, fill = Type)) +
    geom_tile(colour = "grey") +
    ggpubr::theme_pubr() +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = my_colo, na.value = "#F2F2F2") +
    facet_wrap(facets = vars(Panel), ncol = 3, scales = "free_y")

my_targets <- interesting_sites_long %>%
    dplyr::filter(., name == "In_actsite" & value == TRUE) %>%
    .[["Site"]] %>%
    unique(.)

toplot <- interesting_sites_long %>%
    dplyr::filter(., Site %in% my_targets)

toplot_panel <- split(x = unique(toplot$Name), ceiling(seq_along(unique(toplot$Name))/27)) %>%
    plyr::ldply(., data.table::data.table, .id = "Panel") %>%
    set_colnames(c("Panel", "Name"))

toplot %<>%
    dplyr::left_join(x = ., y = toplot_panel, by = "Name")

my_plot[["heatmap_actsite"]] <- ggplot(
    toplot,
    aes(x = name, y = Site, fill = Type)) +
    geom_tile(colour = "grey") +
    ggpubr::theme_pubr() +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = my_colo, na.value = "#F2F2F2") +
    facet_wrap(facets = vars(Panel), ncol = 4, scales = "free_y")

pdf("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sites_of_interest/Site_heatmap.pdf", width = 10, height = 9)
my_plot
dev.off()


