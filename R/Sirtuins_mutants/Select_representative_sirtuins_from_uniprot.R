


rm(list = ls())

library(magrittr)
library(ggplot2)

my_plot <- list()

my_annot_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25.tsv"

my_fasta_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25.fasta"

my_interpro_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25_interpro_simplified.txt"

my_outlier_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/Tree_2024-06-26/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25_selected_msa_fasttree.outliers"

my_annot <- data.table::fread(input = my_annot_f, sep = "\t", header = T)

my_lineage <- unique(my_annot$`Common taxon ID`) %>%
    split(., ceiling(seq_along(.)/100)) %>%
    lapply(., function(x) {
        tryCatch(
            taxize::classification(sci_id = x, db = "ncbi"),
            error = function(err) NA)
    })

my_lineage_df <- my_lineage %>%
    unlist(., recursive = FALSE) %>%
    plyr::ldply(., data.table::data.table, .id = "Common taxon ID") %>%
    dplyr::mutate(., `Common taxon ID` = sub(".+\\.", "", `Common taxon ID`)) %>%
    dplyr::select(., `Common taxon ID`, name, rank, id)

if (any(is.na(my_lineage_df$name))) {
    my_lineage_df <- unique(my_lineage_df[is.na(my_lineage_df$name),][["Common taxon ID"]]) %>%
        taxize::classification(sci_id = ., db = "ncbi") %>%
        plyr::ldply(., data.table::data.table, .id = "Common taxon ID") %>%
        dplyr::bind_rows(my_lineage_df[!is.na(my_lineage_df$name),], .)
}

my_lineage_df %<>%
    unique(.) %>%
    dplyr::mutate(., id = as.integer(id))

data.table::fwrite(
    x = my_lineage_df, file = sub(".tsv", "_taxon.txt", my_annot_f),
    append = F, quote = F, sep = "\t",
    row.names = F, col.names = T)

my_interpro <- data.table::fread(
    input = my_interpro_f, sep = "\t", header = T, quote = "")

my_sirtuin_ids <- my_interpro %>%
    dplyr::filter(
        ., grepl("(deacylase|deacetylase|SIR|Sirtuin|cobyrinic acid a)", `Signature description`) |
            grepl("deacylase|deacetylase|SIR|Sirtuin|cobyrinic acid a", `InterPro description`)) %>%
    .[["Protein accession"]] %>%
    unique(.)

length(unique(my_interpro$`Protein accession`))

my_ignore <- my_interpro %>%
    dplyr::filter(., !`Protein accession` %in% my_sirtuin_ids) %>%
    .[["Signature description"]] %>%
    table(.) %>%
    data.table::as.data.table(.) %>%
    dplyr::arrange(., dplyr::desc(N))

my_kept <- my_interpro %>%
    dplyr::filter(., `Protein accession` %in% my_sirtuin_ids) %>%
    .[["Signature description"]] %>%
    table(.) %>%
    data.table::as.data.table(.) %>%
    dplyr::arrange(., dplyr::desc(N))

my_interpro_gather <- my_interpro %>%
    dplyr::filter(., `Protein accession` %in% my_sirtuin_ids) %>%
    dplyr::select(
        ., `Protein accession`, `Signature description`,
        `InterPro description`) %>%
    tidyr::pivot_longer(
        data = ., cols = c(`Signature description`, `InterPro description`),
        names_to = "DB", values_to = "Domain") %>%
    dplyr::group_by(., `Protein accession`) %>%
    dplyr::summarise(
        ., Domain = paste0(sort(unique(Domain)), collapse = ";")) %>%
    dplyr::ungroup(.)

my_annot_interpro <- my_annot %>%
    dplyr::filter(., `Cluster ID` %in% my_sirtuin_ids) %>%
    dplyr::left_join(
    x = ., y = my_interpro_gather,
    by = c("Cluster ID" = "Protein accession")) %>%
    dplyr::filter(., !grepl(
        "active regulator", `Cluster Name`, ignore.case = T)) %>%
    dplyr::mutate(., Type = dplyr::case_when(
        grepl("NAD-dependent deacetylase, Class IV", Domain) | grepl("(SIRT|sirtuin-)(6|7)", `Cluster Name`) ~ "Sirtuin, class IV",
        grepl("Sirtuin, class III", Domain) ~ "Sirtuin, class III",
        grepl("Sirtuin, class II", Domain) ~ "Sirtuin, class II",
        grepl("Sirtuin, class I", Domain) ~ "Sirtuin, class I",
        grepl("Cobyrinic acid a,c-diamide synthase CbiA", Domain) ~ "CobB",
        grepl("Sirtuin, class U", Domain) ~ "Sirtuin, class U",
        TRUE ~ "Sirtuin family"
    ))

my_plot[["count_total"]] <- ggplot(my_annot_interpro, aes(x = Type)) +
    geom_bar() +
    ggpubr::theme_pubr() +
    coord_flip()

my_plot[["length_total"]] <- ggplot(my_annot_interpro, aes(x = Type, y = Length)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    coord_flip()

my_lineage_df_wide <- my_lineage_df %>%
    dplyr::select(., -id) %>%
    dplyr::filter(., rank %in% c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
    tidyr::pivot_wider(data = ., names_from = "rank", values_from = "name") %>%
    dplyr::mutate(., `Common taxon ID` = as.integer(`Common taxon ID`))

my_annot_interpro_tax <- my_lineage_df_wide %>%
    dplyr::left_join(x = my_annot_interpro, y = ., by = "Common taxon ID") %>%
    dplyr::group_by(., family, Type) %>%
    dplyr::mutate(., Priority = 1:dplyr::n()) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Priority = ifelse(!grepl("Cluster", my_annot_interpro$`Cluster Name`), 0, Priority))

my_outliers <- data.table::fread(input = my_outlier_f, sep = "\t", header = T)

my_annot_interpro_tax_filt <- my_annot_interpro_tax %>%
    dplyr::filter(., !`Cluster ID` %in% my_outliers[my_outliers$Distance > 6.5, ][["UniProtID"]]) %>%
    dplyr::filter(., is.na(Length) | (Length < 900 & Length > 200)) %>%
    dplyr::filter(., Priority <= 1)

my_plot[["count_selected"]] <- ggplot(my_annot_interpro_tax_filt, aes(x = Type)) +
    geom_bar() +
    ggpubr::theme_pubr() +
    coord_flip()

my_plot[["length_selected"]] <- ggplot(my_annot_interpro_tax_filt, aes(x = Type, y = Length)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    coord_flip()

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA", as.string = T)

my_fasta_filt <- my_fasta[
    which(names(my_fasta) %in% my_annot_interpro_tax_filt$`Cluster ID`)]

seqinr::write.fasta(
    sequences = my_fasta_filt, names = names(my_fasta_filt),
    file.out = sub(".fasta", "_selected2.fasta", my_fasta_f),
    open = "w", nbchar = 60)

data.table::fwrite(
    x = my_annot_interpro_tax_filt, file = sub(".fasta", "_selected2.txt", my_fasta_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

pdf(sub(".fasta", "_2.pdf", my_fasta_f), 10, 10)
my_plot
dev.off()

