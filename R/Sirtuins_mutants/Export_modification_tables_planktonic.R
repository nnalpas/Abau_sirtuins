


rm(list = ls())

library(magrittr)
library(ggplot2)
library(Easy)

myplots <- list()

my_colo <- c("#2B2B2B", "#F58225") %>%
    set_names(c("Planktonic", "Pellicle"))

my_data_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Acetyl_formatted_per_sites.txt", "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Kentache-2016_acetylation_formatted_per_sites.txt")
pep_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Planktonic/Kentache_Acetylation_anticorps.txt"

my_fasta_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/ATCC_17978_plasmides.fasta")

my_annot_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Annotation/2023-05-16/Acinetobacter_baumannii_ATCC_17978_full_annotation_2023-05-16.txt"

my_data <- lapply(my_data_f, function(x) {
    my_res <- data.table::fread(
        input = x, sep = "\t", quote = "", header = T)
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = NULL) %>%
    tidyr::separate_rows(data = ., Condition, sep = ";") %>%
    dplyr::mutate(
        ., Condition = ifelse(is.na(Condition), "Planktonic", Condition),
        AA = ifelse(is.na(AA), "K", AA))

my_pep <- data.table::fread(
    input = pep_f,
    sep = "\t", quote = "", header = T) %>%
    .[["Sequences"]]

my_data %<>%
    tidyr::separate_rows(data = ., Sequences, sep = ";") %>%
    dplyr::filter(., is.na(Sequences) | toupper(Sequences) %in% toupper(my_pep))

if (any(grepl("^(Pel_WT|Planktonic)", my_data$Condition))) {
    my_data %<>%
        dplyr::filter(., grepl("^(Pel_WT|Planktonic)", Condition))
}

my_data %<>%
    dplyr::mutate(., Condition = sub("^Pel_.*", "Pellicle", Condition)) %>%
    dplyr::select(., -Modifications, -`Modifications protein`, -Conserved, -ID, -`Accessions A1S`, -Sequences) %>%
    unique(.) %>%
    dplyr::mutate(., PTM = sub("^(.).+", "\\1", PTM))

my_data$Condition <- factor(x = my_data$Condition, levels = names(my_colo), ordered = TRUE)


my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA",
    as.string = T, whole.header = F)
names(my_fasta) <- sub("\\|.*", "", names(my_fasta))

my_k_count <- data.table::data.table(
    `Accessions ABYAL` = names(my_fasta),
    Nmb_Lysine = stringr::str_count(string = my_fasta, pattern = "K"))

my_windows <- my_data %>%
    dplyr::select(., `Accessions ABYAL`, Position) %>%
    unique(.) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(
        ., `Sequence window` = site_sequence_window(
            position = Position,
            protein = as.character(my_fasta[[`Accessions ABYAL`]])))

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "", header = T)

my_annot %<>%
    dplyr::select(
        ., `Locus Tag`, `GenBank_accession`, `GenBank_gene_accession`,
        `Gene Name`, `Product Name`,
        `KEGG Pathway Name`, `GOBP Term`,
        `GOCC Term`, `GOMF Term`, COG_function,
        `Subcellular Localization [b2g]`) %>%
    dplyr::left_join(
        x = my_k_count,
        y = ., by = c("Accessions ABYAL" = "Locus Tag"))

for (x in unique(my_data$PTM)) {
    
    my_sites_explained <- my_data %>%
        dplyr::filter(., PTM == x) %>%
        dplyr::left_join(x = ., y = my_windows) %>%
        tidyr::unite(data = ., col = "Modification", Position, PTM, sep = "")
    
    my_sites_explained %<>%
        dplyr::mutate(., value = TRUE)
    
    my_prot_explained <- my_sites_explained %>%
        dplyr::group_by(., `Accessions ABYAL`, Condition) %>%
        dplyr::summarise(
            ., Nmb_modification = dplyr::n_distinct(Modification),
            value = TRUE) %>%
        dplyr::ungroup(.) %>%
        tidyr::pivot_wider(
            data = ., names_from = "Condition", values_from = "value")
    my_prot_explained[is.na(my_prot_explained)] <- FALSE
    
    my_sites_explained %<>%
        tidyr::pivot_wider(
            data = ., names_from = "Condition", values_from = "value")
    my_sites_explained[is.na(my_sites_explained)] <- FALSE
    
    my_data_wide <- my_data %>%
        dplyr::filter(., PTM == x) %>%
        tidyr::unite(
            data = ., col = "Modification", Position, PTM, sep = "") %>%
        dplyr::group_by(., `Accessions ABYAL`, Condition) %>%
        dplyr::summarise(
            ., Modification = paste0(sort(Modification), collapse = ";")) %>%
        tidyr::pivot_wider(
            data = ., names_from = Condition, values_from = Modification) %>%
        dplyr::left_join(x = ., y = my_annot)
    
    data.table::fwrite(
        x = my_data_wide,
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/", x, "_explain.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
    data.table::fwrite(
        x = my_sites_explained,
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/", x, "_windows_for_OA.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
    data.table::fwrite(
        x = my_prot_explained,
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/", x, "_proteins_for_OA.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
}


