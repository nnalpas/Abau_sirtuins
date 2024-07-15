


rm(list = ls())

library(magrittr)
library(ggplot2)
library(Easy)

myplots <- list()

my_colo <- c("#5757f9ff", "#f94040ff", "#00c000ff", "#fdd61aff") %>%
    set_names(c("WT", "\u0394Sir2-Ab17", "\u0394CobB", "\u0394Sir2-Ab17\u0394CobB"))

my_data_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Acetyl_formatted_per_sites.txt", "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Succinyl_formatted_per_sites.txt")

my_fasta_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/ATCC_17978_plasmides.fasta")

my_annot_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Annotation/2023-05-16/Acinetobacter_baumannii_ATCC_17978_full_annotation_2024-03-12_manual_review_withOperon_withSir-targets.txt"

my_sig <- "D:/PXD004167/combined/txt/Statistics_2023-07-05/Statistics_report.RData"

my_actsites_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/InterPro/PTM_on_activesite_2024-01-10.txt"

my_domain_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/InterPro/PTM_on_domain_2024-01-10.txt"

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
    dplyr::rename(., Sample = Condition) %>%
    dplyr::mutate(
        ., Condition = sub("_(Acet|Succi)_Rep.", "", Sample)) %>%
    dplyr::mutate(., Condition = dplyr::case_when(
        Condition %in% c("Pel_WT", "WT", "BF_WT") ~ names(my_colo)[1],
        Condition %in% c("Pel_dKDAC", "dKDAC (=Ab17Sir2)", "BF_dKDAC") ~ names(my_colo)[2],
        Condition %in% c("Pel_dNpdA", "dNpdA (=CobB)", "BF_dNpdA") ~ names(my_colo)[3],
        Condition %in% c("Pel_dNpdA_dKDAC", "dKDAC_NpdA (=double mutant)", "BF_dNpdA_dKDAC") ~ names(my_colo)[4])) %>%
    dplyr::select(., -Modifications, -`Modifications protein`, -Conserved, -ID, -`Accessions A1S`) %>%
    unique(.) %>%
    dplyr::mutate(
        ., PTM = sub("^(.).+", "\\1", PTM),
        Replicate = sub("Pel_", "", Sample) %>% sub("(Acet|Succi)_", "", .))

my_data$Condition <- factor(
    x = my_data$Condition, levels = names(my_colo), ordered = TRUE)


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
    input = my_annot_f, sep = "\t", header = T)

my_annot %<>%
    dplyr::select(
        ., `Locus Tag`, `GenBank_accession`, `GenBank_gene_accession`,
        `Gene Name [combined]`, `Product Name`,
        `KEGG Pathway Name`, `GOBP Term`,
        `GOCC Term`, `GOMF Term`, COG_function,
        `Subcellular Localization [b2g]`, COG_function_review,
        OperonID, essentiality, `Known sirtuin target`, `RefSeq ID`) %>%
    dplyr::left_join(
        x = my_k_count,
        y = ., by = c("Accessions ABYAL" = "Locus Tag"))


my_session <- new.env()
load(my_sig, envir = my_session)
my_stats <- my_session$my_stats_final %>%
    dplyr::filter(
        ., `ANOVA BH corrected p-value` <= 0.05 &
            !is.na(`Tukey adj.p.value`) & `Tukey adj.p.value` <= 0.05) %>%
    dplyr::select(
        ., `Protein IDs`, `ANOVA BH corrected p-value`,
        contrast, `fold change`) %>%
    tidyr::pivot_wider(
        data = ., names_from = "contrast", values_from = "fold change") %>%
    tidyr::separate_rows(data = ., `Protein IDs`, sep = ";")

my_annot %<>%
    dplyr::left_join(
        x = .,
        y = my_stats, by = c("RefSeq ID" = "Protein IDs"))


my_domain <- data.table::fread(
    input = my_domain_f, sep = "\t", quote = "", header = T) %>%
    dplyr::filter(., grepl("Robin", Study))

my_domain_sum <- my_domain %>%
    dplyr::filter(
        ., !is.na(`signature.description`) & `signature.description` != "") %>%
    dplyr::group_by(., `Accessions ABYAL`, Position) %>%
    dplyr::summarise(
        ., `Functional domain` = paste0(
            unique(`signature.description`), collapse = ";")) %>%
    dplyr::ungroup(.)


my_actsite <- data.table::fread(
    input = my_actsites_f, sep = "\t", quote = "", header = T) %>%
    dplyr::filter(., grepl("Robin", Study))

my_actsite_sum <- my_actsite %>%
    dplyr::group_by(., `Accessions ABYAL`, Position) %>%
    dplyr::summarise(., `Active sites` = paste0(unique(`locs.description`), collapse = ";")) %>%
    dplyr::ungroup(.)


for (x in unique(my_data$PTM)) {
    
    my_sites_explained <- my_data %>%
        dplyr::filter(., PTM == x) %>%
        dplyr::left_join(x = ., y = my_windows) %>%
        tidyr::unite(data = ., col = "Modification", Position, PTM, sep = "") %>%
        dplyr::group_by(., `Accessions ABYAL`, Modification, `Sequence window`) %>%
        dplyr::summarise(
            ., Conditions = paste0(sort(unique(Condition)), collapse = ";"),
            Replicates = paste0(sort(unique(Replicate)), collapse = ";")) %>%
        #dplyr::mutate(., Comment = dplyr::case_when(
        #    Conditions == "ΔSir2-Ab17;ΔSir2-Ab17ΔCobB" ~ "Only by Ab17Sir2",
        #    Conditions == "ΔSir2-Ab17;ΔCobB;ΔSir2-Ab17ΔCobB" ~ "By Ab17Sir2 & CobB",
        #    Conditions == "ΔCobB;ΔSir2-Ab17ΔCobB" ~ "Only by CobB",
        #    Conditions == "WT;ΔSir2-Ab17;ΔCobB;ΔSir2-Ab17ΔCobB" ~ "Common",
        #    Conditions == "WT" ~ "By other enzyme or chemical",
        #    TRUE ~ "Unconclusive"
        #))
        dplyr::mutate(., Comment = dplyr::case_when(
            Conditions %in% c("ΔSir2-Ab17", "ΔSir2-Ab17;ΔSir2-Ab17ΔCobB") ~ "Only by Ab17Sir2",
            Conditions %in% c("ΔSir2-Ab17ΔCobB", "ΔCobB;ΔSir2-Ab17;ΔSir2-Ab17ΔCobB", "ΔSir2-Ab17;ΔCobB;ΔSir2-Ab17ΔCobB") ~ "By Ab17Sir2 & CobB",
            Conditions %in% c("ΔCobB", "ΔCobB;ΔSir2-Ab17ΔCobB") ~ "Only by CobB",
            Conditions %in% c("WT;ΔCobB;ΔSir2-Ab17;ΔSir2-Ab17ΔCobB", "WT;ΔSir2-Ab17;ΔCobB;ΔSir2-Ab17ΔCobB") ~ "Common",
            TRUE ~ "Unconclusive"
        ))
    
    my_data_explained <- my_sites_explained %>%
        dplyr::group_by(., `Accessions ABYAL`, Comment) %>%
        dplyr::summarise(
            ., Modification = paste0(sort(unique(Modification)), collapse = ";")) %>%
        tidyr::pivot_wider(
            data = ., names_from = Comment, values_from = Modification)
    
    my_sites_explained %<>%
        tidyr::pivot_longer(data = ., cols = c(Conditions, Comment), values_to = "Explanation") %>%
        dplyr::mutate(., value = TRUE)
    
    my_prot_explained <- my_sites_explained %>%
        dplyr::group_by(., `Accessions ABYAL`, name) %>%
        dplyr::summarise(
            ., Nmb_modification = dplyr::n_distinct(Modification),
            Explanation = paste0(unique(Explanation), collapse = ";")) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(
            ., Explanation = ifelse(
                grepl(";", Explanation), "Multiple explanation", Explanation),
            value = TRUE) %>%
        tidyr::unite(
            data = ., col = "Comment", name, Explanation, sep = ": ") %>%
        dplyr::arrange(., Comment) %>%
        tidyr::pivot_wider(
            data = ., names_from = "Comment", values_from = "value")
    my_prot_explained[is.na(my_prot_explained)] <- FALSE
    
    my_sites_explained %<>%
        tidyr::separate_rows(data = ., Explanation, sep = ";") %>%
        tidyr::unite(
            data = ., col = "Comment", name, Explanation, sep = ": ") %>%
        dplyr::arrange(., Comment) %>%
        tidyr::pivot_wider(
            data = ., names_from = "Comment", values_from = "value")
    my_sites_explained[is.na(my_sites_explained)] <- FALSE
    
    my_sites_explained_final <- my_sites_explained %>%
        dplyr::mutate(., Position = as.integer(sub("[A-Z]+", "", Modification))) %>%
        dplyr::left_join(x = ., y = my_domain_sum) %>%
        dplyr::left_join(x = ., y = my_actsite_sum) %>%
        dplyr::left_join(x = ., y = my_annot[, c(
            "Accessions ABYAL", "Gene Name [combined]",
            "GenBank_accession", "GenBank_gene_accession", "RefSeq ID")])
    
    my_data_wide <- my_data %>%
        dplyr::filter(., PTM == x) %>%
        tidyr::unite(
            data = ., col = "Modification", Position, PTM, sep = "") %>%
        dplyr::group_by(., `Accessions ABYAL`, Condition) %>%
        dplyr::summarise(
            ., Modification = paste0(sort(Modification), collapse = ";")) %>%
        tidyr::pivot_wider(
            data = ., names_from = Condition, values_from = Modification) %>%
        dplyr::left_join(x = ., y = my_data_explained) %>%
        dplyr::left_join(x = ., y = my_annot)
    
    data.table::fwrite(
        x = my_data_wide,
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Condition_explanation/", x, "_explain_2024-07-15.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
    data.table::fwrite(
        x = my_sites_explained_final,
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Condition_explanation/", x, "_sites_explain_2024-07-15.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
    data.table::fwrite(
        x = my_sites_explained,
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Condition_explanation/", x, "_windows_for_OA_2024-07-15.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
    data.table::fwrite(
        x = my_prot_explained,
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Condition_explanation/", x, "_proteins_for_OA_2024-07-15.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
}

my_data_wide <- my_data %>%
    dplyr::group_by(., `Accessions ABYAL`, Position, Condition) %>%
    dplyr::summarise(
        ., PTM = paste0(sort(unique(PTM)), collapse = "+")) %>%
    tidyr::unite(data = ., col = "Modification", Position, PTM, sep = "") %>%
    dplyr::group_by(., `Accessions ABYAL`, Condition) %>%
    dplyr::summarise(
        ., Modification = paste0(sort(Modification), collapse = ";")) %>%
    tidyr::pivot_wider(
        data = ., names_from = Condition, values_from = Modification) %>%
    dplyr::left_join(x = ., y = my_annot)

my_domain_sum_gather <- my_domain_sum %>%
    dplyr::mutate(., `Functional domain` = paste(Position, " [", `Functional domain`, "]", sep = "")) %>%
    dplyr::group_by(., `Accessions ABYAL`) %>%
    dplyr::summarise(., `Functional domain` = paste0(`Functional domain`, collapse = ";")) %>%
    dplyr::ungroup(.)

my_actsite_sum_gather <- my_actsite_sum %>%
    dplyr::mutate(., `Active sites` = paste(Position, " [", `Active sites`, "]", sep = "")) %>%
    dplyr::group_by(., `Accessions ABYAL`) %>%
    dplyr::summarise(., `Active sites` = paste0(`Active sites`, collapse = ";")) %>%
    dplyr::ungroup(.)

my_data_wide %<>%
    dplyr::left_join(x = ., y = my_domain_sum_gather) %>%
    dplyr::left_join(x = ., y = my_actsite_sum_gather)

data.table::fwrite(
    x = my_data_wide,
    file = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Condition_explanation/Acetylation_Succinylation_sirtuins_mutants_2024-07-15.txt",
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)


