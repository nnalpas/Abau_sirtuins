


rm(list = ls())

library(magrittr)
library(ggplot2)

myplots <- list()

my_colo <- c("#5757f9ff", "#f94040ff", "#00c000ff", "#fdd61aff") %>%
    set_names(c("WT", "\u0394Sir2-Ab17", "\u0394CobB", "\u0394Sir2-Ab17\u0394CobB"))

my_data_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Acetyl_formatted_per_sites.txt", "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Succinyl_formatted_per_sites.txt")

#my_data_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin-2023_sirtuins_formatted_per_sites.txt")

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
    dplyr::mutate(., Condition = sub("_(Acet|Succi)_Rep.", "", Condition)) %>%
    dplyr::mutate(., Condition = dplyr::case_when(
        Condition %in% c("Pel_WT", "WT", "BF_WT") ~ names(my_colo)[1],
        Condition %in% c("Pel_dKDAC", "dKDAC (=Ab17Sir2)", "BF_dKDAC") ~ names(my_colo)[2],
        Condition %in% c("Pel_dNpdA", "dNpdA (=CobB)", "BF_dNpdA") ~ names(my_colo)[3],
        Condition %in% c("Pel_dNpdA_dKDAC", "dKDAC_NpdA (=double mutant)", "BF_dNpdA_dKDAC") ~ names(my_colo)[4]
    )) %>%
    unique(.)

my_data$Condition <- factor(x = my_data$Condition, levels = names(my_colo), ordered = TRUE)

my_data_wide <- my_data %>%
    dplyr::mutate(., value = TRUE) %>%
    tidyr::pivot_wider(
        data = ., names_from = Condition,
        values_from = value, values_fill = FALSE)

my_data_count <- my_data %>%
    dplyr::select(., `Accessions ABYAL`, Position, PTM, Condition) %>%
    unique(.) %>%
    dplyr::group_by(., PTM, Condition) %>%
    dplyr::summarise(
        ., Sites = dplyr::n(),
        Proteins = dplyr::n_distinct(`Accessions ABYAL`)) %>%
    dplyr::ungroup(.) %>%
    tidyr::pivot_longer(data = ., cols = c(Sites, Proteins)) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(., Label = paste(Condition, name))

my_data_count$Label <- factor(
    x = my_data_count$Label,
    levels = paste(rep(levels(my_data_count$Condition), each = 2), rep(c("Proteins", "Sites"), times = 4)),
    ordered = T)

my_data_multi <- my_data %>%
    dplyr::select(., `Accessions ABYAL`, Position, PTM, Condition) %>%
    unique(.) %>%
    dplyr::group_by(., PTM, Condition, `Accessions ABYAL`) %>%
    dplyr::summarise(
        ., Sites = dplyr::n()) %>%
    dplyr::mutate(., Sites = ifelse(Sites > 4, "5 and more", Sites)) %>%
    dplyr::group_by(., PTM, Condition, Sites) %>%
    dplyr::summarise(
        ., Count = dplyr::n_distinct(`Accessions ABYAL`)) %>%
    dplyr::mutate(., Percentage = Count * 100 / sum(Count)) %>%
    dplyr::ungroup(.)

my_data_pep <- lapply(sub("per_sites", "peptides", my_data_f), function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = T)
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = NULL) %>%
    tidyr::separate_rows(data = ., Condition, sep = ";")

if (any(grepl("^Pel", my_data_pep$Condition))) {
    my_data_pep %<>%
        dplyr::filter(., grepl("^Pel", Condition))
}

my_data_pep %<>%
    dplyr::mutate(., Condition = sub("_(Acet|Succi)_Rep.", "", Condition)) %>%
    dplyr::mutate(., Condition = dplyr::case_when(
        Condition %in% c("Pel_WT", "WT") ~ names(my_colo)[1],
        Condition %in% c("Pel_dKDAC", "dKDAC (=Ab17Sir2)") ~ names(my_colo)[2],
        Condition %in% c("Pel_dNpdA", "dNpdA (=CobB)") ~ names(my_colo)[3],
        Condition %in% c("Pel_dNpdA_dKDAC", "dKDAC_NpdA (=double mutant)") ~ names(my_colo)[4]
    )) %>%
    unique(.)

my_data_pep_count <- my_data_pep %>%
    dplyr::select(., `Accessions ABYAL`, Position, PTM, Condition, Sequences) %>%
    unique(.) %>%
    dplyr::group_by(., PTM, Condition) %>%
    dplyr::summarise(
        ., Redundant = dplyr::n(),
        `Non-redundant` = dplyr::n_distinct(toupper(Sequences))) %>%
    dplyr::ungroup(.) %>%
    tidyr::pivot_longer(data = ., cols = c(Redundant, `Non-redundant`)) %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(., Label = paste(Condition, name))

my_data_pep_count$Label <- factor(
    x = my_data_pep_count$Label,
    levels = paste(rep(levels(my_data_count$Condition), each = 2), rep(c("Redundant", "Non-redundant"), times = 4)),
    ordered = T)

for (x in unique(my_data_wide$PTM)) {
    
    myplots[[paste0(x, "_venn")]] <- ggvenn::ggvenn(
        data = my_data_wide %>% dplyr::filter(., PTM == x),
        columns = levels(my_data$Condition), show_percentage = T,
        fill_color = unname(my_colo)) +
        ggtitle(x)
    
    myplots[[paste0(x, "_count")]] <- ggplot(
        my_data_count %>% dplyr::filter(., PTM == x),
        aes(x = paste(Condition, name), y = value, group = Condition, fill = Condition, label = value)) +
        geom_bar(stat = "identity", position = "dodge", colour = "black") +
        geom_text(stat = "identity", position = position_dodge(width = 0.9), vjust = -0.3) +
        ggpubr::theme_pubr() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "right") +
        xlab("Features") +
        ylab(paste0("Number of ", x, "ation")) +
        scale_fill_manual(values = my_colo)
    
    myplots[[paste0(x, "_multi")]] <- ggplot(
        my_data_multi %>% dplyr::filter(., PTM == x),
        aes(x = Sites, y = Percentage, fill = Condition, label = Percentage)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black", width = 0.8) +
        #geom_text(stat = "identity", position = position_dodge(width = 0.9), vjust = -0.3) +
        ggpubr::theme_pubr() +
        #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "right") +
        xlab(paste0("K-", x, "ated sites per protein")) +
        ylab(paste0(x, "ated proteins")) +
        scale_fill_manual(values = my_colo)
    
    myplots[[paste0(x, "_pep_count")]] <- ggplot(
        my_data_pep_count %>% dplyr::filter(., PTM == x),
        aes(x = paste(Condition, name), y = value, group = Condition, fill = Condition, label = value)) +
        geom_bar(stat = "identity", position = "dodge", colour = "black") +
        geom_text(stat = "identity", position = position_dodge(width = 0.9), vjust = -0.3) +
        ggpubr::theme_pubr() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "right") +
        xlab("Features") +
        ylab(paste0("Number of ", x, "ated peptides")) +
        scale_fill_manual(values = my_colo)
    
}

data.table::fwrite(
    x = my_data_wide,
    file = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Identified_sites/Pel_BF_Overlap_detected_sites.txt",
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

#pdf("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Identified_sites/Number_and_overlap_sites.pdf", 8, 8)
#myplots
#dev.off()

cairo_pdf(filename = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Identified_sites/Pel_BF_Number_and_overlap_sites_cairo.pdf", width = 8, height = 8, onefile = T)
for (x in names(myplots)) {
    print(myplots[[x]])
    #grid::grid.newpage()
}
dev.off()

save.image("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Identified_sites/Pel_BF_Number_and_overlap_sites.RData")


