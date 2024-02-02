


rm(list = ls())

library(magrittr)
library(ggplot2)

myplots <- list()

my_colo <- c("#2B2B2B", "#F58225") %>%
    set_names(c("Planktonic", "Pellicle"))

my_data_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Acetyl_formatted_per_sites.txt", "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Kentache-2016_acetylation_formatted_per_sites.txt")
pep_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Planktonic/Kentache_Acetylation_anticorps.txt"

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
    levels = paste(rep(levels(my_data_count$Condition), each = 2), rep(c("Proteins", "Sites"), times = 2)),
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

my_data_wide_venn <- my_data %>%
    dplyr::group_by(., `Accessions ABYAL`, Position, PTM, AA) %>%
    dplyr::summarise(., Condition = paste0(Condition, collapse = "__")) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., value = TRUE) %>%
    tidyr::pivot_wider(
        data = ., names_from = Condition,
        values_from = value, values_fill = FALSE)

my_data_wide_prot <- my_data %>%
    dplyr::select(., `Accessions ABYAL`, Condition, PTM) %>%
    dplyr::mutate(., value = TRUE) %>%
    unique(.) %>%
    tidyr::pivot_wider(
        data = ., names_from = Condition,
        values_from = value, values_fill = FALSE)

for (x in unique(my_data_wide$PTM)) {
    
    myplots[[paste0(x, "_venn")]] <- ggvenn::ggvenn(
        data = my_data_wide %>% dplyr::filter(., PTM == x),
        columns = levels(my_data$Condition), show_percentage = T,
        fill_color = unname(my_colo)) +
        ggtitle(x)
    
    myplots[[paste0(x, "_count")]] <- ggplot(
        my_data_count %>% dplyr::filter(., PTM == x),
        aes(x = Label, y = value, group = Condition, fill = Condition, label = value)) +
        geom_bar(stat = "identity", position = "dodge", colour = "black") +
        geom_text(stat = "identity", position = position_dodge(width = 0.9), vjust = -0.3) +
        ggpubr::theme_pubr() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "right") +
        xlab("Features") +
        ylab(paste0("Number of ", x, "ation")) +
        scale_fill_manual(values = my_colo)
    
    myplots[[paste0(x, "_prot_venn")]] <- ggvenn::ggvenn(
        data = my_data_wide_prot %>% dplyr::filter(., PTM == x),
        columns = levels(my_data$Condition), show_percentage = T,
        fill_color = unname(my_colo)) +
        ggtitle(x)
    
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
    
    data.table::fwrite(
        x = my_data_wide %>% dplyr::filter(., PTM == x),
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/Condition_", x, "_for_OA.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
    data.table::fwrite(
        x = my_data_wide_venn %>% dplyr::filter(., PTM == x),
        file = paste0("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/Venn_", x, "_for_OA.txt"),
        append = F, quote = F, sep = "\t", row.names = F, col.names = T)
    
}

data.table::fwrite(
    x = my_data_wide,
    file = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/Number_and_overlap_sites.txt",
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

#pdf("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/Number_and_overlap_sites.pdf", 8, 8)
#myplots
#dev.off()

cairo_pdf(filename = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/Number_and_overlap_sites_cairo.pdf", width = 8, height = 8, onefile = T)
for (x in names(myplots)) {
    print(myplots[[x]])
    #grid::grid.newpage()
}
dev.off()

save.image("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/Number_and_overlap_sites.RData")


