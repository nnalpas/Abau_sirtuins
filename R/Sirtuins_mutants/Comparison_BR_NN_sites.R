


rm(list = ls())

library(magrittr)

myplots <- list()

my_br_old_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin-2023_sirtuins_formatted_per_sites.txt"

my_br_new_f <- c("C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Acetyl_formatted_per_sites.txt", "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/Robin_et_al._Sirtuins_mutants_Succinyl_formatted_per_sites.txt")

my_br_old <- data.table::fread(
    input = my_br_old_f, sep = "\t", quote = "", header = T) %>%
    dplyr::mutate(., Analysis = "Old")

my_br_new <- lapply(my_br_new_f, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = T)
    }) %>%
    plyr::ldply(., dplyr::bind_rows, .id = NULL) %>%
    dplyr::mutate(., Analysis = "New") %>%
    dplyr::filter(., grepl("Pel", Condition))

my_data <- dplyr::bind_rows(my_br_old, my_br_new)

my_data_no_cond_check <- my_data %>%
    dplyr::select(., `Accessions ABYAL`, Position, PTM, `Accessions A1S`, Analysis) %>%
    unique(.) %>%
    dplyr::mutate(., value = TRUE) %>%
    tidyr::pivot_wider(
        data = ., names_from = Analysis, values_from = value, values_fill = FALSE)

my_data_no_cond_check[!my_data_no_cond_check$Old | !my_data_no_cond_check$New, ]

myplots[["All_venn"]] <- ggvenn::ggvenn(
    data = my_data_no_cond_check,
    columns = c("Old", "New"), show_percentage = T) +
    ggplot2::ggtitle("All sites")

my_data_cond_check <- my_data %>%
    dplyr::select(., `Accessions ABYAL`, Position, PTM, `Accessions A1S`, Analysis, Condition) %>%
    tidyr::separate_rows(data = ., Condition, sep = ";") %>%
    dplyr::mutate(
        ., Condition = dplyr::case_when(
            Condition == "dKDAC_NpdA (=double mutant)" ~ "Pel_dNpdA_dKDAC",
            Condition == "dNpdA (=CobB)" ~ "Pel_dNpdA",
            Condition == "dKDAC (=Ab17Sir2)" ~ "Pel_dKDAC",
            Condition == "WT" ~ "Pel_WT",
            TRUE ~ sub("_(Succi|Acet)_Rep.$", "", Condition)
        ),
        value = TRUE) %>%
    unique(.) %>%
    tidyr::pivot_wider(
        data = ., names_from = Analysis, values_from = value, values_fill = FALSE) %>%
    tidyr::unite(data = ., col = Name, PTM, Condition, sep = "_")

for (x in unique(my_data_cond_check$Name)) {
    
    my_data_cond_check[my_data_cond_check$Name == x & (!my_data_cond_check$Old | !my_data_cond_check$New), ]
    
    myplots[[paste0(x, "_venn")]] <- ggvenn::ggvenn(
        data = my_data_cond_check %>% dplyr::filter(., Name == x),
        columns = c("Old", "New"), show_percentage = T) +
        ggplot2::ggtitle(x)
    
}

data.table::fwrite(
    x = my_data_no_cond_check,
    file = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparison_BR_NN/All_detected_sites.txt",
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

data.table::fwrite(
    x = my_data_cond_check,
    file = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparison_BR_NN/Per_condition_detected_sites.txt",
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

pdf("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparison_BR_NN/Detected_sites.pdf", 8, 8)
myplots
dev.off()


