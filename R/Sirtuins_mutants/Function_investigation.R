


rm(list = ls())

library(magrittr)
library(ggplot2)

myplots <- list()

my_colo <- c("#5757f9ff", "#f94040ff", "#00c000ff", "#fdd61aff") %>%
    set_names(c("WT", "\u0394Sir2-Ab17", "\u0394CobB", "\u0394Sir2-Ab17\u0394CobB"))

my_annot_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Annotation/2023-05-16/Acinetobacter_baumannii_ATCC_17978_full_annotation_2023-05-16.txt"

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
    dplyr::mutate(
        ., Replicate = sub(".+_(Rep.)", "\\1", Condition),
        Condition = sub("_(Acet|Succi)_Rep.", "", Condition)) %>%
    dplyr::mutate(., Condition = dplyr::case_when(
        Condition %in% c("Pel_WT", "WT", "BF_WT") ~ names(my_colo)[1],
        Condition %in% c("Pel_dKDAC", "dKDAC (=Ab17Sir2)", "BF_dKDAC") ~ names(my_colo)[2],
        Condition %in% c("Pel_dNpdA", "dNpdA (=CobB)", "BF_dNpdA") ~ names(my_colo)[3],
        Condition %in% c("Pel_dNpdA_dKDAC", "dKDAC_NpdA (=double mutant)", "BF_dNpdA_dKDAC") ~ names(my_colo)[4]
    )) %>%
    unique(.)

my_data$Condition <- factor(
    x = my_data$Condition, levels = names(my_colo), ordered = TRUE)

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "", header = T)

my_data_annot <- my_annot %>%
    dplyr::select(
        ., `Locus Tag`, `KEGG Pathway Name`, `GOBP Term`,
        `GOCC Term`, `GOMF Term`, COG_function,
        `Subcellular Localization [b2g]`) %>%
    dplyr::left_join(
        x = unique(my_data[, c("Accessions ABYAL", "Condition", "PTM")]),
        y = ., by = c("Accessions ABYAL" = "Locus Tag")) %>%
    tidyr::pivot_longer(
        data = ., cols = c(
            `KEGG Pathway Name`, `GOBP Term`,
            `GOCC Term`, `GOMF Term`, COG_function,
            `Subcellular Localization [b2g]`),
        names_to = "DB", values_to = "Category") %>%
    dplyr::mutate(
        ., Category_no_multi = dplyr::case_when(
            is.na(Category) | Category == "" ~ "Not found",
            grepl(";", Category) ~ "Multiple functions",
            TRUE ~ Category),
        Category = dplyr::case_when(
            is.na(Category) | Category == "" ~ "Not found",
            TRUE ~ Category))

my_data_annot %<>%
    tidyr::separate_rows(data = ., Category, sep = ";") %>%
    unique(.)


    #dplyr::group_by(., PTM, DB, Category, Condition) %>%
    #dplyr::summarise(., Count = dplyr::n_distinct(`Accessions ABYAL`)) %>%
    #dplyr::ungroup(.)

for (p in unique(my_data_annot$PTM)) {
    for (d in unique(my_data_annot$DB)) {
        
        my_data_annot_count <- my_data_annot %>%
            dplyr::filter(., PTM == p & DB == d) %>%
            plyr::ddply(
                .data = ., .variables = c("Category", "Condition"),
                .fun = dplyr::summarise,
                Count = dplyr::n_distinct(`Accessions ABYAL`),
                .drop = F)
        
        my_data_annot_count$Category <- factor(
            x = my_data_annot_count$Category,
            levels = rev(sort(unique(my_data_annot_count$Category))),
            ordered = T)
        
        myplots[[paste(p, d, "count", sep = "_")]] <- ggplot(
            my_data_annot_count,
            aes(x = Category, y = Count, group = Condition, fill = Condition)) +
            geom_bar(stat = "identity", position = "dodge", colour = "black") +
            ggpubr::theme_pubr() +
            xlab(d) +
            ylab(paste0("Number of K-", p, "ated proteins")) +
            scale_fill_manual(values = my_colo) +
            coord_flip()
        
    }
}

cairo_pdf(filename = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Functional_annotation/Functional_overview_cairo.pdf", width = 10, height = 20, onefile = T)
for (x in names(myplots)) {
    print(myplots[[x]])
    #grid::grid.newpage()
}
dev.off()

