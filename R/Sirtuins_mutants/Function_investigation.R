


rm(list = ls())

library(magrittr)
library(ggplot2)

myplots <- list()
myplots2 <- list()

max_entries <- 100

my_colo <- c("#f94040ff", "#00c000ff", "#fdd61aff") %>%
    set_names(c("Only by Ab17Sir2", "Only by CobB", "By Ab17Sir2 & CobB"))

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
        ., Replicates = sub(".+_Rep(.)", "\\1", Condition),
        Condition = sub("_(Acet|Succi)_Rep.", "", Condition)) %>%
    dplyr::mutate(., Condition = dplyr::case_when(
        Condition %in% c("Pel_WT", "WT", "BF_WT") ~ "WT",
        Condition %in% c("Pel_dKDAC", "dKDAC (=Ab17Sir2)", "BF_dKDAC") ~ "ΔSir2-Ab17",
        Condition %in% c("Pel_dNpdA", "dNpdA (=CobB)", "BF_dNpdA") ~ "ΔCobB",
        Condition %in% c("Pel_dNpdA_dKDAC", "dKDAC_NpdA (=double mutant)", "BF_dNpdA_dKDAC") ~ "ΔSir2-Ab17ΔCobB")) %>%
    dplyr::select(., -Modifications, -`Modifications protein`, -Conserved, -ID, -`Accessions A1S`) %>%
    unique(.)

#my_data_hq <- my_data %>%
#    dplyr::group_by(., `Accessions ABYAL`, Position, PTM, Condition, AA) %>%
#    dplyr::summarise(., Rep_count = dplyr::n_distinct(Replicates)) %>%
#    dplyr::ungroup(.) %>%
#    dplyr::group_by(., )

my_data %<>%
    dplyr::group_by(., `Accessions ABYAL`, Position, PTM, AA) %>%
    dplyr::summarise(
        ., Conditions = paste0(sort(unique(Condition)), collapse = ";")) %>%
    #dplyr::mutate(., Comment = dplyr::case_when(
    #    Conditions == "ΔSir2-Ab17;ΔSir2-Ab17ΔCobB" ~ "Only by Ab17Sir2",
    #    Conditions == "ΔCobB;ΔSir2-Ab17;ΔSir2-Ab17ΔCobB" ~ "By Ab17Sir2 & CobB",
    #    Conditions == "ΔCobB;ΔSir2-Ab17ΔCobB" ~ "Only by CobB",
    #    Conditions == "WT;ΔCobB;ΔSir2-Ab17;ΔSir2-Ab17ΔCobB" ~ "Common",
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

my_data$Comment <- factor(
    x = my_data$Comment, levels = names(my_colo), ordered = TRUE)

my_annot <- data.table::fread(
    input = my_annot_f, sep = "\t", quote = "", header = T)

my_data_annot <- my_annot %>%
    dplyr::select(
        ., `Locus Tag`, `KEGG Pathway Name`, `GOBP Term`,
        `GOCC Term`, `GOMF Term`, COG_function,
        `Subcellular Localization [b2g]`) %>%
    dplyr::left_join(
        x = unique(my_data[, c("Accessions ABYAL", "Comment", "PTM")]),
        y = ., by = c("Accessions ABYAL" = "Locus Tag")) %>%
    tidyr::pivot_longer(
        data = ., cols = c(
            `KEGG Pathway Name`, `GOBP Term`,
            `GOCC Term`, `GOMF Term`, COG_function,
            `Subcellular Localization [b2g]`),
        names_to = "DB", values_to = "Category") %>%
    dplyr::mutate(
        ., Category_no_multi = dplyr::case_when(
            is.na(Category) | Category == "" ~ "Function unknown",
            grepl(";", Category) ~ "Multiple functions",
            TRUE ~ Category),
        Category = dplyr::case_when(
            is.na(Category) | Category == "" ~ "Function unknown",
            TRUE ~ Category))

my_data_annot %<>%
    tidyr::separate_rows(data = ., Category, sep = ";") %>%
    unique(.)

pdf("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Functional_annotation/Functional_overview.pdf", width = 12, height = 12)

for (p in unique(my_data_annot$PTM)) {
    for (d in unique(my_data_annot$DB)) {
        
        my_data_annot_total <- my_data_annot %>%
            dplyr::filter(., PTM == p & DB == d & Comment %in% names(my_colo)) %>%
            plyr::ddply(
                .data = ., .variables = "Comment",
                .fun = dplyr::summarise,
                Total = dplyr::n_distinct(`Accessions ABYAL`),
                .drop = F)
        
        my_data_annot_count <- my_data_annot %>%
            dplyr::filter(., PTM == p & DB == d & Comment %in% names(my_colo)) %>%
            plyr::ddply(
                .data = ., .variables = c("Category", "Comment"),
                .fun = dplyr::summarise,
                Count = dplyr::n_distinct(`Accessions ABYAL`),
                .drop = F)
        
        my_data_annot_count %<>%
            dplyr::full_join(x = ., y = my_data_annot_total) %>%
            dplyr::group_by(., Comment) %>%
            dplyr::mutate(., Percentage = Count / Total * 100) %>%
            dplyr::ungroup(.)
        
        facet_min_size <- floor(max_entries/length(unique(my_data_annot_count$Comment)))
        panels_df <- split(
            x = sort(unique(my_data_annot_count$Category)),
            f = ceiling(seq_along(
                sort(unique(my_data_annot_count$Category))) / facet_min_size)) %>%
            plyr::ldply(., data.table::data.table) %>%
            set_colnames(c("Panel", "Category"))
        
        my_data_annot_count %<>%
            dplyr::left_join(x = ., y = panels_df)
        
        my_data_annot_count$Category <- factor(
            x = my_data_annot_count$Category,
            levels = rev(sort(unique(my_data_annot_count$Category))),
            ordered = T)
        
        for (i in unique(my_data_annot_count$Panel)) {
            
            myplots[[paste(p, d, "count", i, sep = "_")]] <- ggplot(
                my_data_annot_count %>% dplyr::filter(., Panel == i),
                aes(x = Category, y = Count, group = Comment, fill = Comment)) +
                geom_bar(stat = "identity", position = "dodge", colour = "black") +
                ggpubr::theme_pubr() +
                xlab(d) +
                ylab(paste0("Number of K-", p, "ated proteins")) +
                scale_fill_manual(values = my_colo) +
                coord_flip() +
                ggtitle(paste("Part", i))
            print(myplots[[paste(p, d, "count", i, sep = "_")]])
            
            myplots[[paste(p, d, "percent", i, sep = "_")]] <- ggplot(
                my_data_annot_count %>% dplyr::filter(., Panel == i),
                aes(x = Category, y = Percentage, group = Comment, fill = Comment)) +
                geom_bar(stat = "identity", position = "dodge", colour = "black") +
                ggpubr::theme_pubr() +
                xlab(d) +
                ylab(paste0("Percentage of K-", p, "ated proteins")) +
                scale_fill_manual(values = my_colo) +
                coord_flip() +
                ggtitle(paste("Part", i))
            print(myplots[[paste(p, d, "percent", i, sep = "_")]])
            
        }
        
    }
}

dev.off()

cairo_pdf(filename = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Functional_annotation/Functional_overview_COG.pdf", width = 12, height = 12, onefile = T)
for (x in grep('COG', names(myplots), value = T)) {
    print(myplots[[x]])
    #grid::grid.newpage()
}
dev.off()


