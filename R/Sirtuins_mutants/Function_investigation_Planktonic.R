


rm(list = ls())

library(magrittr)
library(ggplot2)

myplots <- list()

max_entries <- 100

my_colo <- c("#2B2B2B", "#F58225") %>%
    set_names(c("Planktonic", "Pellicle"))

my_annot_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/Annotation/2023-05-16/Acinetobacter_baumannii_ATCC_17978_full_annotation_2023-05-16.txt"

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
            is.na(Category) | Category == "" ~ "Function unknown",
            grepl(";", Category) ~ "Multiple functions",
            TRUE ~ Category),
        Category = dplyr::case_when(
            is.na(Category) | Category == "" ~ "Function unknown",
            TRUE ~ Category))

my_data_annot %<>%
    tidyr::separate_rows(data = ., Category, sep = ";") %>%
    unique(.)

pdf("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/Functional_overview.pdf", width = 12, height = 12)

for (p in unique(my_data_annot$PTM)) {
    for (d in unique(my_data_annot$DB)) {
        
        my_data_annot_total <- my_data_annot %>%
            dplyr::filter(., PTM == p & DB == d) %>%
            plyr::ddply(
                .data = ., .variables = "Condition",
                .fun = dplyr::summarise,
                Total = dplyr::n_distinct(`Accessions ABYAL`),
                .drop = F)
        
        my_data_annot_count <- my_data_annot %>%
            dplyr::filter(., PTM == p & DB == d) %>%
            plyr::ddply(
                .data = ., .variables = c("Category", "Condition"),
                .fun = dplyr::summarise,
                Count = dplyr::n_distinct(`Accessions ABYAL`),
                .drop = F)
        
        my_data_annot_count %<>%
            dplyr::full_join(x = ., y = my_data_annot_total) %>%
            dplyr::group_by(., Condition) %>%
            dplyr::mutate(., Percentage = Count / Total * 100) %>%
            dplyr::ungroup(.)
        
        facet_min_size <- floor(max_entries/length(unique(my_data_annot_count$Condition)))
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
                aes(x = Category, y = Count, group = Condition, fill = Condition)) +
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
                aes(x = Category, y = Percentage, group = Condition, fill = Condition)) +
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

cairo_pdf(filename = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Comparaison_Pellicle_Planktonic/Functional_overview_cairo.pdf", width = 12, height = 12, onefile = T)
for (x in names(myplots)) {
    print(myplots[[x]])
    #grid::grid.newpage()
}
dev.off()


