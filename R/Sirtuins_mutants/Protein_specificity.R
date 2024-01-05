


library(magrittr)
library(ggplot2)

my_plot <- list()

my_colo <- c("#5757f9ff", "#f94040ff", "#00c000ff", "#fdd61aff") %>%
    set_names(c("WT", "\u0394Sir2-Ab17", "\u0394CobB", "\u0394Sir2-Ab17\u0394CobB"))

my_path <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Condition_explanation"

my_f <- list.files(
    path = my_path, pattern = "_proteins_for_OA.txt", full.names = T) %>%
    set_names(sub("_proteins_for_OA.txt", "", basename(.)))

my_data <- lapply(my_f, function(x) {
    data.table::fread(input = x, sep = "\t", header = T)
}) %>%
    plyr::ldply(., data.table::data.table, .id = "PTM")

my_data_long <- my_data %>%
    tidyr::pivot_longer(data = ., cols = tidyselect::starts_with(c("Comment", "Conditions"))) %>%
    dplyr::filter(., value)

toplot <- my_data_long %>%
    dplyr::filter(., grepl("Comment", name)) %>%
    plyr::ddply(
        .data = ., .variables = c("PTM", "name", "Nmb_modification"),
        .fun = dplyr::summarise,
        Count = dplyr::n_distinct(`Accessions ABYAL`),
        ID = paste0(unique(`Accessions ABYAL`), collapse = ";"), .drop = F)

for (p in unique(my_data_long$PTM)) {
    my_plot[[paste(p, "comment_all", sep = "_")]] <- ggplot(
        toplot %>% dplyr::filter(., PTM == p),
        aes(x = Nmb_modification, y = Count, fill = name)) +
        geom_bar(stat = "identity", position = "dodge") +
        ggpubr::theme_pubr() +
        ggtitle(p)
    my_plot[[paste(p, "comment_filt", sep = "_")]] <- ggplot(
        toplot %>% dplyr::filter(., PTM == p & grepl("Ab17|CobB", name)),
        aes(x = Nmb_modification, y = Count, fill = name)) +
        geom_bar(stat = "identity", position = "dodge") +
        ggpubr::theme_pubr() +
        ggtitle(p)
}

pdf("C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Protein_specificity/Sirtuins_protein_specificity.pdf", 10, 10)
my_plot
dev.off()

data.table::fwrite(
    x = toplot, file = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Protein_specificity/Sirtuins_protein_specificity.txt",
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)


