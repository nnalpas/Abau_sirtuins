


library(magrittr)
library(ggplot2)

my_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/A1S_1281_Abau_isolates/A1S_1281_presence_table.txt"

my_data <- data.table::fread(
    input = my_f, sep = "\t", quote = "", header = T, fill = T)

pdf()

ggplot(my_data, aes(x = `% identity`)) +
    geom_histogram(binwidth = 3, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 3, vjust = -0.5) +
    ylab("Number of isolates")

ggplot(my_data, aes(x = `Alignment length`)) +
    geom_histogram(binwidth = 20, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 20, vjust = -0.5) +
    ylab("Number of isolates")

ggplot(
    my_data %>% dplyr::filter(., `% identity` > 90),
    aes(x = factor(x = `ST (MLST (Pasteur))`, levels = sort(unique(`ST (MLST (Pasteur))`)), ordered = T))) +
    geom_bar(stat = "count", position = "dodge", fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_count(
        aes(y = ..count.., label = ..count..), geom = "text",
        vjust = -0.5) +
    xlab("MLST (Pasteur)") +
    ylab("Number of isolates") +
    ggtitle("Isolates with > 90% identity")

toplot <- my_data %>%
    dplyr::group_by(., `ST (MLST (Pasteur))`) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(`Isolate id`)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., `ST (MLST (Pasteur))` = ifelse(Count > 50, `ST (MLST (Pasteur))`, "Others")) %>%
    dplyr::group_by(., `ST (MLST (Pasteur))`) %>%
    dplyr::summarise(., Count = sum(Count)) %>%
    dplyr::ungroup(.)
toplot$`ST (MLST (Pasteur))` <- factor(
    x = toplot$`ST (MLST (Pasteur))`,
    levels = sort(unique(toplot$`ST (MLST (Pasteur))`)),
    ordered = T)

ggplot(
    toplot,
    aes(x = `ST (MLST (Pasteur))`, y = Count, label = Count)) +
    geom_bar(stat = "identity", position = "dodge", fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    geom_text(stat = "identity", position = position_dodge(width = 0.9), vjust = -0.5) +
    xlab("MLST (Pasteur)") +
    ylab("Number of isolates") +
    ggtitle("All isolates")

dev.off()
