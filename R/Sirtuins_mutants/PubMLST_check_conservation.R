


rm(list = ls())

library(magrittr)
library(ggplot2)

length_prot <- 473
length_prot <- 232

my_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/pubMLST_Blast/cobB_presence_table.txt"

my_data <- data.table::fread(
    input = my_f, sep = "\t", quote = "", header = T, fill = T) %>%
    dplyr::filter(., species == "Acinetobacter baumannii")

pdf(sub(".txt", ".pdf", my_f), 10, 10)

toplot <- my_data %>%
    dplyr::group_by(., `ST (MLST (Pasteur))`) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(`Isolate id`)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., `ST (MLST (Pasteur))` = ifelse(Count > 10, `ST (MLST (Pasteur))`, "Others")) %>%
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
    geom_text(
        stat = "identity", position = position_dodge(width = 0.9), hjust = -0.1, angle = 90) +
    xlab("MLST (Pasteur)") +
    ylab("Number of isolates") +
    ggtitle("All isolates") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

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

ggplot(my_data, aes(x = `E-value`)) +
    geom_histogram(binwidth = 0.1, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 0.1, vjust = -0.5) +
    ylab("Number of isolates")

ggplot(my_data, aes(x = `Bit score`)) +
    geom_histogram(binwidth = 20, fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_bin(
        aes(y=..count.., label=..count..), geom = "text",
        binwidth = 20, vjust = -0.5) +
    ylab("Number of isolates")

ggplot(
    my_data %>% dplyr::filter(., `% identity` > 90 & `Alignment length`/length_prot > 0.90 & `E-value` < 0.1),
    aes(x = factor(x = `ST (MLST (Pasteur))`, levels = sort(unique(`ST (MLST (Pasteur))`)), ordered = T))) +
    geom_bar(stat = "count", position = "dodge", fill = "darkgrey", colour = "black") +
    ggpubr::theme_pubr() +
    stat_count(
        aes(y = ..count.., label = ..count..), geom = "text",
        vjust = -0.5) +
    xlab("MLST (Pasteur)") +
    ylab("Number of isolates") +
    ggtitle("Isolates with > 90% identity over 90% sequence and E-value < 0.1")

toplot <- my_data %>%
    dplyr::group_by(., `ST (MLST (Pasteur))`) %>%
    dplyr::summarise(
        ., Total = dplyr::n_distinct(`Isolate id`),
        Conserved = sum(`% identity` > 90 & `Alignment length`/length_prot > 0.90 & `E-value` < 0.1),
        Percentage = Conserved / Total * 100) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., dplyr::desc(Conserved))

ggplot(
    toplot %>% dplyr::filter(., `Percentage` > 0),
    aes(x = factor(x = `ST (MLST (Pasteur))`, levels = sort(unique(`ST (MLST (Pasteur))`)), ordered = T), y = Conserved, label = paste(round(Percentage, digit = 2), "%"))) +
    geom_bar(stat = "identity", position = "dodge", fill = "darkgrey", colour = "black") +
    geom_text(stat = "identity", position = position_dodge(width = 0.9), vjust = -0.3) +
    ggpubr::theme_pubr() +
    xlab("MLST (Pasteur)") +
    ylab("Number of isolates") +
    ggtitle("Isolates with > 90% identity over 90% sequence and E-value < 0.1")

dev.off()


