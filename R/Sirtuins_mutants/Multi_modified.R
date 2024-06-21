


rm(list = ls())

library(magrittr)
library(ggplot2)

my_other_ptm_f <- "C:/Users/nalpanic/SynologyDrive/Work/Abaumannii_trimeth/dbPTM/"

my_other_ptm_files <- list.files(path = my_other_ptm_f, pattern = "formatted_per_sites.txt", recursive = F, full.names = T) %>%
    set_names(sub("_formatted_per_sites.txt", "", basename(.))) %>%
    grep("Robin", ., value = T)

my_other_ptm <- lapply(my_other_ptm_files, function(x) {
    data.table::fread(input = x, sep = "\t", header = T)
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "PTMDB") %>%
    dplyr::mutate(., AA = sub("(.).+", "\\1", `Modifications protein`)) %>%
    dplyr::filter(., !is.na(Position)) %>%
    dplyr::filter(., grepl("Pel", Condition))

my_multi <- my_other_ptm %>%
    dplyr::group_by(., `Accessions ABYAL`, Position) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(PTM)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(., Condition = "Total")

my_multi <- my_other_ptm %>%
    tidyr::separate_rows(data = ., Condition, sep = ";") %>%
    dplyr::filter(., grepl("^Pel", Condition)) %>%
    dplyr::mutate(., Condition = sub("(_Succi_|_Acet_).+", "", Condition)) %>%
    dplyr::group_by(., `Accessions ABYAL`, Condition, Position) %>%
    dplyr::summarise(., Count = dplyr::n_distinct(PTM)) %>%
    dplyr::ungroup(.) %>%
    dplyr::bind_rows(my_multi, .)

ggplot(my_multi %>% dplyr::filter(., Count == 2), aes(x = Condition, fill = Condition, colour = Condition)) +
    geom_bar(stat = "count", position = "dodge", alpha = 0.7) +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.3) +
    ggpubr::theme_pubr()

