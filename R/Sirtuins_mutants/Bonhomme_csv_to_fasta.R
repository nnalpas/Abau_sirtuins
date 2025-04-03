


rm(list = ls())

library(magrittr)

target_round <- "Round 5"

my_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/SIRim/media-2.csv"

my_data <- data.table::fread(input = my_f, sep = ",", quote = '"', header = T, colClasses = "character")

my_data_filt <- my_data %>%
    dplyr::filter(., Round == target_round) %>%
    dplyr::mutate(
        ., id = paste0("seq", 1:dplyr::n()), Header = paste(id, `Hit id`))

if (length(my_data_filt$`Hit id`) != length(unique(sub(" .*", "", my_data_filt$`Hit id`)))) {
    warning("The identifiers are not unique!")
} else {
    seqinr::write.fasta(
        sequences = as.list(my_data_filt$`Hit sequence`),
        names = my_data_filt$Header,
        file.out = paste0(dirname(my_f), "/Bonhomme_SIR2_subfamily.fasta"),
        as.string = T, open = "w")
}


