


library(magrittr)

my_aln_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf_optimised/Sir2_optimised_Consurf_Outputs_1716881395/msa_fasta.aln"

my_aln <- seqinr::read.alignment(file = my_aln_f, format = "fasta")

my_dist <- seqinr::dist.alignment(x = my_aln, "identity")

my_dist_df <- as.matrix(my_dist) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "ID") %>%
    data.table::data.table(.)

my_dist_target <- my_dist_df %>%
    dplyr::filter(., ID == "Input_seq") %>%
    dplyr::select(., tidyselect::contains("klebsiella_pneumoniae"))


