


library(magrittr)

my_path <- 'C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf'

my_f <- list.files(path = my_path, pattern = "msa_fasta.aln", full.names = T, recursive = T) %>%
    set_names(sub("_Consurf_Outputs_.+", "", basename(dirname(.))))

my_fasta <- lapply(my_f, function(x) {
    seqinr::read.fasta(file = x, seqtype = "AA", as.string = T, whole.header = T)
}) %>%
    unlist(.) %>%
    gsub("-", "", .)

my_fasta_format <- my_fasta %>%
    set_names(sub(".+?\\|", "", names(.))) %>%
    set_names(sub("\\|.+", "", names(.)))

my_fasta_format <- my_fasta_format[!duplicated(names(my_fasta_format))]

seqinr::write.fasta(sequences = as.list(my_fasta_format), names = names(my_fasta_format), file.out = paste(my_path, "All_consurf_sequences.fasta", sep = "/"), open = "w")

writeLines(text = names(my_fasta), con = paste(my_path, "All_consurf_sequences_header.txt", sep = "/"))


