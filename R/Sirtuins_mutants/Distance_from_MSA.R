


rm(list = ls())

library(magrittr)
library(foreach)
library(doParallel)
library(tcltk)

cl <- makeCluster(12)
registerDoParallel(cl)

#my_aln_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/ConSurf_optimised/Sir2_optimised_Consurf_Outputs_1716881395/msa_fasta.aln"
my_aln_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/Tree_2024-06-27/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25_selected2_msa.faa"

my_aln <- seqinr::read.alignment(file = my_aln_f, format = "fasta")

my_dist <- seqinr::dist.alignment(x = my_aln, "identity", gap = FALSE)

my_dist_df <- as.matrix(my_dist) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(.data = ., var = "ID") %>%
    data.table::data.table(.)

#my_dist_target <- my_dist_df %>%
#    dplyr::filter(., ID == "Input_seq") %>%
#    dplyr::select(., tidyselect::contains("klebsiella_pneumoniae"))


identity <- function(a, b, gap = c("-", ".")) {
    a <- unlist(strsplit(x = a, split = ""))
    b <- unlist(strsplit(x = b, split = ""))
    if (length(a) != length(b)) {
        stop(paste0(
            "Lengths of the two sequences (",
            length(a), " vs. ", length(b),
            ") are not equal!"))
    }
    identity <- 0
    ignore <- 0
    for (i in 1:length(a)) {
        identity <- identity + (a[i] == b[i])
        ignore <- ignore + all(c(a[i], b[i]) %in% gap)
    }
    list(
        identity = (identity - ignore),
        identity_percentage = (identity - ignore) / (length(a) - ignore))
}

identity_par <- function(x) {
    identity(a = x[["a"]][[1]], b = x[["b"]][[1]], gap = "-")
}



a <- my_aln$seq %>%
    set_names(my_aln$nam)

dt <- tidyr::crossing(a = a, b = a)

identity_df <- foreach(
    x = iterators::iter(dt, by = 'row'),
    .inorder = T, .packages = "tcltk", .combine = "rbind") %dopar%
    {
        #if(!exists("pb")) {
        #    pb <- tkProgressBar("Parallel task", min = 1, max = 10000)
        #}
        #setTkProgressBar(pb, value = "v")
        identity_par(x)
    }

stopCluster(cl)

identity_df %<>%
    cbind(tidyr::crossing(a_name = names(a), b_name = names(a)), .)

data.table::fwrite(
    x = identity_df, file = sub("\\.faa", "_identity.txt", my_aln_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)


