


library(magrittr)
library(ggplot2)

my_annot_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/uniref_SIRT1_OR_SIRT2_OR_SIRT3_OR_SI_2024_06_19.tsv"

my_fasta_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/uniref_SIRT1_OR_SIRT2_OR_SIRT3_OR_SI_2024_06_19.fasta"

my_annot <- data.table::fread(input = my_annot_f, sep = "\t", header = T) %>%
    dplyr::mutate(
        ., Name = sub("Cluster: ", "", `Cluster Name`) %>%
            sub("NAD-dependent", "", ., ignore.case = T) %>%
            sub("domain-containing protein", "", ., ignore.case = T) %>%
            sub("family", "", ., ignore.case = T) %>%
            sub("DEACETYLASE", "", ., ignore.case = T) %>%
            sub("DEACYLASE", "", ., ignore.case = T) %>%
            sub("LIPOAMIDASE", "", ., ignore.case = T) %>%
            sub("Histone", "", ., ignore.case = T) %>%
            sub("protein", "", ., ignore.case = T) %>%
            sub("putative", "", ., ignore.case = T) %>%
            sub("SILENT INFORMATION REGULATOR", "", ., ignore.case = T) %>%
            sub("REGULATOR(Y)?", "", ., ignore.case = T) %>%
            sub("TRANSCRIPTIONAL", "", ., ignore.case = T) %>%
            sub("Domain", "", ., ignore.case = T) %>%
            sub("ACTIVE", "", ., ignore.case = T) %>%
            sub("cobalamin(e)?", "", ., ignore.case = T) %>%
            sub("biosynthesis", "", ., ignore.case = T) %>%
            sub("nucleotide", "", ., ignore.case = T) %>%
            sub("binding", "", ., ignore.case = T) %>%
            sub("GLUTAMINE", "", ., ignore.case = T) %>%
            sub("AMIDOTRANSFERASE", "", ., ignore.case = T) %>%
            sub("NAD\\+-DEPENDENT", "", ., ignore.case = T) %>%
            sub("CHROMATIN", "", ., ignore.case = T) %>%
            sub("NADDEPENDENT", "", ., ignore.case = T) %>%
            sub("PROBABLE", "", ., ignore.case = T) %>%
            sub("RELATED TO", "", ., ignore.case = T) %>%
            sub("(-LIKE) ?-LIKE", "\\1", ., ignore.case = T) %>%
            sub(" ([0-9]+)", "-\\1", .) %>%
            sub(" of ", "", ., ignore.case = T) %>%
            sub(",", "", .) %>%
            sub("  +", " ", .) %>%
            sub(" +$", "", .) %>%
            sub("^ +", "", .) %>%
            sub("^-", "", .) %>%
            sub("-$", "", .) %>%
            toupper(.))

my_annot_count <- my_annot %>%
    dplyr::group_by(., `Name`) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(., Name)

data.table::fwrite(
    x = my_annot_count,
    file = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/All_sirtuins_annotation.txt",
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

my_annot_label <- data.table::fread(
    input = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/All_sirtuins_annotation.txt",
    sep = "\t", header = T) %>%
    dplyr::select(., Name, Label) %>%
    unique(.) %>%
    dplyr::left_join(x = my_annot %>% dplyr::mutate(., ID = 1:dplyr::n()), y = ., by = "Name")

ggplot(my_annot_label, aes(x = Label, y = Length)) +
    geom_boxplot() +
    ggpubr::theme_pubr() +
    coord_flip()

my_annot_filt <- my_annot_label_count %>%
    dplyr::filter(., !is.na(`Label`) & `Label` != "") %>%
    dplyr::group_by(., `Label`, `Common taxon ID`, `Common taxon`) %>%
    dplyr::filter(., Length == max(Length, na.rm = T)) %>%
    dplyr::ungroup(.)

my_annot_filt_count <- my_annot_filt %>%
    dplyr::group_by(., `Label`) %>%
    dplyr::summarise(., Count = dplyr::n()) %>%
    dplyr::ungroup(.)

my_annot_filt_split <- split(my_annot_filt, my_annot_filt$Label)

my_annot_filt_select <- lapply(names(my_annot_filt_split), function(x) {
    my_samp <- floor(my_annot_filt_count[my_annot_filt_count$Label == x, ][["Count"]]*0.15) %>%
        ifelse(. < 1, 1, .)
    my_annot_filt_split[[x]][sample(x = nrow(my_annot_filt_split[[x]]), size = my_samp, replace = FALSE), ]
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "Label")

my_lineage <- unique(my_annot_filt_select$`Common taxon ID`) %>%
    split(., ceiling(seq_along(.)/1)) %>%
    lapply(., function(x) {
        tryCatch(
            taxize::classification(sci_id = x, db = "ncbi"),
            error = function(err) NA)
    })

my_lineage_df <- my_lineage %>%
    unlist(., recursive = FALSE) %>%
    plyr::ldply(., data.table::data.table, .id = "Common taxon ID") %>%
    dplyr::mutate(., `Common taxon ID` = sub(".+\\.", "", `Common taxon ID`)) %>%
    dplyr::select(., -V1)

my_lineage_df <- unique(my_lineage_df[is.na(my_lineage_df$name),][["Common taxon ID"]]) %>%
    taxize::classification(sci_id = ., db = "ncbi") %>%
    plyr::ldply(., data.table::data.table, .id = "Common taxon ID") %>%
    dplyr::bind_rows(my_lineage_df[!is.na(my_lineage_df$name),], .) %>%
    unique(.) %>%
    dplyr::mutate(., id = as.integer(id))

my_annot_filt_select_annot <- my_annot_filt_select %>%
    dplyr::left_join(
        x = ., y = my_lineage_df %>% dplyr::select(., rank, id) %>% unique(.),
        by = c("Common taxon ID" = "id"))

my_annot_filt_select_annot_species <- my_annot_filt_select_annot %>%
    dplyr::filter(., rank %in% c("species"))

my_fasta <- seqinr::read.fasta(
    file = my_fasta_f, seqtype = "AA", as.string = T)

my_fasta_filt <- my_fasta[
    which(names(my_fasta) %in% my_annot_filt_select_annot_species$`Cluster ID`)]

seqinr::write.fasta(
    sequences = my_fasta_filt, names = names(my_fasta_filt),
    file.out = sub(".fasta", "_selected.fasta", my_fasta_f),
    open = "w", nbchar = 60)

data.table::fwrite(
    x = my_annot_label, file = sub(".fasta", "_annotated.txt", my_fasta_f),
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)


