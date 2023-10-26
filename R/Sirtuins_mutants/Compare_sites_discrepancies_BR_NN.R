


library(magrittr)

my_f <- c(
    #excel = "E:/Brandon_Robin_Sirtuins_mutants/ATCC_Pel_delta_A1S_1281_Intra_Acet_n1_01_ROBIN_20191029_SS.txt",
    excel_R1 = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Formatted/Pel_dKDAC_Acet_Rep1.csv",
    msf_R1 = "E:/Brandon_Robin_Sirtuins_mutants/ATCC_Pel_delta_A1S_1281_Intra_Acet_n1_01_ROBIN_20191029_SS_psms.txt",
    excel_R2 = "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Formatted/Pel_dKDAC_Acet_Rep2.csv",
    msf_R2 = "E:/Brandon_Robin_Sirtuins_mutants/ATCC_DeltaKDAC_intra_Pel_R2_Acet_01_psms.txt")

my_data <- lapply(my_f, function(x) {
    my_sep <- "\t"
    if (sub(".+\\.(.+?$)", "\\1", x) == "csv") { my_sep <- "," }
    my_res <- data.table::fread(
        input = x, sep = my_sep, header = T,
        stringsAsFactors = F, colClasses = "character")
    if ("A2" %in% colnames(my_res)) {
        my_res %<>%
            dplyr::rename(., `Confidence Level` = `A2`)
    }
    if (!"A2" %in% colnames(my_res)) {
        my_res %<>%
            dplyr::group_by(
                ., Sequence, `Protein Group Accessions`, Modifications) %>%
            dplyr::mutate(., `A2` = dplyr::case_when(
                `Confidence Level` == "High" ~ 1,
                `Confidence Level` == "Middle" ~ 2,
                `Confidence Level` == "Low" ~ 3,
            )) %>%
            dplyr::arrange(., `A2`, dplyr::desc(IonScore)) %>%
            dplyr::mutate(., `# PSMs` = dplyr::n()) %>%
            dplyr::slice(., 1) %>%
            dplyr::ungroup(.)
    }
    my_res
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "Origin") %>%
    dplyr::select(
        ., Origin, `Confidence Level`, Sequence, `# PSMs`,
        `Protein Group Accessions`, Modifications,
        `q-Value`, PEP, IonScore, `Exp Value`,
        Charge, `MH+ [Da]`, `m/z [Da]`,
        `RT [min]`, `# Missed Cleavages`)# %>%
    #tidyr::separate(data = ., col = "Origin", into = c("File", "Replicate"), sep = "_")

my_data_filt <- my_data %>%
    dplyr::filter(., grepl("Acetyl", Modifications))

table(my_data_filt$Origin)

my_data_filt %<>%
    dplyr::filter(., `Confidence Level` == "High")

table(my_data_filt$Origin)

my_data_filt %<>%
    dplyr::filter(., IonScore >= 14)

table(my_data_filt$Origin)

my_data_merge <- my_data_filt %>%
    dplyr::select(., Origin, Sequence, `Protein Group Accessions`, Modifications) %>%
    dplyr::group_by(., Sequence, `Protein Group Accessions`, Modifications) %>%
    dplyr::summarise(., Origin = paste0(sort(unique(Origin)), collapse = ";")) %>%
    dplyr::ungroup(.) %>%
    dplyr::left_join(x = ., y = my_data_filt)

data.table::fwrite(
    x = my_data_merge, file = "E:/Brandon_Robin_Sirtuins_mutants/ATCC_Pel_delta_A1S_1281_Intra_Acet_excel_msf_merged.txt",
    append = F, quote = F, sep = "\t", row.names = F, col.names = T)

my_data_psm <- lapply(my_f, function(x) {
    my_res <- data.table::fread(
        input = x, sep = "\t", header = T,
        stringsAsFactors = F, colClasses = "character")
    if ("A2" %in% colnames(my_res)) {
        my_res %<>%
            dplyr::rename(., `Confidence Level` = `A2`)
    }
    my_res
}) %>%
    plyr::ldply(., dplyr::bind_rows, .id = "Origin") %>%
    dplyr::select(
        ., Origin, `Confidence Level`, Sequence, `# PSMs`,
        `Protein Group Accessions`, Modifications,
        `q-Value`, PEP, IonScore, `Exp Value`,
        Charge, `MH+ [Da]`, `m/z [Da]`,
        `RT [min]`, `# Missed Cleavages`)
