
library(magrittr)

my_data <- data.table::fread(
    input = "E:/Annotation_Brassica.tab", sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, integer64 = "double", data.table = FALSE)

my_data_format <- my_data %>%
    dplyr::select(., -`Gene ontology IDs`, -`Gene ontology (GO)`) %>%
    dplyr::mutate_all(.funs = gsub, pattern = " +$", replacement = "") %>%
    dplyr::mutate_all(.funs = gsub, pattern = ";+$", replacement = "") %>%
    dplyr::mutate_all(.funs = gsub, pattern = "; +", replacement = ";") %>%
    unique(.)

my_data_go <- my_data_format %>%
    dplyr::select(., Entry, dplyr::matches("Gene ontology \\("))

my_data_nogo <- my_data_format %>%
    dplyr::select(., -dplyr::matches("Gene ontology \\("))


my_data_final <- my_data_go %>%
    dplyr::select(., Entry, `Gene ontology (biological process)`) %>%
    splitstackshape::cSplit(
        indt = ., splitCols = "Gene ontology (biological process)",
        sep = ";", direction = "long", fixed = TRUE, drop = FALSE,
        type.convert = FALSE) %>%
    tidyr::separate(
        data = ., col = "Gene ontology (biological process)",
        into = c("GOBP Name", "GOBP ID"), sep = " \\[GO:",
        remove = TRUE, convert = FALSE) %>%
    dplyr::mutate(., `GOBP ID` = sub("^", "GO:", `GOBP ID`)) %>%
    dplyr::mutate(., `GOBP ID` = sub("\\]", "", `GOBP ID`)) %>%
    dplyr::group_by(., Entry) %>%
    dplyr::summarise(
        ., `GOBP Name` = paste(unique(`GOBP Name`), collapse = ";"),
        `GOBP ID` = paste(unique(`GOBP ID`), collapse = ";")) %>%
    dplyr::left_join(x = my_data_nogo, y = ., by = "Entry")


my_data_final <- my_data_go %>%
    dplyr::select(., Entry, `Gene ontology (cellular component)`) %>%
    splitstackshape::cSplit(
        indt = ., splitCols = "Gene ontology (cellular component)",
        sep = ";", direction = "long", fixed = TRUE, drop = FALSE,
        type.convert = FALSE) %>%
    tidyr::separate(
        data = ., col = "Gene ontology (cellular component)",
        into = c("GOCC Name", "GOCC ID"), sep = " \\[GO:",
        remove = TRUE, convert = FALSE) %>%
    dplyr::mutate(., `GOCC ID` = sub("^", "GO:", `GOCC ID`)) %>%
    dplyr::mutate(., `GOCC ID` = sub("\\]", "", `GOCC ID`)) %>%
    dplyr::group_by(., Entry) %>%
    dplyr::summarise(
        ., `GOCC Name` = paste(unique(`GOCC Name`), collapse = ";"),
        `GOCC ID` = paste(unique(`GOCC ID`), collapse = ";")) %>%
    dplyr::left_join(x = my_data_final, y = ., by = "Entry")



my_data_final <- my_data_go %>%
    dplyr::select(., Entry, `Gene ontology (molecular function)`) %>%
    splitstackshape::cSplit(
        indt = ., splitCols = "Gene ontology (molecular function)",
        sep = ";", direction = "long", fixed = TRUE, drop = FALSE,
        type.convert = FALSE) %>%
    tidyr::separate(
        data = ., col = "Gene ontology (molecular function)",
        into = c("GOMF Name", "GOMF ID"), sep = " \\[GO:",
        remove = TRUE, convert = FALSE) %>%
    dplyr::mutate(., `GOMF ID` = sub("^", "GO:", `GOMF ID`)) %>%
    dplyr::mutate(., `GOMF ID` = sub("\\]", "", `GOMF ID`)) %>%
    dplyr::group_by(., Entry) %>%
    dplyr::summarise(
        ., `GOMF Name` = paste(unique(`GOMF Name`), collapse = ";"),
        `GOMF ID` = paste(unique(`GOMF ID`), collapse = ";")) %>%
    dplyr::left_join(x = my_data_final, y = ., by = "Entry")


write.table(
    x = my_data_final, file = "E:/Annotation_Brassica.txt",
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


