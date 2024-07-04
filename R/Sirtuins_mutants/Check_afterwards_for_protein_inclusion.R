
my_f <- "C:/Users/nalpanic/SynologyDrive/Work/Colleagues shared work/Brandon_Robin/Abaumannii_mutants/Analysis/Sirtuin_conservation/UniProt/Tree_2024-06-27/Class_U_ids.txt"

my_classu <- data.table::fread(input = my_f, sep = "\t", header = T)

my_classu[my_classu$Accession %in% my_gene$UniProtID, ]


my_gene %<>%
    dplyr::mutate(., Target_type = dplyr::case_when(
        `Cluster Name` %in% unique(my_gene$`Cluster Name`)[1:13] ~ `Cluster Name`,
        UniProtID == my_classu$Accession[my_classu$Accession %in% my_gene$UniProtID] ~ "Sir2",
        TRUE ~ NA_character_
    ))

label_col <- "Target_type"

table(my_gene[[label_col]])

tocolour <- my_gene %>%
    dplyr::filter(., !is.na(!!as.name(label_col))) %>%
    dplyr::select(., !!as.name(label_col)) %>%
    unique(.) %>%
    dplyr::bind_cols(., Colour = my_big_cols[1:nrow(.)])

colors <- my_gene %>%
    dplyr::filter(., UniProtID %in% my_ids$UniProtID) %>%
    dplyr::left_join(
        x = ., y = tocolour,
        by = label_col) %>%
    #dplyr::mutate(
    #    ., Colour = tidyr::replace_na(
    #        data = Colour, replace = "#404040")) %>%
    dplyr::arrange(., match(UniProtID, my_ids$UniProtID))

tipcol <- colors[["Colour"]]

legend <- colors %>%
    dplyr::select(., !!as.name(label_col), Colour)
legend %<>%
    unique(.)

plot(x = my_tree, "unrooted", main = "NJ", tip.color = tipcol, show.tip.label = F)
tiplabels(pch = 16, col = tipcol, cex = 0.7)
i <- which(my_tree$tip.label %in% my_gene$UniProtID[!is.na(my_gene$Target_type)])
tiplabels(my_tree$tip.label[i], i, col = tipcol, adj = 0.5, frame = "none")
legend("bottomright", legend = legend[[label_col]],
       title = "Legend",
       col = legend[["Colour"]], pch = 16, cex = 0.8)
add.scale.bar(cex = 1, font = 2, col = "black")


