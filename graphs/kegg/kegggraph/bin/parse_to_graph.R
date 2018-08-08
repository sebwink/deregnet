library(KEGGgraph)

parse_to_dataframe <- function(kgml_dir) {
    kgmls <- list.files(path=kgml_dir, pattern="*.xml", full.names=TRUE, recursive=FALSE)
    dfs <- lapply(kgmls, parseKGML2DataFrame, expandGenes=TRUE)
    return(do.call(rbind, dfs))
}

df_to_sif <- function(df, path) {
    df <- df[,c(1,3,2)]
    write.table(df, path, row.names=FALSE, col.names=FALSE, fileEncoding="utf8", sep="\t", quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
kgml_dir = args[1]
graph_path = args[2]

df <- parse_to_dataframe(kgml_dir)
df_to_sif(df, graph_path)
