# import_fs <- function(ns, libs = c(), incl = c(), excl = c()) {
#   tmp <- sapply(libs, library, character.only = TRUE); rm(tmp)
#   if (length(incl) != 0 || length(excl) != 0) {
#     import_list <- getNamespaceImports(ns)
#     if (length(incl) == 0)
#       import_list[names(import_list) %in% c("base", "methods", "stats", "utils", libs, excl)] <- NULL
#     else
#       import_list <- import_list[names(import_list) %in% incl]
#     import_names <- names(import_list)
#
#     for (i in seq_len(length(import_list))) {
#       fun <- import_list[[i]]
#       lib <- import_names[[i]]
#       eval(parse(text = paste0("import::from(",lib,", '",paste0(fun,collapse="', '"),"')")))
#     }
#   }
#   invisible()
# }
#
# if (!"package:intermahpr" %in% search())
#   import_fs("intermahpr", incl = c("magrittr","ggplot2","tidyr","dplyr","broom"))
#
