library("intermahpr")

PC <- data.table::fread("data-raw/impc.csv")
RR <- intermahpr::rr_default

microbenchmark::microbenchmark(compute_total_aafs(PC = PC, RR = RR))

params <- derive_params_from_pc(PC)
curves <- deduce_relative_risk_curves_from_rr(RR)


aaf_arguments <- lapply(curves,
                        function(C) apply(params[Gender == C[[1]][[3]], ],
                                          1,
                                          function(P) combine(C,P)))

aaf_list <- lapply(aaf_arguments, function(x) lapply(x, compute_aaf))

lapply(aaf_list, function(x) print(length(x)))

aaf_frames <- lapply(aaf_list,
                     function(x) data.frame(matrix(unlist(x),
                                                   nrow = length(x),
                                                   byrow = T)))

aaf_frame <- do.call("rbind", aaf_frames)

aaf_frame_rename <- plyr::rename(aaf_frame, c("X1" = "Region",
                          "X2" = "Year",
                          "X3" = "Gender",
                          "X4" = "Age_Group",
                          "X5" = "IM",
                          "X6" = "Condition",
                          "X7" = "Outcome",
                          "X8" = "AAF_Total"))

View(aaf_frame_rename)

