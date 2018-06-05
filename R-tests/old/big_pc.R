library("intermahpr")

PC <- intermahpr::pc_default
  # data.table::fread("data-raw/impc.csv")
RR <- intermahpr::rr_default

betas <- data.matrix(rep(0,16))
betas

curves <- intermahpr::deduce_relative_risk_curves_from_rr()

plot_rr <- function(rr_item) {
  fn_rr <- rr_item[[3]]
  pts <- seq(0.03,250,0.1)
  data <- data.frame(x = pts, y = fn_rr(pts))
  print(ggplot(data = data, mapping = aes(x = x, y = y)) +
    geom_line()
  )
}

curves[1:10]

L <- lapply(curves[1:10], plot_rr)


ggplot(data = data.frame(x=0), mapping = aes(x=x)) +
  stat_function(fun = curves[[1]][[3]]) +
  xlim(-250,250)

curves[[27]][[3]](13)

aaf <- compute_all_aafs(PC, RR)



microbenchmark::microbenchmark(compute_all_aafs(PC, RR))

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

