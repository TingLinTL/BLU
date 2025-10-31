R <- 100                       # number of replicates
source("main_sim_one_rep.R")   # loads run_one_rep()

for (r in 1:R) {
  cat(">>> Replicate", r, "starting...\n")
  res_r <- run_one_rep(r)
  saveRDS(res_r, sprintf("Results_rep%03d.rds", r))
  cat(">>> Replicate", r, "finished.\n")
}

# (optional) combine at the end:
all_results <- do.call(cbind, lapply(1:R, function(r) readRDS(sprintf("Results_rep%03d.rds", r))))
saveRDS(all_results, "Results_ALL.rds")

true_spce_aft <- -0.1384194

summary_results <- apply(all_results, 2, function(col) {
  mean_val <- mean(col)
  sd_val <- sd(col)
  CI <- as.numeric(quantile(col, c(0.025, 0.975)))  # remove names
  lo <- CI[1]
  hi <- CI[2]
  cover <- ifelse(true_spce_aft >= lo & true_spce_aft <= hi, "Yes", "No")
  c(mean = mean_val, sd = sd_val, lower = lo, upper = hi, cover = cover)
})

summary_results <- t(summary_results)
summary_results
write.csv(summary_results, file = "summary_results.csv", row.names = FALSE)
mean(as.numeric(summary_results[,1]))
