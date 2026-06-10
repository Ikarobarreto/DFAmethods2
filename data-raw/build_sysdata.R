# data-raw/build_sysdata.R
#
# Builds R/sysdata.rda from the theoretical kappa table shipped with the paper
# "Auditing the assumptions of DFA-based multiple regression" (Barreto et al.
# 2026, Table A.1). Run interactively when the table is updated; the resulting
# R/sysdata.rda is what ships in the package.

raw <- read.csv("data-raw/kappa_theoretical.csv")
stopifnot(all(c("s", "H", "kappa_th") %in% names(raw)))

kappa_th_table <- raw[, c("s", "H", "kappa_th")]
kappa_th_table <- kappa_th_table[order(kappa_th_table$H, kappa_th_table$s), ]
rownames(kappa_th_table) <- NULL

save(kappa_th_table, file = "R/sysdata.rda", compress = "xz",
     compression_level = 9)
message("Wrote R/sysdata.rda with ", nrow(kappa_th_table),
        " rows (s in [", min(kappa_th_table$s), ", ", max(kappa_th_table$s),
        "], H in [", min(kappa_th_table$H), ", ", max(kappa_th_table$H), "]).")
