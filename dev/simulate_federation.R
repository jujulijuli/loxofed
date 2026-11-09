# Federacion simulada de 7 sitios con heterogeneidad conocida: prueba de punta a punta y de recuperacion de la verdad.
# Uso: Rscript dev/simulate_federation.R [carpeta_salida] [tau]
suppressMessages({library(loxofed); library(data.table)})
a <- commandArgs(TRUE); root <- if (length(a) >= 1) a[1] else file.path(tempdir(), "fed7"); tau <- if (length(a) >= 2) as.numeric(a[2]) else 0.25
dir.create(root, recursive = TRUE, showWarnings = FALSE)
sites <- data.table(id = c("NAV", "ARA", "AND", "CAT", "GAL", "VAL", "MAD"), n = c(7000, 5000, 9000, 6000, 4500, 8000, 5500),
                    res = c("day", "day", "month", "day", "month", "day", "day"), cap = c(NA, NA, 4000, NA, NA, NA, NA))
truth <- list(); cfgs <- list()
for (i in seq_len(nrow(sites))) {
  s <- sites[i]; d <- file.path(root, s$id, "data")
  tr <- attr(lox_simulate_cdm(s$n, d, seed = 100 + i, date_res = s$res, tau = tau, het_seed = 900 + i), "truth")
  tr[, site := s$id]; truth[[i]] <- tr
  cfgs[[i]] <- lox_config(site_id = s$id, paths = list(pacientes = file.path(d, "pacientes.csv"), enfermedades = file.path(d, "enfermedades.csv"),
      output_dir = file.path(root, s$id, "out")), cdm = list(sep = "|"), cohort = list(sample_n = if (is.na(s$cap)) NA else s$cap),
      structure = list(boot_R = 20, seed = i), model = list(min_events_node = 40, min_events_hazard = 10))
}
truth <- rbindlist(truth)
for (cf in cfgs) suppressMessages(lox_run_round1(cf))
cen <- file.path(root, "central")
suppressMessages(lox_run_central_consensus(file.path(root, sites$id, "out", "round1"), cen, cfgs[[1]]))
for (cf in cfgs) suppressMessages(lox_run_round2(cf, file.path(cen, "consensus")))
mt <- suppressMessages(lox_run_central_meta(file.path(root, sites$id, "out", "round2"), cen, cfgs[[1]]))
# --- recuperacion: efecto combinado vs media de las verdades de los sitios; tau verdadera vs tau estimada
tm <- truth[, .(beta_true = mean(beta), tau_true = sd(beta), beta_base = beta_base[1]), by = .(from, to)]
x <- merge(mt[capa == "incidencia" & status == "ok", .(to = nodo, from = sub("_t0$", "", termino), beta_re, se_hk, tau2, I2, k)], tm, by = c("from", "to"))
x[, `:=`(err = beta_re - beta_true, tau_est = sqrt(tau2))]
cat("\nArcos verdaderos recuperados:", nrow(x), "de", nrow(tm), "\n"); print(x[order(-beta_true), .(from, to, k, beta_true = round(beta_true, 2), beta_re = round(beta_re, 2), tau_true = round(tau_true, 2), tau_est = round(tau_est, 2), I2 = round(I2, 2))])
cat(sprintf("\nError medio %.3f; |error| medio %.3f; correlacion verdad-estimacion %.3f\n", mean(x$err), mean(abs(x$err)), cor(x$beta_true, x$beta_re)))
cat(sprintf("Cobertura del IC95%% HK sobre la media verdadera: %d de %d\n", sum(abs(x$err) < qt(0.975, x$k - 1) * x$se_hk), nrow(x)))
invisible(list(root = root, mt = mt, truth = truth))
