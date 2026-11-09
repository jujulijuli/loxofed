# Validacion estadistica del meta-analisis: cobertura con k sitios, para distintos tau.
# Uso: Rscript dev/validate_meta.R   (requiere el paquete instalado)
suppressMessages({library(loxofed); library(data.table)})
set.seed(1); R <- 3000; k <- 7; mu <- 0.5
res <- rbindlist(lapply(c(0, 0.1, 0.3, 0.5), function(tau) {
  r <- rbindlist(lapply(seq_len(R), function(i) {
    se <- runif(k, 0.03, 0.25); b <- rnorm(k, mu, sqrt(tau^2 + se^2)); new <- rnorm(1, mu, tau)   # efecto verdadero de un sitio nuevo
    m <- loxofed:::.lox_meta1(b, se, 3L); d <- loxofed:::.lox_dl(b, se)
    tc <- qt(0.975, k - 1)
    data.table(dl_norm = abs(d$beta_re - mu) < 1.96 * d$se_re,
               reml_norm = abs(m$beta_re - mu) < 1.96 * m$se_re,
               reml_hk = abs(m$beta_re - mu) < tc * m$se_hk,
               pi = new > m$beta_re - qt(0.975, k - 2) * sqrt(m$tau2 + m$se_hk^2) & new < m$beta_re + qt(0.975, k - 2) * sqrt(m$tau2 + m$se_hk^2),
               tau2_dl = d$tau2, tau2_reml = m$tau2)
  }))
  data.table(tau = tau, cobertura_DL = mean(r$dl_norm), cobertura_REML_normal = mean(r$reml_norm),
             cobertura_REML_HK = mean(r$reml_hk), cobertura_PI = mean(r$pi),
             tau2_medio_DL = mean(r$tau2_dl), tau2_medio_REML = mean(r$tau2_reml), tau2_verdadero = tau^2)
}))
print(res, digits = 3)
