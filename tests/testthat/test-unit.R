test_that("el parseo de fechas detecta la resolucion e imputa", {
  r <- loxofed:::.lox_parse_date(c("2012", "2012-03", "2012-03-05", "", "basura"))
  expect_equal(r$res, c("year", "month", "day", NA, NA))
  expect_equal(format(r$date[1:3]), c("2012-07-01", "2012-03-15", "2012-03-05"))
  expect_true(all(is.na(r$date[4:5])))
})

test_that("la huella cambia con lo que debe coincidir y no con lo local", {
  a <- lox_config(); fp <- lox_config_fingerprint(a)
  expect_equal(fp, lox_config_fingerprint(lox_config(site_id = "OTRO", structure = list(seed = 99),
                                                      paths = list(output_dir = "/x"))))
  expect_false(fp == lox_config_fingerprint(lox_config(cohort = list(age_min = 65))))
  expect_false(fp == lox_config_fingerprint(lox_config(model = list(age_knots = c(72, 82, 92)))))
})

test_that("la configuracion se valida", {
  expect_error(lox_config(cohort = list(year_first = 2020, year_last = 2012)), "year_last")
  expect_error(lox_config(model = list(age_knots = c(50, 80, 90))), "age_knots")
  expect_warning(lox_config(disclosure = list(min_cell = 3)), "min_cell")
})

test_that("la base de spline depende solo de los nudos fijos, no de los datos", {
  cfg <- lox_config()
  b1 <- lox_age_basis(c(60, 70, 85, 95), cfg); b2 <- lox_age_basis(c(70, 100, 62, 90, 75), cfg)
  expect_equal(names(b1), paste0("edad_ns", 1:4))
  expect_equal(unlist(b1[2]), unlist(b2[1]))        # la edad 70 se codifica igual en cualquier sitio
})

test_that("la supresion anula celdas pequenas y el total complementario", {
  x <- data.table(a = c(5L, 50L, 0L), b = c(30L, 40L, 20L), tot = c(35L, 90L, 20L))
  s <- lox_suppress(x, c("a", "b"), 10, total = "tot")
  expect_equal(s$suppressed, c(TRUE, FALSE, FALSE))
  expect_true(is.na(s$a[1]) && is.na(s$tot[1]))
  expect_equal(s$tot[2:3], c(90L, 20L))             # 0 no es celda pequena
})

test_that("meta-analisis REML + Hartung-Knapp: casos con solucion cerrada", {
  mk <- function(b, se) data.table(site = LETTERS[seq_along(b)], capa = "x", nodo = "n", termino = "t", beta = b, se = se)
  # homogeneo, 3 sitios: tau2 = 0 y efecto = media
  m <- lox_meta(mk(c(0.5, 0.5, 0.5), 0.1)); expect_equal(m$beta_re, 0.5); expect_equal(m$tau2, 0); expect_equal(m$status, "ok")
  # varianzas iguales: tau2_REML = var(b) - se^2 ; EE_HK = sqrt(var(b)/k) ; PI con t(k-2)
  b <- c(0.1, 0.5, 0.9); se <- 0.05; m <- lox_meta(mk(b, se))
  expect_equal(m$tau2, var(b) - se^2, tolerance = 1e-6); expect_equal(m$beta_re, 0.5, tolerance = 1e-8)
  expect_gt(m$I2, 0.9); expect_equal(m$lo, exp(0.5 - qt(0.975, 2) * sqrt(var(b) / 3)), tolerance = 1e-6)
  expect_equal(m$pi_hi, exp(0.5 + qt(0.975, 1) * sqrt(var(b) - se^2 + var(b) / 3)), tolerance = 1e-6)
  expect_true(m$lo < m$exp_re && m$exp_re < m$hi && m$pi_lo < m$lo && m$pi_hi > m$hi)
  expect_equal(m$beta_dl, 0.5, tolerance = 1e-8); expect_equal(m$tau2_dl, var(b) - se^2, tolerance = 1e-6)  # DL = REML con varianzas iguales
  # el IC de Hartung-Knapp nunca es mas estrecho que el normal de REML
  mh <- lox_meta(mk(c(0.50, 0.51, 0.49, 0.50), c(0.05, 0.06, 0.05, 0.04)))
  expect_gte(log(mh$hi) - log(mh$lo), 2 * 1.96 * mh$se_re - 1e-9)
  # k insuficiente: no hay resultado combinado, pero se informa k y el efecto de efectos fijos
  m2 <- lox_meta(mk(c(0.4, 0.6), c(0.1, 0.1))); expect_match(m2$status, "insuficiente"); expect_true(is.na(m2$beta_re)); expect_equal(m2$k, 2L)
  expect_equal(lox_meta(mk(c(0.4, 0.6), c(0.1, 0.1)), min_k = 2L)$status, "ok")
  expect_equal(nrow(lox_meta(mk(0.5, 0.1))), 1)                # un solo sitio no rompe
  # dejar uno fuera: el sitio atipico mueve el rango
  mo <- lox_meta(mk(c(0.5, 0.5, 0.5, 2), 0.1)); expect_true(mo$loo_min < mo$exp_re && mo$loo_max > mo$loo_min)
  expect_equal(mo$loo_min, exp(0.5), tolerance = 1e-6)
  # sitios con se invalido se descartan
  expect_equal(lox_meta(rbind(mk(c(0.5, 0.5, 0.5), 0.1), data.table(site = "Z", capa = "x", nodo = "n", termino = "t", beta = 9, se = NA)))$k, 3L)
})

test_that("la config de ejemplo carga, deja age_max sin limite y coincide con los valores comunes", {
  ex <- system.file("config", "config.example.yml", package = "loxofed")
  if (!nzchar(ex)) ex <- file.path("..", "..", "inst", "config", "config.example.yml")
  cfg <- lox_config(ex)
  expect_true(is.na(cfg$cohort$age_max)); expect_equal(cfg$cohort$age_min, 60)
  expect_equal(cfg$site_id, "NAV"); expect_equal(cfg$cdm$cols$baja, "F_Baja")
  expect_equal(lox_config_fingerprint(cfg), lox_config_fingerprint(lox_config()))
})

test_that("sample_n limita personas, conserva sus enfermedades y no cambia la huella", {
  d <- file.path(tempdir(), "samp"); lox_simulate_cdm(n = 6000, out_dir = d, seed = 1)
  cfg <- lox_default_config(); cfg$cdm$sep <- "|"
  cfg$paths$pacientes <- file.path(d, "pacientes.csv"); cfg$paths$enfermedades <- file.path(d, "enfermedades.csv")
  full <- lox_build_cohort(lox_read_cdm(cfg), cfg)
  cfg2 <- cfg; cfg2$cohort$sample_n <- 1000
  s <- lox_build_cohort(lox_read_cdm(cfg2), cfg2)
  expect_equal(nrow(s$cohorte), 1000)
  expect_equal(s$n_eligible, nrow(full$cohorte))
  expect_true(all(s$enf$id %in% s$cohorte$id))
  expect_equal(nrow(s$enf), nrow(full$enf[id %in% s$cohorte$id]))
  s2 <- lox_build_cohort(lox_read_cdm(cfg2), cfg2)
  expect_identical(s$cohorte$id, s2$cohorte$id)
  expect_identical(lox_config_fingerprint(cfg), lox_config_fingerprint(cfg2))
  cfg3 <- cfg; cfg3$cohort$sample_n <- 10
  expect_error(loxofed:::.lox_validate_config(cfg3))
})

test_that("la configuracion del calendario se valida y entra en la huella", {
  expect_error(lox_config(model = list(calendar = list(enabled = TRUE, knots = c(2005, 2015), pandemic_years = 2020L))), "knots")
  a <- lox_config(); b <- lox_config(model = list(calendar = list(enabled = FALSE, knots = c(2015, 2018), pandemic_years = c(2020L, 2021L))))
  expect_false(identical(lox_config_fingerprint(a), lox_config_fingerprint(b)))
  B <- lox_calendar_basis(c(2012.5, 2016.2, 2022.9), a); expect_equal(names(B), c("anio_ns1", "anio_ns2", "anio_ns3"))
})

test_that("QC de baja frente a muerte: cuenta las muertes perdidas por la regla de censura", {
  library(data.table)
  n <- 400; D <- as.Date
  pac <- data.table(id = sprintf("P%03d", 1:n), birth_year = 1940L, sex = "M", death_date = D(NA), death_res = "day", baja_date = D(NA))
  pac[1:100,  `:=`(death_date = D("2015-06-01"), baja_date = D("2015-03-01"))]   # muerte 92 dias despues de la baja: perdida
  pac[101:130, `:=`(death_date = D("2016-06-01"), baja_date = D("2016-05-15"))]  # 17 dias despues: perdida, dentro de 30 dias
  pac[131:200, `:=`(death_date = D("2017-06-01"), baja_date = D("2017-06-01"))]  # igual: cuenta
  pac[201:300, `:=`(death_date = D("2018-06-01"), baja_date = D("2019-01-01"))]  # baja posterior: cuenta
  enf <- data.table(id = pac$id, code = "HYPERTENSION", dx_date = D("2013-01-01"), dx_res = "day")
  cfg <- lox_config(paths = list(pacientes = "x", enfermedades = "y", output_dir = tempdir()))
  q <- lox_qc(list(pac = pac, enf = enf), cfg); v <- function(ch) q$value[q$check == ch]
  expect_equal(v("baja_y_exitus"), "300")
  expect_equal(v("exitus_con_baja_anterior"), "130")
  expect_equal(v("exitus_con_baja_igual"), "70")
  expect_equal(v("exitus_con_baja_posterior"), "100")
  expect_equal(v("exitus_perdidos_por_baja_anterior"), "130")
  expect_equal(v("exitus_perdidos_dentro_de_30_dias_de_la_baja"), "30")
  expect_equal(q[check == "exitus_perdidos_por_baja_anterior", flag], "WARN")
  # sin fechas de baja no se anaden filas de baja
  pac2 <- copy(pac)[, baja_date := D(NA)]
  expect_false(any(grepl("baja", lox_qc(list(pac = pac2, enf = enf), cfg)$check)))
})


test_that("SE robusto por conglomerado coincide con sandwich; contraste de edad; tablas por arco; metricas de red", {
  library(data.table)
  set.seed(1); n <- 1500; d <- data.table(id = rep(seq_len(n), each = 4)); d[, x := rbinom(.N, 1, .3)]; u <- rep(rnorm(n, 0, .8), each = 4)
  d[, y := rbinom(.N, 1, plogis(-1 + .5 * x + u))]; f <- y ~ x; fit <- glm(f, binomial(), data = d)
  V <- loxofed:::.lox_cluster_vcov(fit, f, d, d$id)
  expect_gt(sqrt(V["x", "x"]), 0); expect_gt(sqrt(V["(Intercept)", "(Intercept)"]), sqrt(vcov(fit)["(Intercept)", "(Intercept)"]))   # correlacion intra-persona: el intercepto se ensancha
  if (requireNamespace("sandwich", quietly = TRUE))
    expect_equal(sqrt(diag(V)), sqrt(diag(sandwich::vcovCL(fit, cluster = d$id, type = "HC1"))), tolerance = 0.01)
  cfg <- lox_config(paths = list(pacientes = "x", enfermedades = "y", output_dir = tempdir()))
  B <- as.matrix(lox_age_basis(c(70, 80), cfg)); b <- setNames(c(0.3, -0.1, 0.2, 0.05), colnames(B))
  Vd <- diag(0.01, 4); dimnames(Vd) <- list(names(b), names(b))
  ac <- loxofed:::.lox_age_contrast(b, Vd, cfg, "incidencia", "X", 100, 10)
  expect_equal(ac$beta, sum((B[2, ] - B[1, ]) * b)); expect_equal(ac$termino, "edad_70a80"); expect_gt(ac$se, 0)
  d2 <- data.table(A_t0 = rep(c(1L, 0L), c(300, 700)), B_t0 = 0L, B_t1 = c(rep(1L, 60), rep(0L, 240), rep(1L, 35), rep(0L, 665)))
  setnames(d2, c("A_t0", "B_t0", "B_t1"), c("A_t0", "B_t0", "B_t1"))
  at <- lox_arc_tables(d2, data.table(from = "A", to = "B"), min_cell = 10L)
  expect_equal(unlist(at[, .(a1_b1, a1_b0, a0_b1, a0_b0)]), c(a1_b1 = 60, a1_b0 = 240, a0_b1 = 35, a0_b0 = 665))
  at2 <- lox_arc_tables(d2[c(1:300, 701:1000)], data.table(from = "A", to = "B"), min_cell = 100L); expect_true(is.na(at2$a0_b1) || at2$a0_b1 == 0L)
  # metricas de red: cadena A -> B -> C y A -> C
  m <- lox_network_metrics(data.frame(from = c("A", "B", "A"), to = c("B", "C", "C"), strength = c(1, 1, 0.5)))
  expect_equal(m[node == "A", out_degree], 2L); expect_equal(m[node == "C", in_degree], 2L); expect_equal(m[node == "B", betweenness], 0)
  m2 <- lox_network_metrics(data.frame(from = c("A", "B"), to = c("B", "C")))
  expect_equal(m2[node == "B", betweenness], 1); expect_gt(m2[node == "A", harmonic_out], m2[node == "C", harmonic_out])
  # comunidades y Rand ajustado
  if (requireNamespace("igraph", quietly = TRUE)) {
    E <- data.frame(from = c("a1", "a2", "a3", "b1", "b2", "b3", "a1"), to = c("a2", "a3", "a1", "b2", "b3", "b1", "b1"), strength = c(1, 1, 1, 1, 1, 1, 0.1))
    cm <- lox_communities(E, "louvain"); expect_equal(uniqueN(cm$community), 2L)
    expect_equal(cm$community[cm$node == "a1"], cm$community[cm$node == "a2"])
  }
  expect_equal(loxofed:::.lox_ari(c(1, 1, 2, 2), c(5, 5, 9, 9)), 1)
  expect_lt(loxofed:::.lox_ari(c(1, 1, 2, 2), c(1, 2, 1, 2)), 0.1)
})

test_that("F_Baja es la columna por defecto y `~` la desactiva", {
  expect_equal(lox_config(paths = list(pacientes = "x", enfermedades = "y"))$cdm$cols$baja, "F_Baja")
  tf <- tempfile(fileext = ".yml"); writeLines(c("cdm:", "  cols:", "    baja: ~"), tf)
  expect_null(lox_config(tf, paths = list(pacientes = "x", enfermedades = "y"))$cdm$cols$baja)
})

test_that("0.7.6: nodo de un solo sexo se ajusta solo en ese sexo y sin termino Sexo", {
  set.seed(3); n <- 6000
  d2 <- data.table::data.table(id = seq_len(n), Sexo = factor(sample(c("Hombre", "Mujer"), n, TRUE)), edad_t0 = runif(n, 60, 95))
  d2[, A_t0 := rbinom(n, 1, .3)]
  d2[, P_t0 := 0L]
  d2[, P_t1 := as.integer(Sexo == "Hombre" & runif(n) < plogis(-3 + .03 * (edad_t0 - 70)))]
  d2[, A_t1 := as.integer(A_t0 == 1L | runif(n) < .1)]
  cfg <- lox_config(paths = list(pacientes = "x", enfermedades = "y", output_dir = tempdir()), model = list(min_events_node = 20))
  edges <- data.table::data.table(from = "A", to = "P")
  r <- lox_fit_layer1(d2, c("A", "P"), edges, cfg)$coefs
  expect_false(any(grepl("^Sexo", r[nodo == "P"]$termino)))
  expect_true(all(r[nodo == "P"]$restriccion == "solo_Hombre"))
  expect_true(any(grepl("^Sexo", r[nodo == "A"]$termino)))
  expect_true(all(abs(r$beta) < 10))
})
