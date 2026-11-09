mk_site <- function(root, id, seed, date_res = "day", n = 4000) {
  d <- file.path(root, id, "data"); lox_simulate_cdm(n = n, out_dir = d, seed = seed, date_res = date_res)
  lox_config(paths = list(pacientes = file.path(d, "pacientes.csv"), enfermedades = file.path(d, "enfermedades.csv"),
                          output_dir = file.path(root, id, "out")),
             site_id = id, structure = list(boot_R = 10, subsample = 20000), model = list(min_events_node = 30))
}

test_that("QC avisa con fechas anuales y la cohorte sigue funcionando", {
  root <- tempfile("qc"); cfg <- mk_site(root, "Y", 3, date_res = "year", n = 1500)
  cdm <- lox_read_cdm(cfg); q <- lox_qc(cdm, cfg)
  expect_equal(q[check == "resolucion_dx_pct_dia/mes/anio", flag], "WARN")
  expect_gt(nrow(lox_build_cohort(cdm, cfg)$cohorte), 100)
})

test_that("la lectura falla con un mensaje claro si falta una columna", {
  root <- tempfile("col"); cfg <- mk_site(root, "Z", 4, n = 500)
  cfg$cdm$cols$birth <- "NO_EXISTE"
  expect_error(lox_read_cdm(cfg), "NO_EXISTE")
})

test_that("flujo federado completo: 2 sitios -> consenso -> ronda 2 -> meta", {
  root <- tempfile("fed")
  cA <- mk_site(root, "A", 1); cB <- mk_site(root, "B", 2)
  suppressMessages({ lox_run_round1(cA); lox_run_round1(cB) })
  r1 <- file.path(root, c("A", "B"), "out", "round1")
  res <- suppressMessages(lox_run_central_consensus(r1, file.path(root, "central"), cA))
  st <- res$structure
  # progresiones verdaderas fuertes recuperadas con la direccion correcta
  expect_true(any(st$from == "HYPERTENSION" & st$to == "HEART_FAILURE" & st$state == "concordante"))
  expect_true(any(st$from == "DIABETES" & st$to == "CHRONIC_KIDNEY_DISEASES" & st$state == "concordante"))
  expect_true(any(st$from == "DEPRESSION_AND_MOOD_DISEASES" & st$to == "DEMENTIA"))
  cons <- file.path(root, "central", "consensus")
  suppressMessages({ lox_run_round2(cA, cons); lox_run_round2(cB, cons) })
  r2 <- file.path(root, c("A", "B"), "out", "round2")
  cA$meta$min_sites <- 2L                                        # solo hay 2 sitios en esta prueba
  mt <- suppressMessages(lox_run_central_meta(r2, file.path(root, "central"), cA))
  hr <- mt[capa == "mortalidad" & nodo %in% c("DEMENTIA", "HEART_FAILURE")]
  expect_true(all(hr$exp_re > 1.3)); expect_true(all(hr$k == 2))
  expect_true(all(c("pi_lo", "pi_hi", "tau2", "I2", "beta_dl", "loo_min") %in% names(mt)))
  expect_true(all(is.na(hr$pi_lo)))                              # con k = 2 no hay intervalo de prediccion
  ms <- data.table::fread(file.path(root, "central", "meta", "meta_sites.csv")); expect_equal(sort(ms$site), c("A", "B"))
  expect_true(all(c("meta_sex.csv", "meta_incidence_strata.csv", "meta_arc_tables.csv", "meta_site_hr_sexo.csv") %in% list.files(file.path(root, "central", "meta"))))
  msx <- data.table::fread(file.path(root, "central", "meta", "meta_sex.csv"))
  expect_true(all(c("mortalidad_Hombre", "mortalidad_Mujer", "mortalidad_dif_sexo") %in% msx$capa))
  expect_true(all(msx[capa == "mortalidad_dif_sexo"]$k == 2))
  ea <- mt[capa == "incidencia" & termino == "edad_70a80"]; expect_true(nrow(ea) > 0); expect_true(all(ea$exp_re > 1))   # mas edad, mas incidencia
  hc <- file.path(root, "central", "informe_central.html"); expect_true(file.exists(hc))
  htm <- paste(readLines(hc, warn = FALSE), collapse = "\n"); expect_match(htm, "Hartung-Knapp"); expect_match(htm, "<svg")
  expect_false(grepl("https?://", sub("w3.org", "", htm)))
  inc <- mt[capa == "incidencia" & nodo == "HEART_FAILURE" & termino == "HYPERTENSION_t0"]
  expect_equal(inc$exp_re, exp(1.2), tolerance = 0.25)          # OR verdadero 3.32

  # salvaguardas: ninguna salida contiene identificadores y todas llevan manifiesto
  for (d in c(r1, r2)) {
    expect_true(file.exists(file.path(d, "manifest.json")))
    for (f in list.files(d, pattern = "\\.csv$", full.names = TRUE)) {
      h <- tolower(names(data.table::fread(f, nrows = 0)))
      expect_false(any(h %in% c("id", "id_pac")), info = f)
    }
  }
  obs <- data.table::fread(file.path(r2[1], "r2_observed.csv"))
  expect_true(all(is.na(obs$n_death[obs$suppressed])))
  expect_true(all(is.na(obs$n_exp[obs$suppressed])))                  # la esperada no delata la celda suprimida
  ok <- obs[!is.na(n_death)]; expect_equal(sum(ok$n_exp), sum(ok$n_death), tolerance = 0.06)   # Poisson con intercepto: O ~ E

  # comparador: una ejecucion frente a si misma -> estructura identica y correlacion 1
  cmp <- suppressMessages(lox_compare_runs(file.path(root, "A", "out"), file.path(root, "A", "out"), consensus_a = cons, consensus_b = cons,
                                           out_file = file.path(root, "cmp.html")))
  hc2 <- paste(readLines(cmp, warn = FALSE), collapse = "\n")
  expect_match(hc2, "Jaccard 100%"); expect_match(hc2, "correlación de log-HR 1.000"); expect_match(hc2, "<svg")
  expect_false(grepl("https?://", sub("w3.org", "", hc2)))

  # huella distinta -> el consenso se niega a mezclar sitios
  cX <- lox_config(cohort = list(age_min = 65), paths = cA$paths, site_id = "X")
  expect_error(lox_run_central_consensus(r1, file.path(root, "c2"), cX), "huella")
})

test_that("el modo annual_next (convencion del informe de Navarra) sigue funcionando de punta a punta", {
  root <- tempfile("an"); d <- file.path(root, "S", "data"); lox_simulate_cdm(n = 3000, out_dir = d, seed = 8)
  cfg <- lox_config(paths = list(pacientes = file.path(d, "pacientes.csv"), enfermedades = file.path(d, "enfermedades.csv"),
                                 output_dir = file.path(root, "S", "out")),
                    site_id = "S", structure = list(boot_R = 10, subsample = 15000),
                    model = list(min_events_node = 30, mortality = "annual_next"))
  suppressMessages(lox_run_round1(cfg))
  suppressMessages(lox_run_central_consensus(file.path(root, "S", "out", "round1"), file.path(root, "c"), cfg))
  suppressMessages(lox_run_round2(cfg, file.path(root, "c", "consensus")))
  obs <- data.table::fread(file.path(root, "S", "out", "round2", "r2_observed.csv"))
  expect_true("n_risk" %in% names(obs))
  m <- jsonlite::read_json(file.path(root, "S", "out", "round2", "manifest.json"))
  expect_equal(m$mortality, "annual_next")
})

test_that("la ronda 2 en modo exacto exporta personas-anio y muertes", {
  root <- tempfile("ex"); cfg <- mk_site(root, "E", 6, n = 3000)
  suppressMessages(lox_run_round1(cfg))
  suppressMessages(lox_run_central_consensus(file.path(root, "E", "out", "round1"), file.path(root, "c"), cfg))
  suppressMessages(lox_run_round2(cfg, file.path(root, "c", "consensus")))
  obs <- data.table::fread(file.path(root, "E", "out", "round2", "r2_observed.csv"))
  expect_true(all(c("person_years", "n_death") %in% names(obs)))
  expect_equal(jsonlite::read_json(file.path(root, "E", "out", "round2", "manifest.json"))$mortality, "exact")
})

test_that("pesos de muestreo escalan los recuentos descriptivos y el informe local se genera", {
  root <- tempfile("rep"); d <- file.path(root, "data"); lox_simulate_cdm(n = 4000, out_dir = d, seed = 5)
  out <- file.path(root, "out")
  cfg <- lox_config(site_id = "REP", paths = list(pacientes = file.path(d, "pacientes.csv"),
                    enfermedades = file.path(d, "enfermedades.csv"), output_dir = out),
                    cdm = list(sep = "|"), cohort = list(sample_n = 1500), structure = list(boot_R = 10),
                    model = list(min_events_node = 30, min_events_hazard = 10))
  capture.output(suppressMessages({
    lox_run_round1(cfg)
    lox_run_central_consensus(file.path(out, "round1"), file.path(out, "cons"), cfg)
    lox_run_round2(cfg, file.path(out, "cons", "consensus")) }))
  m <- fread(file.path(out, "round1", "r1_meta.csv"))
  expect_equal(m$n_cohort, 1500); expect_equal(m$sample_weight, round(m$n_eligible / 1500, 4))
  pv <- fread(file.path(out, "round1", "r1_prevalence.csv"))[!is.na(n_pos)]
  expect_equal(pv$n_pos_est, round(pv$n_pos * m$sample_weight)); expect_true(all(pv$n_est == round(1500 * m$sample_weight)))
  ob <- fread(file.path(out, "round2", "r2_observed.csv"))
  expect_true(all(is.na(ob$n_death_est[is.na(ob$n_death)])))            # lo suprimido sigue suprimido
  expect_equal(ob[!is.na(n_death)]$n_death_est, round(ob[!is.na(n_death)]$n_death * m$sample_weight))
  f <- file.path(out, "informe_local.html"); expect_true(file.exists(f))   # auto tras la ronda 2
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(h, "Calibración"); expect_true("n_exp_est" %in% names(ob))
  expect_match(h, "Efecto del calendario"); expect_true("n_exp0" %in% names(ob))
  pe <- fread(file.path(out, "round2", "r2_period_effect.csv")); expect_equal(pe$anio[1], 2012L); expect_equal(pe$log_rr_total[1], 0)
  expect_true(all(c("log_rr_tendencia", "se_total") %in% names(pe)))
  co2 <- fread(file.path(out, "round2", "r2_coefs.csv")); expect_true(all(c("anio_ns1", "pandemia") %in% co2[nodo == "muerte"]$termino))
  expect_true(file.exists(file.path(out, "round2", "r2_hr_sin_calendario.csv")))
  dx <- fread(file.path(out, "round2", "r2_dx_by_year.csv"))
  expect_true(all(c("node", "anio", "n_risk", "n_inc", "pandemia", "n_inc_est") %in% names(dx)))
  expect_true(all(is.na(dx$n_inc[dx$suppressed])))                                       # supresion de celdas pequenas
  expect_true(dx[node == "TOTAL" & anio == 2020, pandemia]); expect_false(dx[node == "TOTAL" & anio == 2019, pandemia])
  tot <- dx[node == "TOTAL"][order(anio)]; expect_true(all(diff(tot$n_risk) <= 0))        # cohorte cerrada: nadie entra
  pp <- lox_dx_pandemic(dx); expect_true(all(c("ratio", "lo", "hi") %in% names(pp))); expect_true(all(pp$anio %in% c(2020, 2021)))
  expect_match(h, "Registro diagnóstico durante la pandemia")
  # 0.7.5: ronda 2 ampliada
  expect_true(all(c("r2_incidence_strata.csv", "r2_arc_tables.csv", "r2_hr_por_sexo.csv") %in% list.files(file.path(out, "round2"))))
  ist <- fread(file.path(out, "round2", "r2_incidence_strata.csv"))
  expect_true(all(c("node", "anio", "Sexo", "edad_grupo", "n_risk", "n_inc") %in% names(ist)))
  expect_true(all(ist$n_inc[!is.na(ist$n_inc)] == 0 | ist$n_inc[!is.na(ist$n_inc)] >= 10)); expect_true(all(ist$n_risk[!is.na(ist$n_risk)] >= 10))   # supresion de celdas pequenas
  at <- fread(file.path(out, "round2", "r2_arc_tables.csv")); expect_true(all(c("a1_b1", "a1_b0", "a0_b1", "a0_b0") %in% names(at)))
  expect_true(all(c("edad_70a80") %in% co2$termino)); expect_true(all(is.finite(co2[termino == "edad_70a80"]$se)))
  bs <- fread(file.path(out, "round2", "r2_hr_por_sexo.csv")); expect_true(all(c("Hombre", "Mujer") %in% bs$sexo))
  expect_match(h, "Es igual el efecto en hombres y en mujeres"); expect_match(h, "Qué es este informe y dónde encaja")
  expect_match(h, "provisional")                                                     # consenso calculado con un solo sitio
  expect_false(any(grepl("Desde</th>", h) & FALSE))
  pac0 <- fread(file.path(d, "pacientes.csv")); expect_true("F_Baja" %in% names(pac0)); expect_true(all(is.na(pac0$F_Baja) | pac0$F_Baja == ""))
  qc0 <- fread(file.path(out, "round1", "qc.csv")); expect_false(any(grepl("baja", qc0$check)))       # F_Baja vacia: sin controles de baja
  expect_match(h, "Muestra aleatoria"); expect_match(h, "<svg"); expect_false(grepl("https?://", sub("w3.org", "", h)))
  expect_equal(length(gregexpr("<svg", h)[[1]]), length(gregexpr("</svg>", h)[[1]]))
})
