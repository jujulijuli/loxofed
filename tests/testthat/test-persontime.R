# Cohorte construida a mano para comprobar el tiempo-persona caso por caso -------------------------------
mk_cdm <- function() {
  D <- as.Date
  pac <- data.table(
    id = paste0("p", 1:7), birth_year = 1940L, sex_raw = "Mujer", sex = "Mujer",
    death_date = D(c("2012-01-01", "2015-06-15", NA, "2018-02-01", "2022-12-31", "2023-03-01", NA)),
    death_res = "day",
    baja_date = D(c(NA, NA, "2016-05-01", "2016-05-01", NA, NA, NA)))
  # todas con 2 cronicas prevalentes; p2 desarrolla IC el 2014-03-10
  base <- CJ(id = pac$id, code = c("HYPERTENSION", "DIABETES"))[, `:=`(dx_date = D("2005-01-01"), dx_res = "day")]
  inc  <- data.table(id = "p2", code = "HEART_FAILURE", dx_date = D("2014-03-10"), dx_res = "day")
  list(pac = pac, enf = rbind(base, inc))
}
cfg0 <- function(model = list()) lox_config(model = c(list(min_events_hazard = 1), model))

test_that("tiempo-persona: casos a mano (muerte basal, baja, fin de seguimiento, diagnostico incidente)", {
  cfg <- cfg0(); co <- lox_build_cohort(mk_cdm(), cfg)
  nodes <- c("HYPERTENSION", "DIABETES", "HEART_FAILURE")
  pt <- lox_build_persontime(co, nodes, cfg)
  idp <- function(i) which(co$cohorte$id == i)
  y <- 365.2425; b0 <- as.numeric(as.Date("2012-01-01"))
  one <- function(i) pt[pid == idp(i)][order(t0)]

  p1 <- one("p1")                                  # muerte el propio dia basal: no se pierde
  expect_equal(sum(p1$event), 1); expect_gt(sum(p1$exposure), 0)

  p2 <- one("p2")                                  # IC el 2014-03-10, muerte 2015-06-15
  expect_equal(nrow(p2), 5)                        # cortes: 2013-01-01, 2014-01-01, dx, 2015-01-01
  expect_equal(sum(p2$exposure), (as.numeric(as.Date("2015-06-15")) - b0) / y)
  expect_equal(p2$event, c(0, 0, 0, 0, 1))
  expect_equal(p2$HEART_FAILURE_t0, c(0, 0, 0, 1, 1))
  expect_equal(p2$t0[4], as.numeric(as.Date("2014-03-10")))

  p3 <- one("p3")                                  # baja 2016-05-01: censura, no muerte
  expect_equal(sum(p3$event), 0)
  expect_equal(max(p3$t1), as.numeric(as.Date("2016-05-01")))

  expect_equal(sum(one("p4")$event), 0)            # muere DESPUES de la baja: no se observa
  expect_equal(sum(one("p5")$event), 1)            # 2022-12-31 sigue dentro del seguimiento
  expect_equal(sum(one("p6")$event), 0)            # 2023-03-01: fuera, se censura a 2023-01-01
  expect_equal(max(one("p6")$t1), as.numeric(as.Date("2023-01-01")))
})

test_that("tiempo-persona: cuenta todas las muertes, es contiguo y los estados son absorbentes", {
  d <- tempfile("pt"); lox_simulate_cdm(n = 3000, out_dir = d, seed = 5)
  cfg <- lox_config(paths = list(pacientes = file.path(d, "pacientes.csv"),
                                 enfermedades = file.path(d, "enfermedades.csv")))
  cdm <- lox_read_cdm(cfg); co <- lox_build_cohort(cdm, cfg)
  nodes <- c("HYPERTENSION", "HEART_FAILURE", "DEMENTIA", "DIABETES")
  pt <- lox_build_persontime(co, nodes, cfg); setorder(pt, pid, t0)

  n_muertes <- co$cohorte[!is.na(death_date) & death_date < as.Date("2023-01-01"), .N]
  expect_equal(sum(pt$event), n_muertes)                              # ninguna muerte sin contar
  fin <- pt[, .(first = t0[1], last = t1[.N], py = sum(exposure)), by = pid]
  expect_true(all(fin$first == as.numeric(as.Date("2012-01-01"))))
  expect_equal(sum(fin$py), sum(fin$last - fin$first) / 365.2425)
  pt[, gap := t0 - shift(t1), by = pid]; expect_true(all(is.na(pt$gap) | pt$gap == 0))   # sin huecos ni solapes
  for (v in paste0(nodes, "_t0")) {
    pt[, dif := get(v) - shift(get(v)), by = pid]; expect_true(all(is.na(pt$dif) | pt$dif >= 0), info = v)
  }
  # el tiempo expuesto solo crece al partir mas fino (misma exposicion total)
  c6 <- lox_config(paths = cfg$paths, model = list(time_split_months = 6))
  expect_equal(sum(lox_build_persontime(co, nodes, c6)$exposure), sum(pt$exposure))
})

test_that("dx_lag_days retrasa la exposicion; con un retardo enorme solo cuentan los prevalentes", {
  cfg <- cfg0(list(dx_lag_days = 100000)); co <- lox_build_cohort(mk_cdm(), cfg)
  pt <- lox_build_persontime(co, c("HEART_FAILURE"), cfg)
  expect_true(all(pt$HEART_FAILURE_t0 == 0))
  cfg1 <- cfg0(list(dx_lag_days = 30)); pt1 <- lox_build_persontime(co, "HEART_FAILURE", cfg1)
  p2 <- pt1[pid == which(co$cohorte$id == "p2")][order(t0)]
  expect_equal(p2$t0[p2$HEART_FAILURE_t0 == 1][1], as.numeric(as.Date("2014-03-10")) + 30)
})

test_that("la mortalidad exacta recupera HR mayores que 1 y menos atenuados que la convencion anual", {
  d <- tempfile("ex"); lox_simulate_cdm(n = 9000, out_dir = d, seed = 21)
  cfg <- lox_config(paths = list(pacientes = file.path(d, "pacientes.csv"),
                                 enfermedades = file.path(d, "enfermedades.csv")))
  cdm <- lox_read_cdm(cfg); co <- lox_build_cohort(cdm, cfg); arr <- lox_build_states(co, cfg)
  nodes <- lox_prevalence(arr)[prev >= 0.01, node]
  pt <- lox_build_persontime(co, nodes, cfg); d2 <- lox_build_transitions(arr, co, nodes, cfg)
  hx <- lox_fit_layer2_exact(pt, nodes, cfg); ha <- lox_fit_layer2(d2, nodes, cfg)
  g <- function(h, estado) exp(h$coefs[nodo == estado, beta])
  expect_true(hx$converged)
  for (e in c("HEART_FAILURE", "DEMENTIA")) { expect_gt(g(hx, e), 1.5); expect_gt(g(hx, e), g(ha, e)) }
  expect_gt(sum(pt$event), sum(d2$death_next))                         # la convencion anual deja muertes fuera
})
