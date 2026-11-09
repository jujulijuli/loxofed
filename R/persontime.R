# Mortalidad con fechas exactas: tiempo-persona con covariables que cambian en el tiempo ---------------

#' Tabla de tiempo-persona (una fila por tramo) con estados que cambian en cada diagnostico
#'
#' El seguimiento de cada persona va del 1 de enero de `year_first` hasta la primera de: muerte,
#' baja administrativa y 1 de enero de `year_last + 1`. Se parte en cada diagnostico incidente
#' (de las condiciones `nodes`) y cada `time_split_months` meses desde el origen (para que la edad
#' avance). El estado de cada tramo es el vigente al inicio del tramo; el desenlace es la muerte
#' al final del tramo. Asi se cuentan TODAS las muertes del seguimiento y no hay cortes anuales
#' en la exposicion.
#'
#' `dx_lag_days` retrasa el momento en que un diagnostico cuenta como exposicion (sensibilidad
#' para diagnosticos codificados junto a la muerte). Los diagnosticos prevalentes
#' (<= 1 de enero de `year_first`) cuentan siempre.
#'
#' @return data.table: `pid`, `Sexo`, `t0`, `t1` (dias desde 1970), `exposure` (anios), `event`,
#'   `edad_t0` (edad en el punto medio del tramo), `anio` (anio del inicio del tramo) y `<nodo>_t0`.
#' @export
lox_build_persontime <- function(cohort, nodes, cfg) {
  yf <- cfg$cohort$year_first; yl <- cfg$cohort$year_last; m <- cfg$model
  base <- as.Date(sprintf("%d-01-01", yf)); fu_end <- as.Date(sprintf("%d-01-01", yl + 1))
  b0 <- as.numeric(base); b1 <- as.numeric(fu_end)

  P <- copy(cohort$cohorte[, .(id, Sexo, edad_base, death_date, baja_date)])
  P[, pid := .I]
  P[, dth := fifelse(is.na(death_date), Inf, as.numeric(death_date))]
  P[, bj  := fifelse(is.na(baja_date),  Inf, as.numeric(baja_date))]
  P[, event := as.integer(dth < b1 & dth <= bj)]                 # muerte dentro del seguimiento y antes de la baja
  P[, end := pmin(dth, bj, b1)]
  P[event == 1L & end <= b0, end := b0 + 0.5]                    # muerte el propio dia basal: tramo minimo, no se pierde
  P <- P[end > b0]

  # cortes calendario (la edad avanza) y cortes por diagnostico incidente
  cc  <- as.numeric(seq(base, fu_end, by = sprintf("%d months", m$time_split_months)))
  cc  <- cc[cc > b0 & cc < b1]
  cal <- CJ(pid = P$pid, t = cc)
  cal <- cal[P[, .(pid, end_i = end)], on = "pid", nomatch = NULL][t < end_i, .(pid, t)]

  e <- cohort$enf[code %in% nodes, .(id, code, dx_date)]
  e[, pid := P$pid[match(id, P$id)]]
  e <- e[!is.na(pid)]
  e[, teff := fifelse(dx_date <= base, -Inf, as.numeric(dx_date) + m$dx_lag_days)]
  e[P, end_i := i.end, on = "pid"]
  dxc <- e[teff > b0 & teff < end_i, .(pid, t = teff)]

  cu <- unique(rbind(cal, dxc, P[, .(pid, t = b0)], P[, .(pid, t = end)]))
  setorder(cu, pid, t)
  cu[, t1 := shift(t, 1L, type = "lead"), by = pid]
  iv <- cu[!is.na(t1)]; setnames(iv, "t", "t0")
  rm(cu, cal, dxc); invisible(gc())

  iv[P, `:=`(event_p = i.event, end_p = i.end, Sexo = i.Sexo, edad_base = i.edad_base), on = "pid"]
  iv[, event    := as.integer(event_p == 1L & t1 == end_p)]
  iv[, exposure := (t1 - t0) / 365.2425]
  # edad en el punto medio del tramo, suponiendo nacimiento a mitad del anio de nacimiento
  iv[, edad_t0  := 1970 + ((t0 + t1) / 2) / 365.2425 - (yf - edad_base + 0.5)]
  iv[, anio     := data.table::year(data.table::as.IDate(t0))]
  iv[, cal_t0   := 1970 + ((t0 + t1) / 2) / 365.2425]            # tiempo calendario en el punto medio del tramo (anio decimal)

  nP <- max(P$pid)
  for (nd in nodes) {
    w <- rep(NA_real_, nP); en <- e[code == nd]; w[en$pid] <- en$teff
    set(iv, j = paste0(nd, "_t0"), value = as.integer(!is.na(w[iv$pid]) & w[iv$pid] <= iv$t0))
  }
  iv[, c("event_p", "end_p", "edad_base") := NULL]
  iv[]
}

#' Base de spline natural del tiempo calendario con nudos FIJOS (comun a todos los sitios)
#' @return data.table con columnas `anio_ns1..anio_nsK`
#' @export
lox_calendar_basis <- function(cal, cfg) {
  B <- splines::ns(cal, knots = cfg$model$calendar$knots, Boundary.knots = c(cfg$cohort$year_first, cfg$cohort$year_last + 1))
  B <- as.data.table(unclass(B)); setnames(B, paste0("anio_ns", seq_len(ncol(B)))); B
}

#' Capa 2 con fechas exactas: modelo de Poisson con offset sobre el tiempo-persona (exp(beta) = HR)
#'
#' Mismas covariables que `lox_fit_layer2` (spline de edad con nudos fijos, sexo y estados con
#' >= `model$min_events_hazard` muertes entre expuestos), pero con estados que cambian en el tiempo
#' y exposicion exacta. Los coeficientes son log-HR de tasas instantaneas y son meta-analizables.
#' @export
lox_fit_layer2_exact <- function(pt, nodes, cfg) {
  B <- lox_age_basis(pt$edad_t0, cfg); bn <- names(B)
  v0 <- paste0(nodes, "_t0")
  ev <- vapply(v0, function(v) sum(pt$event[pt[[v]] == 1L]), numeric(1))
  keep <- v0[ev >= cfg$model$min_events_hazard]
  cal <- cfg$model$calendar; use_cal <- !is.null(cal) && isTRUE(cal$enabled) && "cal_t0" %in% names(pt)
  dat <- cbind(pt[, c("event", "exposure", "Sexo", keep), with = FALSE], B)
  dat[, loff := log(exposure)]
  cn <- character(); pand <- NULL
  if (use_cal) {
    Bc <- lox_calendar_basis(pt$cal_t0, cfg); cn <- names(Bc); dat <- cbind(dat, Bc)
    py <- as.integer(cal$pandemic_years)
    if (length(py) && any(pt$anio %in% py)) { dat[, pandemia := as.integer(pt$anio %in% py)]; pand <- "pandemia" }
  }
  ft <- function(terms) stats::glm(stats::reformulate(c(terms, "offset(loff)"), "event"), family = stats::poisson(), data = dat,
                                   model = FALSE, y = FALSE, control = stats::glm.control(maxit = 100))
  fit <- ft(c(bn, "Sexo", keep, cn, pand))
  r <- .lox_glm_coefs(fit, "mortalidad", "muerte", nrow(dat), sum(dat$event))
  r <- r[!is.na(beta) & termino != "(Intercept)"]
  r[, nodo := ifelse(termino %in% keep, sub("_t0$", "", termino), "muerte")]
  ac <- tryCatch(.lox_age_contrast(stats::coef(fit), stats::vcov(fit), cfg, "mortalidad", "muerte", nrow(dat), sum(dat$event)), error = function(e) NULL)
  if (!is.null(ac)) r <- rbind(r, ac, fill = TRUE)
  out <- list(coefs = r, excluded = sub("_t0$", "", setdiff(v0, keep)), converged = isTRUE(fit$converged),
              expected = as.numeric(fit$fitted.values))   # muertes esperadas por tramo (misma fila que `pt`)
  # --- HR por sexo: mismo modelo ajustado por separado en hombres y mujeres (para el contraste por sexo)
  bs <- lapply(unique(dat$Sexo), function(sx) {
    d <- dat[Sexo == sx]; ks <- keep[vapply(keep, function(v) sum(d$event[d[[v]] == 1L]) >= cfg$model$min_events_hazard, logical(1))]
    if (!length(ks)) return(NULL)
    f <- tryCatch(suppressWarnings(stats::glm(stats::reformulate(c(bn, ks, cn, pand), "event"), family = stats::poisson(), data = d, offset = loff,
                                              model = FALSE, y = FALSE, control = stats::glm.control(maxit = 100))), error = function(e) NULL)
    if (is.null(f) || !isTRUE(f$converged)) return(NULL)
    rr <- .lox_glm_coefs(f, "mortalidad", "muerte", nrow(d), sum(d$event))[termino %in% ks]
    rr[, `:=`(nodo = sub("_t0$", "", termino), sexo = sx, n_events = vapply(termino, function(v) sum(d$event[d[[v]] == 1L]), numeric(1)))]
    rr[, .(sexo, nodo, beta, se, n_events)]
  })
  out$by_sex <- rbindlist(bs)
  if (use_cal) {
    # --- curva del efecto calendario (tasa relativa al primer anio) con IC por el metodo delta
    V <- stats::vcov(fit); yrs <- seq(cfg$cohort$year_first, cfg$cohort$year_last)
    Bg <- as.matrix(lox_calendar_basis(c(cfg$cohort$year_first + 0.5, yrs + 0.5), cfg)); L <- sweep(Bg[-1, , drop = FALSE], 2, Bg[1, ], "-")
    cf <- stats::coef(fit); ok <- cn[cn %in% names(cf) & !is.na(cf[cn])]
    Lt <- L[, match(ok, cn), drop = FALSE]; bt <- cf[ok]; Vt <- V[ok, ok, drop = FALSE]
    tend <- drop(Lt %*% bt); se_t <- sqrt(rowSums((Lt %*% Vt) * Lt))
    has_p <- !is.null(pand) && !is.na(cf[pand]); ip <- as.numeric(yrs %in% as.integer(cal$pandemic_years))
    Lall <- if (has_p) cbind(Lt, ip) else Lt; nm <- c(ok, if (has_p) pand)
    ball <- cf[nm]; Vall <- V[nm, nm, drop = FALSE]
    tot <- drop(Lall %*% ball); se_a <- sqrt(rowSums((Lall %*% Vall) * Lall))
    out$period <- data.table(anio = yrs, log_rr_tendencia = tend, se_tendencia = se_t, log_rr_total = tot, se_total = se_a)
    # --- modelo SIN calendario: cuanto cambian esperadas y HR al no ajustar por periodo
    fit0 <- ft(c(bn, "Sexo", keep))
    out$expected0 <- as.numeric(fit0$fitted.values)
    r0 <- .lox_glm_coefs(fit0, "mortalidad", "muerte", nrow(dat), sum(dat$event))
    r0 <- r0[!is.na(beta) & termino %in% keep]; r0[, nodo := sub("_t0$", "", termino)]
    out$coefs0 <- r0[, .(nodo, beta, se)]
  }
  out
}

# Personas-anio y muertes observadas por sexo x grupo de edad x anio (agregado) a partir del tiempo-persona
.lox_observed_exact <- function(pt, expected = NULL, expected0 = NULL) {
  if (!is.null(expected) && length(expected) == nrow(pt)) pt <- cbind(pt[, .(Sexo, edad_t0, anio, exposure, event)], expected = expected)
  if (!is.null(expected0) && "expected" %in% names(pt)) pt[, expected0 := expected0]
  g <- cut(pt$edad_t0, c(0, 50, 60, 70, 80, 90, Inf), right = FALSE,
           labels = c("<50", "50-59", "60-69", "70-79", "80-89", "90+"))
  o <- if ("expected" %in% names(pt))
    pt[, c(list(person_years = sum(exposure), n_death = sum(event), n_exp = sum(expected)),
           if ("expected0" %in% names(pt)) list(n_exp0 = sum(expected0))),
       by = .(Sexo, grupo_edad = g, anio = anio)][order(Sexo, grupo_edad, anio)]
  else pt[, .(person_years = sum(exposure), n_death = sum(event)),
          by = .(Sexo, grupo_edad = g, anio = anio)][order(Sexo, grupo_edad, anio)]
  o[, person_years := round(person_years, 1)]
  if ("n_exp" %in% names(o)) o[, n_exp := round(n_exp, 1)]
  if ("n_exp0" %in% names(o)) o[, n_exp0 := round(n_exp0, 1)]
  o[]
}
