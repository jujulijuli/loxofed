# Diagnosticos incidentes por anio natural (registro diagnostico) ---------------------------------------
# Para cada enfermedad (nodo) y anio: personas en riesgo al 1 de enero (vivas, en seguimiento y SIN la enfermedad)
# y diagnosticos incidentes del anio. Sirve para ver el infrarregistro de la pandemia (confinamiento, colapso
# asistencial, menor acceso) y para contextualizar el efecto calendario de la mortalidad. Solo agregados.

#' Diagnosticos incidentes y personas en riesgo por enfermedad y anio
#' @return data.table `node`, `anio`, `n_risk`, `n_inc`; la fila `TOTAL` suma todas las enfermedades
#' @export
lox_dx_by_year <- function(cohort, nodes, cfg) {
  yf <- cfg$cohort$year_first; yl <- cfg$cohort$year_last
  P <- cohort$cohorte[, .(id, dth = as.numeric(death_date), bj = as.numeric(baja_date))]
  b1 <- as.numeric(as.Date(sprintf("%d-01-01", yl + 1)))
  P[, end := pmin(fifelse(is.na(dth), Inf, dth), fifelse(is.na(bj), Inf, bj), b1)]
  e <- cohort$enf[code %in% nodes, .(id, code, d = as.numeric(dx_date))]
  e[P, end := i.end, on = "id"]; e <- e[!is.na(end)]
  out <- vector("list", yl - yf + 1)
  for (y in yf:yl) {
    j1 <- as.numeric(as.Date(sprintf("%d-01-01", y))); j2 <- as.numeric(as.Date(sprintf("%d-01-01", y + 1)))
    alive <- sum(P$end > j1)
    prev <- if (y == yf) e[d <= j1 & end > j1, .N, by = code] else e[d < j1 & end > j1, .N, by = code]
    inc  <- if (y == yf) e[d > j1 & d < j2 & d <= end, .N, by = code] else e[d >= j1 & d < j2 & d <= end, .N, by = code]
    o <- data.table(node = nodes, anio = y)
    o[, n_risk := alive - prev$N[match(node, prev$code)]]; o[is.na(n_risk), n_risk := alive]
    o[, n_inc := inc$N[match(node, inc$code)]]; o[is.na(n_inc), n_inc := 0L]
    out[[y - yf + 1]] <- rbind(o, data.table(node = "TOTAL", anio = y, n_risk = sum(o$n_risk), n_inc = sum(o$n_inc)))
  }
  rbindlist(out)[order(node, anio)]
}

#' Diagnosticos observados frente a los esperados por la tendencia previa a la pandemia
#'
#' Por enfermedad (y `TOTAL`): Poisson con offset `log(n_risk)` y tendencia lineal del anio, ajustado SOLO con los anios
#' previos a los anios pandemicos (`pandemia == FALSE` y anterior al primero), extrapolado a cada anio pandemico.
#' Cociente observados / esperados con IC95 % (Poisson de lo observado + incertidumbre de la extrapolacion).
#' @param dxy tabla de `lox_dx_by_year()` con la columna logica `pandemia`
#' @export
lox_dx_pandemic <- function(dxy, min_obs = 30L) {
  dxy <- as.data.table(dxy); if (!"pandemia" %in% names(dxy) || !any(dxy$pandemia)) return(NULL)
  y0 <- min(dxy$anio[dxy$pandemia]); res <- list()
  for (nd in unique(dxy$node)) {
    d <- dxy[node == nd & !is.na(n_inc)]; pre <- d[anio < y0]; pan <- d[pandemia == TRUE]
    if (nrow(pre) < 3 || !nrow(pan) || sum(pre$n_inc) < 20) next
    fit <- tryCatch(suppressWarnings(stats::glm(n_inc ~ anio + offset(log(n_risk)), family = stats::poisson(), data = pre)), error = function(e) NULL)
    if (is.null(fit)) next
    pr <- stats::predict(fit, newdata = pan, se.fit = TRUE)
    mu <- exp(pr$fit)                                              # predict() suma el offset log(n_risk) de `newdata`
    O <- pan$n_inc
    r <- data.table(node = nd, anio = pan$anio, O = O, E = mu, ratio = O / mu)
    r[, se := sqrt(1 / pmax(O, 1) + pr$se.fit^2)]; r[, `:=`(lo = ratio * exp(-1.96 * se), hi = ratio * exp(1.96 * se))]
    res[[nd]] <- r[O >= 0]
  }
  if (!length(res)) return(NULL)
  r <- rbindlist(res); r[, E := round(E, 1)]; r[]
}
