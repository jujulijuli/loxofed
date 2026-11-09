# Ronda 2: ajuste de las dos capas sobre la estructura consenso --------------------------------------

.lox_glm_coefs <- function(fit, capa, nodo, n, n_events, V = NULL) {
  s <- summary(fit)$coefficients
  se <- s[, 2]
  if (!is.null(V)) { sr <- sqrt(diag(V)); se[names(sr)] <- sr[names(sr)] }       # error estandar robusto por persona
  data.table(capa = capa, nodo = nodo, termino = rownames(s), beta = s[, 1], se = se,
             n = n, n_events = n_events, converged = isTRUE(fit$converged),
             flag = if (!isTRUE(fit$converged)) "no_convergio" else if (any(s[, 2] > 20)) "posible_separacion" else "")
}

#' Capa 1: CPD logistica de incidencia de cada nodo dado sus padres validados + edad + sexo
#'
#' `P(B_t1 = 1 | padres_t0, spline(edad), sexo)` entre quienes NO tienen B en t0 (la persistencia
#' es determinista por absorbencia). Nodos con menos de `model$min_events_node` incidentes no se
#' ajustan (regla en frecuencias absolutas, no en prevalencia). Si casi todos los casos (>= 99% o
#' todos salvo 5) son de un solo sexo, se ajusta solo en ese sexo y no se estima el efecto del sexo
#' (columna `restriccion` = `solo_<sexo>`).
#' @param edges data.table(from, to): arcos de la red causal consenso
#' @export
lox_fit_layer1 <- function(d2, nodes, edges, cfg) {
  B <- lox_age_basis(d2$edad_t0, cfg); bn <- names(B)
  res <- list(); skipped <- list()
  for (nd in nodes) {
    c0 <- paste0(nd, "_t0"); c1 <- paste0(nd, "_t1")
    idx <- which(d2[[c0]] == 0L)
    ev <- sum(d2[[c1]][idx])
    if (ev < cfg$model$min_events_node) { skipped[[nd]] <- data.table(nodo = nd, motivo = "pocos_eventos", n_events = ev); next }
    pa <- intersect(edges$from[edges$to == nd], nodes); pa <- setdiff(pa, nd)
    pcols <- if (length(pa)) paste0(pa, "_t0") else character(0)   # paste0(character(0), "_t0") = "_t0"
    # 0.7.6: enfermedades plausiblemente de un solo sexo (p. ej. prostata) -> se ajusta solo en ese sexo y
    # se omite el termino Sexo (con 0 casos en el otro sexo habria separacion perfecta, beta ~ 18)
    sx <- d2$Sexo[idx]; ev_s <- tapply(d2[[c1]][idx], sx, sum); ev_s[is.na(ev_s)] <- 0
    solo <- NULL
    if (length(ev_s) >= 2 && ev > 0 && min(ev_s) <= max(5, 0.01 * ev)) solo <- names(ev_s)[which.max(ev_s)]
    if (!is.null(solo)) idx <- idx[sx == solo]
    sub <- cbind(d2[idx, c(c1, "Sexo", pcols), with = FALSE], B[idx]); ids <- d2$id[idx]
    pcols <- pcols[vapply(pcols, function(p) uniqueN(sub[[p]]) > 1, logical(1))]
    f <- stats::reformulate(c(bn, if (is.null(solo)) "Sexo", pcols), c1)
    fit <- tryCatch(suppressWarnings(stats::glm(f, family = stats::binomial(), data = sub,
                                                model = FALSE, y = FALSE)), error = function(e) NULL)
    if (is.null(fit)) { skipped[[nd]] <- data.table(nodo = nd, motivo = "error_ajuste", n_events = ev); next }
    V <- tryCatch(.lox_cluster_vcov(fit, f, sub, ids), error = function(e) NULL)    # sandwich agrupado por persona (varias filas/persona)
    r <- .lox_glm_coefs(fit, "incidencia", nd, nrow(sub), ev, V)
    r <- r[!is.na(beta) & termino != "(Intercept)"]
    r[, restriccion := if (is.null(solo)) "" else paste0("solo_", solo)]
    ac <- tryCatch(.lox_age_contrast(stats::coef(fit), if (is.null(V)) summary(fit)$cov.unscaled else V, cfg, "incidencia", nd, nrow(sub), ev), error = function(e) NULL)
    res[[nd]] <- rbind(r, ac, fill = TRUE); res[[nd]][, restriccion := if (is.null(solo)) "" else paste0("solo_", solo)]
  }
  sk <- rbindlist(skipped, fill = TRUE)
  if (!nrow(sk)) sk <- data.table(nodo = character(), motivo = character(), n_events = integer())
  list(coefs = rbindlist(res), skipped = sk)
}

#' Capa 2: hazard discreto cloglog de mortalidad (exp(beta) = HR)
#'
#' Edad con spline de nudos fijos comunes, sexo y estados t0 con >= `model$min_events_hazard`
#' muertes entre portadores (el resto se excluye para evitar separacion).
#' @export
lox_fit_layer2 <- function(d2, nodes, cfg) {
  B <- lox_age_basis(d2$edad_t0, cfg); bn <- names(B)
  v0 <- paste0(nodes, "_t0")
  ev <- vapply(v0, function(v) sum(d2$death_next[d2[[v]] == 1L]), numeric(1))
  keep <- v0[ev >= cfg$model$min_events_hazard]
  dat <- cbind(d2[, c("death_next", "Sexo", keep), with = FALSE], B)
  f <- stats::reformulate(c(bn, "Sexo", keep), "death_next")
  fit <- stats::glm(f, family = stats::binomial("cloglog"), data = dat, model = FALSE, y = FALSE,
                    control = stats::glm.control(maxit = 100))
  r <- .lox_glm_coefs(fit, "mortalidad", "muerte", nrow(dat), sum(dat$death_next))
  r <- r[!is.na(beta) & termino != "(Intercept)"]
  r[, nodo := ifelse(termino %in% keep, sub("_t0$", "", termino), "muerte")]
  list(coefs = r, excluded = sub("_t0$", "", setdiff(v0, keep)), converged = isTRUE(fit$converged))
}

# Persona-anio y muertes observadas por sexo x grupo de edad x anio (agregado, para calibracion central)
.lox_observed <- function(d2, cfg) {
  g <- cut(d2$edad_t0, c(0, 50, 60, 70, 80, 90, Inf), right = FALSE,
           labels = c("<50", "50-59", "60-69", "70-79", "80-89", "90+"))
  d2[, .(n_risk = .N, n_death = sum(death_next)),
     by = .(Sexo, grupo_edad = g, anio_t0 = anio)][order(Sexo, grupo_edad, anio_t0)]
}
