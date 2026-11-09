# Nodo central: estructura consenso y meta-analisis --------------------------------------------------

#' Construye nodos y red causal consenso a partir de los resultados de la Ronda 1
#'
#' @param prev   lista nombrada (por sitio) de `r1_prevalence`
#' @param edges  lista nombrada de `r1_edges`
#' @param nodes  lista nombrada de `r1_nodes` (nodos usados en el aprendizaje local)
#' @param dirs   lista nombrada de `r1_direction`
#' @param w      pesos por sitio (transiciones usadas en el aprendizaje), nombrados
#' @return list(nodes, structure, pairs)
#' @export
lox_consensus <- function(prev, edges, nodes, dirs, w, cfg) {
  sites <- names(prev)
  # 1) nodos: prevalencia basal agrupada >= suelo (celdas suprimidas = 0)
  pv <- rbindlist(lapply(sites, function(s) prev[[s]][, .(node, n_pos = fifelse(is.na(n_pos), 0L, n_pos), n)][, site := s]))
  nd <- pv[, .(n_pos = sum(n_pos), n = sum(n)), by = node][, prev := n_pos / n]
  cons_nodes <- nd[prev >= cfg$model$prev_min][order(-prev)]
  keep <- cons_nodes$node
  # 2) fuerza de arco agrupada, solo entre sitios donde AMBOS nodos fueron nodos locales
  cand <- CJ(from = keep, to = keep)[from != to]
  rows <- lapply(sites, function(s) {
    ok <- nodes[[s]]$node
    x <- merge(cand[from %in% ok & to %in% ok], edges[[s]][, .(from, to, strength)], by = c("from", "to"), all.x = TRUE)
    x[is.na(strength), strength := 0]; x[, site := s][, w := w[[s]]]
  })
  er <- rbindlist(rows)
  thr <- cfg$consensus$threshold
  ag <- er[, .(strength = sum(w * strength) / sum(w), frac_sites = mean(strength >= thr), n_sites = .N), by = .(from, to)]
  arcs <- ag[strength >= thr & frac_sites >= cfg$consensus$min_site_frac]
  if (!nrow(arcs)) return(list(nodes = cons_nodes, structure = data.table(), pairs = ag))
  # 3) una direccion por par (la de mayor fuerza) y validacion por precedencia de fechas agrupada
  arcs[, `:=`(pa = pmin(from, to), pb = pmax(from, to))]
  arcs <- arcs[order(-strength), .SD[1], by = .(pa, pb)]
  dsum <- rbindlist(dirs, idcol = "site")[, lapply(.SD, sum, na.rm = TRUE), by = .(a, b),
                                           .SDcols = c("n_ab", "n_ba", "n_tie", "n_both")]
  arcs <- merge(arcs, dsum, by.x = c("pa", "pb"), by.y = c("a", "b"), all.x = TRUE)
  arcs[, tot := n_ab + n_ba]
  arcs[, n_apoya := fifelse(from == pa, n_ab, n_ba)]
  arcs[, frac := fifelse(!is.na(tot) & tot > 0, n_apoya / tot, NA_real_)]
  dm <- cfg$structure$dir_min
  arcs[, state := fifelse(is.na(frac), "sin datos", fifelse(frac >= dm, "concordante",
                  fifelse(frac <= 1 - dm, "reorientada", "incierta")))]
  arcs[, `:=`(dbn_from = from, dbn_to = to)]
  arcs[state == "reorientada", `:=`(from = dbn_to, to = dbn_from)]
  arcs[, pct_dir_final := fifelse(is.na(frac), NA_real_, round(100 * pmax(frac, 1 - frac)))]
  st <- arcs[, .(from, to, strength = round(strength, 3), state, pct_dir_final, n_orden = tot,
                 n_sites, frac_sites = round(frac_sites, 2), dbn_from, dbn_to)][order(-strength)]
  list(nodes = cons_nodes, structure = st, pairs = ag)
}

# Meta-analisis de efectos aleatorios (DerSimonian-Laird) por (capa, nodo, termino) ------------------
# --- Meta-analisis de efectos aleatorios ---------------------------------------------------------------
# Principal: REML para tau^2 + IC de Hartung-Knapp (t con k-1 gl, truncado para no ser mas estrecho que el IC
# normal) + intervalo de prediccion (t con k-2 gl). Sensibilidad: DerSimonian-Laird y rango "dejar un sitio fuera".
# Con k ~ 7 sitios tau^2 es muy incierto y DL tiende a subestimarlo: por eso REML+HK es el principal.

.lox_dl <- function(b, se) {
  k <- length(b); w <- 1 / se^2; fe <- sum(w * b) / sum(w)
  Q <- sum(w * (b - fe)^2); df <- k - 1; C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (k > 1 && C > 0) max(0, (Q - df) / C) else 0
  ws <- 1 / (se^2 + tau2); re <- sum(ws * b) / sum(ws); sre <- sqrt(1 / sum(ws))
  list(k = k, beta_fe = fe, se_fe = sqrt(1 / sum(w)), beta_re = re, se_re = sre,
       z = re / sre, p = 2 * stats::pnorm(-abs(re / sre)), tau2 = tau2,
       I2 = if (k > 1 && Q > 0) max(0, (Q - df) / Q) else 0, Q = Q)
}

# tau^2 por REML (iteracion de la ecuacion de estimacion; parte de DL)
.lox_reml_tau2 <- function(y, v, tol = 1e-10, maxit = 2000L) {
  k <- length(y); if (k < 2) return(0)
  tau2 <- .lox_dl(y, sqrt(v))$tau2
  for (i in seq_len(maxit)) {
    w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w)
    new <- max(0, sum(w^2 * ((y - mu)^2 - v)) / sum(w^2) + 1 / sum(w))
    if (abs(new - tau2) < tol) { tau2 <- new; break }
    tau2 <- new
  }
  tau2
}

.lox_reml_mu <- function(y, v) { t2 <- .lox_reml_tau2(y, v); w <- 1 / (v + t2); sum(w * y) / sum(w) }

.lox_meta1 <- function(b, se, min_k = 3L) {
  k <- length(b); v <- se^2; na <- NA_real_
  out <- list(k = k, status = if (k >= min_k) "ok" else sprintf("insuficiente (k<%d)", min_k),
              beta_re = na, se_re = na, se_hk = na, lo_b = na, hi_b = na, p = na, tau2 = na, I2 = na, Q = na, p_Q = na,
              pi_lo_b = na, pi_hi_b = na, beta_dl = na, se_dl = na, tau2_dl = na, loo_min_b = na, loo_max_b = na,
              beta_fe = if (k) sum(b / v) / sum(1 / v) else na, se_fe = if (k) sqrt(1 / sum(1 / v)) else na)
  if (k < min_k) return(out)
  t2 <- .lox_reml_tau2(b, v); w <- 1 / (v + t2); mu <- sum(w * b) / sum(w); se_re <- sqrt(1 / sum(w))
  wf <- 1 / v; fe <- sum(wf * b) / sum(wf); Q <- sum(wf * (b - fe)^2)
  sv <- (k - 1) * sum(wf) / (sum(wf)^2 - sum(wf^2))                  # varianza "tipica" (Higgins)
  q_hk <- sum(w * (b - mu)^2) / (k - 1); se_hk <- se_re * sqrt(max(1, q_hk))
  tc <- stats::qt(0.975, k - 1)
  out[c("beta_re", "se_re", "se_hk", "tau2", "Q")] <- list(mu, se_re, se_hk, t2, Q)
  out$lo_b <- mu - tc * se_hk; out$hi_b <- mu + tc * se_hk
  out$p <- 2 * stats::pt(-abs(mu / se_hk), k - 1); out$p_Q <- stats::pchisq(Q, k - 1, lower.tail = FALSE)
  out$I2 <- if (t2 + sv > 0) t2 / (t2 + sv) else 0
  if (k >= 3) { pw <- stats::qt(0.975, k - 2) * sqrt(t2 + se_hk^2); out$pi_lo_b <- mu - pw; out$pi_hi_b <- mu + pw }
  d <- .lox_dl(b, se); out$beta_dl <- d$beta_re; out$se_dl <- d$se_re; out$tau2_dl <- d$tau2
  if (k - 1 >= min_k) {
    lo <- vapply(seq_len(k), function(i) .lox_reml_mu(b[-i], v[-i]), 0)
    out$loo_min_b <- min(lo); out$loo_max_b <- max(lo)
  }
  out
}

#' Meta-analisis de coeficientes entre sitios (REML + Hartung-Knapp)
#'
#' @param coefs data.table con `site, capa, nodo, termino, beta, se`
#' @param min_k numero minimo de sitios para dar resultado combinado (por debajo: `status = "insuficiente"`)
#' @return una fila por (capa, nodo, termino): `k`, `beta_re` (REML), IC de Hartung-Knapp (`lo`, `hi`, en escala
#'   exp), intervalo de prediccion (`pi_lo`, `pi_hi`), `tau2`, `I2`, `Q`/`p_Q`, sensibilidad DerSimonian-Laird
#'   (`beta_dl`, `tau2_dl`) y rango dejando un sitio fuera (`loo_min`, `loo_max`, en escala exp)
#' @export
lox_meta <- function(coefs, min_k = 3L) {
  x <- coefs[is.finite(beta) & is.finite(se) & se > 0]
  m <- x[, .lox_meta1(beta, se, min_k), by = .(capa, nodo, termino)]
  m[, `:=`(exp_re = exp(beta_re), lo = exp(lo_b), hi = exp(hi_b), pi_lo = exp(pi_lo_b), pi_hi = exp(pi_hi_b),
           loo_min = exp(loo_min_b), loo_max = exp(loo_max_b))]
  m[, c("lo_b", "hi_b", "pi_lo_b", "pi_hi_b", "loo_min_b", "loo_max_b") := NULL]
  m[order(capa, nodo, termino)]
}
