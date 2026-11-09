# Utilidades de la ronda 2 ampliada (0.7.5): errores estandar robustos, contraste de edad, tablas por arco,
# incidencia por estratos y metricas de red del nodo central. Todo agregado; nada individual sale del sitio.

# Matriz de covarianza robusta por conglomerado (sandwich agrupado por persona) para un GLM binomial/Poisson.
# Las filas de `data` deben ser las mismas (y en el mismo orden) que las del ajuste.
.lox_cluster_vcov <- function(fit, formula, data, id) {
  if (!is.null(fit$na.action)) return(NULL)
  b <- stats::coef(fit); ok <- !is.na(b)
  X <- stats::model.matrix(formula, data)[, ok, drop = FALSE]
  y <- data[[all.vars(formula)[1]]]
  if (nrow(X) != length(y) || nrow(X) != length(id)) return(NULL)
  U <- X * (y - fit$fitted.values)                       # scores: x_i (y_i - mu_i), valido para logit y para Poisson con log
  gi <- match(id, unique(id)); Ug <- rowsum(U, gi, reorder = FALSE); G <- nrow(Ug)
  if (G < 3) return(NULL)
  bread <- summary(fit)$cov.unscaled[ok, ok, drop = FALSE]
  V <- G / (G - 1) * bread %*% crossprod(Ug) %*% bread
  dimnames(V) <- dimnames(bread); V
}

# Contraste de edad: log-OR / log-HR de pasar de `from` a `to` anios (efecto de la spline resumido en UN numero meta-analizable)
.lox_age_contrast <- function(beta, V, cfg, capa, nodo, n, n_events, from = 70, to = 80) {
  B <- as.matrix(lox_age_basis(c(from, to), cfg)); L <- B[2, ] - B[1, ]
  nm <- names(L); nm <- nm[nm %in% names(beta) & !is.na(beta[nm]) & nm %in% rownames(V)]
  if (length(nm) != length(L)) return(NULL)
  l <- L[nm]; est <- sum(l * beta[nm]); se <- sqrt(drop(t(l) %*% V[nm, nm, drop = FALSE] %*% l))
  data.table(capa = capa, nodo = nodo, termino = sprintf("edad_%da%d", from, to), beta = est, se = se,
             n = n, n_events = n_events, converged = TRUE, flag = "")
}

#' Tablas 2x2 por arco del consenso: enfermedad A en t0 (si/no) frente a B en t1 entre quienes NO tenian B en t0
#' @param d2 transiciones (`lox_build_transitions`)
#' @param edges data.table(from, to) del consenso
#' @return data.table con `from`, `to`, `a1_b1`, `a1_b0`, `a0_b1`, `a0_b0` (celdas < `min_cell` suprimidas)
#' @export
lox_arc_tables <- function(d2, edges, min_cell = 10L) {
  cols <- c("a1_b1", "a1_b0", "a0_b1", "a0_b0")
  if (is.null(edges) || !nrow(edges)) return(data.table(from = character(), to = character(), a1_b1 = integer(), a1_b0 = integer(), a0_b1 = integer(), a0_b0 = integer()))
  r <- rbindlist(lapply(seq_len(nrow(edges)), function(i) {
    a <- paste0(edges$from[i], "_t0"); b0 <- paste0(edges$to[i], "_t0"); b1 <- paste0(edges$to[i], "_t1")
    if (!all(c(a, b0, b1) %in% names(d2))) return(NULL)
    k <- d2[[b0]] == 0L; A <- d2[[a]][k]; B <- d2[[b1]][k]
    data.table(from = edges$from[i], to = edges$to[i], a1_b1 = sum(A == 1L & B == 1L), a1_b0 = sum(A == 1L & B == 0L),
               a0_b1 = sum(A == 0L & B == 1L), a0_b0 = sum(A == 0L & B == 0L))
  }))
  for (cl in cols) r[get(cl) > 0L & get(cl) < min_cell, (cl) := NA_integer_]
  r[]
}

#' Personas en riesgo e incidencia por enfermedad, anio, sexo y grupo de edad (al 1 de enero)
#' @return data.table `node`, `anio`, `Sexo`, `edad_grupo`, `n_risk`, `n_inc` (n_inc < `min_cell` -> NA)
#' @export
lox_dx_strata <- function(cohort, nodes, cfg, min_cell = 10L) {
  yf <- cfg$cohort$year_first; yl <- cfg$cohort$year_last
  P <- cohort$cohorte[, .(id, Sexo, edad_base, dth = as.numeric(death_date), bj = as.numeric(baja_date))]
  b1 <- as.numeric(as.Date(sprintf("%d-01-01", yl + 1)))
  P[, end := pmin(fifelse(is.na(dth), Inf, dth), fifelse(is.na(bj), Inf, bj), b1)]
  e <- cohort$enf[code %in% nodes, .(id, code, d = as.numeric(dx_date))]
  e[P, `:=`(end = i.end, Sexo = i.Sexo, eb = i.edad_base), on = "id"]; e <- e[!is.na(end)]
  grp <- function(a) cut(a, c(-Inf, 70, 80, 90, Inf), labels = c("<70", "70-79", "80-89", "90+"), right = FALSE)
  out <- vector("list", yl - yf + 1)
  for (y in yf:yl) {
    j1 <- as.numeric(as.Date(sprintf("%d-01-01", y))); j2 <- as.numeric(as.Date(sprintf("%d-01-01", y + 1)))
    al <- P[end > j1]; al[, edad_grupo := grp(edad_base + (y - yf))]
    nal <- al[, .(alive = .N), by = .(Sexo, edad_grupo)]
    e[, edad_grupo := grp(eb + (y - yf))]
    prev <- (if (y == yf) e[d <= j1 & end > j1] else e[d < j1 & end > j1])[, .(prev = .N), by = .(code, Sexo, edad_grupo)]
    inc  <- (if (y == yf) e[d > j1 & d < j2 & d <= end] else e[d >= j1 & d < j2 & d <= end])[, .(n_inc = .N), by = .(code, Sexo, edad_grupo)]
    o <- CJ(node = nodes, Sexo = unique(P$Sexo), edad_grupo = factor(levels(grp(0)), levels = levels(grp(0))))
    o[nal, alive := i.alive, on = .(Sexo, edad_grupo)]; o[is.na(alive), alive := 0L]
    o[prev, prev := i.prev, on = .(node = code, Sexo, edad_grupo)]; o[is.na(prev), prev := 0L]
    o[inc, n_inc := i.n_inc, on = .(node = code, Sexo, edad_grupo)]; o[is.na(n_inc), n_inc := 0L]
    o[, `:=`(n_risk = alive - prev, anio = y)][, c("alive", "prev") := NULL]
    out[[y - yf + 1]] <- o[n_risk > 0]
  }
  r <- rbindlist(out); r[, edad_grupo := as.character(edad_grupo)]
  r[n_inc > 0L & n_inc < min_cell, n_inc := NA_integer_]
  r[, n_risk := ifelse(n_risk < min_cell, NA_integer_, n_risk)]
  r[order(node, anio, Sexo, edad_grupo)][]
}

#' Metricas de centralidad de la red (pura R, sin dependencias): grado, fuerza, centralidad armonica e intermediacion
#'
#' `out_*` miden la capacidad de un nodo de "disparar" a otros (disparadores); `in_*`, de recibir (receptores).
#' La centralidad armonica usa distancias dirigidas; la intermediacion, caminos mas cortos no ponderados (Brandes).
#' @param edges data.frame(from, to, [strength])
#' @export
lox_network_metrics <- function(edges) {
  edges <- as.data.frame(edges); nd <- sort(unique(c(edges$from, edges$to))); n <- length(nd)
  if (!n) return(data.table(node = character()))
  w <- if ("strength" %in% names(edges)) edges$strength else rep(1, nrow(edges))
  ix <- function(v) match(v, nd); from <- ix(edges$from); to <- ix(edges$to)
  adj <- lapply(seq_len(n), function(i) to[from == i])
  dist_from <- function(s) { d <- rep(Inf, n); d[s] <- 0; q <- s
    while (length(q)) { u <- q[1]; q <- q[-1]; for (v in adj[[u]]) if (is.infinite(d[v])) { d[v] <- d[u] + 1; q <- c(q, v) } }; d }
  D <- t(vapply(seq_len(n), dist_from, numeric(n)))            # D[i, j] = distancia i -> j
  H <- function(m) { m[is.infinite(m) | m == 0] <- NA; rowSums(1 / m, na.rm = TRUE) / max(n - 1, 1) }
  bet <- numeric(n)                                            # Brandes, grafo dirigido no ponderado
  for (s in seq_len(n)) {
    S <- integer(); P <- vector("list", n); sg <- numeric(n); sg[s] <- 1; d <- rep(-1, n); d[s] <- 0; q <- s
    while (length(q)) { v <- q[1]; q <- q[-1]; S <- c(S, v)
      for (wv in adj[[v]]) { if (d[wv] < 0) { q <- c(q, wv); d[wv] <- d[v] + 1 }
        if (d[wv] == d[v] + 1) { sg[wv] <- sg[wv] + sg[v]; P[[wv]] <- c(P[[wv]], v) } } }
    dl <- numeric(n)
    for (wv in rev(S)) { for (v in P[[wv]]) dl[v] <- dl[v] + sg[v] / sg[wv] * (1 + dl[wv]); if (wv != s) bet[wv] <- bet[wv] + dl[wv] }
  }
  data.table(node = nd, out_degree = tabulate(from, n), in_degree = tabulate(to, n),
             out_strength = vapply(seq_len(n), function(i) sum(w[from == i]), 0), in_strength = vapply(seq_len(n), function(i) sum(w[to == i]), 0),
             harmonic_out = H(D), harmonic_in = H(t(D)), betweenness = bet)
}

#' Comunidades (Louvain o Leiden) sobre la red no dirigida y ponderada por fuerza; requiere `igraph`
#' @param resolution parametro de resolucion (mas alto = mas comunidades)
#' @export
lox_communities <- function(edges, method = c("louvain", "leiden"), resolution = 1, seed = 1L) {
  method <- match.arg(method)
  if (!requireNamespace("igraph", quietly = TRUE)) return(NULL)
  edges <- as.data.frame(edges); if (!nrow(edges)) return(NULL)
  g <- igraph::graph_from_data_frame(data.frame(from = edges$from, to = edges$to, weight = if ("strength" %in% names(edges)) edges$strength else 1), directed = FALSE)
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum"))
  set.seed(seed)
  cl <- if (method == "louvain") igraph::cluster_louvain(g, resolution = resolution)
        else igraph::cluster_leiden(g, objective_function = "modularity", resolution = resolution)
  data.table(node = igraph::V(g)$name, community = as.integer(igraph::membership(cl)), modularity = igraph::modularity(cl))
}

# Indice de Rand ajustado entre dos particiones de los mismos nodos
.lox_ari <- function(a, b) {
  tb <- table(a, b); n <- sum(tb); c2 <- function(x) x * (x - 1) / 2
  si <- sum(c2(tb)); sa <- sum(c2(rowSums(tb))); sb <- sum(c2(colSums(tb))); ex <- sa * sb / c2(n)
  if (isTRUE(all.equal(0.5 * (sa + sb), ex))) return(1)
  (si - ex) / (0.5 * (sa + sb) - ex)
}
