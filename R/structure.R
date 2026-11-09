# Ronda 1: aprendizaje de estructura y estadisticos agregados --------------------------------------

#' Prevalencia basal (1 de enero del primer anio) de las 60 condiciones
#' @return data.table(node, n_pos, n, prev)
#' @export
lox_prevalence <- function(arr) {
  base <- arr[, , 1, drop = TRUE]
  data.table(node = colnames(base), n_pos = as.integer(colSums(base)), n = nrow(base),
             prev = colSums(base) / nrow(base))
}

#' Precedencia empirica por fechas para todos los pares de nodos
#'
#' Para cada par no ordenado (a < b alfabeticamente) cuenta en cuantos pacientes con ambos
#' diagnosticos `a` precede a `b` (`n_ab`), `b` a `a` (`n_ba`) o coinciden (`n_tie`). Con fechas
#' anuales los empates son frecuentes; los empates NO cuentan como orden (igual que el informe).
#' @export
lox_direction_counts <- function(cohort, nodes) {
  nodes <- sort(unique(nodes)); ids <- cohort$cohorte$id
  M <- matrix(NA_real_, length(ids), length(nodes), dimnames = list(NULL, nodes))
  e <- cohort$enf[code %in% nodes]
  M[cbind(match(e$id, ids), match(e$code, nodes))] <- as.numeric(e$dx_date)
  pr <- t(utils::combn(nodes, 2))
  res <- lapply(seq_len(nrow(pr)), function(i) {
    x <- M[, pr[i, 1]]; y <- M[, pr[i, 2]]; ok <- !is.na(x) & !is.na(y)
    x <- x[ok]; y <- y[ok]
    c(n_ab = sum(x < y), n_ba = sum(y < x), n_tie = sum(x == y), n_both = length(x))
  })
  data.table(a = pr[, 1], b = pr[, 2], do.call(rbind, res))
}

#' Aprende la estructura 2-TBN (boot.strength sobre una submuestra) y devuelve las fuerzas de arco
#'
#' Replica el informe: raices = grupo de edad + sexo; lista negra que impide (i) arcos hacia las
#' raices, (ii) arcos hacia el pasado (t0) y (iii) arcos t1->t1 (co-ocurrencia, no sucesion).
#' Se exportan TODAS las fuerzas de progresion A_t0 -> B_t1 (no solo las que pasan el umbral),
#' para que el nodo central decida el consenso.
#' @return list(edges = data.table(from, to, strength, direction), n_fit)
#' @export
lox_learn_structure <- function(d2, nodes, cfg) {
  sc <- cfg$structure; set.seed(sc$seed)
  sub <- if (!is.na(sc$subsample) && nrow(d2) > sc$subsample) d2[sample(.N, sc$subsample)] else d2
  v0 <- paste0(nodes, "_t0"); v1 <- paste0(nodes, "_t1")
  d1 <- data.frame(grupo_edad = cut(sub$edad_t0, c(0, 50, 60, 70, 80, 90, Inf), right = FALSE,
                                    labels = c("<50", "50-59", "60-69", "70-79", "80-89", "90+")),
                   Sexo = sub$Sexo)
  for (v in c(v0, v1)) d1[[v]] <- factor(sub[[v]], levels = c(0, 1))
  d1 <- droplevels(d1)
  raices <- c("grupo_edad", "Sexo")
  ng <- names(d1)[vapply(d1, nlevels, integer(1)) >= 2]
  bl <- rbind(expand.grid(from = ng, to = intersect(raices, ng), stringsAsFactors = FALSE),
              expand.grid(from = ng, to = intersect(v0, ng),     stringsAsFactors = FALSE),
              expand.grid(from = intersect(v1, ng), to = intersect(v1, ng), stringsAsFactors = FALSE))
  bl <- unique(bl[bl$from != bl$to, ])
  st <- bnlearn::boot.strength(d1[, ng, drop = FALSE], R = sc$boot_R, algorithm = sc$algorithm,
                               algorithm.args = list(blacklist = bl, score = sc$score, iss = sc$iss))
  st <- as.data.table(st)
  ed <- st[grepl("_t0$", from) & grepl("_t1$", to) & sub("_t0$", "", from) != sub("_t1$", "", to)]
  ed <- ed[, .(from = sub("_t0$", "", from), to = sub("_t1$", "", to), strength, direction)]
  list(edges = ed[strength > 0], n_fit = nrow(sub),
       nodes_fit = unique(sub("_t[01]$", "", ng[grepl("_t[01]$", ng)])))
}
