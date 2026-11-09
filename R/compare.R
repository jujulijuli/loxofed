# Comparador de dos ejecuciones (anterior / actual) ---------------------------------------------------
# Compara las salidas AGREGADAS de dos ejecuciones de un mismo sitio (o de dos analisis centrales): cohorte,
# estructura consenso, HR de mortalidad, OR de incidencia y calibracion O/E. Sirve para documentar el efecto de
# un filtro de datos, de un tope de muestreo, de otra version del paquete o de un cambio de configuracion.
# HTML autocontenido en R puro; solo lee tablas que ya se comparten entre sitios.

.lox_find_consensus <- function(d) {
  for (s in c("consensus", file.path("central_tmp", "consensus"), file.path("central", "consensus")))
    if (file.exists(file.path(d, s, "consensus_structure.csv"))) return(file.path(d, s))
  NULL
}

# SVG: dispersion A (x) frente a B (y) con diagonal y etiquetas de los puntos destacados
.lox_svg_scatter <- function(x, y, lab, hl = rep(FALSE, length(x)), xlab = "", ylab = "", log = TRUE, w = 460) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; lab <- lab[ok]; hl <- hl[ok]
  if (length(x) < 2) return("<p class='muted'>Sin datos suficientes.</p>")
  ml <- 52; mb <- 40; mt <- 10; mr <- 12; h <- w
  tr <- if (log) base::log else identity; r <- range(tr(c(x, y))); r <- r + c(-1, 1) * 0.06 * diff(r)
  X <- function(v) ml + (w - ml - mr) * (tr(v) - r[1]) / diff(r); Y <- function(v) mt + (h - mt - mb) * (1 - (tr(v) - r[1]) / diff(r))
  tk <- if (log) c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 8, 12) else pretty(r, 5)
  tk <- tk[tr(tk) >= r[1] & tr(tk) <= r[2]]
  grid <- paste0(paste0(sprintf("<line x1='%.1f' x2='%.1f' y1='%d' y2='%d' class='grid'/><text x='%.1f' y='%d' text-anchor='middle' class='lbl'>%s</text>",
                                X(tk), X(tk), mt, h - mb, X(tk), h - mb + 14, format(tk)), collapse = ""),
                 paste0(sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' class='grid'/><text x='%d' y='%.1f' text-anchor='end' class='lbl'>%s</text>",
                                ml, w - mr, Y(tk), Y(tk), ml - 6, Y(tk) + 4, format(tk)), collapse = ""))
  diag <- sprintf("<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' class='ref' stroke-dasharray='4 3'/>", X(exp(r[1])), Y(exp(r[1])), X(exp(r[2])), Y(exp(r[2])))
  if (!log) diag <- sprintf("<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' class='ref' stroke-dasharray='4 3'/>", X(r[1]), Y(r[1]), X(r[2]), Y(r[2]))
  pts <- paste0(sprintf("<circle cx='%.1f' cy='%.1f' r='%s' class='%s'/>", X(x), Y(y), ifelse(hl, "4.5", "3"), ifelse(hl, "pt-up", "pt-ns")), collapse = "")
  tx <- if (any(hl)) paste0(sprintf("<text x='%.1f' y='%.1f' text-anchor='%s' class='lbl halo'>%s</text>", X(x[hl]) + ifelse(X(x[hl]) > 0.62 * w, -6, 6), Y(y[hl]) - 5,
                                               ifelse(X(x[hl]) > 0.62 * w, "end", "start"), .lox_esc(lab[hl])), collapse = "") else ""
  sprintf("<svg viewBox='0 0 %d %d' width='100%%' style='max-width:%dpx' role='img'>%s%s%s%s<text x='%.1f' y='%d' text-anchor='middle' class='lbl'>%s</text><text x='12' y='%.1f' transform='rotate(-90 12 %.1f)' text-anchor='middle' class='lbl'>%s</text></svg>",
          w, h, w, grid, diag, pts, tx, (ml + w - mr) / 2, h - 4, .lox_esc(xlab), (mt + h - mb) / 2, (mt + h - mb) / 2, .lox_esc(ylab))
}

# cambio relativo en % con signo (sin "-0")
.lox_chg <- function(r) { v <- as.integer(round(100 * (r - 1))); ifelse(v == 0L, "0", sprintf("%+d", v)) }

#' Informe comparativo de dos ejecuciones
#'
#' @param dir_a,dir_b carpetas de salida (cada una con `round1/` y/o `round2/`); A = anterior, B = actual
#' @param labels etiquetas de las dos ejecuciones
#' @param consensus_a,consensus_b carpetas `consensus/` (si no se dan se buscan en `consensus/`,
#'   `central_tmp/consensus/` o `central/consensus/`)
#' @param out_file fichero HTML de salida
#' @export
lox_compare_runs <- function(dir_a, dir_b, labels = c("anterior", "actual"), consensus_a = NULL, consensus_b = NULL,
                             out_file = file.path(dir_b, "comparacion.html")) {
  rd <- function(d, sub, f) { p <- file.path(d, sub, f); if (file.exists(p)) data.table::fread(p, na.strings = c("NA", "")) else NULL }
  mj <- function(d, sub) { p <- file.path(d, sub, "manifest.json"); if (file.exists(p)) jsonlite::read_json(p) else NULL }
  La <- labels[1]; Lb <- labels[2]
  if (is.null(consensus_a)) consensus_a <- .lox_find_consensus(dir_a)
  if (is.null(consensus_b)) consensus_b <- .lox_find_consensus(dir_b)
  A <- list(meta = rd(dir_a, "round1", "r1_meta.csv"), co = rd(dir_a, "round2", "r2_coefs.csv"), obs = rd(dir_a, "round2", "r2_observed.csv"),
            m2 = mj(dir_a, "round2"), m1 = mj(dir_a, "round1"),
            cs = if (!is.null(consensus_a)) data.table::fread(file.path(consensus_a, "consensus_structure.csv")) else NULL)
  B <- list(meta = rd(dir_b, "round1", "r1_meta.csv"), co = rd(dir_b, "round2", "r2_coefs.csv"), obs = rd(dir_b, "round2", "r2_observed.csv"),
            m2 = mj(dir_b, "round2"), m1 = mj(dir_b, "round1"),
            cs = if (!is.null(consensus_b)) data.table::fread(file.path(consensus_b, "consensus_structure.csv")) else NULL)
  if (all(vapply(c(A, B), is.null, TRUE))) stop("No se encontraron salidas en ", dir_a, " ni ", dir_b, call. = FALSE)
  H <- character(); add <- function(...) H[length(H) + 1] <<- paste0(...)
  sm <- character()                                                    # frases del resumen
  pct <- function(x) sprintf("%.0f%%", 100 * x)
  vers <- function(X) { m <- if (!is.null(X$m2)) X$m2 else X$m1; if (is.null(m)) NA else m$package_version }

  add("<h1>Comparación de ejecuciones · ", .lox_esc(La), " frente a ", .lox_esc(Lb), "</h1><p class='muted'>loxofed ",
      .lox_esc(lox_version()), " · generado ", .lox_esc(format(Sys.time(), "%Y-%m-%d %H:%M")), ". Todas las cifras proceden de tablas agregadas ya compartibles.</p>")

  # ---- 1. cohorte
  sec_coh <- ""
  if (!is.null(A$meta) && !is.null(B$meta)) {
    g <- function(X, v) if (v %in% names(X$meta)) X$meta[[v]][1] else NA
    py <- function(X) if (!is.null(X$obs) && "person_years" %in% names(X$obs)) sum(X$obs$person_years, na.rm = TRUE) else NA
    dd <- function(X) if (!is.null(X$co)) suppressWarnings(max(X$co[capa == "mortalidad", n_events], na.rm = TRUE)) else NA
    tb <- data.frame(Medida = c("Versión de loxofed", "Personas elegibles", "Personas en el análisis", "Peso de muestreo", "Transiciones (capa 1)",
                                "Muertes en el ajuste", "Personas-año de seguimiento (observadas)"),
                     a = c(vers(A), .lox_big(g(A, "n_eligible")), .lox_big(g(A, "n_cohort")), .lox_fmt(g(A, "sample_weight"), 2), .lox_big(g(A, "n_transitions")), .lox_big(dd(A)), .lox_big(py(A))),
                     b = c(vers(B), .lox_big(g(B, "n_eligible")), .lox_big(g(B, "n_cohort")), .lox_fmt(g(B, "sample_weight"), 2), .lox_big(g(B, "n_transitions")), .lox_big(dd(B)), .lox_big(py(B))))
    names(tb)[2:3] <- c(La, Lb); sec_coh <- .lox_html_table(tb, num = c(La, Lb))
    dn <- g(B, "n_cohort") - g(A, "n_cohort")
    sm <- c(sm, sprintf("<b>Cohorte:</b> %s → %s personas en el análisis (%s%s).", .lox_big(g(A, "n_cohort")), .lox_big(g(B, "n_cohort")),
                        ifelse(dn >= 0, "+", "−"), .lox_big(abs(dn))))
    if (!is.na(vers(A)) && !is.na(vers(B)) && vers(A) != vers(B)) sm <- c(sm, sprintf("<b>Versión:</b> %s frente a %s; el centro no combina versiones distintas.", vers(A), vers(B)))
  }

  # ---- 2. estructura
  sec_st <- ""
  if (!is.null(A$cs) && !is.null(B$cs)) {
    act <- function(cs) cs[state %in% c("concordante", "reorientada")]
    ka <- act(A$cs)[, paste(dbn_from, dbn_to, sep = " → ")]; kb <- act(B$cs)[, paste(dbn_from, dbn_to, sep = " → ")]
    both <- intersect(ka, kb); jac <- length(both) / length(union(ka, kb))
    sm <- c(sm, sprintf("<b>Estructura:</b> %d arcos causales en %s y %d en %s; comunes %d (acuerdo de Jaccard %s). Entran %d y salen %d.",
                        length(ka), La, length(kb), Lb, length(both), pct(jac), length(setdiff(kb, ka)), length(setdiff(ka, kb))))
    M <- merge(A$cs[, .(arco = paste(dbn_from, dbn_to, sep = " → "), fa = strength, ea = state)],
               B$cs[, .(arco = paste(dbn_from, dbn_to, sep = " → "), fb = strength, eb = state)], by = "arco", all = TRUE)
    M[, `:=`(ea = ifelse(is.na(ea), "—", ea), eb = ifelse(is.na(eb), "—", eb))]
    M[, `:=`(activo_a = ea %in% c("concordante", "reorientada"), activo_b = eb %in% c("concordante", "reorientada"))]
    ch <- M[activo_a != activo_b | (activo_a & activo_b & abs(fa - fb) >= 0.1) | (ea != eb & (activo_a | activo_b))]
    if (nrow(ch)) {
      ch[, `Cambio` := ifelse(!activo_a & activo_b, "entra", ifelse(activo_a & !activo_b, "sale", ifelse(ea != eb, "cambia estado", "cambia fuerza")))]
      ch <- ch[order(factor(Cambio, c("entra", "sale", "cambia estado", "cambia fuerza")), -pmax(fa, fb, na.rm = TRUE))]
      tbc <- ch[, .(Arco = .lox_short(sub(" → .*", "", arco)), Hasta = .lox_short(sub(".* → ", "", arco)), Cambio,
                    fuerza_a = .lox_fmt(fa, 2), estado_a = ea, fuerza_b = .lox_fmt(fb, 2), estado_b = eb)]
      setnames(tbc, c("Desde", "Hasta", "Cambio", paste("Fuerza", La), paste("Estado", La), paste("Fuerza", Lb), paste("Estado", Lb)))
      sec_st <- paste0("<p>Arcos que entran o salen del consenso, cambian de estado o cambian de fuerza en 0,10 o más (de ", nrow(M), " pares candidatos).</p>",
                       .lox_html_table(tbc, num = c(paste("Fuerza", La), paste("Fuerza", Lb))))
    } else sec_st <- "<p>La estructura consenso es idéntica en ambas ejecuciones (mismos arcos activos, fuerzas dentro de 0,10).</p>"
    pm <- M[!is.na(fa) & !is.na(fb)]
    if (nrow(pm) > 2) sec_st <- paste0(.lox_svg_scatter(pm$fa, pm$fb, .lox_short(pm$arco), hl = abs(pm$fa - pm$fb) >= 0.1, xlab = paste("Fuerza", La), ylab = paste("Fuerza", Lb), log = FALSE),
                                       "<p class='muted'>Cada punto es un par candidato. Los naranjas difieren en 0,10 o más.</p>", sec_st)
  } else sec_st <- "<p class='muted'>Falta el consenso en una de las dos ejecuciones (indica <code>consensus_a</code> / <code>consensus_b</code>).</p>"

  # ---- 3. HR de mortalidad y 4. OR de incidencia
  sec_hr <- ""; sec_or <- ""
  if (!is.null(A$co) && !is.null(B$co)) {
    j <- merge(A$co[, .(capa, nodo, termino, ba = beta, sa = se)], B$co[, .(capa, nodo, termino, bb = beta, sb = se)], by = c("capa", "nodo", "termino"))
    j[, `:=`(d = bb - ba, z = (bb - ba) / sqrt(sa^2 + sb^2))]
    hz <- j[capa == "mortalidad" & nodo != "muerte"]
    if (nrow(hz) > 2) {
      hz[, `:=`(ra = rank(-ba), rb = rank(-bb))]
      rho <- suppressWarnings(stats::cor(hz$ba, hz$bb, method = "spearman")); r <- stats::cor(hz$ba, hz$bb)
      sm <- c(sm, sprintf("<b>HR de mortalidad:</b> %d estados comparables; correlación de log-HR %.3f (Spearman %.3f); mayor cambio %.0f%% en el HR de %s.",
                          nrow(hz), r, rho, 100 * (exp(max(abs(hz$d))) - 1), .lox_short(hz$nodo[which.max(abs(hz$d))])))
      hz[, `:=`(hra = exp(ba), hrb = exp(bb))]
      top <- hz[order(-abs(d))]
      tb <- top[, .(Estado = .lox_short(nodo), a = .lox_fmt(hra, 2), b = .lox_fmt(hrb, 2), `Cambio %` = .lox_chg(hrb / hra),
                    `Puesto a` = ra, `Puesto b` = rb)]
      setnames(tb, c("Estado", paste("HR", La), paste("HR", Lb), "Cambio %", paste("Puesto", La), paste("Puesto", Lb)))
      sec_hr <- paste0(.lox_svg_scatter(hz$hra, hz$hrb, .lox_short(hz$nodo), hl = rank(-abs(hz$d)) <= 5, xlab = paste("HR", La), ylab = paste("HR", Lb)),
                       "<p class='muted'>Escala logarítmica; los cinco estados con mayor cambio están etiquetados. La línea discontinua es la igualdad.</p>",
                       "<details open><summary>Tabla (ordenada por tamaño del cambio)</summary>", .lox_html_table(tb, num = names(tb)[-1]), "</details>")
    }
    inc <- j[capa == "incidencia" & !grepl("^edad_|^Sexo", termino)]
    if (nrow(inc) > 2) {
      r <- stats::cor(inc$ba, inc$bb)
      sm <- c(sm, sprintf("<b>OR de incidencia:</b> %d arcos comparables; correlación de log-OR %.3f; |Δ log-OR| máximo %.2f; %s de los arcos difieren menos de 2 errores estándar combinados.",
                          nrow(inc), r, max(abs(inc$d)), pct(mean(abs(inc$z) < 2))))
      inc[, `:=`(ora = exp(ba), orb = exp(bb))]
      top <- inc[order(-abs(d))][seq_len(min(15, .N))]
      tb <- top[, .(Desde = .lox_short(sub("_t0$", "", termino)), Hasta = .lox_short(nodo), a = .lox_fmt(ora, 2), b = .lox_fmt(orb, 2),
                    `Cambio %` = .lox_chg(orb / ora))]
      setnames(tb, c("Desde", "Hasta", paste("OR", La), paste("OR", Lb), "Cambio %"))
      sec_or <- paste0(.lox_svg_scatter(inc$ora, inc$orb, paste(.lox_short(sub("_t0$", "", inc$termino)), "→", .lox_short(inc$nodo)),
                                        hl = rank(-abs(inc$d)) <= 3, xlab = paste("OR", La), ylab = paste("OR", Lb)),
                       "<details><summary>Los 15 arcos con mayor cambio</summary>", .lox_html_table(tb, num = names(tb)[-(1:2)]), "</details>")
    }
  } else sec_hr <- "<p class='muted'>Falta <code>round2/r2_coefs.csv</code> en una de las dos ejecuciones (hay que haber completado la ronda 2 en ambas).</p>"

  # ---- 5. calibracion
  sec_oe <- ""
  oa <- .lox_oe(A$obs, "grupo_edad"); ob <- .lox_oe(B$obs, "grupo_edad")
  if (!is.null(oa) && !is.null(ob)) {
    lv <- c("<50", "50-59", "60-69", "70-79", "80-89", "90+"); m <- merge(oa, ob, by = "grupo_edad", suffixes = c("_a", "_b")); m <- m[order(match(grupo_edad, lv))]
    tb <- m[, .(`Grupo de edad` = grupo_edad, a = .lox_fmt(oe_a, 2), b = .lox_fmt(oe_b, 2))]; setnames(tb, c("Grupo de edad", paste("O/E", La), paste("O/E", Lb)))
    sec_oe <- paste0("<p>Cociente muertes observadas / esperadas por grupo de edad; 1 = el modelo reproduce la mortalidad del grupo.</p>", .lox_html_table(tb, num = names(tb)[-1]))
    sm <- c(sm, sprintf("<b>Calibración:</b> el O/E por grupo de edad va de %.2f a %.2f en %s y de %.2f a %.2f en %s.", min(m$oe_a), max(m$oe_a), La, min(m$oe_b), max(m$oe_b), Lb))
  } else sec_oe <- "<p class='muted'>La calibración requiere <code>n_exp</code> en <code>r2_observed.csv</code> de ambas ejecuciones (loxofed &gt;= 0.6.0 en la ronda 2).</p>"

  add("<h2>Resumen</h2>"); add(if (length(sm)) paste0("<ul>", paste0("<li>", sm, "</li>", collapse = ""), "</ul>") else "<p class='muted'>Sin datos comparables.</p>")
  add("<h2>1. Cohorte</h2>", sec_coh); add("<h2>2. Estructura consenso</h2>", sec_st)
  add("<h2>3. HR de mortalidad</h2>", sec_hr); add("<h2>4. OR de incidencia</h2>", sec_or); add("<h2>5. Calibración O/E</h2>", sec_oe)
  add("<div class='note'><b>Cómo leerlo.</b> Si las dos ejecuciones usan los mismos datos de base (p. ej. muestra frente a cohorte completa), no son independientes: las diferencias son orientativas y el valor de Δ/EE no es una prueba. ",
      "Si el filtro o el cambio de versión no altera la estructura ni los efectos, las correlaciones se acercan a 1 y los arcos que entran o salen son los de fuerza intermedia.</div>")
  html <- paste0("<!doctype html><html lang='es'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Comparación ",
                 .lox_esc(La), " vs ", .lox_esc(Lb), "</title><style>", .lox_css, "</style></head><body><main>", paste(H, collapse = "\n"), "</main></body></html>")
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE); writeLines(html, out_file, useBytes = TRUE)
  .lox_msg("Comparación -> %s", out_file); invisible(out_file)
}
