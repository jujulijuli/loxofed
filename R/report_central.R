# Informe central en HTML (R puro, autocontenido): consenso + meta-analisis entre sitios ------------------

# SVG: un arco/coeficiente -> estimaciones por sitio (IC95% normal) + efecto combinado (diamante) + prediccion
.lox_svg_forest_arc <- function(site, b, se, pooled, w = 470) {
  k <- length(b); rowh <- 19; mt <- 22; lw <- 70; mr <- 118; h <- mt + rowh * (k + 1) + 12
  lo <- b - 1.96 * se; hi <- b + 1.96 * se
  allv <- c(lo, hi, 0, if (!is.na(pooled$pi_lo)) c(log(pooled$pi_lo), log(pooled$pi_hi)))
  r <- range(allv[is.finite(allv)]); r <- r + c(-1, 1) * 0.05 * max(diff(r), 0.2)
  X <- function(v) lw + (w - lw - mr) * (v - r[1]) / diff(r)
  tk <- pretty(exp(r), 4); tk <- tk[tk > 0 & log(tk) >= r[1] & log(tk) <= r[2]]
  g <- paste0(paste0(sprintf("<line x1='%.1f' x2='%.1f' y1='%d' y2='%d' class='%s'/><text x='%.1f' y='%d' text-anchor='middle' class='lbl mu'>%s</text>",
                      X(log(tk)), X(log(tk)), mt - 6, h - 8, "grid", X(log(tk)), mt - 9, format(tk)), collapse = ""),
              sprintf("<line x1='%.1f' x2='%.1f' y1='%d' y2='%d' class='ref'/>", X(0), X(0), mt - 6, h - 8))
  y <- mt + rowh * (seq_len(k) - 0.5)
  rows <- paste0(sprintf("<text x='%d' y='%.1f' text-anchor='end' class='lbl'>%s</text><line x1='%.1f' x2='%.1f' y1='%.1f' y2='%.1f' class='ci'/><rect x='%.1f' y='%.1f' width='6' height='6' class='sq'/><text x='%d' y='%.1f' class='lbl mu'>%.2f (%.2f–%.2f)</text>",
                         lw - 6, y + 4, .lox_esc(site), X(lo), X(hi), y, y, X(b) - 3, y - 3, w - mr + 6, y + 4, exp(b), exp(lo), exp(hi)), collapse = "")
  yp <- mt + rowh * (k + 0.5)
  pool <- ""
  if (!is.na(pooled$beta)) {
    if (!is.na(pooled$pi_lo)) pool <- sprintf("<line x1='%.1f' x2='%.1f' y1='%.1f' y2='%.1f' class='pi'/>", X(log(pooled$pi_lo)), X(log(pooled$pi_hi)), yp, yp)
    pool <- paste0(pool, sprintf("<polygon points='%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f' class='dia'/><text x='%d' y='%.1f' class='lbl' font-weight='600'>%.2f (%.2f–%.2f)</text>",
                    X(log(pooled$lo)), yp, X(pooled$beta), yp - 6, X(log(pooled$hi)), yp, X(pooled$beta), yp + 6, w - mr + 6, yp + 4,
                    exp(pooled$beta), pooled$lo, pooled$hi))
  }
  sprintf("<svg viewBox='0 0 %d %d' width='100%%' role='img'>%s%s<text x='%d' y='%.1f' text-anchor='end' class='lbl' font-weight='600'>Combinado</text>%s</svg>",
          w, h, g, rows, lw - 6, yp + 4, pool)
}

.lox_css_central <- ".grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:14px}.card{border:1px solid var(--ln);border-radius:8px;padding:8px 10px;background:var(--card)}.card h4{margin:.1em 0 .3em;font-size:.88rem}
.card .sub{color:var(--mu);font-size:.78rem;margin-bottom:2px}svg .sq{fill:var(--mu)}svg .dia{fill:var(--c2)}svg .pi{stroke:var(--c2);stroke-width:1.4;stroke-dasharray:4 3}.hi{color:var(--warn);font-weight:600}"

#' Informe central en HTML (consenso + meta-analisis)
#'
#' @param central_dir carpeta de salida del central: contiene `consensus/` y `meta/`
#' @param out_file ruta del HTML (por defecto `<central_dir>/informe_central.html`)
#' @param top numero de arcos de incidencia con diagrama de bosque desplegado (el resto, plegado)
#' @return la ruta del fichero (invisible)
#' @export
lox_report_central <- function(central_dir, out_file = file.path(central_dir, "informe_central.html"), top = 12L) {
  cd <- file.path(central_dir, "consensus"); md <- file.path(central_dir, "meta")
  if (!file.exists(file.path(md, "meta_coefs.csv"))) stop("Falta meta/meta_coefs.csv en ", central_dir, " (ejecuta central-meta antes)", call. = FALSE)
  rd <- function(d, f) { p <- file.path(d, f); if (file.exists(p)) data.table::fread(p, na.strings = c("NA", "")) else NULL }
  mt <- rd(md, "meta_coefs.csv"); sc <- rd(md, "meta_site_coefs.csv"); si <- rd(md, "meta_sites.csv")
  cs <- rd(cd, "consensus_structure.csv"); mm <- if (file.exists(file.path(md, "manifest.json"))) jsonlite::read_json(file.path(md, "manifest.json")) else list()
  H <- character(); add <- function(...) H[length(H) + 1] <<- paste0(...)
  fx <- .lox_fmt; ci <- function(a, b, c) sprintf("%s (%s–%s)", fx(a, 2), fx(b, 2), fx(c, 2))
  nsite <- if (!is.null(si)) nrow(si) else NA
  min_ok <- mt[status == "ok"]

  add("<h1>Informe central · loxofed</h1><p class='muted'>Meta-análisis federado · generado ", .lox_esc(format(Sys.time(), "%Y-%m-%d %H:%M")),
      if (!is.null(mm$package_version)) paste0(" · loxofed ", .lox_esc(mm$package_version)), if (!is.null(mm$fingerprint)) paste0(" · huella <code>", .lox_esc(substr(mm$fingerprint, 1, 8)), "</code>"), "</p>")
  kp <- function(v, l) sprintf("<div class='kpi'><b>%s</b><span>%s</span></div>", .lox_esc(v), .lox_esc(l))
  add("<div class='kpis'>", kp(nsite, "sitios"), kp(.lox_big(nrow(mt)), "coeficientes meta-analizados"),
      kp(.lox_big(nrow(min_ok)), "con resultado combinado"), kp(.lox_big(nrow(mt) - nrow(min_ok)), "con k insuficiente"), "</div>")
  add("<div class='note'><b>Método.</b> Efectos aleatorios con τ² por REML, IC95% de Hartung-Knapp (t con k−1 gl), intervalo de predicción (t con k−2 gl, k≥3) ",
      "y, como sensibilidad, DerSimonian-Laird y el rango dejando un sitio fuera. Con pocos sitios, el IC normal es demasiado estrecho; por eso el principal es Hartung-Knapp. ",
      "Los OR (incidencia) y HR (mortalidad) se combinan en escala logar\u00edtmica. Con menos sitios que el m\u00ednimo configurado (<code>meta.min_sites</code>, 3 por defecto) no se da resultado combinado.</div>")

  # --- sitios
  add("<h2>1. Sitios participantes</h2>")
  if (!is.null(si)) {
    t <- copy(si); for (cc in c("n_eligible", "n_cohort")) if (cc %in% names(t)) t[[cc]] <- .lox_big(t[[cc]])
    if ("sample_weight" %in% names(t)) t$sample_weight <- fx(t$sample_weight, 2)
    add(.lox_html_table(t, num = intersect(c("n_eligible", "n_cohort", "sample_weight", "n_coefs"), names(t))))
    if (!is.null(si$version) && length(unique(si$version)) > 1) add("<div class='note warn'><b>Versiones distintas</b> del paquete entre sitios.</div>")
    if ("sample_weight" %in% names(si) && any(si$sample_weight > 1.0001, na.rm = TRUE))
      add("<p class='muted'>Los sitios con <i>peso de muestreo</i> &gt; 1 trabajaron con una muestra aleatoria de personas: sus errores estándar son mayores y pesan menos, sin sesgo.</p>")
  }

  # --- consenso
  add("<h2>2. Estructura consenso</h2>")
  if (!is.null(cs) && nrow(cs)) {
    ok <- cs[state %in% c("concordante", "reorientada")]
    add("<p>", nrow(ok), " arcos causales de ", nrow(cs), " candidatos (concordantes + reorientados por la precedencia de fechas; los inciertos no entran en el ajuste).</p>")
    add(.lox_svg_network(data.frame(from = ok$dbn_from, to = ok$dbn_to, strength = ok$strength)))
    tb <- cs[order(-frac_sites, -strength)][, .(Desde = .lox_short(dbn_from), Hasta = .lox_short(dbn_to), Fuerza = fx(strength, 2), Estado = state,
                                                `% precedencia` = pct_dir_final, `Sitios` = n_sites, `Fracción` = fx(frac_sites, 2))]
    add("<details><summary>Arcos y validación de dirección</summary>", .lox_html_table(tb, num = c("Fuerza", "% precedencia", "Sitios", "Fracción")), "</details>")
  } else add("<p class='muted'>Sin consenso disponible en <code>consensus/</code>.</p>")

  # --- mortalidad
  hz <- mt[capa == "mortalidad" & nodo != "muerte"]
  add("<h2>3. Mortalidad: efecto de cada enfermedad (HR combinado)</h2>")
  hk <- hz[status == "ok"][order(-exp_re)]
  if (nrow(hk)) {
    add("<p class='muted'>Punto = HR combinado; barra = IC95% de Hartung-Knapp. Naranja: IC por encima de 1; azul: por debajo; gris: incluye 1.</p>")
    add(.lox_svg_forest(.lox_short(hk$nodo), hk$exp_re, hk$lo, hk$hi))
    tb <- hk[, .(Estado = .lox_short(nodo), k, HR = fx(exp_re, 2), `IC95% HK` = sprintf("%s–%s", fx(lo, 2), fx(hi, 2)),
                 `Predicción 95%` = ifelse(is.na(pi_lo), NA, sprintf("%s–%s", fx(pi_lo, 2), fx(pi_hi, 2))), `τ²` = fx(tau2, 3), `I² %` = fx(100 * I2, 0),
                 `Sin 1 sitio` = ifelse(is.na(loo_min), NA, sprintf("%s–%s", fx(loo_min, 2), fx(loo_max, 2))))]
    add("<details><summary>Tabla con heterogeneidad y sensibilidad</summary>", .lox_html_table(tb, num = c("k", "HR", "τ²", "I² %")), "</details>")
  } else add("<p class='muted'>Sin estados con k suficiente.</p>")

  # --- incidencia
  inc <- mt[capa == "incidencia" & !grepl("^edad_|^Sexo", termino)]
  add("<h2>4. Incidencia: efecto de las enfermedades previas (OR combinado)</h2>")
  ik <- inc[status == "ok"][order(-exp_re)]
  if (nrow(ik)) {
    hetq <- ik[I2 >= 0.75]
    add("<p>", nrow(ik), " arcos con resultado combinado; <b>", nrow(hetq), "</b> con heterogeneidad alta (I² ≥ 75%), donde el efecto medio debe leerse junto con el intervalo de predicción.</p>")
    tb <- ik[, .(Desde = .lox_short(sub("_t0$", "", termino)), Hasta = .lox_short(nodo), k, OR = fx(exp_re, 2), `IC95% HK` = sprintf("%s–%s", fx(lo, 2), fx(hi, 2)),
                 `Predicción 95%` = ifelse(is.na(pi_lo), NA, sprintf("%s–%s", fx(pi_lo, 2), fx(pi_hi, 2))), `τ²` = fx(tau2, 3), `I² %` = fx(100 * I2, 0))]
    add(.lox_html_table(tb, num = c("k", "OR", "τ²", "I² %")))
    # heterogeneidad
    hh <- ik[order(-I2)][seq_len(min(.N, 20))]
    add("<h3>Heterogeneidad (I²) de los 20 arcos más heterogéneos</h3>")
    add(.lox_svg_bars(paste(.lox_short(sub("_t0$", "", hh$termino)), "→", .lox_short(hh$nodo)), 100 * hh$I2, fmt = function(v) sprintf("%.0f%%", v)))
    # fragilidad
    fr <- ik[!is.na(loo_min) & ((lo > 1 & loo_min < 1) | (hi < 1 & loo_max > 1))]
    add("<h3>Resultados frágiles</h3>")
    if (nrow(fr)) add("<p>El efecto deja de ser distinto de 1 al excluir un solo sitio en:</p>", .lox_html_table(fr[, .(Desde = .lox_short(sub("_t0$", "", termino)), Hasta = .lox_short(nodo),
         OR = fx(exp_re, 2), `Sin 1 sitio` = sprintf("%s–%s", fx(loo_min, 2), fx(loo_max, 2)))], num = "OR"))
    else add("<p class='muted'>Ningún arco con resultado significativo cambia de lado al dejar un sitio fuera.</p>")
    # sensibilidad DL
    dd <- ik[abs(beta_dl - beta_re) > 0.05 | (tau2_dl < tau2 / 2 & tau2 > 0)]
    add("<p class='muted'>Sensibilidad DerSimonian-Laird: ", nrow(ik) - nrow(dd), " de ", nrow(ik), " arcos dan el mismo efecto (Δβ &lt; 0,05); DL estima en media τ² = ", fx(mean(ik$tau2_dl), 3), " frente a ", fx(mean(ik$tau2), 3), " por REML.</p>")

    # forest por arco
    add("<h3>Diagramas de bosque por arco</h3><p class='muted'>Cada cuadrado es un sitio (IC95%); el diamante es el efecto combinado y la línea discontinua el intervalo de predicción. Ordenados por heterogeneidad.</p>")
    ord <- ik[order(-I2, -abs(beta_re))]
    card <- function(i) {
      r <- ord[i]; s <- sc[capa == "incidencia" & nodo == r$nodo & termino == r$termino & is.finite(beta) & is.finite(se) & se > 0][order(site)]
      sprintf("<div class='card'><h4>%s → %s</h4><div class='sub'>k = %d · I² = %s%% · τ² = %s</div>%s</div>", .lox_esc(.lox_short(sub("_t0$", "", r$termino))), .lox_esc(.lox_short(r$nodo)),
              r$k, fx(100 * r$I2, 0), fx(r$tau2, 3),
              .lox_svg_forest_arc(s$site, s$beta, s$se, list(beta = r$beta_re, lo = r$lo, hi = r$hi, pi_lo = r$pi_lo, pi_hi = r$pi_hi)))
    }
    n1 <- min(top, nrow(ord)); add("<div class='grid2'>", paste(vapply(seq_len(n1), card, ""), collapse = ""), "</div>")
    if (nrow(ord) > n1) add("<details><summary>Resto de arcos (", nrow(ord) - n1, ")</summary><div class='grid2'>", paste(vapply(seq(n1 + 1, nrow(ord)), card, ""), collapse = ""), "</div></details>")
  } else add("<p class='muted'>Sin arcos con k suficiente.</p>")

  ins <- mt[status != "ok"]
  if (nrow(ins)) add("<h2>5. Coeficientes sin resultado combinado</h2><p>", nrow(ins), " coeficientes con menos sitios que el mínimo; se muestran sus estimaciones por sitio en <code>meta_site_coefs.csv</code>.</p>")
  add("<footer>Informe generado a partir de tablas agregadas de los sitios (celdas &lt; 10 suprimidas). No contiene datos individuales.</footer>")
  html <- paste0("<!doctype html><html lang='es'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Informe central loxofed</title><style>",
                 .lox_css, .lox_css_central, "</style></head><body><main>", paste(H, collapse = "\n"), "</main></body></html>")
  writeLines(html, out_file, useBytes = TRUE); .lox_msg("Informe central -> %s", out_file); invisible(out_file)
}
