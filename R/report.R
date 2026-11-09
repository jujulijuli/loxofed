# Informe local en HTML autocontenido (R puro: sin Quarto/pandoc, sin red) -----------------------------
# Lee SOLO las tablas agregadas que ya escribieron las rondas (con supresion de celdas pequenas), asi que el
# informe se puede compartir igual que las carpetas round1/round2.

.lox_esc <- function(x) { x <- as.character(x); x[is.na(x)] <- "—"
  x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }

.lox_html_table <- function(df, num = character(), class = "") {
  if (is.null(df) || !nrow(df)) return("<p class='muted'>Sin datos.</p>")
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  th <- paste0("<th", ifelse(names(df) %in% num, " class='n'", ""), ">", .lox_esc(names(df)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(df)), function(i)
    paste0("<tr>", paste0("<td", ifelse(names(df) %in% num, " class='n'", ""), ">", .lox_esc(unlist(df[i, ], use.names = FALSE)), "</td>", collapse = ""), "</tr>"), "")
  sprintf("<div class='scroll'><table class='%s'><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>", class, th, paste(rows, collapse = ""))
}

.lox_fmt <- function(x, d = 2) ifelse(is.na(x), NA, formatC(x, format = "f", digits = d, big.mark = ""))
.lox_big <- function(x) ifelse(is.na(x), NA, format(round(x), big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE))


# Cuadro didactico: explicacion en lenguaje llano junto al contenido tecnico
.lox_plain <- function(...) {
  x <- c(...); x <- x[!is.na(x) & nzchar(x)]
  paste0("<div class='plain'><b>Cómo leerlo</b>", paste0("<p>", x, "</p>", collapse = ""), "</div>")
}
.lox_vs <- function(r, pos = "más", neg = "menos") { v <- round(100 * abs(r - 1)); if (v == 0) "prácticamente igual" else sprintf("un %d%% %s", v, if (r < 1) neg else pos) }
.lox_vs_dx <- function(r) { v <- round(100 * abs(r - 1)); if (v == 0) "prácticamente los diagnósticos esperados" else sprintf("un %d%% %s de diagnósticos que los esperados", v, if (r < 1) "menos" else "más") }
.lox_dec <- function(x, d = 2) sub("\\.", ",", .lox_fmt(x, d))      # decimal con coma en las frases

# Observadas / esperadas (O/E) con IC95% de Poisson exacto (Garwood) sobre las muertes observadas
.lox_oe <- function(obs, by, e = "n_exp") {
  if (is.null(obs) || !all(c("n_death", e) %in% names(obs))) return(NULL)
  obs <- data.table::copy(obs); obs[, ee := get(e)]
  o <- obs[!is.na(n_death) & !is.na(ee), .(O = sum(n_death), E = sum(ee)), by = by]
  o <- o[E >= 5]                       # grupos con menos de 5 esperadas no son interpretables
  if (!nrow(o)) return(NULL)
  o[, `:=`(oe = O / E, lo = ifelse(O > 0, stats::qchisq(.025, 2 * O) / 2, 0) / E, hi = stats::qchisq(.975, 2 * (O + 1)) / 2 / E)]
  o[]
}

# Etiquetas cortas: misma regla que `corto()` del informe de Navarra (sin colisiones en las 60 del catalogo)
.lox_short <- function(x) {
  x <- gsub("_AND_.*|_EMPHYSEMA.*|_VISUAL.*|_AND_OTHER.*|_RELATED.*", "", x)
  x <- gsub("CHRONIC_|_DISEASES?$|_DISEASE$", "", x)
  x <- gsub("_", " ", x); x <- tools::toTitleCase(tolower(x))
  sub("^Copd$", "COPD", x)
}

# SVG: bosque (estimacion + IC95%, escala log, linea de referencia en 1)
.lox_svg_forest <- function(lab, est, lo, hi, w = 640) {
  n <- length(est); if (!n) return("")
  rowh <- 20; mt <- 26; mb <- 24; lw <- 190; mr <- 120; h <- mt + rowh * n + mb
  lo <- pmax(lo, 0.05); hi <- pmin(hi, 50)
  r <- range(log(c(lo, hi, 1))); r <- r + c(-1, 1) * 0.04 * diff(r)
  X <- function(v) lw + (w - lw - mr) * (log(v) - r[1]) / diff(r)
  tk0 <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12); tk <- tk0[tk0 >= exp(r[1]) & tk0 <= exp(r[2])]
  g <- paste0(sprintf("<line x1='%.1f' x2='%.1f' y1='%d' y2='%d' class='%s'/><text x='%.1f' y='%d' text-anchor='middle' class='lbl'>%s</text>",
                      X(tk), X(tk), mt - 6, h - mb, ifelse(tk == 1, "ref", "grid"), X(tk), mt - 10, format(tk)), collapse = "")
  y <- mt + rowh * (seq_len(n) - 0.5)
  rows <- paste0(sprintf("<text x='%d' y='%.1f' text-anchor='end' class='lbl'>%s</text><line x1='%.1f' x2='%.1f' y1='%.1f' y2='%.1f' class='ci'/><circle cx='%.1f' cy='%.1f' r='3.6' class='%s'/><text x='%d' y='%.1f' class='lbl mu'>%.2f (%.2f\u2013%.2f)</text>",
                         lw - 8, y + 4, .lox_esc(lab), X(lo), X(hi), y, y, X(est), y, ifelse(lo > 1, "pt-up", ifelse(hi < 1, "pt-dn", "pt-ns")),
                         w - mr + 8, y + 4, est, lo, hi), collapse = "")
  ax <- sprintf("<text x='%.1f' y='%d' text-anchor='middle' class='lbl mu'>escala logar\u00edtmica (1 = sin diferencia)</text>", (lw + w - mr) / 2, h - 4)
  sprintf("<svg viewBox='0 0 %d %d' width='100%%' role='img'>%s%s%s</svg>", w, h, g, rows, ax)
}

# SVG: barras horizontales
.lox_svg_bars <- function(lab, val, fmt = function(v) sprintf("%.1f%%", v), w = 640) {
  n <- length(val); if (!n) return("")
  h <- 22 * n + 8; lw <- 290; lab <- ifelse(nchar(lab) > 44, paste0(substr(lab, 1, 43), "\u2026"), lab); mx <- max(val, 1e-9)
  bars <- vapply(seq_len(n), function(i) {
    y <- 4 + 22 * (i - 1); bw <- (w - lw - 70) * val[i] / mx
    sprintf("<text x='%d' y='%d' text-anchor='end' class='lbl'>%s</text><rect x='%d' y='%d' width='%.1f' height='14' rx='2' class='bar'/><text x='%.1f' y='%d' class='lbl'>%s</text>",
            lw - 6, y + 12, .lox_esc(lab[i]), lw, y + 1, bw, lw + bw + 5, y + 12, fmt(val[i]))
  }, "")
  sprintf("<svg viewBox='0 0 %d %d' width='100%%' role='img'>%s</svg>", w, h, paste(bars, collapse = ""))
}

# SVG: lineas (tasa por grupo de edad, una serie por sexo)
.lox_svg_lines <- function(x, series, ylab = "", w = 640, h = 280, ymin = 0, ref = NULL) {
  ml <- 56; mb <- 34; mt <- 14; mr <- 110
  ys <- unlist(lapply(series, `[[`, "y")); ys <- ys[is.finite(ys)]
  if (!length(ys)) return("<p class='muted'>Sin datos suficientes.</p>")
  ymax <- if (ymin == 0) max(ys) * 1.08 else max(c(ys, ref)) + 0.05 * diff(range(c(ys, ref, ymin)))
  X <- function(i) ml + (w - ml - mr) * (i - 1) / max(length(x) - 1, 1)
  Y <- function(v) mt + (h - mt - mb) * (1 - (v - ymin) / (ymax - ymin))
  tk <- pretty(c(ymin, ymax), 4); tk <- tk[tk <= ymax & tk >= ymin]
  grid <- paste0(sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' class='grid'/><text x='%d' y='%.1f' text-anchor='end' class='lbl'>%s</text>",
                         ml, w - mr, Y(tk), Y(tk), ml - 6, Y(tk) + 4, format(tk)), collapse = "")
  if (!is.null(ref)) grid <- paste0(grid, sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#c0392b' stroke-dasharray='4 3'/>", ml, w - mr, Y(ref), Y(ref)))
  xl <- paste0(sprintf("<text x='%.1f' y='%d' text-anchor='middle' class='lbl'>%s</text>", X(seq_along(x)), h - 12, .lox_esc(x)), collapse = "")
  cl <- c("c1", "c2", "c3", "c4"); out <- ""
  for (k in seq_along(series)) {
    s <- series[[k]]; ok <- is.finite(s$y)
    if (sum(ok) > 1) out <- paste0(out, sprintf("<polyline points='%s' class='ln %s' style='fill:none'/>",
                                                paste(sprintf("%.1f,%.1f", X(seq_along(x))[ok], Y(s$y[ok])), collapse = " "), cl[k]))
    out <- paste0(out, paste0(sprintf("<circle cx='%.1f' cy='%.1f' r='3' class='pt %s'/>", X(seq_along(x))[ok], Y(s$y[ok]), cl[k]), collapse = ""),
                  sprintf("<text x='%d' y='%d' class='lbl %s'>%s</text>", w - mr + 8, 30 + 16 * k, cl[k], .lox_esc(s$name)))
  }
  sprintf("<svg viewBox='0 0 %d %d' width='100%%' role='img'>%s%s%s<text x='12' y='%d' class='lbl' transform='rotate(-90 12 %d)' text-anchor='middle'>%s</text></svg>",
          w, h, grid, xl, out, h / 2, h / 2, .lox_esc(ylab))
}

# SVG: red en capas (izquierda -> derecha segun la cadena de progresion; arcos from -> to)
.lox_svg_network <- function(edges, w = 980) {
  if (is.null(edges) || !nrow(edges)) return("<p class='muted'>Sin arcos.</p>")
  nd <- sort(unique(c(edges$from, edges$to))); n <- length(nd)
  ly <- setNames(rep(0L, n), nd)
  for (it in seq_len(n)) { ch <- FALSE                                 # capa = camino mas largo desde las fuentes (con tope por ciclos)
    for (i in seq_len(nrow(edges))) { v <- min(ly[[edges$from[i]]] + 1L, n)
      if (v > ly[[edges$to[i]]]) { ly[[edges$to[i]]] <- v; ch <- TRUE } }
    if (!ch) break }
  L <- sort(unique(ly)); ly <- setNames(match(ly, L), nd); nl <- max(ly)
  pos <- list(); ord <- list()
  for (k in seq_len(nl)) { m <- nd[ly == k]
    if (k > 1) { bc <- vapply(m, function(v) { f <- edges$from[edges$to == v]; mean(match(f, unlist(ord[seq_len(k - 1)]), nomatch = 0)) }, 0); m <- m[order(bc, m)] }
    ord[[k]] <- m }
  rowh <- 24; mx <- max(lengths(ord)); h <- 30 + rowh * mx
  colw <- (w - 150) / max(nl - 1, 1)
  px <- setNames(numeric(n), nd); py <- px
  for (k in seq_len(nl)) { m <- ord[[k]]; off <- (mx - length(m)) * rowh / 2
    px[m] <- 20 + colw * (k - 1) * (if (nl > 1) 1 else 0); py[m] <- 24 + off + rowh * (seq_along(m) - 1) }
  ar <- vapply(seq_len(nrow(edges)), function(i) {
    a <- edges$from[i]; b <- edges$to[i]; dx <- px[b] - px[a]; dy <- py[b] - py[a]; d <- sqrt(dx^2 + dy^2); if (d == 0) return("")
    sprintf("<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' class='edge' stroke-width='%.1f' marker-end='url(#ah)'/>",
            px[a] + dx / d * 6, py[a] + dy / d * 6, px[b] - dx / d * 9, py[b] - dy / d * 9,
            0.8 + 5 * pmax(0, ((if (is.null(edges$strength)) 1 else edges$strength[i]) - 0.6) / 0.4))
  }, "")
  lab <- vapply(seq_len(n), function(i)
    sprintf("<circle cx='%.1f' cy='%.1f' r='5' class='node'/><text x='%.1f' y='%.1f' class='lbl halo'>%s</text>",
            px[i], py[i], px[i] + 9, py[i] - 7, .lox_esc(.lox_short(nd[i]))), "")
  sprintf("<svg viewBox='0 0 %d %.0f' width='100%%' role='img'><defs><marker id='ah' viewBox='0 0 10 10' refX='9' refY='5' markerUnits='userSpaceOnUse' markerWidth='9' markerHeight='9' orient='auto'><path d='M0,0 L10,5 L0,10 z' class='arrow'/></marker></defs>%s%s</svg>",
          w, h, paste(ar, collapse = ""), paste(lab, collapse = ""))
}

.lox_css <- "
:root{--bg:#fff;--fg:#1c2733;--mu:#667585;--ln:#d9e0e7;--ac:#1f6fb2;--c1:#1f6fb2;--c2:#d9772b;--c3:#2e8b57;--c4:#8a4fbf;--card:#f6f8fa;--warn:#b26a00;--err:#b3261e}
@media(prefers-color-scheme:dark){:root{--bg:#14181d;--fg:#e4e9ee;--mu:#98a5b3;--ln:#2d3640;--ac:#5aa9e6;--c1:#5aa9e6;--c2:#f0a35e;--c3:#6cc58f;--c4:#b992e6;--card:#1c2229;--warn:#e0a040;--err:#f2867f}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
main{max-width:960px;margin:0 auto;padding:24px 16px 64px}h1{font-size:1.6rem;margin:.2em 0}h2{font-size:1.2rem;margin:2em 0 .5em;border-bottom:1px solid var(--ln);padding-bottom:.25em}h3{font-size:1rem;margin:1.4em 0 .4em}
.muted{color:var(--mu)}.gl dt{font-weight:600;margin-top:.6em}.gl dd{margin:.1em 0 0 1em;color:var(--mu)}.plain{background:rgba(31,111,178,.07);border-left:4px solid var(--c1);border-radius:6px;padding:8px 16px;margin:14px 0}.plain>b{display:block;color:var(--c1);font-size:.85rem;text-transform:uppercase;letter-spacing:.04em}.plain p{margin:.45em 0}.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0}.kpi{background:var(--card);border:1px solid var(--ln);border-radius:8px;padding:10px 12px}.kpi b{display:block;font-size:1.3rem}.kpi span{color:var(--mu);font-size:.82rem}
.note{background:var(--card);border-left:3px solid var(--ac);padding:8px 12px;margin:12px 0;border-radius:4px}.note.warn{border-color:var(--warn)}
.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:.88rem}th,td{padding:4px 8px;border-bottom:1px solid var(--ln);text-align:left;white-space:nowrap}th{color:var(--mu);font-weight:600}.n{text-align:right;font-variant-numeric:tabular-nums}
.OK{color:var(--c3)}.WARN{color:var(--warn);font-weight:600}.ERROR{color:var(--err);font-weight:700}
svg .lbl{font-size:11px;fill:var(--fg)}svg .bar{fill:var(--ac)}svg .grid{stroke:var(--ln);stroke-width:1}svg .ln{fill:none;stroke-width:2}svg .c1{stroke:var(--c1);fill:var(--c1)}svg .c2{stroke:var(--c2);fill:var(--c2)}svg .c3{stroke:var(--c3);fill:var(--c3)}svg .c4{stroke:var(--c4);fill:var(--c4)}svg .pt{stroke:none}svg text.c1,svg text.c2,svg text.c3,svg text.c4{stroke:none}
svg .ref{stroke:var(--fg);stroke-width:1.2}svg .ci{stroke:var(--mu);stroke-width:1.6}svg .pt-up{fill:var(--c2)}svg .pt-dn{fill:var(--c1)}svg .pt-ns{fill:var(--mu)}svg .mu{fill:var(--mu)}
svg .halo{paint-order:stroke;stroke:var(--bg);stroke-width:3px;stroke-linejoin:round}svg .edge{stroke:var(--mu);fill:none;opacity:.8}svg .arrow{fill:var(--mu)}svg .node{fill:var(--ac)}
footer{margin-top:3em;color:var(--mu);font-size:.8rem}"

#' Informe local en HTML (R puro, autocontenido, sin red)
#'
#' Resume lo que produjeron las rondas en este sitio: QC, cohorte (con nota de muestreo), prevalencias, estructura
#' aprendida (y consenso si se da), efectos de las enfermedades sobre la incidencia y la mortalidad, y mortalidad
#' observada. Usa solo las tablas agregadas ya escritas (con supresion de celdas pequenas), asi que es tan
#' compartible como `round1/` y `round2/`.
#' @param output_dir carpeta de salida del sitio (la que contiene `round1/` y `round2/`)
#' @param consensus_dir opcional: carpeta `consensus/` recibida del nodo central
#' @param out_file ruta del HTML; por defecto `<output_dir>/informe_local.html`
#' @return la ruta del fichero (invisible)
#' @export
lox_report_local <- function(output_dir, consensus_dir = NULL, out_file = file.path(output_dir, "informe_local.html")) {
  r1 <- file.path(output_dir, "round1"); r2 <- file.path(output_dir, "round2")
  rd <- function(d, f) { p <- file.path(d, f); if (file.exists(p)) data.table::fread(p, na.strings = c("NA", "")) else NULL }
  m1 <- if (file.exists(file.path(r1, "manifest.json"))) jsonlite::read_json(file.path(r1, "manifest.json")) else NULL
  m2 <- if (file.exists(file.path(r2, "manifest.json"))) jsonlite::read_json(file.path(r2, "manifest.json")) else NULL
  if (is.null(m1) && is.null(m2)) stop("No hay round1/ ni round2/ en ", output_dir, call. = FALSE)
  mf <- if (!is.null(m2)) m2 else m1
  qc <- rd(if (!is.null(m1)) r1 else r2, "qc.csv"); meta <- rd(r1, "r1_meta.csv"); prev <- rd(r1, "r1_prevalence.csv")
  ed1 <- rd(r1, "r1_edges.csv"); co <- rd(r2, "r2_coefs.csv"); obs <- rd(r2, "r2_observed.csv")
  cs <- if (!is.null(consensus_dir)) rd(consensus_dir, "consensus_structure.csv") else NULL
  H <- character()
  add <- function(...) H[length(H) + 1] <<- paste0(...)

  sid <- mf$site_id; wgt <- if (!is.null(meta) && "sample_weight" %in% names(meta)) meta$sample_weight[1] else 1
  nel <- if (!is.null(meta)) meta$n_eligible[1] else if (!is.null(m2$n_eligible)) m2$n_eligible else NA
  nco <- if (!is.null(meta)) meta$n_cohort[1] else if (!is.null(m2$n_cohort)) m2$n_cohort else NA
  add("<h1>Informe local · ", .lox_esc(sid), "</h1><p class='muted'>loxofed ", .lox_esc(mf$package_version),
      " · generado ", .lox_esc(format(Sys.time(), "%Y-%m-%d %H:%M")), " · huella de configuracion <code>",
      .lox_esc(substr(mf$fingerprint, 1, 8)), "</code>")
  mort <- if (!is.null(m2$mortality)) m2$mortality else NULL
  if (!is.null(mort)) add(" · mortalidad: <b>", if (mort == "exact") "fechas exactas (tiempo-persona)" else "anual t→t+1", "</b>")
  add("</p>")

  # --- que es este informe y donde encaja
  fase <- if (!is.null(m2)) "tras la <b>ronda 2</b>" else if (!is.null(cs)) "tras recibir el <b>consenso</b>, antes de la ronda 2" else "tras la <b>ronda 1</b>, antes del consenso"
  if (!is.null(cs) && "n_sites" %in% names(cs) && nrow(cs) && max(cs$n_sites, na.rm = TRUE) <= 1)
    fase <- paste0(fase, ", con una estructura <b>provisional</b> (el consenso se ha calculado solo con este sitio)")
  add("<h2>Qué es este informe y dónde encaja</h2><div class='plain'><p>Este análisis es <b>federado</b>: cada comunidad analiza sus propios datos en su servidor y solo se comparten tablas agregadas. Se hace en varios pasos:</p><ol>",
      "<li><b>Ronda 1 (en cada sitio).</b> Se comprueba la calidad de los datos y se aprende qué relaciones entre enfermedades se repiten. Se guardan tablas agregadas en la carpeta <code>round1/</code>, que se envía al centro.</li>",
      "<li><b>Consenso (en el centro).</b> Con las <code>round1/</code> de todos los sitios se construye una única estructura común de relaciones, que se devuelve a cada sitio.</li>",
      "<li><b>Ronda 2 (en cada sitio).</b> Con la estructura común se estiman los efectos de las enfermedades sobre la mortalidad y la incidencia, y se evalúa la calidad del modelo. Se guardan en <code>round2/</code>, que se envía al centro.</li>",
      "<li><b>Metaanálisis (en el centro).</b> Se combinan las estimaciones de los sitios y se contrasta si los efectos son consistentes entre comunidades.</li></ol>",
      "<p><b>Este informe se ha generado ", fase, ".</b> Es de uso interno del sitio: sirve para revisar los resultados propios y detectar problemas antes de enviar nada. ",
      "Lo que sale del sitio son solo las carpetas <code>round1/</code> y <code>round2/</code>, con celdas de menos de 10 personas suprimidas y ningún dato individual. ",
      "Los resultados de un único sitio no son todavía el resultado del proyecto: lo que se contrastará en el metaanálisis es si se repiten entre comunidades.</p></div>")

  # --- resumen
  k <- function(v, l) sprintf("<div class='kpi'><b>%s</b><span>%s</span></div>", .lox_esc(v), .lox_esc(l))
  nm <- if (!is.null(co)) co[capa == "mortalidad", max(n_events, na.rm = TRUE)] else NA
  add("<div class='kpis'>", k(.lox_big(nco), "personas en el análisis"), k(.lox_big(nel), "personas elegibles"),
      if (!is.null(prev)) k(.lox_big(nrow(prev[prev >= 0.01])), "nodos con prevalencia ≥ 1%") else "",
      if (!is.null(co)) k(.lox_big(nm), "muertes en el ajuste") else "", "</div>")
  if (isTRUE(wgt > 1.0001))
    add("<div class='note'><b>Muestra aleatoria:</b> se analizaron ", .lox_big(nco), " de ", .lox_big(nel), " personas elegibles (peso de muestreo ",
        .lox_fmt(wgt, 2), "). Los coeficientes se estiman <b>sin ponderar</b> (con pesos constantes no cambian y los errores estándar serían falsos). ",
        "Los recuentos marcados <i>est.</i> están reescalados al tamaño del sitio multiplicando por el peso.</div>")

  # --- resumen para lectores no tecnicos
  res <- character()
  yr <- if (!is.null(obs) && nrow(obs)) range(obs$anio) else NULL
  res <- c(res, sprintf("Se ha seguido a <b>%s personas</b>%s%s.", .lox_big(nco),
                        if (!is.null(yr)) sprintf(" entre %d y %d", yr[1], yr[2]) else "",
                        if (!is.null(co) && is.finite(nm)) sprintf(", entre las que se registraron <b>%s muertes</b>", .lox_big(nm)) else ""))
  if (!is.null(prev)) { t3 <- head(prev[!is.na(prev) & prev > 0][order(-prev)], 3)
    res <- c(res, sprintf("Las enfermedades más frecuentes al empezar eran <b>%s</b>.", paste0(.lox_short(t3$node), " (", .lox_fmt(100 * t3$prev, 0), "%)", collapse = ", "))) }
  if (!is.null(co)) { hz0 <- co[capa == "mortalidad" & nodo != "muerte" & is.finite(beta)][order(-beta)]
    if (nrow(hz0) >= 3) res <- c(res, sprintf("A igual edad y sexo, las enfermedades que más se asocian a morir antes son <b>%s</b> (el número es cuántas veces mayor es el riesgo).",
                                                paste0(.lox_short(hz0$nodo[1:3]), " (×", .lox_dec(exp(hz0$beta[1:3]), 1), ")", collapse = ", ")))
    inc0 <- co[capa == "incidencia" & !grepl("^edad_|^Sexo", termino)]; if (!is.null(cs) && nrow(cs)) res <- c(res,
      sprintf("Se han identificado <b>%d relaciones estables</b> entre enfermedades, en las que una suele preceder a la otra (sección 3).", nrow(cs[state %in% c("concordante", "reorientada")]))) }
  oe0 <- .lox_oe(obs, "grupo_edad")
  if (!is.null(oe0)) { dv <- max(abs(oe0$oe - 1))
    res <- c(res, sprintf("El modelo de mortalidad %s por grupos de edad: la mayor desviación entre las muertes reales y las que predice es de <b>%.0f%%</b> (sección 7).",
                          if (dv <= 0.08) "acierta bien" else if (dv <= 0.2) "se ajusta con desviaciones moderadas" else "se desvía de forma importante", 100 * dv)) }
  if (!is.null(co)) { pd0 <- co[nodo == "muerte" & termino == "pandemia"]
    if (nrow(pd0)) res <- c(res, sprintf("En 2020-2021 la mortalidad fue <b>%s</b> de lo que cabía esperar por la tendencia previa (IC95%% %s–%s) (sección 8).",
                                         .lox_vs(exp(pd0$beta), "más alta", "más baja"), .lox_dec(exp(pd0$beta - 1.96 * pd0$se), 2), .lox_dec(exp(pd0$beta + 1.96 * pd0$se), 2))) }
  dxy0 <- rd(r2, "r2_dx_by_year.csv"); pp0 <- if (!is.null(dxy0) && "pandemia" %in% names(dxy0)) lox_dx_pandemic(dxy0) else NULL
  if (!is.null(pp0) && nrow(pp0[node == "TOTAL"])) { tt0 <- pp0[node == "TOTAL"][order(anio)]
    res <- c(res, sprintf("Durante la pandemia se <b>diagnosticaron menos enfermedades de las esperadas</b>: %s (sección 9).",
                          paste0(tt0$anio, ": ", vapply(tt0$ratio, .lox_vs, ""), collapse = "; "))) }
  if (length(res)) add("<h2>Resumen para lectores no especializados</h2><div class='plain'><b>Lo esencial</b><ul>", paste0("<li>", res, "</li>", collapse = ""),
                       "</ul><p class='muted'>Más abajo, cada sección incluye el detalle técnico y un recuadro «Cómo leerlo» que explica cómo leerla. Al final hay un glosario.</p></div>")

  # --- QC
  if (!is.null(qc)) {
    add("<h2>1. Control de calidad de las entradas</h2>")
    add(.lox_plain("Antes de analizar nada, el programa comprueba que los ficheros están completos y bien formados: sin personas duplicadas, con fechas legibles, con sexos reconocidos y con códigos de enfermedad del catálogo.",
                   "<b>OK</b> = correcto. <b>WARN</b> = aviso que conviene revisar pero que no impide el análisis. <b>ERROR</b> = hay que corregirlo antes de seguir."))
    if ("exitus_perdidos_por_baja_anterior" %in% qc$check) {
      g <- function(ch) { v <- qc$value[qc$check == ch]; if (length(v)) v[1] else "NA" }
      add("<div class='note'><b>Fecha de baja y fecha de muerte.</b> Este sitio aporta fecha de baja. La regla del análisis es la misma en todos los sitios: ",
          "la baja <b>censura</b> el seguimiento (se deja de observar a la persona en esa fecha) y una muerte solo se cuenta si ocurre en la fecha de baja o antes. ",
          "Hay <b>", .lox_esc(g("baja_y_exitus")), "</b> personas con ambas fechas. En <b>", .lox_esc(g("exitus_con_baja_anterior")), "</b> la muerte es posterior a la baja y, por la regla, ",
          "no se cuenta; dentro del periodo de seguimiento son <b>", .lox_esc(g("exitus_perdidos_por_baja_anterior")), "</b> muertes (", .lox_esc(g("exitus_perdidos_pct_de_exitus_seguimiento")), "% de las muertes del seguimiento), ",
          "de las cuales <b>", .lox_esc(g("exitus_perdidos_dentro_de_30_dias_de_la_baja")), "</b> ocurren en los 30 días siguientes a la baja (posible retraso en registrar la baja por defunción). ",
          "En <b>", .lox_esc(g("exitus_con_baja_posterior")), "</b> casos la baja es posterior a la muerte; ahí la muerte sí cuenta. ",
          "Si el porcentaje de muertes perdidas es alto, la mortalidad observada del sitio estará infraestimada y conviene revisar el origen de ambas fechas.</div>")
    }
    nw <- sum(qc$flag != "OK")
    add("<p>", if (nw) sprintf("%d de %d comprobaciones requieren atención.", nw, nrow(qc)) else "Todas las comprobaciones correctas.", "</p>")
    q <- as.data.frame(qc); q$flag <- sprintf("<span class='%s'>%s</span>", q$flag, q$flag)
    t <- .lox_html_table(transform(q, flag = "@@"), num = "value"); # marcador para no escapar el span
    for (i in seq_len(nrow(q))) t <- sub("@@", q$flag[i], t, fixed = TRUE)
    t <- gsub("&lt;span class='(OK|WARN|ERROR)'&gt;(.*?)&lt;/span&gt;", "<span class='\\1'>\\2</span>", t)
    add(t)
  }

  # --- prevalencia
  if (!is.null(prev)) {
    add("<h2>2. Prevalencia inicial de las enfermedades</h2>")
    pp <- prev[!is.na(prev) & prev > 0][order(-prev)]; top <- head(pp, 25)
    add("<p class='muted'>Prevalencia en la cohorte analizada (a 1 de enero del primer año). Se muestran las 25 más frecuentes.</p>")
    add(.lox_svg_bars(.lox_short(top$node), 100 * top$prev))
    tb <- pp[, .(Enfermedad = .lox_short(node), `Prevalencia %` = .lox_fmt(100 * prev, 1), n = .lox_big(n_pos),
                 `n est.` = if ("n_pos_est" %in% names(pp)) .lox_big(n_pos_est) else NA)]
    add("<details><summary>Tabla completa</summary>", .lox_html_table(tb, num = c("Prevalencia %", "n", "n est.")), "</details>")
  }

  # --- estructura
  add("<h2>3. Estructura de progresión</h2>")
  add(.lox_plain("Cada <b>flecha</b> va de una enfermedad a otra que suele aparecer <i>después</i> en las mismas personas (por ejemplo, de la hipertensión a la insuficiencia cardiaca). Una relación se conserva solo si cumple dos condiciones:",
                 "<b>1. Es estable.</b> El programa repite el aprendizaje 100 veces con muestras aleatorias distintas de las personas. La <b>fuerza</b> de una flecha es la proporción de repeticiones en las que aparece (0,60 = 60 de 100). Se exige al menos 0,60 y que se cumpla en al menos la mitad de los sitios.",
                 "<b>2. Tiene un orden temporal coherente.</b> Se miran las fechas de diagnóstico de las personas que tienen las dos enfermedades. Si en al menos el 55% de los casos la primera se diagnostica antes que la segunda, la flecha se mantiene; si la mayoría va al revés, se invierte; si no hay una mayoría clara, se descarta.",
                 "Cuanto más gruesa es la flecha, mayor es su fuerza. Una flecha indica una asociación con orden temporal; <b>no demuestra</b> por sí sola que una enfermedad cause la otra."))
  if (!is.null(cs) && nrow(cs)) {
    ok <- cs[state %in% c("concordante", "reorientada")]
    nsit <- if ("n_sites" %in% names(cs) && nrow(cs)) max(cs$n_sites, na.rm = TRUE) else NA
    add("<p>", if (isTRUE(nsit <= 1)) "Estructura <b>provisional</b> (consenso calculado solo con este sitio; cuando se incorporen los demás sitios puede cambiar y habrá que repetir la ronda 2): "
              else sprintf("Estructura <b>consenso</b> (común a %s sitios): ", if (is.finite(nsit)) nsit else "todos los"), nrow(ok), " arcos (relaciones con orden temporal) de ", nrow(cs), " candidatos. ",
        "<span class='muted'>El grosor es la fuerza del arco.</span></p>")
    add(.lox_svg_network(data.frame(from = ok$dbn_from, to = ok$dbn_to, strength = ok$strength)))
    tb <- cs[order(-strength)][, .(Desde = .lox_short(dbn_from), Hasta = .lox_short(dbn_to), Fuerza = .lox_fmt(strength, 2),
                                    Estado = state, `% precedencia` = pct_dir_final, `n pares` = n_orden)]
    add("<details><summary>Arcos del consenso y validación de dirección</summary>", .lox_html_table(tb, num = c("Fuerza", "% precedencia", "n pares")), "</details>")
  } else if (!is.null(ed1) && nrow(ed1)) {
    thr <- 0.6; ok <- ed1[strength >= thr]
    add("<p>Estructura aprendida <b>localmente</b> (fuerza ≥ ", thr, "): ", nrow(ok), " arcos. ",
        "Aún no se ha recibido el consenso; pasa <code>consensus_dir</code> a <code>lox_report_local()</code> para verlo.</p>")
    add(.lox_svg_network(ok))
    add("<details><summary>Todos los arcos locales</summary>",
        .lox_html_table(ed1[order(-strength)][, .(Desde = .lox_short(from), Hasta = .lox_short(to), Fuerza = .lox_fmt(strength, 2))], num = "Fuerza"), "</details>")
  } else add("<p class='muted'>Sin estructura disponible.</p>")

  # --- coeficientes
  if (!is.null(co)) {
    co[, `:=`(or = exp(beta), lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se))]
    est <- function(d) d[, .(Efecto = .lox_fmt(or, 2), `IC95%` = sprintf("%s–%s", .lox_fmt(lo, 2), .lox_fmt(hi, 2)))]
    hz <- co[capa == "mortalidad" & nodo != "muerte"][order(-or)]
    add("<h2>4. Efecto de cada enfermedad sobre la mortalidad</h2>")
    add(.lox_plain("El <b>HR</b> (cociente de tasas) compara el riesgo de morir de una persona que tiene la enfermedad con el de otra de la misma edad y sexo que no la tiene, teniendo en cuenta las demás enfermedades. <b>HR = 1</b>: sin diferencia. <b>HR = 2</b>: el doble de riesgo. <b>HR = 0,8</b>: un 20% menos.",
                   "La barra es el <b>intervalo de confianza al 95%</b>: el rango de valores compatible con los datos. Si la barra cruza el 1, no se puede afirmar que haya diferencia.",
                   if (nrow(hz) >= 3) sprintf("En este sitio, las enfermedades con mayor riesgo son %s.", paste0("<b>", .lox_short(hz$nodo[1:3]), "</b> (HR ", .lox_dec(hz$or[1:3], 2), ")", collapse = ", "))))
    add("<p>", if (identical(mort, "exact")) "Cociente de tasas (HR) ajustado por edad y sexo, estado de enfermedad variable en el tiempo."
        else "Cociente de riesgos ajustado por edad y sexo (convención anual).", " Solo estados con eventos suficientes.</p>")
    if (nrow(hz)) {
      add("<p class='muted'>Naranja: IC95% por encima de 1; azul: por debajo; gris: incluye 1. Escala logar\u00edtmica.</p>")
      add(.lox_svg_forest(.lox_short(hz$nodo), hz$or, hz$lo, hz$hi))
      add("<details><summary>Tabla</summary>", .lox_html_table(cbind(Estado = .lox_short(hz$nodo), est(hz), `Muertes (modelo)` = .lox_big(hz$n_events)),
                                                           num = c("Efecto", "Muertes (modelo)")), "</details>")
    } else add("<p class='muted'>Sin estados ajustados.</p>")
    bsx <- rd(r2, "r2_hr_por_sexo.csv")
    if (!is.null(bsx) && nrow(bsx)) {
      w <- data.table::dcast(bsx, nodo ~ sexo, value.var = c("beta", "se"))
      if (all(c("beta_Hombre", "beta_Mujer", "se_Hombre", "se_Mujer") %in% names(w))) {
        w <- w[is.finite(beta_Hombre) & is.finite(beta_Mujer)]
        w[, `:=`(d = beta_Hombre - beta_Mujer, sed = sqrt(se_Hombre^2 + se_Mujer^2))]; w[, `:=`(rhr = exp(d), lo = exp(d - 1.96 * sed), hi = exp(d + 1.96 * sed))]
        w <- w[order(-rhr)]
        add("<h3>¿Es igual el efecto en hombres y en mujeres?</h3>")
        add(.lox_plain("Se ha repetido el modelo de mortalidad por separado en hombres y en mujeres. El <b>cociente de HR (hombres / mujeres)</b> compara el efecto relativo de cada enfermedad entre sexos: <b>1</b> = el mismo efecto; <b>mayor que 1</b>, la enfermedad se asocia a una mortalidad relativa más alta en hombres; <b>menor que 1</b>, en mujeres. Si la barra cruza el 1, no se puede afirmar que el efecto difiera.",
                       "No compara cuántos hombres y mujeres mueren (eso lo recoge el término de sexo del modelo), sino si <i>la misma enfermedad</i> pesa más en uno de los dos."))
        add(.lox_svg_forest(.lox_short(w$nodo), w$rhr, w$lo, w$hi))
        add("<details><summary>Tabla por sexo</summary>", .lox_html_table(w[, .(Estado = .lox_short(nodo), `HR hombres` = .lox_fmt(exp(beta_Hombre), 2), `HR mujeres` = .lox_fmt(exp(beta_Mujer), 2),
                                                                             `Cociente H/M` = .lox_fmt(rhr, 2), `IC95%` = sprintf("%s–%s", .lox_fmt(lo, 2), .lox_fmt(hi, 2)))],
                                                        num = c("HR hombres", "HR mujeres", "Cociente H/M")), "</details>")
      }
    }
    add("<h2>5. Incidencia: efecto de las enfermedades previas</h2>")
    add(.lox_plain("Pregunta: <i>si una persona ya tiene la enfermedad A, ¿cuánto más probable es que desarrolle la B en el año siguiente?</i> El <b>odds ratio (OR)</b> lo mide, comparando con personas de la misma edad y sexo que no tienen A. <b>OR = 1</b>: ninguna diferencia. <b>OR = 2</b>: las probabilidades se duplican. <b>OR = 0,5</b>: se reducen a la mitad.",
                   "Para enfermedades poco frecuentes el OR se parece mucho al aumento del riesgo; para las muy frecuentes tiende a exagerarlo."))
    inc <- co[capa == "incidencia" & !grepl("^edad_|^Sexo", termino)]
    add("<p>Odds ratio de que aparezca la enfermedad (<i>Hasta</i>) en el año siguiente si la persona ya tiene la enfermedad <i>Desde</i>, ajustado por edad y sexo. ",
        "Solo arcos de la estructura consenso.</p>")
    if (nrow(inc)) {
      inc <- inc[order(-or)]
      add(.lox_html_table(cbind(Desde = .lox_short(sub("_t0$", "", inc$termino)), Hasta = .lox_short(inc$nodo), est(inc),
                                `Eventos` = .lox_big(inc$n_events)), num = c("Efecto", "Eventos")))
    } else add("<p class='muted'>Sin arcos ajustados (falta el consenso o no hay padres).</p>")
    at <- rd(r2, "r2_arc_tables.csv")
    if (!is.null(at) && nrow(at) && !anyNA(at[, .(a1_b1, a1_b0, a0_b1, a0_b0)])) {
      at[, `:=`(r1 = a1_b1 / (a1_b1 + a1_b0), r0 = a0_b1 / (a0_b1 + a0_b0))]; at[, orc := (as.numeric(a1_b1) * a0_b0) / (as.numeric(a1_b0) * a0_b1)]
      add("<details><summary>Tablas 2×2 por arco (cifras sin ajustar)</summary><p class='muted'>Entre las personas que no tenían la enfermedad <i>Hasta</i>: proporción que la desarrolla en el año siguiente según tengan o no la enfermedad <i>Desde</i>. Es la comprobación más simple del OR del modelo.</p>",
          .lox_html_table(at[order(-orc), .(Desde = .lox_short(from), Hasta = .lox_short(to), `Con Desde %` = .lox_fmt(100 * r1, 1), `Sin Desde %` = .lox_fmt(100 * r0, 1), `OR sin ajustar` = .lox_fmt(orc, 2),
                                           `Personas-año con Desde` = .lox_big(a1_b1 + a1_b0))], num = c("Con Desde %", "Sin Desde %", "OR sin ajustar", "Personas-año con Desde")), "</details>")
    }
    if ("restriccion" %in% names(co)) {
      rs <- unique(co[!is.na(restriccion) & restriccion != "", .(nodo, restriccion)])
      if (nrow(rs)) add("<div class='note'><b>Enfermedades de un solo sexo:</b> en ", paste0(.lox_esc(rs$nodo), " (", sub("^solo_", "solo ", .lox_esc(rs$restriccion)), ")", collapse = ", "),
                        " casi todos los casos son de un sexo, así que el modelo de incidencia se ajusta solo en ese sexo y no estima el efecto del sexo.</div>")
    }
    bad <- co[converged %in% c(FALSE, "FALSE") | (!is.na(flag) & flag != "")]
    if (nrow(bad)) add("<div class='note warn'><b>Avisos del ajuste:</b> ", nrow(bad), " coeficientes con aviso o sin convergencia (ver <code>r2_coefs.csv</code>).</div>")
  }

  # --- mortalidad observada
  if (!is.null(obs) && nrow(obs)) {
    add("<h2>6. Mortalidad observada</h2>")
    add(.lox_plain("Una <b>persona-año</b> es una persona observada durante un año: la unidad de «tiempo vivido bajo observación». 100 personas observadas durante 6 meses suman 50 personas-año; 10 personas durante un año, 10 personas-año. La tasa es el número de muertes por cada 1.000 personas-año. Aquí se muestra la mortalidad <b>real</b>, sin ningún ajuste: sube mucho con la edad y suele ser mayor en hombres."))
    if ("person_years" %in% names(obs)) {
      g <- obs[!is.na(n_death), .(py = sum(person_years), d = sum(n_death)), by = .(Sexo, grupo_edad)]
      lv <- c("<50", "50-59", "60-69", "70-79", "80-89", "90+"); lv <- lv[lv %in% g$grupo_edad]
      ser <- lapply(unique(g$Sexo), function(s) { x <- g[Sexo == s]; y <- 1000 * x$d[match(lv, x$grupo_edad)] / x$py[match(lv, x$grupo_edad)]
                                                    list(name = as.character(s), y = y) })
      add("<p>Tasa bruta de mortalidad por 1.000 personas-año, por grupo de edad (cada grupo sin las celdas suprimidas).</p>")
      add(.lox_svg_lines(lv, ser, ylab = "muertes / 1.000 p-a"))
      tot <- obs[, .(`Personas-año` = sum(person_years), `Muertes` = if (anyNA(n_death)) NA else sum(n_death)), by = .(Sexo)]
      if ("person_years_est" %in% names(obs)) tot[, `Personas-año est.` := obs[, sum(person_years_est), by = Sexo]$V1]
      tot[, `Personas-año` := .lox_big(`Personas-año`)]
      if ("Personas-año est." %in% names(tot)) tot[, `Personas-año est.` := .lox_big(`Personas-año est.`)]
      add("<details><summary>Totales por sexo</summary>", .lox_html_table(tot, num = setdiff(names(tot), "Sexo")), "</details>")
      add("<p class='muted'>Observada, sin estandarizar por edad. La comparación con lo esperado por el modelo está en la sección 7.</p>")
    } else {
      g <- obs[!is.na(n_death), .(n = sum(n_risk), d = sum(n_death)), by = .(Sexo, grupo_edad)]
      lv <- c("<50", "50-59", "60-69", "70-79", "80-89", "90+"); lv <- lv[lv %in% g$grupo_edad]
      ser <- lapply(unique(g$Sexo), function(s) { x <- g[Sexo == s]; list(name = as.character(s), y = 1000 * x$d[match(lv, x$grupo_edad)] / x$n[match(lv, x$grupo_edad)]) })
      add("<p>Riesgo anual de muerte por 1.000 personas, por grupo de edad.</p>"); add(.lox_svg_lines(lv, ser, ylab = "muertes / 1.000"))
    }
  }


  # --- calibracion: observadas frente a esperadas
  oe_age <- .lox_oe(obs, "grupo_edad")
  if (!is.null(oe_age)) {
    lv <- c("<50", "50-59", "60-69", "70-79", "80-89", "90+"); oe_age <- oe_age[match(lv[lv %in% grupo_edad], grupo_edad)]
    add("<h2>7. Calibración: muertes observadas frente a esperadas</h2>")
    oy <- .lox_oe(obs, "anio"); bad_y <- if (!is.null(oy)) oy[abs(oe - 1) > 0.10][order(-abs(oe - 1))] else NULL
    add(.lox_plain("<b>Qué modelo es.</b> Es un modelo de Poisson que estima la <b>tasa de mortalidad</b> de cada persona en cada tramo de tiempo a partir de su edad (con una curva flexible), su sexo y las enfermedades que tiene en ese momento (que cambian cuando se las diagnostican). Se ajusta una sola vez con todas las personas del sitio.",
                   "<b>Cómo se obtienen las muertes esperadas.</b> Para cada grupo (por ejemplo, mujeres de 70-79 años en un año concreto) se suman las tasas del modelo por el tiempo que sus personas estuvieron en observación. Esa suma son las muertes que el modelo <b>esperaría</b>.",
                   "<b>Qué se compara.</b> El cociente <b>O/E</b> compara las muertes reales (<b>O</b>bservadas) con las esperadas (<b>E</b>): <b>O/E = 1</b> significa que el modelo acierta; <b>1,10</b>, que hubo un 10% más de muertes de las previstas; <b>0,90</b>, un 10% menos.",
                   "Sirve para saber si podemos fiarnos del modelo en cada grupo de edad y en cada año. Más abajo se comparan dos versiones del modelo: <b>sin calendario</b> (ignora en qué año ocurre cada cosa) y <b>con calendario</b> (añade una curva por año y un efecto específico de 2020-2021). Si las dos dan lo mismo, el calendario no cambia nada; si difieren, se ve en qué años.",
                   if (!is.null(bad_y) && nrow(bad_y)) sprintf("Los años en que el modelo se aleja más de la realidad (más de un 10%%) son %s: ahí sus predicciones son menos fiables.",
                        paste0("<b>", bad_y$anio, "</b> (", .lox_dec(bad_y$oe, 2), ")", collapse = ", "))
                   else "En todos los años el modelo se mantiene a menos de un 10% de la realidad."))
    add("<p>Muertes <b>esperadas</b> por el modelo de mortalidad ajustado (edad, sexo y enfermedades) frente a las <b>observadas</b>, ",
        "por grupo de edad. Un cociente O/E cercano a 1 indica que el modelo reproduce la mortalidad de ese grupo; ",
        "por construcción el cociente global es 1, así que lo informativo son los grupos de edad y los años.</p>")
    add(.lox_svg_forest(oe_age$grupo_edad, oe_age$oe, oe_age$lo, oe_age$hi))
    add(.lox_html_table(oe_age[, .(`Grupo de edad` = grupo_edad, Observadas = .lox_big(O), Esperadas = .lox_big(E),
                                    `O/E` = .lox_fmt(oe, 2), `IC95%` = sprintf("%s–%s", .lox_fmt(lo, 2), .lox_fmt(hi, 2)))],
                        num = c("Observadas", "Esperadas", "O/E", "IC95%")))
    oe_y <- .lox_oe(obs, "anio")
    if (!is.null(oe_y) && nrow(oe_y) > 2) {
      oe_y <- oe_y[order(anio)]
      oe_y0 <- .lox_oe(obs, "anio", "n_exp0")
      if (!is.null(oe_y0)) {
        oe_y0 <- oe_y0[match(oe_y$anio, anio)]
        add("<p>Por <b>año natural</b>, con el modelo <b>sin</b> término calendario (rojo) y con él (azul). La diferencia entre ambas curvas es lo que absorbe el ajuste por periodo (sección 8).</p>")
        add(.lox_svg_lines(as.character(oe_y$anio), list(list(name = "con calendario", y = oe_y$oe), list(name = "sin calendario", y = oe_y0$oe)),
                           ylab = "observadas / esperadas", ymin = min(0.8, floor(10 * min(c(oe_y$oe, oe_y0$oe), na.rm = TRUE)) / 10), ref = 1))
      } else {
        add("<p>Por <b>año natural</b> (el modelo no incluye tendencia temporal: una desviación sostenida o un pico, como el de 2020, aparece aquí).</p>")
        add(.lox_svg_lines(as.character(oe_y$anio), list(list(name = "O/E", y = oe_y$oe)), ylab = "observadas / esperadas",
                           ymin = min(0.8, floor(10 * min(oe_y$oe)) / 10), ref = 1))
      }
      add("<details><summary>Tabla por año</summary>", .lox_html_table(oe_y[, .(Año = anio, Observadas = .lox_big(O), Esperadas = .lox_big(E),
          `O/E` = .lox_fmt(oe, 2), `IC95%` = sprintf("%s–%s", .lox_fmt(lo, 2), .lox_fmt(hi, 2)))],
          num = c("Observadas", "Esperadas", "O/E", "IC95%")), "</details>")
    }
    oe_s <- .lox_oe(obs, c("Sexo", "grupo_edad"))
    if (!is.null(oe_s)) add("<details><summary>Por sexo y grupo de edad</summary>", .lox_html_table(oe_s[order(Sexo, match(grupo_edad, lv))][,
          .(Sexo, `Grupo de edad` = grupo_edad, Observadas = .lox_big(O), Esperadas = .lox_big(E), `O/E` = .lox_fmt(oe, 2),
            `IC95%` = sprintf("%s–%s", .lox_fmt(lo, 2), .lox_fmt(hi, 2)))], num = c("Observadas", "Esperadas", "O/E", "IC95%")), "</details>")
  }

  # --- efecto calendario
  pe <- rd(r2, "r2_period_effect.csv"); h0 <- rd(r2, "r2_hr_sin_calendario.csv")
  if (!is.null(pe) && nrow(pe)) {
    add("<h2>8. Efecto del calendario en la mortalidad</h2>")
    pd8 <- if (!is.null(co)) co[nodo == "muerte" & termino == "pandemia"] else NULL
    add(.lox_plain("La mortalidad no es igual todos los años: cambian los tratamientos, el sistema sanitario y, en 2020, llegó la pandemia de covid-19. Si el modelo no lo tuviera en cuenta, esos cambios acabarían atribuyéndose, por error, a las enfermedades. Por eso ahora el modelo <b>tiene en cuenta el año</b>.",
                   "La <b>tasa relativa</b> dice cuánto mayor o menor fue la mortalidad de cada año respecto al primero, para personas con la misma edad, sexo y enfermedades (1,0 = igual que el primer año; 0,8 = un 20% menos).",
                   if (!is.null(pd8) && nrow(pd8)) sprintf("El <b>HR pandémico</b> mide cuánto más alta fue la mortalidad en 2020-2021 de lo que predecía la tendencia: aquí fue %s (IC95%% %s–%s).",
                        .lox_vs(exp(pd8$beta), "más alta", "más baja"), .lox_dec(exp(pd8$beta - 1.96 * pd8$se), 2), .lox_dec(exp(pd8$beta + 1.96 * pd8$se), 2)),
                   "<b>Precaución:</b> la línea de tendencia es un ajuste estadístico, no un hallazgo. Al seguir a las mismas personas durante diez años no se puede distinguir si la mortalidad baja de verdad o si simplemente van quedando las personas más robustas. Lo que sí es fiable es el salto de los años pandémicos."))
    add("<p><b>Qué se ha hecho.</b> El modelo de mortalidad incluye, además de edad, sexo y enfermedades, (i) un <b>spline natural del año natural</b> ",
        "con nudos fijos y comunes a todos los sitios y (ii) un <b>indicador de los años pandémicos</b>. Sin ellos, toda la década se resume en un único nivel de riesgo ",
        "y las muertes en exceso de un año (p. ej. 2020) o un descenso sostenido de la mortalidad se reparten entre los coeficientes de las enfermedades. ",
        "El informe muestra el efecto estimado y cuánto cambian los resultados al incluirlo.</p>")
    pe[, `:=`(rr_t = exp(log_rr_tendencia), rr = exp(log_rr_total), lo = exp(log_rr_total - 1.96 * se_total), hi = exp(log_rr_total + 1.96 * se_total))]
    add("<p><b>Tasa de mortalidad relativa por año</b> (a igual edad, sexo y enfermedades; referencia = primer año). Línea azul: tendencia suave; ",
        "naranja: tendencia más el efecto de los años pandémicos.</p>")
    add(.lox_svg_lines(as.character(pe$anio), list(list(name = "tendencia", y = pe$rr_t), list(name = "con pandemia", y = pe$rr)),
                       ylab = "tasa relativa", ymin = min(0.6, floor(10 * min(pe$rr, pe$rr_t)) / 10), ref = 1))
    pd <- if (!is.null(co)) co[nodo == "muerte" & termino == "pandemia"] else NULL
    if (!is.null(pd) && nrow(pd))
      add("<div class='kpis'>", sprintf("<div class='kpi'><b>%s</b><span>HR de los años pandémicos (IC95%% %s–%s), sobre la tendencia</span></div>",
          .lox_fmt(exp(pd$beta), 2), .lox_fmt(exp(pd$beta - 1.96 * pd$se), 2), .lox_fmt(exp(pd$beta + 1.96 * pd$se), 2)), "</div>")
    add("<details><summary>Tabla por año</summary>", .lox_html_table(pe[, .(Año = anio, `Tasa relativa` = .lox_fmt(rr, 2),
        `IC95%` = sprintf("%s–%s", .lox_fmt(lo, 2), .lox_fmt(hi, 2)), `Solo tendencia` = .lox_fmt(rr_t, 2))], num = c("Tasa relativa", "IC95%", "Solo tendencia")), "</details>")
    if (!is.null(h0) && !is.null(co)) {
      m <- merge(h0[, .(nodo, b0 = beta)], co[capa == "mortalidad" & nodo != "muerte", .(nodo, b1 = beta)], by = "nodo")
      if (nrow(m) > 2) {
        m[, `:=`(hr0 = exp(b0), hr1 = exp(b1), d = b1 - b0)]
        add("<p><b>Efecto sobre los HR de las enfermedades.</b> Cada punto es un estado de enfermedad: HR sin término calendario (eje x) frente a HR con él (eje y). ",
            "Si el ajuste por periodo no importa, los puntos caen sobre la diagonal (HR igual en los dos modelos). <b>Por encima</b> de la diagonal, al tener en cuenta el año la enfermedad resulta <i>más</i> peligrosa que sin tenerlo; <b>por debajo</b>, <i>menos</i>. Los puntos más alejados son las enfermedades cuya estimación depende más de cómo se trate el calendario.</p>")
        add(.lox_svg_scatter(m$hr0, m$hr1, .lox_short(m$nodo), hl = rank(-abs(m$d)) <= 5, xlab = "HR sin calendario", ylab = "HR con calendario"))
        add(sprintf("<p class='muted'>Correlación de log-HR %.3f; cambio máximo %.0f%% (%s); mediana del cambio absoluto %.1f%%.</p>",
                    stats::cor(m$b0, m$b1), 100 * (exp(max(abs(m$d))) - 1), .lox_short(m$nodo[which.max(abs(m$d))]), 100 * (exp(stats::median(abs(m$d))) - 1)))
        tb <- m[order(-abs(d))][seq_len(min(8, .N)), .(Estado = .lox_short(nodo), `HR sin calendario` = .lox_fmt(hr0, 2), `HR con calendario` = .lox_fmt(hr1, 2),
                                                       `Cambio %` = .lox_chg(hr1 / hr0))]
        add("<details><summary>Los 8 estados que más cambian</summary>", .lox_html_table(tb, num = names(tb)[-1]), "</details>")
      }
    }
    add("<div class='note warn'><b>Cómo interpretar la tendencia.</b> En una cohorte cerrada (todos entran en el primer año) el año natural y el tiempo transcurrido desde el inicio del seguimiento ",
        "van juntos: la tendencia suave mezcla un posible cambio secular de la mortalidad con la selección de supervivientes y con un registro diagnóstico más o menos completo según el año. ",
        "Se incluye como <i>ajuste</i> y no como hallazgo. Los años de los extremos (primero y último) son los de peor ajuste, con o sin término de calendario (ver sección 7): léanse con cautela. El efecto de los años pandémicos, en cambio, es una desviación brusca respecto a esa tendencia y está mejor identificado.</div>")
  }

  # --- registro diagnostico y pandemia
  dxy <- rd(r2, "r2_dx_by_year.csv")
  pp <- if (!is.null(dxy) && "pandemia" %in% names(dxy)) lox_dx_pandemic(dxy) else NULL
  if (!is.null(pp) && nrow(pp[node == "TOTAL"])) {
    add("<h2>9. Registro diagnóstico durante la pandemia</h2>")
    t9 <- pp[node == "TOTAL"][order(anio)]; d9 <- pp[node != "TOTAL" & anio == min(anio) & pmin(O, E) >= 30][order(ratio)]
    add(.lox_plain("Durante la pandemia se <b>diagnosticaron menos enfermedades de lo habitual</b>: consultas aplazadas, confinamiento, centros de salud saturados. Eso no significa que la gente enfermara menos: muchas personas enfermaron y <i>no fueron diagnosticadas</i> (o lo fueron más tarde). Es lo que se llama <b>infrarregistro</b>.",
                   "Aquí se cuentan los diagnósticos nuevos de cada año y se comparan con los que habría habido si hubiera continuado la tendencia de años anteriores. Un valor de <b>0,84</b> significa un 16% menos de diagnósticos de los esperados; <b>1,00</b>, los esperados.",
                   sprintf("En este sitio: %s.", paste0("<b>", t9$anio, "</b>: ", vapply(t9$ratio, .lox_vs_dx, ""), collapse = "; ")),
                   if (nrow(d9) >= 3) sprintf("Las enfermedades con mayor déficit en %d fueron %s; las que menos cambiaron, %s. Suele ser mayor en enfermedades que se detectan en consultas o revisiones programadas, y menor en las de presentación urgente.",
                        min(pp$anio), paste0("<b>", .lox_short(d9$node[1:3]), "</b> (", .lox_dec(d9$ratio[1:3], 2), ")", collapse = ", "),
                        paste0(.lox_short(rev(d9$node)[1:3]), " (", .lox_dec(rev(d9$ratio)[1:3], 2), ")", collapse = ", "))))
    add("<p>Cuántos diagnósticos <b>nuevos</b> se registraron cada año, frente a los que cabía esperar según la <b>tendencia de los años previos</b> a la pandemia ",
        "(Poisson con tendencia lineal del año, ajustado solo con los años anteriores y extrapolado). Un déficit de diagnósticos es el <b>infrarregistro</b> descrito en la literatura ",
        "durante la pandemia (confinamiento, saturación asistencial, menor acceso): personas que enferman pero no figuran como enfermas hasta más tarde, o nunca.</p>")
    tt <- dxy[node == "TOTAL"][order(anio)]; tt[, tasa := 1000 * n_inc / n_risk]
    pt <- pp[node == "TOTAL"]; ex <- rep(NA_real_, nrow(tt)); ex[match(pt$anio, tt$anio)] <- 1000 * pt$E / tt$n_risk[match(pt$anio, tt$anio)]
    add("<p><b>Todas las enfermedades juntas</b>: diagnósticos incidentes por 1.000 personas en riesgo (sin la enfermedad al 1 de enero). La tasa baja con los años por el envejecimiento ",
        "de la cohorte cerrada; lo informativo es la distancia a la línea esperada en los años pandémicos.</p>")
    add(.lox_svg_lines(as.character(tt$anio), list(list(name = "observada", y = tt$tasa), list(name = "esperada", y = ex)), ylab = "diagnósticos / 1.000 en riesgo",
                       ymin = max(0, floor(0.8 * min(tt$tasa, na.rm = TRUE)))))
    add("<div class='kpis'>", paste0(sprintf("<div class='kpi'><b>%s</b><span>%d: observados / esperados (IC95%% %s–%s)</span></div>",
        .lox_fmt(pt$ratio, 2), pt$anio, .lox_fmt(pt$lo, 2), .lox_fmt(pt$hi, 2)), collapse = ""), "</div>")
    y0 <- min(pp$anio); f <- pp[node != "TOTAL" & anio == y0 & pmin(O, E) >= 30][order(ratio)]
    if (nrow(f) > 2) {
      add("<p><b>Por enfermedad, ", y0, "</b>: cociente de diagnósticos observados / esperados. Por debajo de 1 = menos diagnósticos de los esperados. Solo enfermedades con al menos 30 diagnósticos observados y esperados.</p>")
      add(.lox_svg_forest(.lox_short(f$node), f$ratio, f$lo, f$hi))
      w <- dcast(pp[node != "TOTAL"], node ~ anio, value.var = "ratio")
      tb <- w[order(w[[2]])][, c(list(Enfermedad = .lox_short(node)), lapply(.SD, .lox_fmt, d = 2)), .SDcols = setdiff(names(w), "node")]
      add("<details><summary>Tabla por enfermedad y año pandémico (todas)</summary>", .lox_html_table(tb, num = setdiff(names(tb), "Enfermedad")), "</details>")
    }
    add("<div class='note warn'><b>Relación con la mortalidad.</b> Si los diagnósticos se registraron menos en estos años, parte de quienes murieron figuraban como <i>sin</i> la enfermedad, ",
        "y el indicador pandémico del modelo de mortalidad (sección 8) absorbe parte de ese efecto. Es un patrón de <b>registro</b>, no necesariamente de incidencia real. ",
        "La esperada es una extrapolación lineal y es sensible al número de años previos disponibles.</div>")
  }

  add("<h2>Glosario</h2><details open><summary>Términos usados en este informe</summary><dl class='gl'>",
      "<dt>Cohorte</dt><dd>Grupo de personas a las que se sigue en el tiempo. Aquí es <i>cerrada</i>: todas entran al inicio y no se incorporan nuevas.</dd>",
      "<dt>Prevalencia</dt><dd>Porcentaje de personas que tienen una enfermedad en un momento dado.</dd>",
      "<dt>Incidencia</dt><dd>Aparición de casos <i>nuevos</i> de una enfermedad en un periodo.</dd>",
      "<dt>Persona-año</dt><dd>Una persona observada durante un año; es la unidad con la que se calculan las tasas.</dd>",
      "<dt>HR (cociente de tasas)</dt><dd>Cuántas veces mayor es la tasa de muerte en un grupo que en otro. 1 = igual; 2 = el doble.</dd>",
      "<dt>OR (odds ratio)</dt><dd>Cuánto cambian las probabilidades de que ocurra algo. Se parece al aumento del riesgo cuando el suceso es poco frecuente.</dd>",
      "<dt>Intervalo de confianza al 95% (IC95%)</dt><dd>Rango de valores compatibles con los datos. Si incluye el 1 (en HR u OR), no se puede afirmar una diferencia.</dd>",
      "<dt>O/E (observadas / esperadas)</dt><dd>Muertes reales frente a las que predice el modelo. 1 = el modelo acierta.</dd>",
      "<dt>Spline</dt><dd>Curva flexible que permite que un efecto (la edad, el año) cambie de forma suave en lugar de ser una recta.</dd>",
      "<dt>Arco / flecha</dt><dd>Relación estable en la que una enfermedad suele aparecer antes que otra. Indica orden temporal, no causalidad demostrada.</dd>",
      "<dt>Supresión de celdas</dt><dd>Para proteger la privacidad, los recuentos de menos de 10 personas no se muestran.</dd>",
      "<dt>Infrarregistro</dt><dd>Enfermedades que existen pero no quedan registradas (o se registran más tarde) en la historia clínica.</dd>",
      "</dl></details>")

  add("<footer>Informe generado a partir de tablas agregadas con supresión de celdas &lt; ", .lox_esc(mf$min_cell),
      ". No contiene datos individuales. R ", .lox_esc(sub("R version ", "", mf$r_version)), ".</footer>")
  html <- paste0("<!doctype html><html lang='es'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
                 "<title>Informe local ", .lox_esc(sid), "</title><style>", .lox_css, "</style></head><body><main>",
                 paste(H, collapse = "\n"), "</main></body></html>")
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_file, useBytes = TRUE)
  .lox_msg("Informe local -> %s", out_file)
  invisible(out_file)
}
