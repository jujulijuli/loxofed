# Orquestadores: Ronda 1 y 2 (en cada sitio) y pasos centrales ---------------------------------------

.lox_cfg <- function(cfg) if (is.character(cfg)) lox_config(cfg) else cfg

.lox_prepare <- function(cfg, dir, env) {
  cdm <- lox_read_cdm(cfg)
  qc <- lox_qc(cdm, cfg)
  .lox_write(qc, "qc.csv", dir, env)
  if (any(qc$flag == "ERROR")) {
    print(qc[flag != "OK"]); stop("El control de calidad encontro ERRORES (ver qc.csv). Corrigelos antes de seguir.", call. = FALSE)
  }
  if (any(qc$flag == "WARN")) { .lox_msg("QC con avisos:"); print(qc[flag == "WARN"]) }
  cohort <- lox_build_cohort(cdm, cfg)
  if (!nrow(cohort$cohorte)) stop("La cohorte esta vacia con los criterios actuales", call. = FALSE)
  if (cohort$n_eligible > nrow(cohort$cohorte))
    .lox_msg("Cohorte: %d personas (muestra aleatoria de %d elegibles, sample_n = %d)",
             nrow(cohort$cohorte), cohort$n_eligible, cfg$cohort$sample_n)
  else .lox_msg("Cohorte: %d personas", nrow(cohort$cohorte))
  list(cdm_resolution = qc[check == "resolucion_dx_pct_dia/mes/anio", value], cohort = cohort,
       arr = lox_build_states(cohort, cfg))
}

#' Ronda 1 (en el sitio): QC, cohorte, estructura, prevalencia y precedencia por fechas
#'
#' Escribe en `<output_dir>/round1/` solo tablas agregadas y un `manifest.json`.
#' @param cfg lista de configuracion o ruta al YAML
#' @export
lox_run_round1 <- function(cfg) {
  cfg <- .lox_cfg(cfg); mc <- cfg$disclosure$min_cell
  dir <- file.path(cfg$paths$output_dir, "round1"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  env <- .lox_manifest_new()
  p <- .lox_prepare(cfg, dir, env); arr <- p$arr; cohort <- p$cohort

  prev <- lox_prevalence(arr)
  nodes <- prev[prev >= cfg$model$prev_min, node]
  .lox_msg("Nodos locales (prev >= %.1f%%): %d", 100 * cfg$model$prev_min, length(nodes))
  dcnt <- lox_direction_counts(cohort, nodes)
  d2 <- lox_build_transitions(arr, cohort, nodes, cfg); rm(arr); invisible(gc())
  .lox_msg("Transiciones: %d | aprendiendo estructura (R=%d)", nrow(d2), cfg$structure$boot_R)
  st <- lox_learn_structure(d2, nodes, cfg)

  w <- cohort$sample_weight       # N elegibles / n muestreados (1 si no hay tope); solo para recuentos descriptivos
  pv <- lox_suppress(prev, c("n_pos"), mc)
  pv[, `:=`(n_pos_est = round(n_pos * w), n_est = round(n * w))]   # tras suprimir: no revela celdas pequenas
  .lox_write(pv, "r1_prevalence.csv", dir, env)
  .lox_write(data.table(node = nodes), "r1_nodes.csv", dir, env)
  .lox_write(st$edges, "r1_edges.csv", dir, env)
  .lox_write(lox_suppress(dcnt, c("n_ab", "n_ba", "n_tie"), mc, total = "n_both"), "r1_direction.csv", dir, env)
  meta <- data.table(site_id = cfg$site_id, n_eligible = cohort$n_eligible, n_cohort = nrow(cohort$cohorte),
                     sample_weight = round(w, 4), n_transitions = nrow(d2),
                     n_fit = st$n_fit, boot_R = cfg$structure$boot_R, dx_resolution_pct = p$cdm_resolution,
                     bnlearn_version = as.character(utils::packageVersion("bnlearn")))
  .lox_write(lox_suppress(meta, c("n_eligible", "n_cohort"), mc), "r1_meta.csv", dir, env)
  .lox_manifest_write(env, dir, "round1", cfg)
  .lox_msg("Ronda 1 completa -> %s (comparte esta carpeta con el nodo central)", dir)
  invisible(dir)
}

#' Paso central A: construye la estructura consenso a partir de las carpetas round1 de cada sitio
#' @param site_dirs vector de rutas a carpetas `round1` (una por sitio)
#' @param out_dir carpeta de salida (se crea `consensus/`)
#' @param cfg configuracion comun (lista o YAML); define umbrales y huella
#' @export
lox_run_central_consensus <- function(site_dirs, out_dir, cfg) {
  cfg <- .lox_cfg(cfg); fp <- lox_config_fingerprint(cfg)
  mf <- lapply(site_dirs, .lox_manifest_read)
  for (i in seq_along(mf)) {
    if (mf[[i]]$round != "round1") stop("No es una carpeta round1: ", site_dirs[i], call. = FALSE)
    if (!identical(mf[[i]]$fingerprint, fp))
      stop("La configuracion del sitio ", mf[[i]]$site_id, " no coincide con la comun (huella distinta): ",
           "cohorte/modelo/estructura deben ser identicos en todos los sitios", call. = FALSE)
  }
  vv <- unique(vapply(mf, function(m) m$package_version, ""))
  if (length(vv) > 1) stop("Versiones del paquete distintas entre sitios: ", paste(vv, collapse = ", "), call. = FALSE)
  ids <- vapply(mf, function(m) m$site_id, ""); if (anyDuplicated(ids)) stop("site_id repetido", call. = FALSE)
  rd <- function(f) setNames(lapply(site_dirs, function(d) fread(file.path(d, f))), ids)
  meta <- rd("r1_meta.csv")
  w <- setNames(vapply(meta, function(m) as.numeric(m$n_fit), numeric(1)), ids)
  res <- lox_consensus(prev = rd("r1_prevalence.csv"), edges = rd("r1_edges.csv"),
                       nodes = rd("r1_nodes.csv"), dirs = rd("r1_direction.csv"), w = as.list(w), cfg = cfg)
  dir <- file.path(out_dir, "consensus"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  env <- .lox_manifest_new()
  .lox_write(res$nodes, "consensus_nodes.csv", dir, env)
  .lox_write(res$structure, "consensus_structure.csv", dir, env)
  .lox_write(res$pairs, "consensus_pairs_all.csv", dir, env)
  .lox_manifest_write(env, dir, "consensus", cfg, extra = list(sites = as.list(ids)))
  # Descriptivos por sitio para el informe central. Van en `descriptive/` (NO en consensus/): solo `consensus/`
  # se devuelve a los sitios; asi cada sitio no recibe las tablas de los demas.
  dd <- file.path(out_dir, "descriptive"); dir.create(dd, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(rbindlist(lapply(ids, function(s) rd("r1_prevalence.csv")[[s]][, .(site = s, node, n_pos, n, prev)])), file.path(dd, "prevalence_by_site.csv"))
  data.table::fwrite(rbindlist(lapply(ids, function(s) rd("r1_edges.csv")[[s]][, .(site = s, from, to, strength)])), file.path(dd, "edges_by_site.csv"))
  n_c <- if (nrow(res$structure)) res$structure[state %in% c("concordante", "reorientada"), .N] else 0L
  .lox_msg("Consenso: %d nodos, %d arcos causales (%d sitios) -> %s", nrow(res$nodes), n_c, length(ids), dir)
  invisible(res)
}

#' Ronda 2 (en el sitio): ajusta capa 1 y capa 2 sobre la estructura consenso
#' @param consensus_dir carpeta `consensus/` recibida del nodo central
#' @export
lox_run_round2 <- function(cfg, consensus_dir) {
  cfg <- .lox_cfg(cfg); mc <- cfg$disclosure$min_cell
  cm <- .lox_manifest_read(consensus_dir)
  if (!identical(cm$fingerprint, lox_config_fingerprint(cfg)))
    stop("La configuracion local no coincide con la que genero el consenso (huella distinta).", call. = FALSE)
  nodes <- fread(file.path(consensus_dir, "consensus_nodes.csv"))$node
  st <- fread(file.path(consensus_dir, "consensus_structure.csv"))
  edges <- if (nrow(st)) st[state %in% c("concordante", "reorientada"), .(from, to)] else data.table(from = character(), to = character())
  dir <- file.path(cfg$paths$output_dir, "round2"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  env <- .lox_manifest_new()
  p <- .lox_prepare(cfg, dir, env)
  n_elig <- p$cohort$n_eligible; n_coh <- nrow(p$cohort$cohorte); w <- p$cohort$sample_weight
  d2 <- lox_build_transitions(p$arr, p$cohort, nodes, cfg); p$arr <- NULL; invisible(gc())
  dxy <- lox_dx_by_year(p$cohort, nodes, cfg)                    # diagnosticos incidentes por anio (agregado)
  dxs <- lox_dx_strata(p$cohort, nodes, cfg, mc)                  # incidencia por anio x sexo x grupo de edad (agregado)
  arct <- lox_arc_tables(d2, edges, mc)                          # tablas 2x2 por arco del consenso (agregado)
  .lox_msg("Ajustando capa 1 (%d nodos)", length(nodes))
  l1 <- lox_fit_layer1(d2, nodes, edges, cfg)
  if (identical(cfg$model$mortality, "exact")) {
    pt <- lox_build_persontime(p$cohort, nodes, cfg)
    .lox_msg("Capa 2 exacta: %d tramos persona, %d muertes, %.0f personas-anio",
             nrow(pt), sum(pt$event), sum(pt$exposure))
    l2 <- lox_fit_layer2_exact(pt, nodes, cfg); obs <- .lox_observed_exact(pt, l2$expected, l2$expected0); rm(pt)
    l2$expected <- NULL; l2$expected0 <- NULL
    obs <- lox_suppress(obs, "n_death", mc)
    for (cl in intersect(c("n_exp", "n_exp0"), names(obs))) obs[suppressed == TRUE, (cl) := NA_real_]   # la esperada delataria la celda suprimida
  } else {
    .lox_msg("Capa 2 anual (convencion t -> t+1 del informe de Navarra)")
    l2 <- lox_fit_layer2(d2, nodes, cfg)
    obs <- lox_suppress(.lox_observed(d2, cfg), c("n_risk", "n_death"), mc)
  }
  rm(p); invisible(gc())
  if ("person_years" %in% names(obs)) obs[, person_years_est := round(person_years * w, 1)]
  if ("n_risk" %in% names(obs)) obs[, n_risk_est := round(n_risk * w)]
  if ("n_exp" %in% names(obs)) obs[, n_exp_est := round(n_exp * w, 1)]
  if ("n_exp0" %in% names(obs)) obs[, n_exp0_est := round(n_exp0 * w, 1)]
  obs[, n_death_est := round(n_death * w)]                      # tras suprimir (NA se mantiene)
  co <- rbindlist(list(l1$coefs, l2$coefs), fill = TRUE); co[, site := cfg$site_id]
  setcolorder(co, c("site", "capa", "nodo", "termino", "beta", "se"))
  .lox_write(co, "r2_coefs.csv", dir, env)
  .lox_write(l1$skipped, "r2_layer1_skipped.csv", dir, env)
  .lox_write(data.table(estado_excluido_hazard = as.character(l2$excluded)), "r2_layer2_excluded.csv", dir, env)
  .lox_write(obs, "r2_observed.csv", dir, env)
  dxy[, pandemia := anio %in% as.integer(cfg$model$calendar$pandemic_years)]
  dxy <- lox_suppress(dxy, "n_inc", mc); dxy[, n_inc_est := round(n_inc * w)]; dxy[, n_risk_est := round(n_risk * w)]
  .lox_write(dxy, "r2_dx_by_year.csv", dir, env)
  .lox_write(dxs, "r2_incidence_strata.csv", dir, env)
  .lox_write(arct, "r2_arc_tables.csv", dir, env)
  if (!is.null(l2$by_sex) && nrow(l2$by_sex)) .lox_write(l2$by_sex, "r2_hr_por_sexo.csv", dir, env)
  if (!is.null(l2$period)) {                                    # efecto calendario: curva y HR sin ajuste (para el informe)
    .lox_write(l2$period, "r2_period_effect.csv", dir, env)
    .lox_write(l2$coefs0, "r2_hr_sin_calendario.csv", dir, env)
  }
  .lox_manifest_write(env, dir, "round2", cfg, extra = list(layer2_converged = l2$converged,
                                                           mortality = cfg$model$mortality,
                                                           n_eligible = n_elig, n_cohort = n_coh,
                                                           sample_weight = w))
  .lox_msg("Ronda 2 completa -> %s", dir)
  if (isTRUE(cfg$report$auto)) tryCatch(lox_report_local(cfg$paths$output_dir, consensus_dir = consensus_dir),
                                        error = function(e) .lox_msg("Informe local no generado: %s", conditionMessage(e)))
  invisible(dir)
}

#' Paso central B: meta-analisis de los coeficientes de la Ronda 2
#' @param site_dirs rutas a carpetas `round2`
#' @export
lox_run_central_meta <- function(site_dirs, out_dir, cfg) {
  cfg <- .lox_cfg(cfg); fp <- lox_config_fingerprint(cfg)
  mf <- lapply(site_dirs, .lox_manifest_read)
  for (m in mf) if (!identical(m$fingerprint, fp)) stop("Huella distinta en el sitio ", m$site_id, call. = FALSE)
  co <- rbindlist(lapply(site_dirs, function(d) fread(file.path(d, "r2_coefs.csv"))))
  mt <- lox_meta(co, min_k = cfg$meta$min_sites)
  dir <- file.path(out_dir, "meta"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  env <- .lox_manifest_new()
  .lox_write(mt, "meta_coefs.csv", dir, env)
  .lox_write(mt[capa == "mortalidad" & nodo != "muerte", .(estado = nodo, k, status, HR = exp_re, lo, hi, pi_lo, pi_hi, p, tau2, I2)][order(-HR)],
             "meta_hazard_HR.csv", dir, env)
  .lox_write(co, "meta_site_coefs.csv", dir, env)            # coeficientes por sitio (para los diagramas de bosque)
  ob <- rbindlist(lapply(seq_along(site_dirs), function(i) { f <- file.path(site_dirs[i], "r2_observed.csv")
    if (file.exists(f)) fread(f)[, site := mf[[i]]$site_id] }), fill = TRUE)
  if (nrow(ob)) .lox_write(ob, "meta_observed.csv", dir, env)  # mortalidad observada por sitio (celdas ya suprimidas)
  pe <- rbindlist(lapply(seq_along(site_dirs), function(i) { f <- file.path(site_dirs[i], "r2_period_effect.csv")
    if (file.exists(f)) fread(f)[, site := mf[[i]]$site_id] }), fill = TRUE)
  if (nrow(pe)) .lox_write(pe, "meta_period.csv", dir, env)   # efecto calendario por sitio
  dx <- rbindlist(lapply(seq_along(site_dirs), function(i) { f <- file.path(site_dirs[i], "r2_dx_by_year.csv")
    if (file.exists(f)) fread(f)[, site := mf[[i]]$site_id] }), fill = TRUE)
  if (nrow(dx)) .lox_write(dx, "meta_dx_by_year.csv", dir, env)   # diagnosticos incidentes por anio y sitio
  rdsite <- function(f) rbindlist(lapply(seq_along(site_dirs), function(i) { ff <- file.path(site_dirs[i], f)
    if (file.exists(ff)) fread(ff)[, site := mf[[i]]$site_id] }), fill = TRUE)
  bsx <- rdsite("r2_hr_por_sexo.csv")
  if (nrow(bsx)) {                                              # HR por sexo y contraste hombre - mujer (estratos independientes)
    .lox_write(bsx, "meta_site_hr_sexo.csv", dir, env)
    ms <- lox_meta(bsx[, .(site, capa = paste0("mortalidad_", sexo), nodo, termino = "estado", beta, se)], min_k = cfg$meta$min_sites)
    wd <- dcast(bsx, site + nodo ~ sexo, value.var = c("beta", "se"))
    if (all(c("beta_Hombre", "beta_Mujer") %in% names(wd))) {
      wd <- wd[is.finite(beta_Hombre) & is.finite(beta_Mujer)]
      wd[, `:=`(capa = "mortalidad_dif_sexo", termino = "hombre_menos_mujer", beta = beta_Hombre - beta_Mujer, se = sqrt(se_Hombre^2 + se_Mujer^2))]
      ms <- rbind(ms, lox_meta(wd[, .(site, capa, nodo, termino, beta, se)], min_k = cfg$meta$min_sites), fill = TRUE)
    }
    .lox_write(ms, "meta_sex.csv", dir, env)
  }
  ist <- rdsite("r2_incidence_strata.csv"); if (nrow(ist)) .lox_write(ist, "meta_incidence_strata.csv", dir, env)
  at <- rdsite("r2_arc_tables.csv"); if (nrow(at)) .lox_write(at, "meta_arc_tables.csv", dir, env)
  g <- function(m, k) if (is.null(m[[k]])) NA else m[[k]]
  si <- rbindlist(lapply(mf, function(m) data.table(site = m$site_id, version = m$package_version, mortality = g(m, "mortality"),
        n_eligible = g(m, "n_eligible"), n_cohort = g(m, "n_cohort"), sample_weight = g(m, "sample_weight"))), fill = TRUE)
  si <- merge(si, co[, .(n_coefs = .N, converged = all(converged)), by = site], by = "site", all.x = TRUE)
  .lox_write(si, "meta_sites.csv", dir, env)
  .lox_manifest_write(env, dir, "meta", cfg, extra = list(sites = lapply(mf, function(m) m$site_id)))
  .lox_msg("Meta-analisis: %d coeficientes, %d sitios -> %s", nrow(mt), length(site_dirs), dir)
  if (isTRUE(cfg$report$auto) && file.exists(file.path(out_dir, "consensus", "consensus_structure.csv")))
    tryCatch(lox_report_central(out_dir), error = function(e) .lox_msg("Informe central no generado: %s", conditionMessage(e)))
  invisible(mt)
}
