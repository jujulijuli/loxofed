# Lectura del modelo comun de datos (CDM) y control de calidad ------------------------------------

#' Lee las dos tablas del CDM (personas y enfermedades) y las normaliza
#'
#' Devuelve `pac` (id, birth_year, sex, death_date, baja_date) y `enf` (id, code, dx_date, dx_res).
#' Las fechas se parsean con resolucion variable (dia, mes o anio); ver `.lox_parse_date`.
#' @export
lox_read_cdm <- function(cfg) {
  .lox_need(cfg$paths$pacientes, "paths$pacientes"); .lox_need(cfg$paths$enfermedades, "paths$enfermedades")
  cc <- cfg$cdm$cols
  rd <- function(path, cols, what) {
    if (!file.exists(path)) stop("No existe ", what, ": ", path, call. = FALSE)
    hdr  <- names(data.table::fread(path, sep = cfg$cdm$sep, nrows = 0, colClasses = "character"))
    miss <- setdiff(unname(unlist(cols)), hdr)
    if (length(miss)) stop(sprintf("%s: faltan columnas %s (disponibles: %s)%s", what,
                                   paste(miss, collapse = ", "), paste(hdr, collapse = ", "),
                                   if (!is.null(cc$baja) && cc$baja %in% miss) sprintf(". Si el sitio no tiene la fecha de baja, pon `baja: ~` en cdm$cols del config; si se llama distinto, indica su nombre.", cc$baja) else ""), call. = FALSE)
    data.table::fread(path, sep = cfg$cdm$sep, select = unname(unlist(cols)), colClasses = "character",
                      na.strings = c("", "NA"), showProgress = FALSE)
  }
  pc <- Filter(Negate(is.null), cc[c("id", "birth", "sex", "death", "baja")])
  ec <- cc[c("id", "dx_code", "dx_date")]
  p <- rd(cfg$paths$pacientes, pc, "tabla de pacientes")
  e <- rd(cfg$paths$enfermedades, ec, "tabla de enfermedades")

  pd <- .lox_parse_date(p[[cc$death]])
  pac <- data.table(id = as.character(p[[cc$id]]), birth_year = .lox_year(p[[cc$birth]]),
                    sex_raw = as.character(p[[cc$sex]]),
                    death_date = pd$date, death_res = pd$res,
                    baja_date = if (!is.null(cc$baja)) .lox_parse_date(p[[cc$baja]])$date else as.Date(NA))
  sm <- cfg$cdm$sex_map
  pac[, sex := fifelse(sex_raw %in% sm$female, "Mujer", fifelse(sex_raw %in% sm$male, "Hombre", NA_character_))]

  ed <- .lox_parse_date(e[[cc$dx_date]])
  enf <- data.table(id = as.character(e[[cc$id]]), code = as.character(e[[cc$dx_code]]),
                    dx_date = ed$date, dx_res = ed$res)
  list(pac = pac, enf = enf)
}

#' Control de calidad de las entradas (solo agregados; seguro de compartir)
#'
#' Devuelve una tabla `check / value / flag`. `flag` es "OK", "WARN" o "ERROR". Los conteos
#' pequenos se suprimen con `cfg$disclosure$min_cell`.
#' @export
lox_qc <- function(cdm, cfg) {
  pac <- cdm$pac; enf <- cdm$enf; mc <- cfg$disclosure$min_cell
  out <- list()
  add <- function(check, value, flag = "OK") out[[length(out) + 1L]] <<-
    data.table(check = check, value = as.character(value), flag = flag)
  sup <- function(n) if (is.na(n)) NA_character_ else if (n > 0 && n < mc) paste0("<", mc) else as.character(n)

  add("n_pacientes", nrow(pac))
  nd <- sum(duplicated(pac$id)); add("ids_duplicados_pacientes", sup(nd), if (nd > 0) "ERROR" else "OK")
  nb <- sum(is.na(pac$birth_year)); add("sin_anio_nacimiento", sup(nb), if (nb > 0) "WARN" else "OK")
  br <- range(pac$birth_year, na.rm = TRUE)
  add("rango_anio_nacimiento", paste(br, collapse = "-"), if (br[1] < 1880 || br[2] > 2025) "WARN" else "OK")
  ns <- sum(is.na(pac$sex)); add("sexo_sin_mapear", sup(ns), if (ns > 0) "ERROR" else "OK")
  add("n_exitus", sum(!is.na(pac$death_date)))
  dr <- table(factor(pac$death_res, levels = c("day", "month", "year")))
  thr_d <- if (identical(cfg$model$mortality, "annual_next")) 0.2 else 0.05   # en modo exacto la muerte necesita mes/dia
  add("resolucion_exitus_pct_dia/mes/anio", paste(round(100 * prop.table(dr), 1), collapse = "/"),
      if (sum(dr) > 0 && dr[["year"]] / sum(dr) > thr_d) "WARN" else "OK")
  dbb <- sum(!is.na(pac$death_date) & !is.na(pac$birth_year) &
             as.integer(format(pac$death_date, "%Y")) < pac$birth_year)
  add("exitus_antes_de_nacer", sup(dbb), if (dbb > 0) "ERROR" else "OK")

  add("n_filas_enfermedades", nrow(enf))
  unk <- enf[!code %in% LOXO_DIS60, .N]
  add("filas_con_codigo_fuera_de_catalogo", sup(unk), if (unk > 0) "WARN" else "OK")
  add("codigos_del_catalogo_presentes", uniqueN(enf[code %in% LOXO_DIS60, code]),
      if (uniqueN(enf[code %in% LOXO_DIS60, code]) < 50) "WARN" else "OK")
  nf <- enf[is.na(dx_date), .N]; add("filas_con_fecha_no_parseable", sup(nf), if (nf > 0) "WARN" else "OK")
  rr <- table(factor(enf$dx_res, levels = c("day", "month", "year")))
  fy <- if (sum(rr) > 0) rr[["year"]] / sum(rr) else NA
  add("resolucion_dx_pct_dia/mes/anio", paste(round(100 * prop.table(rr), 1), collapse = "/"),
      if (!is.na(fy) && fy > 0.2) "WARN" else "OK")
  add("dx_huerfanos_sin_paciente", sup(enf[!id %in% pac$id, .N]),
      if (enf[!id %in% pac$id, .N] > 0) "WARN" else "OK")
  m <- merge(enf[!is.na(dx_date)], pac[!is.na(death_date), .(id, death_date)], by = "id")
  da <- m[dx_date > death_date, .N]; add("dx_posterior_al_exitus", sup(da), if (da > 0) "WARN" else "OK")
  # --- fecha de baja frente a fecha de muerte (regla del analisis: la baja CENSURA; la muerte solo cuenta si ocurre en la fecha de baja o antes)
  if (any(!is.na(pac$baja_date))) {
    yf <- cfg$cohort$year_first; yl <- cfg$cohort$year_last
    d0 <- as.Date(sprintf("%d-01-01", yf)); d1 <- as.Date(sprintf("%d-01-01", yl + 1L))
    nbj <- sum(!is.na(pac$baja_date)); add("n_baja", sup(nbj))
    both <- pac[!is.na(death_date) & !is.na(baja_date)]
    add("baja_y_exitus", sup(nrow(both)))
    add("exitus_con_baja_anterior", sup(both[death_date > baja_date, .N]))          # la muerte se pierde (censura en la baja)
    add("exitus_con_baja_igual", sup(both[death_date == baja_date, .N]))
    add("exitus_con_baja_posterior", sup(both[death_date < baja_date, .N]))         # la muerte cuenta; la baja es posterior
    lost <- both[death_date > baja_date & death_date >= d0 & death_date < d1]
    nde <- pac[!is.na(death_date) & death_date >= d0 & death_date < d1, .N]
    add("exitus_perdidos_por_baja_anterior", sup(nrow(lost)), if (nde > 0 && nrow(lost) / nde > 0.01) "WARN" else "OK")
    add("exitus_perdidos_pct_de_exitus_seguimiento", if (nde > 0) sprintf("%.1f", 100 * nrow(lost) / nde) else NA_character_)
    add("exitus_perdidos_dentro_de_30_dias_de_la_baja", sup(lost[as.numeric(death_date - baja_date) <= 30, .N]))
    mb <- merge(enf[!is.na(dx_date)], pac[!is.na(baja_date), .(id, baja_date)], by = "id")
    add("dx_posterior_a_baja", sup(mb[dx_date > baja_date, .N]))                    # tras la baja no hay seguimiento
  }
  yr <- as.integer(format(range(enf$dx_date, na.rm = TRUE), "%Y"))
  add("rango_anio_dx", paste(yr, collapse = "-"), if (yr[1] < 1900 || yr[2] > 2030) "WARN" else "OK")
  rbindlist(out)
}
