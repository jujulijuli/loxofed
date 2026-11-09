# Control de divulgacion estadistica -----------------------------------------------------------------

#' Suprime celdas pequenas (0 < n < min_cell) en las columnas de conteo
#' @param dt data.table
#' @param cols columnas de conteo
#' @param min_cell umbral
#' @param total columna que es la SUMA de `cols` (p.ej. `n_both`). Si en una fila se suprime algun
#'   componente, tambien se anula el total (supresion complementaria): si no, el componente
#'   suprimido se despeja por diferencia.
#' @return dt con las celdas pequenas a NA y una columna logica `suppressed`
#' @export
lox_suppress <- function(dt, cols, min_cell, total = NULL) {
  dt <- copy(as.data.table(dt)); dt[, suppressed := FALSE]
  for (cl in cols) {
    small <- !is.na(dt[[cl]]) & dt[[cl]] > 0 & dt[[cl]] < min_cell
    if (any(small)) { dt[small, suppressed := TRUE]; set(dt, which(small), cl, NA_integer_) }
  }
  if (!is.null(total) && any(dt$suppressed)) set(dt, which(dt$suppressed), total, NA_integer_)
  dt[]
}

# Escribe una tabla de salida con salvaguardas: sin identificadores, tamano razonable
.lox_write <- function(dt, name, dir, env) {
  bad <- intersect(tolower(names(dt)), c("id", "id_pac", "idpac", "patient_id", "paciente"))
  if (length(bad)) stop("Salida bloqueada: la tabla ", name, " contiene identificadores: ", bad, call. = FALSE)
  if (nrow(dt) > 2e6) stop("Salida bloqueada: ", name, " tiene demasiadas filas (", nrow(dt),
                           "); parece microdato", call. = FALSE)
  f <- file.path(dir, name); data.table::fwrite(dt, f)
  env$files[[name]] <- list(rows = nrow(dt), cols = ncol(dt), md5 = unname(tools::md5sum(f)))
  invisible(f)
}

.lox_manifest_new <- function() { e <- new.env(); e$files <- list(); e }

.lox_manifest_write <- function(env, dir, round, cfg, extra = list()) {
  m <- c(list(package = "loxofed", package_version = lox_version(), r_version = R.version.string,
              site_id = cfg$site_id, round = round, created = format(Sys.time(), tz = "UTC", usetz = TRUE),
              fingerprint = lox_config_fingerprint(cfg), min_cell = cfg$disclosure$min_cell,
              files = env$files), extra)
  jsonlite::write_json(m, file.path(dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(m)
}

.lox_manifest_read <- function(dir) {
  f <- file.path(dir, "manifest.json")
  if (!file.exists(f)) stop("No hay manifest.json en ", dir, call. = FALSE)
  jsonlite::read_json(f, simplifyVector = FALSE)
}
