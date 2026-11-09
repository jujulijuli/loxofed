.lox_usage <- "loxofed %s

Uso:  loxofed <orden> [opciones]

En cada sitio:
  qc                 --config FICHERO.yml                    Control de calidad de las entradas
  round1             --config FICHERO.yml                    Estructura + agregados (-> output/round1)
  round2             --config FICHERO.yml --consensus DIR    Ajuste capa 1 y 2 (-> output/round2)
  report             --config FICHERO.yml [--consensus DIR]  Informe local HTML (-> output/informe_local.html)
  simulate           --out DIR [--n 5000] [--seed 1] [--date-res day|month|year]

En el nodo central:
  central-consensus  --config FICHERO.yml --sites DIR1,DIR2,... --out DIR
  central-meta       --config FICHERO.yml --sites DIR1,DIR2,... --out DIR   (genera informe_central.html)
  central-report     --out DIR                                Regenera informe_central.html
  compare            --a DIR --b DIR [--labels A,B] [--out F]  Compara dos ejecuciones de un sitio (anterior/actual)

  version            Version del paquete
"

.lox_parse_args <- function(a) {
  out <- list(); i <- 1
  while (i <= length(a)) {
    if (!startsWith(a[i], "--")) stop("Argumento inesperado: ", a[i], call. = FALSE)
    k <- sub("^--", "", a[i]); if (i == length(a)) stop("Falta el valor de --", k, call. = FALSE)
    out[[k]] <- a[i + 1]; i <- i + 2
  }
  out
}

#' Interfaz de linea de comandos (la usa el contenedor Docker)
#' @export
lox_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || args[1] %in% c("-h", "--help", "help")) { cat(sprintf(.lox_usage, lox_version())); return(invisible()) }
  cmd <- args[1]; o <- .lox_parse_args(args[-1])
  need <- function(k) if (is.null(o[[k]])) stop("Falta --", k, call. = FALSE) else o[[k]]
  sites <- function() strsplit(need("sites"), ",")[[1]]
  switch(cmd,
    version = cat(lox_version(), "\n"),
    simulate = lox_simulate_cdm(n = as.integer(if (is.null(o$n)) 5000 else o$n), out_dir = need("out"),
                                seed = as.integer(if (is.null(o$seed)) 1 else o$seed),
                                date_res = if (is.null(o[["date-res"]])) "day" else o[["date-res"]]),
    qc = { cfg <- lox_config(need("config")); cdm <- lox_read_cdm(cfg); q <- lox_qc(cdm, cfg)
           print(q, nrows = 100); dir.create(cfg$paths$output_dir, recursive = TRUE, showWarnings = FALSE)
           data.table::fwrite(q, file.path(cfg$paths$output_dir, "qc.csv")) },
    report = { cfg <- lox_config(need("config")); lox_report_local(cfg$paths$output_dir, consensus_dir = o$consensus) },
    round1 = lox_run_round1(need("config")),
    round2 = lox_run_round2(need("config"), need("consensus")),
    `central-consensus` = lox_run_central_consensus(sites(), need("out"), need("config")),
    `central-meta`      = lox_run_central_meta(sites(), need("out"), need("config")),
    `central-report`    = lox_report_central(need("out")),
    compare = lox_compare_runs(need("a"), need("b"), labels = if (is.null(o$labels)) c("anterior", "actual") else strsplit(o$labels, ",")[[1]],
                               out_file = if (is.null(o$out)) file.path(need("b"), "comparacion.html") else o$out),
    stop("Orden desconocida: ", cmd, ". Usa 'loxofed help'.", call. = FALSE))
  invisible()
}
