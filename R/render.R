#' Renderiza el informe central en Quarto (en el nodo central, en local)
#'
#' Usa la plantilla `inst/quarto/informe_central.qmd`. Requiere Quarto instalado en la maquina del central
#' (los sitios NO lo necesitan: su informe local es HTML en R puro) y los paquetes de R que lista la plantilla.
#' @param central_dir carpeta de salida del central (con `consensus/`, `descriptive/` y `meta/`)
#' @param out_file nombre del HTML (por defecto `informe_central_quarto.html` dentro de `central_dir`)
#' @param quarto ruta al ejecutable de Quarto (por defecto el que haya en el PATH)
#' @export
lox_render_central <- function(central_dir, out_file = "informe_central_quarto.html", quarto = Sys.which("quarto")) {
  if (!nzchar(quarto)) stop("No se encuentra Quarto (https://quarto.org). Instalalo en la maquina del nodo central.", call. = FALSE)
  central_dir <- normalizePath(central_dir, mustWork = TRUE)
  if (!file.exists(file.path(central_dir, "meta", "meta_coefs.csv"))) stop("Falta meta/meta_coefs.csv: ejecuta central-meta antes", call. = FALSE)
  tpl <- system.file("quarto", "informe_central.qmd", package = "loxofed")
  if (!nzchar(tpl)) stop("No se encuentra la plantilla del informe en el paquete", call. = FALSE)
  work <- file.path(central_dir, "informe_quarto"); dir.create(work, showWarnings = FALSE)   # el .qmd se copia aqui (editable)
  qmd <- file.path(work, "informe_central.qmd"); file.copy(tpl, qmd, overwrite = TRUE)
  old <- setwd(work); on.exit(setwd(old), add = TRUE)            # Quarto resuelve las rutas relativas desde aqui
  st <- system2(quarto, c("render", "informe_central.qmd", "-P", shQuote(paste0("central_dir:", central_dir))))
  if (st != 0) stop("quarto render fallo (codigo ", st, ")", call. = FALSE)
  out <- file.path(central_dir, out_file)
  file.copy(file.path(work, "informe_central.html"), out, overwrite = TRUE)
  .lox_msg("Informe Quarto -> %s", out); invisible(out)
}
