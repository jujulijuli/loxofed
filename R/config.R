#' Configuracion por defecto
#'
#' Todo lo que debe ser IGUAL en los 7 sitios (cohorte, modelo, estructura) entra en
#' `lox_config_fingerprint()`; el nodo central rechaza combinar sitios con huellas distintas.
#' Lo que es propio de cada sitio (rutas, nombres de columna, mapa de sexo, semilla) no.
#' @export
lox_default_config <- function() list(
  site_id = "SITE",
  paths   = list(pacientes = NULL, enfermedades = NULL, output_dir = "/output"),
  cdm = list(
    sep  = "|",
    cols = list(id = "Id_pac", birth = "F_Nac", sex = "Sexo", death = "F_Exitus",
                baja = "F_Baja", dx_date = "F_EnfC", dx_code = "EnfC"),
    # Recodificacion de Sexo a {Mujer, Hombre}; "Mujer" es la referencia en todos los modelos.
    # El modelo de datos comun codifica Sexo como H / M (H = hombre, M = mujer). Es deliberadamente estricto:
    # cualquier otro valor (p. ej. "F", "V", "1") da ERROR en `qc` en lugar de asignarse en silencio.
    sex_map = list(female = c("M", "m", "Mujer", "mujer"),
                   male   = c("H", "h", "Hombre", "hombre"))),
  cohort = list(year_first = 2012, year_last = 2022, age_min = 60, age_max = NA,
                min_chronic = 2, exclude_baja_baseline = TRUE,
                # Tope de personas: muestra aleatoria de la cohorte ELEGIBLE (ya con >= min_chronic) y de sus
                # enfermedades. NA = sin tope. Pensado para sitios muy grandes (memoria/tiempo). Es propio del
                # sitio: NO entra en la huella, y el tamano real queda registrado en los manifiestos.
                sample_n = NA, sample_seed = 2026),
  model = list(prev_min = 0.01,            # suelo de prevalencia basal para ser nodo
               min_events_node = 100,      # incidentes minimos para ajustar la CPD de un nodo
               min_events_hazard = 20,     # muertes minimas entre portadores para entrar al hazard
               # Spline de edad con nudos FIJOS y comunes: sin esto los coeficientes de edad
               # no son comparables ni meta-analizables entre sitios.
               age_knots = c(70, 80, 90), age_boundary = c(60, 105),
               # Mortalidad (capa 2):
               #  "exact"       tiempo-persona con fechas exactas: el seguimiento se parte en cada
               #                diagnostico y cada `time_split_months`; estado variable en el tiempo;
               #                modelo de Poisson con offset (exp(beta) = HR). Cuenta TODAS las muertes.
               #  "annual_next" replica el informe de Navarra (foto a 1-ene, muerte en t+1).
               mortality = "exact", time_split_months = 12,
               # Retardo (dias) hasta que un diagnostico cuenta como exposicion; 0 = inmediato.
               # Util como sensibilidad ante diagnosticos codificados junto a la muerte.
               dx_lag_days = 0,
               # Efecto calendario en la mortalidad (solo `mortality = "exact"`): spline natural del anio con
               # nudos FIJOS y comunes (sin esto la mortalidad de 2012 y la de 2022 se mezclan en un unico nivel)
               # mas un indicador de los anios pandemicos. `enabled = FALSE` reproduce el modelo sin calendario.
               # Se ajusta ademas el modelo sin calendario para mostrar en el informe cuanto cambia todo.
               calendar = list(enabled = TRUE, knots = c(2015, 2018), pandemic_years = c(2020L, 2021L))),
  structure = list(threshold = 0.60, boot_R = 100, subsample = 60000, dir_min = 0.55,
                   algorithm = "hc", score = "bde", iss = 1, seed = 2026),
  consensus = list(threshold = 0.60, min_site_frac = 0.5),
  meta = list(min_sites = 3),               # nodo central: sitios minimos para dar resultado combinado
  report = list(auto = TRUE),              # generar informe_local.html al acabar la ronda 2
  disclosure = list(min_cell = 10)
)

#' Lee la configuracion YAML y la fusiona con los valores por defecto
#' @param path ruta a un YAML (opcional)
#' @param ... sobrescrituras como lista anidada, p.ej. `cohort = list(age_min = 65)`
#' @export
lox_config <- function(path = NULL, ...) {
  cfg <- lox_default_config()
  if (!is.null(path)) {
    if (!file.exists(path)) stop("No existe el fichero de configuracion: ", path, call. = FALSE)
    cfg <- utils::modifyList(cfg, yaml::read_yaml(path))
  }
  ov <- list(...); if (length(ov)) cfg <- utils::modifyList(cfg, ov)
  .lox_validate_config(cfg)
}

.lox_validate_config <- function(cfg) {
  # YAML `~` se lee como NULL y modifyList() lo ELIMINA: se restituye como NA ("sin limite")
  for (k in c("age_min", "age_max", "sample_n")) if (is.null(cfg$cohort[[k]])) cfg$cohort[k] <- list(NA)
  co <- cfg$cohort
  if (!(co$year_last > co$year_first)) stop("cohort: year_last debe ser > year_first", call. = FALSE)
  if (!is.na(co$age_min) && !is.na(co$age_max) && co$age_min >= co$age_max)
    stop("cohort: age_min debe ser < age_max", call. = FALSE)
  m <- cfg$model
  if (length(m$age_boundary) != 2 || m$age_boundary[1] >= m$age_boundary[2])
    stop("model$age_boundary debe ser c(min, max)", call. = FALSE)
  if (any(m$age_knots <= m$age_boundary[1] | m$age_knots >= m$age_boundary[2]))
    stop("model$age_knots debe quedar dentro de age_boundary", call. = FALSE)
  if (!cfg$model$mortality %in% c("exact", "annual_next"))
    stop("model$mortality debe ser 'exact' o 'annual_next'", call. = FALSE)
  if (!cfg$model$time_split_months %in% c(1, 2, 3, 4, 6, 12))
    stop("model$time_split_months debe dividir 12 (1, 2, 3, 4, 6 o 12)", call. = FALSE)
  cal <- m$calendar
  if (!is.null(cal) && isTRUE(cal$enabled)) {
    if (any(cal$knots <= co$year_first | cal$knots >= co$year_last + 1))
      stop("model$calendar$knots debe quedar dentro de (year_first, year_last + 1)", call. = FALSE)
    if (any(!is.na(cal$pandemic_years) & (cal$pandemic_years < co$year_first | cal$pandemic_years > co$year_last)))
      warning("model$calendar$pandemic_years fuera del periodo de seguimiento; esos anios no aportan", call. = FALSE)
  }
  if (cfg$model$dx_lag_days < 0) stop("model$dx_lag_days debe ser >= 0", call. = FALSE)
  if (!is.na(co$sample_n) && co$sample_n < 1000)
    stop("cohort$sample_n debe ser NA (sin tope) o >= 1000", call. = FALSE)
  if (cfg$disclosure$min_cell < 1) stop("disclosure$min_cell debe ser >= 1", call. = FALSE)
  if (cfg$disclosure$min_cell < 5)
    warning("disclosure$min_cell < 5: umbral inferior a las practicas habituales", call. = FALSE)
  cfg
}

#' Huella de los parametros que deben coincidir entre sitios
#' @export
lox_config_fingerprint <- function(cfg) {
  x <- cfg[c("cohort", "model", "structure", "consensus")]
  x$structure$seed <- NULL                      # la semilla puede diferir entre sitios
  x$cohort$sample_n <- NULL; x$cohort$sample_seed <- NULL   # el tope de personas es decision del sitio
  canon <- function(v) if (is.list(v)) lapply(v[order(names(v))], canon) else v   # independiente del orden de claves
  digest::digest(jsonlite::toJSON(canon(x), auto_unbox = TRUE, digits = NA, null = "null"), algo = "md5")
}
