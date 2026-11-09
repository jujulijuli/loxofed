# Utilidades comunes ---------------------------------------------------------------------------

utils::globalVariables(c(".", ".N", ".SD", ":=", "id", "code", "dx_date", "anio", "muerto",
                         "anio_muerte", "death_next", "Sexo", "edad2012", "edad_base",
                         "n_cronicas", "dx_res", "birth_year", "sex", "strength", "direction",
                         "site", "capa", "nodo", "termino", "beta", "se", "k", "from", "to"))

lox_version <- function() as.character(utils::packageVersion("loxofed"))

.lox_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " [loxofed] ", sprintf(...))

# Parseo tolerante de fechas con resolución variable. Devuelve la fecha imputada y la resolución
# detectada de cada elemento ("year" -> 1 jul, "month" -> día 15, "day" -> exacta).
# La resolución importa para la validación de dirección y la latencia (ver README).
.lox_parse_date <- function(x) {
  x   <- trimws(as.character(x)); x[x %in% c("", "NA")] <- NA_character_
  res <- ifelse(grepl("^[0-9]{4}$", x), "year",
         ifelse(grepl("^[0-9]{4}-[0-9]{2}$", x), "month",
         ifelse(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}", x), "day", NA_character_)))
  txt <- ifelse(res == "year",  paste0(x, "-07-01"),
         ifelse(res == "month", paste0(x, "-15"),
         ifelse(res == "day",   substr(x, 1, 10), NA_character_)))
  list(date = suppressWarnings(as.Date(txt, format = "%Y-%m-%d")), res = res)
}

.lox_year <- function(x) suppressWarnings(as.integer(substr(trimws(as.character(x)), 1, 4)))

.lox_need <- function(x, what) if (is.null(x)) stop("Falta en la configuración: ", what, call. = FALSE)
