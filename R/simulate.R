#' Genera un CDM sintetico (personas + enfermedades) con una verdad conocida
#'
#' Sirve para probar la instalacion del contenedor en un sitio sin tocar datos reales y para los
#' tests del paquete. Las progresiones verdaderas estan en el atributo `truth`.
#' @param n numero de personas
#' @param out_dir carpeta donde escribir `pacientes.csv` y `enfermedades.csv`
#' @param date_res resolucion de las fechas escritas: "day", "month" o "year"
#' @param tau desviacion tipica (escala logit/log) de la heterogeneidad entre sitios: cada sitio desplaza sus
#'   progresiones verdaderas y los efectos de las enfermedades sobre la mortalidad con N(0, tau). 0 = sin heterogeneidad.
#' @param het_seed semilla de los desplazamientos (distinta por sitio); si es NULL usa `seed`
#' @export
lox_simulate_cdm <- function(n = 5000, out_dir, seed = 1, date_res = c("day", "month", "year"),
                             tau = 0, het_seed = NULL) {
  date_res <- match.arg(date_res); set.seed(if (is.null(het_seed)) seed else het_seed)
  shift_inc <- if (tau > 0) stats::rnorm(8, 0, tau) else rep(0, 8)
  shift_mort <- if (tau > 0) stats::rnorm(5, 0, tau) else rep(0, 5)
  set.seed(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dis <- c("HYPERTENSION", "DIABETES", "HEART_FAILURE", "CHRONIC_KIDNEY_DISEASES", "ATRIAL_FIBRILLATION",
           "CEREBROVASCULAR_DISEASE", "DEMENTIA", "DEPRESSION_AND_MOOD_DISEASES",
           "OSTEOARTHRITIS_AND_OTHER_DEGENERATIVE_JOINT_DISEASES", "DYSLIPIDEMIA")
  prev0 <- c(0.60, 0.30, 0.06, 0.10, 0.10, 0.07, 0.04, 0.20, 0.50, 0.50)
  a_inc <- c(-3.0, -3.5, -4.5, -4.2, -4.0, -4.3, -4.8, -3.5, -3.0, -3.0)
  s_age <- c(0.03, 0.03, 0.05, 0.05, 0.06, 0.05, 0.12, 0.01, 0.03, 0.01)
  names(prev0) <- names(a_inc) <- names(s_age) <- dis
  truth <- data.table(from = c("HYPERTENSION", "DIABETES", "HEART_FAILURE", "ATRIAL_FIBRILLATION",
                               "HEART_FAILURE", "DEPRESSION_AND_MOOD_DISEASES", "HYPERTENSION", "DIABETES"),
                      to   = c("HEART_FAILURE", "CHRONIC_KIDNEY_DISEASES", "CHRONIC_KIDNEY_DISEASES",
                               "CEREBROVASCULAR_DISEASE", "ATRIAL_FIBRILLATION", "DEMENTIA",
                               "CEREBROVASCULAR_DISEASE", "HEART_FAILURE"),
                      beta = c(1.2, 1.0, 0.8, 1.0, 0.7, 0.8, 0.5, 0.6))
  truth[, beta_base := beta][, beta := beta_base + shift_inc]
  mort <- c(HEART_FAILURE = 0.7, DEMENTIA = 0.8, CHRONIC_KIDNEY_DISEASES = 0.4, CEREBROVASCULAR_DISEASE = 0.3,
            ATRIAL_FIBRILLATION = 0.2)
  mort_site <- mort + shift_mort
  birth <- sample(1920:1952, n, TRUE); male <- stats::rbinom(n, 1, 0.45)
  age12 <- 2012 - birth
  P <- sapply(dis, function(d) stats::rbinom(n, 1, stats::plogis(stats::qlogis(prev0[d]) + 0.03 * (age12 - 70)))) == 1
  colnames(P) <- dis
  fmt <- function(d) switch(date_res, day = format(d, "%Y-%m-%d"), month = format(d, "%Y-%m"), year = format(d, "%Y"))
  ev <- list()
  for (d in dis) {
    w <- which(P[, d])
    if (length(w)) ev[[length(ev) + 1L]] <- data.table(i = w, code = d,
        date = as.Date("2000-01-01") + sample(0:(365 * 12 - 2), length(w), TRUE))
  }
  # Simulacion mensual (casi tiempo continuo): incidencia y muerte se evaluan cada mes con el estado
  # vigente, de modo que el estado cambia en la fecha del diagnostico y la verdad es de tasa instantanea.
  alive <- rep(TRUE, n); death <- rep(as.Date(NA), n)
  pm <- function(p) 1 - (1 - p)^(1 / 12)                              # probabilidad anual -> mensual
  for (Y in 2012:2022) for (mo in 1:12) {
    m0 <- as.Date(sprintf("%d-%02d-01", Y, mo)); age <- Y + (mo - 0.5) / 12 - birth - 0.5; newP <- P
    for (d in dis) {
      eta <- a_inc[d] + s_age[d] * (age - 70)
      for (r in which(truth$to == d)) eta <- eta + truth$beta[r] * P[, truth$from[r]]
      w <- which(alive & !P[, d] & (stats::runif(n) < pm(stats::plogis(eta))))
      if (length(w)) { ev[[length(ev) + 1L]] <- data.table(i = w, code = d, date = m0 + sample(0:27, length(w), TRUE)); newP[w, d] <- TRUE }
    }
    eta <- -4.7 + 0.09 * (age - 70) + 0.2 * male + mort_site[["HEART_FAILURE"]] * P[, "HEART_FAILURE"] +
           mort_site[["DEMENTIA"]] * P[, "DEMENTIA"] + mort_site[["CHRONIC_KIDNEY_DISEASES"]] * P[, "CHRONIC_KIDNEY_DISEASES"] +
           mort_site[["CEREBROVASCULAR_DISEASE"]] * P[, "CEREBROVASCULAR_DISEASE"] +
           mort_site[["ATRIAL_FIBRILLATION"]] * P[, "ATRIAL_FIBRILLATION"]
    die <- alive & (stats::runif(n) < pm(stats::plogis(eta)))
    death[die] <- m0 + sample(0:27, sum(die), TRUE); alive[die] <- FALSE
    P <- newP
  }
  e <- rbindlist(ev); ids <- sprintf("P%06d", seq_len(n)); e[, id := ids[i]]
  e <- e[is.na(death[i]) | date <= death[i]]                       # sin diagnosticos posteriores a la muerte
  pac <- data.table(Id_pac = ids, F_Nac = birth, Sexo = ifelse(male == 1, "Hombre", "Mujer"),
                    F_Exitus = ifelse(is.na(death), "", fmt(death)), F_Baja = "", CCAA = "SIM")
  enf <- e[, .(Id_pac = id, F_EnfC = fmt(date), EnfC = code)]
  data.table::fwrite(pac, file.path(out_dir, "pacientes.csv"), sep = "|")
  data.table::fwrite(enf, file.path(out_dir, "enfermedades.csv"), sep = "|")
  invisible(structure(list(dir = out_dir), truth = truth, mortality_truth = mort_site, tau = tau))
}
