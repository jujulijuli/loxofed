# Cohorte, array de estados y tabla de transiciones persona-anio --------------------------------

#' Construye la cohorte LOXO y reduce las enfermedades a su primer diagnostico
#'
#' Replica el informe de Navarra: viva a `year_first`, banda de edad, >= `min_chronic`
#' condiciones del catalogo con diagnostico <= 1 de enero de `year_first`.
#' @export
lox_build_cohort <- function(cdm, cfg) {
  co <- cfg$cohort; by <- co$year_first; base <- as.Date(sprintf("%d-01-01", by))
  pac <- copy(cdm$pac)
  pac <- pac[!is.na(birth_year) & !is.na(sex) & !duplicated(id)]
  pac[, edad_base := by - birth_year]
  pac[, anio_muerte := as.integer(format(death_date, "%Y"))]
  pac[, anio_baja   := as.integer(format(baja_date, "%Y"))]
  pac <- pac[is.na(anio_muerte) | anio_muerte >= by]                       # viva al basal
  if (isTRUE(co$exclude_baja_baseline)) pac <- pac[is.na(anio_baja) | anio_baja > by]
  if (!is.na(co$age_min)) pac <- pac[edad_base >= co$age_min]
  if (!is.na(co$age_max)) pac <- pac[edad_base <= co$age_max]

  enf <- cdm$enf[!is.na(dx_date) & code %in% LOXO_DIS60 & id %in% pac$id]
  setorder(enf, id, code, dx_date)
  enf <- unique(enf, by = c("id", "code"))                                  # primer diagnostico
  basal <- enf[dx_date <= base, .N, by = id][N >= co$min_chronic]
  pac <- pac[id %in% basal$id]
  enf <- enf[id %in% pac$id]
  pac[, Sexo := factor(sex, levels = c("Mujer", "Hombre"))]
  n_elig <- nrow(pac)
  if (!is.na(co$sample_n) && n_elig > co$sample_n) {            # muestra aleatoria reproducible de personas
    set.seed(co$sample_seed)
    pac <- pac[sort(sample.int(n_elig, co$sample_n))]
    enf <- enf[id %in% pac$id]                                  # y SUS enfermedades
  }
  list(n_eligible = n_elig, sample_weight = n_elig / nrow(pac),
       cohorte = pac[, .(id, Sexo, edad_base, anio_muerte, anio_baja, death_date, baja_date)],
       enf = enf[, .(id, code, dx_date, dx_res)])
}

#' Array de estados (paciente x enfermedad x anio): foto a 1 de enero, enfermedades absorbentes
#'
#' Vectorizado como en el informe (`findInterval` + `sequence`): sin bucles por fila.
#' @export
lox_build_states <- function(cohort, cfg) {
  yrs <- cfg$cohort$year_first:cfg$cohort$year_last; K <- length(yrs)
  refs <- as.numeric(as.Date(sprintf("%d-01-01", yrs)))
  ids <- cohort$cohorte$id; enf <- cohort$enf
  idx <- setNames(seq_along(ids), ids); posd <- setNames(seq_along(LOXO_DIS60), LOXO_DIS60)
  i_vec <- idx[enf$id]; j_vec <- posd[enf$code]
  k0 <- findInterval(as.numeric(enf$dx_date) - 0.5, refs) + 1L              # primer anio presente
  ok <- !is.na(i_vec) & !is.na(j_vec) & k0 <= K
  reps <- K - k0[ok] + 1L
  arr <- array(0L, dim = c(length(ids), length(LOXO_DIS60), K),
               dimnames = list(ids, LOXO_DIS60, as.character(yrs)))
  arr[cbind(rep(i_vec[ok], reps), rep(j_vec[ok], reps), sequence(reps, from = k0[ok]))] <- 1L
  arr
}

#' Tabla de transiciones persona-anio (data.table)
#'
#' Una fila por (paciente, anio t) vivo en t y con t+1 observado. Columnas: `<nodo>_t0`,
#' `<nodo>_t1`, `death_next`, `edad_t0`, `Sexo`. La convencion de `death_next` es la del informe de
#' Navarra: muerte en el anio t+1 entre quienes sobreviven el anio t (ver README, "Convencion").
#' @param nodes nodos a incluir (por defecto los de prevalencia basal >= `model$prev_min`)
#' @export
lox_build_transitions <- function(arr, cohort, nodes, cfg) {
  yrs <- cfg$cohort$year_first:cfg$cohort$year_last
  panel <- rbindlist(lapply(seq_along(yrs), function(k) {
    m <- arr[, nodes, k]
    if (!is.matrix(m)) m <- matrix(m, ncol = length(nodes), dimnames = list(dimnames(arr)[[1]], nodes))
    dt <- as.data.table(m, keep.rownames = "id"); dt[, anio := yrs[k]]; dt
  }))
  panel <- merge(panel, cohort$cohorte[, .(id, Sexo, edad_base, anio_muerte, anio_baja)], by = "id")
  panel[, muerto := as.integer(!is.na(anio_muerte) & anio == anio_muerte)]
  panel <- panel[is.na(anio_muerte) | anio <= anio_muerte]
  setorder(panel, id, anio)
  v0 <- paste0(nodes, "_t0"); v1 <- paste0(nodes, "_t1")
  panel[, death_next := shift(muerto, 1L, type = "lead"), by = id]
  panel[, (v1) := shift(.SD, 1L, type = "lead"), by = id, .SDcols = nodes]
  tr <- panel[muerto == 0 & !is.na(death_next)]
  # Censura por baja administrativa (migracion, cambio de sistema): el desenlace del anio t+1 solo
  # es observable si la baja es posterior a t+1, o si en t+1 hubo muerte registrada (muerte <= baja).
  tr <- tr[is.na(anio_baja) | (anio + 1L) < anio_baja | (death_next == 1L & (anio + 1L) <= anio_baja)]
  tr[, edad_t0 := edad_base + (anio - cfg$cohort$year_first)]
  setnames(tr, nodes, v0)
  tr[, count_t0 := rowSums(.SD), .SDcols = v0]
  tr[, c("anio_muerte", "anio_baja", "muerto", "edad_base") := NULL]
  tr[]
}

#' Base de spline natural con nudos FIJOS (comun a todos los sitios)
#' @return data.table con columnas `edad_ns1..edad_nsK`
#' @export
lox_age_basis <- function(edad, cfg) {
  B <- splines::ns(edad, knots = cfg$model$age_knots, Boundary.knots = cfg$model$age_boundary)
  B <- as.data.table(unclass(B)); setnames(B, paste0("edad_ns", seq_len(ncol(B)))); B
}
