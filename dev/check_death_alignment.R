# Contraste, en datos sinteticos con verdad conocida, de tres formas de alinear la mortalidad:
#   anual t -> t+1   convencion del informe de Navarra (estado a 1-ene de t, muerte en t+1 entre quienes sobreviven t)
#   mismo anio       estado a 1-ene de t, muerte DENTRO de t (cortes anuales; el estado no cambia en el anio)
#   exacto           tiempo-persona con fechas exactas (el estado cambia en cada diagnostico)
# La "verdad" (columna `OR_sim`) es exp(beta) en la escala logit anual del simulador; en escala de tasa
# instantanea el HR verdadero es algo menor cuando la mortalidad basal es alta.
#   Rscript dev/check_death_alignment.R
suppressMessages({ library(loxofed); library(data.table) })
d <- file.path(tempdir(), "sim"); lox_simulate_cdm(n = 40000, out_dir = d, seed = 7)
cfg <- lox_config(paths = list(pacientes = file.path(d, "pacientes.csv"), enfermedades = file.path(d, "enfermedades.csv"),
                               output_dir = file.path(tempdir(), "out")))
cdm <- lox_read_cdm(cfg); co <- lox_build_cohort(cdm, cfg); arr <- lox_build_states(co, cfg)
nodes <- lox_prevalence(arr)[prev >= 0.01, node]

d2 <- lox_build_transitions(arr, co, nodes, cfg)
anual <- lox_fit_layer2(d2, nodes, cfg)$coefs[nodo != "muerte", .(nodo, HR_anual_t_t1 = exp(beta))]

yrs <- 2012:2022
pn <- rbindlist(lapply(seq_along(yrs), function(k) { dt <- as.data.table(arr[, nodes, k], keep.rownames = "id"); dt[, anio := yrs[k]]; dt }))
pn <- merge(pn, co$cohorte, by = "id"); pn <- pn[is.na(anio_muerte) | anio_muerte >= anio]
pn[, death_next := as.integer(!is.na(anio_muerte) & anio_muerte == anio)]
pn[, edad_t0 := edad_base + anio - 2012]; setnames(pn, nodes, paste0(nodes, "_t0"))
mismo <- lox_fit_layer2(pn, nodes, cfg)$coefs[nodo != "muerte", .(nodo, HR_mismo_anio = exp(beta))]

pt <- lox_build_persontime(co, nodes, cfg)
exacto <- lox_fit_layer2_exact(pt, nodes, cfg)$coefs[nodo != "muerte", .(nodo, HR_exacto = exp(beta))]

truth <- data.table(nodo = c("HEART_FAILURE", "DEMENTIA", "CHRONIC_KIDNEY_DISEASES", "CEREBROVASCULAR_DISEASE", "ATRIAL_FIBRILLATION"),
                    OR_sim = exp(c(.7, .8, .4, .3, .2)))
r <- Reduce(function(x, y) merge(x, y, by = "nodo"), list(truth, anual, mismo, exacto))
r[, (2:5) := lapply(.SD, round, 2), .SDcols = 2:5]; print(r)
muertes <- co$cohorte[!is.na(death_date) & death_date < as.Date("2023-01-01"), .N]
cat(sprintf("Muertes en la ventana: %d | contadas: anual t->t+1 = %d, mismo anio = %d, exacto = %d\n",
            muertes, sum(d2$death_next), sum(pn$death_next), sum(pt$event)))
