# ==========================================================
# OBJETIVO 1
# ANÁLISIS MULTINIVEL DE MORTALIDAD HOSPITALARIA
# Pacientes de Neurología y Neurocirugía atendidos en hospitales de India
# ==========================================================

# Objetivo general:
# Evaluar si la mortalidad hospitalaria en pacientes neurológicos y
# neuroquirúrgicos se explica únicamente por características individuales,
# sociodemográficas y clínicas, o también por diferencias entre hospitales.
#
# Decisiones metodológicas centrales:
# - nombre_categoria se usa para restringir la muestra a Neurología y Neurocirugía.
# - especialidad_neuro se conserva para describir la muestra.
# - procedimiento_modelo es la variable clínica principal.
# - nombre_casta se evalúa como determinante sociodemográfico de nivel paciente.
# - nombre_casta NO se modela como efecto aleatorio.
# - tipo_hospital se evalúa como variable contextual, pero no se conserva en el
#   modelo final si no mejora el ajuste.
# - monto_reclamado se analiza solo como sensibilidad.
# - hospital es el nivel contextual principal.
# - distrito se evalúa solo como estructura candidata.
# - el ranking hospitalario se interpreta como desviación ajustada, no como
#   indicador directo de calidad asistencial.


# ==========================================================
# 0. PAQUETES
# ==========================================================

library(tidyverse)
library(janitor)
library(lubridate)
library(arrow)
library(lme4)
library(car)
library(pROC)
library(broom.mixed)
library(DHARMa)
library(splines)


# ==========================================================
# 1. FUNCIONES AUXILIARES
# ==========================================================

calcular_icc_logistico <- function(modelo) {

  var_random <- as.data.frame(VarCorr(modelo)) %>%
    select(grp, vcov)

  var_total_random <- sum(var_random$vcov)
  var_logistica <- pi^2 / 3

  tabla_icc <- var_random %>%
    mutate(
      ICC = vcov / (var_total_random + var_logistica),
      ICC_pct = ICC * 100
    )

  list(
    varianzas = var_random,
    ICC_total = var_total_random / (var_total_random + var_logistica),
    ICC_total_pct = 100 * var_total_random / (var_total_random + var_logistica),
    ICC_por_nivel = tabla_icc
  )
}

extraer_or_wald <- function(modelo) {
  exp(cbind(
    OR = fixef(modelo),
    confint(modelo, parm = "beta_", method = "Wald")
  ))
}

evaluar_modelo <- function(modelo) {
  list(
    singular = isSingular(modelo, tol = 1e-4),
    convergencia = modelo@optinfo$conv$lme4$messages
  )
}

tabla_modelos <- function(...) {

  modelos <- list(...)
  nombres <- names(modelos)

  tibble(
    modelo = nombres,
    AIC = map_dbl(modelos, AIC),
    BIC = map_dbl(modelos, BIC),
    logLik = map_dbl(modelos, ~ as.numeric(logLik(.x))),
    devianza = map_dbl(modelos, deviance),
    singular = map_lgl(modelos, isSingular)
  )
}


# ==========================================================
# 2. CARGA, LIMPIEZA Y RECODIFICACIÓN
# ==========================================================

datos <- read_parquet(
  "~/UNAB/taller big data y salud/pacientes/ntrarogyaseva.parquet"
) %>%
  clean_names() %>%
  rename(
    edad = age,
    sexo = sex,
    nombre_casta = caste_name,
    nombre_categoria = category_name,
    nombre_cirugia = surgery,
    nombre_hospital = hosp_name,
    tipo_hospital = hosp_type,
    distrito_hospital = hosp_district,
    fecha_preautorizacion = preauth_date,
    monto_preautorizado = preauth_amt,
    fecha_reclamo = claim_date,
    monto_reclamado = claim_amount,
    fecha_cirugia = surgery_date,
    fecha_egreso = discharge_date,
    mortalidad = mortality_y_n
  ) %>%
  mutate(
    across(
      c(
        sexo,
        mortalidad,
        nombre_casta,
        nombre_categoria,
        nombre_cirugia,
        nombre_hospital,
        tipo_hospital,
        distrito_hospital
      ),
      ~ str_squish(str_to_lower(as.character(.x)))
    ),

    sexo = factor(
      sexo,
      levels = c("male", "female"),
      labels = c("Masculino", "Femenino")
    ),

    tipo_hospital = factor(
      tipo_hospital,
      levels = c("c", "g"),
      labels = c("Privado", "Público")
    ),

    mortalidad = factor(
      mortalidad,
      levels = c("no", "yes")
    ),

    mortalidad_binaria = case_when(
      mortalidad == "yes" ~ 1,
      mortalidad == "no"  ~ 0,
      TRUE ~ NA_real_
    ),

    edad = ifelse(edad < 0 | edad > 120, NA, edad),

    across(
      c(
        nombre_casta,
        nombre_categoria,
        nombre_cirugia,
        nombre_hospital,
        distrito_hospital
      ),
      as.factor
    )
  )


# ==========================================================
# 3. RESTRICCIÓN A NEUROLOGÍA Y NEUROCIRUGÍA
# ==========================================================

categorias_neuro <- c("neurology", "neurosurgery")

datos_neuro <- datos %>%
  filter(as.character(nombre_categoria) %in% categorias_neuro) %>%
  mutate(
    especialidad_neuro = case_when(
      as.character(nombre_categoria) == "neurology" ~ "Neurología",
      as.character(nombre_categoria) == "neurosurgery" ~ "Neurocirugía",
      TRUE ~ NA_character_
    ),
    especialidad_neuro = factor(
      especialidad_neuro,
      levels = c("Neurología", "Neurocirugía")
    )
  ) %>%
  droplevels()

datos_neuro %>%
  count(especialidad_neuro)


# ==========================================================
# 4. VARIABLES DERIVADAS
# ==========================================================

datos_neuro <- datos_neuro %>%
  mutate(
    fecha_cirugia = as_date(fecha_cirugia),
    fecha_egreso = as_date(fecha_egreso),

    dias_internacion = as.numeric(fecha_egreso - fecha_cirugia),
    dias_internacion = ifelse(dias_internacion < 0, NA, dias_internacion),

    monto_reclamado = ifelse(monto_reclamado < 0, NA, monto_reclamado),
    monto_preautorizado = ifelse(monto_preautorizado < 0, NA, monto_preautorizado),

    diferencia_montos = monto_reclamado - monto_preautorizado,

    ratio_montos = monto_reclamado / monto_preautorizado,
    ratio_montos = ifelse(is.infinite(ratio_montos), NA, ratio_montos),

    log_monto_reclamado = log1p(monto_reclamado),
    log_monto_preautorizado = log1p(monto_preautorizado)
  )


# ==========================================================
# 5. BASE ANALÍTICA PRINCIPAL
# ==========================================================

datos_modelo <- datos_neuro %>%
  select(
    mortalidad_binaria,
    edad,
    sexo,
    nombre_casta,
    especialidad_neuro,
    nombre_cirugia,
    tipo_hospital,
    distrito_hospital,
    nombre_hospital,
    dias_internacion,
    monto_reclamado,
    monto_preautorizado,
    log_monto_reclamado,
    log_monto_preautorizado
  ) %>%
  drop_na(
    mortalidad_binaria,
    edad,
    sexo,
    nombre_casta,
    especialidad_neuro,
    nombre_cirugia,
    tipo_hospital,
    distrito_hospital,
    nombre_hospital
  ) %>%
  mutate(
    edad_z = as.numeric(scale(edad)),
    log_monto_reclamado_z = as.numeric(scale(log_monto_reclamado)),
    log_monto_preautorizado_z = as.numeric(scale(log_monto_preautorizado))
  ) %>%
  droplevels()

resumen_base_inicial <- datos_modelo %>%
  summarise(
    pacientes = n(),
    hospitales = n_distinct(nombre_hospital),
    distritos = n_distinct(distrito_hospital),
    castas = n_distinct(nombre_casta),
    procedimientos_originales = n_distinct(nombre_cirugia),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100
  )

resumen_base_inicial


# ==========================================================
# 6. DESCRIPTIVOS CLÍNICOS Y SOCIODEMOGRÁFICOS
# ==========================================================

resumen_especialidad <- datos_modelo %>%
  group_by(especialidad_neuro) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    procedimientos = n_distinct(nombre_cirugia),
    hospitales = n_distinct(nombre_hospital),
    .groups = "drop"
  )

resumen_especialidad

resumen_casta <- datos_modelo %>%
  group_by(nombre_casta) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    hospitales = n_distinct(nombre_hospital),
    procedimientos = n_distinct(nombre_cirugia),
    .groups = "drop"
  ) %>%
  arrange(desc(pacientes))

resumen_casta

resumen_procedimiento_original <- datos_modelo %>%
  group_by(nombre_cirugia) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    especialidades = n_distinct(especialidad_neuro),
    hospitales = n_distinct(nombre_hospital),
    .groups = "drop"
  ) %>%
  arrange(pacientes)

resumen_procedimiento_original

summary(resumen_procedimiento_original$pacientes)
summary(resumen_procedimiento_original$muertes)


# ==========================================================
# 7. AGRUPAMIENTO DE PROCEDIMIENTOS
# ==========================================================

min_pacientes_proc <- 100
min_eventos_proc <- 5

procedimientos_inestables <- resumen_procedimiento_original %>%
  filter(
    pacientes < min_pacientes_proc |
      muertes < min_eventos_proc
  ) %>%
  pull(nombre_cirugia)

datos_modelo <- datos_modelo %>%
  mutate(
    procedimiento_modelo = if_else(
      nombre_cirugia %in% procedimientos_inestables,
      "Otros procedimientos",
      as.character(nombre_cirugia)
    ),
    procedimiento_modelo = factor(procedimiento_modelo)
  ) %>%
  droplevels()

resumen_procedimiento_modelo <- datos_modelo %>%
  group_by(procedimiento_modelo) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    procedimientos_originales = n_distinct(nombre_cirugia),
    hospitales = n_distinct(nombre_hospital),
    .groups = "drop"
  ) %>%
  arrange(pacientes)

resumen_procedimiento_modelo


# ==========================================================
# 8. DEFINICIÓN DE HOSPITALES ESTABLES
# ==========================================================

resumen_hospital_todos <- datos_modelo %>%
  group_by(nombre_hospital) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    tipo_hospital = first(tipo_hospital),
    distrito_hospital = first(distrito_hospital),
    .groups = "drop"
  ) %>%
  arrange(pacientes)

resumen_hospital_todos

summary(resumen_hospital_todos$pacientes)
summary(resumen_hospital_todos$muertes)

hospitales_estables <- resumen_hospital_todos %>%
  filter(
    pacientes >= 30,
    muertes >= 5
  ) %>%
  pull(nombre_hospital)

datos_modelo_hosp_estables <- datos_modelo %>%
  filter(nombre_hospital %in% hospitales_estables) %>%
  droplevels()

resumen_base_estable <- datos_modelo_hosp_estables %>%
  summarise(
    pacientes = n(),
    hospitales = n_distinct(nombre_hospital),
    distritos = n_distinct(distrito_hospital),
    castas = n_distinct(nombre_casta),
    procedimientos = n_distinct(procedimiento_modelo),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100
  )

resumen_base_estable

resumen_procedimiento_estable <- datos_modelo_hosp_estables %>%
  group_by(procedimiento_modelo) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    .groups = "drop"
  ) %>%
  arrange(pacientes)

resumen_procedimiento_estable


# ==========================================================
# 9. ESTRUCTURA JERÁRQUICA DESCRIPTIVA
# ==========================================================

estructura_distrito <- datos_modelo_hosp_estables %>%
  group_by(distrito_hospital) %>%
  summarise(
    pacientes = n(),
    hospitales = n_distinct(nombre_hospital),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mortalidad_pct))

estructura_distrito

estructura_hospital <- datos_modelo_hosp_estables %>%
  group_by(nombre_hospital) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_cruda_pct = mean(mortalidad_binaria) * 100,
    tipo_hospital = first(tipo_hospital),
    distrito_hospital = first(distrito_hospital),
    .groups = "drop"
  ) %>%
  arrange(desc(mortalidad_cruda_pct))

estructura_hospital

summary(estructura_hospital$pacientes)
summary(estructura_hospital$muertes)
summary(estructura_hospital$mortalidad_cruda_pct)


# ==========================================================
# 10. MODELOS LOGÍSTICOS CLÁSICOS DE REFERENCIA
# ==========================================================

glm_referencia_basico <- glm(
  mortalidad_binaria ~
    edad_z +
    sexo,
  data = datos_modelo_hosp_estables,
  family = binomial()
)

summary(glm_referencia_basico)
exp(cbind(OR = coef(glm_referencia_basico), confint.default(glm_referencia_basico)))

glm_referencia_sociodemografico <- glm(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_casta,
  data = datos_modelo_hosp_estables,
  family = binomial()
)

summary(glm_referencia_sociodemografico)
exp(cbind(OR = coef(glm_referencia_sociodemografico), confint.default(glm_referencia_sociodemografico)))
car::vif(glm_referencia_sociodemografico)

glm_referencia_clinico <- glm(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo,
  data = datos_modelo_hosp_estables,
  family = binomial()
)

summary(glm_referencia_clinico)
exp(cbind(OR = coef(glm_referencia_clinico), confint.default(glm_referencia_clinico)))
car::vif(glm_referencia_clinico)

glm_referencia_contextual <- glm(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo +
    tipo_hospital,
  data = datos_modelo_hosp_estables,
  family = binomial()
)

summary(glm_referencia_contextual)
exp(cbind(OR = coef(glm_referencia_contextual), confint.default(glm_referencia_contextual)))
car::vif(glm_referencia_contextual)


# Sensibilidad clásica con monto reclamado

datos_modelo_sens_monto <- datos_modelo_hosp_estables %>%
  drop_na(log_monto_reclamado_z)

glm_referencia_monto <- glm(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo +
    log_monto_reclamado_z,
  data = datos_modelo_sens_monto,
  family = binomial()
)

summary(glm_referencia_monto)
exp(cbind(OR = coef(glm_referencia_monto), confint.default(glm_referencia_monto)))
car::vif(glm_referencia_monto)


# ==========================================================
# 11. MODELOS NULOS MULTINIVEL
# ==========================================================

# ----------------------------------------------------------
# 11.1 Modelo nulo hospitalario
# ----------------------------------------------------------

glmm_nulo_hospital <- glmer(
  mortalidad_binaria ~
    1 +
    (1 | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(glmm_nulo_hospital)

icc_nulo_hospital <- calcular_icc_logistico(glmm_nulo_hospital)

icc_nulo_hospital$varianzas
icc_nulo_hospital$ICC_total_pct
icc_nulo_hospital$ICC_por_nivel

evaluar_modelo(glmm_nulo_hospital)


# ----------------------------------------------------------
# 11.2 Modelo nulo distrito/hospital
# ----------------------------------------------------------

glmm_nulo_distrito_hospital <- glmer(
  mortalidad_binaria ~
    1 +
    (1 | distrito_hospital / nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(glmm_nulo_distrito_hospital)

icc_nulo_distrito_hospital <- calcular_icc_logistico(
  glmm_nulo_distrito_hospital
)

icc_nulo_distrito_hospital$ICC_por_nivel

anova(
  glmm_nulo_hospital,
  glmm_nulo_distrito_hospital
)

AIC(
  glmm_nulo_hospital,
  glmm_nulo_distrito_hospital
)

BIC(
  glmm_nulo_hospital,
  glmm_nulo_distrito_hospital
)

evaluar_modelo(glmm_nulo_distrito_hospital)


# ==========================================================
# 12. CONSTRUCCIÓN PROGRESIVA DEL MODELO MULTINIVEL
# ==========================================================

# ----------------------------------------------------------
# 12.1 Modelo demográfico básico: edad + sexo
# ----------------------------------------------------------

glmm_demografico <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    (1 | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(glmm_demografico)

icc_demografico <- calcular_icc_logistico(glmm_demografico)

icc_demografico$ICC_total_pct
icc_demografico$ICC_por_nivel

anova(
  glmm_nulo_hospital,
  glmm_demografico
)

AIC(
  glmm_nulo_hospital,
  glmm_demografico
)

BIC(
  glmm_nulo_hospital,
  glmm_demografico
)


# ----------------------------------------------------------
# 12.2 Modelo sociodemográfico ampliado: edad + sexo + casta
# ----------------------------------------------------------

glmm_sociodemografico <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_casta +
    (1 | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 3e5)
  )
)

summary(glmm_sociodemografico)

icc_sociodemografico <- calcular_icc_logistico(glmm_sociodemografico)

icc_sociodemografico$ICC_total_pct
icc_sociodemografico$ICC_por_nivel

anova(
  glmm_demografico,
  glmm_sociodemografico
)

AIC(
  glmm_demografico,
  glmm_sociodemografico
)

BIC(
  glmm_demografico,
  glmm_sociodemografico
)

or_glmm_sociodemografico <- extraer_or_wald(glmm_sociodemografico)
or_glmm_sociodemografico


# ----------------------------------------------------------
# 12.3 Modelo clínico base: edad + sexo + procedimiento
# ----------------------------------------------------------

glmm_clinico <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo +
    (1 | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 3e5)
  )
)

summary(glmm_clinico)

icc_clinico <- calcular_icc_logistico(glmm_clinico)

icc_clinico$ICC_total_pct
icc_clinico$ICC_por_nivel

anova(
  glmm_demografico,
  glmm_clinico
)

AIC(
  glmm_demografico,
  glmm_clinico
)

BIC(
  glmm_demografico,
  glmm_clinico
)


# ----------------------------------------------------------
# 12.4 Modelo clínico sociodemográfico: edad + sexo + casta + procedimiento
# ----------------------------------------------------------

glmm_clinico_sociodemografico <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_casta +
    procedimiento_modelo +
    (1 | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 3e5)
  )
)

summary(glmm_clinico_sociodemografico)

icc_clinico_sociodemografico <- calcular_icc_logistico(
  glmm_clinico_sociodemografico
)

icc_clinico_sociodemografico$ICC_total_pct
icc_clinico_sociodemografico$ICC_por_nivel

anova(
  glmm_clinico,
  glmm_clinico_sociodemografico
)

AIC(
  glmm_clinico,
  glmm_clinico_sociodemografico
)

BIC(
  glmm_clinico,
  glmm_clinico_sociodemografico
)

or_glmm_clinico_sociodemografico <- extraer_or_wald(
  glmm_clinico_sociodemografico
)

or_glmm_clinico_sociodemografico


# ----------------------------------------------------------
# 12.5 Modelo contextual: + tipo_hospital
# ----------------------------------------------------------

glmm_contextual <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo +
    tipo_hospital +
    (1 | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 3e5)
  )
)

summary(glmm_contextual)

icc_contextual <- calcular_icc_logistico(glmm_contextual)

icc_contextual$ICC_total_pct
icc_contextual$ICC_por_nivel

anova(
  glmm_clinico,
  glmm_contextual
)

AIC(
  glmm_clinico,
  glmm_contextual
)

BIC(
  glmm_clinico,
  glmm_contextual
)


# ----------------------------------------------------------
# 12.6 Modelo contextual ampliado: casta + procedimiento + tipo hospital
# ----------------------------------------------------------

glmm_contextual_ampliado <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_casta +
    procedimiento_modelo +
    tipo_hospital +
    (1 | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 3e5)
  )
)

summary(glmm_contextual_ampliado)

icc_contextual_ampliado <- calcular_icc_logistico(glmm_contextual_ampliado)

icc_contextual_ampliado$ICC_total_pct
icc_contextual_ampliado$ICC_por_nivel

anova(
  glmm_clinico_sociodemografico,
  glmm_contextual_ampliado
)

AIC(
  glmm_clinico_sociodemografico,
  glmm_contextual_ampliado
)

BIC(
  glmm_clinico_sociodemografico,
  glmm_contextual_ampliado
)


# ----------------------------------------------------------
# 12.7 Pendiente aleatoria para edad sobre modelo clínico parsimonioso
# ----------------------------------------------------------

glmm_edad_pendiente_aleatoria <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo +
    (1 + edad_z | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 5e5)
  )
)

summary(glmm_edad_pendiente_aleatoria)

icc_edad_pendiente_aleatoria <- calcular_icc_logistico(
  glmm_edad_pendiente_aleatoria
)

icc_edad_pendiente_aleatoria$ICC_total_pct
icc_edad_pendiente_aleatoria$ICC_por_nivel

evaluar_modelo(glmm_edad_pendiente_aleatoria)

anova(
  glmm_clinico,
  glmm_edad_pendiente_aleatoria
)

AIC(
  glmm_clinico,
  glmm_edad_pendiente_aleatoria
)

BIC(
  glmm_clinico,
  glmm_edad_pendiente_aleatoria
)


# ----------------------------------------------------------
# 12.8 Pendiente aleatoria para edad sobre modelo clínico sociodemográfico
# Se conserva solo si casta aporta al modelo clínico.
# ----------------------------------------------------------

glmm_edad_pendiente_aleatoria_casta <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_casta +
    procedimiento_modelo +
    (1 + edad_z | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 5e5)
  )
)

summary(glmm_edad_pendiente_aleatoria_casta)

icc_edad_pendiente_aleatoria_casta <- calcular_icc_logistico(
  glmm_edad_pendiente_aleatoria_casta
)

icc_edad_pendiente_aleatoria_casta$ICC_total_pct
icc_edad_pendiente_aleatoria_casta$ICC_por_nivel

evaluar_modelo(glmm_edad_pendiente_aleatoria_casta)

anova(
  glmm_edad_pendiente_aleatoria,
  glmm_edad_pendiente_aleatoria_casta
)

AIC(
  glmm_edad_pendiente_aleatoria,
  glmm_edad_pendiente_aleatoria_casta
)

BIC(
  glmm_edad_pendiente_aleatoria,
  glmm_edad_pendiente_aleatoria_casta
)


# ----------------------------------------------------------
# 12.9 Pendiente aleatoria exploratoria para sexo
# Se evalúa como sensibilidad. No se retiene salvo mejora clara y estabilidad.
# ----------------------------------------------------------

glmm_sexo_pendiente_aleatoria <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo +
    (1 + sexo | nombre_hospital),
  data = datos_modelo_hosp_estables,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 5e5)
  )
)

summary(glmm_sexo_pendiente_aleatoria)

icc_sexo_pendiente_aleatoria <- calcular_icc_logistico(
  glmm_sexo_pendiente_aleatoria
)

icc_sexo_pendiente_aleatoria$ICC_total_pct
icc_sexo_pendiente_aleatoria$ICC_por_nivel

evaluar_modelo(glmm_sexo_pendiente_aleatoria)

anova(
  glmm_clinico,
  glmm_sexo_pendiente_aleatoria
)

AIC(
  glmm_clinico,
  glmm_sexo_pendiente_aleatoria
)

BIC(
  glmm_clinico,
  glmm_sexo_pendiente_aleatoria
)


# ----------------------------------------------------------
# 12.10 Sensibilidad económica: + log_monto_reclamado_z
# ----------------------------------------------------------

glmm_clinico_sens_base <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo +
    (1 | nombre_hospital),
  data = datos_modelo_sens_monto,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 3e5)
  )
)

glmm_sensibilidad_monto <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    procedimiento_modelo +
    log_monto_reclamado_z +
    (1 | nombre_hospital),
  data = datos_modelo_sens_monto,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 3e5)
  )
)

summary(glmm_sensibilidad_monto)

icc_sensibilidad_monto <- calcular_icc_logistico(glmm_sensibilidad_monto)

icc_sensibilidad_monto$ICC_total_pct
icc_sensibilidad_monto$ICC_por_nivel

anova(
  glmm_clinico_sens_base,
  glmm_sensibilidad_monto
)

AIC(
  glmm_clinico_sens_base,
  glmm_sensibilidad_monto
)

BIC(
  glmm_clinico_sens_base,
  glmm_sensibilidad_monto
)


# ==========================================================
# 13. COMPARACIÓN OBJETIVA DE MODELOS
# ==========================================================

comparacion_modelos_principales <- tabla_modelos(
  "Nulo hospital" = glmm_nulo_hospital,
  "Demográfico" = glmm_demografico,
  "Sociodemográfico: casta" = glmm_sociodemografico,
  "Clínico: procedimiento" = glmm_clinico,
  "Clínico sociodemográfico" = glmm_clinico_sociodemografico,
  "Contextual: tipo hospital" = glmm_contextual,
  "Contextual ampliado" = glmm_contextual_ampliado,
  "Pendiente aleatoria edad" = glmm_edad_pendiente_aleatoria,
  "Pendiente aleatoria edad + casta" = glmm_edad_pendiente_aleatoria_casta,
  "Pendiente aleatoria sexo" = glmm_sexo_pendiente_aleatoria
)

comparacion_modelos_principales

comparacion_icc <- tibble(
  modelo = c(
    "Nulo hospital",
    "Demográfico",
    "Sociodemográfico: casta",
    "Clínico: procedimiento",
    "Clínico sociodemográfico",
    "Contextual: tipo hospital",
    "Contextual ampliado",
    "Pendiente aleatoria edad",
    "Pendiente aleatoria edad + casta",
    "Pendiente aleatoria sexo"
  ),
  ICC_pct = c(
    icc_nulo_hospital$ICC_total_pct,
    icc_demografico$ICC_total_pct,
    icc_sociodemografico$ICC_total_pct,
    icc_clinico$ICC_total_pct,
    icc_clinico_sociodemografico$ICC_total_pct,
    icc_contextual$ICC_total_pct,
    icc_contextual_ampliado$ICC_total_pct,
    icc_edad_pendiente_aleatoria$ICC_total_pct,
    icc_edad_pendiente_aleatoria_casta$ICC_total_pct,
    icc_sexo_pendiente_aleatoria$ICC_total_pct
  )
)

comparacion_icc

# Comparaciones secuenciales principales

anova(
  glmm_nulo_hospital,
  glmm_demografico,
  glmm_sociodemografico
)

anova(
  glmm_demografico,
  glmm_clinico
)

anova(
  glmm_clinico,
  glmm_clinico_sociodemografico
)

anova(
  glmm_clinico,
  glmm_contextual
)

anova(
  glmm_clinico_sociodemografico,
  glmm_contextual_ampliado
)

anova(
  glmm_clinico,
  glmm_edad_pendiente_aleatoria
)

anova(
  glmm_edad_pendiente_aleatoria,
  glmm_edad_pendiente_aleatoria_casta
)

anova(
  glmm_clinico,
  glmm_sexo_pendiente_aleatoria
)

# Sensibilidad con monto

comparacion_sensibilidad_monto <- tabla_modelos(
  "Clínico sin monto - base sensibilidad" = glmm_clinico_sens_base,
  "Clínico con monto - sensibilidad" = glmm_sensibilidad_monto
)

comparacion_sensibilidad_monto


# ==========================================================
# 14. SELECCIÓN DEL MODELO FINAL
# ==========================================================

# Regla de decisión:
# 1. El modelo clínico con procedimiento es el modelo fijo principal.
# 2. Casta se conserva solo si mejora el ajuste del modelo clínico y tiene
#    relevancia epidemiológica clara.
# 3. Tipo_hospital se descarta del modelo final si no mejora AIC/BIC/LRT.
# 4. La pendiente aleatoria de edad se conserva si converge, no es singular y
#    mejora AIC/BIC/LRT.
# 5. La pendiente aleatoria de sexo queda como sensibilidad si no mejora
#    claramente el ajuste.
# 6. El monto reclamado queda como sensibilidad, no como modelo principal.

modelo_final <- glmm_edad_pendiente_aleatoria

# Si casta mejora claramente el modelo con pendiente aleatoria, puede cambiarse:
if (
  !isSingular(glmm_edad_pendiente_aleatoria_casta) &&
  AIC(glmm_edad_pendiente_aleatoria_casta) < AIC(glmm_edad_pendiente_aleatoria) &&
  BIC(glmm_edad_pendiente_aleatoria_casta) < BIC(glmm_edad_pendiente_aleatoria)
) {
  modelo_final <- glmm_edad_pendiente_aleatoria_casta
}

summary(modelo_final)

or_modelo_final <- extraer_or_wald(modelo_final)

or_modelo_final

evaluar_modelo(modelo_final)


# ==========================================================
# 15. DIAGNÓSTICOS DEL MODELO FINAL
# ==========================================================

evaluar_modelo(modelo_final)

set.seed(123)

residuos_modelo_final <- simulateResiduals(
  fittedModel = modelo_final,
  n = 1000
)

plot(residuos_modelo_final)

testUniformity(residuos_modelo_final)
testDispersion(residuos_modelo_final)
testOutliers(residuos_modelo_final)

datos_modelo_hosp_estables$prob_modelo_final <- predict(
  modelo_final,
  type = "response"
)

summary(datos_modelo_hosp_estables$prob_modelo_final)


# ==========================================================
# 16. EFECTOS ALEATORIOS HOSPITALARIOS
# ==========================================================

efectos_hospital <- ranef(
  modelo_final,
  condVar = TRUE
)$nombre_hospital

efectos_hospital <- efectos_hospital %>%
  as.data.frame() %>%
  rownames_to_column("nombre_hospital")

# En modelos con pendiente aleatoria, ranef puede devolver más de una columna.
# Se conserva el intercepto hospitalario ajustado como desviación principal.

efectos_hospital <- efectos_hospital %>%
  rename(
    efecto_aleatorio_intercepto = `(Intercept)`
  )

resumen_hospital_ajustado <- datos_modelo_hosp_estables %>%
  group_by(nombre_hospital) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_cruda_pct = mean(mortalidad_binaria) * 100,
    tipo_hospital = first(tipo_hospital),
    distrito_hospital = first(distrito_hospital),
    .groups = "drop"
  ) %>%
  left_join(
    efectos_hospital,
    by = "nombre_hospital"
  ) %>%
  mutate(
    odds_ratio_hospital_ajustado = exp(efecto_aleatorio_intercepto)
  ) %>%
  arrange(desc(efecto_aleatorio_intercepto))

resumen_hospital_ajustado

hospitales_superior_esperado <- resumen_hospital_ajustado %>%
  filter(efecto_aleatorio_intercepto > 0) %>%
  arrange(desc(efecto_aleatorio_intercepto))

hospitales_inferior_esperado <- resumen_hospital_ajustado %>%
  filter(efecto_aleatorio_intercepto < 0) %>%
  arrange(efecto_aleatorio_intercepto)

hospitales_superior_esperado
hospitales_inferior_esperado

# Gráfico de interceptos aleatorios hospitalarios

resumen_hospital_ajustado %>%
  mutate(
    nombre_hospital = fct_reorder(nombre_hospital, efecto_aleatorio_intercepto)
  ) %>%
  ggplot(
    aes(
      x = nombre_hospital,
      y = efecto_aleatorio_intercepto
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  coord_flip() +
  labs(
    x = "Hospital",
    y = "Efecto aleatorio hospitalario ajustado",
    title = "Diferencias hospitalarias ajustadas en mortalidad"
  )

# Gráfico opcional de pendientes aleatorias de edad, si el modelo final las contiene

if ("edad_z" %in% names(efectos_hospital)) {

  resumen_hospital_ajustado %>%
    mutate(
      nombre_hospital = fct_reorder(nombre_hospital, edad_z)
    ) %>%
    ggplot(
      aes(
        x = nombre_hospital,
        y = edad_z
      )
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point() +
    coord_flip() +
    labs(
      x = "Hospital",
      y = "Pendiente aleatoria de edad",
      title = "Variación hospitalaria del efecto de la edad"
    )
}


# ==========================================================
# 17. AUC COMO ANÁLISIS COMPLEMENTARIO
# ==========================================================

prob_modelo_final <- predict(
  modelo_final,
  type = "response"
)

roc_modelo_final <- roc(
  response = datos_modelo_hosp_estables$mortalidad_binaria,
  predictor = prob_modelo_final,
  quiet = TRUE
)

auc_modelo_final <- auc(roc_modelo_final)

auc_modelo_final

plot(roc_modelo_final)

# Comparación opcional con GLM clásico clínico

prob_glm_referencia_clinico <- predict(
  glm_referencia_clinico,
  type = "response"
)

roc_glm_referencia_clinico <- roc(
  response = datos_modelo_hosp_estables$mortalidad_binaria,
  predictor = prob_glm_referencia_clinico,
  quiet = TRUE
)

auc_glm_referencia_clinico <- auc(roc_glm_referencia_clinico)

auc_glm_referencia_clinico


# ==========================================================
# 18. CONCLUSIÓN METODOLÓGICA DEL OBJETIVO 1
# ==========================================================

# El análisis confirmó la existencia de dependencia entre pacientes atendidos
# en el mismo hospital, justificando el uso de modelos multinivel.
#
# La estructura distrito/hospital fue evaluada, pero el distrito se conserva
# únicamente si aporta varianza adicional y mejora AIC/BIC/LRT. En caso contrario,
# el hospital se mantiene como principal nivel contextual.
#
# La edad, el sexo y el procedimiento específico constituyen el núcleo del modelo
# clínico. El procedimiento representa el principal componente de case-mix en la
# cohorte neurológica y neuroquirúrgica.
#
# La casta se evaluó como determinante sociodemográfico de nivel paciente mediante
# efecto fijo. No se modeló como efecto aleatorio porque sus categorías no
# representan una muestra aleatoria de una población mayor de niveles, sino
# categorías sociales sustantivas cuyo efecto debe estimarse explícitamente.
#
# El tipo de hospital se evaluó como variable contextual, pero no se conserva en
# el modelo final si no mejora el ajuste luego de controlar por edad, sexo,
# procedimiento y hospital.
#
# La pendiente aleatoria de edad permite evaluar si el efecto de la edad sobre
# la mortalidad varía entre hospitales. Si converge, no es singular y mejora los
# criterios de ajuste, representa una especificación final epidemiológicamente
# defendible.
#
# El monto reclamado se conserva como análisis de sensibilidad, dado que puede
# reflejar duración, complejidad, complicaciones o desenlace del episodio y no
# una característica basal del paciente.
#
# Los efectos aleatorios hospitalarios deben interpretarse como heterogeneidad
# residual ajustada y no como ranking directo de calidad asistencial.
