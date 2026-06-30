sink("output.txt")

# ==========================================================
# ANÁLISIS MULTINIVEL DE MORTALIDAD
# Pacientes atendidos en hospitales de India
# ==========================================================

# ==========================================================
# 0. PAQUETES
# ==========================================================

library(tidyverse)
library(janitor)
library(lubridate)
library(arrow)
library(lme4)
library(performance)
library(car)
library(pROC)

# ==========================================================
# 0.5. SET DIRECTORIO DEL TP
# ==========================================================

setwd("C:/Github/NTR-Arogyaseva-analysis")

# ==========================================================
# 1. FUNCIÓN PARA ICC LOGÍSTICO MULTINIVEL
# ==========================================================

calcular_icc_logistico <- function(modelo) {
  
  var_random <- as.data.frame(VarCorr(modelo)) %>%
    select(grp, vcov)
  
  var_total_random <- sum(var_random$vcov)
  
  tabla_icc <- var_random %>%
    mutate(
      ICC = vcov / (var_total_random + (pi^2 / 3)),
      ICC_pct = ICC * 100
    )
  
  list(
    varianzas = var_random,
    ICC_total = var_total_random / (var_total_random + (pi^2 / 3)),
    ICC_total_pct = 100 * var_total_random / (var_total_random + (pi^2 / 3)),
    ICC_por_nivel = tabla_icc
  )
}

# ==========================================================
# 2. CARGA Y LIMPIEZA INICIAL
# ==========================================================

datos <- read_parquet(
  "C:/Github/NTR-Arogyaseva-analysis/data/clean/ntrarogyaseva.parquet"
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
    
    mortalidad_binaria = ifelse(mortalidad == "yes", 1, 0),
    
    across(
      c(
        nombre_hospital,
        distrito_hospital,
        nombre_casta,
        nombre_categoria,
        nombre_cirugia
      ),
      as.factor
    )
  )

# ==========================================================
# 3. VARIABLES DERIVADAS Y TRANSFORMACIONES
# ==========================================================

datos <- datos %>%
  mutate(
    dias_internacion = as.numeric(as_date(fecha_egreso) - as_date(fecha_cirugia)),
    dias_reclamo_cirugia = as.numeric(as_date(fecha_reclamo) - as_date(fecha_cirugia)),
    diferencia_montos = monto_reclamado - monto_preautorizado,
    ratio_montos = monto_reclamado / monto_preautorizado,
    
    edad = ifelse(edad < 0 | edad > 120, NA, edad),
    dias_internacion = ifelse(dias_internacion < 0, NA, dias_internacion),
#    dias_reclamo_cirugia = ifelse(dias_reclamo_cirugia < 0, NA, dias_reclamo_cirugia),
    ratio_montos = ifelse(is.infinite(ratio_montos), NA, ratio_montos),
    
    log_monto_reclamado = log1p(monto_reclamado),
    log_diferencia_montos = log1p(abs(diferencia_montos))
  )

# Interpretación:
# Las transformaciones logarítmicas reducen la asimetría positiva de las
# variables económicas. log_diferencia_montos se conserva para descriptivos
# y sensibilidad, pero no se usará en el modelo principal porque ocurre
# durante/después del proceso asistencial.

# ==========================================================
# 4. BASE ANALÍTICA INICIAL
# ==========================================================

datos_modelo_inicial <- datos %>%
  select(
    mortalidad_binaria,
    edad,
    sexo,
    nombre_casta,
    nombre_categoria,
    nombre_cirugia,
    tipo_hospital,
    distrito_hospital,
    nombre_hospital,
    log_monto_reclamado,
    log_diferencia_montos
  ) %>%
  drop_na(
    mortalidad_binaria,
    edad,
    sexo,
    nombre_casta,
    nombre_categoria,
    nombre_cirugia,
    tipo_hospital,
    distrito_hospital,
    nombre_hospital,
    log_monto_reclamado,
    log_diferencia_montos
  ) %>%
  mutate(
    edad_z = as.numeric(scale(edad)),
    log_monto_reclamado_z = as.numeric(scale(log_monto_reclamado)),
    log_diferencia_montos_z = as.numeric(scale(log_diferencia_montos))
  )

# ==========================================================
# 5. AGRUPAMIENTO DE CATEGORÍAS CLÍNICAS INESTABLES
# ==========================================================

resumen_categoria <- datos_modelo_inicial %>%
  group_by(nombre_categoria) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    .groups = "drop"
  ) %>%
  arrange(pacientes)

categorias_inestables <- resumen_categoria %>%
  filter(pacientes < 100 | muertes < 50) %>%
  pull(nombre_categoria)

datos_modelo <- datos_modelo_inicial %>%
  mutate(
    nombre_categoria_modelo = if_else(
      nombre_categoria %in% categorias_inestables,
      "Other",
      as.character(nombre_categoria)
    ),
    nombre_categoria_modelo = factor(nombre_categoria_modelo)
  )

resumen_categoria_modelo <- datos_modelo %>%
  group_by(nombre_categoria_modelo) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    .groups = "drop"
  ) %>%
  arrange(pacientes)

resumen_categoria_modelo

# Interpretación:
# nombre_categoria_modelo será la variable clínica usada en los modelos.
# La variable original nombre_categoria queda solo para descriptivos.

# ==========================================================
# 6. DEFINICIÓN DE HOSPITALES ESTABLES
# ==========================================================

resumen_hospital_todos <- datos_modelo %>%
  group_by(nombre_hospital) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    .groups = "drop"
  )

hospitales_estables <- resumen_hospital_todos %>%
  filter(pacientes >= 40, muertes >= 10) %>%
  pull(nombre_hospital)

datos_modelo_hosp_estables <- datos_modelo %>%
  filter(nombre_hospital %in% hospitales_estables) %>%
  droplevels()

resumen_base_estable <- datos_modelo_hosp_estables %>%
  summarise(
    pacientes = n(),
    hospitales = n_distinct(nombre_hospital),
    distritos = n_distinct(distrito_hospital),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100
  )

resumen_base_estable

# Interpretación:
# Esta será la base principal del análisis multinivel.
# Se prioriza estabilidad de estimaciones hospitalarias y conservación
# de eventos de mortalidad.

# ==========================================================
# 7. ESTRUCTURA JERÁRQUICA DESCRIPTIVA
# ==========================================================

estructura_distrito_estables <- datos_modelo_hosp_estables %>%
  group_by(distrito_hospital) %>%
  summarise(
    pacientes = n(),
    hospitales = n_distinct(nombre_hospital),
    muertes = sum(mortalidad_binaria),
    mortalidad_pct = mean(mortalidad_binaria) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mortalidad_pct))

estructura_distrito_estables

resumen_hospital_estables <- datos_modelo_hosp_estables %>%
  group_by(nombre_hospital) %>%
  summarise(
    pacientes = n(),
    muertes = sum(mortalidad_binaria),
    mortalidad_cruda_pct = mean(mortalidad_binaria) * 100,
    tipo_hospital = first(tipo_hospital),
    distrito_hospital = first(distrito_hospital),
    .groups = "drop"
  )

summary(resumen_hospital_estables$pacientes)
summary(resumen_hospital_estables$muertes)
summary(resumen_hospital_estables$mortalidad_cruda_pct)

# Interpretación:
# Las tasas crudas por hospital no deben usarse para identificar desempeño.
# La identificación de hospitales atípicos se hará con efectos aleatorios
# ajustados.

# ==========================================================
# INTERPRETACIÓN DE LOS BLOQUES 0 A 7
# ==========================================================
#
# La carga, limpieza y preparación inicial de la base se ejecutaron sin
# errores. La variable de resultado quedó codificada como mortalidad_binaria
# en formato 0/1, y las variables categóricas principales fueron definidas
# como factores, lo cual es adecuado para los modelos logísticos y
# multinivel posteriores.
#
# Se construyeron variables económicas transformadas mediante log1p().
# Esta decisión es apropiada porque los montos sanitarios suelen presentar
# fuerte asimetría positiva. La variable log_monto_reclamado_z se conserva
# como predictor principal, mientras que log_diferencia_montos_z se conserva
# solo para análisis descriptivos o de sensibilidad, dado que representa
# información producida durante o después del episodio asistencial.
#
# Luego de aplicar los criterios de inclusión y eliminar casos con valores
# perdidos en variables centrales, se conformó la base datos_modelo para
# el análisis.
#
# La agrupación de categorías clínicas inestables redujo la variable
# nombre_categoria a 18 categorías modelables mediante
# nombre_categoria_modelo. Esta decisión es metodológicamente adecuada
# porque evita categorías con muy pocos pacientes o pocos eventos de
# mortalidad, que podrían generar separación completa, errores estándar
# elevados o estimaciones inestables.
#
# La categoría "Other" concentra 88.886 pacientes y 94 muertes, con una
# mortalidad cruda baja. Esta categoría debe interpretarse como una
# categoría técnica de estabilización del modelo y no como una entidad
# clínica homogénea.
#
# La base principal restringida a hospitales estables quedó conformada por
# 340.736 pacientes, 120 hospitales, 14 distritos y 9.751 muertes. La
# mortalidad cruda global en esta base fue de 2,86%.
#
# La restricción a hospitales estables es adecuada para el análisis
# principal porque conserva un número muy elevado de pacientes y eventos,
# y reduce la influencia de hospitales con bajo volumen asistencial o pocos
# eventos, cuyas tasas de mortalidad serían altamente inestables.
#
# La estructura jerárquica descriptiva confirma que los pacientes están
# agrupados en hospitales y que estos hospitales se distribuyen en
# distritos. Los distritos presentan entre 2 y 17 hospitales, con
# mortalidad cruda entre 1,53% y 4,30%.
#
# A nivel hospitalario, los hospitales estables tienen entre 198 y 21.854
# pacientes, con una mediana de 1.820 pacientes. El número de muertes por
# hospital varía entre 10 y 770, con mediana de 37,5 eventos.
#
# La mortalidad cruda hospitalaria presenta variabilidad relevante, con
# valores entre 0,22% y 11,00%. Sin embargo, estas tasas no deben
# interpretarse directamente como desempeño hospitalario porque no están
# ajustadas por edad, sexo, categoría clínica ni complejidad del caso.
#
# Estos resultados justifican avanzar con modelos multinivel. El siguiente
# paso metodológico es estimar modelos nulos para cuantificar formalmente
# la proporción de variabilidad atribuible al hospital y evaluar si el
# distrito aporta variabilidad adicional.

# ==========================================================
# 8. MODELO LOGÍSTICO CLÁSICO DE REFERENCIA
# ==========================================================

glm_referencia <- glm(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_categoria_modelo +
    log_monto_reclamado_z +
    tipo_hospital,
  data = datos_modelo_hosp_estables,
  family = binomial()
)

summary(glm_referencia)

or_glm_referencia <- exp(cbind(
  OR = coef(glm_referencia),
  confint.default(glm_referencia)
))

or_glm_referencia

car::vif(glm_referencia)

# ==========================================================
# INTERPRETACIÓN DEL MODELO LOGÍSTICO DE REFERENCIA
# ==========================================================
#
# El modelo logístico clásico se estimó como referencia para identificar
# los principales factores asociados con la mortalidad antes de considerar
# la estructura jerárquica de los datos.
#
# La edad mostró una asociación positiva con la mortalidad. Por cada
# incremento de un desvío estándar en la edad, las odds de fallecimiento
# aumentaron aproximadamente un 29% (OR = 1.29; IC95%: 1.25–1.33).
#
# El monto reclamado presentó una asociación inversa con la mortalidad
# (OR = 0.39; IC95%: 0.38–0.40), indicando que los episodios con mayores
# montos reclamados tendieron a presentar menor probabilidad de muerte.
# Dado que esta variable se genera durante el proceso asistencial, su
# interpretación debe entenderse como un factor asociado al episodio de
# atención y no como un predictor basal del paciente.
#
# La categoría clínica fue el predictor con mayor capacidad explicativa,
# observándose importantes diferencias en el riesgo de mortalidad entre
# especialidades. Esto confirma que el case-mix clínico constituye un
# determinante fundamental de la mortalidad y justifica su incorporación
# en los modelos multinivel posteriores.
#
# El tipo de hospital también mostró una asociación estadísticamente
# significativa. En este modelo, los hospitales públicos presentaron un
# 18% mayor odds de mortalidad que los privados (OR = 1.18; IC95%:
# 1.13–1.24). Sin embargo, este efecto aún no considera la dependencia
# entre pacientes atendidos en un mismo hospital.
#
# El sexo femenino presentó una asociación estadísticamente significativa
# (OR = 1.05), aunque con una magnitud de efecto muy pequeña. Debido al
# gran tamaño muestral, esta diferencia probablemente refleja un efecto
# de escasa relevancia clínica. La variable se conservará en los modelos
# posteriores por su importancia epidemiológica como variable de ajuste.
#
# Los factores de inflación de la varianza fueron bajos
# (GVIF^(1/(2Df)) < 1.5 para todas las variables), indicando ausencia de
# problemas relevantes de colinealidad entre los predictores incluidos.
#
# En conjunto, este modelo confirma la importancia de las características
# individuales y clínicas para explicar la mortalidad. Sin embargo, como
# supone independencia entre observaciones, sus coeficientes aún pueden
# estar sesgados por la agrupación de pacientes dentro de hospitales.
# Por ello, el siguiente paso consiste en estimar modelos multinivel para
# cuantificar formalmente el efecto hospitalario mediante el cálculo del
# coeficiente de correlación intraclase (ICC).

# ==========================================================
# 9. MODELOS NULOS MULTINIVEL
# ==========================================================

# ----------------------------------------------------------
# 9.1 Modelo nulo principal: paciente dentro de hospital
# ----------------------------------------------------------

glmm_nulo_hospital <- glmer(
  mortalidad_binaria ~ 1 + (1 | nombre_hospital),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

summary(glmm_nulo_hospital)

icc_nulo_hospital <- calcular_icc_logistico(glmm_nulo_hospital)

icc_nulo_hospital$varianzas
icc_nulo_hospital$ICC_total_pct

# Interpretación:
# Este ICC responde qué proporción de la variabilidad latente de mortalidad
# se atribuye al hospital antes de ajustar por características del paciente.

# ----------------------------------------------------------
# 9.2 Modelo nulo de tres niveles: distrito / hospital
# ----------------------------------------------------------

glmm_nulo_distrito_hospital <- glmer(
  mortalidad_binaria ~ 1 + (1 | distrito_hospital / nombre_hospital),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

summary(glmm_nulo_distrito_hospital)

icc_nulo_distrito_hospital <- calcular_icc_logistico(
  glmm_nulo_distrito_hospital
)

icc_nulo_distrito_hospital$ICC_por_nivel

anova(glmm_nulo_hospital, glmm_nulo_distrito_hospital)

AIC(glmm_nulo_hospital, glmm_nulo_distrito_hospital)
BIC(glmm_nulo_hospital, glmm_nulo_distrito_hospital)

# Interpretación:
# Si la varianza del distrito es cero o el modelo no mejora el ajuste,
# se conserva el modelo de dos niveles por parsimonia.
# ==========================================================
# INTERPRETACIÓN DE LOS MODELOS NULOS MULTINIVEL
# ==========================================================
#
# El modelo nulo de dos niveles permitió cuantificar la variabilidad de la
# mortalidad atribuible al hospital antes de incorporar cualquier variable
# explicativa.
#
# La varianza del intercepto aleatorio fue de 0.534, lo que corresponde a
# un ICC de 13.96%. Esto indica que aproximadamente el 14% de la
# variabilidad latente de la mortalidad se debe a diferencias entre
# hospitales, mientras que el 86% restante corresponde a diferencias entre
# pacientes atendidos dentro de los hospitales.
#
# Este resultado demuestra que los pacientes atendidos en un mismo hospital
# presentan probabilidades de fallecimiento más similares entre sí que
# pacientes atendidos en hospitales diferentes, justificando el uso de un
# modelo multinivel en lugar de una regresión logística convencional.
#
# Posteriormente se evaluó un modelo de tres niveles incorporando el
# distrito hospitalario como un nivel jerárquico adicional.
#
# La varianza estimada para el distrito fue exactamente igual a cero
# (ICC = 0%), indicando que, una vez considerado el hospital, no existe
# variabilidad adicional atribuible al distrito.
#
# Esta conclusión fue consistente con la comparación formal entre modelos.
# La prueba de razón de verosimilitud no mostró mejoras al incorporar el
# nivel distrito (χ² = 0.00; gl = 1; p = 0.9999), mientras que tanto el
# AIC como el BIC fueron ligeramente mayores para el modelo de tres
# niveles.
#
# En consecuencia, no existe evidencia empírica que justifique mantener
# un tercer nivel jerárquico. El modelo multinivel de dos niveles
# (pacientes anidados en hospitales) representa la estructura más
# parsimoniosa y será utilizado en el resto del análisis.
#
# Preguntas de investigación respondidas:
#
# ✔ ¿Los pacientes atendidos en un mismo hospital presentan probabilidades
#   de fallecimiento más similares entre sí?
#   Sí. Aproximadamente el 14% de la variabilidad de la mortalidad se
#   atribuye al hospital.
#
# ✔ ¿Existe efecto territorial del distrito?
#   No. El distrito no explicó variabilidad adicional una vez considerado
#   el hospital.
#
# ✔ ¿Es necesario utilizar modelos multinivel?
#   Sí. El ICC hospitalario obtenido justifica modelar explícitamente la
#   dependencia entre pacientes atendidos en un mismo hospital mediante un
#   modelo multinivel de dos niveles.

# ==========================================================
# 10. MODELOS MULTINIVEL EXPLICATIVOS
# ==========================================================

var_hosp_nulo <- icc_nulo_hospital$varianzas$vcov[
  icc_nulo_hospital$varianzas$grp == "nombre_hospital"
]

# ----------------------------------------------------------
# 10.1 Modelo demográfico
# ----------------------------------------------------------

glmm_demografico <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    (1 | nombre_hospital),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

icc_demografico <- calcular_icc_logistico(glmm_demografico)

icc_demografico$varianzas
icc_demografico$ICC_total_pct


# ----------------------------------------------------------
# 10.2 Modelo clínico del paciente
# ----------------------------------------------------------

glmm_paciente_clinico <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_categoria_modelo +
    (1 | nombre_hospital),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

icc_paciente_clinico <- calcular_icc_logistico(glmm_paciente_clinico)

icc_paciente_clinico$varianzas
icc_paciente_clinico$ICC_total_pct

# ----------------------------------------------------------
# 10.3 Modelo paciente completo principal
# ----------------------------------------------------------

glmm_paciente_principal <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_categoria_modelo +
    log_monto_reclamado_z +
    (1 | nombre_hospital),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(glmm_paciente_principal)

or_glmm_paciente_principal <- exp(cbind(
  OR = fixef(glmm_paciente_principal),
  confint(glmm_paciente_principal, parm = "beta_", method = "Wald")
))

or_glmm_paciente_principal

icc_paciente_principal <- calcular_icc_logistico(glmm_paciente_principal)

# ----------------------------------------------------------
# 10.4 Modelo contextual con tipo de hospital
# ----------------------------------------------------------

glmm_contextual_hospital <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_categoria_modelo +
    log_monto_reclamado_z +
    tipo_hospital +
    (1 | nombre_hospital),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(glmm_contextual_hospital)

or_glmm_contextual_hospital <- exp(cbind(
  OR = fixef(glmm_contextual_hospital),
  confint(glmm_contextual_hospital, parm = "beta_", method = "Wald")
))

or_glmm_contextual_hospital

icc_contextual_hospital <- calcular_icc_logistico(glmm_contextual_hospital)

# ----------------------------------------------------------
# 10.5 Modelo de sensibilidad con casta
# ----------------------------------------------------------

glmm_sensibilidad_casta <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_casta +
    nombre_categoria_modelo +
    log_monto_reclamado_z +
    tipo_hospital +
    (1 | nombre_hospital),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(glmm_sensibilidad_casta)

icc_sensibilidad_casta <- calcular_icc_logistico(glmm_sensibilidad_casta)

# Interpretación:
# Casta no forma parte del modelo principal por parsimonia y bajo aporte
# previo, pero se conserva como análisis de sensibilidad por relevancia
# epidemiológica en desigualdades sociales.

# ==========================================================
# 11. COMPARACIÓN DE MODELOS PRINCIPALES
# ==========================================================

comparacion_modelos_principales <- tibble(
  modelo = c(
    "Nulo hospital",
    "Demográfico",
    "Paciente clínico",
    "Paciente principal",
    "Paciente + tipo hospital",
    "Sensibilidad con casta"
  ),
  objeto = c(
    "glmm_nulo_hospital",
    "glmm_demografico",
    "glmm_paciente_clinico",
    "glmm_paciente_principal",
    "glmm_contextual_hospital",
    "glmm_sensibilidad_casta"
  ),
  AIC = c(
    AIC(glmm_nulo_hospital),
    AIC(glmm_demografico),
    AIC(glmm_paciente_clinico),
    AIC(glmm_paciente_principal),
    AIC(glmm_contextual_hospital),
    AIC(glmm_sensibilidad_casta)
  ),
  BIC = c(
    BIC(glmm_nulo_hospital),
    BIC(glmm_demografico),
    BIC(glmm_paciente_clinico),
    BIC(glmm_paciente_principal),
    BIC(glmm_contextual_hospital),
    BIC(glmm_sensibilidad_casta)
  ),
  varianza_hospital = c(
    icc_nulo_hospital$varianzas$vcov[1],
    icc_demografico$varianzas$vcov[1],
    icc_paciente_clinico$varianzas$vcov[1],
    icc_paciente_principal$varianzas$vcov[1],
    icc_contextual_hospital$varianzas$vcov[1],
    icc_sensibilidad_casta$varianzas$vcov[1]
  ),
  ICC_pct = c(
    icc_nulo_hospital$ICC_total_pct,
    icc_demografico$ICC_total_pct,
    icc_paciente_clinico$ICC_total_pct,
    icc_paciente_principal$ICC_total_pct,
    icc_contextual_hospital$ICC_total_pct,
    icc_sensibilidad_casta$ICC_total_pct
  )
) %>%
  mutate(
    PCV_hospital_pct = 100 *
      (var_hosp_nulo - varianza_hospital) / var_hosp_nulo
  )

comparacion_modelos_principales

anova(
  glmm_nulo_hospital,
  glmm_demografico,
  glmm_paciente_clinico,
  glmm_paciente_principal,
  glmm_contextual_hospital
)

anova(
  glmm_contextual_hospital,
  glmm_sensibilidad_casta
)

# Interpretación:
# El modelo final no se elige solo por AIC.
# Se evalúa simultáneamente:
# - mejora del ajuste,
# - reducción de la varianza hospitalaria,
# - parsimonia,
# - coherencia epidemiológica,
# - estabilidad del modelo.

# ==========================================================
# 11.5. COMPARACIÓN DE MODELOS POSTERIORES
# ==========================================================

# utilizamos datos iniciales, sin filtrado de categorías ni hospitales
glmm_paciente_clinico_nuevo <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    (1 | nombre_categoria) +
    (1 | nombre_hospital),
  data = datos_modelo_inicial,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

icc_paciente_clinico_nuevo <- calcular_icc_logistico(glmm_paciente_clinico_nuevo)

icc_paciente_clinico_nuevo$varianzas
icc_paciente_clinico_nuevo$ICC_por_nivel


glmm_paciente_clinico_otro <- glmer(
  mortalidad_binaria ~
    edad_z +
    sexo +
    nombre_categoria +
    (1 | nombre_hospital),
  data = datos_modelo_inicial,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

icc_paciente_clinico_otro <- calcular_icc_logistico(glmm_paciente_clinico_otro)

icc_paciente_clinico_otro$varianzas
icc_paciente_clinico_otro$ICC_total_pct

anova(
  glmm_paciente_clinico_otro,
  glmm_paciente_clinico_nuevo
)

comparacion_modelos_posteriores <- tibble(
  modelo = c(
    "Paciente clínico (efecto aleatorio)",
    "Paciente clínico (efecto fijo)"
  ),
  objeto = c(
    "glmm_paciente_clinico_nuevo",
    "glmm_paciente_clinico_otro"
  ),
  AIC = c(
    AIC(glmm_paciente_clinico_nuevo),
    AIC(glmm_paciente_clinico_otro)
  ),
  BIC = c(
    BIC(glmm_paciente_clinico_nuevo),
    BIC(glmm_paciente_clinico_otro)
  ),
  varianza_hospital = c(
    icc_paciente_clinico_nuevo$varianzas$vcov[1],
    icc_paciente_clinico_otro$varianzas$vcov[1]
  ),
  ICC_pct = c(
    icc_paciente_clinico_nuevo$ICC_total_pct,
    icc_paciente_clinico_otro$ICC_total_pct
  )
) %>%
  mutate(
    PCV_hospital_pct = 100 *
      (var_hosp_nulo - varianza_hospital) / var_hosp_nulo
  )

comparacion_modelos_posteriores

# ==========================================================
# INTERPRETACIÓN DE LOS MODELOS MULTINIVEL EXPLICATIVOS
# ==========================================================
#
# Los modelos multinivel se construyeron de forma secuencial para evaluar
# el aporte incremental de las características individuales y contextuales
# en la explicación de la mortalidad, así como su capacidad para explicar
# las diferencias observadas entre hospitales.
#
# La incorporación de las variables demográficas (edad y sexo) mejoró
# significativamente el ajuste respecto del modelo nulo
# (LRT = 429.10; p < 0.001). Sin embargo, la varianza hospitalaria aumentó
# levemente y el ICC pasó de 13.96% a 15.56%.
#
# Este incremento no debe interpretarse como un empeoramiento del modelo.
# En modelos logísticos multinivel el ICC puede aumentar tras incorporar
# covariables debido a cambios en la escala latente. Por este motivo, la
# reducción de la varianza hospitalaria (PCV) constituye un indicador más
# apropiado para evaluar cuánto explican las covariables sobre las
# diferencias entre hospitales.
#
# La incorporación de la categoría clínica produjo la mayor mejora del
# ajuste de todo el proceso de modelización (LRT = 6080.44; p < 0.001).
# Además, redujo la varianza hospitalaria de 0.606 a 0.283, lo que indica
# que una proporción importante de las diferencias entre hospitales se
# explica por el perfil clínico de los pacientes atendidos (case-mix).
#
# La incorporación del monto reclamado mejoró nuevamente el ajuste del
# modelo (LRT = 3131.66; p < 0.001). Sin embargo, la varianza hospitalaria
# aumentó moderadamente (0.329), lo que sugiere que esta variable explica
# principalmente diferencias individuales de mortalidad y no diferencias
# sistemáticas entre hospitales.
#
# Posteriormente se evaluó si el tipo de hospital aportaba información
# adicional una vez ajustado el perfil de los pacientes. El efecto del tipo
# de hospital dejó de ser estadísticamente significativo (OR = 1.17;
# IC95%: 0.87–1.57; p = 0.295) y la comparación entre modelos no mostró una
# mejora significativa del ajuste (LRT = 1.09; p = 0.296).
#
# Esto indica que las diferencias inicialmente observadas entre hospitales
# públicos y privados en el modelo logístico clásico pueden explicarse por
# diferencias en las características de los pacientes atendidos y por la
# variabilidad propia de cada hospital, más que por el tipo de gestión del
# establecimiento.
#
# Finalmente, la incorporación de la variable casta tampoco produjo una
# mejora significativa del ajuste (LRT = 4.89; gl = 5; p = 0.429), ni
# modificó la varianza hospitalaria. Por lo tanto, se decidió conservarla
# únicamente como análisis de sensibilidad y no como parte del modelo
# principal, siguiendo el principio de parsimonia.
#
# En conjunto, los resultados muestran que la mayor parte de la capacidad
# explicativa proviene de las características clínicas de los pacientes,
# mientras que las variables contextuales evaluadas (tipo de hospital y
# casta) no aportaron información adicional una vez controlado el case-mix.
#
# Preguntas de investigación respondidas:
#
# ✔ ¿Qué variables individuales aumentan el riesgo de mortalidad?
#   La edad, la categoría clínica y el monto reclamado mostraron asociación
#   independiente con la mortalidad. El sexo no presentó asociación
#   estadísticamente significativa luego del ajuste multinivel.
#
# ✔ ¿Las diferencias entre hospitales permanecen luego del ajuste?
#   Sí. Aunque la varianza hospitalaria disminuyó respecto del modelo nulo,
#   persiste heterogeneidad residual entre hospitales, lo que indica que
#   continúan existiendo diferencias no explicadas por las variables
#   individuales incluidas en el modelo.
#
# ✔ ¿El tipo de hospital explica parte de la variabilidad hospitalaria?
#   No. Una vez ajustadas las características de los pacientes, el tipo de
#   hospital no mostró una asociación independiente con la mortalidad ni
#   mejoró el ajuste del modelo.
#
# ✔ ¿La casta mejora el modelo?
#   No. Su incorporación no produjo mejoras significativas y se mantuvo
#   únicamente como análisis de sensibilidad.
# ==========================================================
# 12. SELECCIÓN DEL MODELO FINAL
# ==========================================================

# Objetivo:
# Seleccionar el modelo que mejor representa la mortalidad ajustada
# considerando simultáneamente ajuste, parsimonia y coherencia
# epidemiológica.

# Criterios utilizados:
# - Prueba de razón de verosimilitud (LRT)
# - AIC y BIC
# - Reducción de la varianza hospitalaria
# - Principio de parsimonia

# ----------------------------------------------------------
# Modelo seleccionado
# ----------------------------------------------------------

modelo_final <- glmm_paciente_principal

icc_modelo_final <- icc_paciente_principal

or_modelo_final <- or_glmm_paciente_principal

summary(modelo_final)

or_modelo_final

icc_modelo_final$ICC_total_pct

icc_modelo_final$ICC_por_nivel

# ==========================================================
# INTERPRETACIÓN DEL MODELO MULTINIVEL FINAL
# ==========================================================
#
# El modelo paciente principal fue seleccionado como modelo final del
# estudio por presentar el mejor equilibrio entre ajuste, parsimonia y
# capacidad explicativa.
#
# La incorporación del tipo de hospital no mejoró significativamente el
# ajuste del modelo, no redujo la variabilidad hospitalaria residual y
# aumentó ligeramente los criterios de información (AIC y BIC). En
# consecuencia, el tipo de hospital no fue retenido en el modelo final y
# se informó únicamente como un análisis contextual complementario.
#
# Luego del ajuste simultáneo por todas las covariables incluidas, la edad
# continuó asociándose de forma independiente con la mortalidad. Por cada
# incremento de un desvío estándar en la edad, las odds de fallecimiento
# aumentaron aproximadamente un 34% (OR = 1.34; IC95%: 1.30–1.38).
#
# La categoría clínica representó el principal determinante individual de
# la mortalidad, observándose importantes diferencias de riesgo entre las
# distintas especialidades incluso luego del ajuste por el resto de las
# variables.
#
# El monto reclamado mantuvo una asociación inversa con la mortalidad
# (OR = 0.39; IC95%: 0.37–0.40). Dado que esta variable se genera durante
# el proceso asistencial, debe interpretarse como un marcador del episodio
# de atención y no como un predictor basal del paciente.
#
# El sexo no mostró una asociación independiente con la mortalidad
# (p = 0.131), indicando que las diferencias inicialmente observadas se
# explican por el perfil clínico y las demás covariables incorporadas al
# modelo.
#
# A pesar del ajuste, la varianza hospitalaria permaneció distinta de
# cero (σ² = 0.329), con un ICC residual de 9.08%. Esto indica que
# aproximadamente el 9% de la variabilidad latente residual de la
# mortalidad continúa siendo atribuible al hospital donde fue atendido el
# paciente.
#
# En comparación con el modelo nulo, el ICC disminuyó desde 13.96% hasta
# 9.08%, lo que demuestra que una parte importante de las diferencias
# entre hospitales se explica por las características de los pacientes.
# Sin embargo, persiste heterogeneidad hospitalaria residual que no puede
# atribuirse únicamente a la edad, el sexo, la categoría clínica o el
# monto reclamado.
#
# Estos resultados sugieren la existencia de factores hospitalarios no
# medidos —como diferencias en la organización asistencial, disponibilidad
# de recursos, experiencia de los equipos o calidad de la atención— que
# podrían contribuir a explicar la variabilidad residual observada entre
# hospitales.
#
# Preguntas de investigación respondidas:
#
# ✔ ¿Qué variables individuales explican la mortalidad?
#   Principalmente la edad y el perfil clínico del paciente, mientras que
#   el sexo no mostró un efecto independiente luego del ajuste.
#
# ✔ ¿Las diferencias entre hospitales permanecen luego del ajuste?
#   Sí. Aunque disminuyen considerablemente, persiste un efecto
#   hospitalario residual que justifica la estimación de efectos aleatorios
#   para identificar hospitales con mortalidad ajustada mayor o menor que
#   la esperada.

#insight:Las variables individuales incorporadas en el modelo explicaron aproximadamente el 38% de la variabilidad entre hospitales observada en el modelo nulo. 

# ==========================================================
# 13. EFECTOS ALEATORIOS HOSPITALARIOS
# Mortalidad ajustada superior o inferior a la esperada
# ==========================================================

# Objetivo:
# Identificar hospitales con mortalidad ajustada superior o inferior a la
# esperada, utilizando los efectos aleatorios del modelo final.

# Nota metodológica:
# Los efectos aleatorios hospitalarios no son tasas crudas de mortalidad.
# Representan desviaciones ajustadas respecto del promedio hospitalario,
# luego de controlar por las covariables incluidas en el modelo final.

# ----------------------------------------------------------
# 13.1 Extraer efectos aleatorios hospitalarios
# ----------------------------------------------------------

efectos_hospital <- ranef(
  modelo_final,
  condVar = TRUE
)$nombre_hospital %>%
  rownames_to_column("nombre_hospital") %>%
  rename(efecto_hospital = `(Intercept)`) %>%
  left_join(
    resumen_hospital_estables,
    by = "nombre_hospital"
  ) %>%
  arrange(desc(efecto_hospital))

# Interpretación:
# efecto_hospital > 0 indica mortalidad ajustada superior a la esperada.
# efecto_hospital < 0 indica mortalidad ajustada inferior a la esperada.

# ----------------------------------------------------------
# 13.2 Hospitales con mortalidad ajustada superior a la esperada
# ----------------------------------------------------------

hospitales_sobre_esperado <- efectos_hospital %>%
  slice_max(
    efecto_hospital,
    n = 20
  )

hospitales_sobre_esperado

# Interpretación:
# Estos hospitales presentan los mayores efectos aleatorios positivos.
# Esto sugiere una mortalidad ajustada superior a la esperada según el
# perfil de pacientes atendidos.
# Deben interpretarse como candidatos para auditoría o análisis posterior,
# no como evidencia directa de mala calidad asistencial.

# ----------------------------------------------------------
# 13.3 Hospitales con mortalidad ajustada inferior a la esperada
# ----------------------------------------------------------

hospitales_bajo_esperado <- efectos_hospital %>%
  slice_min(
    efecto_hospital,
    n = 20
  )

hospitales_bajo_esperado

# Interpretación:
# Estos hospitales presentan los efectos aleatorios más negativos.
# Esto sugiere una mortalidad ajustada inferior a la esperada según el
# perfil de pacientes atendidos.
# Pueden orientar estudios de buenas prácticas, pero no prueban por sí
# solos mejor calidad asistencial.

# ==========================================================
# INTERPRETACIÓN GENERAL DEL BLOQUE
# ==========================================================
#
# A partir del modelo multinivel final se estimaron los efectos aleatorios
# hospitalarios, también llamados BLUPs.
#
# Estos efectos representan cuánto se desvía cada hospital de la mortalidad
# esperada promedio, luego de ajustar por las características incluidas en
# el modelo final.
#
# Los valores positivos indican mortalidad ajustada superior a la esperada,
# mientras que los valores negativos indican mortalidad ajustada inferior
# a la esperada.
#
# Estos resultados no deben interpretarse como un ranking directo de
# calidad hospitalaria. Los BLUPs permiten identificar hospitales que
# requieren análisis adicional, pero no prueban por sí solos buena o mala
# calidad asistencial.
#
# La ventaja de utilizar efectos aleatorios en lugar de tasas crudas es que
# el modelo incorpora contracción estadística (shrinkage), reduciendo la
# influencia de hospitales con menor volumen de pacientes o menor cantidad
# de eventos.
#
# En consecuencia, este bloque permite identificar hospitales candidatos
# para auditorías clínicas, revisión de procesos asistenciales o análisis
# institucionales más detallados.



# ==========================================================
# 14. AUC DEL MODELO FINAL
# ==========================================================

prob_modelo_final <- predict(
  modelo_final,
  type = "response"
)

roc_modelo_final <- roc(
  datos_modelo_hosp_estables$mortalidad_binaria,
  prob_modelo_final,
  quiet = TRUE
)

auc(roc_modelo_final)
plot(roc_modelo_final)

sink("log.txt", split = TRUE)

# ==========================================================
# INTERPRETACIÓN DEL DESEMPEÑO DISCRIMINATIVO
# ==========================================================
#
# La curva ROC evalúa la capacidad discriminativa del modelo final, es
# decir, su habilidad para distinguir entre pacientes que fallecieron y
# pacientes que sobrevivieron.
#
# El modelo obtuvo un AUC de 0.806, lo que indica una buena capacidad de
# discriminación. En promedio, existe aproximadamente un 81% de
# probabilidad de que el modelo asigne un mayor riesgo a un paciente que
# falleció que a otro que sobrevivió.
#
# Este resultado confirma que las variables incorporadas en el modelo
# permiten discriminar adecuadamente entre pacientes con distinto riesgo
# de mortalidad.
#
# No obstante, el AUC constituye una medida complementaria del desempeño
# predictivo y no reemplaza la interpretación epidemiológica del modelo
# multinivel.
#
# El objetivo principal del estudio fue explicar la variabilidad de la
# mortalidad entre hospitales. Para ello, los indicadores fundamentales
# continúan siendo el coeficiente de correlación intraclase (ICC), la
# reducción de la varianza hospitalaria y los efectos aleatorios
# hospitalarios (BLUPs).
#
# En conjunto, los resultados muestran que el modelo presenta una buena
# capacidad discriminativa y que, aun después del ajuste por las
# características individuales de los pacientes, persisten diferencias
# sistemáticas entre hospitales que justifican futuras investigaciones
# sobre factores organizacionales y de calidad asistencial.



