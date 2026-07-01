
# ==========================================================
# ANÁLISIS MULTINIVEL DE mortality_y_n
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

# Destinamos los outputs a un archivo .txt
sink("R_scripts/outputs/02_output.txt")

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
  mutate(
    sex = factor(
      sex,
      levels = c("male", "female"),
      labels = c("Masculino", "Femenino")
    ),
    
    hosp_type = factor(
      hosp_type,
      levels = c("c", "g"),
      labels = c("Privado", "Público")
    ),
    
    mortality_y_n = factor(
      mortality_y_n,
      levels = c("no", "yes")
    ),
    
    mortality = ifelse(mortality_y_n == "yes", 1, 0),
    
    across(
      c(
        hosp_name,
        hosp_district,
        caste_name,
        category_name,
        surgery,
        district_name
      ),
      as.factor
    )
  )

# ==========================================================
# 3. VARIABLES DERIVADAS Y TRANSFORMACIONES
# ==========================================================

datos <- datos %>%
  mutate(
#    dias_internacion = as.numeric(as_date(fecha_egreso) - as_date(fecha_cirugia)),
    claim_surgery_days = as.numeric(as_date(claim_date) - as_date(surgery_date)),
    amount_diff = claim_amount - preauth_amt,
    amount_diff_ratio = amount_diff / preauth_amt,
    
    age = ifelse(age < 0 | age > 120, NA, age),
#    dias_internacion = ifelse(dias_internacion < 0, NA, dias_internacion),
    claim_surgery_days = ifelse(claim_surgery_days < 0, NA, claim_surgery_days),
    amount_diff_ratio = ifelse(is.infinite(amount_diff_ratio), NA, amount_diff_ratio),
    
    log_claim_amount = log1p(claim_amount),
    log_amount_diff = log1p(abs(amount_diff))
  )

# ==========================================================
# 4. BASE ANALÍTICA INICIAL
# ==========================================================

datos_modelo_inicial <- datos %>%
  select(
    mortality,
    age,
    sex,
    caste_name,
    category_name,
    surgery,
    hosp_type,
    hosp_district,
    hosp_name,
    log_claim_amount,
    log_amount_diff
  ) %>%
  drop_na(
    mortality,
    age,
    sex,
    caste_name,
    category_name,
    surgery,
    hosp_type,
    hosp_district,
    hosp_name,
    log_claim_amount,
    log_amount_diff
  ) %>%
  mutate(
    age_z = as.numeric(scale(age)),
    log_claim_amount_z = as.numeric(scale(log_claim_amount)),
    log_amount_diff_z = as.numeric(scale(log_amount_diff))
  )


# ==========================================================
# 5. FILTRADO DE CATEGORÍAS DE INTERÉS (NEUROLOGÍA Y NEUROCIRUGÍA)
# ==========================================================

datos_modelo <- datos_modelo_inicial %>%
  filter(
    category_name %in% c("neurology", "neurosurgery")
  ) %>%
  droplevels()


# ==========================================================
# 6. FILTRADO POR HOSPITALES CON AL MENOS 10 PACIENTES
# ==========================================================

datos_modelo <- datos_modelo %>%
  group_by(hosp_name) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  droplevels()


# ==========================================================
# 7. MODELO LOGÍSTICO CLÁSICO DE REFERENCIA
# ==========================================================

glm_referencia <- glm(
  mortality ~
    age_z +
    sex +
    category_name +
    log_claim_amount_z +
    hosp_type,
  data = datos_modelo,
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
# 8. MODELO NULO MULTINIVEL
# ==========================================================

glmm_nulo_hospital <- glmer(
  mortality ~ 1 + (1 | hosp_name),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)

summary(glmm_nulo_hospital)

icc_nulo_hospital <- calcular_icc_logistico(glmm_nulo_hospital)

icc_nulo_hospital$varianzas
icc_nulo_hospital$ICC_total_pct

# ==========================================================
# 10. MODELOS MULTINIVEL
# ==========================================================

var_hosp_nulo <- icc_nulo_hospital$varianzas$vcov[
  icc_nulo_hospital$varianzas$grp == "hosp_name"
]

# ----------------------------------------------------------
# 10.1 Modelo demográfico
# ----------------------------------------------------------

glmm_demografico <- glmer(
  mortality ~
    age_z +
    sex +
    (1 | hosp_name),
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
  mortality ~
    age_z +
    sex +
    category_name +
    (1 | hosp_name),
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
  mortality ~
    age_z +
    sex +
    category_name +
    log_claim_amount_z +
    (1 | hosp_name),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

icc_paciente_principal <- calcular_icc_logistico(glmm_paciente_principal)

icc_paciente_principal$varianzas
icc_paciente_principal$ICC_total_pct

# ----------------------------------------------------------
# 10.4 Modelo contextual con tipo de hospital
# ----------------------------------------------------------

glmm_contextual_hospital <- glmer(
  mortality ~
    age_z +
    sex +
    category_name +
    log_claim_amount_z +
    hosp_type +
    (1 | hosp_name),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

icc_contextual_hospital <- calcular_icc_logistico(glmm_contextual_hospital)

icc_contextual_hospital$varianzas
icc_contextual_hospital$ICC_total_pct


# ----------------------------------------------------------
# 10.5 Modelo de sensibilidad con casta
# ----------------------------------------------------------

glmm_sensibilidad_casta <- glmer(
  mortality ~
    age_z +
    sex +
    caste_name +
    category_name +
    log_claim_amount_z +
    hosp_type +
    (1 | hosp_name),
  data = datos_modelo,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

icc_sensibilidad_casta <- calcular_icc_logistico(glmm_sensibilidad_casta)

icc_sensibilidad_casta$varianzas
icc_sensibilidad_casta$ICC_total_pct

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


# ==========================================================
# 12. COMPARACIÓN DE MODELOS EN TEST
# ==========================================================

set.seed(123)  # para reproducibilidad
n <- nrow(datos_modelo)

# índices de entrenamiento (80%)
train_index <- sample(1:n, size = 0.8 * n)

# crear los conjuntos
train <- datos_modelo[train_index, ]
test  <- datos_modelo[-train_index, ]


glmm_paciente_principal <- glmer(
  mortality ~
    age_z +
    sex +
    category_name +
    log_claim_amount_z +
    (1 | hosp_name),
  data = train,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

y_true <- test$mortality
y_pred_prob <- predict(glmm_paciente_principal, newdata = test, type = "response")

roc_curve <- roc(y_true, y_pred_prob)

roc_auc <- auc(roc_curve)

roc_auc

confusion_matrix <- table(
  Actual = y_true,
  Predicted = ifelse(y_pred_prob > 0.5, 1, 0)
)

confusion_matrix
