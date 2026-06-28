# load required libraries 
library(readr)
library(glmmTMB)
library(broom.mixed)
library(dplyr)

# load in data files

model_df <- read_csv("model_fits_per_segment.csv")
determ_df <- read_csv("Deterministic-Combined-Processed.csv")
prob_df <- read_csv("Probabilistic-Combined-Processed.csv")

## DRL TRAINING
model_df_determ_train <- model_df[
  model_df$session_type == "MOUSE DRL" &
    model_df$session_id %in% 1:10 &
    model_df$block_prob == "1",
]
determ_df_train <- determ_df[
  determ_df$Task_Stage == "DRL Training" &
    determ_df$Session_Number %in% 1:10,
]

model_df_determ_train$genotype <- relevel(factor(model_df_determ_train$genotype), ref = "WT")
determ_df_train$Genotype <- relevel(factor(determ_df_train$Genotype), ref = "WT")

# create overdispersion confirmation function

check_overdispersion <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  list(ratio = ratio, pval = pval)
}

## STABLE DRL
model_df_determ_stable <- model_df[
  model_df$session_type == "MOUSE DRL Overtraining" &
    model_df$session_id %in% 1:6 &
    model_df$block_prob == "1",
]
determ_df_stable <- determ_df[
  determ_df$Task_Stage == "DRL OT" &
    determ_df$Session_Number %in% 1:6,
]

model_df_determ_stable$genotype <- relevel(factor(model_df_determ_stable$genotype), ref = "WT")
determ_df_stable$Genotype <- relevel(factor(determ_df_stable$Genotype), ref = "WT")

## PRL training
model_df_prob_train <- model_df[
  model_df$session_type == "MOUSE PRL" &
    model_df$session_id %in% 1:10 &
    model_df$block_prob == "0.8",
]
prob_df_train <- prob_df[
  prob_df$Task_Stage == "PRL" &
    prob_df$Session_Number %in% 1:10 &
    prob_df$Block_Prob == "0.8",
]

model_df_prob_train$genotype <- relevel(factor(model_df_prob_train$genotype), ref = "WT")
prob_df_train$Genotype <- relevel(factor(prob_df_train$Genotype), ref = "WT")

## Stable PRL
model_df_prob_stable <- model_df[
  model_df$session_type == "MOUSE PRL" &
    model_df$session_id %in% 11:16 &
    model_df$block_prob == "0.8",
]
prob_df_stable <- prob_df[
  prob_df$Task_Stage == "PRL" &
    prob_df$Session_Number %in% 11:16 &
    prob_df$Block_Prob == "0.8",
]

model_df_prob_stable$genotype <- relevel(factor(model_df_prob_stable$genotype), ref = "WT")
prob_df_stable$Genotype <- relevel(factor(prob_df_stable$Genotype), ref = "WT")

## PRL baselining
model_df_prob_base <- model_df[
  model_df$session_type == "MOUSE PRL" &
    model_df$session_id %in% c(21:23, 27:29) &
    model_df$block_prob == "0.8",
]
prob_df_base <- prob_df[
  prob_df$Task_Stage == "PRL" &
    prob_df$Session_Number %in% c(21:23, 27:29) &
    prob_df$Block_Prob == "0.8",
]

model_df_prob_base$genotype <- relevel(factor(model_df_prob_base$genotype), ref = "WT")
prob_df_base$Genotype <- relevel(factor(prob_df_base$Genotype), ref = "WT")

## PRL Probe 1
model_df_prob_probe_1 <- model_df[
  model_df$session_type == "MOUSE PRL" &
    model_df$session_id %in% 18:20 &
    model_df$block_prob == "0.7",
]
prob_df_probe_1 <- prob_df[
  prob_df$Task_Stage == "PRL" &
    prob_df$Session_Number %in% 18:20 &
    prob_df$Block_Prob == "0.7",
]

model_df_prob_probe_1$genotype <- relevel(factor(model_df_prob_probe_1$genotype), ref = "WT")
prob_df_probe_1$Genotype <- relevel(factor(prob_df_probe_1$Genotype), ref = "WT")

## PRL Probe 2
model_df_prob_probe_2 <- model_df[
  model_df$session_type == "MOUSE PRL" &
    model_df$session_id %in% 24:26 &
    model_df$block_prob == "0.6",
]
prob_df_probe_2 <- prob_df[
  prob_df$Task_Stage == "PRL" &
    prob_df$Session_Number %in% 24:26 &
    prob_df$Block_Prob == "0.6",
]

model_df_prob_probe_2$genotype <- relevel(factor(model_df_prob_probe_2$genotype), ref = "WT")
prob_df_probe_2$Genotype <- relevel(factor(prob_df_probe_2$Genotype), ref = "WT")

## PRL Probe 3
model_df_prob_probe_3 <- model_df[
  model_df$session_type == "MOUSE PRL 7 to Reversal" &
    model_df$session_id %in% 1:3 &
    model_df$block_prob == "0.8",
]
prob_df_probe_3 <- prob_df[
  prob_df$Task_Stage == "PRL 7 to Reversal" &
    prob_df$Session_Number %in% 1:3 &
    prob_df$Block_Prob == "0.8",
]

model_df_prob_probe_3$genotype <- relevel(factor(model_df_prob_probe_3$genotype), ref = "WT")
prob_df_probe_3$Genotype <- relevel(factor(prob_df_probe_3$Genotype), ref = "WT")

# Run models

## DRL TRAINING
determ_df_train <- determ_df_train %>%
  filter(!is.na(Reversals))

determ_train_poisson_model <- glmmTMB(
  Reversals ~ Genotype + (1 | Animal_ID),
  data = determ_df_train,
  family = poisson(link = "log")
)

determ_train_od <- check_overdispersion(determ_train_poisson_model)

if (determ_train_od$ratio > 1.2 && determ_train_od$pval < 0.05) {
  message("Overdispersion detected — refitting with Negative Binomial")
  reversals_model_determ_train <- glmmTMB(
    Reversals ~ Genotype + (1 | Animal_ID),
    data = determ_df_train,
    family = nbinom2(link = "log")
  )
} else {
  message("No substantial overdispersion — keeping Poisson")
  reversals_model_determ_train <- determ_train_poisson_model
}

determ_train_model_summary <- tidy(reversals_model_determ_train, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "DRL",
    Metric = "Number of Reversals",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

## STABLE DRL
determ_df_stable <- determ_df_stable %>%
  filter(!is.na(Reversals))

determ_stable_poisson_model <- glmmTMB(
  Reversals ~ Genotype + (1 | Animal_ID),
  data = determ_df_stable,
  family = poisson(link = "log")
)

determ_stable_od <- check_overdispersion(determ_stable_poisson_model)

if (determ_stable_od$ratio > 1.2 && determ_stable_od$pval < 0.05) {
  message("Overdispersion detected — refitting with Negative Binomial")
  reversals_model_determ_stable <- glmmTMB(
    Reversals ~ Genotype + (1 | Animal_ID),
    data = determ_df_stable,
    family = nbinom2(link = "log")
  )
} else {
  message("No substantial overdispersion — keeping Poisson")
  reversals_model_determ_stable <- determ_stable_poisson_model
}

determ_stable_model_summary <- tidy(reversals_model_determ_stable, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "DRL OT",
    Metric = "Number of Reversals",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )


## PRL TRAINING
prob_df_train <- prob_df_train %>%
  filter(!is.na(Reversals))

prob_train_poisson_model <- glmmTMB(
  Reversals ~ Genotype + (1 | Animal_ID),
  data = prob_df_train,
  family = poisson(link = "log")
)

prob_train_od <- check_overdispersion(prob_train_poisson_model)

if (prob_train_od$ratio > 1.2 && prob_train_od$pval < 0.05) {
  message("Overdispersion detected — refitting with Negative Binomial")
  reversals_model_prob_train <- glmmTMB(
    Reversals ~ Genotype + (1 | Animal_ID),
    data = prob_df_train,
    family = nbinom2(link = "log")
  )
} else {
  message("No substantial overdispersion — keeping Poisson")
  reversals_model_prob_train <- prob_train_poisson_model
}

prob_train_model_summary <- tidy(reversals_model_prob_train, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL",
    Metric = "Number of Reversals",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )


## STABLE PRL
prob_df_stable <- prob_df_stable %>%
  filter(!is.na(Reversals))

prob_stable_poisson_model <- glmmTMB(
  Reversals ~ Genotype + (1 | Animal_ID),
  data = prob_df_stable,
  family = poisson(link = "log")
)

prob_stable_od <- check_overdispersion(prob_stable_poisson_model)

if (prob_stable_od$ratio > 1.2 && prob_stable_od$pval < 0.05) {
  message("Overdispersion detected — refitting with Negative Binomial")
  reversals_model_prob_stable <- glmmTMB(
    Reversals ~ Genotype + (1 | Animal_ID),
    data = prob_df_stable,
    family = nbinom2(link = "log")
  )
} else {
  message("No substantial overdispersion — keeping Poisson")
  reversals_model_prob_stable <- prob_stable_poisson_model
}

prob_stable_model_summary <- tidy(reversals_model_prob_stable, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL OT",
    Metric = "Number of Reversals",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )


## PRL BASELINING
prob_df_base <- prob_df_base %>%
  filter(!is.na(Reversals))

prob_base_poisson_model <- glmmTMB(
  Reversals ~ Genotype + (1 | Animal_ID),
  data = prob_df_base,
  family = poisson(link = "log")
)

prob_base_od <- check_overdispersion(prob_base_poisson_model)

if (prob_base_od$ratio > 1.2 && prob_base_od$pval < 0.05) {
  message("Overdispersion detected — refitting with Negative Binomial")
  reversals_model_prob_base <- glmmTMB(
    Reversals ~ Genotype + (1 | Animal_ID),
    data = prob_df_base,
    family = nbinom2(link = "log")
  )
} else {
  message("No substantial overdispersion — keeping Poisson")
  reversals_model_prob_base <- prob_base_poisson_model
}

prob_base_model_summary <- tidy(reversals_model_prob_base, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL Baselining",
    Metric = "Number of Reversals",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )


## PRL PROBE 1
prob_df_probe_1 <- prob_df_probe_1 %>%
  filter(!is.na(Reversals))

prob_probe_1_poisson_model <- glmmTMB(
  Reversals ~ Genotype + (1 | Animal_ID),
  data = prob_df_probe_1,
  family = poisson(link = "log")
)

prob_probe_1_od <- check_overdispersion(prob_probe_1_poisson_model)

if (prob_probe_1_od$ratio > 1.2 && prob_probe_1_od$pval < 0.05) {
  message("Overdispersion detected — refitting with Negative Binomial")
  reversals_model_prob_probe_1 <- glmmTMB(
    Reversals ~ Genotype + (1 | Animal_ID),
    data = prob_df_probe_1,
    family = nbinom2(link = "log")
  )
} else {
  message("No substantial overdispersion — keeping Poisson")
  reversals_model_prob_probe_1 <- prob_probe_1_poisson_model
}

prob_probe_1_model_summary <- tidy(reversals_model_prob_probe_1, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL Probe 1",
    Metric = "Number of Reversals",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )


## PRL PROBE 2
prob_df_probe_2 <- prob_df_probe_2 %>%
  filter(!is.na(Reversals))

prob_probe_2_poisson_model <- glmmTMB(
  Reversals ~ Genotype + (1 | Animal_ID),
  data = prob_df_probe_2,
  family = poisson(link = "log")
)

prob_probe_2_od <- check_overdispersion(prob_probe_2_poisson_model)

if (prob_probe_2_od$ratio > 1.2 && prob_probe_2_od$pval < 0.05) {
  message("Overdispersion detected — refitting with Negative Binomial")
  reversals_model_prob_probe_2 <- glmmTMB(
    Reversals ~ Genotype + (1 | Animal_ID),
    data = prob_df_probe_2,
    family = nbinom2(link = "log")
  )
} else {
  message("No substantial overdispersion — keeping Poisson")
  reversals_model_prob_probe_2 <- prob_probe_2_poisson_model
}

prob_probe_2_model_summary <- tidy(reversals_model_prob_probe_2, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL Probe 2",
    Metric = "Number of Reversals",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )


## PRL PROBE 3
prob_df_probe_3 <- prob_df_probe_3 %>%
  filter(!is.na(Reversals))

prob_probe_3_poisson_model <- glmmTMB(
  Reversals ~ Genotype + (1 | Animal_ID),
  data = prob_df_probe_3,
  family = poisson(link = "log")
)

prob_probe_3_od <- check_overdispersion(prob_probe_3_poisson_model)

if (prob_probe_3_od$ratio > 1.2 && prob_probe_3_od$pval < 0.05) {
  message("Overdispersion detected — refitting with Negative Binomial")
  reversals_model_prob_probe_3 <- glmmTMB(
    Reversals ~ Genotype + (1 | Animal_ID),
    data = prob_df_probe_3,
    family = nbinom2(link = "log")
  )
} else {
  message("No substantial overdispersion — keeping Poisson")
  reversals_model_prob_probe_3 <- prob_probe_3_poisson_model
}

prob_probe_3_model_summary <- tidy(reversals_model_prob_probe_3, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL Probe 3",
    Metric = "Number of Reversals",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

# Combine outputs from models
glmm_summaries <- bind_rows(determ_train_model_summary,
                            determ_stable_model_summary,
                            prob_train_model_summary,
                            prob_stable_model_summary,
                            prob_base_model_summary,
                            prob_probe_1_model_summary,
                            prob_probe_2_model_summary,
                            prob_probe_3_model_summary)

# Export to CSV
write.csv(glmm_summaries, "glmm_summaries_number_of_reversals.csv", row.names = FALSE)