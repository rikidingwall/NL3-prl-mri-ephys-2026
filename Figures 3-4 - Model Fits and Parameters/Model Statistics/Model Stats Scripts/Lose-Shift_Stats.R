# load required libraries 
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
lose_shift_model_determ_train <- glmmTMB(
  p_loss_switch ~ genotype + (1 | mouse_id),
  data = model_df_determ_train,
  family = beta_family(link = "logit")
)

determ_train_model_summary <- tidy(lose_shift_model_determ_train, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "DRL",
    Metric = "Lose-Shift",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

## STABLE DRL
lose_shift_model_determ_stable <- glmmTMB(
  p_loss_switch ~ genotype + (1 | mouse_id),
  data = model_df_determ_stable,
  family = beta_family(link = "logit")
)

determ_stable_model_summary <- tidy(lose_shift_model_determ_stable, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "DRL OT",
    Metric = "Lose-Shift",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

## PRL TRAINING
lose_shift_model_prob_train <- glmmTMB(
  p_loss_switch ~ genotype + (1 | mouse_id),
  data = model_df_prob_train,
  family = beta_family(link = "logit")
)

prob_train_model_summary <- tidy(lose_shift_model_prob_train, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL",
    Metric = "Lose-Shift",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

## STABLE PRL
lose_shift_model_prob_stable <- glmmTMB(
  p_loss_switch ~ genotype + (1 | mouse_id),
  data = model_df_prob_stable,
  family = beta_family(link = "logit")
)

prob_stable_model_summary <- tidy(lose_shift_model_prob_stable, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL OT",
    Metric = "Lose-Shift",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

## PRL BASELINING
lose_shift_model_prob_base <- glmmTMB(
  p_loss_switch ~ genotype + (1 | mouse_id),
  data = model_df_prob_base,
  family = beta_family(link = "logit")
)

prob_base_model_summary <- tidy(lose_shift_model_prob_base, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL Baselining",
    Metric = "Lose-Shift",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

## PRL PROBE 1
lose_shift_model_prob_probe_1 <- glmmTMB(
  p_loss_switch ~ genotype + (1 | mouse_id),
  data = model_df_prob_probe_1,
  family = beta_family(link = "logit")
)

prob_probe_1_model_summary <- tidy(lose_shift_model_prob_probe_1, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL Probe 1",
    Metric = "Lose-Shift",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

## PRL PROBE 2
lose_shift_model_prob_probe_2 <- glmmTMB(
  p_loss_switch ~ genotype + (1 | mouse_id),
  data = model_df_prob_probe_2,
  family = beta_family(link = "logit")
)

prob_probe_2_model_summary <- tidy(lose_shift_model_prob_probe_2, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL Probe 2",
    Metric = "Lose-Shift",
    term = ifelse(term == "(Intercept)", "ReferenceWT", term)
  )

## PRL PROBE 3
lose_shift_model_prob_probe_3 <- glmmTMB(
  p_loss_switch ~ genotype + (1 | mouse_id),
  data = model_df_prob_probe_3,
  family = beta_family(link = "logit")
)

prob_probe_3_model_summary <- tidy(lose_shift_model_prob_probe_3, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Task_Stage = "PRL Probe 3",
    Metric = "Lose-Shift",
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
write.csv(glmm_summaries, "glmm_summaries_lose_shift.csv", row.names = FALSE)
