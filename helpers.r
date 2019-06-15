library(readstata13)
library(tidyverse)
library(outliers)

# Analyze The Missing Values
na_table <- function(x) {
  na_table <- data.frame()
  for (i in 1:ncol(x)) {
    n_na <- nrow(x[is.na(x[,i]),])
    na_ratio <- n_na / nrow(x)
    na_table[i, 1] <- colnames(x)[[i]]
    na_table[i, 2] <- n_na
    na_table[i, 3] <- na_ratio
    colnames(na_table) <- c("covariate", "n", "ratio")
  }
  return(na_table)
}

# Load And pre-process the dataset
load_endline1 <- function() {
  endlines <- read.dta13("data/2013-0533_data_endlines1and2.dta",
                         convert.factors = FALSE,
                         generate.factors = TRUE)

  endline1 <- endlines %>%
    filter(sample1 == 1) %>%
    select(colnames(endlines)[1:16],
           contains("_1"),
           -c(w, w1, w2, sample1, sample2, visitday_1, visitmonth_1, visityear_1),
           -starts_with("area_"),
           -ends_with("_mo_1"),
           -ends_with("_annual_1")) %>%
    mutate(old_biz = ifelse(any_old_biz == 0 | is.na(any_old_biz) == TRUE,
                            0,
                            old_biz),
           total_biz_1 = ifelse(any_biz_1 == 0 | is.na(any_biz_1) == TRUE,
                                0,
                                total_biz_1),
           newbiz_1 = ifelse(any_new_biz_1 == 0 | is.na(any_new_biz_1) == TRUE,
                             0,
                             newbiz_1)) %>%
    select(-c(any_old_biz, any_biz_1, any_new_biz_1,
              hhsize_1,
              anymfi_1, anymfi_amt_1,
              hours_week_1, hours_headspouse_week_1, hours_child1620_week_1,
              total_exp_mo_pc_1))

  na_delete_threshold <- 0.1
  na_delete_col <- (na_table(endline1) %>% filter(ratio > na_delete_threshold))[,1]
  # delete those variables
  for (col in na_delete_col) {
    endline1[,col] <- NULL
  }

  endline1 <- endline1 %>%
    mutate(bizassets_1 = ifelse(total_biz_1 == 0 | is.na(total_biz_1),
                                0,
                                bizassets_1),
           bizinvestment_1 = ifelse(total_biz_1 == 0 | is.na(total_biz_1),
                                    0,
                                    bizinvestment_1),
           bizrev_1 = ifelse(total_biz_1 == 0 | is.na(total_biz_1),
                             0,
                             bizrev_1),
           bizexpense_1 = ifelse(total_biz_1 == 0 | is.na(total_biz_1),
                                 0,
                                 bizexpense_1),
           bizprofit_1 = ifelse(total_biz_1 == 0 | is.na(total_biz_1),
                                0,
                                bizprofit_1),
           bizemployees_1 = ifelse(total_biz_1 == 0 | is.na(total_biz_1),
                                   0,
                                   bizemployees_1))

  covariates_name <- endline1 %>%
    select(-contains("index")) %>%
    colnames()
  for (covar in covariates_name) {
    endline1[is.na(endline1[, covar]), covar] <-
      median(endline1[, covar], na.rm = TRUE)
  }

  endline1 <- na.omit(endline1)

  exp_col <- endline1 %>%
    select(contains("exp_mo_pc")) %>%
    colnames()
  for (covar in exp_col) {
    covar_outlier <- scores(x = endline1[, covar], type = "iqr", lim = 5)
    endline1 <- endline1[!covar_outlier, ]
  }

  return(endline1)
}


# Two-models Approach (with Sorted Groups ATE)
## This function contains two parts.
tm_gates <- function(target, treatment, data,
                     split_ratio=0.5, cluster=0, num.iter=100,
                     ml_method="rf") {

  # a list to store the results
  results <- list()

  # TODO implement user specified num.groups
  num.groups <- 5

  if (cluster != 0) {
    strati_target <- cluster
  } else {
    strati_target <- treatment
  }

  # In the Chernozhukov paper, they suggested that we run the model multiple times
  # to mitigate the sample selection bias (because we split the dataset).
  # And take the medium of the final results
  for (i in 1:num.iter) {
    # set seed for reproduction
    set.seed(i)

    # seperate auxi and main sample
    auxi_index <- createDataPartition(data[,strati_target],
                                      p = split_ratio,
                                      list = FALSE)
    auxi <- data[auxi_index, ]
    main <- data[-auxi_index,]

    # seperate treatment & control in the auxiliary sample
    auxi_treat_index <- which(auxi[,treatment] == 1)
    auxi_treat <- auxi[auxi_treat_index, ]
    auxi_contr <- auxi[-auxi_treat_index,]

    # use the specified machine learning method to predict the conditional treatment effect
    if (ml_method == "rf") { # use two RANDOM FOREST to predict the CTE
      # fit a random forest on auxi_treat and auxi_contr
      auxi_formula <- as.formula(paste(target, " ~ .", "-", treatment))
      # we use only the control group to fit a model that predicts the baseline value
      auxi_yi0 <- randomForest(auxi_formula,
                               data = auxi_contr,
                               ntree = 3000,
                               mtry = 3,
                               replace = TRUE,
                               type = "regression")
      # and use only the treatment group to fit a model that predict the treated value
      auxi_yi1 <- randomForest(auxi_formula,
                               data = auxi_treat,
                               ntree = 3000,
                               mtry = 3,
                               replace = TRUE,
                               type = "regression")
      # predict the baseline and treated value on main sample
      main_yi0 <- predict(auxi_yi0, newdata = main)
      main_yi1 <- predict(auxi_yi1, newdata = main)
      main$baseline <- main_yi0
      # we assumed the difference between baseline value and treated value is
      # conditional treatment effect (CTE)
      main$cte <- (main_yi1 - main_yi0)
    } else if (ml_method == "crf") { # use CAUSAL RANFOM FOREST to predict the CTE
      # fit a causal random forest on auxi sample
      auxi_X <- auxi %>%
        select(everything(), -target, -treatment)
      auxi_Y <- auxi[, target]
      auxi_W <- auxi[, treatment]
      auxi_crf <- causal_forest(X = auxi_X,
                                Y = auxi_Y,
                                W = auxi_W,
                                honesty = TRUE,
                                mtry = 3,
                                num.trees = 3000)
      # predict the conditional treatment effect on main sample
      auxi_crf_pred <- predict(auxi_crf, newdata = main)
      main$cte <- auxi_crf_pred$predictions
    }

    # 1. TWO-MODELS APPROACH
    # Fit regression on conditional treatment effect
    tm_exclude_col <- c(target, treatment, cluster,
                        "baseline", "cte")
    data_col <- names(main)
    tm_formula <- as.formula(
      paste(
        "cte", "~",
        paste(data_col[!data_col %in% tm_exclude_col], collapse = " + ")))

    tm_model <- lm(tm_formula, data = main)

    results[[i]] <- tm_model

    # 2. SORTED GROUP AVERAGE TREATMENT EFFECT
    # calculate propensity score (treated/all)
    # TODO implement option to use non-randomized treatment assignment
    prop_score <- nrow(data[data$treatment == 1, ])/nrow(data)
    main$prop_score <- prop_score

    # divide observations based on their predicted conditional treatment effect
    breaks <- quantile(main$cte, seq(0,1, 1/num.groups), include.lowest = TRUE)
    breaks[1] <- breaks[1] - 0.001
    breaks[6] <- breaks[6] + 0.001
    main$treat_group <- cut(main$cte, breaks = breaks)

    # calculate the propensity score offset for each observation in main sample
    main$prop_offset <- main$treatment - main$prop_score

    # construct matrix from each observation's group factor
    SGX <- model.matrix(~-1+main$treat_group)
    # construct D-p(X)*1(G_k) and weight for each observation
    DSG <- data.frame(main$prop_offset*SGX)
    colnames(DSG) <- c("G1", "G2", "G3", "G4", "G5")
    main[,c("G1", "G2", "G3", "G4", "G5", "weight")] <- cbind(
      DSG$G1, DSG$G2, DSG$G3, DSG$G4, DSG$G5,
      1/prop_score*(1-prop_score))

    # fit weighted ols
    if (ml_method == "rf") {
      gates_formula <- as.formula(paste(target,
                                        "~",
                                        "-1+baseline+cte+G1+G2+G3+G4+G5",
                                        "|0|0|",
                                        cluster))
    } else if (ml_method == "crf") {
      gates_formula <- as.formula(paste(target,
                                        "~",
                                        "cte+G1+G2+G3+G4+G5",
                                        "|0|0|",
                                        cluster))
    }

    # use Weighted OLS
    gates_model <- felm(gates_formula,
                        data = main,
                        weights = main$weight)

    results[[num.iter+i]] <- gates_model
  }
  return(results)
}

# Plot The Variable Importance Plot for CRF
var_imp_plot <- function(forest, decay.exponent = 2L, max.depth = 4L) {

  # Calculate variable importance of all features
  # (from print.R)
  split.freq <- split_frequencies(forest, max.depth)
  split.freq <- split.freq / pmax(1L, rowSums(split.freq))
  weight <- seq_len(nrow(split.freq)) ^ -decay.exponent
  var.importance <- t(split.freq) %*% weight / sum(weight)

  # Format data frame
  p <- ncol(forest$X.orig)

  var.names <- colnames(forest$X.orig)[seq_len(p)]
  if (is.null(var.names)) {
    var.names <- paste0('x', seq_len(p))
  }
  df <- tibble(Variable = var.names,
               Importance = as.numeric(var.importance)) %>%
    arrange(Importance) %>%
    mutate(Variable = factor(Variable, levels = unique(Variable)))

  # Plot results
  p <- ggplot(df, aes(Variable, Importance)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    ggtitle('Variable Importance') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}

# Plot the Trend of Variables to CATE
trend_plots <- function(crf, test) {
  # Get the variable importance table
  var_imp <- crf %>%
    variable_importance() %>%
    as.data.frame() %>%
    mutate(variable = colnames(crf$X.orig)) %>%
    arrange(desc(V1))
  # for the first four most important variable
  # create a plot that shows if there are trend of correlation
  p1 <- ggplot(test, aes(x = test[, var_imp$variable[1]], y = preds)) +
    geom_point() +
    geom_smooth(method = "loess", span = 1) +
    theme_light() +
    labs(x = var_imp$variable[1], y = "pred. CTE")
  p2 <- ggplot(test, aes(x = test[, var_imp$variable[2]], y = preds)) +
    geom_point() +
    geom_smooth(method = "loess", span = 1) +
    theme_light() +
    labs(x = var_imp$variable[2], y = "pred. CTE")
  p3 <- ggplot(test, aes(x = test[, var_imp$variable[3]], y = preds)) +
    geom_point() +
    geom_smooth(method = "loess", span = 1) +
    theme_light() +
    labs(x = var_imp$variable[3], y = "pred. CTE")
  p4 <- ggplot(test, aes(x = test[, var_imp$variable[4]], y = preds)) +
    geom_point() +
    geom_smooth(method = "loess", span = 1) +
    theme_light() +
    labs(x = var_imp$variable[4], y = "pred. CTE")

  # combine those plots
  cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)
}
