# 02_transform.R ####
# Author: Simon Dedman
# Date: 2026-01-28
# Purpose: Test assumptions and transform data as needed
#
# Description:
#   1. Test assumptions (homogeneity, normality)
#   2. Diagnostic plots (residuals, QQ)
#   3. Distributional transformations (log, sqrt, reciprocal, power, arcsin)
#   4. Outlier handling (winsorize)
#   5. Scaling for modeling (stdize, range_stdize)
#   6. Encoding categorical variables (one-hot)
#
# Inputs:
#   - Cleaned data from 01_tidy-data.R
#
# Outputs:
#   - Diagnostic plots in outputs/02_transform/
#   - Transformed data objects




# Functions ####
## Assumption test functions ####

#' Test homogeneity of variance
#' Runs Bartlett and Levene tests
#' @param formula Formula: response ~ group
#' @param data Data frame
#' @return List with test results
#'
test_homogeneity <- function(formula, data) {
  bartlett_result <- bartlett.test(formula, data = data)
  levene_result <- car::leveneTest(formula, data = data)
  list(
    bartlett = bartlett_result,
    levene = levene_result
  )
}

#' Test normality
#' Runs Shapiro-Wilk (n < 50) or Kolmogorov-Smirnov (n >= 50)
#'
#' @param x Numeric vector (typically residuals)
#' @return Test result
#'
test_normality <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)

  if (n < 3) {
    warning("Need at least 3 observations for normality test")
    return(NULL)
  }

  if (n < 50) {
    shapiro.test(x)
  } else {
    ks.test(x, "pnorm", mean(x), sd(x))
  }
}


## Diagnostic plot functions ####

#' Plot residuals vs fitted values
#' @param model Fitted model object
#' @param output_path Optional path to save plot
#' @return ggplot object
#'
plot_residuals_fitted <- function(model, output_path = NULL) {
  df <- data.frame(
    fitted = fitted(model),
    residuals = residuals(model)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = fitted, y = residuals)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "blue") +
    ggplot2::labs(
      title = "Residuals vs Fitted Values",
      x = "Fitted Values",
      y = "Residuals"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(output_path, plot = p, width = 8, height = 6)
  }
  return(p)
}


#' Plot histogram of residuals
#' @param residuals Numeric vector of residuals
#' @param output_path Optional path to save plot
#' @return ggplot object
#'
plot_residuals_histogram <- function(residuals, output_path = NULL) {
  df <- data.frame(residuals = residuals)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = residuals)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    ggplot2::labs(
      title = "Histogram of Residuals",
      x = "Residuals",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(output_path, plot = p, width = 8, height = 6)
  }
  return(p)
}

#' Plot QQ plot of residuals
#' @param residuals Numeric vector of residuals
#' @param output_path Optional path to save plot
#' @return ggplot object
#'
plot_qq <- function(residuals, output_path = NULL) {
  df <- data.frame(residuals = residuals)

  p <- ggplot2::ggplot(df, ggplot2::aes(sample = residuals)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = "red") +
    ggplot2::labs(
      title = "Normal Q-Q Plot",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(output_path, plot = p, width = 8, height = 6)
  }
  return(p)
}


## Distributional transformation functions ####

#' Log transformation (safe for zeros)
#' @param x Numeric vector
#' @param base Log base (default: natural log)
#' @param add_constant Value to add before log if zeros present (default: 1)
#' @return Transformed vector with attributes for reversal
#'
transform_log <- function(x, base = exp(1), add_constant = 1) {
  has_zeros <- any(x == 0, na.rm = TRUE)
  has_negatives <- any(x < 0, na.rm = TRUE)

  if (has_negatives) stop("Log transformation not valid for negative values")

  if (has_zeros) {
    x_trans <- log(x + add_constant, base = base)
    attr(x_trans, "add_constant") <- add_constant
  } else {
    x_trans <- log(x, base = base)
    attr(x_trans, "add_constant") <- 0
  }
  attr(x_trans, "base") <- base
  return(x_trans)
}


#' Reverse log transformation
#' @param x Log-transformed vector
#' @return Original scale values
#'
untransform_log <- function(x) {
  base <- attr(x, "base")
  add_constant <- attr(x, "add_constant")
  if (is.null(base)) stop("Missing transformation attributes")
  x_orig <- base^x - add_constant
  return(x_orig)
}


#' Square root transformation
#' @param x Numeric vector
#' @param add_constant Value to add if zeros/negatives (default: 0)
#' @return Transformed vector with attributes
#'
transform_sqrt <- function(x, add_constant = 0) {
  min_val <- min(x, na.rm = TRUE)

  if (min_val < 0) {
    add_constant <- abs(min_val) + 0.001
    warning(paste("Negative values detected. Adding", add_constant, "before sqrt"))
  } else if (min_val == 0 && add_constant == 0) {
    add_constant <- 0.001
  }
  x_trans <- sqrt(x + add_constant)
  attr(x_trans, "add_constant") <- add_constant
  return(x_trans)
}


#' Reverse square root transformation
#' @param x Sqrt-transformed vector
#' @return Original scale values
#'
untransform_sqrt <- function(x) {
  add_constant <- attr(x, "add_constant")
  if (is.null(add_constant)) stop("Missing transformation attributes")
  x_orig <- x^2 - add_constant
  return(x_orig)
}


#' Reciprocal transformation
#' @param x Numeric vector
#' @param add_constant Value to add if zeros (default: 1)
#' @return Transformed vector with attributes
#'
transform_reciprocal <- function(x, add_constant = 1) {
  has_zeros <- any(x == 0, na.rm = TRUE)

  if (has_zeros) {
    x_trans <- 1 / (x + add_constant)
    attr(x_trans, "add_constant") <- add_constant
  } else {
    x_trans <- 1 / x
    attr(x_trans, "add_constant") <- 0
  }
  return(x_trans)
}


#' Reverse reciprocal transformation
#' @param x Reciprocal-transformed vector
#' @return Original scale values
#'
untransform_reciprocal <- function(x) {
  add_constant <- attr(x, "add_constant")
  if (is.null(add_constant)) stop("Missing transformation attributes")

  x_orig <- (1 / x) - add_constant
  return(x_orig)
}


#' Power transformation
#' @param x Numeric vector
#' @param power Power to raise to (default: 2)
#' @return Transformed vector with attributes
#'
transform_power <- function(x, power = 2) {
  x_trans <- x^power
  attr(x_trans, "power") <- power

  return(x_trans)
}


#' Reverse power transformation
#' @param x Power-transformed vector
#' @return Original scale values
#'
untransform_power <- function(x) {
  power <- attr(x, "power")
  if (is.null(power)) stop("Missing transformation attributes")

  x_orig <- x^(1 / power)
  return(x_orig)
}


#' Arcsin (angular) transformation for proportions
#' @param x Numeric vector of proportions (0-1)
#' @param n Sample size (for edge corrections)
#' @return Transformed vector
#'
transform_arcsin <- function(x, n = NULL) {
  if (any(x < 0 | x > 1, na.rm = TRUE)) {
    stop("Arcsin transformation requires values between 0 and 1")
  }
  # Edge corrections for 0 and 1 values
  if (!is.null(n)) {
    x <- ifelse(x == 0, 1 / (4 * n), x)
    x <- ifelse(x == 1, 1 - 1 / (4 * n), x)
  }
  x_trans <- asin(sqrt(x))
  return(x_trans)
}


#' Reverse arcsin transformation
#' @param x Arcsin-transformed vector
#' @return Original proportions
#'
untransform_arcsin <- function(x) {
  x_orig <- sin(x)^2
  return(x_orig)
}


## Outlier functions ####

#' Winsorize outliers
#' Replace extreme values with specified percentile values
#'
#' @param x Numeric vector
#' @param probs Lower and upper probability thresholds (default: c(0.05, 0.95))
#' @return Winsorized vector with attributes
#'
winsorize <- function(x, probs = c(0.05, 0.95)) {
  limits <- quantile(x, probs = probs, na.rm = TRUE)
  x_wins <- x
  x_wins[x < limits[1]] <- limits[1]
  x_wins[x > limits[2]] <- limits[2]
  attr(x_wins, "lower_limit") <- limits[1]
  attr(x_wins, "upper_limit") <- limits[2]
  attr(x_wins, "n_winsorized_lower") <- sum(x < limits[1], na.rm = TRUE)
  attr(x_wins, "n_winsorized_upper") <- sum(x > limits[2], na.rm = TRUE)
  return(x_wins)
}


## Scaling functions ####

#' Standardize to 0-1 scale (centre, scale to SD, rescale)
#' Centres around mean, scales to SD, then rescales to 0:1
#' Stores attributes for reverse transformation
#' @param x Numeric vector to standardize
#' @return Standardized vector with attributes for reversal
#'
stdize <- function(x) {
  orig_mean <- mean(x, na.rm = TRUE)
  orig_sd <- sd(x, na.rm = TRUE)

  # Center & scale to SD
  x_scaled <- (x - orig_mean) / orig_sd

  # Shift so min is 0
  scaled_min <- min(x_scaled, na.rm = TRUE)
  x_shifted <- x_scaled - scaled_min
  # Scale so max is 1
  shifted_max <- max(x_shifted, na.rm = TRUE)
  x_final <- x_shifted / shifted_max

  # Store params for reversal
  attr(x_final, "orig_mean") <- orig_mean
  attr(x_final, "orig_sd") <- orig_sd
  attr(x_final, "scaled_min") <- scaled_min
  attr(x_final, "shifted_max") <- shifted_max

  return(x_final)
}

#' Reverse stdize transformation
#' @param x Standardized vector (with attributes from stdize)
#' @return Original scale values
#'
unstdize <- function(x) {
  orig_mean <- attr(x, "orig_mean")
  orig_sd <- attr(x, "orig_sd")
  scaled_min <- attr(x, "scaled_min")
  shifted_max <- attr(x, "shifted_max")
  if (is.null(orig_mean)) stop("Missing transformation attributes. Use stdize() first.")

  # Reverse: multiply by max, add min, multiply by sd, add mean
  x_shifted <- x * shifted_max
  x_scaled <- x_shifted + scaled_min
  x_orig <- (x_scaled * orig_sd) + orig_mean

  return(x_orig)
}

#' Range standardization (Steinley 2004, 2006)
#' zi = (xi - min(x)) / (max(x) - min(x))
#' Stores attributes for reverse transformation
#'
#' @param x Numeric vector to standardize
#' @return Range-standardized vector (0-1) with attributes for reversal
#'
range_stdize <- function(x) {
  orig_min <- min(x, na.rm = TRUE)
  orig_max <- max(x, na.rm = TRUE)
  x_final <- (x - orig_min) / (orig_max - orig_min)
  attr(x_final, "orig_min") <- orig_min
  attr(x_final, "orig_max") <- orig_max
  return(x_final)
}


#' Reverse range standardization
#' @param x Range-standardized vector (with attributes from range_stdize)
#' @return Original scale values
#'
unrange_stdize <- function(x) {
  orig_min <- attr(x, "orig_min")
  orig_max <- attr(x, "orig_max")
  if (is.null(orig_min)) stop("Missing transformation attributes. Use range_stdize() first.")
  x_orig <- x * (orig_max - orig_min) + orig_min
  return(x_orig)
}


## Encoding functions ####

#' One-hot encode a categorical column
#' Creates dummy variables using tidyverse map/unnest approach
#'
#' @param df Data frame
#' @param col Name of column to encode (unquoted or string)
#' @param drop_original Drop original column after encoding (default: TRUE)
#' @return Data frame with one-hot encoded columns
#'
#' @examples
#' iris |> one_hot_encode(Species)
#' df |> one_hot_encode(category, drop_original = FALSE)
#'
one_hot_encode <- function(df, col, drop_original = TRUE) {
  col_name <- rlang::ensym(col)
  col_str <- rlang::as_string(col_name)

  # Ensure column is a factor
  if (!is.factor(df[[col_str]])) {
    df[[col_str]] <- as.factor(df[[col_str]])
  }

  lvls <- levels(df[[col_str]])

  df <- df |>
    dplyr::mutate(
      "{col_str}_one_hot" := purrr::map(
        .data[[col_str]],
        ~ rlang::set_names(lvls == .x, lvls)
      )
    ) |>
    tidyr::unnest_wider(tidyselect::all_of(paste0(col_str, "_one_hot")))

  if (drop_original) {
    df[[col_str]] <- NULL
  }

  return(df)
}




# Load data ####
# Load cleaned data from previous script
# source(here::here("R", "01_tidy-data.R"))
# Or load saved data:
# load(here::here("data", "cleaned_data.rda"))




# Test assumptions ####
## Example usage ####
# Fit a model first
# model <- lm(response ~ predictor, data = df)

# Test homogeneity of variance
# test_homogeneity(response ~ group, data = df)

# Test normality of residuals
# test_normality(residuals(model))

# Plot diagnostics
# plot_residuals_fitted(model, here::here("outputs", "02_transform", "residuals_fitted.png"))
# plot_residuals_histogram(residuals(model), here::here("outputs", "02_transform", "residuals_hist.png"))
# plot_qq(residuals(model), here::here("outputs", "02_transform", "qq_plot.png"))




# Apply distributional transformations ####
## Example usage ####
# Log transform
# df$var_log <- transform_log(df$var)
# df$var_original <- untransform_log(df$var_log)

# Square root transform
# df$var_sqrt <- transform_sqrt(df$var)
# df$var_original <- untransform_sqrt(df$var_sqrt)




# Handle outliers ####
## Example usage ####
# df$var_wins <- winsorize(df$var, probs = c(0.05, 0.95))




# Scale for modeling ####
## Example usage ####
# Standardize (0-1 with SD scaling)
# df$var_std <- stdize(df$var)
# df$var_original <- unstdize(df$var_std)

# Range standardize (simple 0-1)
# df$var_range <- range_stdize(df$var)
# df$var_original <- unrange_stdize(df$var_range)




# Encode categorical variables ####
## Example usage ####
# df_encoded <- df |>
#   one_hot_encode(category1) |>
#   one_hot_encode(category2)
