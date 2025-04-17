get_response <- function(formula) {
  if (class(formula) != "formula") {
    stop("Not formula object")
  }
  all.vars(formula[[length(formula) - 1]])
}

get_covs <- function(formula) {
  if (class(formula) != "formula") {
    stop("Not formula object")
  }
  all.vars(formula[[length(formula)]])
}

robust_glm <- function(formula, data, subset = NULL, family = "quasibinomial",
                       robust_weights = T, sandwich = T, add_ci = T, p = NULL,
                       ...) {
  #' Fit a Generalized Linear Model with Robust Weights and Errors
  #'
  #' This function fits a generalized linear model (GLM) with robust weights per
  #'  observation, to down-weight outliers, and robust sandwich standard errors,
  #'   to account for the lack of independance or heteroscedasticity.
  #'
  #' @param formula A formula specifying the model.
  #' @param data A data frame containing the variables in the model.
  #' @param subset An optional character vector specifying a subset of observations
  #'               to be used in the model.
  #' @param family A character string or function (see \code{lm()}) specifying
  #' the distribution family in the GLM (default is "binomial").
  #' @param robust_weights Logical indicating whether to compute robust model weights
  #'                       (default is TRUE).
  #' @param sandwich Logical indicating whether to compute sandwich standard errors
  #'                 (default is TRUE).
  #' @param add_ci Logical indicating whether to add confidence intervals to the model
  #'               coefficients (default is TRUE).
  #' @param p An optional progressor object to monitor progress (default is NULL).
  #' @param ... Other arguments passed on to \code{stats::glm}.
  #'
  #' @return An object of class \code{"glm"} with additional attributes such as
  #'         confidence intervals, sandwich standard errors, and collinear terms.
  #'
  #' @examples
  #' df <- data.frame(y = runif(50, 0, 1),
  #'                  x1 = rep(1:2, 50),
  #'                  x2 = runif(50, .5, 1),
  #'                  x3 = runif(50, 0, .5))
  #' fit <- robust_glm(y ~ x1 + x2, df, family = "quasibinomial")
  #' summary(fit)
  #'
  #' @importFrom MASS rlm
  #' @importFrom bestNormalize bestNormalize
  #' @importFrom sandwich vcovHC
  #' @importFrom lmtest coeftest coefci
  #' @export

  # Input validation
  if (missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' arguments must be provided.")
  }
  if (!inherits(formula, "formula") && !is.character(formula)) {
    stop("Argument 'formula' must be either a character formula or a formula object.")
  }
  if (!inherits(formula, "formula") && is.character(formula)) {
    formula <- as.formula(formula)
  }
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame.")
  }
  if (!is.null(subset) && !is.character(subset)) {
    stop("Argument 'subset' must be a character vector.")
  }

  # Subset data
  if (!is.null(subset)) {
    data <- subset(data, eval(parse(text = subset)))
  }

  # Extract model terms
  y <- get_response(formula)
  covs <- get_covs(formula)

  # Check if terms are present in the data
  missing_terms <- c(y, covs)[!(c(y, covs) %in% names(data))]
  if (length(missing_terms) > 0) {
    stop("Variable(s) '", paste(missing_terms, collapse = "', '"),
         "' not found in the data.")
  }

  # Clean data
  data <- data[, c(y, covs)]
  data <- na.omit(data)

  # Extract the response variable values
  prop <- data[[y]]

  # Set up the progressor progress bar
  if (!is.null(p) && inherits(p, "progressor")) {
    p()
  }

  # Check for collinear terms
  x <- model.matrix(formula, data)
  ncovs <- ncol(x)
  QR <- qr(x)
  if (QR$rank < ncovs) {
    collinear_terms <- colnames(x)[QR$pivot[(QR$rank + 1):ncovs]]
    if (all(collinear_terms %in% covs)) {
      formula <- reformulate(covs[!covs %in% collinear_terms], response = y)
    }
    collinear_terms <- paste(collinear_terms, collapse = ", ")

    warning("Dropped collinear terms: ", collinear_terms)
  } else {
    collinear_terms <- NULL
  }

  # Compute robust model weights if requested
  if (robust_weights && is.numeric(data[[y]]) &&
      length(unique(data[[y]])) > 2) {
    # From [0, 1] to (-Inf, +Inf)
    data$prop_norm <- bestNormalize::bestNormalize(prop, loo = T, quiet = T)$x.t
    robust_formula <- update.formula(formula, prop_norm ~ .)
    rweights <- MASS::rlm(robust_formula, data = data)$w

    if (length(rweights) != length(prop)) {
      stop("Weight and response have different lengths. Any NA maybe?")
    }

    data$rweights <- rweights

    # Fit a generalized linear model with robust weights
    fit <- stats::glm(formula, data = data, family = family, weights = rweights,
                      ...)
  } else {
    fit <- stats::glm(formula, data = data, family = family, ...)
  }

  if (add_ci) {
    fit$ci <- suppressMessages(confint(fit))
    colnames(fit$ci) <- c("conf.low", "conf.high")
  }

  # Add sandwich errors
  if (sandwich) {
    vcov_mat <- sandwich::vcovHC(fit, type = "HC3")
    fit$sandwich <- lmtest::coeftest(fit, vcov. = vcov_mat)
    sandwich_ci <- lmtest::coefci(fit, vcov. = vcov_mat)
    fit$sandwich <- cbind(fit$sandwich, sandwich_ci)
    colnames(fit$sandwich) <- c("estimate", "std.error", "statistic", "p.value",
                                "conf.low", "conf.high")
  }

  # Add collinear terms
  fit$collinear_terms <- collinear_terms

  # Add class robust_glm
  class(fit) <- c(class(fit), "robust_glm")

  return(fit)
}

tidy_terms <- function(model, id = NULL, exponentiate = F) {
  #' Tidy Model Terms
  #'
  #' This function tidies model terms, including coefficients, confidence intervals,
  #' and sandwich standard errors, from a list of model objects or a single model object.
  #'
  #' @param model A list of or a single robust_glm object.
  #' @param id A character string specifying the identifier column name
  #'  (default is "id").
  #' @param exponentiate Logical. Exponentiate estimate and confidence
  #'  intervals.
  #'
  #' @return A data frame containing tidied model terms with columns:
  #' \describe{
  #'   \item{\code{id}}{Identifier column name.}
  #'   \item{\code{term}}{Model term names.}
  #'   \item{\code{estimate}}{Estimate of coefficients.}
  #'   \item{\code{std.error}}{Standard error of coefficients.}
  #'   \item{\code{statistic}}{Value of test statistics.}
  #'   \item{\code{p.value}}{p-value of test statistics.}
  #'   \item{\code{conf.low}}{Lower bound of confidence interval.}
  #'   \item{\code{conf.high}}{Upper bound of confidence interval.}
  #'   \item{\code{std.error_hc}}{Sandwich standard error of coefficients.}
  #'   \item{\code{statistic_hc}}{Value of test statistics with sandwich
  #'    standard errors.}
  #'   \item{\code{p.value_hc}}{p-value of test statistics with sandwich
  #'    standard errors.}
  #'   \item{\code{conf.low_hc}}{Lower bound of confidence interval with
  #'    sandwich standard errors.}
  #'   \item{\code{conf.high_hc}}{Upper bound of confidence interval with
  #'    sandwich standard errors.}
  #'    \item{\code{estimate_exp}}{Exponentiated estimate.}
  #'   \item{\code{conf.low_exp}}{Lower bound of exponentiated confidence
  #'    interval.}
  #'   \item{\code{conf.high_exp}}{Upper bound of exponentiated confidence
  #'    interval.}
  #'   \item{\code{std.error_exp}}{Exponentiated standard error.}
  #'   \item{\code{std.error_hc_exp}}{Exponentiated sandwich standard error.}
  #' }
  #'
  #' @importFrom broom tidy
  #' @importFrom tibble rownames_to_column
  #' @importFrom dplyr left_join select rename_with
  #' @importFrom purrr map_dfr reduce
  #' @export

  # If lm or glm supplied, convert to a list
  if (inherits(model, c("lm", "glm"))) {
    model <- list(model)
  }

  # Check if model contains at least one model object
  if (length(model) == 0) {
    stop("model must contain at least one model object.")
  }

  # Tidy terms for each model object and combine results
  terms_df <- purrr::map_dfr(model, function(mod) {
    # Check if mod is a valid model object
    if (!inherits(mod, c("lm", "glm"))) {
      stop("model must be a valid 'lm' or 'glm' object.")
    }
    if (!inherits(mod, "robust_glm")) {
      stop("model must be of class robust_glm.")
    }

    # Check if mod contains ci and sandwich, and they are not NA
    ci_df <- if (!is.null(mod$ci) && !all(is.na(mod$ci))) {
      mod$ci %>%
        as.data.frame() %>%
        tibble::rownames_to_column("term")
    } else {
      data.frame(term = names(coef(mod)), conf.low = NA, conf.high = NA)
    }

    sandwich_df <- if (!is.null(mod$sandwich) && !all(is.na(mod$sandwich))) {
      mod$sandwich %>%
        as.data.frame() %>%
        dplyr::select(-estimate) %>%
        dplyr::rename_with(~ paste0(.x, "_hc")) %>%
        tibble::rownames_to_column("term")
    } else {
      data.frame(term = names(coef(mod)), std.error_hc = NA, statistic_hc = NA,
                 p.value_hc = NA, conf.low_hc = NA, conf.high_hc = NA)
    }

    # Tidy model terms and combine with ci and sandwich
    df <- purrr::reduce(list(broom::tidy(mod), ci_df, sandwich_df),
                        dplyr::left_join, by = "term")

    # Exponentiate estimate and confidence intervals
    if (exponentiate) {
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::matches("conf|estimate"), exp,
                                    .names = "{.col}_exp"))
    }

    return(df)
  }, .id = id)

  return(terms_df)
}

tidy_model <- function(model, id = NULL, ...) {
  #' Tidy Model Summary Statistics
  #'
  #' Tidies summary statistics of model objects, such as convergence status,
  #' dispersion, and performance metrics.
  #'
  #' @param model A list of or a single robust_glm object.
  #' @param id A character string specifying the identifier column name
  #'  (default is "id").
  #' @param ... Arguments passed to \code{broom::tidy()}.
  #'
  #' @return A data frame containing tidied model summary statistics with
  #'  columns:
  #' \describe{
  #'   \item{\code{id}}{Identifier column name.}
  #'   \item{\code{converged}}{Convergence status of the model.}
  #'   \item{\code{dispersion}}{Dispersion value of the model.}
  #'   \item{\code{glance}}{Summary statistics of the model.}
  #'   \item{\code{rmse}}{Root mean square error.}
  #'   \item{\code{mse}}{Mean square error.}
  #'   \item{\code{hosmer_chisq}}{Hosmer goodness-of-fit test statistic.}
  #'   \item{\code{hosmer_df}}{Degrees of freedom for Hosmer goodness-of-fit
  #'    test.}
  #'   \item{\code{hosmer_p.value}}{p-value of Hosmer goodness-of-fit test.}
  #'   \item{\code{R2_Tjur}}{Tjur's R-squared.}
  #'   \item{\code{Log_loss}}{Logarithmic loss.}
  #' }
  #'
  #' @importFrom purrr map_dfr
  #' @importFrom broom glance
  #' @importFrom performance performance_rmse performance_mse performance_hosmer
  #' @importFrom performance r2_tjur performance_logloss
  #' @importFrom tibble as_tibble
  #' @importFrom dplyr bind_cols
  #' @export

  # If lm or glm supplied, convert to a list
  if (inherits(model, c("lm", "glm"))) {
    model <- list(model)
  }

  # Check if model contains at least one model object
  if (length(model) == 0) {
    stop("Model must contain at least one model object.")
  }

  # Tidy summary statistics for each model object and combine results
  mods_df <- purrr::map_dfr(model, function(mod) {
    # Check if mod is a valid model object
    if (!inherits(mod, "lm") && !inherits(mod, "glm")) {
      stop("Model must be a valid 'lm' or 'glm' object.")
    }
    if (!inherits(mod, "robust_glm")) {
      stop("model must be of class robust_glm.")
    }

    # Extract summary statistics from the model object
    summary_stats <- list(
      converged = tryCatch(mod$converged, error = function(e) NA),
      dispersion = tryCatch(
        if ("dispersion" %in% names(summary(mod))) summary(mod)$dispersion else NA,
        error = function(e) NA
      ),
      rmse = tryCatch(performance::performance_rmse(mod),
                      error = function(e) NA),
      mse = tryCatch(performance::performance_mse(mod),
                     error = function(e) NA),
      hosmer_chisq = tryCatch(performance::performance_hosmer(mod)$chisq,
                              error = function(e) NA),
      hosmer_df = tryCatch(performance::performance_hosmer(mod)$df,
                           error = function(e) NA),
      hosmer_p.value = tryCatch(performance::performance_hosmer(mod)$p.value,
                                error = function(e) NA),
      R2_Tjur = tryCatch(performance::r2_tjur(mod)[[1]],
                         error = function(e) NA),
      Log_loss = tryCatch(performance::performance_logloss(mod)[[1]],
                          error = function(e) NA)
    )
    glance_df <- tryCatch(broom::glance(mod, ...), error = function(e) NA)
    if (!is.null(glance_df)) summary_stats <- c(as.list(glance_df),
                                                summary_stats)

    # Combine summary statistics into a data frame
    tibble::as_tibble(summary_stats)
  }, .id = id)

  return(mods_df)
}