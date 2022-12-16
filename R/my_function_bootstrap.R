#' Final Project Parametric Bootstrap
#'
#' @param data the data set for the example
#' @param example the name of example, in our case we focus on "tortoise"
#' @param B sample size of bootstrap
#'
#' @return A list with summary statistics including point estimates of coefficients, standard errors, confidence intervals and p-value for hypothesis test.
#' @export
#'
#' @examples
#' ##### Tortoise #####
#' ## Load data
#' tortoise <- read_csv("gopher_tortoise.csv")
#'
#' ## Choose the desired sample size for bootstrap
#' B <- 1000
#'
#' ## Apply the para_boot function to the data set with selected B
#' para_boot(tortoise, "tortoise", B)
para_boot <- function(data, example, B){
  ## Number of rows of data
  n <- nrow(data)

  ## Obtain the results we will use after fitting the model from run_model function
  # suppressMessages is used to avoid the "isSingular" message
  model_summary <- suppressMessages(run_model(data, example))
  beta_0 <- model_summary[[1]][[1]]
  beta_1 <- model_summary[[1]][[2]]
  year_factor <- c(0, model_summary[[1]][[3]], model_summary[[1]][[4]])
  sigma <- sqrt(model_summary[[2]])

  ## The steps of parametric bootstrap are as follows:
  ## 1. Generate random effect Z* from N(0, sigma^2)
  Z <- rnorm(10, mean = 0, sd = sigma)
  Z <- rep(Z, each = 3)

  ## 2. Compute mu_hat = exp(beta_0 + log(x_area) + year_factor + beta_1 * x_prev + Z
  # x_area
  x_area <- data %>%
    select(Area)
  # x_prev
  x_prev <- data %>%
    select(prev)
  # year_factor
  year_factor <- rep(year_factor , times = 10)
  # mu_hat
  mu_hat <- exp(beta_0 + log(x_area) + year_factor + beta_1 * x_prev + Z)

  ## 3. Resample (Y_1, ..., Y_B) from Poisson(mu_hat) for each i, j
  Y_btsp <- tibble(B = 1:B) %>%
    crossing(data) %>%
    group_by(B) %>%
    mutate(shells = rpois(n, mu_hat[1:n,]))

  Y_btsp <- tibble(B = 1:B) %>%
    crossing(data) %>%
    mutate(shells = Y_btsp$shells) %>%
    nest_by(B)

  ## 4. Refit with the sampled (Y_1, ..., Y_B) to obtain estimated beta_1 and test statistics
  beta_btsp <- tibble(t = 1:B) %>%
    group_by(t) %>%
    summarize(beta_1_star = suppressMessages(run_model(Y_btsp$data[[t]], example))[[1]][[2]],
              test_statistic = suppressMessages(run_model(Y_btsp$data[[t]], example))[[4]][[2]])

  ## Bootstrap standard deviation
  beta_sd <- sd(beta_btsp$beta_1_star)

  ## Bootstrap 95% confidence interval
  beta_ci <- quantile(beta_btsp$beta_1_star,c(.025,.975))

  ## P-value for hypothesis test, where H0: beta_1 = 0
  p <- mean(abs(beta_btsp$beta_1_star) > abs(beta_btsp$test_statistic))

  ## Return values
  # Histogram of bootstrap sampled beta_1
  hist(beta_btsp$beta_1_star,
       breaks = 20,
       main = "Histogram of sampled beta_1",
       xlab = "beta_1")

  # Summary statistics
  list(point_estimation = beta_1,
       bootstrapped_samples = beta_btsp$beta_1_star,
       standard_deviation = beta_sd,
       `95% confidence interval` = beta_ci,
       p_value = p)
}
