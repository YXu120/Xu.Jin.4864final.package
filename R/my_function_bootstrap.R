#' Final Project Parametric Bootstrap
#'
#' @param data the data set for the example
#' @param example the name of example, in our case we focus on "tortoise"
#' @param B sample size of bootstrapping
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
  beta_0 <- run_model(tortoise, "tortoise")[[1]][[1]]
  beta_1 <- run_model(tortoise, "tortoise")[[1]][[2]]
  year_factor <- c(0, run_model(tortoise, "tortoise")[[1]][[3]], run_model(tortoise, "tortoise")[[1]][[4]])
  sigma <- sqrt(run_model(tortoise, "tortoise")[[2]])

  ## The steps of parametric bootstrap are as follows:
  ## 1. Generate random effect Z* from N(0, sigma^2)
  Z <- rnorm(n, mean = 0, sd = sigma)
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
    summarize(shells = rpois(30, mu_hat[1:30,]), .groups = "drop") %>%
    nest_by(B)

  ## 4. Refit with the sampled (Y_1, ..., Y_B) to approximate betas
  beta_1_star <- array(data = NA, dim = B)
  for(i in 1:B){
    data$shells <- unlist(Y_btsp$data[[i]])
    beta_1_star[i] <- run_model(data, "tortoise")[[1]][[2]]
  }

  ## Bootstrap standard deviation
  beta_sd <- sd(beta_1_star)

  ## Bootstrap 95% confidence interval
  beta_ci <- quantile(beta_1_star,c(.025,.975))

  ## Test statistics and p-values


  ## Return values
  list(point_estimation = beta_1,
       standard_deviation = beta_sd,
       `95%_confidence_interval` = beta_ci,
       p_value = test_stat)
}
