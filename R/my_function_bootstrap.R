#' Final Project Parametric Bootstrap
#'
#' @param data the data set for the example
#' @param example the name of example: one of "culcita", "ctsib", "epilepsy", or "tortoise"
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

  ## Obtain the fitted values of coefficients, and sigma by using run_model function
  beta_0 <- run_model(tortoise, "tortoise")[[1]][[1]]
  beta_1 <- run_model(tortoise, "tortoise")[[1]][[2]]
  sigma <- sqrare(run_model(tortoise, "tortoise")[[2]])

  ## Generate random effect Z* from N(0, sigma^2)
  bootstrap_Z <- tibble(B = 1:B) %>%
    group_by(B) %>%
    summarize(d = rnorm(n, mean = 0, sd = sigma), .groups = "keep")

  ## Bootstrap standard deviation
  sd(bootstrap1$lambda_hat)

  ## Bootstrap confidence interval
  quantile(bootstrap1$lambda_hat,c(.025,.975))

  ## Bootstrap bias
  mean(bootstrap1$lambda_hat) - lambda_hat
}
