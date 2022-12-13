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
  sigma <- sqrare(run_model(tortoise, "tortoise")[[2]])

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

  ## 3. Generate B data (Y_1, ..., Y_B) from Poisson(mu_hat)
  Y_btsp <- tibble(B = 1:B) %>%
    group_by(B) %>%
    summarize(Y_hat = rpois(30, mu_hat[1:30,]), .groups = "drop") %>%
    nest_by(B) %>%
    summarize(glmer(Y_hat~prev+offset(log(Area))+factor(year)+(1|Site),
                   family=poisson,
                   data=data), .groups = "drop")





  Y_btsp



  ## Bootstrap standard deviation
  sd(bootstrap1$lambda_hat)

  ## Bootstrap confidence interval
  quantile(bootstrap1$lambda_hat,c(.025,.975))

  ## Bootstrap bias
  mean(bootstrap1$lambda_hat) - lambda_hat
}
