# Main function sample
review_two_d_sim_two_node_gp_sum_each_tree <- function(n, # Number of observations
                                                       mu1 =  c(-10,0,10), # Mean of the first terminal node
                                                       mu2  = c(5,20,-15) , # Mean of the second terminal nodde
                                                       nu = 0.1, # Getting the \nu parameter
                                                       phi = 0.1, # Getting the \phi parameter
                                                       tau = 10,
                                                       unif_sample = FALSE,
                                                       seed = NULL # Setting the seed
){
  #
  # Setting the seed
  set.seed(seed)

  # Defining the kernel function structure
  omega_function <- function(x, x_star = NULL, nu, phi) {

    # Calculating the square matrix
    if (is.null(x_star)) {
      kernel_matrix <- (nu^-1) * exp(-1 / (2 * phi^2) * as.matrix(stats::dist(x))^2)  +
        diag(1e-8, nrow = nrow(x))
    } else {
      kernel_matrix <- (nu^-1) * exp(-1 / (2 * phi^2) * as.matrix(stats::dist(x, x_star) )^2)
    }

    # Getting the kernel matrix
    return(kernel_matrix)
  }


  # Generating the x axis
  if(unif_sample){
    x <- expand.grid(lat = stats::runif(n = round(sqrt(n)),min = -10,max = 10),
                     lon = stats::runif(n = round(sqrt(n)),min = -10,max = 10))
  } else {
    x <- expand.grid(lat = seq(-10,10,length.out = round(sqrt(n))),
                     lon = seq(-10,10,length.out = round(sqrt(n))))

  }
  # Creating the true response
  y <- numeric(nrow(x))

  # Getting observation from different tree rules
  tree_one_n1 <- length(which(x[,1]< x[,2]))
  tree_one_n2 <- length(which(x[,1]>= x[,2]))

  tree_two_n1 <- length(which(x[,1] < -x[,2]))
  tree_two_n2 <- length(which(x[,1] >= - x[,2]))

  tree_three_n1 <- length(which(x[,1] < 0))
  tree_three_n2 <- length(which(x[,1] >= 0))


  # ==== Sampling for the first treee ======

  # -- Sampling for the first node

  # Defining the variance
  var_tree_one_node_one <- omega_function(x = x[x[,1] < x[,2], ,drop =FALSE],
                                          nu = nu,
                                          phi = phi)
  # Defining the mean
  mean_tree_one_node_one <- rep(mu1[1],tree_one_n1)


  # Sampling from the first node
  y [ x[,1] < x[,2] ] <-  y [ x[,1] < x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_one_node_one,
                                                              Sigma = var_tree_one_node_one)



  # -- Sampling for the second node

  # Defining the variance
  var_tree_one_node_two <-  omega_function(x = x[x[,1] >= x[,2], ,drop =FALSE],
                                           nu = nu,
                                           phi = phi)
  # Defining the mean
  mean_tree_one_node_two <- rep(mu2[1],tree_one_n2)

  y [ x[,1] >= x[,2] ] <-  y [ x[,1] >= x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_one_node_two,
                                                                Sigma = var_tree_one_node_two)


  # ==== Sampling for the second treee ======

  # -- Sampling for the first node

  # Defining the variance
  var_tree_two_node_one <-  omega_function(x = x[x[,1] < - x[,2], ,drop =FALSE],
                                           nu = nu,
                                           phi = phi)
  # Defining the mean
  mean_tree_two_node_one <- rep(mu1[2],tree_two_n1)


  # Sampling from the first node
  y [ x[,1] < -x[,2] ] <-  y [ x[,1] < -x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_two_node_one,
                                                                Sigma = var_tree_two_node_one)



  # -- Sampling for the second node

  # Defining the variance
  var_tree_two_node_two <-  omega_function(x = x[x[,1] >= -x[,2], ,drop =FALSE],
                                           nu = nu,
                                           phi = phi)
  mean_tree_two_node_two <- rep(mu2[2],tree_two_n2)

  y [ x[,1] >= -x[,2] ] <-  y [ x[,1] >= -x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_two_node_two,
                                                                  Sigma = var_tree_two_node_two)


  # ==== Sampling for the third treee ======

  # -- Sampling for the first node

  # Defining the variance
  var_tree_three_node_one <-  omega_function(x = x[x[,1] < 0, ,drop =FALSE],
                                             nu = nu,
                                             phi = phi)
  # Defining the mean
  mean_tree_three_node_one <- rep(mu1[3],tree_three_n1)


  # Sampling from the first node
  y [ x[,1] < 0 ] <-  y [ x[,1] < 0 ] + MASS::mvrnorm(n = 1,mu = mean_tree_three_node_one,
                                                      Sigma = var_tree_three_node_one)



  # -- Sampling for the second node

  # Defining the variance
  var_tree_three_node_two <-  omega_function(x = x[x[,1] >= 0, ,drop =FALSE],
                                             nu = nu,
                                             phi = phi)
  # Defining the mean
  mean_tree_three_node_two <- rep(mu2[3],tree_three_n2)

  y [ x[,1] >= 0 ] <-  y [ x[,1] >= 0 ] + MASS::mvrnorm(n = 1,mu = mean_tree_three_node_two,
                                                        Sigma = var_tree_three_node_two)

  y_true <-y
  y <- y + stats::rnorm(n = length(y),mean = 0,sd = tau^(-1/2))

  return(as.matrix(data.frame(x, y,y_true)))

}

# New simulation scenario
sim_2d_two_node <- function(n, seed, sd_value){


  # Setting a seed

  # Creating x and y
  x <- replicate(stats::runif(n = n,min = -1,max = 1),n = 2)
  colnames(x) <- c("lon","lat")
  y <- numeric(n)

  # Getting the left and right nodes
  left_index <- which(x[,"lon"]<0)
  right_index <- which(x[,"lon"]>=0)

  y[left_index] <- exp( (-0.5*(x[left_index,"lon"]^2+x[left_index,"lat"]^2)))
  y[right_index] <- -exp( (-0.5*(x[right_index,"lon"]^2+x[right_index,"lat"]^2)))
  f <- y
  y <- f + stats::rnorm(sd = sd_value,n = n)


  sim_data <- cbind(x,y,f)
  return(sim_data)
}

