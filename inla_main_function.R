# Sourcing the library
library(INLA)
devtools::install_github("MateusMaiaDS/gpbart") # Install if necessary
setwd(this.path::this.dir()) # Just a quick way to set the working directory
source("sim_function.R")
source("kfold_sim.R")

# Setting the cross-validation parameterisation
N <- 100
seed <- 42
tau <- 10
nu <- 0.1
phi <- 3

# Simulating the data
sim_dataset <- review_two_d_sim_two_node_gp_sum_each_tree(
    n = N,
    mu1 = c(-10,0,10),
    mu2 = c(5,20,-15),
    tau = tau,
    nu = nu,
    phi = phi,
    seed = seed,unif_sample = FALSE)[,1:3,drop = FALSE]

##----------------------------------------------------------------------
library(tidyverse)
sim_dataset2 <- as_tibble(sim_dataset)
ggplot(sim_dataset2, aes(x = lon, y = lat, size = y)) +
    geom_point(alpha = 0.7)
##----------------------------------------------------------------------

# Getting the cross-validation object
cross_validation_object <- k_fold(
    data = sim_dataset,
    dependent_variable = "y",
    k_partitions = 5,seed = 42)



# ==== This is usually commented, just used here to test te function
idx <- 5
x.train <- cross_validation_object[[idx]]$x_train
y.train <- cross_validation_object[[idx]]$y_train
x.test <- cross_validation_object[[idx]]$x_test
y.test <- cross_validation_object[[idx]]$y_test

inla_mod <- function(x.train,
                     y.train,
                     x.test,
                     y.test){
    ## Getting the train.df
    train.df <- cbind(x.train,y.train)
    test.df <- cbind(x.test,y.test)
    ## Scaling the data
    y_mean <- train.df[,"y"] %>% mean
    y_sd <- train.df[,"y"] %>% sd
    train.df[,"y"] <- (train.df[,"y"]-y_mean)/y_sd
    ## Getting the coordinates
    coo <- cbind(train.df[,"lat"],train.df[,"lon"])
    ## Creating the mesh
    ## Avoid mesh problems
    ## https://github.com/hrue/r-inla/issues/39
    n <- nrow(coo)
    coo2 <- coo + matrix(runif(2 * n) * 1e-9, n, 2)
    mesh <- inla.mesh.2d(loc = coo2, max.edge = 0.2)
    cat(">>> MESH DONE\n")
    ## Creating the SPDE
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
    ## Creating the indexs
    indexs <- inla.spde.make.index("s", spde$n.spde)
    ## Creating the Projector Matrix
    A <- inla.spde.make.A(mesh = mesh, loc = coo)
    ## Getting the new observations
    coop <- test.df[,c("lat","lon")]
    ## Creating the Projection Matrix for the new one
    Ap <- inla.spde.make.A(mesh = mesh,loc = coop)
    ## Creating the stacking
    stk.e <- inla.stack(
        tag = "est",
        data = list(y = train.df[,"y"]),
        A = list(1,A),
        effects = list(data.frame(
            b0 = 1,
            train.df[,!(colnames(train.df) %in% c("lat","lon","y"))]),
            s = indexs)
    )
    ## For the test data now we would have
    stk.p <- inla.stack(
        tag = "pred",
        data = list(y = NA),
        A = list(1,Ap),
        effects = list(data.frame(
            b0 = 1,
            test.df[,!(colnames(test.df) %in% c("lat","lon","y"))]),
            s = indexs)
    )
    ## Stacking both
    stk_full <- inla.stack(stk.e,stk.p)
    ## Creating the formula
    ## if(colna)
    cov_formula <- paste(colnames(
        train.df[,!(colnames(train.df) %in% c("lat","lon","y") )]),
        collapse = " + ")
    ## Restricting the case where there only spatial covariates
    if(cov_formula==""){
        formula <- y ~ 0 + b0 + f(s, model = spde)
    } else {
        formula <- as.formula(
            paste0("y ~ 0 + b0 +",cov_formula," + f(s,model = spde)")
        )
    }
    ## Calling the modelling function
    cat(sprintf(">>> Running model with formula `%s`\n",
                format(formula)))
    res <- inla(formula,
                family = "gaussian",
                data = inla.stack.data(stk_full),
                control.predictor = list(
                    compute = TRUE,
                    A = inla.stack.A(stk_full)),
                quantiles = c(0.25,0.75)
                )
    cat(">>> Model done\n")
    ## Getting the results from the test.prediction set we would have
    index <- inla.stack.index(stack = stk_full, tag = "pred")$data
    ## Getting the mean value
    prev_mean <- (res$summary.fitted.values[index,"mean"]*y_sd)+y_mean
    prev_ll <- res$summary.fitted.values[index, "0.25quant"]*y_sd + y_mean
    prev_ul <- res$summary.fitted.values[index, "0.75quant"]*y_sd + y_mean
    ## Prediction object
    prediction_object <- list(prev_mean = prev_mean,
                              prev_ll = prev_ll,
                              prev_ul = prev_ul)
    rmse <- sqrt(mean((prev_mean-test.df[,"y"])^2))
    pi <- mean((prev_ll <= test.df[,"y"]) & (prev_ul >= test.df[,"y"]))
    crps <- gpbart::crps(
                        y = y.test,
                        means = res$summary.fitted.values[index,"mean"],
                        sds = res$summary.fitted.values[index,"sd"])$CRPS
    ## Returning the list
    return(list(prediction_object = prediction_object,
                rmse = rmse,
                pi = pi,
                crps = crps))
}


##======================================================================
## Testing

cross_validation_object <- k_fold(
    data = sim_dataset,
    dependent_variable = "y",
    k_partitions = 5,seed = 42 + 1)

final <- vector(mode = "list", 5)

for(i in seq_along(final)) {
    x.train <- cross_validation_object[[i]]$x_train
    y.train <- cross_validation_object[[i]]$y_train
    x.test <- cross_validation_object[[i]]$x_test
    y.test <- cross_validation_object[[i]]$y_test
    cat(sprintf(">>>> Partition %s <<<<\n", i))
    final[[i]] <- try(inla_mod(x.train, y.train, x.test, y.test))
    ## stopifnot(class(final[[i]]) != "try-error")
}

any(sapply(final, class) == "try-error")
