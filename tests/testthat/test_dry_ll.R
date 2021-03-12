cla_secsse_loglik_rhs <- function(t, y, parameter) {
  ly <- length(y)
  d <- ly/2
  Es <- y[1:d]
  Ds <- y[(d + 1):ly]
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  Q <- parameter[[3]]
  diag(Q) <- 0

  all_states <- cbind(Ds, Es)
  a <- cbind(all_states[, 2], all_states[, 1])
  b <- t(all_states)
  cross_D_E <- a %*% b

  dD <- -((unlist(lapply(lambdas, sum))) +
            mus + Q %*% (rep(1, d))) *
    Ds + (Q %*% Ds) + unlist(lapply(lapply(lambdas, "*", cross_D_E), sum))
  dE <- -((unlist(lapply(lambdas, sum))) + mus + Q %*% (rep(1, d))) *
    Es +
    (Q %*% Es) +
    mus +
    unlist(lapply(lapply(lambdas, "*", Es %*% t(Es)), sum))

  return(list(c(dE, dD)))
}
phy <- NULL; rm(phy);
utils::data('example_phy_GeoSSE', package = 'secsse');
traits <- as.numeric(phy$tip.state)
lambdas <- list()
lambdas[[1]] <- matrix(0,ncol = 9,nrow = 9,byrow = TRUE)
lambdas[[1]][2,1] <- 1.5
lambdas[[1]][3,1] <- 0.5
lambdas[[1]][3,2] <- 1
for (i in 2:9) {
  lambdas[[i]] <- lambdas[[1]]
}
mus <- rep(0,9)
Q <- matrix(stats::runif(81),ncol = 9,nrow = 9,byrow = TRUE)
#diag(Q) <- NA
parameter <- list()
parameter[[1]] <- lambdas
parameter[[2]] <- mus
parameter[[3]] <- Q

num_concealed_states <- 3
sampling_fraction <- c(1,1,1)


lambdas <- parameter[[1]]
mus <- parameter[[2]]
parameter[[3]][is.na(parameter[[3]])] <- 0
Q <- parameter[[3]]

func <- cla_secsse_loglik_rhs

setting_calculation = NULL
root_state_weight = "maddison_weights"
is_complete_tree = FALSE
penalty = 0
see_ancestral_states = FALSE


if (is.null(setting_calculation)) {
  check_input(traits,
              phy,
              sampling_fraction,
              root_state_weight,
              is_complete_tree)
  setting_calculation <- build_initStates_time(phy,
                                               traits,
                                               num_concealed_states,
                                               sampling_fraction,
                                               is_complete_tree,
                                               mus)
}


states <- setting_calculation$states
forTime <- setting_calculation$forTime
ances <- setting_calculation$ances

write.table(states, file = "/Users/janzen/states_cla.txt", quote = F, col.names = F, row.names = F)
write.table(forTime, file = "/Users/janzen/forTime_cla.txt", quote = F, col.names = F, row.names = F)
write.table(Q, file = "/Users/janzen/Q_cla.txt", quote = F, col.names = F, row.names = F)
write.table(ances, file = "/Users/janzen/ances_cla.txt", quote = F, col.names = F, row.names = F)

if (num_concealed_states != round(num_concealed_states)) {
  # for testing
  d <- ncol(states)/2
  new_states <- states[, c(1:sqrt(d), (d + 1):((d + 1) + sqrt(d) - 1))]
  new_states <- states[, c(1, 2, 3, 10, 11, 12)]
  states <- new_states
}

loglik <- 0
ly <- ncol(states)
d <- ncol(states)/2

logliks <- c()

use_fortran <- FALSE
methode = "ode45"
# ances <- ances_vec
for (i in 1:length(ances)) {
  loglik <- 0
  calcul <- cla_calThruNodes(ances[i],
                             states,
                             loglik,
                             forTime,
                             parameter,
                             use_fortran = use_fortran,
                             methode = methode,
                             phy = phy,
                             func = func)
  states <- calcul$states
  loglik <- calcul$loglik
  nodeN <- calcul$nodeN
  logliks[i] <- loglik
}

logliks
sum(logliks)
mergeBranch <- calcul$mergeBranch
nodeM <- calcul$nodeM

cla_secsse_loglik_cpp(parameter = parameter,
                      phy = phy,
                      traits = traits,
                      num_concealed_states = num_concealed_states,
                      use_fortran = FALSE,
                      methode = "ode45",
                      cond = "maddison_cond",
                      root_state_weight = "maddison_weights",
                      sampling_fraction = sampling_fraction,
                      run_parallel = FALSE,
                      setting_calculation = NULL,
                      setting_parallel = NULL,
                      see_ancestral_states = FALSE,
                      loglik_penalty = 0,
                      is_complete_tree = FALSE)


lambdas <- parameter[[1]]
mus <- parameter[[2]]
parameter[[3]][is.na(parameter[[3]])] <- 0
Q <- parameter[[3]]


if (is.null(setting_calculation)) {
  check_input(traits,
              phy,
              sampling_fraction,
              root_state_weight,
              is_complete_tree)
  setting_calculation <- build_initStates_time(phy,
                                               traits,
                                               num_concealed_states,
                                               sampling_fraction,
                                               is_complete_tree,
                                               mus)
}


states <- setting_calculation$states
forTime <- setting_calculation$forTime
ances <- setting_calculation$ances

if (num_concealed_states != round(num_concealed_states)) {
  # for testing
  d <- ncol(states)/2
  new_states <- states[, c(1:sqrt(d), (d + 1):((d + 1) + sqrt(d) - 1))]
  new_states <- states[, c(1, 2, 3, 10, 11, 12)]
  states <- new_states
}

loglik <- 0
ly <- ncol(states)
d <- ncol(states)/2

calcul <- cla_calThruNodes_cpp(ances,
                               states,
                               forTime,
                               lambdas,
                               mus,
                               Q)

mergeBranch <- calcul$mergeBranch
nodeM <- calcul$nodeM




logliks <- c()

use_fortran <- FALSE
methode = "ode45"
# ances <- ances_vec
for (i in 1:length(ances)) {
  loglik <- 0
  calcul <- cla_calThruNodes(ances[i],
                             states,
                             loglik,
                             forTime,
                             parameter,
                             use_fortran = use_fortran,
                             methode = methode,
                             phy = phy,
                             func = func)
  states <- calcul$states
  loglik <- calcul$loglik
  nodeN <- calcul$nodeN
  logliks[i] <- loglik
}

logliks
sum(logliks)
mergeBranch2 <- calcul$mergeBranch
nodeM2 <- calcul$nodeM

all.equal(nodeM2, nodeM)
all.equal(mergeBranch2, mergeBranch)
