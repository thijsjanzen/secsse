
#' Logikelihood calculation for the SecSSE model given a set of parameters and data
#' @title Likelihood for SecSSE model
#' @param parameter list where first vector represents lambdas, the second mus and the third transition rates.
#' @param phy phylogenetic tree of class phylo, ultrametric, fully-resolved, rooted and with branch lengths.
#' @param traits vector with trait states, order of states must be the same as tree tips, for help, see vignette.
#' @param num_concealed_states number of concealed states, generally equivalent to number of examined states.
#' @param use_fortran Should the Fortran code for numerical integration be called? Default is TRUE.
#' @param methode Solver for the set of equations, default is "ode45".
#' @param cond condition on the existence of a node root: "maddison_cond","proper_cond"(default). For details, see vignette.
#' @param root_state_weight the method to weigh the states:"maddison_weights","proper_weights"(default) or "equal_weights". It can also be specified the root state:the vector c(1,0,0) indicates state 1 was the root state.
#' @param sampling_fraction vector that states the sampling proportion per trait state. It must have as many elements as trait states.
#' @param run_parallel should the routine to run in parallel be called?
#' @param setting_calculation argument used internally to speed up calculation. It should be left blank (default : setting_calculation = NULL)
#' @param setting_parallel argument used internally to set a parallel calculation. It should be left blank (default : setting_parallel = NULL)
#' @param see_ancestral_states should the ancestral states be shown? Deafault FALSE
#' @param loglik_penalty the size of the penalty for all parameters; default is 0 (no penalty)
#' @param is_complete_tree whether or not a tree with all its extinct species is provided
#' @param func function to be used in solving the ODE system Currently only for testing purposes.
#' @note To run in parallel it is needed to load the following libraries when windows: apTreeshape, doparallel and foreach. When unix, it requires: apTreeshape, doparallel, foreach and doMC
#' @return The loglikelihood of the data given the parameters
#' @examples
#' rm(list = ls(all = TRUE))
#' library(secsse)
#' library(DDD)
#' library(deSolve)
#' library(apTreeshape)
#' library(foreach)
#' set.seed(13)
#' phylotree <- ape::rcoal(31, tip.label = 1:31)
#' traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace = TRUE)
#' num_concealed_states <- 2
#' use_fortran <- TRUE
#' cond <- "proper_cond"
#' methode <- "ode45"
#' root_state_weight <- "proper_weights"
#' sampling_fraction <- c(1,1,1)
#' run_parallel <- FALSE
#' drill <- id_paramPos(traits,num_concealed_states)
#' drill[[1]][] <- c(0.12,0.01,0.2,0.21,0.31,0.23)
#' drill[[2]][] <- 0
#' drill[[3]][,] <- 0.1
#' diag(drill[[3]]) <- NA
#' secsse_loglik(parameter = drill,phylotree,traits,num_concealed_states,
#'    use_fortran,methode,cond,root_state_weight,sampling_fraction,see_ancestral_states = FALSE)
#'
#' #[1] -113.1018
#' @export
secsse_loglik_cpp <- function(parameter,
                          phy,
                          traits,
                          num_concealed_states,
                          methode = "ode45",
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          sampling_fraction,
                          setting_calculation = NULL,
                          setting_parallel = NULL,
                          see_ancestral_states = FALSE,
                          loglik_penalty = 0,
                          is_complete_tree = FALSE) {
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]

  if (is.null(setting_calculation)) {
    check_input(traits,phy,sampling_fraction,root_state_weight,is_complete_tree)
    setting_calculation <- build_initStates_time(phy,traits,num_concealed_states,sampling_fraction,is_complete_tree,mus)
  }

  states <- setting_calculation$states
  forTime <- setting_calculation$forTime
  ances <- setting_calculation$ances

  if (num_concealed_states != round(num_concealed_states)) { # for testing
    d <- ncol(states) / 2
    new_states <- states[,c(1:sqrt(d),(d + 1):((d + 1) + sqrt(d) - 1))]
    new_states <- states[,c(1,2,3,10,11,12)]
    states <- new_states
  }

  ly <- ncol(states)
  d <- ncol(states) / 2

  calcul <- calThruNodes_cpp(ances,
                             states,
                             forTime,
                             lambdas,
                             mus,
                             Q)

  mergeBranch <- calcul$mergeBranch
  nodeM <- calcul$nodeM
  loglik <- calcul$loglik

  ## At the root
  mergeBranch2 <- (mergeBranch)
  if (is.numeric(root_state_weight)) {
    giveWeights <- root_state_weight / num_concealed_states
    weightStates <- rep(giveWeights,num_concealed_states)
  } else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }

    if (root_state_weight == "proper_weights") {
      weightStates <- (mergeBranch2/(lambdas * (1 - nodeM[1:d]) ^ 2))/sum((mergeBranch2 / (lambdas * (1 - nodeM[1:d]) ^ 2)))
    }

    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1/length(mergeBranch2),length(mergeBranch2))
    }
  }

  if (cond == "maddison_cond") {
    mergeBranch2 <-
      mergeBranch2 / sum(weightStates * lambdas * (1 - nodeM[1:d]) ^ 2)
  }

  if (cond == "proper_cond") {
    mergeBranch2 <- mergeBranch2 / (lambdas * (1 - nodeM[1:d]) ^ 2)
  }

  atRoot <- ((mergeBranch2) * (weightStates))

  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - penalty(pars = parameter,loglik_penalty = loglik_penalty)

  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states),]
    ancestral_states <- ancestral_states[,-(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states,LL = LL))
  } else {
    return(LL)
  }
}
