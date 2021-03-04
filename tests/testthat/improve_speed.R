context("improve speed")


test_that("secsse new gives same result as secsse old", {

  set.seed(42)
  out <- DDD::dd_sim(pars = c(0.4,0.1,40), age = 15)
  phy <- out$tes
  traits <- sample(c(0,1),ape::Ntip(phy),replace = T)
  b <- c(0.04,0.04)  # lambda
  d <- rep(0,2)
  userTransRate <- 0.2 # transition rate among trait states
  num_concealed_states <- 2
  sampling_fraction <- c(1,1)
  toCheck <- secsse::id_paramPos(traits,num_concealed_states)
  toCheck[[1]][] <- b
  toCheck[[2]][] <- d
  toCheck[[3]][,] <- userTransRate
  diag(toCheck[[3]]) <- NA
  root_state_weight <- "maddison_weights"
  use_fortran <- TRUE
  methode <- "ode45"
  cond <- "noCondit"

  loglik1 <- as.numeric(secsse::secsse_loglik(parameter = toCheck,
                                              phy = phy,
                                              traits = traits,
                                              num_concealed_states = num_concealed_states,
                                              use_fortran = TRUE,
                                              methode = "ode45",
                                              cond = cond,
                                              root_state_weight = root_state_weight,
                                              sampling_fraction = sampling_fraction,
                                              is_complete_tree = TRUE,
                                              func = "secsse_runmod_ct")
  )
  loglik2 <- as.numeric(secsse::secsse_loglik(parameter = toCheck,
                                              phy = phy,
                                              traits = traits,
                                              num_concealed_states = num_concealed_states,
                                              use_fortran = TRUE,
                                              methode = "ode45",
                                              cond = cond,
                                              root_state_weight = root_state_weight,
                                              sampling_fraction = sampling_fraction)
  )
  # check that the likelihood for a specifically complete tree without extinct lineages with 0 extinction
  # is equal to the likelihood for a tree with extant species only and 0 extinction rate
  testthat::expect_equal(loglik1, loglik2)
})