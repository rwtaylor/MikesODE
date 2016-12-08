library(shiny)
library(ggplot2)
library(dplyr)
library(deSolve)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  ## Calculation of heterozygosity trajectories
  # Function that returns (dHs/dt, dHt/dt)
  # Assumptions: infinite alleles per site.
  # asserions: two lineages can coalesce only if they are in the same deme.
  # Derivation: backwards-coalescence frame of thought. 
  # Within-demes:
  # H(t+1) = Ht(1-1/2N-2m-2Miu) + 2mHgen(t) + 2miu(*), whose reordering gives dH/timestep,
  # As apears in the function below. [note that the case of coalescence doesn't have a
  # term of its own, because its contribution to H is zero]
  # * : If we assume an infinite allele model, any mutation always leads to a contribution 
  # to H. If we assume a 2-allele model with back mutation, we need to add the terms:
  # 2miu(1-Ht) - 2miu*Ht, where the first is a mutation in a pair that were identical, 
  # and the second a back-mutation in an initially non-identical pair.
  # Hgeneral (HG): we assume no two events occur in a single timestep, so ignore the option
  # of migration and coalescence or mig+mutation, and we even ignore the option that both
  # individuals happen to be from the same deme and one of them had a mutation (although
  # this last thing isn't completely justified, I think).
  # HG(t+1) = HG(t)(1-1/L-2mu) + 1/L * H(t) + 2mu
  # HG(t+1)-HG(t) = 1/L(H(t)-HG(t)) + 2mu(1-HG(t))
  # Which looks, admittedly, quite different from Mike's derivation which is based on a
  # forward-in-time-in-his-head intuition/derivation-that-I-didn't-follow.

  heterozygosities <- function(times, y, parms){
    library(deSolve)
    L = parms[1]; N = parms[2]; mu = parms[3]; m = parms[4];
    dH_within      = -y[1]/(2*N) + 2*m*(y[2]-y[1]) + 2*mu*(1-y[1])
    dH_generalMike = -y[1]/(2*N*L) + 2*mu*(1-y[2])
    dH_withinOren  = -y[3]/(2*N) + 2*m*(y[4]-y[3]) + 2*mu*(1-y[3])
    dH_generalOren = (y[3]-y[4])/L + 2*mu*(1-y[4])
    # list(c(dH_within,dH_generalMike))
    list(c(dH_within,dH_generalMike,dH_withinOren,dH_generalOren))
  }

  trajectories <- function(generations = 2000, L = 5, N = 50, mu = 1e-8, Ntot0 = 500*100, m = 0, deriv="mike") {
    library(dplyr)
    H0 = 4 * input$Ntot0 * input$mu
    yini <- c(H0,H0,H0,H0)
    out_m <- ode (times = 1:generations, y = yini, func = heterozygosities, parms = c(L, N, mu, m))

    if(deriv == "mike"){
      out_df <- data_frame(gens = rep(1:generations, 2), type = rep(c("within", "total"), each = generations),
                           heterozygocity = c(out_m[, 2], out_m[, 3]))
    } else {
      out_df <- data_frame(gens = rep(1:generations, 2), type = rep(c("within", "total"), each = generations),
                           heterozygocity = c(out_m[, 4], out_m[, 5]))
    }
    return(out_df)
  }


  #input <- data_frame(L = 5, N = 50, Ntot0 = 500*100, mu = 1e-8, generations = 2000)
  output$distPlot <- renderPlot({
    plot_data <- trajectories(input$generations, input$L, input$N, input$mu, input$Ntot0)
    ggplot(plot_data, aes(x = gens, y = heterozygocity, color = type)) + geom_line() + xlab("Generation") + ylab("Heterozygocity")
  })
})