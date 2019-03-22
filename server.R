library(tidyverse)
library(shiny)
library(rsconnect)
library(ggthemes)
library(shinyBS)

## Notation
# HT     = Total or global Heterozygosity
# HS     = Sub-population, local, or within-deme heterozygosity
# R      = Fraction of alleles that are purged during a genetic rescue event
# L      = Sub-population size
# mu     = Mutation rate
# N      = Current total population size
# Ntot0  = Total population size before anthropogenic bottleneck
# Time   = Generation or generations
# G_t    = Total number of generations to simulate
# lambda = Probability of random rare rescues per deme per gen (expected << 1)
colors <- c("No migration"    = rgb(0,0,255, maxColorValue = 255),
            "Full admixture"  = rgb(0,0,0, maxColorValue = 255),
            "One mig per gen" = rgb(0,255,0, maxColorValue = 255),
            "Scheme 1"        = rgb(255, 165, 0, maxColorValue = 255),
            "Scheme 2"        = rgb(160, 32, 240, maxColorValue = 255),
            "Random Rescue"   = rgb(255,192,203, maxColorValue = 255))


# "Standard model" = finite island model with mutation, migration, and drift,
# plus rare and random rescues with probability lambda per deme per generation.
# Assumes lambda << 1.
trajectory_standard <- function(scenario, Ntot0, G_t, L, N, mu, m, lambda = 0,
    R = 0.2) {
    #browser()
    # Note --- the approximation isn't technically valid for m~1 in the
    # discrete-time formulation, and here becomes numerically unstable for
    # m > 0.99.
    #
    # In general, need to consider the numerical stability issues in the
    # discrete sim's---I think has to do with HS becoming > than HT.
    
    theta0 = 4 * Ntot0 * mu
    HS0 = HT0 = theta0 / (1 + theta0)
    
    # Record results in a data frame that we allocate ahead of time. Note that
    # Time runs from 0 to G_t, with row i corresponding to Time = i - 1
    traj <- tibble(Time = seq(0, G_t), HS = NA, HT = NA, scenario = scenario)
    # Set H's at Time = 0 with the provided initial conditions
    traj[1,c("HS", "HT")] <- c(HS0, HT0)
    # Iterate over generations, calculating value in next generation from the
    # previous generation
    for (i in seq(G_t)) {
        HS <- traj[[i, "HS"]]
        HT <- traj[[i, "HT"]]
        # Standard WF island model w/ mutation, migration, and drift
        # + rare and random rescues
        dHS <- 2 * mu * (1 - HS) + 2 * m * (HT - HS) - HS / (2 * N) + 
            2 * lambda * R * (HT - HS) - lambda * R^2 * (2 * HT - HS)
        dHT <- 2 * mu * (1 - HT) - HS / (2*N*L) - lambda * R^2 * (2*HT - HS) / L 
        traj[i + 1, "HS"] <- HS + dHS
        traj[i + 1, "HT"] <- HT + dHT
    }
    traj
}

# Scheme 1 simulation, implemented as follows. Finite island model with
# mutation, migration, and drift. At start of each generation, check if HS <
# threshold. If so, do a simultaneous rescue on all demes. The extent to which
# rescues affect diversity is controlled by R.
trajectory_scheme1 <- function(scenario, Ntot0, G_t, L, N, mu, m = 0, R = 0.2, 
    threshold = 0.1) {
    theta0 = 4 * Ntot0 * mu
    HS0 = HT0 = theta0 / (1 + theta0)
    traj <- tibble(Time = seq(0, G_t), HS = NA, HT = NA, scenario = scenario)
    traj[1,c("HS", "HT")] <- c(HS0, HT0)
    for (i in seq(G_t)) {
        HS <- traj[[i, "HS"]]
        HT <- traj[[i, "HT"]]
        # Simultaneously rescue all demes if HS below the threshold
        if (HS < threshold) {
            dHS <- 2 * R * (HT - HS) - R^2 * (2*HT - HS)
            dHT <- - R^2 * HT / L
            HS <- HS + dHS
            HT <- HT + dHT
        }
        # Standard WF island model w/ mutation, migration, and drift
        dHS <- 2 * mu * (1 - HS) + 2 * m * (HT - HS) - HS / (2 * N)
        dHT <- 2 * mu * (1 - HT) - HS / (2*N*L)
        traj[i + 1, "HS"] <- HS + dHS
        traj[i + 1, "HT"] <- HT + dHT
    }
    traj
}

# Scheme 2 simulation, implemented as follows. Finite island model with
# mutation, migration, and drift. After update_interval, check if HS <
# threshold. If so, set the migration rate to a value designed to keep HS above
# the threshold. Do not raise migration rate above m_max.
trajectory_scheme2 <- function(scenario, Ntot0, G_t, L, N, mu, m = 0,
    threshold = 0.2, update_interval = 10, m_mult = 2, m_max = 0.1) {
    theta0 = 4 * Ntot0 * mu
    HS0 = HT0 = theta0 / (1 + theta0)
    traj <- tibble(Time = seq(0, G_t), HS = NA, HT = NA, scenario = scenario)
    traj[1,c("HS", "HT")] <- c(HS0, HT0)
#    migration_rates <- c(1/(6*N), 1/(4*N), 1/(2*N), 1/N, 2/N, 3/N)
    migration_rates <- c(1/(20*N), 1/(10*N), 1/(6*N), 1/(4*N), 1/(2*N), 1/N, 2/N, 3/N)
    m_itt <- 0
    for (i in seq(G_t)) {
        HS <- traj[[i, "HS"]]
        HT <- traj[[i, "HT"]]
        # Periodically check HS against the threshold and update the migration
        # rate.
        if (((i %% update_interval) == 0) & (HS < threshold) & m_itt < 8) {
            if (HT > threshold) {
                m_itt <- m_itt + 1
                # print(m_itt)
                #m <- m_mult / (4*N) * threshold / (HT - threshold)
                m <- migration_rates[m_itt]
                m <- min(m, m_max)
            }
        }
        # Standard WF island model w/ mutation, migration, and drift
        dHS <- 2 * mu * (1 - HS) + 2 * m * (HT - HS) - HS / (2 * N)
        dHT <- 2 * mu * (1 - HT) - HS / (2*N*L)
        traj[i + 1, "HS"] <- HS + dHS
        traj[i + 1, "HT"] <- HT + dHT
    }
    traj
}

shinyServer( function(input, output){

  scenario_switcher <- function(x){
    if(x == "no_migration") {
      trajectory_standard(scenario = "No migration", input$Ntot0, input$G_t, L = input$L, N = input$N, mu = input$mu, m = 0, lambda = 0)
    } else if (x == "full_admixture") {
      trajectory_standard(scenario = "Full admixture", input$Ntot0, input$G_t, L = input$L, N = input$N, mu = input$mu, m = 0.99, lambda = 0)
    } else if (x == "ompg") {
      trajectory_standard(scenario = "One mig per gen", input$Ntot0, input$G_t, L = input$L, N = input$N, mu = input$mu, m = 1/input$N, lambda = 0)
    } else if (x == "scheme_1") {
      trajectory_scheme1(scenario = "Scheme 1", input$Ntot0, input$G_t, L = input$L, N = input$N, mu = input$mu, R = input$R, threshold = input$h_scheme_1)
    } else if (x == "scheme_2") {
      trajectory_scheme2(scenario = "Scheme 2", input$Ntot0, input$G_t, L = input$L, N = input$N, mu = input$mu, threshold = input$h_scheme_2)
    } else if (x == "randomrescue") {
      trajectory_standard(scenario = "Random Rescue", input$Ntot0, input$G_t, L = input$L, N = input$N, mu = input$mu, m = 0, lambda = input$lambda)
    }
  }

  get_plot_data <- reactive({
    #browser()
    scenario_list <- input$scenarios
    if(input$randomrescue == TRUE){scenario_list = c(scenario_list, "randomrescue")}
    out <- map_dfr(scenario_list, scenario_switcher)
    out$scenario <- factor(out$scenario, levels = c("No migration","Full admixture", "One mig per gen", "Scheme 1", "Scheme 2", "Random Rescue"))
    if(input$which_het == 1){
      out <- out %>% filter(type == "S")
    } else if(input$which_het == 2){
      out <- out %>% filter(type == "T")
    }
    out
  })

  make_plot <- reactive({
    plot_data <- get_plot_data()
    plot_data <- gather(plot_data, key = "type", value = "H", HS:HT, -Time, -scenario)
    plot_data$type[plot_data$type == "HS"] <- "Within-deme H"
    plot_data$type[plot_data$type == "HT"] <- "Global H"

    p <- ggplot(plot_data, aes(x = Time, y = H, color = scenario, linetype = type)) +
      geom_line(alpha = 0.9, size = 1) +
      xlab("Generation") +
      ylab("Heterozygosity") + 
      theme_minimal(base_size = 16) +
      scale_color_manual(values = colors, name = "Migration") + 
      scale_linetype_manual(
          values = c("Global H" = "solid","Within-deme H" = "dashed"),
          name = "Heterozygosity"
          ) +
      theme(
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key.size = unit(16, "pt"),
        legend.box = "horizontal"
        ) + ylim(input$ylims)

    if (any(input$scenarios == "scheme_1")) {
      p <- p + geom_hline(yintercept = input$h_scheme_1,
        color = 'red', linetype = 'dotted')
    }
    if (any(input$scenarios == "scheme_2")) {
      p <- p + geom_hline(yintercept = input$h_scheme_2,
        color = 'red', linetype = 'dotted')
    }
    p
  })

  output$distPlot <- renderPlot({ make_plot() })

  output$downloadPlot <- downloadHandler(
      filename = function(){input$plot_filename},
      content = function(file) {
          ggsave(file, plot = make_plot(), device = "pdf", width = 16, height = 10)
      }
  ) # downloadHandler
}) #shinyServer
