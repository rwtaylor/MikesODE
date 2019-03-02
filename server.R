library(shiny)
library(ggplot2)
library(plotly)
library(tidyverse)
library(deSolve)
library(foreach)
library(rsconnect)
library(ggthemes)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  heterozygosities <- function(times, y, parms){
    L = parms[1]; N = parms[2]; mu = parms[3]; m = parms[4]; rrate = parms[5]; rfrac = parms[6];
    dH_within  = -y[1]/(2*N) + 2*mu*(1-y[1]) + 2*m*(y[2]-y[1]) + 2*rrate*rfrac*(y[2]-y[1]) - rrate*(rfrac^2)*(2*y[2]-y[1]);
    dH_general = -y[1]/(2*N*L) + 2*mu*(1-y[2]) - (2*y[2]-y[1]) * rrate * (rfrac^2) / L
    list(c(dH_within,dH_general))
  }

  trajectories <- function(generations = 1:input$generations, L = input$L,
    N = input$N, mu = input$mu, Ntot0 = input$Ntot0, m = 0,
    rrate = input$r_rate, rfrac = input$r_frac, yini) {
    times <- generations - min(generations) + 1

    out_m <- ode(
      times = times,
               y = yini,
               func = heterozygosities,
               parms = c(L, N, mu, m, rrate, rfrac)
             )
    out_df <- tibble(
                generation = rep(generations, 2),
                type = rep(c("within", "total"), each = max(times)),
                heterozygosity = c(out_m[, 2], out_m[, 3])
              )
    return(out_df)
  }

  trajectories2 <- function(generations = 1:input$generations, L = input$L,
    N = input$N, mu = input$mu, Ntot0 = input$Ntot0, m = 0,
    rrate = input$r_rate, rfrac = input$r_frac, yini) {
    browser()
    times <- generations - min(generations) + 1
    out_m <- ode(
      times = times,
               y = yini,
               func = heterozygosities,
               parms = c(L, N, mu, m, rrate, rfrac)
             )
    out_df <- tibble(
                generation = rep(generations, 2),
                type = rep(c("within", "total"), each = max(times)),
                heterozygosity = c(out_m[, 2], out_m[, 3])
              )
    return(out_df)
  }


  colors <- c("No migration" = rgb(7,54,248, maxColorValue = 255),
              "Full admixture" = rgb(0,0,0, maxColorValue = 255),
              "One mig per gen" = rgb(45,249,66, maxColorValue = 255),
              "Het. threshold" = rgb(254, 173, 67, maxColorValue = 255))

  migration_calculator <- function(m, N ) {
    # Calculates migration rate from input text
    case_when(
      m == "No migration"    ~ 0,
      m == "Full admixture"  ~ 1,
      m == "One mig per gen" ~ 1/N
      )
  }
  # Scenario functions.

  s_no_migration <- reactive({
    teta = 4 * input$Ntot0 * input$mu
    H0 = teta / (1+teta)

    out <- trajectories(m = 0, rrate = 0, yini = c(H0, H0))
    out$scenario = "No migration"
    out
  })

  s_full_admix <- reactive({
    teta = 4 * input$Ntot0 * input$mu
    H0 = teta / (1+teta)
    out <- trajectories(m = 1, rrate = 0, yini = c(H0, H0))
    out$scenario = "Full admixture"
    out
  })

  s_ompg <- reactive({
    teta = 4 * input$Ntot0 * input$mu
    H0 = teta / (1+teta)
    out <- trajectories(m = 1/input$N, rrate = 0, yini = c(H0, H0))
    out$scenario = "One mig per gen"
    out
  })

  last_hets <- function(x) {
    x <- x %>% filter(generation == max(generation))
    hl <- x %>% filter(type == "within") %>% pull(heterozygosity)
    hg <- x %>% filter(type == "total") %>% pull(heterozygosity)
    return(c(hl, hg))
  }

  s_het_thresh <- reactive({

    teta = 4 * input$Ntot0 * input$mu
    H0 = teta / (1+teta)
    out_data <- trajectories(m = 0, rrate = 0, yini = c(H0, H0))
    gen_h_crit <- out_data %>% filter(type == "within", heterozygosity > input$h_critical) %>% pull(generation) %>% max()
    out_data <- out_data %>% filter(generation < gen_h_crit)

    while(max(out_data$generation) < input$generations){
      rescue_gens <- max(out_data$generation):(max(out_data$generation) + 1)
      rescue_data <- trajectories(generations = rescue_gens, m = 0, rrate = input$r_rate, rfrac = input$r_frac, yini = last_hets(out_data))
      rescue_data <- rescue_data %>% filter(generation == max(generation))
      decline_data <- trajectories(generations = max(out_data$generation)+1:input$generations, m = 0, rrate = 0, yini = last_hets(rescue_data))
      gen_h_crit <- decline_data %>% filter(type == "within", heterozygosity > input$h_critical) %>% pull(generation) %>% max()
      decline_data <- decline_data %>% filter(generation < gen_h_crit)
      out_data <- rbind(out_data, rescue_data, decline_data)
    }

    out_data$scenario = "Het. threshold"
    out_data
  })

  scenario_switcher <- function(x){
    if(x == "No migration") {
      s_no_migration()
    } else if (x == "Full admixture") {
      s_full_admix()
    } else if (x == "One migrant per generation") {
      s_ompg()
    } else if (x == "Heterozygosity threshold") {
      s_het_thresh()
    }
  }

  get_plot_data <- reactive({
    H0 = 4 * input$Ntot0 * input$mu / (1 + 4 * input$Ntot0 * input$mu)
    out <- bind_rows(lapply(input$scenarios, scenario_switcher))
    out$scenario <- factor(out$scenario, levels = c("No migration","Full admixture","One mig per gen", "Het. threshold", "Pop. decline"))
    if(input$which_het == 1){
      out <- out %>% filter(type == "within")
    } else if(input$which_het == 2){
      out <- out %>% filter(type == "total")
    }
    out
  })


output$distPlot <- renderPlot({
  plot_data <- get_plot_data()

  p <- ggplot(plot_data, aes(x = generation, y = heterozygosity, color = scenario, linetype = type)) + geom_line(alpha = 0.9, size = 1) + xlab("Generation") + ylab("Heterozygosity")

  font_size = 24
  p <- p + theme_minimal(base_size = font_size)
  p <- p + scale_color_manual(values = colors, name = "Migration") + scale_linetype_manual(values = c("total" = "solid", "within" = "dashed"), name = "Heterozygosity") + scale_size_manual(values = c("total" = 3, "within" = 5), name = "Heterozygosity")
  if (any(input$plot_thresholds == "h_critical")) {
    p <- p + geom_hline(yintercept = input$h_critical, color = 'red', linetype = 'dotted')
  }
  if (any(input$plot_thresholds == "h_pop_decline")) {
    p <- p + geom_hline(yintercept = input$h_pop_decline, color = 'red', linetype = 'dotted')
  }

  p

})

output$downloadPlot <- downloadHandler(
    filename = function(){input$plot_filename},
    content = function(file) {
        ggsave(file, plot = makePlot(), device = "pdf", width = 16, height = 10)
    }
  )
})
