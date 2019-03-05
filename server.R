library(tidyverse)
library(shiny)
library(ggplot2)
library(deSolve)
library(rsconnect)
library(ggthemes)
library(shinyBS)

# Define server logic required to draw a histogram
shinyServer( function(input, output) {

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
                type = rep(c("local", "global"), each = max(times)),
                heterozygosity = c(out_m[, 2], out_m[, 3])
              )
    return(out_df)
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
    hl <- x %>% filter(type == "local") %>% pull(heterozygosity)
    hg <- x %>% filter(type == "global") %>% pull(heterozygosity)
    return(c(hl, hg))
  }


  s_scheme_1 <- reactive({
    itt = 1
    teta = 4 * input$Ntot0 * input$mu
    H0 = teta / (1+teta)
    out_data <- trajectories(m = 0, rrate = 0, yini = c(H0, H0))
    gen_h_crit <- out_data %>% filter(type == "local", heterozygosity > input$h_scheme_1) %>% pull(generation) %>% max()
    out_data <- out_data %>% filter(generation < gen_h_crit)

    while(max(out_data$generation) < input$generations){
      itt = itt + 1
      rescue_gens <- max(out_data$generation):(max(out_data$generation) + 1)
      rescue_data <- trajectories(generations = rescue_gens, m = 0, rrate = input$r_rate, rfrac = input$r_frac, yini = last_hets(out_data))
      rescue_data <- rescue_data %>% filter(generation == max(generation))
      decline_data <- trajectories(
        generations = max(rescue_data$generation)+1:input$generations,
        m = 0, rrate = 0, yini = last_hets(rescue_data))
      rescue_local_het <- rescue_data %>% filter(type == "local") %>% pull(heterozygosity)
      if(rescue_local_het - 0.01 > input$h_scheme_1){
        gen_h_crit <- decline_data %>% filter(type == "local", heterozygosity > input$h_scheme_1) %>% pull(generation) %>% max()
        decline_data <- decline_data %>% filter(generation < gen_h_crit)
      }
      out_data <- rbind(out_data, rescue_data, decline_data)
    }
    out_data <- out_data %>% filter(generation <= input$generations)
    print(paste("scheme 1 itterations:", itt))
    out_data$scenario = "Scheme 1"
    out_data
  })
  
  s_scheme_2 <- reactive({
    teta = 4 * input$Ntot0 * input$mu
    H0 = teta / (1+teta)
    out_data <- trajectories(m = 0, rrate = 0, yini = c(H0, H0))
    gen_h_crit <- out_data %>% filter(type == "local", heterozygosity > input$h_scheme_2) %>% pull(generation) %>% max()
    out_data <- out_data %>% filter(generation < gen_h_crit)
    itt = 1
    while(max(out_data$generation) < input$generations ){
      itt <- itt + 1
      N <- input$N
      migration_rates <- c(0,1/(6*N),1/(4*N),1/(2*N),1/N,2/N,3/N,4/N)
      migration_period <- trajectories(generations = max(out_data$generation)+1:input$generations, m = migration_rates[itt], rrate = 0, yini = last_hets(out_data))
      max_het <- migration_period %>% filter(type == "local") %>% summarize(heterozygosity = max(heterozygosity))
      if(max_het - 0.01 > input$h_scheme_2 & itt < 8){
        gen_h_crit <- migration_period %>% filter(type == "local", heterozygosity > input$h_scheme_2) %>% pull(generation) %>% max()
        migration_period <- migration_period %>% filter(generation < gen_h_crit)
      }
      out_data <- rbind(out_data, migration_period)
    }
    out_data <- out_data %>% filter(generation <= input$generations)
    print(paste("scheme 2 itterations:", itt))
    out_data$scenario = "Scheme 2"
    out_data
  })


  scenario_switcher <- function(x){
    if(x == "no_migration") {
      s_no_migration()
    } else if (x == "full_admixture") {
      s_full_admix()
    } else if (x == "ompg") {
      s_ompg()
    } else if (x == "scheme_1") {
      s_scheme_1()
    } else if (x == "scheme_2") {
      s_scheme_2()
    }
  }

  get_plot_data <- reactive({
    H0 = 4 * input$Ntot0 * input$mu / (1 + 4 * input$Ntot0 * input$mu)
    out <- bind_rows(lapply(input$scenarios, scenario_switcher))
    out$scenario <- factor(out$scenario, levels = c("No migration","Full admixture","One mig per gen", "Scheme 1", "Scheme 2"))
    if(input$which_het == 1){
      out <- out %>% filter(type == "local")
    } else if(input$which_het == 2){
      out <- out %>% filter(type == "global")
    }
    out
  })

  colors <- c("No migration"    = rgb(0,0,255, maxColorValue = 255),
              "Full admixture"  = rgb(0,0,0, maxColorValue = 255),
              "One mig per gen" = rgb(0,255,0, maxColorValue = 255),
              "Scheme 1"        = rgb(255, 165, 0, maxColorValue = 255),
              "Scheme 2"        = rgb(160, 32, 240, maxColorValue = 255))


  output$distPlot <- renderPlot({
    plot_data <- get_plot_data()

    p <- ggplot(plot_data, aes(x = generation, y = heterozygosity, color = scenario, linetype = type)) + geom_line(alpha = 0.9, size = 1) + xlab("Generation") + ylab("Heterozygosity")

    font_size = 24
    p <- p + theme_minimal(base_size = font_size)
    p <- p + scale_color_manual(values = colors, name = "Migration") + scale_linetype_manual(values = c("global" = "solid", "local" = "dashed"), name = "Heterozygosity") + scale_size_manual(values = c("global" = 3, "local" = 5), name = "Heterozygosity")
    if (any(input$scenarios == "scheme_1")) {
      p <- p + geom_hline(yintercept = input$h_scheme_1, color = 'red', linetype = 'dotted')
    }
    if (any(input$scenarios == "scheme_2")) {
      p <- p + geom_hline(yintercept = input$h_scheme_2, color = 'red', linetype = 'dotted')
    }

    p

  }) #renderPlot

  output$downloadPlot <- downloadHandler(
      filename = function(){input$plot_filename},
      content = function(file) {
          ggsave(file, plot = output$distPlot, device = "pdf", width = 16, height = 10)
      }
  ) # downloadHandler
}) #shinyServer
