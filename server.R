library(tidyverse)
library(shiny)
library(ggplot2)
library(deSolve)
library(rsconnect)
library(ggthemes)
library(shinyBS)

shinyServer( function(input, output){

  ## Notation
  # HT = Total or global Heterozygosity
  # HS = Sub-population, local, or within-deme heterozygosity
  # R   = Fraction of alleles that are purged during a genetic rescue event
  # L   = Sub-population size
  # mu  = Mutation rate
  # N   = Current total population size
  # N0  = Total population size before anthropogenic bottleneck

  rescue_H <- function(RG, L, HS, HT, R){
    delta_HS = 2 * R * (HT - HS) - R^2 * (HT - HS)
    delta_HT = - HT * R^2 / L
    out <- tibble(generation = RG + 1, type = c("S", "T"),
                  H = c(HS + delta_HS, HT + delta_HT))
    return(out)
  }

  model_H <- function(times, y, p){
    dHS  = -y[1]/(2*p$N) + 2*p$mu*(1-y[1]) + 2*p$m*(y[2]-y[1])
    dHT  = -y[1]/(2*p$N*p$L) + 2*p$mu*(1-y[2])
    list(c(dHS,dHT))
  }

  forward_H <- function(generations = 1:input$generations, L = input$L,
    N = input$N, mu = input$mu, N0 = input$N0, m = 0, yini) {
    times <- generations - min(generations) + 1
    plist = list(L = L, N = N, mu = mu, m = m)
    future_H <- ode( times = times, y = yini, func = model_H, parms = plist)
    out <- tibble(generation = rep(generations, 2),
                  type = rep(c("S", "T"), each = max(times)),
                  H = c(future_H[, 2], future_H[, 3]))
    return(out)
  }


  # Scenario functions.
  H0 <- reactive({
    teta = 4 * input$N0 * input$mu
    teta / (1 + teta)
  })

  s_no_migration <- reactive({
    out <- forward_H(m = 0, yini = c(H0(), H0()))
    out$scenario = "No migration"
    out
  })

  s_full_admix <- reactive({
    out <- forward_H(m = 1, yini = c(H0(), H0()))
    out$scenario = "Full admixture"
    out
  })

  s_ompg <- reactive({
    out <- forward_H(m = 1/input$N, yini = c(H0(), H0()))
    out$scenario = "One mig per gen"
    out
  })

  last_H <- function(x) {
    x <- x %>% filter(generation == max(generation))
    HS <- x %>% filter(type == "S") %>% pull(H)
    HT <- x %>% filter(type == "T") %>% pull(H)
    return(c(HS, HT))
  }

  s_scheme_1 <- reactive({
    itt = 1
    teta = 4 * input$N0 * input$mu
    H0 = teta / (1 + teta)
    out_data <- forward_H(m = 0, yini = c(H0(), H0()))
    gen_h_crit <- out_data %>% filter(type == "S", H > input$h_scheme_1) %>% pull(generation) %>% max()
    out_data <- out_data %>% filter(generation < gen_h_crit)

    while(max(out_data$generation) < input$generations){
      itt = itt + 1
      H_R0 <- last_H(out_data) # Heterozygosities prior to rescue
      rescued_H <-   rescue_H(max(out_data$generation), input$L, H_R0[1], H_R0[2], input$R)
      decline_data <- forward_H(
        generations = max(rescued_H$generation)+1:input$generations,
        m = 0, yini = last_H(rescued_H))
      rescue_local_het <- rescued_H %>% filter(type == "S") %>% pull(H)
      if(rescue_local_het - 0.01 > input$h_scheme_1){
        gen_h_crit <- decline_data %>% filter(type == "S", H > input$h_scheme_1) %>% pull(generation) %>% max()
        decline_data <- decline_data %>% filter(generation < gen_h_crit)
      }
      out_data <- rbind(out_data, rescued_H, decline_data)
    }
    out_data <- out_data %>% filter(generation <= input$generations)
    print(paste("scheme 1 itterations:", itt))
    out_data$scenario = "Scheme 1"
    out_data
  })
  
  s_scheme_2 <- reactive({
    teta = 4 * input$N0 * input$mu
    H0 = teta / (1 + teta)
    out_data <- forward_H(m = 0, yini = c(H0(), H0()))
    gen_h_crit <- out_data %>% filter(type == "S", H > input$h_scheme_2) %>% pull(generation) %>% max()
    out_data <- out_data %>% filter(generation < gen_h_crit)
    itt = 1
    while(max(out_data$generation) < input$generations ){
      itt <- itt + 1
      N <- input$N
      migration_rates <- c(0, 1/(6*N), 1/(4*N), 1/(2*N), 1/N, 2/N, 3/N, 4/N)
      migration_period <- forward_H(generations = max(out_data$generation)+1:input$generations, m = migration_rates[itt], yini = last_H(out_data))
      max_het <- migration_period %>% filter(type == "S") %>% summarize(H = max(H))
      if(max_het - 0.01 > input$h_scheme_2 & itt < 8){
        gen_h_crit <- migration_period %>% filter(type == "S", H > input$h_scheme_2) %>% pull(generation) %>% max()
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
    out <- bind_rows(lapply(input$scenarios, scenario_switcher))
    out$scenario <- factor(out$scenario, levels = c("No migration","Full admixture","One mig per gen", "Scheme 1", "Scheme 2"))
    if(input$which_het == 1){
      out <- out %>% filter(type == "S")
    } else if(input$which_het == 2){
      out <- out %>% filter(type == "T")
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

    p <- ggplot(plot_data, aes(x = generation, y = H, color = scenario, linetype = type)) + geom_line(alpha = 0.9, size = 1) + xlab("Generation") + ylab("Heterozygosity")

    font_size = 24
    p <- p + theme_minimal(base_size = font_size)
    p <- p + scale_color_manual(values = colors, name = "Migration") + scale_linetype_manual(values = c("T" = "solid", "S" = "dashed"), name = "Heterozygosity") + scale_size_manual(values = c("T" = 3, "S" = 5), name = "Heterozygosity")
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
