library(tidyverse)
library(shiny)
library(ggplot2)
library(deSolve)
library(rsconnect)
library(ggthemes)
library(shinyBS)

## Notation
# HT = Total or global Heterozygosity
# HS = Sub-population, local, or within-deme heterozygosity
# R   = Fraction of alleles that are purged during a genetic rescue event
# L   = Sub-population size
# mu  = Mutation rate
# N   = Current total population size
# N0  = Total population size before anthropogenic bottleneck
# G   = Generation or generations

colors <- c("No migration"    = rgb(0,0,255, maxColorValue = 255),
            "Full admixture"  = rgb(0,0,0, maxColorValue = 255),
            "One mig per gen" = rgb(0,255,0, maxColorValue = 255),
            "Scheme 1"        = rgb(255, 165, 0, maxColorValue = 255),
            "Scheme 2"        = rgb(160, 32, 240, maxColorValue = 255))

rescue_H <- function(RG, L, HS, HT, R){
  dHS = 2 * R * (HT - HS) - R^2 * (HT - HS)
  dHT = - HT * R^2 / L
  out <- tibble(G = RG + 1, HS = HS + dHS, HT = HT + dHT)
  return(out)
}

model_H <- function(times, y, p){
  dHS  = -y[1] / (2 * p$N) + 2 * p$mu * (1 - y[1]) + 2 * p$m * (y[2] - y[1])
  dHT  = -y[1] / (2 * p$N * p$L) + 2 * p$mu *(1 - y[2])
  list(c(dHS, dHT))
}

forward_H <- function(G, G_t, L, N, mu, m = 0, yini) {
  plist = list(L = L, N = N, mu = mu, m = m)
  times = G:G_t - G
  future_H <- ode(times = times, y = yini, func = model_H, parms = plist)
  out <- tibble(G = G:G_t, HS = future_H[, 2], HT = future_H[, 3])
  return(out)
}

shinyServer( function(input, output){

  # Scenario functions.
  H0 <- reactive({
    teta = 4 * input$N0 * input$mu
    teta / (1 + teta)
  })

  s_no_migration <- reactive({
    out <- forward_H(G = 1, G_t = input$G_total, L = input$L, N = input$N,
      mu = input$mu, m = 0, yini = c(H0(), H0()))
    out$scenario = "No migration"
    out
  })

  s_full_admix <- reactive({
    out <- forward_H(G = 1, G_t = input$G_total, L = input$L, N = input$N,
      mu = input$mu, m = 1, yini = c(H0(), H0()))
    out$scenario = "Full admixture"
    out
  })

  s_ompg <- reactive({
    out <- forward_H(G = 1, G_t = input$G_total, L = input$L, N = input$N,
      mu = input$mu, m = 1/input$N, yini = c(H0(), H0()))
    out$scenario = "One mig per gen"
    out
  })

  s_scheme_1 <- reactive({
    out <- forward_H(G = 1, G_t = input$G_total, L = input$L, N = input$N,
      mu = input$mu, m = 0, yini = c(H0(), H0()))
    out <- out %>% filter(HS > input$h_scheme_1)
    i = 1
    while(max(out$G) < input$G_total){
      i = i + 1
      H_R0 <- tail(out, n = 1) # H prior to rescue
      rescued_H <- rescue_H(H_R0$G, input$L, H_R0$HS, H_R0$HT, input$R)
      forward_data <- forward_H(G = rescued_H$G + 1, G_t = input$G_total,
        L = input$L, N = input$N, mu = input$mu, m = 0,
        yini = c(rescued_H$HS, rescued_H$HT))
      # Only rescue (on next itteration) if rescue would bump HS at least 0.01 above target.
      if(rescued_H$HS - 0.01 > input$h_scheme_1){
        forward_data <- forward_data %>% filter(HS > input$h_scheme_1)
      }
      out <- rbind(out, rescued_H, forward_data)
    }
    print(paste("scheme 1 ierations:", i))
    out$scenario = "Scheme 1"
    out
  })
  
  s_scheme_2 <- reactive({
    out <- forward_H(G = 1, G_t = input$G_total, L = input$L, N = input$N,
      mu = input$mu, m = 0, yini = c(H0(), H0()))
    out <- out %>% filter(HS > input$h_scheme_2)
    N <- input$N
    m_rates <- c(0, 1/(6*N), 1/(4*N), 1/(2*N), 1/N, 2/N, 3/N, 4/N)
    i = 1
    while(max(out$G) < input$G_total){
      i <- i + 1
      H_M0 <- tail(out, n = 1) # H prior to migration period
      m_period <- forward_H(G = H_M0$G + 1, G_t = input$G_total,
        L = input$L, N = input$N, mu = input$mu, m = m_rates[i],
        yini = c(H_M0$HS, H_M0$HT))
      if(max(m_period$HS) - 0.01 > input$h_scheme_2 & i < length(m_rates)){
        m_period <- m_period %>% filter(HS > input$h_scheme_2)
      }
      out <- rbind(out, m_period)
    }
    print(paste("scheme 2 ierations:", i))
    out$scenario = "Scheme 2"
    out
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

  output$distPlot <- renderPlot({
    plot_data <- get_plot_data()
    plot_data <- gather(plot_data, key = "type", value = "H", HS:HT, -G, -scenario)
    plot_data$type[plot_data$type == "HS"] <- "Within-deme H"
    plot_data$type[plot_data$type == "HT"] <- "Global H"

    p <- ggplot(plot_data, aes(x = G, y = H, color = scenario, linetype = type)) + geom_line(alpha = 0.9, size = 1) + xlab("Generation") + ylab("Heterozygosity")

    font_size = 24
    p <- p + theme_minimal(base_size = font_size)
    p <- p + scale_color_manual(values = colors, name = "Migration") + scale_linetype_manual(values = c("Global H" = "solid", "Within-deme H" = "dashed"), name = "Heterozygosity") + scale_size_manual(values = c("T" = 3, "S" = 5), name = "Heterozygosity")
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
