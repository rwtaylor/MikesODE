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
    library(deSolve)
    L = parms[1]; N = parms[2]; mu = parms[3]; m = parms[4];
    dH_local       = 2*mu*(1-y[1]) + 2*m*(y[2]-y[1]) - y[1]/(2*N)
    dH_global      = 2*mu*(1-y[2]) - y[1]/(2*N*L)
    list(c(dH_local, dH_global))
  }

  trajectories <- function(generations = 3000, L = 5, N = 500, mu = 5e-09, Ntot0 = 1e5, m = 1) {
    library(dplyr)
    H0 = 4 * input$Ntot0 * input$mu
    yini <- c(H0,H0)
    out_m <- ode(
               times = 1:generations,
               y = yini,
               func = heterozygosities,
               parms = c(L, N, mu, m)
             )
    out_df <- data_frame(
                gens = rep(1:generations, 2),
                type = rep(c("within", "total"), each = generations),
                heterozygocity = c(out_m[, 2], out_m[, 3])
              )
    return(out_df)
  }

  colors <- c("No migration" = rgb(7,54,248, maxColorValue = 255), "Full admixture" = rgb(0,0,0, maxColorValue = 255), "One mig per gen" = rgb(45,249,66, maxColorValue = 255))

  migration_calculator <- function(m, N ) {
    # Calculates migration rate from input text
    case_when(
      m == "No migration"    ~ 0,
      m == "Full admixture"  ~ 1,
      m == "One mig per gen" ~ 1/N
      )
  }

makePlot <- reactive ({
  })

output$distPlot <- renderPlot({

  migration_rates <- migration_calculator(input$migration, input$N)

  plot_data <- foreach(i = 1:length(migration_rates), .combine = rbind) %do% {
    out <- trajectories(input$generations, input$L, input$N, input$mu, input$Ntot0, migration_rates[i])
    out$migration = input$migration[i]
    out$m = migration_rates[i]
    out
  }

  plot_data$migration <- factor(plot_data$migration, levels = c("No migration","Full admixture","One mig per gen"))
  if(input$which_het == 1){
    plot_data <- plot_data %>% filter(type == "within")
  } else if(input$which_het == 2){
    plot_data <- plot_data %>% filter(type == "total")
  }

  p <- ggplot(plot_data, aes(x = gens, y = heterozygocity, color = migration, linetype = type, size = type)) + geom_line(alpha = 0.9) + xlab("Generation") + ylab("Heterozygosity")

  font_size = 24
  p <- p + theme_minimal(base_size = font_size)
  p <- p + scale_color_manual(values = colors, name = "Migration") + scale_linetype_manual(values = c("total" = "solid", "within" = "dotted"), name = "Heterozygosity") + scale_size_manual(values = c("total" = 3, "within" = 5), name = "Heterozygosity")
  p <- p + geom_hline(yintercept = input$h_critical)
  p <- p + geom_hline(yintercept = input$h_pop_decline)
#  p <- ggplotly(p)
  p
  
})
  
output$downloadPlot <- downloadHandler(
    filename = function(){input$plot_filename},
    content = function(file) {
        ggsave(file, plot = makePlot(), device = "pdf", width = 16, height = 10)
    }
  )  
})

