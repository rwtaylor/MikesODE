library(shiny)
library(ggplot2)
library(dplyr)
library(deSolve)
library(foreach)
library(rsconnect)
library(ggthemes)




heterozygosities <- function(times, y, parms){
  library(deSolve)
  L <- parms[1]; N <- parms[2]; mu <- parms[3]; m <- parms[4];
  dH_within       <- -y[1]/(2*N) + 2*m*(y[2]-y[1]) + 2*mu*(1-y[1])
  dH_generalMike  <- -y[1]/(2*N*L) + 2*mu*(1-y[2])
  # list(c(dH_within,dH_generalMike))
  list(c(dH_within,dH_generalMike))
}

trajectories <- function(generations = 2000, L = 5, N = 50, mu = 1e-8, Ntot0 = 500*100, m = 1, deriv="mike") {
  library(dplyr)
  H0 <- 4 * input$Ntot0 * input$mu
  yini <- c(H0,H0)
  out_m <- ode(times = 1:generations, y = yini, func = heterozygosities, parms = c(L, N, mu, m))
  out_df <- data_frame(gens = rep(1:generations, 2), type = rep(c("within", "total"), each = generations),
                       heterozygosity = c(out_m[, 2], out_m[, 3]))
  return(out_df)
}

colors <- c("0" = "#000000", "1/1000N" = "#E69F00", "1/100N" = "#56B4E9", "1/10N" = "#009E73", "1/N" ="#F0E442", "1/10" = "#0072B2", "1/100" = "#D55E00", "1" = "#CC79A7")
colors <- c("0" = "#000000", "1/1000N" = "#E69F00", "1/100N" = "#F0E442", "1/10N" = "#009E73", "1/N" ="#56B4E9", "1/10" = "#0072B2", "1/100" = "#D55E00", "1" = "#CC79A7")

# Migration rates
input <- list()
input$migration <- c("0", "1/1000N", "1/100N", "1/10N", "1/N", "1")
input$generations <- 2000
input$L <- 15
input$N <- 20
input$mu <- 1e-8
input$Ntot0 <- 50000


plot_data <- foreach(m_chr = input$migration, .combine = rbind) %do% {
  m_temp <- gsub("1/N", "1/(input$N)", m_chr)
  m <- eval(parse(text = gsub("1/([0-9]+)N", "1/(\\1 * input$N)", m_temp)))
  out <- trajectories(input$generations, input$L, input$N, input$mu, input$Ntot0, m)
  out$migration <- m_chr
  out$m_eval <- m
  out
}


plot_gen <- plot_data %>% group_by(type, migration) %>% mutate(prop_het = heterozygosity / max(heterozygosity)) %>% filter(gens == 100)
ggplot(plot_gen, aes(x = log(m_eval) , y = prop_het, color = type)) + geom_point() + geom_line()

















  makePlot <- reactive ({
    
    plot_data$migration <- factor(plot_data$migration, levels = c("0","1/1000N","1/100N","1/10N" ,"1/N" ,"1/10" ,"1/100","1"))
    if(input$which_het == 1){
      plot_data <- plot_data %>% filter(type == "within")
    } else if(input$which_het == 2){
      plot_data <- plot_data %>% filter(type == "total")
    }
    
    p <- ggplot(plot_data, aes(x = gens, y = heterozygosity, color = migration, linetype = type, size = type)) + geom_line(alpha = 0.9) + xlab("Generation") + ylab("Heterozygocity")
    
    font_size = 24
    
    if(input$pstyle == "bw") {
      p <- p + theme_bw(base_size = font_size)
    } else if(input$pstyle == "minimal") {
      p <- p + theme_minimal(base_size = font_size)
    } else if(input$pstyle == "tufte") {
      p <- p + theme_tufte(base_size = font_size)
    } else if(input$pstyle == "base") {
      p <- p + theme_base(base_size = font_size)
    } else if(input$pstyle == "classic") {
      p <- p + theme_classic(base_size = font_size)
    } else if(input$pstyle == "linedraw") {
      p <- p + theme_linedraw(base_size = font_size)
    } else if(input$pstyle == "grey") {
      p <- p + theme_grey(base_size = font_size)
    }
    p + scale_color_manual(values = colors, name = "Migration") + scale_linetype_manual(values = c("total" = "solid", "within" = "dotted"), name = "Heterozygocity") + scale_size_manual(values = c("total" = 3, "within" = 5), name = "Heterozygocity")
  })

