library(shiny)
library(ggplot2)
library(dplyr)
library(deSolve)
library(foreach)
library(rsconnect)
library(ggthemes)

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
    # list(c(dH_within,dH_generalMike))
    list(c(dH_within,dH_generalMike))
  }

  trajectories <- function(generations = 2000, L = 5, N = 50, mu = 1e-8, Ntot0 = 500*100, m = 1, deriv="mike") {
    library(dplyr)
    H0 = 4 * input$Ntot0 * input$mu
    yini <- c(H0,H0)
    out_m <- ode(times = 1:generations, y = yini, func = heterozygosities, parms = c(L, N, mu, m))
    out_df <- data_frame(gens = rep(1:generations, 2), type = rep(c("within", "total"), each = generations),
                         heterozygocity = c(out_m[, 2], out_m[, 3]))
    return(out_df)
  }

  colors <- c("0" = "#000000", "1/1000N" = "#E69F00", "1/100N" = "#56B4E9", "1/10N" = "#009E73", "1/N" ="#F0E442", "1/10" = "#0072B2", "1/100" = "#D55E00", "1" = "#CC79A7")
  colors <- c("0" = "#000000", "1/1000N" = "#E69F00", "1/100N" = "#F0E442", "1/10N" = "#009E73", "1/N" ="#56B4E9", "1/10" = "#0072B2", "1/100" = "#D55E00", "1" = "#CC79A7")

  makePlot <- reactive ({
    plot_data <- foreach(m_chr = input$migration, .combine = rbind) %do% {
      m_temp <- gsub("1/N", "1/(input$N)", m_chr)
      m <- eval(parse(text = gsub("1/([0-9]+)N", "1/(\\1 * input$N)", m_temp)))
      out <- trajectories(input$generations, input$L, input$N, input$mu, input$Ntot0, m, input$deriv)
      out$migration = m_chr
      out
    }
    
    plot_data$migration <- factor(plot_data$migration, levels = c("0","1/1000N","1/100N","1/10N" ,"1/N" ,"1/10" ,"1/100","1"))
    if(input$which_het == 1){
      plot_data <- plot_data %>% filter(type == "within")
    } else if(input$which_het == 2){
      plot_data <- plot_data %>% filter(type == "total")
    }
    
    p <- ggplot(plot_data, aes(x = gens, y = heterozygocity, color = migration, linetype = type, size = type)) + geom_line(alpha = 0.9) + xlab("Generation") + ylab("Heterozygocity")
    
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

#  input <- list(L = 5, N = 50, Ntot0 = 500*100, mu = 1e-8, generations = 2000, migration = c("0", "1/N", "1/10N", "1"))
output$distPlot <- renderPlot({
  makePlot()
})
  
output$downloadPlot <- downloadHandler(
    filename = function(){input$plot_filename},
    content = function(file) {
        ggsave(file, plot = makePlot(), device = "pdf", width = 16, height = 10)
    }
  )  
})

