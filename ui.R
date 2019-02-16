library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Application title
  titlePanel("One Migrant Too Many"),
  # Sidebar with a slider input for the number of bins
  fluidRow(
    column(3,
      tabsetPanel(type = "tabs",
        tabPanel("Demographics",
          sliderInput("L",
            "Number of demes:",
            min   = 1,
            max   = 50,
            value = 10
            ),
          sliderInput("N",
            "Individuals (Ne) per deme:",
            min   = 1,
            max   = 10000,
            value = 500
            ),
          sliderInput("Ntot0",
            "Original population size:",
            min   = 1000,
            max   = 1e6,
            value = 1e5
            ),
          sliderInput("generations",
            "Generations:",
            min   = 10,
            max   = 5000,
            value = 3000
            ),
          numericInput("mu",
            "Mutation rate:",
            step = 1e-9,
            value = 5e-9
            )
        ),
        tabPanel("Scenarios",
          checkboxGroupInput("migration",
            label = "Choose scenarios to plot",
            choices = list(
                        "No migration" = "No migration",
                        "Full admixture" = "Full admixture",
                        "One mig per gen" = "One mig per gen",
                        "Heterozygosity Threshold" = "Heterozygosity Threshold (Not implemented yet)",
                        "Population Decline (Not implemented yet)"
                      ),
            selected = c("No migration", "Full admixture", "One mig per gen")
          ),
          hr(),
          sliderInput("h_critical",
            "Critical threshold of heterozygosity",
            min = 0, max = 1, value = 0.02
          ),
          hr(),
          sliderInput('h_pop_decline',
            'Heterozygosity threshold for population decline',
            min = 0, max = 1, value = 0.01
          ),
          hr(),
          sliderInput('r_proportion',
            'Proportion of alleles replaced during rescue',
            min = 0, max = 1, value = 0.2
          )
        ),
        tabPanel("Plot Options",
          radioButtons("which_het", label = h3("Plot Heterozygosity"),
             choices = list("Local" = 1, "Global" = 2, "Both" = 3), 
             selected = 3
          ),
          checkboxGroupInput("plot_thresholds",
            label = h3("Thresholds"),
            choices = list("Plot Thresholds" = TRUE)
          )
        )
      )
    ),
    column(9,
      plotOutput("distPlot")
    )
  ),
  fluidRow(
    column(3,
      wellPanel(
        h3("Save Plot"),
        textInput('plot_filename', "OM2M.pdf", label = "Filename"),
        downloadButton('downloadPlot', 'Download Plot')
        ))
    )
))