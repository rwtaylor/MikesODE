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
            value = 100
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
            step = 1e-8,
            value = 1e-5
            )
        ),
        tabPanel("Scenarios",
          h3("Choose scenarios to plot"),
          checkboxGroupInput("scenarios",
          label = NULL,
            choices = list(
                        "No migration" = "No migration",
                        "Full admixture" = "Full admixture",
                        "One migrant per generation" = "One migrant per generation",
                        "Heterozygosity threshold" = "Heterozygosity threshold",
                        "Population decline (Not implemented yet)" = "Population decline"
                      ),
            selected = c("No migration", "Full admixture", "One migrant per generation")
          ),
          hr(),
          h3("Scenario parameters"),
          p("These parameters are shared across all applicable scenarios."),
          sliderInput("h_critical",
            "Critical threshold of heterozygosity",
            min = 0, max = 1, value = 0.2
          ),
          sliderInput('h_pop_decline',
            'Heterozygosity threshold for population decline',
            min = 0, max = 1, value = 0.1
          ),
          sliderInput('r_rate',
            'Rescue rate',
            min = 0, max = 1, value = 1
          ),
          sliderInput('r_frac',
            'Proportion of alleles replaced during rescue',
            min = 0, max = 1, value = 0.2
          )
        ), #tabpanel scenarios
        tabPanel("Plot Options",
          radioButtons("which_het", label = h3("Plot Heterozygosity"),
             choices = list("Local" = 1, "Global" = 2, "Both" = 3),
             selected = 3
          ),
          checkboxGroupInput("plot_thresholds",
            label = h3("Thresholds"),
            choices = list("Plot h critical" = "h_critical",
                           "Plot h population decline" = "h_pop_decline")
          )
        ) #tabpanel
      ) #tabsetpanel
    ), #column3
    column(9,
      plotOutput("distPlot", height = "800px")
    )
  ),
  fluidRow(
    column(3,
      wellPanel(
        h3("Save Plot"),
        textInput('plot_filename', "OM2M.pdf", label = "Filename"),
        downloadButton('downloadPlot', 'Download Plot')
        )
      )
    )#fluidRow
  )#fluidPage
)#shinyui
