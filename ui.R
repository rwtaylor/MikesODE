library(shiny)
library(shinyBS)
# popify(bsButton("pointlessButton", "Button", style = "primary", size = "large"),
# +          "A Pointless Button",
# +          "This button is <b>pointless</b>. It does not do <em>anything</em>!"),
#style="display: flex;justify-content:space-between;,"
#  tags$head(
#        tags$style(type="text/css", "label{ display: inline}")
#      ),

info <- function(label = "Label", info_content = "add content") {
  div(
    p(label,
      popify(icon("info-circle"), title = NULL, content = info_content, placement = "right", trigger = "hover")
    )
  )
}

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Application title
  titlePanel("One Migrant Too Many"),
  # Sidebar with a slider input for the number of bins
  fluidRow(
    column(3,
      tabsetPanel(type = "tabs",
        tabPanel("Demographics",
          sliderInput("L", info("Number of demes:", "The number of isolated populations."),
            min   = 1,
            max   = 50,
            value = 5
            ),
          sliderInput("N",
            "Individuals (Ne) per deme:",
            min   = 1,
            max   = 10000,
            value = 100
            ),
          sliderInput("Ntot0",
            info("Original population size:", "The estimated population size prior to anthropogenic population decline."),
            min   = 1000,
            max   = 1e6,
            value = 1e5
            ),
          sliderInput("generations",
            info("Generations:", "The number of generations to run the simulation for"),
            min   = 10,
            max   = 5000,
            value = 3000
            ),
          numericInput("mu",
            info("Mutation rate:", "Microsattelite ~ 1e-5, SNP ~ 1e-8"),
            step = 1e-8,
            value = 1e-5
            )
        ),
        tabPanel("Scenarios",
          h3("Choose scenarios to plot"),
          checkboxGroupInput("scenarios",
          label = NULL,
            choiceNames = list(
              info("No migration", "All populations are isolated in perpetuity."),
              info("Full admixture", "The populations fully mix with each other."),
              info("One migrant per generation", "One individual per genearation migrates into each population."),
              info("Scheme 1","A genetic rescue occurs when the local heterozygosity drops below a set threshold."),
              info("Scheme 2", "Migration is allowed at increasing rates once heterozygosity drops to a set threshold.")
              ),
            choiceValues = list("no_migration", "full_admixture", "ompg",
                                "scheme_1", "scheme_2"),
            selected = c("no_migration", "full_admixture", "ompg",
                         "scheme_1", "scheme_2")
          ),
          hr(),
          h3("Scenario parameters"),
          sliderInput('h_scheme_1',
            'Heterozygosity threshold for Scheme 1',
            min = 0, max = 1, value = 0.1
          ),
          sliderInput("h_scheme_2",
            "Heterozygosity threshold for Scheme 2",
            min = 0, max = 1, value = 0.2
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
    ),
    HTML("<script>
$(document).ready(function(){
  $('[data-toggle=\\\"popover\\\"]').popover(); 
});
</script>")#fluidRow
  )#fluidPage
)#shinyui
