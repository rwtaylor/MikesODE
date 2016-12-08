library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("One Migrant Too Many"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("L",
                  "Number of demes:",
                  min = 1,
                  max = 50,
                  value = 5),
      sliderInput("N",
                  "Individuals per deme:",
                  min = 1,
                  max = 5000,
                  value = 50),
      sliderInput("Ntot0",
                  "Original population size:",
                  min = 1000,
                  max = 1e6,
                  value = 500*100),
      sliderInput("mu",
                  "Mutation rate:",
                  min = 0,
                  max = 1e-8,
                  value = 1e-8),
      sliderInput("generations",
                  "Generations:",
                  min = 10,
                  max = 1e5,
                  value = 2000),
      radioButtons("deriv",
                    label = h3("Derivation"),
                    choices = list("Mike's" = "mike", "Orren's" = "orren"),
                    selected = "mike")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))