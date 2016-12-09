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
      sliderInput("generations",
        "Generations:",
        min = 10,
        max = 1e4,
        value = 2000),
      numericInput("mu",
        "Mutation rate:",
        step = 1e-9,
        value = 1e-8),
      radioButtons("deriv",
        label = h3("Derivation"),
        choices = list("Mike's" = "mike", "Oren's" = "oren"),
        selected = "mike"),
      checkboxGroupInput("migration",
        label = h3("Migration Rates"), 
        choices = list("0" = "0",
                        "1/1000N" = "1/1000N",
                        "1/100N" = "1/100N",
                        "1/10N" = "1/10N",
                        "1/N" = "1/N",
                        "1/10" = "1/10",
                        "1/100" = "1/100",
                        "1" = "1"),
        selected = c("0", "1")
        )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))