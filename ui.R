library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("One Migrant Too Many"),


  # Sidebar with a slider input for the number of bins
  fluidRow(
    column(3,
      wellPanel(
        h3("Demographics"),
        sliderInput("L",
          "Number of demes:",
          min = 1,
          max = 50,
          value = 5),
        sliderInput("N",
          "Individuals (Ne) per deme:",
          min = 1,
          max = 500,
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
        checkboxGroupInput("migration",
          inline = TRUE,
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
          ),
        radioButtons("which_het", label = h3("Plot Heterozygosity"),
           choices = list("Within deme" = 1, "Total" = 2, "Both" = 3), 
           selected = 3
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
        selectInput("pstyle", label = h3("Plot Theme"),
          choices = list("Normal" = "grey", "BW" = "bw", "Minimal" = "minimal", "Tufte" = "tufte", "Base" = "base", "Classic" = "classic", "Linedraw" = "linedraw"),
          selected = "minimal")
        )
      ),
    column(3,
      wellPanel(
        h3("Save Plot"),
        textInput('plot_filename', "OM2M.pdf", label = "Filename"),
        downloadButton('downloadPlot', 'Download Plot')
        ))
    )
))