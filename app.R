library(shiny)
source('functions/opt_vax_rate.R')

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Otimização de aplicação de duas doses"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      numericInput(inputId = "VAX.INITIAL.STORAGE.NUM",
                   label = "Estoque inicial de doses",
                   value = 0),
      numericInput(inputId = "VAX.PRODUCTION.RATE",
                   label = "Produção de novas doses or dia",
                   value = 0),
      numericInput(inputId = "MAX.VAC.RATE",
                   label = "Capacidade máxima de aplicação por dia",
                   value = 0),
      sliderInput(inputId = "VAX.WINDOW.DAYS",
                  label = "Intervalo entre a primeira e a segunda dose, em semanas",
                  min = 3,
                  max = 12,
                  step = 1,
                  value = 3),
      numericInput(inputId = "MAX.TIME.DAYS",
                   label = "Tempo total da simulação, em dias",
                   value = 120),
      numericInput(inputId = "SECOND.VAX.LOSS.FRAC",
                   label = "Fração de pessoas que não tomam a segunda dose",
                   value = 0.1)
    ),
# VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE, MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC, MAX.TIME.DAYS

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "distPlot")

    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  
  output$distPlot <- renderPlot(with(reactiveValuesToList(input), {
    OPT.VAX.RATE <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS,
                                 SECOND.VAX.LOSS.FRAC, MAX.TIME.DAYS)
    
    plot_vac_schedule(OPT.VAX.RATE, VAX.INITIAL.STORAGE.NUM,
                      VAX.PRODUCTION.RATE, MAX.VAC.RATE, VAX.WINDOW.DAYS,
                      SECOND.VAX.LOSS.FRAC, MAX.TIME.DAYS) 

    }))
}

shinyApp(ui, server)

