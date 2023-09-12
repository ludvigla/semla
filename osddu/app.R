library(shiny)
library(osddu)

ui <- fluidPage(
  titlePanel("reactR HTMLWidget Example"),
  osdduOutput('widgetOutput')
)

server <- function(input, output, session) {
  output$widgetOutput <- renderOsddu(
    osddu("Hello world!")
  )
}

shinyApp(ui, server)