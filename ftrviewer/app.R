library(shiny)
library(ftrviewer)

ui <- fluidPage(
  titlePanel("reactR HTMLWidget Example"),
  ftrviewerOutput('widgetOutput')
)

server <- function(input, output, session) {
  output$widgetOutput <- renderFtrviewer(
    ftrviewer("Hello world!")
  )
}

shinyApp(ui, server)