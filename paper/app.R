library(shiny)
library(paper)

ui <- fluidPage(
  titlePanel("reactR HTMLWidget Example"),
  paperOutput('widgetOutput')
)

server <- function(input, output, session) {
  output$widgetOutput <- renderPaper(
    paper("Hello world!")
  )
}

shinyApp(ui, server)