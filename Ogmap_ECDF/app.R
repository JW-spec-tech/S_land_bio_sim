library(shiny)
library(ggplot2)
library(ggplotify)

# Define UI
ui <- fluidPage(
  titlePanel("Quadratic Equation Right Skewed"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      numericInput("a", "Enter a:", 1),
      numericInput("b", "Enter b:", 5),
      numericInput("c", "Enter c:", 150),
      sliderInput("amp", "Sin Amplitude:", min = 1, max = 50, value = 10),
      sliderInput("bw", "Bandwidth:", min = 0.01, max = 1, value = 0.2),
      sliderInput("steps", "Number of Steps:", min = 1, max = 50, value = 10),
      downloadButton("download", "Download JPG"),
    ),
    
    # Main layout with three plots
    mainPanel(
      plotOutput("plot1"),
      plotOutput("plot2"),
      plotOutput("plot3")
    )
  )
)

# Define server
server <- function(input, output) {
  # Create data frame for quadratic equation with sin
  quadratic_eqn <- reactive({
    x <- seq(-10, 10, by = 1)
    y <- input$a * x^2 + input$b * x + input$c + input$amp * sin(x)
    data.frame(x = x, y = y)
  })
  
  # First plot: quadratic equation right skewed with sin
  output$plot1 <- renderPlot({
    ggplot(quadratic_eqn(), aes(x, y)) +
      geom_point() +
      geom_smooth() +
      theme_minimal()
  })
  
  
  # Second plot: ecdf step plot
  output$plot2 <- renderPlot({
    ecdf_df <- data.frame(x = sort(quadratic_eqn()$y), y = ecdf(quadratic_eqn()$y)(sort(quadratic_eqn()$y)))
    ggplot(ecdf_df, aes(x, y)) +
      geom_step(aes(colour = "ECDF"), size = 1.2) +
      theme_minimal() +
      ggtitle("ECDF Step Plot")
  })
  
  # Third plot: ecdf step plot with adjustable bandwidth and number of steps
  output$plot3 <- renderPlot({
    bw <- input$bw
    n_steps <- input$steps
    ecdf_df <- data.frame(x = sort(quadratic_eqn()$y), y = ecdf(quadratic_eqn()$y)(sort(quadratic_eqn()$y)))
    ecdf_smooth <- approxfun(ecdf_df$x, ecdf_df$y, method = "linear")
    x_smooth <- seq(min(ecdf_df$x), max(ecdf_df$x), length.out = n_steps)
    y_smooth <- ecdf_smooth(x_smooth)
    ecdf_smooth_df <- data.frame(x_smooth, y_smooth)
    ggplot(ecdf_df, aes(x, y)) +
      geom_step(data = ecdf_smooth_df, aes(x = x_smooth, y = y_smooth), size = 1.2, colour = "red") +
      theme_minimal() +
      ggtitle("ECDF Step Plot with Adjustable Bandwidth and Steps")
  })
  
  output$download_plots <- downloadHandler(
    filename = function() {
      paste0("plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".jpg")
    },
    content = function(file) {
      p1 <- ggplot(quadratic_eqn(), aes(x, y)) +
        geom_point() +
        geom_smooth() +
        theme_minimal()
      
      p2 <- ggplot(ecdf_df, aes(x, y)) +
        geom_step(aes(colour = "ECDF"), size = 1.2) +
        theme_minimal() +
        ggtitle("ECDF Step Plot")
      
      p3 <- ggplot(ecdf_df, aes(x, y)) +
        geom_step(data = ecdf_smooth_df, aes(x = x_smooth, y = y_smooth), size = 1.2, colour = "red") +
        theme_minimal() +
        ggtitle("ECDF Step Plot with Adjustable Bandwidth and Steps")
      
      # Save plots as jpg files
      ggsave("plot1.jpg", plot = p1, device = "jpeg")
      ggsave("plot2.jpg", plot = p2, device = "jpeg")
      ggsave("plot3.jpg", plot = p3, device = "jpeg")
      
      # Convert jpg files to pdf
      magick::image_convert(c("plot1.jpg", "plot2.jpg", "plot3.jpg"), output = file, format = "jpg", quality = 90)
    },
    contentType = "image/jpeg"
  )
  
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
