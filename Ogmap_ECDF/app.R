library(shiny)
library(ggplot2)
library(ggplotify)
library(ggpubr)

# Define UI
ui <- fluidPage(
  titlePanel("OGMap ECDF Demo with adjustable steps/bandwidth"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      numericInput("a", "Enter a:", 1),
      numericInput("b", "Enter b:", 5),
      numericInput("c", "Enter c:", 150),
      sliderInput("amp", "Sin Amplitude:", min = 1, max = 50, value = 10),
      sliderInput("steps", "Number of Steps:", min = 1, max = 50, value = 10),
      downloadButton("download_plots", "Download JPG"),
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
    y <- input$a * x^2 + input$b * x + input$c + input$amp * sin(x) *10
    data.frame(x = x, y = y)
  })
  
  generate_plots <- function(quadratic_eqn_data) {
    # First plot: quadratic equation right skewed with sin
    p1 <- ggplot(quadratic_eqn_data, aes(x, y)) +
      geom_point() +
      geom_smooth() +
      theme_minimal() +
      ggtitle("Data points for plotting")
    
    # Second plot: ecdf step plot
    ecdf_df <- data.frame(x = sort(quadratic_eqn_data$y), y = ecdf(quadratic_eqn_data$y)(sort(quadratic_eqn_data$y)))
    p2 <- ggplot(ecdf_df, aes(x, y)) +
      geom_step(aes(colour = "ECDF"), size = 1.2) +
      theme_minimal() +
      ggtitle("ECDF Step Plot")
    
    # Third plot: ecdf step plot with adjustable bandwidth and number of steps
    bw <- input$bw
    n_steps <- input$steps
    ecdf_smooth <- approxfun(ecdf_df$x, ecdf_df$y, method = "linear")
    x_smooth <- seq(min(ecdf_df$x), max(ecdf_df$x), length.out = n_steps)
    y_smooth <- ecdf_smooth(x_smooth)
    ecdf_smooth_df <- data.frame(x_smooth, y_smooth)
    p3 <- ggplot(ecdf_df, aes(x, y)) +
      geom_step(data = ecdf_smooth_df, aes(x = x_smooth, y = y_smooth), size = 1.2, colour = "red") +
      theme_minimal() +
      ggtitle("ECDF Step Plot with Adjustable Bandwidth and Steps")
    
    list(p1 = p1, p2 = p2, p3 = p3)
  }
  
  output$plot1 <- renderPlot({
    plot_objects <- generate_plots(quadratic_eqn())
    plot_objects$p1
  })
  
  output$plot2 <- renderPlot({
    plot_objects <- generate_plots(quadratic_eqn())
    plot_objects$p2
  })
  
  output$plot3 <- renderPlot({
    plot_objects <- generate_plots(quadratic_eqn())
    plot_objects$p3
  })
  
  output$download_plots <- downloadHandler(
    filename = function() {
      paste0("plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".jpg")
    },
    content = function(file) {
      plot_objects <- generate_plots(quadratic_eqn())
      p1 <- plot_objects$p1
      p2 <- plot_objects$p2
      p3 <- plot_objects$p3
      
      # Combine the plots into a single plot
      combined_plot <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3)
      
      # Save the combined plot as a jpg file
      ggsave(file, plot = combined_plot, device = "jpeg", width = 7, height = 15, dpi = 300)
    },
    contentType = "image/jpeg"
  )
}

# Run the application
shinyApp(ui = ui, server = server)
