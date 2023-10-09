library(shiny)
library(ellipse)
library(circlize)
library(fields)
library(plotrix)

drawBoard = function(){
  fcts <- 1:20
  circlize::circos.par("gap.degree" = 0, "cell.padding" = c(0, 0, 0, 0),
                       start.degree = 360/40, track.margin = c(0, 0), "clock.wise" = FALSE)
  circlize::circos.initialize(factors = fcts, xlim = c(0, 1))
  circlize::circos.trackPlotRegion(ylim = c(0,1), factors = fcts, bg.col = "black",
                                   track.height = 0.15)
  circlize::circos.trackText(rep(0.5, 20), rep(0.5, 20),  rep(0.5, 20),
                             labels = c(13, 4, 18, 1, 20, 5, 12, 9, 14, 11, 8, 16, 7, 19, 3, 17, 2, 15, 10, 6),
                             factors = fcts, col = "#EEEEEE", font = 2,
                             facing = "downward")
  circlize::circos.trackPlotRegion(ylim = c(0, 1), factors = fcts,
                                   bg.col = rep(c("#df2623", "#11a551"), 10), bg.border = "#EEEEEE", 
                                   track.height = 0.04)
  circlize::circos.trackPlotRegion(ylim = c(0, 1), factors = fcts,
                                   bg.col = rep(c("black", "#e6cda5"), 10), bg.border = "#EEEEEE", 
                                   track.height = 0.275)
  circlize::circos.trackPlotRegion(ylim = c(0, 1), factors = fcts,
                                   bg.col = rep(c("#df2623", "#11a551"), 10), bg.border = "#EEEEEE",
                                   track.height = 0.04)
  circlize::circos.trackPlotRegion(ylim = c(0, 1), factors = fcts, 
                                   bg.col = rep(c("black", "#e6cda5"), 10), bg.border = "#EEEEEE",
                                   track.height = 0.415)
  circlize::draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360,
                        rou1 = 0.16/2, col = "#11a551", border = "#EEEEEE")
  circlize::draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360,
                        rou1 = 0.0635/2, col = "#df2623", border = "#EEEEEE")
  circlize::circos.clear()
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(HTML("<title>DartsOptimizeR</title>")),
  titlePanel(h1("Darts App", align = "center")),
  br(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h3("Throws"),
      textOutput("meanx"),
      textOutput("meany"),
      textOutput("nthrows"),
      # textOutput("score"),
      textOutput("max_score"),
      actionButton("sim", "Calculate optimal target")
    ),
    mainPanel(
      fluidRow(
        column(6,
               h3("Please enter your throws by clickling:"),
               tags$style(type="text/css", "img{display: flex; align-items: center; justify-content: center;}"),
               plotOutput("board", click = "plot_click"),
               br(), br(), br(), br(), br()),
        column(6, 
               h3("Optimal target (Heatmap):"),
               plotOutput("expScore"))
        )
      )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  val <- reactiveValues(clickx = NULL, clicky = NULL, newx = NULL, newy = NULL, score = NULL, obs = NULL, sim = NULL, time = NULL,
                        clickx_sim = NULL, clicky_sim = NULL, obs_sim = NULL, cov.mat = NULL, mus = NULL, E = 0)
  
  observe({
    input$plot_click
    isolate({
      val$clickx <- c(val$clickx, input$plot_click$x)
      val$clicky <- c(val$clicky, input$plot_click$y)
      # val$newx = input$plot_click$x
      # val$newy = input$plot_click$y
      val$clickx_sim <- c(val$clickx_sim, input$plot_click_sim$x)
      val$clicky_sim <- c(val$clicky_sim, input$plot_click_sim$y)    
      val$obs <- rbind(val$obs, c(input$plot_click$x, input$plot_click$y))
      val$obs_sim <- rbind(val$obs_sim, c(input$plot_click_sim$x, input$plot_click_sim$y))
    })
  }) # adding clicks to list
  
  # reactive({
  #   val$score = darts_score(x = input$newx*20, y = input$newy*20)
  # })
  
  output$board <- renderPlot({
    drawBoard()
    input$plot_click
    isolate({
      
      points(val$obs, col = "white", pch = 4, lwd = 3)
      mus <- tryCatch(colMeans(val$obs), error = function(e) NA)
      val$mus = mus
      x.star <- val$obs
      
      if(sum(is.na(mus)) == 0){
        if(nrow(x.star) > 2){
        x.star <- apply(x.star, 2, function(x) x - mean(x))
        sigma.hat <- 1 / (nrow(x.star) - 1) * (t(x.star) %*% x.star)
        
        sigma.hat_rescaled <- 1 / (nrow(x.star) - 1) * (t(x.star*200) %*% x.star*200)
        val$cov.mat <- sigma.hat_rescaled
        
        sds <- sqrt(diag(sigma.hat))
        diag.elements <- diag(sds)
        cor.matrix <- solve(diag.elements) %*% sigma.hat %*% solve(diag.elements)
        
        points(mus[1], mus[2], pch = 16, lwd = 10, col = "white")
        points(ellipse(cor.matrix, scale = sds, centre = mus, level = 0.1), type = "l", col = "white", lwd = 2.5)
        points(ellipse(cor.matrix, scale = sds, centre = mus, level = 0.3), type = "l", col = "white", lwd = 2.5)
        points(ellipse(cor.matrix, scale = sds, centre = mus, level = 0.5), type = "l", col = "white", lwd = 2.5)
        points(ellipse(cor.matrix, scale = sds, centre = mus, level = 0.7), type = "l", col = "white", lwd = 2.5)
        points(ellipse(cor.matrix, scale = sds, centre = mus, level = 0.9), type = "l", col = "white", lwd = 2.5)
        
        obs.so.far <- val$obs
        val$time <- as.data.frame(obs.so.far)
        
        }}
    })
  }, height = 520, width = 520)
  
  output$expScore = renderPlot({
    input$sim
    isolate({
      if(sum(is.na(val$mus)) == 0){
      if(nrow(val$obs) > 2){
        E = darts::generalExpScores(c(val$cov.mat[1,1], val$cov.mat[2,2], val$cov.mat[1,2]))*3
        val$E = E
        R = 170
        fields::image.plot(-R:R, -R:R, E, axes = F, xlim = c(-R, R), xlab = "", ylab = "",
              ylim = c(-R, R), nlevel = 64)
        # col = heat.colors(50)
        draw.circle(0, 0, radius = 170, border = scales::alpha("black", 0.5), lwd = 2.5)
        draw.circle(0, 0, radius = 162, border = scales::alpha("black", 0.5), lwd = 2.5)
        draw.circle(0, 0, radius = 107, border = scales::alpha("black", 0.5), lwd = 2.5)
        draw.circle(0, 0, radius = 99, border = scales::alpha("black", 0.5), lwd = 2.5)
        draw.circle(0, 0, radius = 16, border = scales::alpha("black", 0.5), lwd = 2.5)
        draw.circle(0, 0, radius = 6.35, border = scales::alpha("black", 0.5), lwd = 2.5)
        t0 = pi/2 + 2*pi/40
        points(c(16*cos(t0),170*cos(t0)), c(16*sin(t0),170*sin(t0)), col = scales::alpha("black", 0.5), type="l" , lwd = 2)
        for (i in 1:19){
          t1 = t0 - i*2*pi/20
          points(c(16*cos(t1),170*cos(t1)), c(16*sin(t1),170*sin(t1)), col = scales::alpha("black", 0.5), type="l" , lwd = 2)
        }
      }}
    })
  }, height = 512, width = 530)
  
  output$meanx <- renderText({
    paste0("Horizontal: ", round(mean(val$clickx)*20, digits = 1))
  })
  
  output$meany <- renderText({
    paste0("Vertical: ", round(mean(val$clicky)*20, digits = 1))
  })
  
  output$nthrows = renderText({
    paste("Number of throws:", length(val$clickx))
  })
  
  # output$score = renderText({
  #   paste("Punktzahl:", val$score)
  # })
  
  output$max_score = renderText({
    paste("Max. expected score (three darts):", round(max(val$E),1))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

