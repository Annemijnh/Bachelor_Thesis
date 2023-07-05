## RUN THIS before the ega output
plotClarkeGrid_adj <- function (referenceVals, testVals, title = "Clarke Error Grid", xlab="", ylab="",
                                linesize = 0.5, linetype = "solid", linecolor = "black",
                                linealpha = 0.6, pointsize = 4, pointalpha = 1, zones = NA, unit='gram')
{
  #unit selection and control
  if (unit != "mol" & unit != "gram") {
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  #axis named depending on unit
  if (unit == "mol") {
    n <- 18 #where n is a scaling factor for unit conversion
    if (xlab==""){
      xlab="Reference Glucose Concentration (mmol/L)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mmol/L)"
    }
  } else {
    n <- 1
    if (xlab==""){
      xlab="Reference Glucose Concentration (mg/dL)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mg/dL)"
    }
  }
  
  # use default zone assignment if none is provided
  if (is.na (zones)) {
    zones <- getClarkeZones (referenceVals, testVals, unit)
  }
  tolerance <- 0.2
  
  # create a df for ggplot (NULL to appease CRAN)
  ref <- test <- NULL
  data <- data.frame (ref=referenceVals, test=testVals, zones=zones)
  
  #better solution for scaling axis automatically with some extra space
  maxX <- max (max (data$ref) + 20 / n, 550 / n)
  maxY <- max ( (data$test + 20 / n), (1 + tolerance) * maxX, 650 / n)
  
  #labels with coordinats and colors
  labels <- data.frame (x=c (240 / n, 120 / n, 350 / n, 120 / n, 163 / n,
                             35 / n, 350 / n, 35 / n, 350 / n),
                        y=c (230 / n, 200 / n, 230 / n, 300 / n, 20 / n,
                             130 / n, 130 / n, 300 / n, 35 / n),
                        label=c ("A", "B", "B", "C", "C", "D", "D", "E", "E"),
                        color=c ("blue", "blue", "blue", "blue", "blue",
                                 "red", "red", "red", "red"))
  
  #segment endpoints for borders (NULL to appease CRAN)
  x1 <- y1 <- xend <- yend <- NULL
  border <- data.frame (x1=c (58.3 / n, 70 / n, 70 / n, 70 / n, 0, 240 / n,
                              0, 70 / n, 180 / n, 240 / n, 180 / n, 130 / n),
                        y1=c (70 / n, 56 / n, 180 / n, 83 / n, 180 / n,
                              180 / n, 70 / n, 0, 70 / n, 70 / n, 0, 0),
                        xend=c (maxX, maxX, 550 / n, 70 / n, 70 / n, maxX,
                                58.3 / n, 70 / n, maxX, 240 / n,
                                180 / n, 180 / n),
                        yend=c ((1 + tolerance) * maxX, (1 - tolerance) * maxX,
                                660 / n, maxY, 180 / n, 180 / n, 70 / n,
                                56 / n, 70 / n, 180 / n, 70 / n, 70 / n))
  
  ceg <- ggplot(data, aes(x = ref, y = test)) +
    scale_x_continuous (breaks = c (round (70 / n, digits=1),
                                    round (100 / n, digits=1),
                                    round (150 / n, digits=1),
                                    round (180 / n, digits=1),
                                    round (240 / n, digits=1),
                                    round (300 / n, digits=1),
                                    round (350 / n, digits=1),
                                    round (400 / n, digits=1),
                                    round (450 / n, digits=1),
                                    round (500 / n, digits=1),
                                    round (550 / n, digits=1),
                                    round (600 / n, digits=1),
                                    round (650 / n, digits=1),
                                    round (700 / n, digits=1),
                                    round (750 / n, digits=1),
                                    round (800 / n, digits=1),
                                    round (850 / n, digits=1),
                                    round (900 / n, digits=1),
                                    round (950 / n, digits=1),
                                    round (1000 / n, digits=1)),
                        expand = c (0, 0)) +
    scale_y_continuous (breaks = c (round (70 / n, digits=1),
                                    round (100 / n, digits=1),
                                    round (150 / n, digits=1),
                                    round (180 / n, digits=1),
                                    round (240 / n, digits=1),
                                    round (300 / n, digits=1),
                                    round (350 / n, digits=1),
                                    round (400 / n, digits=1),
                                    round (450 / n, digits=1),
                                    round (500 / n, digits=1),
                                    round (550 / n, digits=1),
                                    round (600 / n, digits=1),
                                    round (650 / n, digits=1),
                                    round (700 / n, digits=1),
                                    round (750 / n, digits=1),
                                    round (800 / n, digits=1),
                                    round (850 / n, digits=1),
                                    round (900 / n, digits=1),
                                    round (950 / n, digits=1),
                                    round (1000 / n, digits=1)),
                        expand = c (0, 0)) +
    geom_point(aes(color = zones), size = pointsize, alpha = pointalpha) +
    scale_color_manual(
      values=c("#e67c73","#ba67c8","#41b375","#7baaf7","#f7cb4d"),
      name="Zones",
      breaks= c("A", "B", "C", "D", "E"),
      labels = c("A", "B", "C", "D", "E"))+
    geom_segment (aes (x = x1, y = y1, xend = xend, yend = yend),
                  data=border, linetype=linetype) +
    annotate (geom="text", x = labels$x, y = labels$y, size = 16,
              label = labels$label, color=labels$color) +
    theme_bw (18) +
    theme (legend.position = "none") +
    ggtitle (title) +
    xlab (xlab) +
    ylab (ylab)+
    coord_cartesian(ylim=c(0,28),xlim=c(0,28))
  # theme(legend.position = "none")
  ceg
}
