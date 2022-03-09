# custom ggplot theme
newggtheme <- theme(
  plot.title = element_text(hjust = 0, vjust = 1, size = rel(2.3)), 
  panel.background = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black"),
  plot.background = element_blank(), 
  text = element_text(color = "gray20", size = 10), 
  axis.text = element_text(size = rel(1)), 
  axis.text.x = element_text(color = "gray20", size = rel(1.2)), 
  axis.text.y = element_text(color = "gray20", size = rel(1.2)), 
  axis.title.x = element_text(size = rel(1.3), vjust = 0), 
  axis.title.y = element_text(size = rel(1.3), vjust = 1), 
  axis.ticks.y = element_blank(), 
  axis.ticks.x = element_blank(), 
  strip.text.x = element_text(size = rel(1.5)),
  legend.position = "none",
  legend.key=element_rect(fill=NA),
  legend.title = element_blank(),
  legend.text=element_text(size=rel(1.5)),
  plot.subtitle = element_text(color = "gray20", size = rel(1.0), face="italic"),
  plot.caption = element_text(color = "dodgerblue", size = rel(1.0))
)



# with legends
newggtheme_with_legends <- theme(
  plot.title = element_text(hjust = 0, vjust = 1, size = rel(2.3)), 
  panel.background = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black"),
  plot.background = element_blank(), 
  text = element_text(color = "gray20", size = 10), 
  axis.text = element_text(size = rel(1)), 
  axis.text.x = element_text(color = "gray20", size = rel(1.2)), 
  axis.text.y = element_text(color = "gray20", size = rel(1.2)), 
  axis.title.x = element_text(size = rel(1.3), vjust = 0), 
  axis.title.y = element_text(size = rel(1.3), vjust = 1), 
  axis.ticks.y = element_blank(), 
  axis.ticks.x = element_blank(), 
  strip.text.x = element_text(size = rel(1.5)),
  legend.position = "top",
  legend.key=element_rect(fill=NA),
  legend.title = element_blank(),
  legend.text=element_text(size=rel(1.5)),
  plot.subtitle = element_text(color = "gray20", size = rel(1.0), face="italic"),
  plot.caption = element_text(color = "dodgerblue", size = rel(1.0))
)


#sociomatrix plot ===========================================


matrixplot <- function(m, mylevs=NULL, lowcolor="white",highcolor="red1"){
  
  #make the df we will use for plotting
  m.dat <- reshape2::melt(m)
  m.dat <- data.frame(m.dat)
  m.dat <- m.dat[complete.cases(m.dat),] #removing NAs
  
  if(is.null(mylevs)) { mylevs = rownames(m)}
  
  #reorder the levels of the y-axis so plots properly
  m.dat$loser <- factor(m.dat$loser, levels=mylevs)
  m.dat$winner <- factor(m.dat$winner, levels = rev(mylevs))
  m.dat[m.dat == 0] <- NA
  
  
  #plot
  p1<-ggplot(m.dat, aes(loser, winner, fill = value)) + 
    geom_tile(colour="black", 
              size=0.5, stat="identity") + 
    geom_text(data=m.dat, aes(loser, winner, label = value), color="black", size=rel(5.5)) +
    scale_fill_gradient(low = lowcolor, high = highcolor, space = "Lab", na.value = "white", guide = "colourbar") +
    scale_x_discrete(expand = c(0, 0), position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    xlab("Loser") + 
    ylab("Winner") +
    theme(axis.text.x = element_text(vjust = 1),
          axis.text.y = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          axis.text = element_text(color="#3C3C3C", size=rel(1.1)),
          legend.position = "none"        
    ) 
  return(p1)
}


# Dichotomized Matrix

matrixplot0 <- function(m,mylevs=NULL,lowcolor="white",highcolor="firebrick3"){
  
  dcs <- m / (m + t(m))
  m.dat_dc <- reshape2::melt(dcs)
  colnames(m.dat_dc)[3]<-"DC"
  
  m1 <- get_di_matrix(m)
  m.dat1 <- reshape2::melt(m1)
  m.dat1 <- data.frame(m.dat1)
  m.dat1 <- merge(m.dat1, m.dat_dc)
  m.dat1$value <- ifelse(is.na(m.dat1$DC), NA, m.dat1$value)
  
  
  #reorder the levels of the y-axis so plots properly
  if(is.null(mylevs)){mylevs = rownames(m)}
  m.dat1$loser <- factor(m.dat1$loser, levels=mylevs)
  m.dat1$winner <- factor(m.dat1$winner, levels = rev(mylevs))
  
  m.dat1$DC[m.dat1$value == 0] <- NA
  m.dat1$value[m.dat1$value == 0] <- NA
  
  
  p2<-ggplot(m.dat1, aes(loser, winner)) + 
    geom_tile(colour="black", aes(fill=DC),
              size=0.5, stat="identity") + 
    geom_text(data=m.dat1, aes(loser, winner, label = value), color="black", size=rel(3.5)) +
    scale_fill_gradient(low = lowcolor, high = highcolor, space = "Lab", na.value = "white", guide = "colourbar") +
    scale_x_discrete(expand = c(0, 0), position='top') +
    scale_y_discrete(expand = c(0, 0)) +
    xlab("Loser") + 
    ylab("Winner") +
    theme(axis.text.x = element_text(vjust = 1),
          axis.text.y = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          axis.text = element_text(color="#3C3C3C", size=rel(1.1)),
          legend.position = "none"        
    ) 
  return(p2)
}




## I want to change font size

matrixplot_fontsize <- function(m, mylevs=NULL, lowcolor="white",highcolor="red1",fontsize=5){
  
  #make the df we will use for plotting
  m.dat <- reshape2::melt(m)
  m.dat <- data.frame(m.dat)
  m.dat <- m.dat[complete.cases(m.dat),] #removing NAs
  
  if(is.null(mylevs)) { mylevs = rownames(m)}
  
  #reorder the levels of the y-axis so plots properly
  m.dat$loser <- factor(m.dat$loser, levels=mylevs)
  m.dat$winner <- factor(m.dat$winner, levels = rev(mylevs))
  m.dat[m.dat == 0] <- NA
  
  
  #plot
  p1<-ggplot(m.dat, aes(loser, winner, fill = value)) + 
    geom_tile(colour="black", 
              size=0.5, stat="identity") + 
    geom_text(data=m.dat, aes(loser, winner, label = value), color="black", size=rel(fontsize)) +
    scale_fill_gradient(low = lowcolor, high = highcolor, space = "Lab", na.value = "white", guide = "colourbar") +
    scale_x_discrete(expand = c(0, 0), position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    xlab("Loser") + 
    ylab("Winner") +
    theme(axis.text.x = element_text(vjust = 1),
          axis.text.y = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          axis.text = element_text(color="#3C3C3C", size=rel(1.1)),
          legend.position = "none"        
    ) 
  return(p1)
}


# Dichotomized Matrix

matrixplot0_fontsize <- function(m,mylevs=NULL,lowcolor="white",highcolor="firebrick3",fontsize=5){
  
  dcs <- m / (m + t(m))
  m.dat_dc <- reshape2::melt(dcs)
  colnames(m.dat_dc)[3]<-"DC"
  
  m1 <- get_di_matrix(m)
  m.dat1 <- reshape2::melt(m1)
  m.dat1 <- data.frame(m.dat1)
  m.dat1 <- merge(m.dat1, m.dat_dc)
  m.dat1$value <- ifelse(is.na(m.dat1$DC), NA, m.dat1$value)
  
  
  #reorder the levels of the y-axis so plots properly
  if(is.null(mylevs)){mylevs = rownames(m)}
  m.dat1$loser <- factor(m.dat1$loser, levels=mylevs)
  m.dat1$winner <- factor(m.dat1$winner, levels = rev(mylevs))
  
  m.dat1$DC[m.dat1$value == 0] <- NA
  m.dat1$value[m.dat1$value == 0] <- NA
  
  
  p2<-ggplot(m.dat1, aes(loser, winner)) + 
    geom_tile(colour="black", aes(fill=DC),
              size=0.5, stat="identity") + 
    geom_text(data=m.dat1, aes(loser, winner, label = value), color="black", size=rel(fontsize)) +
    scale_fill_gradient(low = lowcolor, high = highcolor, space = "Lab", na.value = "white", guide = "colourbar") +
    scale_x_discrete(expand = c(0, 0), position='top') +
    scale_y_discrete(expand = c(0, 0)) +
    xlab("Loser") + 
    ylab("Winner") +
    theme(axis.text.x = element_text(vjust = 1),
          axis.text.y = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          axis.text = element_text(color="#3C3C3C", size=rel(1.1)),
          legend.position = "none"        
    ) 
  return(p2)
}


## Glickoplot==========================

plotglicko <- function(robj, mycolors=c("black", "grey", "orange", "red"), 
                       ltypes=c(1,2,3,1,2,3,1,2,3,1), thetitle="",  linewd=1, ylim1=1000,ylim2=3250,
                       ndays=1){
  
  x<-as.data.frame(unlist(robj$history))  
  n <- length(x)/4
  
  x.ratings<-x[,1:n]
  x.deviations<-x[,(1+n):(n+n)]
  
  #longform the data
  x.ratingsmelt<-reshape2::melt(x.ratings)
  
  ids<-rownames(x.ratings)       #to make id column
  x.ratingsmelt$ids<-rep(ids, n)  #making id column
  
  l.ids<-length(ids)
  x.ratingsmelt$event<-rep(1:n, each=l.ids) 
  
  
  #add ranks
  xrn<-as.data.frame(x.ratings[n])
  colnames(xrn)<-c("finalrating")
  
  x.ratingsmelt$rank<-rank(-xrn$finalrating, ties.method="random")
  
  #make ids1 a factor with levels defined by rank
  x.ratingsmelt1 <- data.frame(ids=unique(x.ratingsmelt$ids),rank=unique(x.ratingsmelt$rank))
  x.ratingsmelt1 <- x.ratingsmelt1[order(x.ratingsmelt1$rank),]
  x.ratingsmelt$ids1 <- factor(x.ratingsmelt$ids,levels=x.ratingsmelt1$ids)
  
  #define color palette, multiple options (see below)
  colourCount <-length(unique(x.ratingsmelt$ids))
  getPalette = colorRampPalette(mycolors)
  
  
  ### now plot using ids1 instead of ids.
  p1<-ggplot(x.ratingsmelt, aes(x = event, y = value, col=ids1, linetype=ids1)) +
    scale_colour_manual(values = getPalette(colourCount)) +
    scale_linetype_manual(values=ltypes) +
    ylab("Glicko Rating") +
    xlab("Event") +
    ggtitle(thetitle)+
    ylim(ylim1,ylim2)+
    geom_line(lwd = linewd) +
    theme(plot.title = element_text(hjust = 0, vjust = 1, size = rel(1.7)), 
          panel.background = element_blank(), 
          plot.background = element_blank(), 
          panel.grid.major.y = element_line(color = "gray75",linetype = 'dotted'), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank(), 
          strip.background  = element_blank(),
          strip.text = element_text(size=rel(1.1)),
          text = element_text(color="gray20", size=10),
          axis.text = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="gray20", size=rel(1.0)),
          axis.text.y = element_text(color="gray20", size=rel(1.0)),
          axis.title.x = element_text(size=rel(1.0), vjust=0),
          axis.title.y = element_text(size=rel(1.0), vjust=1),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  
  return(p1)
}




### Making half box/half point plot ========================================================

# somewhat hackish solution to:
# https://twitter.com/EamonCaddigan/status/646759751242620928
# based mostly on copy/pasting from ggplot2 geom_violin source:
# https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r

library(ggplot2)
library(dplyr)


"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL 
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

######################################################
#source("https://raw.githubusercontent.com/tidyverse/ggplot2/2f3fef72e140d34210daa9d95917c77b19e89669/R/geom-boxplot.r")
library(grid)

ggname <- function (prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

geom_boxplot2 <- function(mapping = NULL, data = NULL,
                          stat = "boxplot", position = "dodge",
                          ...,
                          outlier.colour = NULL,
                          outlier.color = NULL,
                          outlier.fill = NULL,
                          outlier.shape = 19,
                          outlier.size = 1.5,
                          outlier.stroke = 0.5,
                          outlier.alpha = NULL,
                          notch = FALSE,
                          notchwidth = 0.5,
                          varwidth = FALSE,
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE) {
  
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomBoxplot2,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      outlier.colour = outlier.color %||% outlier.colour,
      outlier.fill = outlier.fill,
      outlier.shape = outlier.shape,
      outlier.size = outlier.size,
      outlier.stroke = outlier.stroke,
      outlier.alpha = outlier.alpha,
      notch = notch,
      notchwidth = notchwidth,
      varwidth = varwidth,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomBoxplot2 <- ggproto("GeomBoxplot2", Geom,
                        setup_data = function(data, params) {
                          data$width <- data$width %||%
                            params$width %||% (resolution(data$x, FALSE) * 0.9)
                          
                          if (!is.null(data$outliers)) {
                            suppressWarnings({
                              out_min <- vapply(data$outliers, min, numeric(1))
                              out_max <- vapply(data$outliers, max, numeric(1))
                            })
                            
                            data$ymin_final <- pmin(out_min, data$ymin)
                            data$ymax_final <- pmax(out_max, data$ymax)
                          }
                          
                          # if `varwidth` not requested or not available, don't use it
                          if (is.null(params) || is.null(params$varwidth) || !params$varwidth || is.null(data$relvarwidth)) {
                            data$xmin <- data$x - data$width / 2
                            data$xmin2 <- data$x - data$width / 4
                            data$xmax <- data$x + data$width / 2
                          } else {
                            # make `relvarwidth` relative to the size of the largest group
                            data$relvarwidth <- data$relvarwidth / max(data$relvarwidth)
                            data$xmin <- data$x - data$relvarwidth * data$width / 2
                            data$xmin2 <- data$x - data$relvarwidth * data$width / 4
                            data$xmax <- data$x + data$relvarwidth * data$width / 2
                          }
                          data$width <- NULL
                          if (!is.null(data$relvarwidth)) data$relvarwidth <- NULL
                          
                          data
                        },
                        
                        draw_group = function(data, panel_params, coord, fatten = 2,
                                              outlier.colour = NULL, outlier.fill = NULL,
                                              outlier.shape = 19,
                                              outlier.size = 1.5, outlier.stroke = 0.5,
                                              outlier.alpha = NULL,
                                              notch = FALSE, notchwidth = 0.5, varwidth = FALSE) {
                          
                          common <- data.frame(
                            colour = data$colour,
                            size = data$size,
                            linetype = data$linetype,
                            fill = alpha(data$fill, data$alpha),
                            group = data$group,
                            stringsAsFactors = FALSE
                          )
                          
                          whiskers <- data.frame(
                            x = c(data$x,data$x,data$xmin2,data$xmin2),
                            xend = c(data$x,data$x,data$x,data$x),
                            y = c(data$upper, data$lower,data$ymax,data$ymin),
                            yend = c(data$ymax, data$ymin,data$ymax,data$ymin),
                            alpha = NA,
                            common,
                            stringsAsFactors = FALSE
                          )
                          
                          
                          box <- data.frame(
                            xmin = data$xmin,
                            xmax = data$x,
                            ymin = data$lower,
                            y = data$middle,
                            ymax = data$upper,
                            ynotchlower = ifelse(notch, data$notchlower, NA),
                            ynotchupper = ifelse(notch, data$notchupper, NA),
                            notchwidth = notchwidth,
                            alpha = data$alpha,
                            common,
                            stringsAsFactors = FALSE
                          )
                          
                          if (!is.null(data$outliers) && length(data$outliers[[1]] >= 1)) {
                            outliers <- data.frame(
                              y = data$outliers[[1]],
                              x = data$x[1],
                              colour = outlier.colour %||% data$colour[1],
                              fill = outlier.fill %||% data$fill[1],
                              shape = outlier.shape %||% data$shape[1],
                              size = outlier.size %||% data$size[1],
                              stroke = outlier.stroke %||% data$stroke[1],
                              fill = NA,
                              alpha = outlier.alpha %||% data$alpha[1],
                              stringsAsFactors = FALSE
                            )
                            outliers_grob <- GeomPoint$draw_panel(outliers, panel_params, coord)
                          } else {
                            outliers_grob <- NULL
                          }
                          
                          ggname("geom_boxplot2", grobTree(
                            outliers_grob,
                            GeomSegment$draw_panel(whiskers, panel_params, coord),
                            GeomCrossbar$draw_panel(box, fatten = fatten, panel_params, coord)
                          ))
                        },
                        
                        draw_key = draw_key_boxplot,
                        
                        default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                                          alpha = NA, shape = 19, linetype = "solid"),
                        
                        required_aes = c("x", "lower", "upper", "middle", "ymin", "ymax")
)

#################################################################
stat_boxplot <- function(mapping = NULL, data = NULL,
                         geom = "boxplot", position = "dodge",
                         ...,
                         coef = 1.5,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatBoxplot,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      coef = coef,
      ...
    )
  )
}


#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatBoxplot <- ggproto("StatBoxplot", Stat,
                       required_aes = c("x", "y"),
                       non_missing_aes = "weight",
                       
                       setup_params = function(data, params) {
                         params$width <- params$width %||% (resolution(data$x) * 0.75)
                         
                         if (is.double(data$x) && !has_groups(data) && any(data$x != data$x[1L])) {
                           warning(
                             "Continuous x aesthetic -- did you forget aes(group=...)?",
                             call. = FALSE)
                         }
                         
                         params
                       },
                       
                       compute_group = function(data, scales, width = NULL, na.rm = FALSE, coef = 1.5) {
                         qs <- c(0, 0.25, 0.5, 0.75, 1)
                         
                         if (!is.null(data$weight)) {
                           mod <- quantreg::rq(y ~ 1, weights = weight, data = data, tau = qs)
                           stats <- as.numeric(stats::coef(mod))
                         } else {
                           stats <- as.numeric(stats::quantile(data$y, qs))
                         }
                         names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
                         iqr <- diff(stats[c(2, 4)])
                         
                         outliers <- data$y < (stats[2] - coef * iqr) | data$y > (stats[4] + coef * iqr)
                         if (any(outliers)) {
                           stats[c(1, 5)] <- range(c(stats[2:4], data$y[!outliers]), na.rm = TRUE)
                         }
                         
                         if (length(unique(data$x)) > 1)
                           width <- diff(range(data$x)) * 0.9
                         
                         df <- as.data.frame(as.list(stats))
                         df$outliers <- list(data$y[outliers])
                         
                         if (is.null(data$weight)) {
                           n <- sum(!is.na(data$y))
                         } else {
                           # Sum up weights for non-NA positions of y and weight
                           n <- sum(data$weight[!is.na(data$y) & !is.na(data$weight)])
                         }
                         
                         df$notchupper <- df$middle + 1.58 * iqr / sqrt(n)
                         df$notchlower <- df$middle - 1.58 * iqr / sqrt(n)
                         
                         df$x <- if (is.factor(data$x)) data$x[1] else mean(range(data$x))
                         df$width <- width
                         df$relvarwidth <- sqrt(n)
                         df
                       }
)



