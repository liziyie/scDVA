# Load packages -----------------
getPackage <- function(pkg, check = TRUE, load = TRUE, silent = FALSE, repos = "http://cran.us.r-project.org", github = NULL) {
  if(check){
    if(!suppressMessages(suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE)))){
      if (is.null(github)){
        try(install.packages(pkg, repos = repos), silent = TRUE)
      }
      else{
        try(remotes::install_github(github))
      }
    }
  }
  if(load) suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE))
  if(load & !silent) message("Loaded ", pkg)
}
packages <- c("ggplot2", "png", "ggedit", "cowplot", "ggbeeswarm", "ggrepel", "corrplot", "grid", "RColorBrewer", "dplyr", "egg", "gtable", "scales", "reshape2")
lapply(packages, getPackage)

# Load functions----
source("./dataprepare_utils.R")
source("./data/color_panel.R")

# Embedding plot -----------------
DoEmbeddingPlot <- function(plot.data, genes, 
                            color.by = "Genes", color.panel,
                            embedding.tag, embedding.used, 
                            seperate.plot = TRUE, nrow = 2,
                            dot.size = 2, font.size = 15, 
                            AI.friendly = FALSE) {
  ## Plot the tSNE scatter plot, cells could be colored by the gene expression,
  ## clusters, tissue, patient, treatment or day.
  ##
  ## Args:
  #'  @plot.data: The plot data created from getPlotData function.
  #'  @genes: The genes used in the plot.
  #'  @color.by: The cells should be colored by which, one of "Genes", 
  #'  "Glbal_Cluster", "Sub_Cluster", "Sample", "Tissue", "Treatment", "Day", 
  #'  "nGene" and "nUMI.
  #'  @color.panel: The color palette used when color by genes expression.
  #'  @embedding.tag: The embeddings' tag used, should be "Global" or "Sub".
  #'  @embedding.used: The embeddings used, should be "tSNE" or "UMAP".
  #'  @seperate.plot: Whether to plot all genes seperately.
  #'  @nrow: The row of final plot.
  #'  @dot.size: The size of point, default by 2.
  #'  @font.szie: The size of text, default by 15.
  #'  @AI.friendly: Whether to save an AI.friendly plot.
  ##
  ## Returns:
  ##  A ggplot item.
  x.axis <- sprintf("%s_%s_1", embedding.tag, embedding.used)
  y.axis <- sprintf("%s_%s_2", embedding.tag, embedding.used)
  if(color.by == "Genes"){  # Color by genes expression
    if(color.panel %in% c("RdYlBu", "RdYlGn", "Spectral", "RdBu")){
      myColorPalette <- colorRampPalette(rev(brewer.pal(8, color.panel)))
    }else{
      myColorPalette <- colorRampPalette(brewer.pal(8, color.panel))
    }
    if(seperate.plot){
      plot.data <- melt(plot.data[,c(genes, x.axis, y.axis)], id.vars = c(x.axis, y.axis))
      ggplot(plot.data %>% arrange(value), aes_string(x = x.axis, y = y.axis)) +
        geom_point(aes(color = value), alpha = 0.8, size = dot.size) +
        facet_wrap(~variable, nrow = nrow) +
        labs(x = x.axis, y = y.axis) +
        scale_colour_gradientn("Exp", colors = myColorPalette(100)) +
        theme_cowplot(font_size = font.size) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              strip.background = element_blank(), strip.text = element_text(size = font.size),
              legend.position = "bottom", legend.justification = "center",
              panel.border = element_rect(colour = "black")) +
        guides(color = guide_colourbar(barwidth = 20, barheight = 1)) -> p
      if(AI.friendly){
        p <- ggAIplot.grid(p, "variable")
      }
    }else{
      ggplot(plot.data %>% arrange(Expression), aes_string(x = x.axis, y = y.axis)) +
        geom_point(aes(color = Expression), alpha = 0.8, size = dot.size) +
        labs(x = x.axis, y = y.axis) +
        scale_colour_gradientn("Exp", colors = myColorPalette(100)) +
        theme_cowplot(font_size = font.size) +
        theme(legend.position = "bottom",
              legend.justification = "center") +
        guides(color = guide_colourbar(barwidth = 20, barheight = 1)) -> p
      if(AI.friendly){
        p <- ggAIplot(p)
      }
    }
  }else if(color.by %in% c("nGene","nUMI")){  # Color by nGene or nUMI
    plot.data[,paste0("log2(",color.by,")")] <- log2(plot.data$color.by + 1)
    myColorPalette <- colorRampPalette(rev(brewer.pal(8, "Spectral")))
    ggplot(plot.data, aes_string(x = x.axis, y = y.axis)) +
      geom_point(aes_string(color = paste0("log2(",color.by,")")), alpha = 0.8, size = dot.size) +
      labs(x = x.axis, y = y.axis) +
      scale_colour_gradientn(paste0("log2(",color.by,")"), colors = myColorPalette(100)) +
      theme_cowplot(font_size = font.size) +
      theme(legend.position = "bottom",
            legend.justification = "center") +
      guides(color = guide_colourbar(barwidth = 20, barheight = 1)) -> p
    if(AI.friendly){
      p <- ggAIplot(p)
    }
  }else{  # Color by cluster or other metadata
    if(exists(paste0(color.by, "_color_panel"))){
      color_panel <- get(paste0(color.by, "_color_panel")) 
      if(sum(!unique(plot.data[, color.by]) %in% names(color_panel)) > 0){
        color_panel <- c68
      }
    }else{
      color_panel <- c68
    }
    ggplot(plot.data, aes_string(x = x.axis, y = y.axis)) +
      geom_point(aes_string(color = color.by), size = dot.size, alpha = 0.8) +
      theme_cowplot(font_size = font.size) +
      labs(x = x.axis, y = y.axis) +
      scale_color_manual(values = color_panel) +
      theme(legend.position = "top", legend.justification = "right") +
      guides(color = guide_legend(title = NULL, nrow = 4, byrow = TRUE,
             override.aes = list(size = 4))) -> p
    if(AI.friendly){
      p <- ggAIplot(p)
    }
    if (color.by %in% c("Global_Cluster", "Sub_Cluster")){
      label_pos <-
        aggregate(plot.data[, c(x.axis, y.axis)], list(plot.data[, color.by]),
                  median) %>% data.frame()
      p <- p + geom_text_repel(data = label_pos, 
                               aes_string(x = x.axis, y = y.axis, label = "Group.1"),
                               size = font.size / 3)
    }
  }
  grid.newpage()
  grid.draw(p)
}
  
# Distribution plot -----------------
DoDistributionPlot <- function(plot.data, genes, plot.type,
                               group.by = "Sub_Cluster", color.by, color.panel,
                               python.type = FALSE, seperate.plot = TRUE, 
                               nrow = 2, font.size = 15, dot.size = 0.6, AI.friendly = FALSE) {
  ## Plot the gene distribution plot, cells could be grouped by clusters,
  ## tissue, patient, treatment or day and each group was colored by the mean expression
  ## in the group.
  ##
  ## Args:
  #'  @plot.data: The plot data created from getPlotData function.
  #'  @genes: The genes used in the plot.
  #'  @plot.type: Should be "Box Plot" or "Violin Plot" .
  #'  @group.by: The distribuion plot should be grouped by which column,
  #'  should be "Patient", "Tissue", "Treatment", "Day", "Sub_Cluster" or "Global_Cluster".
  #'  @color.by: The cells should be colored by which, one of "Exp", 
  #'  "Glbal_Cluster", "Sub_Cluster", "Sample", "Tissue", "Treatment" and "Day".
  #'  @color.panel: The color palette used when color by genes expression.
  #'  @python.type: Plot a python-theme violin, only works when group.by is same as color.by.
  #'  @seperate.plot: Whether to plot all genes seperately.
  #'  @dot.size: The size of point, default by 1.
  #'  @font.szie: The size of text, default by 15.
  #'  @AI.friendly: Whether to save an AI.friendly plot.
  ##
  ## Returns:
  ##  A ggplot item with the distribution plot in box or violin type.
  if(color.by == "Exp"){
    plot.data[,"Exp"] <- "A"
    if (color.panel %in% c("RdYlBu", "RdYlGn", "Spectral", "RdBu")) {
      myColorPalette <- colorRampPalette(rev(brewer.pal(8, color.panel)))
    }else{
      myColorPalette <- colorRampPalette(brewer.pal(8, color.panel))
    }
  }else{
    if(exists(paste0(color.by, "_color_panel"))){
      color_panel <- get(paste0(color.by, "_color_panel")) 
      if(sum(!unique(plot.data[, color.by]) %in% names(color_panel)) > 0){
        color_panel <- c68
      }
    }else{
      color_panel <- c68
    }
  }
  plot.data <- plot.data[,unique(c("Expression", genes, group.by, color.by))]
  plot.data <- melt(plot.data, id.vars = unique(c(group.by, color.by)))
  if(group.by != color.by){
    gene_mean <- aggregate(plot.data$value, list(plot.data[,group.by], plot.data[,color.by], plot.data$variable), mean)
    colnames(gene_mean) <- c(group.by, color.by, "variable", "Mean")
    gene_median <- aggregate(plot.data$value, list(plot.data[,group.by], plot.data[,color.by], plot.data$variable), median)
    colnames(gene_median) <- c(group.by, color.by, "variable", "Median")
  }else{
    gene_mean <- aggregate(plot.data$value, list(plot.data[,group.by], plot.data$variable), mean)
    colnames(gene_mean) <- c(group.by, "variable", "Mean")
    gene_median <- aggregate(plot.data$value, list(plot.data[,group.by], plot.data$variable), median)
    colnames(gene_median) <- c(group.by, "variable", "Median") 
  }
  gene_median$Label <- round(gene_median$Median, 2)
  plot.data <- merge(plot.data, gene_mean)
  plot.data$inter <- interaction(plot.data[,group.by], plot.data[,color.by])
  if(seperate.plot){
    ggplot(plot.data %>% filter(variable != "Expression"), aes_string(x = group.by, y = "value")) +
      labs(x = "", y = "Exp") +
      theme_bw() + 
      theme(
        axis.text.x = element_text(size = font.size * 0.8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = font.size),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = font.size),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black")
      ) -> p
    if(plot.type == "Box Plot"){
      if(color.by == "Exp"){
        p <- p +
          geom_boxplot(outlier.size = -1, aes(color = Mean)) +
          geom_quasirandom(cex = dot.size, width = 0.25, alpha = 0.5, aes(col = Mean), groupOnX = TRUE) +
          geom_text(data = gene_median %>% filter(variable != "Expression"), aes_string(label = "Label", x = group.by, y = "Median"), size = 5, vjust = -0.5) +
          scale_color_gradientn("Exp", colors = myColorPalette(100)) +
          facet_wrap(~variable, nrow = nrow, scales = "free") +
          guides(color = guide_colorbar(barwidth = 1, barheight = 12))
      }else if(color.by == group.by){
        p <- p +
          geom_boxplot(outlier.size = -1, aes_string(color = color.by)) +
          geom_quasirandom(cex = dot.size, width = 0.25, alpha = 0.5, aes_string(col = color.by), groupOnX = TRUE) +
          geom_text(data = gene_median %>% filter(variable != "Expression"), aes_string(label = "Label", x = group.by, y = "Median"), size = 5, vjust = -0.5) +
          scale_color_manual(values = color_panel) +
          facet_wrap(~variable, nrow = nrow, scales = "free")
      }else{
        p <- p +
          geom_boxplot(outlier.size = -1, aes_string(color = color.by)) +
          geom_quasirandom(cex = dot.size, alpha = 0.5, aes_string(col = color.by), groupOnX = TRUE, varwidth = TRUE, dodge.width = 0.7) +
          geom_text(data = gene_median %>% filter(variable != "Expression"), aes_string(label = "Label", x = group.by, y = "Median", group = color.by), size = 5, vjust = -0.5, position = position_dodge(.8)) +
          scale_color_manual(values = color_panel) +
          facet_wrap(~variable, nrow = nrow, scales = "free")
      }
    }else if(plot.type == "Violin Plot"){
      if(color.by == "Exp"){
        p <- p +
          geom_violin(scale = "width", aes(fill = Mean), color = "white", kernel = "gaussian") +
          geom_boxplot(outlier.size = -1, width = .1, fill = "white") +
          geom_text(data = gene_median %>% filter(variable != "Expression"), aes_string(label = "Label", x = group.by, y = "Median"), size = 5, vjust = -0.5) +
          scale_fill_gradientn("Exp", colors = myColorPalette(100)) +
          facet_wrap(~variable, nrow = nrow, scales = "free") +
          guides(fill = guide_colorbar(barwidth = 1, barheight = 12))
      }else if(color.by == group.by){
        p <- p +
          geom_violin(scale = "width", aes_string(fill = color.by), color = "white", kernel = "gaussian") +
          geom_boxplot(outlier.size = -1, width = .1, fill = "white") +
          geom_text(data = gene_median %>% filter(variable != "Expression"), aes_string(label = "Label", x = group.by, y = "Median"), size = 5, vjust = -0.5) +
          scale_fill_manual(values = color_panel) +
          facet_wrap(~variable, nrow = nrow, scales = "free")
        if(python.type){
          p <-
            ggplot(plot.data %>% filter(variable != "Expression"), aes(x = variable, y = value)) +
            geom_violin(scale = "width", aes_string(fill = color.by), color = "white", kernel = "gaussian") +
            facet_grid(as.formula(paste(group.by, "~.")), switch = "y", scales = "free_y") +
            scale_fill_manual(values = color_panel) +
            theme_cowplot(font_size = 15) +
            scale_y_continuous(position = "right") +
            theme(legend.position = "null",
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  axis.title = element_blank(),
                  panel.spacing = unit(0.05, "lines"),
                  panel.border = element_rect(colour = "black", fill=NA, size = .75),
                  strip.background =element_rect(fill="white"),
                  strip.text.y = element_text(angle = 180, vjust = 0.5, hjust = 1))
        }
      }else{
        p <- p +
          geom_violin(scale = "width", aes_string(fill = color.by), color = "white", kernel = "gaussian", position = position_dodge(width = 1)) +
          geom_boxplot(aes(group = inter), width = 0.1, fill = "white", position = position_dodge(width = 1), outlier.shape = NA) +
          geom_text(data = gene_median %>% filter(variable != "Expression"), aes_string(label = "Label", x = group.by, y = "Median", group = color.by), size = 5, vjust = -0.5, position = position_dodge(.8)) +
          scale_fill_manual(values = color_panel) +
          facet_wrap(~variable, nrow = nrow, scales = "free")
      }
    }
  }else{
    ggplot(plot.data %>% filter(variable == "Expression"), aes_string(x = group.by, y = "value")) +
      labs(x = "", y = "Exp") +
      theme_bw() + 
      theme(
        axis.text.x = element_text(size = font.size * 0.8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = font.size),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = font.size),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black")
      ) -> p
    if(plot.type == "Box Plot"){
      if(color.by == "Exp"){
        p <- p +
          geom_boxplot(outlier.size = -1, aes(color = Mean)) +
          geom_quasirandom(cex = dot.size, width = 0.25, alpha = 0.5, aes(col = Mean), groupOnX = TRUE) +
          geom_text(data = gene_median %>% filter(variable == "Expression"), aes_string(label = "Label", x = group.by, y = "Median"), size = 5, vjust = -0.5) +
          scale_color_gradientn("Exp", colors = myColorPalette(100)) +
          guides(color = guide_colorbar(barwidth = 1, barheight = 12))
      }else if(color.by == group.by){
        p <- p +
          geom_boxplot(outlier.size = -1, aes_string(color = color.by)) +
          geom_quasirandom(cex = dot.size, width = 0.25, alpha = 0.5, aes_string(col = color.by), groupOnX = TRUE) +
          geom_text(data = gene_median %>% filter(variable == "Expression"), aes_string(label = "Label", x = group.by, y = "Median"), size = 5, vjust = -0.5) +
          scale_color_manual(values = color_panel)
      }else{
        p <- p +
          geom_boxplot(outlier.size = -1, aes_string(color = color.by)) +
          geom_quasirandom(cex = dot.size, alpha = 0.5, aes_string(col = color.by), groupOnX = TRUE, varwidth = TRUE, dodge.width = 0.7) +
          geom_text(data = gene_median %>% filter(variable == "Expression"), aes_string(label = "Label", x = group.by, y = "Median", group = color.by), size = 5, vjust = -0.5, position = position_dodge(.8)) +
          scale_color_manual(values = color_panel)
      }
    }else if(plot.type == "Violin Plot"){
      if(color.by == "Exp"){
        p <- p +
          geom_violin(scale = "width", aes(fill = Mean), color = "white", kernel = "gaussian") +
          geom_boxplot(outlier.size = -1, width = .1, fill = "white") +
          geom_text(data = gene_median %>% filter(variable == "Expression"), aes_string(label = "Label", x = group.by, y = "Median"), size = 5, vjust = -0.5) +
          scale_fill_gradientn("Exp", colors = myColorPalette(100)) +
          guides(fill = guide_colorbar(barwidth = 1, barheight = 12))
      }else if(color.by == group.by){
        p <- p +
          geom_violin(scale = "width", aes_string(fill = color.by), color = "white", kernel = "gaussian") +
          geom_boxplot(outlier.size = -1, width = .1, fill = "white") +
          geom_text(data = gene_median %>% filter(variable == "Expression"), aes_string(label = "Label", x = group.by, y = "Median"), size = 5, vjust = -0.5) +
          scale_fill_manual(values = color_panel)
      }else{
        p <- p +
          geom_violin(scale = "width", aes_string(fill = color.by), color = "white", kernel = "gaussian", position = position_dodge(width = 1)) +
          geom_boxplot(aes(group = inter), width = 0.1, fill = "white", position = position_dodge(width = 1), outlier.shape = NA) +
          geom_text(data = gene_median %>% filter(variable == "Expression"), aes_string(label = "Label", x = group.by, y = "Median", group = color.by), size = 5, vjust = -0.5, position = position_dodge(.8)) +
          scale_fill_manual(values = color_panel)
      }
    }
  }
  return(p)
}

# Significance plot -----------------
DoSignificancePlot <- function(plot.data, group.by, per.cutoff,
                               font.size, sig.level = 0) {
  ## Plot the gene differentially expressed pattern using a correlation plot.
  ##
  ## Args:
  #'  @plot.data: The plot data created from getPlotData function.
  #'  @group.by: The distribuion plot should be grouped by which column,
  #'  should be "Patient", "Tissue", "Treatment", "Day", "Sub_Cluster" or "Global_Cluster".
  #'  @per.cutoff: The cutoff of a gene seen as expressed in the cells.
  #'  @font.size: The character size showing the log2 fold change.
  #'  @sig.level: The p-value cutoff as significant.
  ##
  ## Returns:
  ##  A correlation style plot. 
  CorrList <- getSigData(plot.data, group.by, per.cutoff)
  if(!is.null(CorrList)){
    corrplot(
      corr = CorrList$percentage, p.mat = CorrList$pvalue,
      method = "pie", outline = "white", col = "lightblue", addgrid.col = NA,
      type = "upper", tl.pos = "n", cl.pos = "n", insig = "blank"
    )
    corrplot(
      corr = CorrList$coeff, p.mat = CorrList$pvalue,
      add = TRUE, method = "color", is.corr = F, type = "lower",
      addCoef.col = "black", number.cex = font.size, 
      tl.cex = font.size, tl.pos = "ld", tl.srt = 45, tl.offset = 1, cl.pos = "n",
      diag = F, insig = "blank", sig.level = sig.level,
      col = rev(brewer.pal(n = 8, name = "RdYlBu"))
    )
  }
}

# Heatmap plot -----------------
DoHeatmapPlot <- function(plot.data, genes, group.by, color.panel, font.size = 15){
  ## Do the geneset heatmap plot.
  ##
  ## Args:
  #'  @plot.data: The plot data created from getPlotData function.
  #'  @genes: The genes used in the plot.
  #'  @group.by: The cells should be ordered by or grouped by.
  #'  @color.panel: The color palette used when color by genes expression.
  #'  @font.size: The size of text, default by 15.
  ##
  ## Returns:
  ## A heatmap plot.
  if(length(genes) > 1){
    expression.matrix <- t(plot.data[, genes])
    genes <- genes[rowSums(expression.matrix) != 0]
    expression.matrix <- expression.matrix[rowSums(expression.matrix) != 0, ]
    groupid <- factor(plot.data[, group.by])
    
    if(length(unique(groupid)) > 1){
      # Boxplot in the left panel
      signature.mean.expression <- colMeans(expression.matrix[genes,])
      boxplot.data <- data.frame(Expression = signature.mean.expression, Group = groupid)
      p_boxplot <- 
        ggplot(boxplot.data, aes(x = Group, y = Expression)) + 
        geom_boxplot(color = "#1979B5", fill = "#1979B5", 
                     outlier.colour = "black", outlier.shape = 1) + 
        theme_bw() +
        theme(legend.position = "NULL",
              axis.text.y = element_text(size = font.size * 0.8),
              axis.text.x = element_text(size = font.size),
              axis.title = element_blank()) +
        scale_x_discrete(limits = rev(levels(as.factor(groupid)))) +
        coord_flip()
      dat <- ggplot_build(p_boxplot)$data[[1]]
      dat$xsp <- 1/2 * (dat$xmax + dat$xmin) - 1/4 * (dat$xmax - dat$xmin)
      dat$xep <- 1/2 * (dat$xmax + dat$xmin) + 1/4 * (dat$xmax - dat$xmin)
      p_boxplot <- 
        p_boxplot + 
        geom_segment(data = dat,  aes(x = xmin, xend = xmax, y = middle, yend = middle), colour = "#2BA147", size = 2) +
        geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymin, yend = ymin), colour = "black", size = 2) +
        geom_segment(data = dat,  aes(x = xsp, xend = xep, y = ymax, yend = ymax), colour = "black", size = 2)
      # Heatmap in the right panel
      if (color.panel %in% c("RdYlBu", "RdYlGn", "Spectral", "RdBu")) {
        myColorPalette <- colorRampPalette(rev(brewer.pal(8, color.panel)))
      }else{
        myColorPalette <- colorRampPalette(brewer.pal(8, color.panel))
      }
      gene.median.expression <- aggregate(t(expression.matrix[genes,]), list(Group = groupid), mean)
      gene.median.expression[, 2:(length(genes) + 1)] <- apply(gene.median.expression[, 2:(length(genes) + 1)], 2, function(x){
        (x - mean(x))/sd(x)
      })
      heatmap.data <- reshape2::melt(gene.median.expression, id = "Group")
      heatmap.data$Group <- as.factor(heatmap.data$Group)
      heatmap.data$Group <- factor(heatmap.data$Group, levels = rev(levels(heatmap.data$Group)))
      max_value <- round(max(abs(heatmap.data$value)) + 0.05, 1)
      scale_breaks <- round(max_value * c(-4/5, -2/5, 0, 2/5, 4/5), 1)
      p_heatmap <- ggplot(heatmap.data, aes(x = factor(variable), y = Group)) +
        geom_tile(aes(fill = value), colour = "white") + 
        theme_bw() + 
        theme(panel.border = element_blank(),
              legend.title = element_blank(),
              axis.line = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                         hjust = 1, size = font.size * 0.8),
              axis.title = element_blank(),
              axis.ticks = element_blank()) +
        scale_fill_gradientn(limits = c(-max_value, max_value), breaks = scale_breaks, expand = c(0,0), colours = myColorPalette(100), guide = guide_colourbar(label.theme = element_text(size = font.size * 0.5))) +
        scale_x_discrete(breaks = unique(heatmap.data$variable), labels = unique(heatmap.data$variable)) + 
        scale_y_discrete(expand = c(0,0))
      gt = ggplotGrob(p_heatmap)
      leg = gtable_filter(gt, "guide-box")
      leg[[1]][[1]][[1]][[1]][[1]][[2]]$height = unit(1, "npc")
      pos = unit.c(unit(0.01,"npc"), unit(.25, "npc"), unit(.5, "npc"), unit(.75, "npc"), unit(.99, "npc"))
      leg[[1]][[1]][[1]][[1]][[1]][[5]]$y0 = pos
      leg[[1]][[1]][[1]][[1]][[1]][[5]]$y1 = pos
      leg[[1]][[1]][[1]][[1]]$heights = unit.c(rep(unit(0, "mm"), 3),
                                               unit(1, "npc"),
                                               unit(0, "mm"))
      leg[[1]][[1]]$heights[[3]] = sum(rep(unit(0, "mm"), 3),
                                       unit(1, "npc"),
                                       unit(0, "mm"))
      if(strsplit(as.character(packageVersion("ggplot2")),"[.]")[[1]][1] == "2"){
        # ggplot 2.0
        leg[[1]][[1]][[1]][[1]][[1]][[3]]$y = pos
        p_heatmap = gtable_add_grob(gt, leg, t = 6, l = 8)
      }else if(strsplit(as.character(packageVersion("ggplot2")),"[.]")[[1]][1] == "3"){
        # ggplot 3.0
        leg[[1]][[1]][[1]][[1]][[1]][[3]]$children[[1]]$y = pos
        p_heatmap = gtable_add_grob(gt, leg, t = 7, l = 9)
      }
      # Final plot
      g1 <- ggplotGrob(p_boxplot)
      g2 <- p_heatmap
      fg1 <- gtable_frame(g1, height = unit(10, "null"))
      fg2 <- gtable_frame(g2, height = unit(40, "null"))
      p <- gtable_frame(gtable_cbind(fg1, fg2))
      grid.newpage()
      grid.draw(p)
    } else{
      ggplot(data.frame()) + 
        ggtitle("There should be more than two groups") + 
        theme_bw() + 
        theme(plot.title = element_text(
          size = 26,
          face = "bold",
          hjust = 0.5)
        )
    }
  }else{
    ggplot(data.frame()) +
      ggtitle("There should be more than two genes") +
      theme_bw() +
      theme(plot.title = element_text(
        size = 26,
        face = "bold",
        hjust = 0.5
      ))
  }
}

# In-silico FACS plot ----------------- 
DoISFACSPlot <- function(plot.data, genes, color.by, color.panel, 
                         dot.size = -1, x.cutoff = 1, y.cutoff = 1, 
                         font.size = 15, AI.friendly = FALSE){
  ## Plot a in-silico FACS plot with two genes, users could set the cutoff of
  ## x and y axis (gate) manually.
  ##
  ## Args:
  #'  @plot.data: The plot data created from getPlotData function.
  #'  @genes: The genes used in the plot, must be a vector with only two items.
  #'  @color.by: The point was colored by "None", "Glbal_Cluster", "Sub_Cluster", 
  #'  "Sample", "Tissue", "Treatment" and "Day".
  #'  @color.panel: The color palette used when color by genes expression.
  #'  @dot.size: The size of point, default by -1, not showing the point.
  #'  @x.cutoff: The cutoff of x axis.
  #'  @y.cutoff: The cutoff of y axis.
  #'  @font.szie: The size of text, default by 15.
  #'  @AI.friendly: Whether to save an AI.friendly plot.
  ##
  ## Returns:
  ##  A ggplot item with the in-silico FACS plot.
  if(color.panel %in% c("RdYlBu", "RdYlGn", "Spectral", "RdBu")) {
    myColorPalette <- colorRampPalette(rev(brewer.pal(8, color.panel)))
  }else{
    myColorPalette <- colorRampPalette(brewer.pal(8, color.panel))
  }
  if(color.by != "None"){
    if(exists(paste0(color.by, "_color_panel"))){
      color_panel <- get(paste0(color.by, "_color_panel")) 
      if(sum(!unique(plot.data[, color.by]) %in% names(color_panel)) > 0){
        color_panel <- c68
      }
    }else{
      color_panel <- c68
    }
  }
  if(length(genes) == 2){
    gene.A <- genes[1]
    gene.B <- genes[2]
    G4.num <- c(
      sum(plot.data[,gene.A] <= x.cutoff & plot.data[,gene.B] > y.cutoff),
      sum(plot.data[,gene.A] <= x.cutoff & plot.data[,gene.B] <= y.cutoff),
      sum(plot.data[,gene.A] > x.cutoff & plot.data[,gene.B] > y.cutoff),
      sum(plot.data[,gene.A] > x.cutoff & plot.data[,gene.B] <= y.cutoff)
    )
    G4.per <- round(G4.num / nrow(plot.data) * 100,2)
    G4.text <- data.frame(
      X = rep(c(x.cutoff / 2, (max(plot.data[,gene.A]) + x.cutoff + 0.5) / 2), each = 2),
      Y = rep(c(max(plot.data[,gene.B]) - 0.5, y.cutoff - 0.5),2),
      group = paste0("Grp", 1:4, " : ", gene.A, c("-","-","+","+"), gene.B, c("+","-","+","-")),
      text = paste0(G4.per,"% (",G4.num,"|",nrow(plot.data),")"))
    if(color.by == "None"){
      ggplot(plot.data, aes_string(x = gene.A, y = gene.B)) +
        geom_point(color = "gray", alpha = 0.8, size = dot.size) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.A]) + 0.5)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.B]) + 0.5)) +
        theme_bw() +
        theme(
          axis.text = element_text(size = font.size * 0.8),
          axis.title = element_text(size = font.size),
          legend.text = element_text(size = font.size * 0.8),
          legend.title = element_text(size = font.size * 0.8)
        ) -> pmain
      if(AI.friendly){
        pmain <- ggAIplot(pmain)
      }
      pmain <- pmain +
        stat_density_2d(aes(fill = stat(level), alpha = stat(level)), geom = "polygon") +
        scale_alpha(range = c(0.6,1)) + 
        scale_fill_gradientn(colours = myColorPalette(100)) +
        geom_hline(yintercept = y.cutoff, linetype = "dashed") +  
        geom_vline(xintercept = x.cutoff, linetype = "dashed") +
        geom_text(data = G4.text, aes(x = X, y = Y, label = paste0(group ,"\n",text)), size = font.size * 0.3) +
        guides(alpha = FALSE, fill = FALSE)
      xdens <- axis_canvas(pmain, axis = "x") +
        geom_density(data = plot.data, aes_string(x = gene.A), fill = "gray",
                     color = NA, alpha = 0.7, size = 0.2) 
      p0 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
      ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
        geom_density(data = plot.data, aes_string(x = gene.B), fill = "gray",
                     color = NA, alpha = 0.7, size = 0.2) +
        coord_flip()
      p <- insert_yaxis_grob(p0, ydens, grid::unit(.2, "null"), position = "right")
    }
    else{
      ggplot(plot.data, aes_string(x = gene.A, y = gene.B)) +
        geom_point(aes_string(color = color.by, group = color.by), alpha  = 0.8, size = dot.size) +
        scale_color_manual("", values = color_panel) +
        guides(color = guide_legend(override.aes = list(fill = NA, nrow = 4, byrow = TRUE, size = 4, linetype = 0))) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.A]) + 0.5)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.B]) + 0.5)) +
        theme_bw() +
        theme(
          axis.text = element_text(size = font.size * 0.8),
          axis.title = element_text(size = font.size),
          legend.text = element_text(size = font.size * 0.8),
          legend.title = element_text(size = font.size * 0.8)
        ) -> pmain
      if(AI.friendly){
        pmain <- ggAIplot(pmain)
      }
      pmain <- pmain +
        stat_density_2d(aes(fill = stat(level), alpha = stat(level)), geom = "polygon") +
        scale_alpha(range = c(0.6,1)) + 
        scale_fill_gradientn(colours = myColorPalette(100)) +
        geom_hline(yintercept = y.cutoff, linetype = "dashed") +  
        geom_vline(xintercept = x.cutoff, linetype = "dashed") +
        geom_text(data = G4.text, aes(x = X, y = Y, label = paste0(group ,"\n",text)), size = font.size * 0.3) +
        guides(alpha = FALSE, fill = FALSE)
      xdens <- axis_canvas(pmain, axis = "x") +
        geom_density(data = plot.data, 
                     aes_string(x = gene.A, fill = color.by),
                     color = NA, alpha = 0.7, size = 0.2) +
        scale_fill_manual(values = color_panel)
      p0 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
      ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
        geom_density(data = plot.data, 
                     aes_string(x = gene.B, fill = color.by),
                     color = NA, alpha = 0.7, size = 0.2) +
        coord_flip() +
        scale_fill_manual(values = color_panel)
      p <- insert_yaxis_grob(p0, ydens, grid::unit(.2, "null"), position = "right")
    }
    grid.newpage()
    grid.draw(p)
  }else{
    ggplot(data.frame()) +
      ggtitle("There should be two genes as input!") +
      theme_bw() +
      theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
  }
} 

DoISFACSVolcanoPlot <- function(limma.data, dot.size = 1, adj.P.Val.cutoff = 0.05,
                                logFC.cutoff = 0.5, ngenes.labeled = 25,
                                font.size = 15, AI.friendly = FALSE){
  ## Plot a in-silico FACS plot with two genes, users could set the cutoff of
  ## x and y axis (gate) manually.
  ##
  ## Args:
  #'  @limma.data: The plot data created from getLIMMAData function.
  #'  @dot.size: The size of point, default by 1.
  #'  @adj.P.Val.cutoff: The cutoff of adjusted p-value to define significant genes.
  #'  @logFC.cutoff: The cutoff of log fold change to define significant genes.
  #'  @ngenes.labeled: The number of genes with gene symbol labeled, ordered by the 
  #'  log fold change. 
  #'  @font.szie: The size of text, default by 15.
  #'  @AI.friendly: Whether to save an AI.friendly plot.
  ##
  ## Returns:
  ##  A ggplot item with the volcano plot.
  if(is.null(limma.data$Warning)){
    limma.data$Sig <- abs(limma.data$logFC) >= logFC.cutoff & limma.data$adj.P.Val < adj.P.Val.cutoff
    limma.data.genes_labeled <- c(order(limma.data$logFC)[1:ngenes.labeled], order(-limma.data$logFC)[1:ngenes.labeled])
    limma.data[limma.data.genes_labeled, "label"] <- limma.data[limma.data.genes_labeled, "Symbol"]
    limma.data$adj.P.Val[limma.data$adj.P.Val == 0] <- min(limma.data$adj.P.Val[limma.data$adj.P.Val > 0])
    p <- ggplot(limma.data, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(aes(color = Sig), size = dot.size) +
      scale_color_manual(values = Binomial_color_panel) +
      labs(x = "log2 fold change", y = "-log10 adj p-value") +
      theme_cowplot(font_size = font.size) +
      theme(legend.position = "none")
    if(AI.friendly){
      p <- ggAIplot(p)
    }
    p <- p +
      geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), color = "grey", linetype = "dashed", lwd = 1) +
      geom_text_repel(data = limma.data, aes(label = label), size = font.size / 3)
    return(p)
  }else{
    ggplot(data.frame()) +
      ggtitle(limma.data$Warning) +
      theme_bw() +
      theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
  }
}
  
# Metadata plot -----------------
DoMetadataPlot <- function(plot.data, plot.type, group.by1, color.by, group.by2, quantify = "Proportion", coord.flip = TRUE, show.cutoff = 10, nrow = 2, font.size = 15) {
  ## Plot the cluster, patient or tissue proportion using bar plot.
  ##
  ## Args:
  #'  @plot.data: The plot data created from getPlotData function.
  #'  @plot.type: Should be "Bar Plot" or "Pie Plot".
  #'  @group.by1: The xlab of the metadata plot. One of "Sub_Cluster", 
  #'  "Global_Cluster", "Day", "Treatment", "Patient" or "Tissue".
  #'  @color.by: The color of each part in the metadata plot. One of "Sub_Cluster", 
  #'  "Global_Cluster", "Day", "Treatment", "Patient" or "Tissue".
  #'  @group.by2: The metadata plot should be facet by which. One of "None", 
  #'  "Sub_Cluster", "Global_Cluster", "Day", "Treatment", "Patient" or "Tissue".
  #'  @quantify: "Proportion" or "Count", default by "Proportion".
  #'  @coord.flip: Whether to inverse the x axis and y axis, default by TRUE.
  #'  @show.cutoff: The minimum proportions to be shown in the pie plot.
  #'  @nrow: The row of facet bar plot or pie plot, default by 2.
  #'  @font.size: The font size, default by 15.
  ##
  ## Returns:
  ##  A ggplot item with the metadata plot in bar plot or pie plot. 
  if(exists(paste0(color.by, "_color_panel"))){
    color_panel <- get(paste0(color.by, "_color_panel")) 
    if(sum(!unique(plot.data[, color.by]) %in% names(color_panel)) > 0){
      color_panel <- c68
    }
  }else{
    color_panel <- c68
  } 
  if(plot.type == "Bar Plot"){
    if(quantify == "Proportion"){
      ggplot(plot.data, aes_string(x = group.by1, fill = color.by)) + 
        geom_bar(position = "fill", alpha = .8, color = "white") +
        scale_fill_manual(values = color_panel) +
        labs(y = "Percentage", x = "") +
        guides(fill = guide_legend(ncol = 4)) +
        theme_cowplot(font_size = font.size) +
        theme(legend.position = "bottom", legend.justification = "right") -> p
    }else{
      ggplot(plot.data, aes_string(x = group.by1, fill = color.by)) + 
        geom_histogram(stat = "count", alpha = .8, color = "white") +
        scale_fill_manual(values = color_panel) +
        labs(y = "Number of cells", x = "") +
        guides(fill = guide_legend(ncol = 4)) +
        theme_cowplot(font_size = font.size) +
        theme(legend.position = "bottom", legend.justification = "right") -> p
    }
    if(group.by2 != "None"){
      p <- p + facet_wrap(as.formula(paste0("~", group.by2)), nrow = nrow)
    }
    if(coord.flip){
      p <- p + coord_flip()
    }else{
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }
  }else if(plot.type == "Pie Plot"){
      pie_df <- list()
      for(i in levels(factor(plot.data[,group.by1]))){
        pie_df[[i]] <- plot.data[plot.data[,group.by1] == i,] %>% select(color.by) %>% table() %>% data.frame()
        colnames(pie_df[[i]]) <- c("group","value")
        pie_df[[i]]$value <- round(pie_df[[i]]$value / sum(pie_df[[i]]$value) * 100, 2)
        pie_df[[i]]$facet.by <- i
        pie_df[[i]] <- cbind(pie_df[[i]], y_height = pie_df[[i]]$value/2 + c(0, cumsum(pie_df[[i]]$value)[-length(pie_df[[i]]$value)]))
        pie_df[[i]]$y_height[which(pie_df[[i]]$value <= show.cutoff)] <- NA
      }
      pie_df_all <- do.call(rbind, pie_df)
      ggplot(pie_df_all, aes(x = "", y = value, fill = group)) +
        geom_bar(alpha = .8, width = 1, stat = "identity",
                 position = position_stack(reverse = TRUE)) +
        scale_fill_manual(color.by, values = color_panel) +
        guides(fill = guide_legend(ncol = 3)) +
        coord_polar("y", start = 0) +
        geom_text(aes(y = y_height, label = value), size = font.size * 0.3, na.rm = TRUE) +
        facet_wrap(~facet.by, nrow = 2) +
        theme_void() +
        theme(
          legend.position = "bottom",
          legend.justification = "right",
          strip.text = element_text(size = font.size * 0.8),
          legend.title = element_text(size = font.size * 0.8),
          legend.text = element_text(size = font.size * 0.8)
        ) -> p
  }
  return(p)
}

# AI friendly plot----
gg0point <- function(plot){
  newLayer <- cloneLayer(plot$layers[[1]])
  newLayer$aes_params$size <- -1
  newplot <- plot %>% rgg("point", 1, newLayer)
  return(newplot)
}
 
MyTheme_transparent <- theme(
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
  legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
  axis.line = element_line(colour = "black") # adding a black line for x and y axis
)

ggAIplot <- function(plot, calibration = c(0,0,0,0), width = 10, height = 10, dpi = 300){
  plot.png <- plot + theme_nothing() + MyTheme_transparent
  ggsave(plot.png, filename = "./tmp/temp.png", width = width, height = height, dpi = dpi, bg = "transparent")
  img <- readPNG("./tmp/temp.png")
  file.remove("./tmp/temp.png")
  blank.plot <- gg0point(plot)
  range.values <- c(
    ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = blank.plot)$layout$panel_params[[1]]$y.range
  )
  plot.pdf <- blank.plot +
    annotation_raster(img, 
                       xmin = range.values[1] + calibration[1], 
                       xmax = range.values[2] + calibration[2],
                       ymin = range.values[3] + calibration[3], 
                       ymax = range.values[4] + calibration[4])
  return(plot.pdf)
}
 
annotation_raster.grid <- 
  function (raster, xmin, xmax, ymin, ymax, interpolate = FALSE, data) {
    raster <- grDevices::as.raster(raster)
    layer(data = data, mapping = NULL, stat = StatIdentity, 
          position = PositionIdentity, geom = GeomRasterAnn, inherit.aes = FALSE, 
          params = list(raster = raster, xmin = xmin, xmax = xmax, 
                        ymin = ymin, ymax = ymax, interpolate = interpolate))
  }
 
ggAIplot.grid <- function(plot, facet.by, calibration = c(0,0,0,0), width = 10, height = 10, dpi = 300){
  plot.data <- plot$data
  plot.pdf <- gg0point(plot)
  range.values <- c(
    ggplot_build(plot = plot.pdf)$layout$panel_params[[1]]$x.range,
    ggplot_build(plot = plot.pdf)$layout$panel_params[[1]]$y.range
  )
  for(variable in unique(plot.data[,facet.by])){
    temp.plot.png <- plot + theme_nothing() + MyTheme_transparent
    temp.plot.png$data <- temp.plot.png$data[temp.plot.png$data[,facet.by] == variable,]
    ggsave(temp.plot.png, filename = "./tmp/temp.png", width = width, height = height, dpi = dpi, bg = "transparent")
    img <- readPNG("./tmp/temp.png")
    file.remove("./tmp/temp.png")
    ggimg <- annotation_raster.grid(
      img,
      xmin = range.values[1] + calibration[1],
      xmax = range.values[2] + calibration[2],
      ymin = range.values[3] + calibration[3],
      ymax = range.values[4] + calibration[4],
      data = temp.plot.png$data[1, ]
    )
    plot.pdf <- plot.pdf + ggimg
  }
  return(plot.pdf)
}
