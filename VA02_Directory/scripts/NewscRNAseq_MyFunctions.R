########################## PLOTTING VOLCANO PLOT  ########################## 
# This function creates a volcano plot starting from any dataframe containing logFC 
# and corresponding p values.

plot.volcano <- function(df.stat,
                         pval.col = "FDR",
                         lfc.col,
                         gene.col = "Genes",
                         point.size = 2,
                         y.title,
                         x.title,
                         y.max,
                         x.min,
                         x.max,
                         center.x.axis = T,
                         plot.title = NULL,
                         lfc.colors = c("red", "#CCCCCC", "#0066CC"),
                         lfc.thresh = c(-1,1),
                         pval.color.thresh = 0.05,
                         alpha_transparency = 0.6,
                         show.vlines = T,
                         vlines.linetype = "dashed",
                         vlines.colors = NULL,
                         show.hline = T,
                         hline.linetype = "dashed",
                         hline.color = NULL,
                         lines.thickness = 0.5,
                         label.genes = NULL,
                         label.size = 4,
                         show.gene.nbr = F)
{
  # Check user inputs.
  stopifnot("`df.stat` should be a dataframe." = class(df.stat) == "data.frame")
  stopifnot("`pval.col` should be numerical." = is.numeric(df.stat[,pval.col]))
  stopifnot("`lfc.col` should be a numerical." = is.numeric(df.stat[,lfc.col]))
  stopifnot("`gene.col` should be missing or contains only characters" = is.character(df.stat[,gene.col]) | missing(gene.col))
  stopifnot("`center.x.axis` should be a boolean value." = center.x.axis %in% c(TRUE, FALSE))
  stopifnot("`plot.title` should be either NULL or character." = is.character(plot.title) | is.null(plot.title))
  stopifnot("`lfc.thresh` should be a numerical vector of length 2." = length(lfc.thresh) == 2 & is.numeric(lfc.thresh))
  stopifnot("`pval.color.thresh` should be a single numerical value." = length(pval.color.thresh) == 1 & is.numeric(pval.color.thresh))
  stopifnot("`show.vlines` should be a boolean value." = show.vlines %in% c(TRUE, FALSE))
  stopifnot("`show.hline` should be a boolean value." = show.hline %in% c(TRUE, FALSE))
  stopifnot("`alpha_transparency` needs to be numerical." = is.numeric(alpha_transparency))
  vlines.linetype <- match.arg(vlines.linetype, choices = c("blank", "solid", "dashed", "dotted", "dotdash", "longdash","twodash"))
  hline.linetype <- match.arg(hline.linetype, choices = c("blank", "solid", "dashed", "dotted", "dotdash", "longdash","twodash"))
  stopifnot("`lfc.colors` should be a color vector of length 3." = sum(areColors(lfc.colors)) == 3 & length(lfc.colors) == 3)
  if(is.null(vlines.colors) == T){
  } else {
    stopifnot("`vlines.colors` should be either NULL, or a color vector of length 2." = sum(areColors(vlines.colors)) == 2 & length(lfc.colors) == 2)
  }
  if(is.null(hline.color) == T){
  } else {
    stopifnot("`hline.color` should be either NULL, or a color vector of length 2." = sum(areColors(hline.color)) == 1 & length(lfc.colors) == 1)
  }
  
  # Create tmp df.
  if(missing(gene.col)){
    data <- data.frame(
      lfc = df.stat[,lfc.col],
      p = -log10(df.stat[,pval.col]),
      p_thresh = df.stat[,pval.col])
  } else {
    data <- data.frame(
      lfc = df.stat[,lfc.col],
      p = -log10(df.stat[,pval.col]),
      p_thresh = df.stat[,pval.col],
      Genes = df.stat[,gene.col])
  }
  
  # Create ggplot arguments based on user inputs. 
  if(missing(y.title)){
    ylabel <- ylab(paste0("-log10(", pval.col,")")) 
  } else {ylabel <- ylab(y.title)}
  
  if(missing(x.title)){
    xlabel <- xlab(lfc.col) 
  } else {xlabel <- xlab(x.title)}
  
  if(missing(y.max)){
    ylimit <- c(0, ceiling(max(data[data$p !=0,deparse(substitute(p))])))
  } else {ylimit <- c(0, y.max)}
  
  if(center.x.axis == T){
    xlimit <- c(-max(abs(ceiling(max(data[,deparse(substitute(lfc))]))), abs(floor(min(data[,deparse(substitute(lfc))])))),
                max(abs(ceiling(max(data[,deparse(substitute(lfc))]))), abs(floor(min(data[,deparse(substitute(lfc))])))))
  } else {
    xlimit <- c(-abs(floor(min(data[,deparse(substitute(lfc))]))),
                abs(ceiling(max(data[,deparse(substitute(lfc))]))))
  }
  
  if(show.vlines == T){
    if(is.null(vlines.colors)){
      vline_min <- geom_vline(xintercept = min(lfc.thresh), linetype = vlines.linetype, color = lfc.colors[1], linewidth = lines.thickness) 
      vline_max <- geom_vline(xintercept = max(lfc.thresh), linetype = vlines.linetype, color = lfc.colors[3], linewidth = lines.thickness)
    } else {
      vline_min <- geom_vline(xintercept = min(lfc.thresh), linetype = vlines.linetype, color = vlines.colors[1], linewidth = lines.thickness) 
      vline_max <- geom_vline(xintercept = max(lfc.thresh), linetype = vlines.linetype, color = vlines.colors[2], linewidth = lines.thickness)
    }
  } else {
    vline_min <- NULL
    vline_max <- NULL
  }
  
  if(show.hline == T){
    if(is.null(hline.color)){
      hline <- geom_hline(yintercept = -log10(pval.color.thresh), linetype = hline.linetype, color = lfc.colors[2], linewidth = lines.thickness)
    } else {
      hline <- geom_hline(yintercept = -log10(pval.color.thresh), linetype = hline.linetype, color = hline.color, linewidth = lines.thickness)
    }
  } else {
    hline <- NULL
  }
  
  # Setup colors for the plot. 
  data$color <- ifelse(data[,deparse(substitute(p))] > -log10(pval.color.thresh) & data[,deparse(substitute(lfc))] < min(lfc.thresh),
                       data$color <- lfc.colors[1],
                       ifelse(data[,deparse(substitute(p))] > -log10(pval.color.thresh) & data[,deparse(substitute(lfc))] > max(lfc.thresh),
                              data$color <- lfc.colors[3],
                              data$color <- lfc.colors[2]))
  data$color <- data$color %>% factor(levels = lfc.colors)
  MyColors <- lfc.colors
  names(MyColors) <- levels(data$color)
  colScale <- scale_color_manual(name = "color", values = MyColors)
  
  if(missing(gene.col)){
    geomlabel <- NULL
  } else {
    if(is.null(label.genes)){
      geomlabel <- NULL
    } else {
      if(is.numeric(label.genes) & length(label.genes) == 1){
        df.highlight <- data[data$p_thresh < pval.color.thresh & abs(data$lfc) > label.genes,]
      }
      if(is.numeric(label.genes) & length(label.genes) == 2){
        df.highlight <- data[data$p_thresh < pval.color.thresh & data$lfc < min(label.genes) | data$p_thresh < pval.color.thresh & data$lfc > max(label.genes),]
      }
      if(is.character(label.genes) & sum(label.genes %in% df.stat[,gene.col]) == length(label.genes)){
        df.highlight <- data[data$Genes %in% label.genes,]
      }
      geomlabel <- geom_label_repel(data = df.highlight, size = label.size, aes(label = Genes, colour = color))
    }
  }
  
  if(show.gene.nbr == F){
    lowest.annot  <- NULL
    highest.annot <- NULL
  } else {
    if(center.x.axis == F){
      lowest.annot  <- geom_text(y = ceiling(max(data[,deparse(substitute(p))]))*0.9, 
                                 x = -(((abs(floor(min(data[,deparse(substitute(lfc))])))-abs(min(lfc.thresh)))/2)+abs(min(lfc.thresh))), 
                                 label = dim(data[data$p_thresh < 0.05 & data$lfc < min(lfc.thresh),])[1], 
                                 size = 7, fontface = "bold", color = lfc.colors[1], check_overlap = T)
      highest.annot <- geom_text(y = ceiling(max(data[,deparse(substitute(p))]))*0.9, 
                                 x = ((abs(ceiling(max(data[,deparse(substitute(lfc))])))-abs(max(lfc.thresh)))/2)+abs(max(lfc.thresh)), 
                                 label = dim(data[data$p_thresh < 0.05 & data$lfc > max(lfc.thresh),])[1], 
                                 size = 7, fontface = "bold", color = lfc.colors[3], check_overlap = T)
    } else {
      lowest.annot  <- geom_text(y = ceiling(max(data[,deparse(substitute(p))]))*0.9, 
                                 x = -(ceiling(max(data[,deparse(substitute(lfc))]))+abs(min(lfc.thresh)))/2, 
                                 label = dim(data[data$p_thresh < 0.05 & data$lfc < min(lfc.thresh),])[1], 
                                 size = 7, fontface = "bold", color = lfc.colors[1], check_overlap = T)
      highest.annot <- geom_text(y = ceiling(max(data[,deparse(substitute(p))]))*0.9, 
                                 x = (ceiling(max(data[,deparse(substitute(lfc))]))+abs(max(lfc.thresh)))/2, 
                                 label = dim(data[data$p_thresh < 0.05 & data$lfc > max(lfc.thresh),])[1], 
                                 size = 7, fontface = "bold", color = lfc.colors[3], check_overlap = T)
    }
  }
  
  # Plot Volcano plot. 
  plot <- ggplot(data, aes(x = lfc, y = p)) +
    geom_point(size = point.size, alpha = alpha_transparency, aes(colour = color)) + 
    ylab(ylabel) +
    xlab(xlabel) +
    ylim(ylimit) +
    xlim(xlimit) +
    ggtitle(plot.title) +
    geomlabel +
    colScale +
    vline_min +
    vline_max +
    hline +
    lowest.annot + 
    highest.annot +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 15, colour = "black"),
          axis.line = element_line(linewidth = 0.5),
          axis.ticks.length = unit(2, "mm"),
          axis.title = element_text(size = 20))
  return(plot)
}

########################## CHECK IF COLORS  ########################## 
areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}

########################## PLOTTING VIOLIN PLOT  ########################## 
# The prepare_data function is designed to process data from a Seurat object (tp), identifying differentially 
# expressed genes between predefined groups. It computes adjusted p-values and assigns significance symbols 
# to genes based on these values. The plot_seq_data function generates violin plots for specified genes across 
# different groups defined by a grouping variable (group.by). Finally, the plot_gene_set function orchestrates 
# the plotting of each violin plots, organizing them into a grid layout and saving the resulting plots directly to PDF files.

# Function to prepare data with customizable identifiers
prepare_data <- function(tp, ident1, ident2) {
  df <- FindMarkers(tp, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0, min.pct = 0)
  df$Genes <- rownames(df)
  
  # Add Symbol to the df object
  df$Symbol <- "ns"
  df$Symbol <- ifelse(df$p_val_adj < 0.0001, "****", 
                      ifelse(df$p_val_adj < 0.001, "***",
                             ifelse(df$p_val_adj < 0.01, "**",
                                    ifelse(df$p_val_adj < 0.05, "*", "ns"))))
  return(df)
}

# Create plot_seq_data function with customizable group.by argument
plot_seq_data <- function(gene, df, tp, group.by = "Clone.ID", cols = c("#339900", "#CC9900")) {
  # Determine caption and title from df
  gene_info <- df[df$Genes == gene, ]
  title_text <- str_to_title(gene_info$Genes)
  caption_text <- gene_info$Symbol
  
  # Create the plot
  plot <- VlnPlot(obj = tp, features = gene, group.by = group.by, pt.size = 0.25, alpha = 0.2, ncol = length(gene), cols = cols) & 
    labs(title = title_text, caption = caption_text) & 
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "italic"),
      plot.caption = element_text(hjust = 0.5, vjust = ifelse(caption_text == "ns", 40, 31), face = "plain", size = ifelse(caption_text == "ns", 16, 20)),
      legend.position = 'none',
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.y = element_blank()
    )
  
  return(plot)
}

# Function to plot gene sets and save to file using ggsave
plot_gene_set <- function(gene_set, file_name, df, tp, group.by = "Clone.ID", cols = c("#339900", "#CC9900")) {
  plots <- lapply(gene_set, function(gene) {
    plot_seq_data(gene, df = df, tp = tp, group.by = group.by, cols = cols)
  })
  
  pdf(file_name, height = 2.5, width = 2 * length(plots))
  grid.arrange(grobs = plots, ncol = length(plots))
  dev.off()
}
