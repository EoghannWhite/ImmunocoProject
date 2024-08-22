############## Function to create a volcano plot from a data frame ##############
# Parameters:
#   df.stat: The data frame containing the data to plot (required). Should have columns for p-values and log-fold changes.
#   pval.col: The name of the column containing the p-values (default: "FDR").
#   lfc.col: The name of the column containing the log-fold changes (required).
#   gene.col: The name of the column containing the gene names (default: "Genes").
#   point.size: The size of the points in the plot (default: 2).
#   plot.title: The title of the plot (default: NULL, no title).
#   y.title: The title for the y-axis (default: NULL, uses "-log10(pval.col)").
#   x.title: The title for the x-axis (default: NULL, uses lfc.col).
#   y.max: The maximum value for the y-axis (default: NULL, uses the maximum value in the data).
#   x.min: The minimum value for the x-axis (default: NULL, uses the minimum value in the data).
#   x.max: The maximum value for the x-axis (default: NULL, uses the maximum value in the data).
#   center.x.axis: Whether to center the x-axis around 0 (default: TRUE).
#   lfc.colors: A vector of 3 colors to use for the different log-fold change regions (default: c("red", "#CCCCCC", "#0066CC")).
#   lfc.thresh: A vector of 2 values specifying the log-fold change thresholds (default: c(-1, 1)).
#   pval.color.thresh: The p-value threshold for coloring points (default: 0.05).
#   alpha_transparency: The transparency of the points (default: 0.6).
#   show.vlines: Whether to show vertical lines at the log-fold change thresholds (default: TRUE).
#   vlines.linetype: The line type for the vertical lines (default: "dashed").
#   vlines.colors: The colors for the vertical lines (default: NULL, uses lfc.colors).
#   show.hline: Whether to show a horizontal line at the p-value threshold (default: TRUE).
#   hline.linetype: The line type for the horizontal line (default: "dashed").
#   hline.color: The color for the horizontal line (default: NULL, uses lfc.colors[2]).
#   lines.thickness: The thickness of the lines (default: 0.5).
#   label.genes: A numeric value (number of genes to label), a numeric vector of 2 values (log-fold change range to label), or a character vector of gene names to label (default: NULL, no labeling).
#   label.size: The size of the gene labels (default: 4).
#   show.gene.nbr: Whether to show the number of genes above/below the log-fold change thresholds (default: FALSE).
# Returns:
#   A ggplot object representing the volcano plot.

plot.volcano <- function(
    df.stat, 
    pval.col = "FDR", 
    lfc.col, 
    gene.col = "Genes",
    point.size = 2, 
    plot.title = NULL,
    y.title = NULL, 
    x.title = NULL, 
    y.max = NULL, 
    x.min = NULL, 
    x.max = NULL, 
    center.x.axis = TRUE,
    lfc.colors = c("red", "#CCCCCC", "#0066CC"), 
    lfc.thresh = c(-1, 1),
    pval.color.thresh = 0.05, 
    alpha_transparency = 0.6,
    show.vlines = TRUE, 
    vlines.linetype = "dashed", 
    vlines.colors = NULL,
    show.hline = TRUE, 
    hline.linetype = "dashed", 
    hline.color = NULL,
    lines.thickness = 0.5,
    label.genes = NULL, 
    label.size = 4, 
    show.gene.nbr = FALSE
) {
  # Input validation.
  stopifnot(!is.null(df.stat) && nrow(df.stat) > 0)
  stopifnot(is.data.frame(df.stat), is.numeric(df.stat[[pval.col]]), is.numeric(df.stat[[lfc.col]]))
  stopifnot(all(c(pval.col, lfc.col, gene.col) %in% names(df.stat)))
  stopifnot(is.character(df.stat[[gene.col]]) | missing(gene.col))
  if (!missing(gene.col)) {
    stopifnot(anyDuplicated(df.stat[[gene.col]]) == 0)
  }
  stopifnot(is.logical(center.x.axis), is.null(plot.title) | is.character(plot.title))
  stopifnot(is.numeric(lfc.thresh) & length(lfc.thresh) == 2)
  stopifnot(is.numeric(pval.color.thresh) & length(pval.color.thresh) == 1)
  stopifnot(is.logical(show.vlines), is.logical(show.hline))
  stopifnot(is.numeric(alpha_transparency))
  stopifnot(is.null(label.genes) || is.numeric(label.genes) || is.character(label.genes))
  
  # Ensure linetype is valid.
  valid_linetypes <- c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  vlines.linetype <- match.arg(vlines.linetype, choices = valid_linetypes)
  hline.linetype <- match.arg(hline.linetype, choices = valid_linetypes)
  
  # Color validation.
  stopifnot(sum(areColors(lfc.colors)) == 3 & length(lfc.colors) == 3)
  
  if (!is.null(vlines.colors)) {
    stopifnot(sum(areColors(vlines.colors)) == 2 & length(vlines.colors) == 2)
  }
  
  if (!is.null(hline.color)) {
    stopifnot(sum(areColors(hline.color)) == 1 & length(hline.color) == 1)
  }
  
  # Create temporary data frame.
  data <- data.frame(
    lfc = df.stat[[lfc.col]],
    p = -log10(df.stat[[pval.col]]),
    p_thresh = df.stat[[pval.col]]
  )
  if (!missing(gene.col)) {
    data$Genes <- df.stat[[gene.col]]
  }
  
  # Axis titles and limits.
  xlabel <- x.title %||% lfc.col
  ylabel <- y.title %||% paste0("-log10(", pval.col, ")")
  
  ylimit <- c(0, y.max %||% ceiling(max(data$p[data$p != 0])))
  
  xlimit <- if (center.x.axis) {
    max_abs_lfc <- max(abs(ceiling(max(data$lfc))), abs(floor(min(data$lfc))))
    c(-max_abs_lfc, max_abs_lfc)
  } else {
    c(-abs(floor(min(data$lfc))), abs(ceiling(max(data$lfc))))
  }
  
  # Vertical and horizontal lines.
  vlines <- if (show.vlines) {
    vlines.colors <- vlines.colors %||% c(lfc.colors[1], lfc.colors[3])
    list(
      geom_vline(xintercept = min(lfc.thresh), linetype = vlines.linetype, color = vlines.colors[1], linewidth = lines.thickness),
      geom_vline(xintercept = max(lfc.thresh), linetype = vlines.linetype, color = vlines.colors[2], linewidth = lines.thickness)
    )
  } else {
    NULL
  }
  
  hline <- if (show.hline) {
    geom_hline(yintercept = -log10(pval.color.thresh), linetype = hline.linetype, color = hline.color %||% lfc.colors[2], linewidth = lines.thickness)
  } else {
    NULL
  }
  
  # Color mapping.
  data$color <- factor(
    ifelse(
      data$p > -log10(pval.color.thresh) & data$lfc < min(lfc.thresh), lfc.colors[1],
      ifelse(data$p > -log10(pval.color.thresh) & data$lfc > max(lfc.thresh), lfc.colors[3], lfc.colors[2])
    ),
    levels = lfc.colors
  )
  col_scale <- scale_color_manual(name = "color", values = setNames(lfc.colors, lfc.colors))
  
  # Label genes.
  geomlabel <- if (!is.null(label.genes)) {
    if (is.numeric(label.genes)) {
      df.highlight <- data[data$p_thresh < pval.color.thresh & abs(data$lfc) > label.genes, ]
    } else if (is.character(label.genes)) {
      df.highlight <- data[data$Genes %in% label.genes, ]
    }
    geom_label_repel(data = df.highlight, size = label.size, aes(label = Genes, colour = color))
  } else {
    NULL
  }
  
  # Annotate number of genes.
  annot_text <- function(y_pos, x_pos, label, color) {
    geom_text(
      y = y_pos,
      x = x_pos,
      label = label,
      size = 7, fontface = "plain", color = color, check_overlap = TRUE
    )
  }
  
  lowest.annot <- if (show.gene.nbr) {
    y_pos <- ceiling(max(data$p)) * 0.9
    x_pos <- -(ceiling(max(data$lfc)) + abs(min(lfc.thresh))) / 2
    label <- nrow(data[data$p_thresh < 0.05 & data$lfc < min(lfc.thresh), ])
    annot_text(y_pos, x_pos, label, lfc.colors[1])
  } else {
    NULL
  }
  
  highest.annot <- if (show.gene.nbr) {
    y_pos <- ceiling(max(data$p)) * 0.9
    x_pos <- (ceiling(max(data$lfc)) + abs(max(lfc.thresh))) / 2
    label <- nrow(data[data$p_thresh < 0.05 & data$lfc > max(lfc.thresh), ])
    annot_text(y_pos, x_pos, label, lfc.colors[3])
  } else {
    NULL
  }
  
  # Plot volcano plot. 
  ggplot(data, aes(x = lfc, y = p)) +
    geom_point(size = point.size, alpha = alpha_transparency, aes(colour = color)) +
    ylab(ylabel) +
    xlab(xlabel) +
    ylim(ylimit) +
    xlim(xlimit) +
    ggtitle(plot.title) +
    geomlabel +
    col_scale +
    vlines +
    hline +
    lowest.annot +
    highest.annot +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text = element_text(size = 15, colour = "black"),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks.length = unit(2, "mm"),
      axis.title = element_text(size = 20)
    )
}

############## Data Preparation for Differential Expression Analysis ##############
# The prepare_data function is designed to process data from a Seurat object (tp), identifying differentially 
# expressed genes between predefined groups. It computes adjusted p-values and assigns significance symbols 
# to genes based on these values. 
# Parameters:
#   tp: A Seurat object containing the single-cell expression data (required).
#   ident1: The identifier for the first group to compare (required).
#   ident2: The identifier for the second group to compare (required).
# Returns:
#   A dataframe containing:
#     - Differential expression statistics (log fold change, p-values, etc.)
#     - Gene names
#     - Significance symbols based on adjusted p-values

prepare_data <- function(tp, ident1, ident2) {
  df <- FindMarkers(tp, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0, min.pct = 0)
  df$Genes <- rownames(df)
  df$Symbol <- "ns"
  df$Symbol <- ifelse(df$p_val_adj < 0.0001, "****", 
                      ifelse(df$p_val_adj < 0.001, "***",
                             ifelse(df$p_val_adj < 0.01, "**",
                                    ifelse(df$p_val_adj < 0.05, "*", "ns"))))
  return(df)
}

############## Wrapper around the VlnPlot function from seurat ############## 
# This function will use the dataframe created with the prepare_data function to
# generates violin plots for specified genes across different groups defined by a grouping variable (group.by).
# Parameters:
#   gene: A vector of genes to plot (required).
#   df: A dataframe containing statistical results, including gene names, significance levels, and symbols (required).
#   tp: A Seurat object containing the single-cell expression data (required).
#   group.by: The variable to group the data by (default: "Clone.ID").
#   cols: A vector of colors for the violin plots (default: c("#339900", "#CC9900")).
# Returns:
#   A ggplot object representing the violin plot for the specified genes, with the gene name and
#   significance annotation added as a title and caption, respectively. The function also
#   customizes the appearance of the plot, including removing the legend, x-axis labels and
#   title, and y-axis title.

plot_seq_data <- function(
    gene, 
    df, 
    tp, 
    group.by = "Clone.ID", 
    cols = c("#339900", "#CC9900")
    ) {
  
  # Determine caption and title from df.
  gene_info <- df[df$Genes == gene, ]
  title_text <- str_to_title(gene_info$Genes)
  caption_text <- gene_info$Symbol
  
  # Create the plot.
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

############## Plot genes as violin plot and export using ggsave ############## 
# The plot_gene_set function orchestrates the creation of violin plots for a set of genes,
# organizes them into a grid layout, and saves the resulting plots directly to PDF files.
# Parameters:
#   gene_set: A vector of gene names to be plotted (required).
#   file_name: The name of the output PDF file (required).
#   df: A dataframe containing gene expression data and statistics (required).
#   tp: A Seurat object containing the single-cell expression data (required).
#   group.by: The variable to group the data by in the violin plots (default: "Clone.ID").
#   cols: A vector of colors for the violin plots (default: c("#339900", "#CC9900")).
# Returns:
#   This function does not return a value. Instead, it produces a side effect:
#   creating and saving a PDF file containing the violin plots.

plot_gene_set <- function(gene_set, file_name, df, tp, group.by = "Clone.ID", cols = c("#339900", "#CC9900")) {
  
  plots <- lapply(gene_set, function(gene) {
    plot_seq_data(gene, df = df, tp = tp, group.by = group.by, cols = cols)
  })
  
  pdf(file_name, height = 2.5, width = 2 * length(plots))
  grid.arrange(grobs = plots, ncol = length(plots))
  dev.off()
}

############## Check whether elements of a vector are colors ############## 

areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}