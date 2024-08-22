############## Wrapper around findMarkers() to perform Wilcoxon rank-sum test on a SingleCellExperiment object ############## 
# Parameters:
#   sce.object: The SingleCellExperiment object on which to perform the Wilcoxon test (required).
#   data.type: The assay slot to use for the analysis, e.g., "logcounts" or "normalized_counts" (default: "logcounts").
#   group.by: The column name in `sce.object@colData` that contains the grouping information (required).
#   correction: The multiple testing correction method to use, one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none" (default: "fdr").
#   compareTo: Indicates which group to compare against (1 or 2) (default: 1).
# Returns:
#   A data frame with the results of the Wilcoxon rank-sum test and log-fold changes, including a "Symbol" column to indicate significance levels.

wilcoxon.sce <- function(
    sce.object,
    data.type = "logcounts",
    group.by,
    correction = "fdr",
    compareTo = 1
    ) {
  # Check user-defined arguments.
  if (missing(sce.object)) {
    stop("`sce.object` is missing. Please provide a SingleCellExperiment object.")
  }
  if (missing(group.by)) {
    stop("`group.by` is missing. Please provide the column name in `sce.object@colData` that contains the grouping information.")
  }
  if (compareTo != 1 && compareTo != 2) {
    stop("`compareTo` needs to be either 1 or 2, indicating which group to compare against.")
  }
  
  # Ensure input is valid.
  stopifnot("`sce.object` should be a SingleCellExperiment object." = class(sce.object) == "SingleCellExperiment")
  stopifnot("`data.type` should be a valid assay slot in the SingleCellExperiment object." = data.type %in% names(sce.object@assays))
  stopifnot("`group.by` should be a valid column name in the SingleCellExperiment object." = group.by %in% colnames(sce.object@colData))
  stopifnot("`group.by` should contain only two observations (for the Wilcoxon test)." = length(unique(sce.object[[group.by]])) == 2)
  
  # Check the correction method.
  correction <- match.arg(correction, choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  
  # Format the sce object and remove genes with 0 expression.
  num.expressed.cells <- nexprs(sce.object, byrow = TRUE)
  sce.object <- sce.object[num.expressed.cells > 0, ]
  param <- sce.object[[group.by]] %>% as.factor()
  
  # Run the Wilcoxon test and the t-test (for log-fold changes).
  wilcoxon.results <- findMarkers(sce.object, sce.object[[group.by]], direction = "any", test = "wilcox")
  wilcoxon.top <- wilcoxon.results[[as.character(unique(param)[compareTo])]] %>% as.data.frame()
  wilcoxon.top$Genes <- rownames(wilcoxon.top)
  wilcoxon.top <- wilcoxon.top[, c("Genes", "p.value", "FDR")]
  
  logfc.results <- findMarkers(sce.object, sce.object[[group.by]], direction = "any", test = "t")
  logfc.top <- logfc.results[[as.character(unique(param)[compareTo])]] %>% as.data.frame()
  logfc.top$Genes <- rownames(logfc.top)
  logfc.top <- logfc.top[, grepl("logFC.", names(logfc.top)) | grepl("Genes", names(logfc.top))]
  
  # Merge the Wilcoxon and t-test results.
  merged.results <- left_join(wilcoxon.top, logfc.top, by = "Genes")
  
  # Add a "Symbol" column to indicate significance levels.
  merged.results <- merged.results %>%
    mutate(Symbol = case_when(
      FDR < 0.05 & FDR > 0.01 ~ "*",
      FDR < 0.01 & FDR > 0.001 ~ "**",
      FDR < 0.001 & FDR > 0.0001 ~ "***",
      FDR < 0.0001 ~ "****",
      TRUE ~ format(FDR, digits = 2)
    ))
  merged.results[merged.results$FDR > 0.05, "Symbol"] <- format(merged.results[merged.results$FDR > 0.05, "FDR"], digits = 1)
  
  # Return the results.
  return(merged.results)
}

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

############## Function to create half-violin plot with significance annotations ##############
# Parameters:
#   sce.object: The SingleCellExperiment object containing the expression data (required).
#   gene.vector: A vector of gene names to plot (required).
#   df.wilcox: A data frame containing the results of the Wilcoxon rank-sum test (default: NULL, will be calculated if missing).
#   group.by: The name of the column in the colData of `sce.object` that contains the grouping information (required).
#   plot.title: The title of the plot (default: NULL).
#   title.size: The font size of the plot title (default: 28).
#   title.face: The font face of the plot title (default: "plain").
#   colors: A vector of 2 colors to use for the violin plots (default: c("red", "blue")).
#   show.signif: Whether to show significance annotations (default: TRUE).
#   signif.param: The type of significance annotation to use, one of "stars", "pval", or "starsORpval" (default: "starsORpval").
#   stars.size: The font size of the star significance annotations (default: 10).
#   pval.size: The font size of the p-value significance annotations (default: 6).
#   stars.above.violin: The vertical position of the star significance annotations relative to the top of the violin plot (default: 0.5).
#   pval.above.violin: The vertical position of the p-value significance annotations relative to the top of the violin plot (default: 1.3).
#   show.hlines: Whether to show horizontal lines at fixed expression levels (default: TRUE).
#   hlines.color: The color of the horizontal lines (default: "grey").
#   hlines.linetype: The line type of the horizontal lines (default: "dotted").
#   alpha.transparency: The transparency of the violin plots (default: 0.6).
#   violin.size: The width of the violin plots (default: 0.5).
#   y.scale.limits: The minimum and maximum values for the y-axis (default: c(0, 13)).
#   y.title: The title of the y-axis (default: "Expression Level").
#   y.title.size: The font size of the y-axis title (default: 20).
#   y.text.size: The font size of the y-axis tick labels (default: 12).
#   x.text.size: The font size of the x-axis tick labels (default: 22).
#   x.text.face: The font face of the x-axis tick labels (default: "italic").
#   y.linewidth: The width of the y-axis line (default: 0.5).
# Returns:
#   A ggplot object representing the violin plot with significance annotations.

plot.violin <- function(
    sce.object,
    gene.vector,
    df.wilcox = NULL,
    group.by,
    plot.title = NULL,
    title.size = 28,
    title.face = "plain",
    colors = c("red", "blue"),
    show.signif = TRUE,
    signif.param = "starsORpval",
    stars.size = 10,
    pval.size = 6,
    stars.above.violin = 0.5,
    pval.above.violin = 1.3,
    show.hlines = TRUE,
    hlines.color = "grey",
    hlines.linetype = "dotted",
    alpha.transparency = 0.6,
    violin.size = 0.5,
    y.scale.limits = c(0, 13),
    y.title = "Expression Level",
    y.title.size = 20,
    y.text.size = 12,
    x.text.size = 22,
    x.text.face = "italic",
    y.linewidth = 0.5
) {
  # Input validation.
  stopifnot(is(sce.object, "SingleCellExperiment"))
  stopifnot(is.character(gene.vector) && length(gene.vector) > 0)
  stopifnot(is.null(df.wilcox) || is.data.frame(df.wilcox))
  stopifnot(is.character(group.by) && group.by %in% colnames(colData(sce.object)))
  stopifnot(is.null(plot.title) || is.character(plot.title))
  stopifnot(is.numeric(title.size) && title.size > 0)
  stopifnot(is.character(title.face) && title.face %in% c("plain", "italic", "bold", "bold.italic"))
  stopifnot(is.numeric(stars.size) && stars.size > 0)
  stopifnot(is.numeric(pval.size) && pval.size > 0)
  stopifnot(is.numeric(stars.above.violin) && stars.above.violin >= 0)
  stopifnot(is.numeric(pval.above.violin) && pval.above.violin >= 0)
  stopifnot(is.logical(show.hlines))
  stopifnot(is.character(hlines.color) && length(hlines.color) == 1)
  stopifnot(hlines.linetype %in% c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"))
  stopifnot(is.numeric(alpha.transparency) && alpha.transparency >= 0 && alpha.transparency <= 1)
  stopifnot(is.numeric(violin.size) && violin.size > 0)
  stopifnot(is.numeric(y.scale.limits) && length(y.scale.limits) == 2)
  stopifnot(is.character(y.title) && length(y.title) == 1)
  stopifnot(is.numeric(y.title.size) && y.title.size > 0)
  stopifnot(is.numeric(y.text.size) && y.text.size > 0)
  stopifnot(is.numeric(x.text.size) && x.text.size > 0)
  stopifnot(is.character(x.text.face) && x.text.face %in% c("plain", "italic", "bold", "bold.italic"))
  stopifnot(is.numeric(y.linewidth) && y.linewidth > 0)
  
  # If df.wilcox is not provided, perform a Wilcoxon rank test.
  if (is.null(df.wilcox)) {
    df.wilcox <- wilcoxon.sce(sce.object = sce.object, group.by = group.by, compareTo = 2)
  }
  
  # Prepare the data for plotting.
  tmp <- as.data.frame(assay(sce.object, "logcounts"))
  tmp <- tmp[rownames(tmp) %in% gene.vector, ] %>% t() %>% as.data.frame()
  
  # Create a data frame for plotting.
  data <- data.frame(
    gene   = str_to_title(rep(gene.vector, each = dim(tmp)[1])) %>%
      factor(levels = str_to_title(df.wilcox[order(df.wilcox$Genes %in% gene.vector, decreasing = FALSE), "Genes"])),
    expr   = as.vector(sapply(gene.vector, function(x) tmp[, x])),
    param  = rep(sce.object[[group.by]] %>% as.factor(), dim(tmp)[2]),
    ymax   = rep(as.vector(sapply(gene.vector, function(x) max(tmp[, x]))), each = dim(tmp)[1]),
    FDR    = rep(as.vector(sapply(gene.vector, function(x) df.wilcox[df.wilcox$Genes == x, "FDR"])), each = dim(tmp)[1]),
    Symbol = rep(as.vector(sapply(gene.vector, function(x) df.wilcox[df.wilcox$Genes == x, "Symbol"])), each = dim(tmp)[1])
  )
  data$gene <- droplevels(data$gene)
  
  # Prepare significance labels if requested.
  if (show.signif == TRUE) {
    if (signif.param == "stars") {
      data$Symbol[str_detect(data$Symbol, pattern = "\\*", negate = TRUE)] <- "ns"
    }
    if (signif.param == "pval") {
      data$Symbol[str_detect(data$Symbol, pattern = "\\*")] <- format(data$FDR[str_detect(data$Symbol, pattern = "\\*")], digits = 2)
    }
    if (signif.param == "starsORpval") {
    }
    
    # Adjust significance label size and position.
    data$pvalsize <- ifelse(str_detect(data$Symbol, pattern = "\\*"), stars.size, pval.size)
    data$abovepval <- ifelse(str_detect(data$Symbol, pattern = "\\*"), stars.above.violin, pval.above.violin)
    
    # Create significance label layer.
    signif.label <- geom_text(aes(x = gene, y = ymax + abovepval, label = Symbol),
                              check_overlap = TRUE, size = data$pvalsize)
  } else {
    signif.label <- NULL
  }
  
  # Plot violin plots.
  plot <- ggplot(data = data, aes(x = gene, y = expr, fill = param)) +
    scale_y_continuous(limits = function(r) { c(min(r), (max(r) + 0.1)) },
                       breaks = function(z) seq(0, range(z)[2], by = 2.5)) +
    geom_hline(yintercept = seq(0, ceiling(max(data$expr) + 1), by = 2.5),
               color = hlines.color, linetype = hlines.linetype) +
    geom_split_violin(alpha = alpha.transparency, trim = TRUE, scale = "width", size = violin.size) +
    scale_fill_manual(values = colors) +
    ggtitle(plot.title) +
    ylab(y.title) +
    signif.label +
    theme(
      axis.title.y = element_text(size = y.title.size),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = y.text.size, colour = "black"),
      axis.text.x = element_text(size = x.text.size, face = x.text.face, colour = "black"),
      axis.ticks.length.y = unit(2, "mm"),
      axis.ticks.x = element_blank(),
      axis.line.y = element_line(linewidth = y.linewidth),
      panel.background = element_blank(),
      plot.title = element_text(size = title.size, hjust = 0.5, face = title.face)
    )
  
  return(plot)
}

############## Wrapper around the VlnPlot function from seurat ############## 
# This function will use the df dataframe created before (with the wilcoxon.sce function) to
# add significance annotations to each violin plot.
# Parameters:
#   gene: A vector of genes to plot (required).
# Returns:
#   A ggplot object representing the violin plot for the specified genes, with the gene name and
#   significance annotation added as a title and caption, respectively. The function also
#   customizes the appearance of the plot, including removing the legend, x-axis labels and
#   title, and y-axis title.

plot_seq_data <- function(gene) {
  VlnPlot(obj = data_seurat, features = gene, group.by = "TP", pt.size=0.5, ncol = length(gene), cols = c("red", "blue")) & 
    labs(title = str_to_title(df[df[,"Genes"] %in% gene,]$Genes),
         caption = df[df[,"Genes"] %in% gene,]$Symbol) & 
    theme(plot.title = element_text(hjust= 0.5, size = 13, face = "italic"),
          plot.caption = element_text(hjust = 0.5, vjust = 23, face = "plain", size = 20),
          legend.position = 'none', 
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.line.x=element_blank(),
          axis.title.y = element_blank())
}

############## Check whether elements of a vector are colors ##############

areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}

############## Opposite of the %in% function ##############

`%notin%` <- Negate(`%in%`)
