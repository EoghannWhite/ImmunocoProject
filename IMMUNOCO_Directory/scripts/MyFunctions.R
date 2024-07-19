########################## CHECKING WHETHER ELEMENTS OF A VECTOR ARE COLORS  ########################## 
areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}

########################## OPPOSITE OF THE %in% FUNCTION  ########################## 
`%notin%` <- Negate(`%in%`)

########################## PERFORMING WILCOXON TEST ON SCE OBJECT  ########################## 
# Function to perform wilcoxon test and calculate logFC on a SingleCellExperiment object.

wilcoxon.sce <- function(sce.object,
                         data.type = "logcounts",
                         group.by,
                         correction = "fdr",
                         compareTo = 1)
  {
  # Checking user-defined arguments.
  if(missing(sce.object)){stop("`sce.object` is missing. Please provide a SingleCellExperiment object.")}
  if(missing(group.by)){stop("`group.by` is missing. Please provide it.")}
  if(compareTo != 1 & compareTo !=2){stop("`compareTo` needs to be a numerical value : 1 or 2.")}  
  stopifnot("`sce.object` should be a SingleCellExperiment object." = class(sce.object)=="SingleCellExperiment")
  stopifnot("`data.type` should be a valid assay slot in the SingleCellExperiment object." = data.type %in% names(sce.object@assays))
  stopifnot("`group.by` should be a valid colname in the SingleCellExperiment object." = group.by %in% colnames(sce.object@colData))
  stopifnot("`group.by` should contain only two observations (for the wilcoxon test)." = length(unique(sce.object[[group.by]])) == 2)
  correction <- match.arg(correction, choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
  # Format sce object and remove genes with 0 expression to reduce computation time.
  num.cells <- nexprs(sce.object, byrow=TRUE)
  sce.object <- sce.object[num.cells>0,]
  param <- sce.object[[group.by]] %>% as.factor()
  # Running wilcoxon test with the `findMarkers` function.
  wilc <- findMarkers(sce.object, sce.object[[group.by]], direction = "any", test = "wilcox")
  wilc.top <- wilc[[as.character(unique(param)[compareTo])]] %>% as.data.frame()
  wilc.top$Genes <- rownames(wilc.top)
  wilc.top <- wilc.top[, c("Genes", "p.value", "FDR")]
  lfc <- findMarkers(sce.object, sce.object[[group.by]], direction = "any", test = "t")
  lfc.top <- lfc[[as.character(unique(param)[compareTo])]] %>% as.data.frame()
  lfc.top$Genes <- rownames(lfc.top)
  lfc.top <- lfc.top[, grepl("logFC.", names(lfc.top)) | grepl("Genes", names(lfc.top))]
  merged <- left_join(wilc.top, lfc.top, by = "Genes")
  merged <- merged %>%
    mutate(Symbol = case_when(
      FDR < 0.05 & FDR > 0.01 ~ "*",
      FDR < 0.01 & FDR > 0.001 ~ "**",
      FDR < 0.001 & FDR > 0.0001 ~ "***",
      FDR < 0.0001 ~ "****",
      TRUE ~ format(FDR, digits = 2)
    ))
  merged[merged$FDR > 0.05, "Symbol"] <- format(merged[merged$FDR > 0.05, "FDR"], digits = 1)
  # Return result.
  return(merged)
}

########################## PLOTTING VOLCANO PLOT  ########################## 
# This function creates a volcano plot starting from any dataframe containing logFC 
# and corresponding p values.

plot.volcano <- function(df.stat, pval.col = "FDR", lfc.col, gene.col = "Genes",
                         point.size = 2, plot.title = NULL,
                         y.title, x.title, y.max, x.min, x.max, center.x.axis = T,
                         lfc.colors = c("red", "#CCCCCC", "#0066CC"), lfc.thresh = c(-1,1),
                         pval.color.thresh = 0.05, alpha_transparency = 0.6,
                         show.vlines = T, vlines.linetype = "dashed", vlines.colors = NULL,
                         show.hline = T, hline.linetype = "dashed", hline.color = NULL,
                         lines.thickness = 0.5,
                         label.genes = NULL, label.size = 4, show.gene.nbr = F)
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
                                 size = 7, fontface = "plain", color = lfc.colors[1], check_overlap = T)
      highest.annot <- geom_text(y = ceiling(max(data[,deparse(substitute(p))]))*0.9, 
                                 x = ((abs(ceiling(max(data[,deparse(substitute(lfc))])))-abs(max(lfc.thresh)))/2)+abs(max(lfc.thresh)), 
                                 label = dim(data[data$p_thresh < 0.05 & data$lfc > max(lfc.thresh),])[1], 
                                 size = 7, fontface = "plain", color = lfc.colors[3], check_overlap = T)
    } else {
      lowest.annot  <- geom_text(y = ceiling(max(data[,deparse(substitute(p))]))*0.9, 
                                 x = -(ceiling(max(data[,deparse(substitute(lfc))]))+abs(min(lfc.thresh)))/2, 
                                 label = dim(data[data$p_thresh < 0.05 & data$lfc < min(lfc.thresh),])[1], 
                                 size = 7, fontface = "plain", color = lfc.colors[1], check_overlap = T)
      highest.annot <- geom_text(y = ceiling(max(data[,deparse(substitute(p))]))*0.9, 
                                 x = (ceiling(max(data[,deparse(substitute(lfc))]))+abs(max(lfc.thresh)))/2, 
                                 label = dim(data[data$p_thresh < 0.05 & data$lfc > max(lfc.thresh),])[1], 
                                 size = 7, fontface = "plain", color = lfc.colors[3], check_overlap = T)
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

########################## PLOTTING HALF-VIOLIN PLOTS ########################## 
# This function creates violin plots of a chosen set of genes. 
# The user need to give a character vector of genes that will be retrieved directly from the wilcoxon.sce function output.

plot.violin <- function(sce.object,
                        gene.vector,
                        df.wilcox,
                        group.by,
                        plot.title = NULL,
                        title.size = 28,
                        title.face = "plain",
                        colors = c("red", "blue"),
                        # Significativity
                        show.signif = T,
                        signif.param = "starsORpval",
                        stars.size = 10,
                        pval.size = 6,
                        stars.above.violin = 0.5,
                        pval.above.violin = 1.3,
                        # Background horizontal lines
                        show.hlines = T,
                        hlines.color = "grey",
                        hlines.linetype = "dotted",
                        # violin arguments
                        alpha.transparency = 0.6,
                        violin.size = 0.5,
                        # Theme arguments
                        y.scale.limits = c(0, 13),
                        y.title = "Expression Level",
                        y.title.size = 20,
                        y.text.size = 12,
                        x.text.size = 22,
                        x.text.face = "italic",
                        y.linewidth = 0.5
                        ){
  
  # Perform wilcoxon rank test with wilcoxon.sce if df.wilcox is missing.
  if(missing(df.wilcox)){
    df.wilcox <- wilcoxon.sce(sce.object = sce.object, group.by = group.by, compareTo = 2)
  }
  # Prepare data for plotting.
  tmp <- as.data.frame(assay(sce.object, "logcounts"))
  tmp <- tmp[rownames(tmp) %in% gene.vector,] %>% t() %>% as.data.frame()
  data <- data.frame(
    gene   = str_to_title(rep(gene.vector, each = dim(tmp)[1])) %>% factor(levels = str_to_title(df.wilcox[order(df.wilcox$Genes %in% gene.vector, decreasing = F), "Genes"])),
    expr   = as.vector(sapply(gene.vector, function(x) tmp[,x])),
    param  = rep(sce.object[[group.by]] %>% as.factor(), dim(tmp)[2]),
    ymax   = rep(as.vector(sapply(gene.vector, function(x) max(tmp[,x]))), each = dim(tmp)[1]),
    FDR    = rep(as.vector(sapply(gene.vector, function(x) df.wilcox[df.wilcox$Genes == x, "FDR"])), each = dim(tmp)[1]),
    Symbol = rep(as.vector(sapply(gene.vector, function(x) df.wilcox[df.wilcox$Genes == x, "Symbol"])), each = dim(tmp)[1])
  )
  data$gene <- droplevels(data$gene)
  # Prepare ggplot arguments
  if(show.signif == T){
    if(signif.param == "stars"){
      data$Symbol[str_detect(data$Symbol, pattern = "\\*", negate = T)] <- "ns"
    }
    if(signif.param == "pval"){
      data$Symbol[str_detect(data$Symbol, pattern = "\\*")] <- format(data$FDR[str_detect(data$Symbol, pattern = "\\*")], digits = 2)
    }
    if(signif.param == "starsORpval"){
    }
    data$pvalsize <- data$Symbol
    data$abovepval <- data$Symbol
    data$pvalsize <- ifelse(str_detect(data$Symbol, pattern = "\\*"), data$pvalsize <- stars.size, data$pvalsize <- pval.size)
    data$abovepval <- ifelse(str_detect(data$Symbol, pattern = "\\*"), data$abovepval <- stars.above.violin, data$abovepval <- pval.above.violin)
    signif.label <- geom_text(aes(x = gene, y = ymax + abovepval, label = Symbol), check_overlap = T, size = data$pvalsize)
  } else { 
    signif.label <- NULL
  }
# Plot
  plot <- ggplot(data = data, aes(x = gene, y = expr, fill = param)) +
    scale_y_continuous(limits = function(r){c(min(r),(max(r)+0.1))}, breaks = function(z) seq(0, range(z)[2], by = 2.5)) +
    geom_hline(yintercept = seq(0, ceiling(max(data$expr)+1), by = 2.5), color = hlines.color, linetype = hlines.linetype) +
    geom_split_violin(alpha = alpha.transparency, trim = T, scale = "width", size = violin.size) +
    scale_fill_manual(values = colors) +
    ggtitle(plot.title) +
    ylab(y.title) +
    signif.label +
    theme(axis.title.y = element_text(size = y.title.size),
          axis.title.x = element_blank(), 
          axis.text.y = element_text(size = y.text.size, colour = "black"),
          axis.text.x = element_text(size = x.text.size, face = x.text.face, colour = "black"),
          axis.ticks.length.y = unit(2, "mm"),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(linewidth = y.linewidth),
          panel.background = element_blank(),
          plot.title = element_text(size = title.size, hjust = 0.5, face = title.face))
  return(plot)
}

########################## WRAPPER AROUND THE VlnPlot FUNCTION FROM SEURAT ########################## 
# This function will use the df dataframe created before (with the wilcoxon.sce function) to 
# add significativity to each violin plot. 

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
