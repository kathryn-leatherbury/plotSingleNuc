#' FeaturePlot_scCustom + DotPlot_scCustom composite for Seurat objects
#'
#' Builds a panel of one or more UMAP feature plots (via
#' \code{scCustomize::FeaturePlot_scCustom}) next to a dot plot
#' (\code{scCustomize::DotPlot_scCustom}). You can control the number of
#' UMAP columns, flip/stack the layout, and style axis labels.
#'
#' Note: This function requires \pkg{scCustomize} at runtime.
#' Install with: \code{remotes::install_github("samuel-marsh/scCustomize")}.
#'
#' @param object A \pkg{Seurat} object.
#' @param features Character vector of feature/gene names to plot.
#' @param label Logical; add text labels on the UMAP feature plots.
#' @param repel Logical; use repelled labels for UMAP (if available).
#' @param reduction Character; embedding/reduction name for UMAP plots
#'   (e.g., \code{"umap"}, \code{"umap.wnn"}). Default \code{"umap.wnn"}.
#' @param pt.size Numeric; point size for UMAP feature plots.
#' @param assay Character; assay to use. If \code{NULL}, uses
#'   \code{Seurat::DefaultAssay(object)}.
#' @param group.by Character scalar naming a column in \code{object@meta.data}
#'   (or a factor/character vector of length \code{ncol(object)}) used to group
#'   cells in the dot plot. If \code{NULL}, defaults to \code{Seurat::Idents(object)}.
#' @param split.by Optional character scalar naming a metadata column to split
#'   the dot plot.
#' @param legend Logical; show color legend on the UMAP feature plots.
#' @param umap_axes Logical; show axes on the UMAP feature plots.
#' @param bold_umap_labels Logical; attempt to bold the UMAP text labels
#'   when \code{label = TRUE}.
#' @param featureplot.label.size Numeric; text size for UMAP labels.
#' @param dotplot.label.size Numeric; axis text size for the dot plot.
#' @param bold.y.axis Logical; bold-face the y-axis tick labels in the dot plot.
#' @param bold.x.axis Logical; bold-face the x-axis tick labels in the dot plot.
#' @param colors Character vector of colors for both UMAP and dot-plot gradients.
#'   If \code{NULL}, a pink palette is used. Provide \eqn{\ge} 3 colors for smooth
#'   gradients.
#' @param flip.axes Logical; if \code{TRUE}, the dot plot is rotated (via
#'   \code{Seurat::RotatedAxis()}) and the layout is stacked vertically
#'   (UMAP panel on top, dot plot below). If \code{FALSE}, panels are side-by-side.
#' @param umaps.ncols Integer; number of columns to arrange the UMAP feature
#'   plots panel. Default \code{2}.
#' @param plot.heights Numeric vector of length 2; relative heights
#'   \code{c(UMAPs, DotPlot)} used when \code{flip.axes = TRUE} (stacked layout).
#'   If \code{NULL}, defaults to \code{c(4, 1)}.
#' @param plot.widths Numeric vector of length 2; relative widths
#'   \code{c(UMAPs, DotPlot)} used when \code{flip.axes = FALSE}
#'   (side-by-side layout). If \code{NULL}, defaults to \code{c(4, 1)}.
#'
#' @return A \pkg{patchwork} object combining the UMAP feature plots and the dot plot.
#'
#' @examples
#' \dontrun{
#' p <- FeatureDotPlot(
#'   object = hypo,
#'   features = c("myrf","olig1","gpr17","spi1"),
#'   reduction = "umap.wnn",
#'   label = TRUE, repel = TRUE,
#'   pt.size = 0.3,
#'   assay = "SCT",
#'   group.by = "subclusters_2.1.1_new",
#'   legend = FALSE,
#'   umap_axes = FALSE,
#'   bold_umap_labels = TRUE,
#'   dotplot.label.size = 8,
#'   featureplot.label.size = 3,
#'   bold.y.axis = TRUE,
#'   bold.x.axis = FALSE,
#'   umaps.ncols = 2,                 # arrange UMAPs in 2 columns
#'   flip.axes = FALSE,               # side-by-side layout
#'   plot.widths = c(4, 1)            # relative widths for side-by-side
#' )
#' p
#'
#' # stacked layout (UMAPs on top, dot plot below)
#' p2 <- FeatureDotPlot(
#'   object = hypo,
#'   features = c("myrf","olig1","gpr17","spi1"),
#'   flip.axes = TRUE,
#'   plot.heights = c(3, 2)           # relative heights when stacked
#' )
#' }
#'
#' @importFrom ggplot2 theme margin element_text
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat Idents NoAxes NoLegend DefaultAssay RotatedAxis
#' @export
FeatureDotPlot <- function(
    object,
    features,
    label = TRUE,
    repel = FALSE,
    reduction = "umap.wnn",
    pt.size = 0.3,
    assay = NULL,                 # default to object's DefaultAssay
    group.by = NULL,
    split.by = NULL,
    legend = TRUE,
    umap_axes = TRUE,
    bold_umap_labels = TRUE,
    featureplot.label.size = 3,
    dotplot.label.size = 10,
    bold.y.axis = FALSE,
    bold.x.axis = FALSE,
    colors = NULL,
    flip.axes = FALSE,            # if TRUE → add RotatedAxis() and stack
    umaps.ncols = 2,            # <--- NEW: columns for FeaturePlot_scCustom panel
    plot.heights = NULL,          # <--- NEW: used when flip.axes = TRUE (stacked)
    plot.widths  = NULL           # <--- NEW: used when flip.axes = FALSE (side-by-side)
) {
  suppressWarnings({
    # ---- assay handling ----
    if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
    old_assay <- Seurat::DefaultAssay(object)
    on.exit(Seurat::DefaultAssay(object) <- old_assay, add = TRUE)
    Seurat::DefaultAssay(object) <- assay
    
    # ---- palette ----
    pink_palette <- c("gray87", "pink", "palevioletred1", "deeppink", "deeppink3", "deeppink4")
    if (is.null(colors)) colors <- pink_palette
    if (!is.null(colors) && length(colors) < 3)
      warning("`colors` has < 3 entries; consider providing >= 3 for smoother gradients.")
    
    # ---- features (no cap) ----
    feat_vec <- as.character(features)
    
    # ---- group.by default ----
    if (is.null(group.by)) group.by <- Seurat::Idents(object)
    
    # ---- UMAP FeaturePlot(s) ----
    feat_list <- lapply(feat_vec, function(feature) {
      p <- scCustomize::FeaturePlot_scCustom(
        seurat_object = object,
        features      = feature,
        order         = TRUE,
        label         = label,
        repel         = repel,
        reduction     = reduction,
        colors_use    = colors,
        na_cutoff     = 0.05,
        pt.size       = pt.size,
        label.size    = featureplot.label.size
      )
      if (!umap_axes) p <- p + NoAxes()
      if (!legend)    p <- p + NoLegend()
      if (bold_umap_labels && isTRUE(label)) {
        if (length(p$layers) >= 2 && !is.null(p$layers[[2]]$aes_params)) {
          p$layers[[2]]$aes_params$fontface <- "bold"
        }
      }
      p
    })
    
    # If only one feature and user wants multi-column panel, add a spacer for symmetry
    if (length(feat_list) == 1 && umaps.ncols > 1) {
      spacer <- patchwork::plot_spacer() + ggplot2::theme_void()
      feat_list <- list(feat_list[[1]], spacer)
    }
    
    # Combine feature plots into panel 'a' with requested ncols
    a <- patchwork::wrap_plots(feat_list, ncol = umaps.ncols)
    
    # ---- DotPlot panel ----
    face <- if (isTRUE(bold.y.axis)) "bold" else "plain"
    face <- if (isTRUE(bold.x.axis)) "bold" else "plain"
    b <- scCustomize::DotPlot_scCustom(
      object,
      features   = feat_vec,
      dot.scale  = 5.5,
      col.min    = 0.05,
      dot.min    = 0.01,
      scale      = TRUE,
      colors_use = colors,
      assay      = assay,
      group.by   = group.by,
      split.by   = split.by,
      flip_axes  = flip.axes
    ) + ggplot2::theme(
      plot.margin = ggplot2::margin(t = 0, r = 10, b = 25, l = 10, unit = "pt"),
      axis.text.y = ggplot2::element_text(size = dotplot.label.size, face = face),
      axis.text.x = ggplot2::element_text(size = dotplot.label.size, face = face)
    )
    
    # ---- Final layout with user-controllable sizes ----
    if (isTRUE(flip.axes)) {
      b <- b + RotatedAxis()
      if (is.null(plot.heights)) plot.heights <- c(4, 1)  # default when stacked
      c <- patchwork::wrap_plots(a, b, ncol = 1, heights = plot.heights)
    } else {
      if (is.null(plot.widths)) plot.widths <- c(4, 1)    # default side-by-side
      c <- patchwork::wrap_plots(a, b, ncol = 2, widths = plot.widths)
    }
    
    return(c)
  })
}