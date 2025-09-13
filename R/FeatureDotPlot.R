#' FeaturePlot_scCustomize + DotPlot_scCustomize composite plot for Seurat objects
#'
#' Note: This function requires the \pkg{scCustomize} package at runtime.
#' If it is not installed, an informative error is raised with installation instructions.
#' For scCustomize see: https://github.com/samuel-marsh/scCustomize
#' Marsh SE (2021). scCustomize: Custom Visualizations & Functions for Streamlined Analyses
#' of Single Cell Sequencing. https://doi.org/10.5281/zenodo.5706430. RRID:SCR_024675.
#'
#' @param object Seurat object.
#' @param features character vector of features.
#' @param label logical; pass to FeaturePlot_scCustom for labeling.
#' @param repel logical; repel labels in FeaturePlot_scCustom.
#' @param reduction character; reduction for FeaturePlot_scCustom, e.g. "umap.wnn".
#' @param pt.size numeric; point size for FeaturePlot_scCustom.
#' @param assay character; assay passed to DotPlot_scCustom. If NULL, uses DefaultAssay(object).
#' @param group.by character or factor; metadata/ident field used for DotPlot grouping.
#'   If NULL, defaults to Seurat::Idents(object).
#' @param split.by character; optional metadata to split DotPlot.
#' @param legend logical; show legend on FeaturePlots.
#' @param umap_axes logical; show axes on FeaturePlots.
#' @param bold_umap_labels logical; try to bold UMAP text labels (best effort).
#' @param featureplot.label.size numeric; label size for FeaturePlot_scCustom.
#' @param dotplot.label.size numeric; y-axis label size in the DotPlot panel.
#' @param bold.y.axis logical; bolden y-axis labels in the DotPlot.
#' @param colors optional character vector of colors for FeaturePlot/DotPlot gradients;
#'   if NULL, a pink palette is used. Provide >= 3 colors for smoother gradients.
#'
#' @return A patchwork/ggplot object.
#' @export
#'
#' @importFrom ggplot2 theme margin element_text
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat Idents NoAxes NoLegend DefaultAssay
#'
#' @examples
#' \dontrun{
#'   p <- FeatureDotPlot(
#'     object = hypo,
#'     features = c('myrf','olig1','gpr17','spi1'),
#'     label = TRUE,
#'     repel = TRUE,
#'     reduction = 'umap.wnn',
#'     pt.size = 0.3,
#'     legend = FALSE,
#'     umap_axes = FALSE,
#'     assay = 'SCT',
#'     group.by = 'subclusters_2.1.1_new',
#'     bold_umap_labels = TRUE,
#'     dotplot.label.size = 8,
#'     featureplot.label.size = 3,
#'     bold.y.axis = TRUE
#'   )
#'   p
#' }
FeatureDotPlot <- function(
    object,
    features,
    label = TRUE,
    repel = FALSE,
    reduction = "umap.wnn",
    pt.size = 0.3,
    assay = NULL,
    group.by = NULL,
    split.by = NULL,
    legend = TRUE,
    umap_axes = TRUE,
    bold_umap_labels = TRUE,
    featureplot.label.size = 3,
    dotplot.label.size = 10,
    bold.y.axis = FALSE,
    colors = NULL
) {
  # ---- dependency guard (Suggests) ----
  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    stop(
      "Package 'scCustomize' is required for FeatureDotPlot().\n",
      "Install with: remotes::install_github('samuel-marsh/scCustomize')",
      call. = FALSE
    )
  }

  suppressWarnings({
    # ---- resolve defaults that depend on 'object' ----
    if (is.null(assay))    assay    <- Seurat::DefaultAssay(object)
    if (is.null(group.by)) group.by <- Seurat::Idents(object)

    # ---- palettes ----
    pink_palette <- c("gray87", "pink", "palevioletred1", "deeppink", "deeppink3", "deeppink4")
    if (!is.null(colors) && length(colors) < 3) {
      # swallow warnings globally per your preference; keeping message logic here is harmless
      # (it will be suppressed by suppressWarnings())
      warning("`colors` has < 3 entries; consider providing >=3 colors for smoother gradients.")
    }
    if (is.null(colors)) colors <- pink_palette

    # ---- feature cap ----
    feat_vec <- as.character(features)
    if (length(feat_vec) > 4) feat_vec <- feat_vec[1:4]

    # ---- UMAP feature plot(s) ----
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

      if (!umap_axes) p <- p + Seurat::NoAxes()
      if (!legend)    p <- p + Seurat::NoLegend()

      if (bold_umap_labels && isTRUE(label)) {
        # try to bold a text-geom layer without assuming fixed layer index
        idx <- which(vapply(p$layers, function(ly) inherits(ly$geom, "GeomText"), logical(1)))
        if (length(idx) >= 1) {
          p$layers[[idx[1]]]$aes_params$fontface <- "bold"
        }
      }
      p
    })

    # ---- combine feature plots ----
    a <- if (length(feat_list) == 1) feat_list[[1]] else patchwork::wrap_plots(feat_list, ncol = 2)

    face <- if (isTRUE(bold.y.axis)) "bold" else "plain"

    # ---- Dot plot ----
    b <- scCustomize::DotPlot_scCustom(
      object      = object,
      features    = feat_vec,
      dot.scale   = 5.5,
      col.min     = 0.05,
      dot.min     = 0.01,
      scale       = TRUE,
      colors_use  = colors,
      assay       = assay,
      group.by    = group.by,
      split.by    = split.by
    ) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0, r = 50, b = 20, l = 50, unit = "pt"),
        axis.text.y = ggplot2::element_text(size = dotplot.label.size, face = face)
      )

    # ---- final layout ----
    patchwork::wrap_plots(a, b, ncol = 2, widths = c(4, 1))
  })
}
