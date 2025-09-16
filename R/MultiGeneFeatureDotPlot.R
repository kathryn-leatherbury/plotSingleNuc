#' Multi-gene UMAP(s) + DotPlot composite for Seurat objects
#'
#' Builds either (a) a single composite UMAP with multi-gene color overlays
#' or (b) a grid of one UMAP per gene (controlled by `split.umaps`), and
#' appends a DotPlot summarizing expression across groups. Gene names can be
#' matched to assay features directly or via a homolog mapping table.
#'
#' @param object Seurat object.
#' @param reduction Character. Embedding to use (e.g., "umap", "umap.wnn").
#' @param genes Character vector of genes (mzebra or human; mixed allowed).
#' @param na.cutoff Numeric. Values ≤ this are treated as “off” for plotting.
#' @param max.expression Optional numeric upper cap for color scales.
#' @param scale.legend Logical. If TRUE, force all genes to share the same limits.
#' @param scale.legend.max Logical. With scale.legend=TRUE, use max of per-gene 95th pct.
#' @param pt.size Numeric point size for UMAP.
#' @param alpha Numeric point alpha for UMAP.
#' @param group.by Column in meta.data (or vector/factor) for labels/DotPlot grouping.
#' @param label Logical. Draw cluster labels on UMAP(s).
#' @param repel Logical. Use ggrepel for labels (if available).
#' @param gene.ref.table Optional data.frame for homolog mapping.
#' @param gene.ref.col Column in gene.ref.table with mzebra IDs.
#' @param new.name.col Column with human homolog names.
#' @param secondary.ref.col Logical column indicating true homolog rows.
#' @param convert.names Logical. If TRUE, show human homologs in legends/x-axis.
#' @param assay Assay to use. Defaults to Seurat::DefaultAssay(object).
#' @param legend.ncols Optional integer. If set (and split.umaps=FALSE),
#'   renders a separate multi-column legend panel next to the UMAP.
#' @param plot.title Optional title/subtitle for UMAP.
#' @param split.by Optional column name for DotPlot splitting.
#' @param split.umaps Logical. If TRUE, draw one UMAP per gene.
#' @param umap.ncols Integer. Columns for per-gene UMAP grid (when split.umaps=TRUE).
#' @param flip.axes Logical; if TRUE, DotPlot rotated & stacked under UMAP(s).
#' @param plot.heights Numeric length-2 for stacked layout (UMAPs, DotPlot).
#' @param plot.widths Numeric length-2 for side-by-side layout (UMAPs, DotPlot).
#'
#' @return A patchwork object combining UMAP (single or grid) and a DotPlot.
#'
#' @examples
#' \dontrun{
#' MultiGeneFeatureDotPlot(
#'   object = seurat_obj,
#'   reduction = "umap.wnn",
#'   genes = c("fat2","LOC101479120","cga","fabp7","myrf","pitx2"),
#'   gene.ref.table = homologs, convert.names = TRUE,
#'   label = TRUE, repel = TRUE,
#'   split.umaps = TRUE, umap.ncols = 3
#' )
#' }
#'
#' @importFrom stats setNames
#' @importFrom rlang .data
#' @export
MultiGeneFeatureDotPlot <- function(
    object,
    reduction          = "umap",
    genes              = NULL,          # mz or human symbols (mixed ok)
    na.cutoff          = 0.5,
    max.expression     = NULL,
    scale.legend       = FALSE,
    scale.legend.max   = FALSE,
    pt.size            = 0.5,
    alpha              = 0.5,
    group.by           = NULL,          # default -> Idents(object)
    label              = TRUE,
    repel              = FALSE,
    gene.ref.table     = NULL,
    gene.ref.col       = "mzebra_seurat_gene_id",
    new.name.col       = "human_homolog",
    secondary.ref.col  = "mzebra_to_human_true_homolog",
    convert.names      = FALSE,
    assay              = Seurat::DefaultAssay(object),
    legend.ncols       = NULL,
    plot.title         = NULL,
    split.by           = NULL,
    split.umaps        = FALSE,
    flip.axes          = FALSE,
    umap.ncols         = 3,
    plot.heights       = NULL,
    plot.widths        = NULL
){
  .build <- function(){
    
    ## --- helpers --------------------------------------------------------------
    stop_if <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)
    to_bool <- function(x){
      if (is.logical(x)) return(x)
      tolower(as.character(x)) %in% c("true","t","1","yes","y")
    }
    need_pkgs <- function(pkgs){
      miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
      if (length(miss)) stop("Please install: {", paste(miss, collapse = ", "), "}", call. = FALSE)
    }
    first_or_na <- function(x) if (length(x)) x[1] else NA_character_
    cap_top <- function(x, p = 0.99){
      nz <- x[is.finite(x) & x > 0]
      if (!length(nz)) return(x)
      pmin(x, stats::quantile(nz, p, na.rm = TRUE))
    }
    finite_pos_q <- function(x, p = 0.99){
      x <- x[is.finite(x) & x > 0]
      if (length(x)) stats::quantile(x, p, na.rm = TRUE) else 0
    }
    finite_pos_max <- function(x){
      x <- x[is.finite(x) & x > 0]
      if (length(x)) max(x) else 0
    }
    three_ticks_strict <- function(lo, hi){
      if (!is.finite(hi)) hi <- 1
      if (hi <= lo) hi <- lo + 1e-8
      c(lo, (lo + hi)/2, hi)
    }
    lab_num <- scales::label_number(accuracy = 0.1, trim = TRUE)
    
    ## --- guards ---------------------------------------------------------------
    stop_if(!inherits(object, "Seurat"), "`object` must be a Seurat object.")
    base_needs <- c("ggplot2","ggnewscale","scales","grid","patchwork","scCustomize","cowplot","ggplotify")
    need_pkgs(base_needs)
    if (isTRUE(repel)) need_pkgs("ggrepel")
    
    stop_if(is.null(genes) || length(genes) < 1, "`genes` must be a non-empty character vector.")
    stop_if(!assay %in% names(object@assays),
            paste0("Assay '", assay, "' not found. Available: ",
                   paste(names(object@assays), collapse = ", ")))
    stop_if(!reduction %in% Seurat::Reductions(object),
            paste0("Reduction '", reduction, "' not found. Available: ",
                   paste(Seurat::Reductions(object), collapse = ", ")))
    
    emb2 <- Seurat::Embeddings(object, reduction = reduction)[, 1:2, drop = FALSE]
    stop_if(!is.matrix(emb2) || ncol(emb2) < 2, "Reduction must have at least 2 dimensions.")
    xlab <- colnames(emb2)[1]; ylab <- colnames(emb2)[2]
    
    ## --- resolve requested genes -> assay features ----------------------------
    feats <- rownames(object[[assay]])
    if (!is.null(gene.ref.table) && !is.null(gene.ref.col) && !is.null(new.name.col)) {
      stop_if(!all(c(gene.ref.col, new.name.col) %in% names(gene.ref.table)),
              "`gene.ref.col` and/or `new.name.col` not found in `gene.ref.table`.")
      stop_if(!secondary.ref.col %in% names(gene.ref.table),
              paste0("`secondary.ref.col` ('", secondary.ref.col, "') not found in `gene.ref.table`."))
      
      g_vec <- as.character(gene.ref.table[[gene.ref.col]])
      h_vec <- as.character(gene.ref.table[[new.name.col]])
      g_lc  <- tolower(g_vec)
      h_lc  <- tolower(h_vec)
      tmask <- to_bool(gene.ref.table[[secondary.ref.col]])
      
      map_one <- function(sym){
        s <- tolower(sym)
        hit_g <- which(g_lc == s)
        if (length(hit_g) >= 1) return(g_vec[first_or_na(hit_g)])
        hit_h <- which(h_lc == s)
        if (length(hit_h) >= 1) {
          hit_true <- hit_h[tmask[hit_h]]
          if (length(hit_true) >= 1) return(g_vec[first_or_na(hit_true)])
          return(g_vec[first_or_na(hit_h)])
        }
        return(NA_character_)
      }
      mapped <- vapply(genes, map_one, character(1))
      stop_if(anyNA(mapped),
              paste0("Not found in ref table: ",
                     paste(genes[is.na(mapped)], collapse = ", ")))
      genes_resolved <- mapped
    } else {
      genes_resolved <- genes
    }
    
    resolve_in_assay <- function(g){
      if (g %in% feats) return(g)
      hit <- feats[toupper(feats) == toupper(g)]
      if (length(hit) == 1) return(hit)
      first_or_na(grep(paste0("^", g, "$"), feats, ignore.case = TRUE, value = TRUE))
    }
    resolved <- vapply(genes_resolved, resolve_in_assay, character(1))
    stop_if(anyNA(resolved),
            paste0("Not found in assay '", assay, "': ",
                   paste(genes_resolved[is.na(resolved)], collapse = ", ")))
    genes_resolved <- as.character(genes_resolved)
    
    ## --- fetch data -----------------------------------------------------------
    old_assay <- Seurat::DefaultAssay(object)
    on.exit(Seurat::DefaultAssay(object) <- old_assay, add = TRUE)
    Seurat::DefaultAssay(object) <- assay
    
    expr <- Seurat::FetchData(object, vars = resolved, slot = "data")
    colnames(expr) <- genes_resolved
    df <- data.frame(as.data.frame(emb2), expr, check.names = FALSE)
    
    ## --- constants & scaling --------------------------------------------------
    legend_title_size <- 11
    legend_label_size <- 10
    pt.size_bckgr     <- pt.size
    lower.limit       <- if (na.cutoff > 0) na.cutoff else 0
    
    for (g in genes_resolved) df[[g]] <- cap_top(df[[g]])
    
    per_gene_q99 <- vapply(genes_resolved, function(g) finite_pos_q(df[[g]], 0.99), numeric(1))
    per_gene_q95 <- vapply(genes_resolved, function(g) finite_pos_q(df[[g]], 0.95), numeric(1))
    per_gene_max <- vapply(genes_resolved, function(g) finite_pos_max(df[[g]]), numeric(1))
    
    if (isTRUE(scale.legend)) {
      shared_upper <- if (!is.null(max.expression)) {
        max.expression
      } else if (isTRUE(scale.legend.max)) {
        max(per_gene_q95)
      } else {
        min(per_gene_q99)
      }
      lo_use <- lower.limit
      if (shared_upper <= lo_use) shared_upper <- lo_use + 1e-8
      limits_for_scale <- c(lo_use, shared_upper)
      brks_common <- three_ticks_strict(lo_use, shared_upper)
      lims_list <- setNames(replicate(length(genes_resolved), limits_for_scale, simplify = FALSE), genes_resolved)
      brks_list <- setNames(replicate(length(genes_resolved), brks_common,  simplify = FALSE), genes_resolved)
    } else {
      lims_list <- lapply(seq_along(genes_resolved), function(i){
        lo <- lower.limit
        hi <- if (is.null(max.expression)) per_gene_max[i] else max.expression
        if (hi <= lo) hi <- lo + 1e-8
        c(lo, hi)
      })
      names(lims_list) <- genes_resolved
      brks_list <- lapply(lims_list, function(lg) three_ticks_strict(lg[1], lg[2]))
      names(brks_list) <- genes_resolved
    }
    
    ## --- palettes & legend labels --------------------------------------------
    pals <- list(
      c("gray80","lightskyblue","dodgerblue3","dodgerblue4"),
      c("gray80","plum2","purple2","purple4"),
      c("gray80","lightgreen","forestgreen","darkgreen"),
      c("gray80","orange","indianred3","firebrick"),
      c("gray80","pink1","deeppink2","deeppink3"),
      c("gray80","darkslategray2","lightseagreen","turquoise4")
    )
    na_color <- "gray80"
    
    map_one_homolog <- function(g){
      if (is.null(gene.ref.table) || is.null(gene.ref.col) || is.null(new.name.col)) return(NA_character_)
      if (!all(c(gene.ref.col, new.name.col) %in% names(gene.ref.table))) return(NA_character_)
      gene_vec <- as.character(gene.ref.table[[gene.ref.col]])
      hom_vec  <- as.character(gene.ref.table[[new.name.col]])
      hit <- which(gene_vec == g)
      if (length(hit)) hom_vec[hit[1]] else NA_character_
    }
    legend_labels <- vapply(genes_resolved, function(g){
      lbl <- g
      if (isTRUE(convert.names)) {
        hm <- map_one_homolog(g); if (is.character(hm) && nzchar(hm)) lbl <- hm
      }
      toupper(as.character(lbl))
    }, character(1))
    
    ## --- base UMAP background -------------------------------------------------
    p_base <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = df, ggplot2::aes_string(x = xlab, y = ylab),
        color = na_color, size = pt.size_bckgr
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        legend.position  = if (is.null(legend.ncols)) "right" else "none",
        legend.direction = "vertical",
        legend.box       = "vertical",
        plot.margin      = grid::unit(c(4,4,4,4), "pt")
      ) +
      Seurat::NoAxes()
    
    ## --- cluster labels once --------------------------------------------------
    labs <- NULL
    if (isTRUE(label)) {
      if (is.null(group.by)) {
        labs <- as.character(Seurat::Idents(object))
      } else if (is.character(group.by) && length(group.by) == 1) {
        stop_if(!group.by %in% colnames(object@meta.data),
                paste0("Column `", group.by, "` not found in object@meta.data."))
        labs <- as.character(object@meta.data[[group.by]])
      } else if (is.vector(group.by) || is.factor(group.by)) {
        stop_if(length(group.by) != ncol(object),
                "`group.by` vector length must equal number of cells in the object.")
        labs <- as.character(group.by)
      } else {
        stop("`group.by` must be NULL, a single column name, or a vector/factor of length ncol(object).",
             call. = FALSE)
      }
    }
    tmp  <- if (!is.null(labs)) data.frame(emb2, label = labs, check.names = FALSE) else NULL
    medoids <- if (!is.null(tmp)) do.call(rbind, lapply(split(tmp, tmp$label), function(d){
      cx <- mean(d[[xlab]], na.rm = TRUE); cy <- mean(d[[ylab]], na.rm = TRUE)
      i  <- which.min((d[[xlab]] - cx)^2 + (d[[ylab]] - cy)^2)
      data.frame(label = d$label[i], x = d[[xlab]][i], y = d[[ylab]][i])
    })) else NULL
    
    ## ==================== BUILD UMAP(S) =======================================
    if (isTRUE(split.umaps)) {
      # per-gene panel list
      umap_list <- lapply(seq_along(genes_resolved), function(i){
        g   <- genes_resolved[i]
        pal <- pals[[ ((i-1) %% length(pals)) + 1 ]]
        lo  <- lims_list[[g]][1]; hi <- lims_list[[g]][2]
        pos_seq <- seq(lo, hi, length.out = length(pal))
        vals    <- scales::rescale(pos_seq, to = c(0,1), from = c(lo, hi))
        p_i <- p_base +
          ggplot2::geom_point(
            data = df[df[[g]] > lower.limit & is.finite(df[[g]]), , drop = FALSE],
            ggplot2::aes_string(x = xlab, y = ylab, color = g),
            size = pt.size, alpha = alpha
          ) +
          ggplot2::scale_color_gradientn(
            name    = legend_labels[i],
            colours = pal,
            values  = vals,
            limits  = c(lo, hi),
            breaks  = brks_list[[g]],
            labels  = lab_num,
            oob     = scales::squish,
            na.value = na_color,
            guide   = ggplot2::guide_colourbar(
              nbin = 100,
              barheight = grid::unit(75, "pt"),
              barwidth  = grid::unit(25, "pt"),
              title.theme = ggplot2::element_text(face = "bold", size = legend_title_size),
              label.theme = ggplot2::element_text(size = legend_label_size)
            )
          ) +
          ggplot2::labs(title = legend_labels[i])
        
        if (!is.null(medoids)) {
          if (isTRUE(repel) && requireNamespace("ggrepel", quietly = TRUE)) {
            p_i <- p_i + ggrepel::geom_text_repel(
              data = medoids,
              ggplot2::aes(x = x, y = y, label = label),
              size = 3.8, fontface = "bold", color = "black",
              box.padding = 0.5, point.padding = 0.3,
              segment.color = "grey50", max.overlaps = Inf
            )
          } else {
            p_i <- p_i + ggplot2::geom_text(
              data = medoids,
              ggplot2::aes(x = x, y = y, label = label),
              size = 3.8, fontface = "bold", color = "black"
            )
          }
        }
        p_i
      })
      
      # turn the grid into a *ggplot* so ggplot2::labs() will work:
      p_grid <- cowplot::plot_grid(plotlist = umap_list, ncol = umap.ncols, align = "hv")
      p_umap <- ggplotify::as.ggplot(p_grid)
      
      # attach title/subtitle the SAME way as your original code
      if (!is.null(plot.title)) {
        p_umap <- p_umap +
          ggplot2::labs(title = plot.title) +
          ggplot2::theme(plot.margin = ggplot2::margin(0, 10, 0, 10, unit = "pt"),
                         plot.title    = ggplot2::element_text(face = "bold", size = 15, hjust = 0, vjust = -2))
      }
      
    } else {
      # original single composite UMAP (overlay)
      order_across_genes <- TRUE
      n_order_slices     <- 5
      df_list <- lapply(genes_resolved, function(g){
        d <- df[df[[g]] > lower.limit & is.finite(df[[g]]), , drop = FALSE]
        if (!nrow(d)) return(d)
        r <- rank(d[[g]], ties.method = "average")
        d$.__pct__ <- r / max(r)
        d
      })
      names(df_list) <- genes_resolved
      
      p <- p_base
      if (!order_across_genes) {
        for (i in seq_along(genes_resolved)) {
          if (i > 1) p <- p + ggnewscale::new_scale_color()
          g   <- genes_resolved[i]
          dat <- df_list[[g]]
          if (!nrow(dat)) next
          pal <- pals[[ ((i-1) %% length(pals)) + 1 ]]
          lo <- lims_list[[g]][1]; hi <- lims_list[[g]][2]
          pos_seq <- seq(lo, hi, length.out = length(pal))
          vals    <- scales::rescale(pos_seq, to = c(0,1), from = c(lo, hi))
          p <- p +
            ggplot2::geom_point(
              data = dat, ggplot2::aes_string(x = xlab, y = ylab, color = g),
              size = pt.size, alpha = alpha
            ) +
            ggplot2::scale_color_gradientn(
              name    = legend_labels[i],
              colours = pal,
              values  = vals,
              limits  = c(lo, hi),
              breaks  = brks_list[[g]],
              labels  = lab_num,
              oob     = scales::squish,
              na.value = na_color,
              guide   = ggplot2::guide_colourbar(
                nbin = 100,
                barheight = grid::unit(75, "pt"),
                barwidth  = grid::unit(25, "pt"),
                title.theme = ggplot2::element_text(face = "bold", size = legend_title_size),
                label.theme = ggplot2::element_text(size = legend_label_size)
              )
            )
        }
      } else {
        cuts <- seq(0, 1, length.out = max(2, n_order_slices) + 1)
        for (s in seq_len(length(cuts) - 1)) {
          for (i in seq_along(genes_resolved)) {
            g   <- genes_resolved[i]
            dat <- df_list[[g]]
            if (!nrow(dat)) next
            slice <- dat[dat$.__pct__ > cuts[s] & dat$.__pct__ <= cuts[s+1], , drop = FALSE]
            if (!nrow(slice)) next
            if (s > 1 || i > 1) p <- p + ggnewscale::new_scale_color()
            pal <- pals[[ ((i-1) %% length(pals)) + 1 ]]
            lo <- lims_list[[g]][1]; hi <- lims_list[[g]][2]
            legend_on <- (s == (length(cuts) - 1))
            pos_seq <- seq(lo, hi, length.out = length(pal))
            vals    <- scales::rescale(pos_seq, to = c(0,1), from = c(lo, hi))
            p <- p +
              ggplot2::geom_point(
                data = slice, ggplot2::aes_string(x = xlab, y = ylab, color = g),
                size = pt.size, alpha = alpha, show.legend = legend_on
              ) +
              ggplot2::scale_color_gradientn(
                name    = legend_labels[i],
                colours = pal,
                values  = vals,
                limits  = c(lo, hi),
                breaks  = brks_list[[g]],
                labels  = lab_num,
                oob     = scales::squish,
                na.value = na_color,
                guide   = if (legend_on) ggplot2::guide_colourbar(
                  nbin = 100,
                  barheight = grid::unit(75, "pt"),
                  barwidth  = grid::unit(25, "pt"),
                  title.theme = ggplot2::element_text(face = "bold", size = legend_title_size),
                  label.theme = ggplot2::element_text(size = legend_label_size)
                ) else "none"
              )
          }
        }
      }
      
      if (!is.null(legend.ncols)) {
        build_leg <- function(i){
          g   <- genes_resolved[i]
          pal <- pals[[ ((i-1) %% length(pals)) + 1 ]]
          lo <- lims_list[[g]][1]; hi <- lims_list[[g]][2]
          pos_seq <- seq(lo, hi, length.out = length(pal))
          vals    <- scales::rescale(pos_seq, to = c(0,1), from = c(lo, hi))
          dummy <- ggplot2::ggplot(data.frame(x=0,y=0,z=0), ggplot2::aes(x,y,color=z)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_gradientn(
              name    = legend_labels[i],
              colours = pal,
              values  = vals,
              limits  = c(lo, hi),
              breaks  = three_ticks_strict(lo, hi),
              labels  = lab_num,
              oob     = scales::squish,
              guide   = ggplot2::guide_colourbar(
                nbin = 100,
                barheight = grid::unit(75, "pt"),
                barwidth  = grid::unit(25, "pt"),
                title.theme = ggplot2::element_text(face = "bold", size = legend_title_size),
                label.theme = ggplot2::element_text(size = legend_label_size)
              )
            ) +
            ggplot2::theme_void() +
            ggplot2::theme(
              legend.position    = "right",
              legend.margin      = ggplot2::margin(0,0,0,0),
              legend.box.margin  = ggplot2::margin(0,0,0,0),
              legend.box.spacing = grid::unit(0, "pt"),
              legend.spacing.y   = ggplot2::unit(0, "pt"),
              plot.margin        = grid::unit(c(0,0,0,0), "pt")
            )
          cowplot::get_legend(dummy)
        }
        legs <- lapply(seq_along(genes_resolved), build_leg)
        tiles <- lapply(legs, function(g) cowplot::ggdraw(g) +
                          cowplot::theme_nothing() +
                          ggplot2::theme(plot.margin = grid::unit(c(0, 2, 0, 0), "pt")))
        nrow_legend  <- ceiling(length(tiles) / legend.ncols)
        legend_panel <- cowplot::plot_grid(plotlist = tiles,
                                           ncol = legend.ncols, nrow = nrow_legend,
                                           align = "hv", axis = "tblr", greedy = TRUE)
        p <- ggplotify::as.ggplot(cowplot::plot_grid(p, legend_panel, ncol = 2, rel_widths = c(1, 0.28), align = "h"))
      }
      
      if (!is.null(plot.title)) {
        p <- p + ggplot2::labs(title = plot.title) +
          ggplot2::theme(plot.margin = ggplot2::margin(0, 10, 0, 10, unit = "pt"),
                         plot.title    = ggplot2::element_text(face = "bold", size = 15, hjust = 0, vjust = -2))
      }
      p_umap <- p
    }
    
    ## ==================== DotPlot =============================================
    features_dp <- genes_resolved
    lab_map     <- stats::setNames(legend_labels, genes_resolved)
    pink_palette <- c("gray87","pink","palevioletred1","deeppink","deeppink3","deeppink4")
    
    b <- scCustomize::DotPlot_scCustom(
      seurat_object = object,
      features      = features_dp,
      dot.scale     = 5.5,
      col.min       = 0.05,
      dot.min       = 0.01,
      scale         = TRUE,
      colors_use    = pink_palette,
      assay         = assay,
      group.by      = group.by,
      split.by      = split.by,
      flip_axes     = flip.axes
    ) +
      ggplot2::scale_x_discrete(limits = features_dp,
                                labels = unname(lab_map[features_dp])) +
      ggplot2::labs(x = NULL) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 15, l = 10, unit = "pt"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.text.y = ggplot2::element_text(size = 9)
      )
    
    ## ==================== Final layout ========================================
    if (isTRUE(flip.axes)) {
      b <- b + RotatedAxis()
      if (is.null(plot.heights)) plot.heights <- c(4, 1)
      c <- patchwork::wrap_plots(p_umap, b, ncol = 1, heights = plot.heights)
    } else {
      if (is.null(plot.widths)) plot.widths <- c(4, 1)
      c <- patchwork::wrap_plots(p_umap, b, ncol = 2, widths = plot.widths)
    }
    
    return(c)
  }
  
  suppressMessages(suppressWarnings(.build()))
}