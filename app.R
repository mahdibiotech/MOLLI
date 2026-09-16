# =========================================================
# MOLLI APP
# Multi-Omics Liver Laboratory Interactive App
#
# Modules :
#   1. Résumé
#   2. Interactome
#   3. Régions CRM / TSS
#   4. Exploration globale des gènes
#   5. Transcriptome
#   6. Intégration Interactome × Transcriptome
#   7. Dataset / QC
#   8. Aide
#
# Dépendances :
# shiny, data.table, dplyr, pheatmap, DT, ggplot2,
# colourpicker, shinyWidgets, scales
# =========================================================


# =========================================================
# Packages
# =========================================================

library(shiny)
library(data.table)
library(dplyr)
library(pheatmap)
library(DT)
library(ggplot2)
library(colourpicker)
library(shinyWidgets)
library(scales)


# =========================================================
# Configuration
# =========================================================

BINDING_FILE <- "./pool_MultiBind_ProseqActiveLiverPromoters_Annotated_withColnames.bed"
RNASEQ_FILE  <- "./RNAseq_DE_results.tsv"


# =========================================================
# HELPERS - Lecture et préparation des données
# =========================================================

read_binding_table <- function(path) {

  if (!file.exists(path)) {
    stop("Fichier interactome introuvable : ", path)
  }

  dt <- fread(path)

  required_cols <- c(
    "chrCRM", "startCRM", "endCRM",
    "tssC", "tssS", "tssE",
    "Name", "inter"
  )

  missing_cols <- setdiff(required_cols, colnames(dt))

  if (length(missing_cols) > 0) {
    stop(
      "Colonnes manquantes dans le fichier interactome : ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt[, Name := trimws(as.character(Name))]
  dt[, chrCRM := as.character(chrCRM)]
  dt[, tssC := as.character(tssC)]
  dt[, startCRM := as.numeric(startCRM)]
  dt[, endCRM := as.numeric(endCRM)]
  dt[, tssS := as.numeric(tssS)]
  dt[, tssE := as.numeric(tssE)]

  protein_cols <- setdiff(colnames(dt), required_cols)

  if (length(protein_cols) == 0) {
    stop("Aucune colonne protéique détectée.")
  }

  for (protein in protein_cols) {
    dt[, (protein) := as.numeric(get(protein))]
    dt[is.na(get(protein)), (protein) := 0]
    dt[, (protein) := as.integer(get(protein) > 0)]
  }

  dt <- dt[
    !is.na(Name) &
      Name != ""
  ]

  list(
    data = dt,
    protein_cols = protein_cols
  )
}


read_rnaseq_table <- function(path) {

  if (!file.exists(path)) {
    stop("Fichier RNA-seq introuvable : ", path)
  }

  dt <- fread(path)

  required_cols <- c(
    "gene",
    "StudiesCond",
    "abslogFC",
    "FDR",
    "State"
  )

  missing_cols <- setdiff(required_cols, colnames(dt))

  if (length(missing_cols) > 0) {
    stop(
      "Colonnes RNA-seq manquantes : ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt[, gene := trimws(as.character(gene))]
  dt[, StudiesCond := trimws(as.character(StudiesCond))]
  dt[, abslogFC := as.numeric(abslogFC)]
  dt[, FDR := as.numeric(FDR)]
  dt[, State := tolower(trimws(as.character(State)))]

  dt <- dt[
    !is.na(gene) &
      gene != "" &
      !is.na(StudiesCond) &
      StudiesCond != "" &
      !is.na(abslogFC)
  ]

  dt <- dt[
    State %in% c("up", "down")
  ]

  dt[
    ,
    signed_logFC := fifelse(
      State == "down",
      -abslogFC,
      abslogFC
    )
  ]

  dt
}


# =========================================================
# HELPERS - Interactome
# =========================================================

get_palette <- function(palette_name) {

  switch(
    palette_name,
    "Rouge" = c("white", "red"),
    "Bleu" = c("white", "royalblue"),
    "Vert" = c("white", "forestgreen"),
    "Noir/Blanc" = c("white", "black"),
    "Violet" = c("white", "purple"),
    c("white", "red")
  )
}


aggregate_to_gene_matrix <- function(
    dt,
    genes,
    protein_cols,
    keep_only_bound = TRUE
) {

  sub <- dt[
    Name %in% genes
  ]

  if (nrow(sub) == 0) {
    return(NULL)
  }

  agg <- sub[
    ,
    lapply(
      .SD,
      function(x) {
        as.integer(
          any(
            x == 1,
            na.rm = TRUE
          )
        )
      }
    ),
    by = Name,
    .SDcols = protein_cols
  ]

  mat <- as.matrix(
    agg[
      ,
      ..protein_cols
    ]
  )

  rownames(mat) <- agg$Name
  storage.mode(mat) <- "numeric"

  ordered_genes <- genes[
    genes %in% rownames(mat)
  ]

  mat <- mat[
    ordered_genes,
    ,
    drop = FALSE
  ]

  if (keep_only_bound) {

    keep_cols <- colSums(
      mat,
      na.rm = TRUE
    ) > 0

    mat <- mat[
      ,
      keep_cols,
      drop = FALSE
    ]
  }

  if (ncol(mat) == 0) {
    return(NULL)
  }

  mat
}


make_heatmap <- function(
    mat,
    palette_name = "Rouge",
    fontsize_row = 12,
    fontsize_col = 7,
    main = "Interactions protéines-régions régulatrices"
) {

  if (is.null(mat)) {
    plot.new()
    title("Aucune donnée disponible")
    return(invisible(NULL))
  }

  cluster_rows_value <- nrow(mat) > 1
  cluster_cols_value <- ncol(mat) > 1

  pheatmap(
    mat,
    cluster_rows = cluster_rows_value,
    cluster_cols = cluster_cols_value,
    color = get_palette(palette_name),
    breaks = c(-0.01, 0.5, 1.01),
    border_color = "grey60",
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    legend_breaks = c(0, 1),
    legend_labels = c("Absence", "Présence"),
    main = main
  )
}


make_crm_protein_table <- function(
    dt,
    protein_cols
) {

  crm_dt <- copy(dt)

  crm_dt[
    ,
    CRM_ID := paste0(
      chrCRM,
      ":",
      startCRM,
      "-",
      endCRM
    )
  ]

  crm_dt[
    ,
    TSS_ID := paste0(
      tssC,
      ":",
      tssS
    )
  ]

  crm_dt[
    ,
    CRM_size := endCRM - startCRM
  ]

  crm_long <- melt(
    crm_dt,
    id.vars = c(
      "Name",
      "CRM_ID",
      "chrCRM",
      "startCRM",
      "endCRM",
      "CRM_size",
      "TSS_ID",
      "tssC",
      "tssS",
      "tssE",
      "inter"
    ),
    measure.vars = protein_cols,
    variable.name = "Protein",
    value.name = "Binding"
  )

  crm_long <- crm_long[
    Binding == 1
  ]

  unique(crm_long)
}


compute_jaccard_matrix <- function(
    dt,
    proteins
) {

  proteins <- proteins[
    proteins %in% colnames(dt)
  ]

  if (length(proteins) == 0) {
    return(NULL)
  }

  mat <- as.matrix(
    dt[
      ,
      ..proteins
    ]
  )

  storage.mode(mat) <- "numeric"

  jaccard <- matrix(
    0,
    nrow = length(proteins),
    ncol = length(proteins),
    dimnames = list(
      proteins,
      proteins
    )
  )

  for (i in seq_along(proteins)) {
    for (j in seq_along(proteins)) {

      a <- mat[, i] == 1
      b <- mat[, j] == 1

      union_n <- sum(
        a | b,
        na.rm = TRUE
      )

      inter_n <- sum(
        a & b,
        na.rm = TRUE
      )

      jaccard[i, j] <- if (
        union_n == 0
      ) {
        NA_real_
      } else {
        inter_n / union_n
      }
    }
  }

  jaccard
}


make_upset_like_data <- function(
    dt,
    proteins,
    top_n = 20
) {

  proteins <- proteins[
    proteins %in% colnames(dt)
  ]

  if (length(proteins) == 0) {
    return(data.table())
  }

  if (length(proteins) > 8) {
    proteins <- proteins[1:8]
  }

  tmp <- copy(
    dt[
      ,
      c(
        "Name",
        proteins
      ),
      with = FALSE
    ]
  )

  tmp[
    ,
    Combination := apply(
      .SD,
      1,
      function(x) {

        active <- proteins[
          as.numeric(x) == 1
        ]

        if (length(active) == 0) {
          "Aucune"
        } else {
          paste(
            active,
            collapse = " + "
          )
        }
      }
    ),
    .SDcols = proteins
  ]

  counts <- tmp[
    Combination != "Aucune",
    .N,
    by = Combination
  ]

  setorder(
    counts,
    -N
  )

  head(
    counts,
    top_n
  )
}


# =========================================================
# HELPERS - Transcriptome
# =========================================================

safe_neg_log10 <- function(x) {

  -log10(
    pmax(
      x,
      .Machine$double.xmin
    )
  )
}


# =========================================================
# HELPERS - Intégration
# =========================================================

build_binding_expression_table <- function(
    binding_dt,
    rnaseq_dt,
    protein,
    study
) {

  if (
    is.null(protein) ||
    is.null(study) ||
    length(protein) != 1 ||
    length(study) != 1
  ) {
    return(data.table())
  }

  if (!(protein %in% colnames(binding_dt))) {
    return(data.table())
  }

  bind <- binding_dt[
    ,
    .(
      gene = Name,
      Bound = get(protein),
      chrCRM,
      startCRM,
      endCRM,
      tssC,
      tssS,
      tssE
    )
  ]

  expr <- rnaseq_dt[
    StudiesCond == study,
    .(
      gene,
      StudiesCond,
      abslogFC,
      signed_logFC,
      FDR,
      State
    )
  ]

  merged <- merge(
    expr,
    bind,
    by = "gene",
    all = FALSE
  )

  merged[
    ,
    Binding_status := fifelse(
      Bound == 1,
      "Lié",
      "Non lié"
    )
  ]

  merged
}


compute_fisher_enrichment <- function(
    binding_dt,
    rnaseq_dt,
    protein,
    study,
    fdr_threshold = 0.05
) {

  dat <- build_binding_expression_table(
    binding_dt,
    rnaseq_dt,
    protein,
    study
  )

  if (nrow(dat) == 0) {
    return(
      data.table(
        Protein = protein,
        Study = study,
        OddsRatio = NA_real_,
        Pvalue = NA_real_,
        N_common = 0,
        Bound_DE = 0,
        Bound_nonDE = 0,
        Unbound_DE = 0,
        Unbound_nonDE = 0
      )
    )
  }

  dat[
    ,
    DE := !is.na(FDR) &
      FDR < fdr_threshold
  ]

  a <- sum(
    dat$Bound == 1 &
      dat$DE,
    na.rm = TRUE
  )

  b <- sum(
    dat$Bound == 1 &
      !dat$DE,
    na.rm = TRUE
  )

  c <- sum(
    dat$Bound == 0 &
      dat$DE,
    na.rm = TRUE
  )

  d <- sum(
    dat$Bound == 0 &
      !dat$DE,
    na.rm = TRUE
  )

  test_matrix <- matrix(
    c(
      a, b,
      c, d
    ),
    nrow = 2,
    byrow = TRUE
  )

  ft <- tryCatch(
    fisher.test(
      test_matrix
    ),
    error = function(e) NULL
  )

  odds <- if (is.null(ft)) {
    NA_real_
  } else {
    unname(ft$estimate)
  }

  pval <- if (is.null(ft)) {
    NA_real_
  } else {
    ft$p.value
  }

  data.table(
    Protein = protein,
    Study = study,
    OddsRatio = odds,
    Pvalue = pval,
    N_common = nrow(dat),
    Bound_DE = a,
    Bound_nonDE = b,
    Unbound_DE = c,
    Unbound_nonDE = d
  )
}


compute_enrichment_grid <- function(
    binding_dt,
    rnaseq_dt,
    proteins,
    studies,
    fdr_threshold
) {

  results <- vector(
    "list",
    length(proteins) * length(studies)
  )

  k <- 1

  for (study in studies) {
    for (protein in proteins) {

      results[[k]] <- compute_fisher_enrichment(
        binding_dt = binding_dt,
        rnaseq_dt = rnaseq_dt,
        protein = protein,
        study = study,
        fdr_threshold = fdr_threshold
      )

      k <- k + 1
    }
  }

  out <- rbindlist(
    results,
    fill = TRUE
  )

  out[
    ,
    log2OR := log2(
      OddsRatio
    )
  ]

  out[
    ,
    P_adj_BH := p.adjust(
      Pvalue,
      method = "BH"
    )
  ]

  out[
    ,
    negLog10P := safe_neg_log10(
      P_adj_BH
    )
  ]

  out
}


compute_directional_enrichment <- function(
    binding_dt, rnaseq_dt, protein, study, fdr_threshold = 0.05
) {
  dat <- build_binding_expression_table(binding_dt, rnaseq_dt, protein, study)
  if (nrow(dat) == 0) return(data.table())
  dat[, Significant := !is.na(FDR) & FDR < fdr_threshold]
  one_direction <- function(direction) {
    event <- dat$Significant & dat$State == direction
    tab <- matrix(c(
      sum(dat$Bound == 1 & event, na.rm=TRUE),
      sum(dat$Bound == 1 & !event, na.rm=TRUE),
      sum(dat$Bound == 0 & event, na.rm=TRUE),
      sum(dat$Bound == 0 & !event, na.rm=TRUE)
    ), nrow=2, byrow=TRUE)
    ft <- tryCatch(fisher.test(tab), error=function(e) NULL)
    data.table(Direction=direction, OddsRatio=if(is.null(ft)) NA_real_ else unname(ft$estimate), Pvalue=if(is.null(ft)) NA_real_ else ft$p.value)
  }
  ans <- rbindlist(list(one_direction("up"), one_direction("down")))
  ans[, P_adj_BH := p.adjust(Pvalue, method="BH")]
  ans
}

compute_binding_wilcoxon <- function(binding_dt, rnaseq_dt, protein, study) {
  dat <- build_binding_expression_table(binding_dt, rnaseq_dt, protein, study)
  if (nrow(dat)==0 || uniqueN(dat$Binding_status)<2) return(list(pvalue=NA_real_, median_bound=NA_real_, median_unbound=NA_real_))
  wt <- tryCatch(wilcox.test(signed_logFC ~ Binding_status, data=dat, exact=FALSE), error=function(e) NULL)
  list(
    pvalue=if(is.null(wt)) NA_real_ else wt$p.value,
    median_bound=median(dat[Binding_status=="Lié", signed_logFC], na.rm=TRUE),
    median_unbound=median(dat[Binding_status=="Non lié", signed_logFC], na.rm=TRUE)
  )
}


# =========================================================
# Chargement des données
# =========================================================

binding_obj <- read_binding_table(
  BINDING_FILE
)

binding_dt <- binding_obj$data
protein_cols <- binding_obj$protein_cols

rnaseq_dt <- read_rnaseq_table(
  RNASEQ_FILE
)

crm_protein_dt <- make_crm_protein_table(
  binding_dt,
  protein_cols
)


# =========================================================
# Listes globales et statistiques
# =========================================================

binding_genes <- sort(
  unique(
    binding_dt$Name
  )
)

expression_genes <- sort(
  unique(
    rnaseq_dt$gene
  )
)

all_genes <- sort(
  unique(
    c(
      binding_genes,
      expression_genes
    )
  )
)

all_studies <- sort(
  unique(
    rnaseq_dt$StudiesCond
  )
)

common_genes <- intersect(
  binding_genes,
  expression_genes
)

binding_only_genes <- setdiff(
  binding_genes,
  expression_genes
)

expression_only_genes <- setdiff(
  expression_genes,
  binding_genes
)

default_proteins <- intersect(
  c(
    "HNF4A",
    "PPARA",
    "FOXA2",
    "CTCF",
    "NR1H4",
    "NR5A2",
    "PROX1",
    "RXRA"
  ),
  protein_cols
)

if (length(default_proteins) == 0) {
  default_proteins <- head(
    protein_cols,
    8
  )
}


# =========================================================
# Tableau global d'exploration
# =========================================================

explore_dt <- binding_dt[
  ,
  lapply(
    .SD,
    function(x) {
      as.integer(
        any(
          x == 1,
          na.rm = TRUE
        )
      )
    }
  ),
  by = Name,
  .SDcols = protein_cols
]

explore_dt[
  ,
  n_bound_proteins := rowSums(
    .SD
  ),
  .SDcols = protein_cols
]

explore_summary <- explore_dt[
  ,
  .(
    Gene = Name,
    Bound_proteins = n_bound_proteins
  )
]

setorder(
  explore_summary,
  -Bound_proteins,
  Gene
)

max_bound_proteins <- if (
  nrow(explore_summary) > 0
) {
  max(
    explore_summary$Bound_proteins,
    na.rm = TRUE
  )
} else {
  0
}


# =========================================================
# UI
# =========================================================

ui <- fluidPage(

  tags$head(

    tags$style(
      HTML(
        "
        hr {
          border-top: 1px solid #cccccc;
        }

        .molli-title {
          margin-bottom: 20px;
        }

        .info-box {
          background-color: #f7f7f7;
          border: 1px solid #dddddd;
          border-radius: 5px;
          padding: 12px;
          margin-bottom: 15px;
        }

        .small-note {
          color: #666666;
          font-size: 12px;
        }

        .section-title {
          margin-top: 15px;
        }
        "
      )
    )
  ),

  div(

    class = "molli-title",

    titlePanel(
      "Interactive Multi-Omics Platform — MOLLI App"
    ),

    tags$p(
      paste(
        "Exploration intégrée des interactions protéines-régions régulatrices,",
        "des régions CRM/TSS et de l'expression différentielle RNA-seq."
      )
    )
  ),

  sidebarLayout(

    sidebarPanel(

      width = 3,

      # ---------------------------------------------------
      # Sélection globale des gènes
      # ---------------------------------------------------

      h4(
        "Sélection des gènes"
      ),

      textInput(
        inputId = "gene_to_add",
        label = "Ajouter un gène",
        placeholder = "Exemple : Cyp4a14"
      ),

      fluidRow(

        column(
          width = 6,
          actionButton(
            "add_gene",
            "Ajouter"
          )
        ),

        column(
          width = 6,
          actionButton(
            "clear_genes",
            "Effacer"
          )
        )
      ),

      tags$br(),

      selectizeInput(
        inputId = "selected_genes",
        label = "Gènes sélectionnés",
        choices = NULL,
        selected = character(0),
        multiple = TRUE,
        options = list(
          placeholder = "Tapez ou sélectionnez des gènes",
          plugins = list(
            "remove_button"
          )
        )
      ),

      tags$hr(),

      tags$div(

        class = "info-box",

        tags$strong(
          "Données disponibles"
        ),

        tags$br(),

        paste(
          length(binding_genes),
          "gènes interactome"
        ),

        tags$br(),

        paste(
          length(protein_cols),
          "protéines"
        ),

        tags$br(),

        paste(
          length(expression_genes),
          "gènes transcriptome"
        ),

        tags$br(),

        paste(
          length(all_studies),
          "comparaisons RNA-seq"
        ),

        tags$br(),

        paste(
          length(common_genes),
          "gènes communs"
        )
      ),

      # ---------------------------------------------------
      # Fiche biologique
      # ---------------------------------------------------

      h4("Fiche biologique"),

      selectInput(
        "gene_card_gene",
        "Gène à détailler",
        choices = all_genes,
        selected = if (length(common_genes) > 0) common_genes[1] else all_genes[1]
      ),

      tags$p(
        class = "small-note",
        "CRM/TSS, protéines associées et profil transcriptomique multi-études."
      ),

      tags$hr(),

      # ---------------------------------------------------
      # Options Interactome
      # ---------------------------------------------------

      h4(
        "Options Interactome"
      ),

      selectInput(
        inputId = "interactome_plot_type",
        label = "Type de figure",
        choices = c(
          "Heatmap gènes × protéines" = "gene_heatmap",
          "Co-binding protéines (Jaccard)" = "jaccard",
          "Fréquence des protéines" = "protein_frequency",
          "Combinaisons de protéines" = "upset",
          "Réseau gènes ↔ protéines" = "network"
        ),
        selected = "gene_heatmap"
      ),

      selectizeInput(
        inputId = "interactome_proteins",
        label = "Protéines à considérer",
        choices = protein_cols,
        selected = default_proteins,
        multiple = TRUE,
        options = list(
          plugins = list(
            "remove_button"
          )
        )
      ),

      checkboxInput(
        inputId = "keep_only_bound",
        label = "Masquer les protéines absentes des gènes sélectionnés",
        value = TRUE
      ),

      conditionalPanel(
        condition = "input.interactome_plot_type == 'gene_heatmap'",

        selectInput(
          inputId = "palette_name",
          label = "Palette heatmap",
          choices = c(
            "Rouge",
            "Bleu",
            "Vert",
            "Noir/Blanc",
            "Violet"
          ),
          selected = "Rouge"
        ),

        numericInput(
          inputId = "fontsize_row",
          label = "Taille police gènes",
          value = 12,
          min = 6,
          max = 20,
          step = 1
        ),

        numericInput(
          inputId = "fontsize_col",
          label = "Taille police protéines",
          value = 7,
          min = 4,
          max = 16,
          step = 1
        )
      ),

      conditionalPanel(
        condition = "input.interactome_plot_type == 'protein_frequency'",

        numericInput(
          inputId = "interactome_top_proteins",
          label = "Nombre de protéines à afficher",
          value = 20,
          min = 5,
          max = length(protein_cols),
          step = 1
        )
      ),

      conditionalPanel(
        condition = "input.interactome_plot_type == 'upset'",

        numericInput(
          inputId = "upset_top_combinations",
          label = "Nombre de combinaisons",
          value = 15,
          min = 5,
          max = 30,
          step = 1
        ),

        tags$p(
          class = "small-note",
          "Pour la lisibilité, les 8 premières protéines sélectionnées sont utilisées."
        )
      ),

      conditionalPanel(
        condition = "input.interactome_plot_type == 'network'",

        numericInput(
          inputId = "network_max_genes",
          label = "Nombre maximal de gènes",
          value = 20,
          min = 2,
          max = 40,
          step = 1
        )
      ),

      tags$br(),

      downloadButton(
        "download_interactome_png",
        "Figure Interactome PNG"
      ),

      tags$br(),
      tags$br(),

      downloadButton(
        "download_interactome_pdf",
        "Figure Interactome PDF"
      ),

      tags$hr(),

      # ---------------------------------------------------
      # Options Transcriptome
      # ---------------------------------------------------

      h4(
        "Options Transcriptome"
      ),

      selectInput(
        inputId = "transcriptome_plot_type",
        label = "Type de figure",
        choices = c(
          "Bubble Plot" = "bubble",
          "Heatmap des log2FC" = "heatmap",
          "Volcano Plot" = "volcano",
          "Profil multi-études d'un gène" = "profile",
          "Top gènes différentiellement exprimés" = "topde"
        ),
        selected = "bubble"
      ),

      selectizeInput(
        inputId = "expression_studies",
        label = "Études / comparaisons",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Sélectionner une ou plusieurs études..."
        )
      ),

      sliderTextInput(
        inputId = "FDRthresh",
        label = "Seuil FDR",
        choices = c(
          0.01,
          0.05,
          0.10,
          0.15,
          1
        ),
        selected = 0.05,
        grid = TRUE
      ),

      checkboxInput(
        inputId = "significant_only",
        label = "Masquer les résultats non significatifs",
        value = FALSE
      ),

      fluidRow(

        column(
          width = 6,
          colourInput(
            "colDown",
            "Couleur Down",
            "purple"
          )
        ),

        column(
          width = 6,
          colourInput(
            "colUp",
            "Couleur Up",
            "#FF00BF"
          )
        )
      ),

      conditionalPanel(
        condition = "input.transcriptome_plot_type == 'bubble'",

        sliderTextInput(
          inputId = "maxBubbleSize",
          label = "Taille maximale des bulles",
          choices = seq(
            5,
            30,
            2.5
          ),
          selected = 15,
          grid = TRUE
        )
      ),

      conditionalPanel(
        condition = "input.transcriptome_plot_type == 'volcano' || input.transcriptome_plot_type == 'topde'",

        selectInput(
          inputId = "transcriptome_single_study",
          label = "Comparaison",
          choices = all_studies,
          selected = if (
            length(all_studies) > 0
          ) {
            all_studies[1]
          } else {
            character(0)
          }
        )
      ),

      conditionalPanel(
        condition = "input.transcriptome_plot_type == 'topde'",

        numericInput(
          inputId = "top_de_n",
          label = "Nombre de gènes",
          value = 20,
          min = 5,
          max = 100,
          step = 5
        )
      ),

      tags$hr(),

      # ---------------------------------------------------
      # Options Intégration
      # ---------------------------------------------------

      h4(
        "Options Intégration"
      ),

      selectInput(
        inputId = "integration_plot_type",
        label = "Analyse intégrative",
        choices = c(
          "Études × protéines : enrichissement" = "enrichment_grid",
          "Liés vs non liés : log2FC" = "bound_vs_unbound",
          "Réseau binding + expression" = "integrated_network"
        ),
        selected = "enrichment_grid"
      ),

      selectInput(
        inputId = "integration_study",
        label = "Étude pour l'analyse ciblée",
        choices = all_studies,
        selected = if (
          length(all_studies) > 0
        ) {
          all_studies[1]
        } else {
          character(0)
        }
      ),

      selectInput(
        inputId = "integration_protein",
        label = "Protéine / régulateur",
        choices = protein_cols,
        selected = if (
          "HNF4A" %in% protein_cols
        ) {
          "HNF4A"
        } else {
          protein_cols[1]
        }
      ),

      selectizeInput(
        inputId = "integration_proteins_grid",
        label = "Protéines pour la matrice d'enrichissement",
        choices = protein_cols,
        selected = default_proteins,
        multiple = TRUE
      ),

      numericInput(
        inputId = "integration_network_top",
        label = "Top gènes pour le réseau intégré",
        value = 20,
        min = 5,
        max = 40,
        step = 5
      )
    ),

    # =====================================================
    # MAIN PANEL
    # =====================================================

    mainPanel(

      width = 9,

      tabsetPanel(

        id = "main_tabs",

        # =================================================
        # RÉSUMÉ
        # =================================================

        tabPanel(

          "Résumé",

          tags$br(),

          h3(
            "Résumé de la sélection"
          ),

          verbatimTextOutput(
            "summary_txt"
          ),

          tags$hr(),

          h4(
            "Statut des gènes"
          ),

          tableOutput(
            "genes_status"
          )
        ),

        # =================================================
        # INTERACTOME
        # =================================================

        tabPanel(

          "Interactome",

          tags$br(),

          h3(
            textOutput(
              "interactome_plot_title",
              inline = TRUE
            )
          ),

          uiOutput(
            "interactome_plot_description"
          ),

          plotOutput(
            "interactome_plot",
            height = "760px"
          )
        ),

        # =================================================
        # CRM / TSS
        # =================================================

        tabPanel(

          "Régions CRM / TSS",

          tags$br(),

          h3(
            "Régions régulatrices associées aux gènes sélectionnés"
          ),

          tags$p(
            paste(
              "Cette vue conserve les coordonnées de la région CRM et du TSS",
              "présents dans le fichier d'interactome."
            )
          ),

          tags$p(
            tags$strong(
              "Remarque : "
            ),
            paste(
              "avec le fichier actuellement chargé, chaque gène apparaît une seule fois.",
              "Cette vue ne doit donc pas être interprétée comme une collection",
              "de plusieurs promoteurs ou isoformes par gène."
            )
          ),

          selectizeInput(
            inputId = "crm_proteins",
            label = "Filtrer par protéine",
            choices = protein_cols,
            selected = character(0),
            multiple = TRUE
          ),

          verbatimTextOutput(
            "crm_summary"
          ),

          downloadButton(
            "download_crm",
            "Exporter CRM/TSS TSV"
          ),

          tags$hr(),

          DTOutput(
            "crm_table"
          )
        ),

        # =================================================
        # EXPLORER
        # =================================================

        tabPanel(

          "Explorer les gènes",

          tags$br(),

          h3(
            "Exploration globale de l'interactome"
          ),

          sliderInput(
            inputId = "min_bound_proteins",
            label = "Nombre minimum de protéines liées",
            min = 0,
            max = max_bound_proteins,
            value = 0,
            step = 1
          ),

          tags$p(
            tags$strong(
              "Cliquez sur une ligne pour ajouter automatiquement le gène à la sélection."
            )
          ),

          DTOutput(
            "available_genes_table"
          ),

          tags$hr(),

          h4(
            "Distribution du nombre de protéines liées par gène"
          ),

          plotOutput(
            "bound_protein_hist",
            height = "350px"
          )
        ),

        # =================================================
        # TRANSCRIPTOME
        # =================================================

        tabPanel(

          "Transcriptome",

          tags$br(),

          h3(
            textOutput(
              "transcriptome_plot_title",
              inline = TRUE
            )
          ),

          uiOutput(
            "transcriptome_plot_description"
          ),

          verbatimTextOutput(
            "expression_summary"
          ),

          plotOutput(
            "transcriptome_plot",
            height = "700px"
          ),

          fluidRow(

            column(
              width = 3,
              downloadButton(
                "download_transcriptome_svg",
                "Figure Transcriptome SVG"
              )
            ),

            column(
              width = 3,
              downloadButton(
                "download_expression",
                "Données RNA-seq TSV"
              )
            )
          ),

          tags$hr(),

          h4(
            "Valeurs"
          ),

          DTOutput(
            "expression_table"
          )
        ),

        # =================================================
        # FICHE GÈNE
        # =================================================

        tabPanel(
          "Fiche gène",
          tags$br(),
          h3(textOutput("gene_card_title", inline = TRUE)),
          fluidRow(
            column(
              5,
              h4("Informations régulatrices"),
              verbatimTextOutput("gene_card_regulatory"),
              h4("Protéines associées"),
              DTOutput("gene_card_proteins")
            ),
            column(
              7,
              h4("Profil transcriptomique"),
              plotOutput("gene_card_expression_plot", height="550px")
            )
          ),
          tags$hr(),
          h4("Résultats RNA-seq"),
          DTOutput("gene_card_expression_table")
        ),

        # =================================================
        # INTÉGRATION
        # =================================================

        tabPanel(

          "Intégration",

          tags$br(),

          h3(
            "Intégration Interactome × Transcriptome"
          ),

          tags$p(
            paste(
              "Cette section croise les gènes communs aux deux jeux de données",
              "afin d'explorer l'association entre liaison d'une protéine",
              "et réponse transcriptionnelle."
            )
          ),

          verbatimTextOutput(
            "integration_summary"
          ),

          h4("Enrichissement directionnel Up / Down"),
          tableOutput("directional_enrichment_table"),

          plotOutput(
            "integration_plot",
            height = "720px"
          ),

          tags$hr(),

          h4(
            "Table binding × expression"
          ),

          downloadButton(
            "download_integration",
            "Exporter table intégrée TSV"
          ),

          tags$br(),
          tags$br(),

          DTOutput(
            "integration_table"
          )
        ),

        # =================================================
        # DATASET / QC
        # =================================================

        tabPanel(

          "Dataset / QC",

          tags$br(),

          h3(
            "Résumé des données chargées"
          ),

          verbatimTextOutput(
            "dataset_qc_txt"
          ),

          tags$hr(),

          h4(
            "Couverture des gènes"
          ),

          plotOutput(
            "dataset_overlap_plot",
            height = "350px"
          ),

          tags$hr(),

          h4(
            "Fréquence des protéines dans l'interactome"
          ),

          plotOutput(
            "protein_frequency_qc",
            height = "500px"
          )
        ),

        # =================================================
        # AIDE
        # =================================================

        tabPanel(

          "Aide",

          tags$br(),

          h3(
            "Utilisation de MOLLI App"
          ),

          tags$h4(
            "Interactome"
          ),

          tags$ul(

            tags$li(
              "La heatmap gènes × protéines représente les signaux de liaison binaires."
            ),

            tags$li(
              "La heatmap Jaccard compare les profils de liaison des protéines sur l'ensemble des gènes."
            ),

            tags$li(
              "La vue Combinaisons montre les associations de protéines les plus fréquentes parmi les gènes."
            ),

            tags$li(
              "Le réseau biparti affiche les gènes sélectionnés et les protéines qui leur sont associées."
            )
          ),

          tags$h4(
            "Régions CRM / TSS"
          ),

          tags$ul(

            tags$li(
              "Les coordonnées CRM et TSS proviennent directement du fichier d'interactome."
            ),

            tags$li(
              paste(
                "Dans le fichier actuellement chargé, un gène correspond à une seule ligne.",
                "La vue CRM/TSS n'est donc pas une annotation multi-promoteurs."
              )
            )
          ),

          tags$h4(
            "Transcriptome"
          ),

          tags$ul(

            tags$li(
              "Le Bubble Plot compare les gènes sélectionnés sur plusieurs études."
            ),

            tags$li(
              "La heatmap représente les log2FC signés reconstruits à partir de abslogFC et State."
            ),

            tags$li(
              "Le Volcano Plot utilise l'ensemble des gènes disponibles pour la comparaison choisie."
            ),

            tags$li(
              "Le profil multi-études suit un ou plusieurs gènes dans les comparaisons sélectionnées."
            ),

            tags$li(
              "La vue Top DE affiche les gènes ayant les plus grands |log2FC| parmi les résultats significatifs."
            )
          ),

          tags$h4(
            "Intégration"
          ),

          tags$ul(

            tags$li(
              paste(
                "L'enrichissement Études × protéines utilise un test exact de Fisher",
                "sur les gènes communs aux deux datasets, avec correction BH",
                "dans la matrice globale."
              )
            ),

            tags$li(
              paste(
                "L'analyse ciblée distingue les enrichissements Up et Down",
                "et compare les log2FC des gènes liés et non liés par Wilcoxon."
              )
            ),

            tags$li(
              paste(
                "La comparaison Liés vs non liés montre la distribution des log2FC",
                "pour la protéine et l'étude choisies."
              )
            ),

            tags$li(
              paste(
                "Le réseau intégré associe protéines et gènes significatifs,",
                "avec une taille de nœud proportionnelle au |log2FC|."
              )
            ),

            tags$li(
              paste(
                "Une liaison à une région régulatrice associée à un gène",
                "et une variation d'expression après perturbation",
                "ne démontrent pas à elles seules une relation causale directe."
              )
            )
          )
        )
      )
    )
  )
)


# =========================================================
# SERVER
# =========================================================

server <- function(
    input,
    output,
    session
) {

  # =======================================================
  # Initialisations
  # =======================================================

  updateSelectizeInput(
    session = session,
    inputId = "selected_genes",
    choices = all_genes,
    selected = character(0),
    server = TRUE
  )

  updateSelectizeInput(
    session = session,
    inputId = "expression_studies",
    choices = all_studies,
    selected = if (
      length(all_studies) > 0
    ) {
      all_studies[1]
    } else {
      character(0)
    },
    server = TRUE
  )


  # =======================================================
  # Sélection gènes
  # =======================================================

  observeEvent(
    input$add_gene,
    {

      new_gene <- trimws(
        input$gene_to_add
      )

      if (!nzchar(new_gene)) {
        return(NULL)
      }

      matched_gene <- all_genes[
        tolower(all_genes) ==
          tolower(new_gene)
      ]

      if (
        length(matched_gene) == 0
      ) {

        showNotification(
          paste(
            "Gène non trouvé :",
            new_gene
          ),
          type = "warning"
        )

        updateTextInput(
          session,
          "gene_to_add",
          value = ""
        )

        return(NULL)
      }

      current_genes <- input$selected_genes

      if (is.null(current_genes)) {
        current_genes <- character(0)
      }

      updated_genes <- unique(
        c(
          current_genes,
          matched_gene[1]
        )
      )

      updateSelectizeInput(
        session = session,
        inputId = "selected_genes",
        choices = all_genes,
        selected = updated_genes,
        server = TRUE
      )

      updateTextInput(
        session,
        "gene_to_add",
        value = ""
      )
    }
  )


  observeEvent(
    input$clear_genes,
    {

      updateSelectizeInput(
        session = session,
        inputId = "selected_genes",
        choices = all_genes,
        selected = character(0),
        server = TRUE
      )

      updateTextInput(
        session,
        "gene_to_add",
        value = ""
      )
    }
  )


  # =======================================================
  # Explorer
  # =======================================================

  filtered_explore_summary <- reactive({

    explore_summary[
      Bound_proteins >=
        input$min_bound_proteins
    ]
  })


  output$available_genes_table <- renderDT({

    datatable(
      filtered_explore_summary(),
      selection = "single",
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        order = list(
          list(
            1,
            "desc"
          )
        )
      )
    )
  })


  observeEvent(
    input$available_genes_table_rows_selected,
    {

      selected_row <- input$available_genes_table_rows_selected

      if (length(selected_row) == 0) {
        return(NULL)
      }

      gene_clicked <- filtered_explore_summary()[
        selected_row,
        Gene
      ]

      current_genes <- input$selected_genes

      if (is.null(current_genes)) {
        current_genes <- character(0)
      }

      updated_genes <- unique(
        c(
          current_genes,
          gene_clicked
        )
      )

      updateSelectizeInput(
        session = session,
        inputId = "selected_genes",
        choices = all_genes,
        selected = updated_genes,
        server = TRUE
      )

      showNotification(
        paste(
          gene_clicked,
          "ajouté à la sélection"
        ),
        type = "message",
        duration = 2
      )
    }
  )


  output$bound_protein_hist <- renderPlot({

    hist(
      explore_summary$Bound_proteins,
      breaks = 30,
      main = "Nombre de protéines liées par gène",
      xlab = "Nombre de protéines liées",
      ylab = "Nombre de gènes"
    )
  })


  # =======================================================
  # Résumé général
  # =======================================================

  gene_matrix <- reactive({

    genes <- input$selected_genes

    req(
      length(genes) > 0
    )

    proteins <- input$interactome_proteins

    if (
      is.null(proteins) ||
      length(proteins) == 0
    ) {
      proteins <- protein_cols
    }

    aggregate_to_gene_matrix(
      dt = binding_dt,
      genes = genes,
      protein_cols = proteins,
      keep_only_bound = input$keep_only_bound
    )
  })


  output$summary_txt <- renderText({

    genes <- input$selected_genes

    if (
      is.null(genes) ||
      length(genes) == 0
    ) {
      return("Aucun gène sélectionné.")
    }

    binding_found <- genes[
      genes %in%
        binding_genes
    ]

    expression_found <- genes[
      genes %in%
        expression_genes
    ]

    common_found <- genes[
      genes %in%
        common_genes
    ]

    mat <- tryCatch(
      gene_matrix(),
      error = function(e) NULL
    )

    heatmap_text <- if (
      is.null(mat)
    ) {
      "Aucune matrice d'interactome disponible."
    } else {
      paste0(
        "Dimensions de la matrice courante : ",
        nrow(mat),
        " gènes × ",
        ncol(mat),
        " protéines"
      )
    }

    paste0(
      "Gènes sélectionnés : ",
      length(genes),
      "\n",
      "Présents dans l'interactome : ",
      length(binding_found),
      "\n",
      "Présents dans le transcriptome : ",
      length(expression_found),
      "\n",
      "Présents dans les deux datasets : ",
      length(common_found),
      "\n",
      heatmap_text
    )
  })


  output$genes_status <- renderTable({

    genes <- input$selected_genes

    if (
      is.null(genes) ||
      length(genes) == 0
    ) {
      return(NULL)
    }

    data.frame(
      Gene = genes,
      Interactome = ifelse(
        genes %in% binding_genes,
        "Présent",
        "Absent"
      ),
      Transcriptome = ifelse(
        genes %in% expression_genes,
        "Présent",
        "Absent"
      ),
      Commun = ifelse(
        genes %in% common_genes,
        "Oui",
        "Non"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })


  # =======================================================
  # CRM / TSS
  # =======================================================

  crm_data <- reactive({

    genes <- input$selected_genes

    if (
      is.null(genes) ||
      length(genes) == 0
    ) {
      return(
        crm_protein_dt[0]
      )
    }

    dat <- crm_protein_dt[
      Name %in% genes
    ]

    proteins <- input$crm_proteins

    if (
      !is.null(proteins) &&
      length(proteins) > 0
    ) {
      dat <- dat[
        Protein %in% proteins
      ]
    }

    dat[
      ,
      gene_order := match(
        Name,
        genes
      )
    ]

    setorder(
      dat,
      gene_order,
      Protein
    )

    dat
  })


  output$crm_summary <- renderText({

    dat <- crm_data()
    genes <- input$selected_genes

    if (
      is.null(genes) ||
      length(genes) == 0
    ) {
      return("Aucun gène sélectionné.")
    }

    if (nrow(dat) == 0) {
      return(
        paste0(
          "Gènes sélectionnés : ",
          length(genes),
          "\n",
          "Aucune interaction CRM-protéine trouvée avec les filtres actuels."
        )
      )
    }

    paste0(
      "Gènes avec interaction : ",
      uniqueN(dat$Name),
      "\n",
      "Régions CRM distinctes : ",
      uniqueN(dat$CRM_ID),
      "\n",
      "Protéines distinctes : ",
      uniqueN(dat$Protein),
      "\n",
      "Interactions CRM × protéine : ",
      nrow(dat)
    )
  })


  output$crm_table <- renderDT({

    dat <- crm_data()

    if (nrow(dat) == 0) {
      return(
        datatable(
          data.frame(
            Information = "Aucune donnée CRM/TSS disponible."
          ),
          rownames = FALSE,
          options = list(
            dom = "t"
          )
        )
      )
    }

    display_data <- dat[
      ,
      .(
        Gene = Name,
        Protein,
        CRM = CRM_ID,
        CRM_chr = chrCRM,
        CRM_start = startCRM,
        CRM_end = endCRM,
        CRM_size,
        TSS = TSS_ID,
        TSS_chr = tssC,
        TSS_position = tssS,
        inter
      )
    ]

    datatable(
      display_data,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 15,
        autoWidth = TRUE,
        scrollX = TRUE
      )
    )
  })


  output$download_crm <- downloadHandler(

    filename = function() {
      paste0(
        "MOLLI_CRM_TSS_",
        Sys.Date(),
        ".tsv"
      )
    },

    content = function(file) {

      dat <- crm_data()

      fwrite(
        dat[
          ,
          .(
            Gene = Name,
            Protein,
            CRM = CRM_ID,
            CRM_chr = chrCRM,
            CRM_start = startCRM,
            CRM_end = endCRM,
            CRM_size,
            TSS = TSS_ID,
            TSS_chr = tssC,
            TSS_position = tssS,
            inter
          )
        ],
        file = file,
        sep = "\t"
      )
    }
  )


  # =======================================================
  # INTERACTOME - titres / descriptions
  # =======================================================

  output$interactome_plot_title <- renderText({

    switch(
      input$interactome_plot_type,
      "gene_heatmap" = "Heatmap gènes × protéines",
      "jaccard" = "Co-binding des protéines — indice de Jaccard",
      "protein_frequency" = "Fréquence des protéines dans l'interactome",
      "upset" = "Combinaisons de protéines les plus fréquentes",
      "network" = "Réseau biparti gènes ↔ protéines",
      "Interactome"
    )
  })


  output$interactome_plot_description <- renderUI({

    txt <- switch(
      input$interactome_plot_type,

      "gene_heatmap" = paste(
        "Chaque case indique l'absence ou la présence du signal de liaison",
        "pour un gène et une protéine."
      ),

      "jaccard" = paste(
        "La similarité Jaccard compare les ensembles de gènes associés à deux protéines.",
        "Une valeur proche de 1 indique des profils de liaison fortement partagés."
      ),

      "protein_frequency" = paste(
        "Le barplot représente le nombre de gènes associés à chaque protéine."
      ),

      "upset" = paste(
        "Cette représentation de type UpSet résume les combinaisons de protéines",
        "observées sur les gènes."
      ),

      "network" = paste(
        "Les gènes sont reliés aux protéines dont le signal de liaison est présent.",
        "Le réseau est limité aux gènes sélectionnés pour rester lisible."
      ),

      ""
    )

    tags$p(txt)
  })


  # =======================================================
  # INTERACTOME - préparation des figures
  # =======================================================

  interactome_plot_object <- reactive({

    plot_type <- input$interactome_plot_type

    proteins <- input$interactome_proteins

    if (
      is.null(proteins) ||
      length(proteins) == 0
    ) {
      proteins <- protein_cols
    }

    if (plot_type == "gene_heatmap") {

      mat <- gene_matrix()

      validate(
        need(
          !is.null(mat),
          "Sélectionnez au moins un gène présent dans l'interactome."
        )
      )

      return(
        list(
          type = "pheatmap",
          mat = mat
        )
      )
    }


    if (plot_type == "jaccard") {

      jac <- compute_jaccard_matrix(
        binding_dt,
        proteins
      )

      validate(
        need(
          !is.null(jac),
          "Sélectionnez au moins une protéine."
        )
      )

      jac_dt <- as.data.table(
        as.table(jac)
      )

      setnames(
        jac_dt,
        c(
          "Protein1",
          "Protein2",
          "Jaccard"
        )
      )

      p <- ggplot(
        jac_dt,
        aes(
          x = Protein1,
          y = Protein2,
          fill = Jaccard
        )
      ) +
        geom_tile(
          colour = "white"
        ) +
        scale_fill_gradient(
          low = "white",
          high = "black",
          limits = c(
            0,
            1
          )
        ) +
        labs(
          x = NULL,
          y = NULL,
          fill = "Jaccard"
        ) +
        theme_classic(
          base_size = 11
        ) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          )
        )

      return(
        list(
          type = "ggplot",
          plot = p
        )
      )
    }


    if (plot_type == "protein_frequency") {

      freq <- data.table(
        Protein = proteins,
        Bound_genes = sapply(
          proteins,
          function(protein) {
            sum(
              binding_dt[[protein]] == 1,
              na.rm = TRUE
            )
          }
        )
      )

      setorder(
        freq,
        -Bound_genes
      )

      top_n <- min(
        as.integer(
          input$interactome_top_proteins
        ),
        nrow(freq)
      )

      freq <- head(
        freq,
        top_n
      )

      p <- ggplot(
        freq,
        aes(
          x = reorder(
            Protein,
            Bound_genes
          ),
          y = Bound_genes
        )
      ) +
        geom_col() +
        coord_flip() +
        labs(
          x = "Protéines",
          y = "Nombre de gènes associés"
        ) +
        theme_classic(
          base_size = 12
        )

      return(
        list(
          type = "ggplot",
          plot = p
        )
      )
    }


    if (plot_type == "upset") {

      dat <- make_upset_like_data(
        binding_dt,
        proteins,
        top_n = as.integer(
          input$upset_top_combinations
        )
      )

      validate(
        need(
          nrow(dat) > 0,
          "Aucune combinaison de protéines à afficher."
        )
      )

      p <- ggplot(
        dat,
        aes(
          x = reorder(
            Combination,
            N
          ),
          y = N
        )
      ) +
        geom_col() +
        coord_flip() +
        labs(
          x = "Combinaison de protéines",
          y = "Nombre de gènes"
        ) +
        theme_classic(
          base_size = 11
        )

      return(
        list(
          type = "ggplot",
          plot = p
        )
      )
    }


    # Réseau biparti gènes ↔ protéines

    genes <- input$selected_genes

    validate(
      need(
        !is.null(genes) &&
          length(genes) > 0,
        "Sélectionnez des gènes pour construire le réseau."
      )
    )

    max_genes <- as.integer(
      input$network_max_genes
    )

    genes <- head(
      genes[
        genes %in%
          binding_genes
      ],
      max_genes
    )

    validate(
      need(
        length(genes) > 0,
        "Aucun des gènes sélectionnés n'est présent dans l'interactome."
      )
    )

    net_dt <- binding_dt[
      Name %in% genes,
      c(
        "Name",
        proteins
      ),
      with = FALSE
    ]

    net_long <- melt(
      net_dt,
      id.vars = "Name",
      measure.vars = proteins,
      variable.name = "Protein",
      value.name = "Binding"
    )

    net_long <- net_long[
      Binding == 1
    ]

    validate(
      need(
        nrow(net_long) > 0,
        "Aucune interaction dans cette sélection."
      )
    )

    gene_levels <- unique(
      net_long$Name
    )

    protein_levels <- unique(
      net_long$Protein
    )

    gene_y <- setNames(
      seq_along(gene_levels),
      gene_levels
    )

    protein_y <- setNames(
      seq(
        1,
        max(
          length(gene_levels),
          length(protein_levels)
        ),
        length.out = length(protein_levels)
      ),
      protein_levels
    )

    edges <- net_long[
      ,
      .(
        x = 1,
        y = gene_y[Name],
        xend = 2,
        yend = protein_y[Protein]
      )
    ]

    gene_nodes <- data.table(
      label = gene_levels,
      x = 1,
      y = unname(
        gene_y[
          gene_levels
        ]
      ),
      type = "Gène"
    )

    protein_nodes <- data.table(
      label = protein_levels,
      x = 2,
      y = unname(
        protein_y[
          protein_levels
        ]
      ),
      type = "Protéine"
    )

    nodes <- rbind(
      gene_nodes,
      protein_nodes
    )

    p <- ggplot() +
      geom_segment(
        data = edges,
        aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend
        ),
        alpha = 0.35
      ) +
      geom_point(
        data = nodes,
        aes(
          x = x,
          y = y,
          shape = type
        ),
        size = 4
      ) +
      geom_text(
        data = gene_nodes,
        aes(
          x = x,
          y = y,
          label = label
        ),
        hjust = 1.1,
        size = 3
      ) +
      geom_text(
        data = protein_nodes,
        aes(
          x = x,
          y = y,
          label = label
        ),
        hjust = -0.1,
        size = 3
      ) +
      scale_x_continuous(
        limits = c(
          0.5,
          2.5
        ),
        breaks = c(
          1,
          2
        ),
        labels = c(
          "Gènes",
          "Protéines"
        )
      ) +
      labs(
        x = NULL,
        y = NULL,
        shape = NULL
      ) +
      theme_void() +
      theme(
        legend.position = "bottom"
      )

    list(
      type = "ggplot",
      plot = p
    )
  })


  draw_interactome_plot <- function(obj) {

    if (
      is.null(obj)
    ) {
      plot.new()
      title("Aucune donnée")
      return(
        invisible(NULL)
      )
    }

    if (
      obj$type == "pheatmap"
    ) {

      make_heatmap(
        mat = obj$mat,
        palette_name = input$palette_name,
        fontsize_row = input$fontsize_row,
        fontsize_col = input$fontsize_col
      )

    } else {

      print(
        obj$plot
      )
    }
  }


  output$interactome_plot <- renderPlot({

    draw_interactome_plot(
      interactome_plot_object()
    )
  })


  output$download_interactome_png <- downloadHandler(

    filename = function() {
      paste0(
        "MOLLI_interactome_",
        input$interactome_plot_type,
        "_",
        Sys.Date(),
        ".png"
      )
    },

    content = function(file) {

      png(
        file,
        width = 1800,
        height = 1200,
        res = 150
      )

      draw_interactome_plot(
        interactome_plot_object()
      )

      dev.off()
    }
  )


  output$download_interactome_pdf <- downloadHandler(

    filename = function() {
      paste0(
        "MOLLI_interactome_",
        input$interactome_plot_type,
        "_",
        Sys.Date(),
        ".pdf"
      )
    },

    content = function(file) {

      pdf(
        file,
        width = 14,
        height = 10
      )

      draw_interactome_plot(
        interactome_plot_object()
      )

      dev.off()
    }
  )


  # =======================================================
  # TRANSCRIPTOME - données
  # =======================================================

  expression_data <- reactive({

    genes <- input$selected_genes

    req(
      length(genes) > 0
    )

    studies <- input$expression_studies

    if (
      is.null(studies) ||
      length(studies) == 0
    ) {
      return(
        rnaseq_dt[0]
      )
    }

    dat <- rnaseq_dt[
      gene %in% genes &
        StudiesCond %in% studies
    ]

    if (
      isTRUE(
        input$significant_only
      )
    ) {

      threshold <- as.numeric(
        input$FDRthresh
      )

      dat <- dat[
        !is.na(FDR) &
          FDR < threshold
      ]
    }

    dat[
      ,
      gene_order := match(
        gene,
        genes
      )
    ]

    setorder(
      dat,
      gene_order
    )

    dat
  })


  output$expression_summary <- renderText({

    genes <- input$selected_genes

    if (
      is.null(genes) ||
      length(genes) == 0
    ) {
      return(
        "Aucun gène sélectionné."
      )
    }

    dat <- expression_data()

    if (nrow(dat) == 0) {
      return(
        paste0(
          "Gènes sélectionnés : ",
          length(genes),
          "\n",
          "Aucun résultat d'expression disponible avec les filtres actuels."
        )
      )
    }

    n_significant <- sum(
      !is.na(dat$FDR) &
        dat$FDR <
        as.numeric(
          input$FDRthresh
        )
    )

    paste0(
      "Gènes sélectionnés : ",
      length(genes),
      "\n",
      "Gènes retrouvés : ",
      uniqueN(dat$gene),
      "\n",
      "Comparaisons affichées : ",
      uniqueN(dat$StudiesCond),
      "\n",
      "Observations : ",
      nrow(dat),
      "\n",
      "Up : ",
      sum(
        dat$State == "up",
        na.rm = TRUE
      ),
      "\n",
      "Down : ",
      sum(
        dat$State == "down",
        na.rm = TRUE
      ),
      "\n",
      "FDR < ",
      input$FDRthresh,
      " : ",
      n_significant
    )
  })


  output$transcriptome_plot_title <- renderText({

    switch(
      input$transcriptome_plot_type,
      "bubble" = "Bubble Plot",
      "heatmap" = "Heatmap des log2FC signés",
      "volcano" = "Volcano Plot",
      "profile" = "Profil multi-études",
      "topde" = "Top gènes différentiellement exprimés",
      "Transcriptome"
    )
  })


  output$transcriptome_plot_description <- renderUI({

    txt <- switch(
      input$transcriptome_plot_type,

      "bubble" = paste(
        "Taille = |log2FC| ; couleur = direction de variation ;",
        "transparence = significativité."
      ),

      "heatmap" = paste(
        "La couleur représente le log2FC signé pour les gènes",
        "et comparaisons sélectionnés."
      ),

      "volcano" = paste(
        "Le Volcano Plot utilise tous les gènes disponibles",
        "pour la comparaison sélectionnée."
      ),

      "profile" = paste(
        "Cette vue suit les gènes sélectionnés à travers",
        "les différentes comparaisons RNA-seq."
      ),

      "topde" = paste(
        "Cette vue affiche les gènes significatifs présentant",
        "les plus grands |log2FC| pour une comparaison."
      ),

      ""
    )

    tags$p(txt)
  })


  transcriptome_plot_object <- reactive({

    plot_type <- input$transcriptome_plot_type
    threshold <- as.numeric(
      input$FDRthresh
    )

    if (plot_type == "bubble") {

      dat <- copy(
        expression_data()
      )

      validate(
        need(
          nrow(dat) > 0,
          "Aucune donnée d'expression à afficher."
        )
      )

      dat[
        ,
        FDR_status := fifelse(
          !is.na(FDR) &
            FDR < threshold,
          "Significatif",
          "Non significatif"
        )
      ]

      p <- ggplot(
        dat,
        aes(
          x = gene,
          y = StudiesCond,
          size = abslogFC,
          colour = State,
          alpha = FDR_status
        )
      ) +
        geom_point() +
        scale_colour_manual(
          values = c(
            "up" = input$colUp,
            "down" = input$colDown
          )
        ) +
        scale_alpha_manual(
          values = c(
            "Significatif" = 1,
            "Non significatif" = 0.3
          )
        ) +
        scale_size_area(
          max_size = as.numeric(
            input$maxBubbleSize
          )
        ) +
        labs(
          x = "Gènes",
          y = "Études / comparaisons",
          size = "|log2FC|",
          colour = "État",
          alpha = "FDR"
        ) +
        theme_classic(
          base_size = 12
        ) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          )
        )

      return(p)
    }


    if (plot_type == "heatmap") {

      dat <- copy(
        expression_data()
      )

      validate(
        need(
          nrow(dat) > 0,
          "Aucune donnée d'expression à afficher."
        )
      )

      heat_dat <- dat[
        ,
        .(
          signed_logFC = mean(
            signed_logFC,
            na.rm = TRUE
          )
        ),
        by = .(
          gene,
          StudiesCond
        )
      ]

      p <- ggplot(
        heat_dat,
        aes(
          x = StudiesCond,
          y = gene,
          fill = signed_logFC
        )
      ) +
        geom_tile(
          colour = "grey85"
        ) +
        scale_fill_gradient2(
          low = input$colDown,
          mid = "white",
          high = input$colUp,
          midpoint = 0
        ) +
        labs(
          x = "Études / comparaisons",
          y = "Gènes",
          fill = "log2FC"
        ) +
        theme_classic(
          base_size = 12
        ) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          )
        )

      return(p)
    }


    if (plot_type == "volcano") {

      study <- input$transcriptome_single_study

      dat <- copy(
        rnaseq_dt[
          StudiesCond == study
        ]
      )

      validate(
        need(
          nrow(dat) > 0,
          "Aucune donnée disponible pour cette comparaison."
        )
      )

      dat[
        ,
        significance := fifelse(
          !is.na(FDR) &
            FDR < threshold,
          State,
          "non_significant"
        )
      ]

      dat[
        ,
        neg_log10_FDR := safe_neg_log10(
          FDR
        )
      ]

      p <- ggplot(
        dat,
        aes(
          x = signed_logFC,
          y = neg_log10_FDR,
          colour = significance
        )
      ) +
        geom_point(
          alpha = 0.65
        ) +
        geom_hline(
          yintercept = -log10(
            threshold
          ),
          linetype = "dashed"
        ) +
        geom_vline(
          xintercept = 0,
          linetype = "dashed"
        ) +
        scale_colour_manual(
          values = c(
            "up" = input$colUp,
            "down" = input$colDown,
            "non_significant" = "grey70"
          )
        ) +
        labs(
          title = study,
          x = "log2FC signé",
          y = "-log10(FDR)",
          colour = "État"
        ) +
        theme_classic(
          base_size = 12
        )

      return(p)
    }


    if (plot_type == "profile") {

      dat <- copy(
        expression_data()
      )

      validate(
        need(
          nrow(dat) > 0,
          "Sélectionnez au moins un gène et une étude."
        )
      )

      dat[
        ,
        Significant := !is.na(FDR) &
          FDR < threshold
      ]

      p <- ggplot(
        dat,
        aes(
          x = StudiesCond,
          y = signed_logFC,
          fill = State,
          alpha = Significant
        )
      ) +
        geom_col(
          position = "dodge"
        ) +
        facet_wrap(
          ~ gene,
          scales = "free_x"
        ) +
        scale_fill_manual(
          values = c(
            "up" = input$colUp,
            "down" = input$colDown
          )
        ) +
        scale_alpha_manual(
          values = c(
            "TRUE" = 1,
            "FALSE" = 0.35
          )
        ) +
        geom_hline(
          yintercept = 0
        ) +
        labs(
          x = "Études / comparaisons",
          y = "log2FC signé",
          fill = "État",
          alpha = paste0(
            "FDR < ",
            threshold
          )
        ) +
        theme_classic(
          base_size = 11
        ) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          )
        )

      return(p)
    }


    # Top DE

    study <- input$transcriptome_single_study

    dat <- copy(
      rnaseq_dt[
        StudiesCond == study &
          !is.na(FDR) &
          FDR < threshold
      ]
    )

    validate(
      need(
        nrow(dat) > 0,
        "Aucun gène significatif avec ce seuil FDR."
      )
    )

    setorder(
      dat,
      -abslogFC
    )

    dat <- head(
      dat,
      as.integer(
        input$top_de_n
      )
    )

    p <- ggplot(
      dat,
      aes(
        x = reorder(
          gene,
          abslogFC
        ),
        y = signed_logFC,
        fill = State
      )
    ) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(
        values = c(
          "up" = input$colUp,
          "down" = input$colDown
        )
      ) +
      labs(
        title = study,
        x = "Gènes",
        y = "log2FC signé",
        fill = "État"
      ) +
      theme_classic(
        base_size = 12
      )

    p
  })


  output$transcriptome_plot <- renderPlot({

    print(
      transcriptome_plot_object()
    )
  })


  output$expression_table <- renderDT({

    dat <- expression_data()

    if (nrow(dat) == 0) {

      return(
        datatable(
          data.frame(
            Information = "Aucune donnée d'expression disponible."
          ),
          rownames = FALSE,
          options = list(
            dom = "t"
          )
        )
      )
    }

    display_data <- dat %>%
      arrange(FDR) %>%
      select(
        gene,
        StudiesCond,
        signed_logFC,
        abslogFC,
        FDR,
        State
      )

    datatable(
      display_data,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 10,
        autoWidth = TRUE
      )
    )
  })


  output$download_transcriptome_svg <- downloadHandler(

    filename = function() {
      paste0(
        "MOLLI_transcriptome_",
        input$transcriptome_plot_type,
        "_",
        Sys.Date(),
        ".svg"
      )
    },

    content = function(file) {

      svg(
        filename = file,
        width = 12,
        height = 8
      )

      print(
        transcriptome_plot_object()
      )

      dev.off()
    }
  )


  output$download_expression <- downloadHandler(

    filename = function() {
      paste0(
        "MOLLI_expression_",
        Sys.Date(),
        ".tsv"
      )
    },

    content = function(file) {

      dat <- expression_data()

      fwrite(
        dat[
          ,
          .(
            gene,
            StudiesCond,
            signed_logFC,
            abslogFC,
            FDR,
            State
          )
        ],
        file = file,
        sep = "\t"
      )
    }
  )


  # =======================================================
  # FICHE BIOLOGIQUE PAR GÈNE
  # =======================================================

  observeEvent(input$selected_genes, {
    genes <- input$selected_genes
    if (!is.null(genes) && length(genes) > 0) {
      updateSelectInput(session, "gene_card_gene", choices=all_genes, selected=genes[1])
    }
  }, ignoreInit=TRUE)

  gene_card_binding <- reactive({ binding_dt[Name == input$gene_card_gene] })
  gene_card_expression <- reactive({ rnaseq_dt[gene == input$gene_card_gene] })

  output$gene_card_title <- renderText({ paste("Fiche biologique —", input$gene_card_gene) })

  output$gene_card_regulatory <- renderText({
    dat <- gene_card_binding()
    if (nrow(dat)==0) return("Ce gène n'est pas présent dans le fichier interactome.")
    vals <- unlist(dat[1, ..protein_cols], use.names=FALSE)
    paste0(
      "CRM : ", dat$chrCRM[1], ":", dat$startCRM[1], "-", dat$endCRM[1], "\n",
      "Taille CRM : ", dat$endCRM[1]-dat$startCRM[1], " pb\n",
      "TSS : ", dat$tssC[1], ":", dat$tssS[1], "\n",
      "Nombre de protéines associées : ", sum(as.numeric(vals)==1, na.rm=TRUE)
    )
  })

  output$gene_card_proteins <- renderDT({
    dat <- gene_card_binding()
    if (nrow(dat)==0) return(datatable(data.frame(Information="Aucune donnée interactome."), rownames=FALSE, options=list(dom="t")))
    vals <- unlist(dat[1, ..protein_cols], use.names=FALSE)
    active <- protein_cols[as.numeric(vals)==1]
    datatable(data.frame(Protein=active), rownames=FALSE, options=list(pageLength=12, dom="tip"))
  })

  output$gene_card_expression_plot <- renderPlot({
    dat <- copy(gene_card_expression())
    validate(need(nrow(dat)>0, "Ce gène n'est pas présent dans le transcriptome."))
    threshold <- as.numeric(input$FDRthresh)
    dat[, Significant := !is.na(FDR) & FDR < threshold]
    print(ggplot(dat, aes(x=reorder(StudiesCond, signed_logFC), y=signed_logFC, fill=State, alpha=Significant)) +
      geom_col() + coord_flip() + geom_hline(yintercept=0) +
      scale_fill_manual(values=c("up"=input$colUp, "down"=input$colDown)) +
      scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.3)) +
      labs(x="Études / comparaisons", y="log2FC signé", fill="État", alpha=paste0("FDR < ",threshold)) +
      theme_classic(base_size=11))
  })

  output$gene_card_expression_table <- renderDT({
    dat <- gene_card_expression()
    if (nrow(dat)==0) return(datatable(data.frame(Information="Aucun résultat RNA-seq."), rownames=FALSE, options=list(dom="t")))
    datatable(dat[order(FDR), .(StudiesCond, signed_logFC, abslogFC, FDR, State)], rownames=FALSE, filter="top", options=list(pageLength=13))
  })

  # =======================================================
  # INTÉGRATION
  # =======================================================

  integration_data <- reactive({

    build_binding_expression_table(
      binding_dt = binding_dt,
      rnaseq_dt = rnaseq_dt,
      protein = input$integration_protein,
      study = input$integration_study
    )
  })


  targeted_enrichment <- reactive({

    compute_fisher_enrichment(
      binding_dt = binding_dt,
      rnaseq_dt = rnaseq_dt,
      protein = input$integration_protein,
      study = input$integration_study,
      fdr_threshold = as.numeric(
        input$FDRthresh
      )
    )
  })


  enrichment_grid <- reactive({

    proteins <- input$integration_proteins_grid

    if (
      is.null(proteins) ||
      length(proteins) == 0
    ) {
      proteins <- default_proteins
    }

    compute_enrichment_grid(
      binding_dt = binding_dt,
      rnaseq_dt = rnaseq_dt,
      proteins = proteins,
      studies = all_studies,
      fdr_threshold = as.numeric(
        input$FDRthresh
      )
    )
  })


  output$integration_summary <- renderText({

    enr <- targeted_enrichment()

    if (
      nrow(enr) == 0
    ) {
      return(
        "Aucune donnée intégrée disponible."
      )
    }

    wilcox_res <- compute_binding_wilcoxon(
      binding_dt, rnaseq_dt, input$integration_protein, input$integration_study
    )

    paste0(
      "Protéine : ",
      input$integration_protein,
      "\n",
      "Étude : ",
      input$integration_study,
      "\n",
      "Gènes communs utilisés : ",
      enr$N_common,
      "\n",
      "Liés + DE : ",
      enr$Bound_DE,
      "\n",
      "Liés + non-DE : ",
      enr$Bound_nonDE,
      "\n",
      "Non liés + DE : ",
      enr$Unbound_DE,
      "\n",
      "Non liés + non-DE : ",
      enr$Unbound_nonDE,
      "\n",
      "Odds ratio Fisher : ",
      ifelse(
        is.na(enr$OddsRatio),
        "NA",
        signif(
          enr$OddsRatio,
          4
        )
      ),
      "\n",
      "p-value Fisher : ",
      ifelse(is.na(enr$Pvalue), "NA", format(enr$Pvalue, scientific=TRUE, digits=4)),
      "\nMédiane log2FC — liés : ", signif(wilcox_res$median_bound,4),
      "\nMédiane log2FC — non liés : ", signif(wilcox_res$median_unbound,4),
      "\nWilcoxon p-value : ", ifelse(is.na(wilcox_res$pvalue), "NA", format(wilcox_res$pvalue, scientific=TRUE, digits=4))
    )
  })

  output$directional_enrichment_table <- renderTable({
    dat <- compute_directional_enrichment(binding_dt, rnaseq_dt, input$integration_protein, input$integration_study, as.numeric(input$FDRthresh))
    if (nrow(dat)==0) return(NULL)
    data.frame(Direction=toupper(dat$Direction), Odds_ratio=signif(dat$OddsRatio,4), P_value=signif(dat$Pvalue,4), P_BH=signif(dat$P_adj_BH,4), check.names=FALSE)
  })


  integration_plot_object <- reactive({

    plot_type <- input$integration_plot_type
    threshold <- as.numeric(
      input$FDRthresh
    )

    if (plot_type == "enrichment_grid") {

      dat <- enrichment_grid()

      validate(
        need(
          nrow(dat) > 0,
          "Aucune donnée d'enrichissement."
        )
      )

      dat[
        ,
        log2OR_plot := pmax(
          pmin(
            log2OR,
            4
          ),
          -4
        )
      ]

      p <- ggplot(
        dat,
        aes(
          x = Protein,
          y = Study,
          size = negLog10P,
          colour = log2OR_plot
        )
      ) +
        geom_point(
          alpha = 0.85
        ) +
        scale_colour_gradient2(
          low = "blue",
          mid = "white",
          high = "red",
          midpoint = 0
        ) +
        scale_size_area(
          max_size = 12
        ) +
        labs(
          x = "Protéines",
          y = "Études / comparaisons",
          size = "-log10(p BH)",
          colour = "log2(OR)"
        ) +
        theme_classic(
          base_size = 11
        ) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          )
        )

      return(p)
    }


    if (plot_type == "bound_vs_unbound") {

      dat <- integration_data()

      validate(
        need(
          nrow(dat) > 0,
          "Aucune donnée commune pour cette protéine et cette étude."
        )
      )

      p <- ggplot(
        dat,
        aes(
          x = Binding_status,
          y = signed_logFC
        )
      ) +
        geom_boxplot(
          outlier.alpha = 0.15
        ) +
        geom_jitter(
          width = 0.18,
          alpha = 0.15,
          size = 0.8
        ) +
        geom_hline(
          yintercept = 0,
          linetype = "dashed"
        ) +
        labs(
          title = paste(
            input$integration_protein,
            "—",
            input$integration_study
          ),
          x = "Statut de liaison",
          y = "log2FC signé"
        ) +
        theme_classic(
          base_size = 12
        )

      return(p)
    }


    # Réseau intégré

    dat <- integration_data()

    validate(
      need(
        nrow(dat) > 0,
        "Aucune donnée commune pour le réseau intégré."
      )
    )

    dat <- dat[
      !is.na(FDR) &
        FDR < threshold &
        Bound == 1
    ]

    validate(
      need(
        nrow(dat) > 0,
        "Aucun gène significatif et lié avec le seuil actuel."
      )
    )

    setorder(
      dat,
      -abslogFC
    )

    dat <- head(
      dat,
      as.integer(
        input$integration_network_top
      )
    )

    genes <- dat$gene

    proteins <- input$integration_proteins_grid

    if (
      is.null(proteins) ||
      length(proteins) == 0
    ) {
      proteins <- input$integration_protein
    }

    net_dt <- binding_dt[
      Name %in% genes,
      c(
        "Name",
        proteins
      ),
      with = FALSE
    ]

    net_long <- melt(
      net_dt,
      id.vars = "Name",
      measure.vars = proteins,
      variable.name = "Protein",
      value.name = "Binding"
    )

    net_long <- net_long[
      Binding == 1
    ]

    validate(
      need(
        nrow(net_long) > 0,
        "Aucune liaison disponible pour les protéines choisies."
      )
    )

    gene_info <- dat[
      ,
      .(
        gene,
        abslogFC,
        State
      )
    ]

    gene_info <- unique(
      gene_info,
      by = "gene"
    )

    gene_levels <- unique(
      net_long$Name
    )

    protein_levels <- unique(
      net_long$Protein
    )

    gene_y <- setNames(
      seq_along(gene_levels),
      gene_levels
    )

    protein_y <- setNames(
      seq(
        1,
        max(
          length(gene_levels),
          length(protein_levels)
        ),
        length.out = length(protein_levels)
      ),
      protein_levels
    )

    edges <- net_long[
      ,
      .(
        x = 1,
        y = gene_y[Name],
        xend = 2,
        yend = protein_y[Protein]
      )
    ]

    gene_nodes <- data.table(
      label = gene_levels,
      x = 1,
      y = unname(
        gene_y[
          gene_levels
        ]
      )
    )

    gene_nodes <- merge(
      gene_nodes,
      gene_info,
      by.x = "label",
      by.y = "gene",
      all.x = TRUE
    )

    protein_nodes <- data.table(
      label = protein_levels,
      x = 2,
      y = unname(
        protein_y[
          protein_levels
        ]
      )
    )

    p <- ggplot() +
      geom_segment(
        data = edges,
        aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend
        ),
        alpha = 0.3
      ) +
      geom_point(
        data = gene_nodes,
        aes(
          x = x,
          y = y,
          size = abslogFC,
          colour = State
        )
      ) +
      geom_point(
        data = protein_nodes,
        aes(
          x = x,
          y = y
        ),
        shape = 15,
        size = 4
      ) +
      geom_text(
        data = gene_nodes,
        aes(
          x = x,
          y = y,
          label = label
        ),
        hjust = 1.1,
        size = 3
      ) +
      geom_text(
        data = protein_nodes,
        aes(
          x = x,
          y = y,
          label = label
        ),
        hjust = -0.1,
        size = 3
      ) +
      scale_colour_manual(
        values = c(
          "up" = input$colUp,
          "down" = input$colDown
        )
      ) +
      scale_size_area(
        max_size = 10
      ) +
      scale_x_continuous(
        limits = c(
          0.5,
          2.5
        ),
        breaks = c(
          1,
          2
        ),
        labels = c(
          "Gènes DE liés",
          "Protéines"
        )
      ) +
      labs(
        x = NULL,
        y = NULL,
        colour = "État",
        size = "|log2FC|"
      ) +
      theme_void() +
      theme(
        legend.position = "bottom"
      )

    p
  })


  output$integration_plot <- renderPlot({

    print(
      integration_plot_object()
    )
  })


  output$integration_table <- renderDT({

    dat <- integration_data()

    if (nrow(dat) == 0) {
      return(
        datatable(
          data.frame(
            Information = "Aucune donnée intégrée disponible."
          ),
          rownames = FALSE,
          options = list(
            dom = "t"
          )
        )
      )
    }

    display_data <- dat[
      ,
      .(
        Gene = gene,
        Study = StudiesCond,
        Protein = input$integration_protein,
        Binding = Binding_status,
        signed_logFC,
        abslogFC,
        FDR,
        State,
        chrCRM,
        startCRM,
        endCRM,
        TSS_chr = tssC,
        TSS_position = tssS
      )
    ]

    setorder(
      display_data,
      FDR
    )

    datatable(
      display_data,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 15,
        autoWidth = TRUE,
        scrollX = TRUE
      )
    )
  })


  output$download_integration <- downloadHandler(

    filename = function() {
      paste0(
        "MOLLI_integration_",
        input$integration_protein,
        "_",
        Sys.Date(),
        ".tsv"
      )
    },

    content = function(file) {

      dat <- integration_data()

      fwrite(
        dat[
          ,
          .(
            Gene = gene,
            Study = StudiesCond,
            Protein = input$integration_protein,
            Binding = Binding_status,
            signed_logFC,
            abslogFC,
            FDR,
            State,
            chrCRM,
            startCRM,
            endCRM,
            TSS_chr = tssC,
            TSS_position = tssS
          )
        ],
        file = file,
        sep = "\t"
      )
    }
  )


  # =======================================================
  # DATASET / QC
  # =======================================================

  output$dataset_qc_txt <- renderText({

    duplicated_binding_genes <- sum(
      duplicated(
        binding_dt$Name
      )
    )

    duplicated_rna_pairs <- rnaseq_dt[
      ,
      .N,
      by = .(
        gene,
        StudiesCond
      )
    ][
      N > 1,
      .N
    ]

    missing_binding <- sum(
      is.na(
        binding_dt
      )
    )

    missing_rna <- sum(
      is.na(
        rnaseq_dt[
          ,
          .(
            gene,
            StudiesCond,
            abslogFC,
            FDR,
            State
          )
        ]
      )
    )

    paste0(
      "INTERACTOME\n",
      "Lignes : ",
      nrow(binding_dt),
      "\n",
      "Gènes uniques : ",
      length(binding_genes),
      "\n",
      "Protéines : ",
      length(protein_cols),
      "\n",
      "Gènes dupliqués : ",
      duplicated_binding_genes,
      "\n",
      "Valeurs manquantes : ",
      missing_binding,
      "\n\n",

      "TRANSCRIPTOME\n",
      "Lignes : ",
      nrow(rnaseq_dt),
      "\n",
      "Gènes uniques : ",
      length(expression_genes),
      "\n",
      "Comparaisons : ",
      length(all_studies),
      "\n",
      "Couples gène × comparaison dupliqués : ",
      duplicated_rna_pairs,
      "\n",
      "Valeurs manquantes (colonnes principales) : ",
      missing_rna,
      "\n\n",

      "INTÉGRATION\n",
      "Gènes communs : ",
      length(common_genes),
      "\n",
      "Uniquement interactome : ",
      length(binding_only_genes),
      "\n",
      "Uniquement transcriptome : ",
      length(expression_only_genes)
    )
  })


  output$dataset_overlap_plot <- renderPlot({

    dat <- data.table(
      Category = factor(
        c(
          "Communs",
          "Interactome uniquement",
          "Transcriptome uniquement"
        ),
        levels = c(
          "Communs",
          "Interactome uniquement",
          "Transcriptome uniquement"
        )
      ),
      N = c(
        length(common_genes),
        length(binding_only_genes),
        length(expression_only_genes)
      )
    )

    ggplot(
      dat,
      aes(
        x = Category,
        y = N
      )
    ) +
      geom_col() +
      labs(
        x = NULL,
        y = "Nombre de gènes"
      ) +
      theme_classic(
        base_size = 12
      ) +
      theme(
        axis.text.x = element_text(
          angle = 20,
          hjust = 1
        )
      )
  })


  output$protein_frequency_qc <- renderPlot({

    freq <- data.table(
      Protein = protein_cols,
      Bound_genes = sapply(
        protein_cols,
        function(protein) {
          sum(
            binding_dt[[protein]] == 1,
            na.rm = TRUE
          )
        }
      )
    )

    setorder(
      freq,
      Bound_genes
    )

    ggplot(
      freq,
      aes(
        x = reorder(
          Protein,
          Bound_genes
        ),
        y = Bound_genes
      )
    ) +
      geom_col() +
      coord_flip() +
      labs(
        x = "Protéines",
        y = "Nombre de gènes associés"
      ) +
      theme_classic(
        base_size = 10
      )
  })
}


# =========================================================
# Lancement
# =========================================================

shinyApp(
  ui = ui,
  server = server
)
