# app.R — Two-stage Shiny app
# Stage A: Welcome/Input + filtering (select methods)
# Stage B: Results UI with dynamic tabs (GO/KEGG/Aggrescan/FoldDisco)
#          FoldDisco shows TWO tables: solubility_up.csv and solubility_down.csv
#          Optional filtering: blank input shows ALL rows; no matches → fallback to ALL rows

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(shinyjs)
  library(DT)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ragg)
  library(tibble)
  # library(shinyBS) # uncomment if you want BS tooltips
})

# ----------------- options -----------------
options(shiny.useragg = TRUE)   # crisp ragg plots
Sys.setenv(R_GTK2_MODULES = "")
options(bitmapType = "cairo")

# ----------------- helpers -----------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

load_csv <- function(fname) {
  fpath_data <- file.path("data", fname)
  if (file.exists(fpath_data)) {
    readr::read_csv(fpath_data, show_col_types = FALSE)
  } else if (file.exists(fname)) {
    readr::read_csv(fname, show_col_types = FALSE)
  } else {
    NULL
  }
}

last_col_as_id <- function(df) {
  if (is.null(df)) return(NULL)
  idcol <- names(df)[ncol(df)]
  df$Protein_ID <- df[[idcol]]
  df
}

first_id_col <- function(df) if (!is.null(df)) names(df)[1] else NULL

find_file <- function(fname) {
  fpath_data <- file.path("data", fname)
  if (file.exists(fpath_data)) fpath_data else fname
}

page_css <- "
.page-transition { animation: fadeIn 0.5s; }
@keyframes fadeIn { from { opacity: 0; } to { opacity: 1; } }
"

# ----------------- GLOBAL CSV fallbacks (optional) -----------------
ref_df_fallback      <- last_col_as_id(load_csv("data/reference_cleaned.csv"))
up_df_fallback       <- last_col_as_id(load_csv("data/up_cleaned.csv"))
down_df_fallback     <- last_col_as_id(load_csv("data/down_cleaned.csv"))
go_df_fallback       <- last_col_as_id(load_csv("data/GO_results.csv"))
kegg_leg_df_fallback <- last_col_as_id(load_csv("data/KEGG_Legacy_results.csv"))
kegg_med_df_fallback <- last_col_as_id(load_csv("data/KEGG_Medicus_results.csv"))
agg_df_fallback      <- last_col_as_id(load_csv("data/demo_proteins_waggrescan.csv"))
fold_df_fallback     <- last_col_as_id(load_csv("data/folddisco_cleaned.csv"))

# Split solubility CSVs used on FoldDisco tab
solu_up_fallback     <- last_col_as_id(load_csv("data/solubility_up.csv"))
solu_down_fallback   <- last_col_as_id(load_csv("data/solubility_down.csv"))

# =============================================================================
# MODULE A — Welcome/Input + Filtering
# =============================================================================

method_descriptions <- list(
  GO = "Gene Ontology enrichment identifies biological processes, molecular functions, and cellular components overrepresented in your protein list. Helps understand what biological roles these proteins play.",
  KEGG = "KEGG pathway analysis maps proteins to known biological pathways and molecular interaction networks. Reveals which cellular processes and signaling cascades are affected.",
  Aggrescan = "Aggrescan3D v2 predicts protein aggregation propensity based on 3D structure. High scores indicate regions prone to forming harmful clumps, relevant in aging and neurodegeneration.",
  FoldDisco = "FoldDisco identifies structural motifs and similarities between proteins using 3D structure comparison. Discovers patterns not visible in sequence analysis alone."
)

appA_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    useShinyjs(),
    theme = bs_theme(
      version = 5,
      bootswatch = "flatly",
      base_font = font_google("Poppins")
    ),
    tags$head(
      tags$style(HTML(page_css)),
      tags$link(rel = "stylesheet",
                href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css"),
      tags$style(HTML("
        .info-icon { color: #3498db; cursor: help; margin-left: 8px; font-size: 18px; }
        .method-checkbox { margin: 15px 0; }
        .welcome-container { text-align: center; padding: 50px 20px; }
        .welcome-title { font-size: 3em; font-weight: 600; margin-bottom: 20px; color: #2c3e50; }
        .welcome-subtitle { font-size: 1.3em; color: #7f8c8d; margin-bottom: 40px; }
        .protein-image { max-width: 500px; margin: 30px auto; border-radius: 15px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
        .start-btn { font-size: 1.2em; padding: 15px 50px; margin-top: 30px; }
      "))
    ),
    
    # Welcome
    div(id = ns("welcome_page"), class = "page-transition",
        div(class = "welcome-container",
            h1(class = "welcome-title", "🧬 Protein Analysis Dashboard"),
            p(class = "welcome-subtitle", "Comprehensive proteomics analysis for your BioHackathon project"),
            img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/5/5f/Myoglobin.png/400px-Myoglobin.png",
                class = "protein-image", alt = "Protein Structure"),
            p("Analyze protein motifs, aggregation propensity, pathways, and structural similarities",
              style = "font-size: 1.1em; color: #34495e; margin: 20px 0;"),
            actionButton(ns("start_btn"), "Start Analysis", class = "btn-primary start-btn")
        )
    ),
    
    # Input page
    hidden(
      div(id = ns("input_page"), class = "page-transition",
          titlePanel("Protein Analysis Input"),
          fluidRow(
            column(8,
                   wellPanel(
                     h4("Step 1: Enter Protein IDs"),
                     textAreaInput(ns("protein_ids"),
                                   label = "Paste protein IDs (one per line or comma-separated):",
                                   placeholder = "e.g.,\nP12345\nQ98765\nO43210",
                                   rows = 8, width = "100%"),
                     p("Or upload a file:", style = "margin-top: 15px;"),
                     fileInput(ns("protein_file"), label = NULL,
                               accept = c(".txt", ".csv"),
                               buttonLabel = "Browse...",
                               placeholder = "No file chosen")
                   ),
                   wellPanel(
                     h4("Step 2: Select Analysis Methods"),
                     div(class = "method-checkbox",
                         checkboxInput(ns("method_go"), label = NULL, value = TRUE),
                         tags$span("Gene Ontology (GO) Enrichment"),
                         tags$i(class = "fas fa-info-circle info-icon", id = ns("info_go"))
                     ),
                     div(class = "method-checkbox",
                         checkboxInput(ns("method_kegg"), label = NULL, value = TRUE),
                         tags$span("KEGG Pathway Analysis"),
                         tags$i(class = "fas fa-info-circle info-icon", id = ns("info_kegg"))
                     ),
                     div(class = "method-checkbox",
                         checkboxInput(ns("method_aggrescan"), label = NULL, value = TRUE),
                         tags$span("Aggrescan3D v2 (Aggregation Prediction)"),
                         tags$i(class = "fas fa-info-circle info-icon", id = ns("info_aggrescan"))
                     ),
                     div(class = "method-checkbox",
                         checkboxInput(ns("method_folddisco"), label = NULL, value = TRUE),
                         tags$span("FoldDisco (Structural Motifs)"),
                         tags$i(class = "fas fa-info-circle info-icon", id = ns("info_folddisco"))
                     )
                   ),
                   actionButton(ns("submit_btn"), "Submit Analysis",
                                class = "btn-success btn-lg",
                                style = "width: 100%; margin-top: 20px;")
            ),
            column(4,
                   wellPanel(
                     h4("ℹ️ Instructions"),
                     tags$ol(
                       tags$li("Enter your protein IDs in the text area or upload a file"),
                       tags$li("Select which analysis methods to run"),
                       tags$li("Click Submit to view results"),
                       tags$li("Hover over the ℹ️ icons for method descriptions")
                     ),
                     hr(),
                     h5("Accepted ID formats:"),
                     tags$ul(
                       tags$li("UniProt IDs (e.g., P12345)"),
                       tags$li("Gene symbols (e.g., TP53)"),
                       tags$li("One per line or comma-separated")
                     )
                   )
            )
          )
      )
    )
  )
}

appA_server <- function(id, results_rv, go_to_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$start_btn, { hide("welcome_page"); show("input_page") })
    
    get_protein_ids <- reactive({
      ids <- character()
      if (!is.null(input$protein_ids) && input$protein_ids != "") {
        ids <- c(ids, unlist(strsplit(input$protein_ids, "[,\\n\\s]+")))
      }
      if (!is.null(input$protein_file)) {
        file_content <- readLines(input$protein_file$datapath, warn = FALSE)
        ids <- c(ids, unlist(strsplit(paste(file_content, collapse = " "), "[,\\n\\s]+")))
      }
      ids <- unique(trimws(ids[ids != ""]))
      ids  # empty character() if nothing entered (optional filtering)
    })
    
    build_results <- eventReactive(input$submit_btn, {
      protein_ids <- get_protein_ids()
      
      # Optional filtering logic:
      # - If no IDs: return full df (no filter)
      # - If IDs given but 0 matches: return full df (fallback)
      filter_by_ids <- function(df) {
        if (is.null(df)) return(NULL)
        if (length(protein_ids) == 0) return(df)  # blank input → show all
        idcol <- names(df)[ncol(df)]              # last_col_as_id adds Protein_ID last
        out <- dplyr::filter(df, .data[[idcol]] %in% protein_ids)
        if (nrow(out) == 0) df else out           # no matches → show all
      }
      
      results <- list(
        protein_ids = protein_ids,
        # general up/down/ref (if present)
        ref       = filter_by_ids(ref_df_fallback),
        up        = filter_by_ids(up_df_fallback),
        down      = filter_by_ids(down_df_fallback),
        # per-method datasets (if present)
        go        = filter_by_ids(go_df_fallback),
        kegg_leg  = filter_by_ids(kegg_leg_df_fallback),
        kegg_med  = filter_by_ids(kegg_med_df_fallback),
        aggrescan = filter_by_ids(agg_df_fallback),
        folddisco = filter_by_ids(fold_df_fallback),
        
        # solubility tables for FoldDisco tab
        solubility_up   = filter_by_ids(solu_up_fallback),
        solubility_down = filter_by_ids(solu_down_fallback)
      )
      
      results$methods <- c(
        if (isTRUE(input$method_go)) "GO" else NULL,
        if (isTRUE(input$method_kegg)) "KEGG" else NULL,
        if (isTRUE(input$method_aggrescan)) "Aggrescan" else NULL,
        if (isTRUE(input$method_folddisco)) "FoldDisco" else NULL
      )
      
      results
    })
    
    observeEvent(build_results(), {
      results_rv$data <- build_results()
      go_to_results()
    })
  })
}

# =============================================================================
# MODULE B — Results (dynamic tabs + split FoldDisco tables)
# =============================================================================

resultsB_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme = bs_theme(version = 5, bootswatch = "flatly"),
    useShinyjs(),
    titlePanel("Protein Analysis Dashboard"),
    actionButton(ns("back_btn"), "← Back to Input", class = "btn-secondary",
                 style = "margin-bottom: 16px;"),
    uiOutput(ns("tabs_dyn"))  # tabs are built in server based on toggles
  )
}

serverB <- function(id, results_rv, back_to_input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$back_btn, back_to_input())
    
    # Prefer Stage-A data; fall back to CSVs
    ref_df   <- reactive({ results_rv$data$ref       %||% ref_df_fallback })
    up_df    <- reactive({ results_rv$data$up        %||% up_df_fallback })
    down_df  <- reactive({ results_rv$data$down      %||% down_df_fallback })
    go_df    <- reactive({ results_rv$data$go        %||% go_df_fallback })
    kleg_df  <- reactive({ results_rv$data$kegg_leg  %||% kegg_leg_df_fallback })
    kmed_df  <- reactive({ results_rv$data$kegg_med  %||% kegg_med_df_fallback })
    agg_df   <- reactive({ results_rv$data$aggrescan %||% agg_df_fallback })
    fold_df  <- reactive({ results_rv$data$folddisco %||% fold_df_fallback })
    
    # Split solubility tables (robust fallback if somehow empty)
    solu_up_df <- reactive({
      df <- results_rv$data$solubility_up
      if (is.null(df) || nrow(df) == 0) solu_up_fallback else df
    })
    solu_dn_df <- reactive({
      df <- results_rv$data$solubility_down
      if (is.null(df) || nrow(df) == 0) solu_down_fallback else df
    })
    
    comb_df <- reactive({
      u <- up_df(); d <- down_df()
      if (is.null(u) && is.null(d)) return(NULL)
      u2 <- if (!is.null(u)) mutate(u, Direction = "Up") else NULL
      d2 <- if (!is.null(d)) mutate(d, Direction = "Down") else NULL
      bind_rows(u2, d2)
    })
    
    # ----------------- General tab outputs -----------------
    counts_df <- reactive({
      tibble(
        Set = c("Reference","Up","Down"),
        Count = c(nrow(ref_df() %||% tibble()),
                  nrow(up_df()  %||% tibble()),
                  nrow(down_df()%||% tibble()))
      )
    })
    
    output$counts_tbl <- renderDT({
      datatable(counts_df(), options = list(dom="t"), rownames = FALSE)
    })
    
    output$counts_plot <- renderImage({
      tmp <- tempfile(fileext = ".png")
      ragg::agg_png(tmp, width = 900, height = 420, res = 120); on.exit(dev.off())
      print(
        ggplot(counts_df(), aes(x = Set, y = Count, fill = Set)) +
          geom_col(width=0.6) +
          geom_text(aes(label = Count), vjust=-0.3, size=5) +
          scale_fill_manual(values = c(Reference="#64748b", Up="#ef4444", Down="#3b82f6")) +
          theme_minimal(base_size = 14) +
          labs(title="Up/Down/Reference counts", y="Count")
      )
      list(src=tmp, contentType="image/png", width="100%")
    }, deleteFile=TRUE)

    
    # ----------------- GO -----------------
    output$go_tbl <- renderDT({
      datatable(go_df() %||% tibble(), options=list(pageLength=10, scrollX=TRUE))
    })
    
    output$go_bar <- renderImage({
      df <- go_df() %||% tibble()
      validate(need(nrow(df) > 0, "No GO data."))
      idcol <- if ("ontology_id" %in% names(df)) "ontology_id" else first_id_col(df)
      tall <- df %>% count(.data[[idcol]], name="Count") %>% arrange(desc(Count)) %>% head(20)
      tmp <- tempfile(fileext = ".png")
      ragg::agg_png(tmp, width = 1000, height = 520, res = 120); on.exit(dev.off())
      print(
        ggplot(tall, aes(x=reorder(.data[[idcol]], Count), y=Count)) +
          geom_col(fill="#0ea5e9") +
          geom_text(aes(label=Count), hjust=-0.2, size=4) +
          coord_flip() +
          theme_minimal(base_size=14) +
          labs(title="Top GO Terms", x="GO Term", y="Count")
      )
      list(src=tmp, contentType="image/png", width="100%")
    }, deleteFile=TRUE)
    
    output$go_heat <- renderImage({
      df <- go_df() %||% tibble()
      validate(need(nrow(df) > 0, "No GO data."))
      if (!"p_value" %in% names(df)) df$p_value <- NA_real_
      df$p_value <- suppressWarnings(as.numeric(df$p_value))
      xcol <- if ("element_id_in_ontology" %in% names(df)) "element_id_in_ontology" else first_id_col(df)
      ycol <- if ("ontology_id" %in% names(df)) "ontology_id" else names(df)[2 %||% 1]
      tmp <- tempfile(fileext = ".png")
      ragg::agg_png(tmp, width=1000, height=520, res=120); on.exit(dev.off())
      print(
        ggplot(df, aes(x=.data[[xcol]], y=.data[[ycol]], fill=p_value)) +
          geom_tile() +
          scale_fill_gradient(low="red", high="gray", name="p_value", na.value="gray90") +
          theme_minimal(base_size=12) +
          theme(axis.text.x=element_text(angle=60,hjust=1)) +
          labs(title="GO Heatmap", x="Element", y="GO Term")
      )
      list(src=tmp, contentType="image/png", width="100%")
    }, deleteFile=TRUE)
    
    output$go_dl <- downloadHandler(
      filename = function() "GO_results.csv",
      content = function(file) {
        df <- go_df()
        if (!is.null(df)) write.csv(df, file, row.names = FALSE) else file.copy(find_file("GO_results.csv"), file)
      }
    )
    
    # ----------------- KEGG Legacy -----------------
    output$kleg_tbl <- renderDT({
      datatable(kleg_df() %||% tibble(), options=list(pageLength=10, scrollX=TRUE))
    })
    
    output$kleg_bar <- renderImage({
      df <- kleg_df() %||% tibble(); validate(need(nrow(df) > 0, "No KEGG Legacy data."))
      idcol <- first_id_col(df)
      tall <- df %>% count(.data[[idcol]], name="Count") %>% arrange(desc(Count)) %>% head(20)
      tmp <- tempfile(fileext = ".png")
      ragg::agg_png(tmp, width = 1000, height = 520, res = 120); on.exit(dev.off())
      print(
        ggplot(tall, aes(x=reorder(.data[[idcol]], Count), y=Count)) +
          geom_col(fill="#22c55e") +
          geom_text(aes(label=Count), hjust=-0.2, size=4) +
          coord_flip() +
          theme_minimal(base_size=14) +
          labs(title="Top KEGG Legacy Pathways", x="Pathway", y="Count")
      )
      list(src=tmp, contentType="image/png", width="100%")
    }, deleteFile=TRUE)
    
    output$kleg_lolli <- renderImage({
      df <- kleg_df() %>% mutate(p_value = suppressWarnings(as.numeric(.data[["p_value"]]))) %||% tibble()
      validate(need(nrow(df) > 0, "No KEGG Legacy data."))
      idcol <- first_id_col(df)
      if (!"p_value" %in% names(df) || all(is.na(df$p_value))) {
        df <- df %>% count(.data[[idcol]], name="Count")
        df$neglog10P <- -log10((max(df$Count, na.rm=TRUE) - df$Count + 1) / (max(df$Count, na.rm=TRUE) + 1))
      } else {
        df$neglog10P <- -log10(df$p_value)
      }
      top <- df %>% arrange(desc(neglog10P)) %>% head(20)
      tmp <- tempfile(fileext = ".png")
      ragg::agg_png(tmp, width = 1000, height = 520, res = 120); on.exit(dev.off())
      print(
        ggplot(top, aes(x=neglog10P, y=reorder(.data[[idcol]], neglog10P))) +
          geom_segment(aes(x=0, xend=neglog10P, y=reorder(.data[[idcol]], neglog10P), yend=reorder(.data[[idcol]], neglog10P)), linewidth=1) +
          geom_point(size=3) +
          theme_minimal(base_size=14) +
          labs(title="KEGG Legacy Lollipop", x="-log10(p)", y="Pathway")
      )
      list(src=tmp, contentType="image/png", width="100%")
    }, deleteFile=TRUE)
    
    output$kleg_dl <- downloadHandler(
      filename = function() "KEGG_Legacy_results.csv",
      content = function(file){
        df <- kleg_df()
        if (!is.null(df)) write.csv(df, file, row.names = FALSE) else file.copy(find_file("KEGG_Legacy_results.csv"), file)
      })
    
    # ----------------- KEGG Medicus -----------------
    output$kmed_tbl <- renderDT({
      datatable(kmed_df() %||% tibble(), options=list(pageLength=10, scrollX=TRUE))
    })
    
    output$kmed_bar <- renderImage({
      df <- kmed_df() %||% tibble(); validate(need(nrow(df) > 0, "No KEGG Medicus data."))
      idcol <- first_id_col(df)
      tall <- df %>% count(.data[[idcol]], name="Count") %>% arrange(desc(Count)) %>% head(20)
      tmp <- tempfile(fileext = ".png")
      ragg::agg_png(tmp, width = 1000, height = 520, res = 120); on.exit(dev.off())
      print(
        ggplot(tall, aes(x=reorder(.data[[idcol]], Count), y=Count)) +
          geom_col(fill="#14b8a6") +
          geom_text(aes(label=Count), hjust=-0.2, size=4) +
          coord_flip() +
          theme_minimal(base_size=14) +
          labs(title="Top KEGG Medicus Pathways", x="Pathway", y="Count")
      )
      list(src=tmp, contentType="image/png", width="100%")
    }, deleteFile=TRUE)
    
    output$kmed_lolli <- renderImage({
      df <- kmed_df() %>% mutate(p_value = suppressWarnings(as.numeric(.data[["p_value"]]))) %||% tibble()
      validate(need(nrow(df) > 0, "No KEGG Medicus data."))
      idcol <- first_id_col(df)
      if (!"p_value" %in% names(df) || all(is.na(df$p_value))) {
        df <- df %>% count(.data[[idcol]], name="Count")
        df$neglog10P <- -log10((max(df$Count, na.rm=TRUE) - df$Count + 1) / (max(df$Count, na.rm=TRUE) + 1))
      } else {
        df$neglog10P <- -log10(df$p_value)
      }
      top <- df %>% arrange(desc(neglog10P)) %>% head(20)
      tmp <- tempfile(fileext = ".png")
      ragg::agg_png(tmp, width = 1000, height = 520, res = 120); on.exit(dev.off())
      print(
        ggplot(top, aes(x=neglog10P, y=reorder(.data[[idcol]], neglog10P))) +
          geom_segment(aes(x=0, xend=neglog10P, y=reorder(.data[[idcol]], neglog10P), yend=reorder(.data[[idcol]], neglog10P)), linewidth=1) +
          geom_point(size=3) +
          theme_minimal(base_size=14) +
          labs(title="KEGG Medicus Lollipop", x="-log10(p)", y="Pathway")
      )
      list(src=tmp, contentType="image/png", width="100%")
    }, deleteFile=TRUE)
    
    output$kmed_dl <- downloadHandler(
      filename = function() "KEGG_Medicus_results.csv",
      content = function(file){
        df <- kmed_df()
        if (!is.null(df)) write.csv(df, file, row.names = FALSE) else file.copy(find_file("KEGG_Medicus_results.csv"), file)
      })
    
    # ----------------- Aggrescan -----------------
    output$agg_tbl <- renderDT({
      datatable(agg_df() %||% tibble(), options=list(pageLength=10, scrollX=TRUE))
    })
    
    # >>> CHANGED: tilt x-axis labels 45° for readability
    output$agg_volcano <- renderImage({
      df <- agg_df() %||% tibble()
      validate(need(nrow(df) > 0, "No Aggrescan data."))
      xcol <- if ("log2FC" %in% names(df)) "log2FC" else names(df)[2 %||% 1]
      ycol <- if ("neglog10P" %in% names(df)) "neglog10P" else names(df)[3 %||% 1]
      tmp <- tempfile(fileext = ".png")
      ragg::agg_png(tmp, width = 1000, height = 520, res = 120); on.exit(dev.off())
      print(
        ggplot(df, aes(x=.data[[xcol]], y=.data[[ycol]])) +
          geom_point(alpha=0.7, size=1.8, color="#ef4444", na.rm = TRUE) +
          theme_minimal(base_size=14) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +   # <-- new
          labs(title="Aggrescan Divergence Volcano", x=xcol, y=ycol)
      )
      list(src=tmp, contentType="image/png", width="100%")
    }, deleteFile=TRUE)
    # <<< END CHANGE
    
    output$agg_dl <- downloadHandler(
      filename = function() "demo_proteins_waggrescan.csv",
      content = function(file){
        df <- agg_df()
        if (!is.null(df)) write.csv(df, file, row.names = FALSE) else file.copy(find_file("demo_proteins_waggrescan.csv"), file)
      })
    
    # ----------------- FoldDisco (two tables + bottom downloads) -----------------
    output$fold_up_tbl <- renderDT({
      datatable(solu_up_df() %||% tibble(),
                options = list(pageLength = 10, scrollX = TRUE))
    })
    output$fold_dn_tbl <- renderDT({
      datatable(solu_dn_df() %||% tibble(),
                options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$fold_up_dl <- downloadHandler(
      filename = function() "solubility_up.csv",
      content = function(file) {
        df <- solu_up_df()
        if (!is.null(df)) write.csv(df, file, row.names = FALSE) else file.copy(find_file("solubility_up.csv"), file)
      }
    )
    output$fold_dn_dl <- downloadHandler(
      filename = function() "solubility_down.csv",
      content = function(file) {
        df <- solu_dn_df()
        if (!is.null(df)) write.csv(df, file, row.names = FALSE) else file.copy(find_file("solubility_down.csv"), file)
      }
    )
    
    # ----------------- Dynamic tabs (methods toggles) -----------------
    output$tabs_dyn <- renderUI({
      req(results_rv$data)
      methods <- results_rv$data$methods %||% character()
      
      tabs <- list(
        tabPanel("General",
                 br(),
                 h4("Dataset counts"),
                 DTOutput(ns("counts_tbl")),
                 br(),
                 h4("Up vs Down Bar Plot"),
                 imageOutput(ns("counts_plot"), height = "380px"),
                 br(),
                
        )
      )
      
      if ("GO" %in% methods) {
        tabs <- c(tabs, list(
          tabPanel("GO",
                   br(),
                   DTOutput(ns("go_tbl")),
                   br(),
                   h4("Top GO Terms"),
                   imageOutput(ns("go_bar"), height = "520px"),
                   br(),
                   h4("GO Heatmap"),
                   imageOutput(ns("go_heat"), height = "520px"),
                   tags$hr(),
                   div(style="text-align:right;",
                       downloadButton(ns("go_dl"), "Download GO CSV")
                   )
          )
        ))
      }
      
      if ("KEGG" %in% methods) {
        tabs <- c(tabs, list(
          tabPanel("KEGG Legacy",
                   br(),
                   DTOutput(ns("kleg_tbl")),
                   br(),
                   h4("Top KEGG Legacy Pathways"),
                   imageOutput(ns("kleg_bar"), height = "520px"),
                   br(),
                   h4("Lollipop (KEGG Legacy)"),
                   imageOutput(ns("kleg_lolli"), height = "520px"),
                   tags$hr(),
                   div(style="text-align:right;",
                       downloadButton(ns("kleg_dl"), "Download KEGG Legacy CSV")
                   )
          ),
          tabPanel("KEGG Medicus",
                   br(),
                   DTOutput(ns("kmed_tbl")),
                   br(),
                   h4("Top KEGG Medicus Pathways"),
                   imageOutput(ns("kmed_bar"), height = "520px"),
                   br(),
                   h4("Lollipop (KEGG Medicus)"),
                   imageOutput(ns("kmed_lolli"), height = "520px"),
                   tags$hr(),
                   div(style="text-align:right;",
                       downloadButton(ns("kmed_dl"), "Download KEGG Medicus CSV")
                   )
          )
        ))
      }
      
      if ("Aggrescan" %in% methods) {
        tabs <- c(tabs, list(
          tabPanel("Aggrescan",
                   br(),
                   DTOutput(ns("agg_tbl")),
                   br(),
                   h4("Aggrescan Divergence Volcano"),
                   imageOutput(ns("agg_volcano"), height = "520px"),
                   tags$hr(),
                   div(style="text-align:right;",
                       downloadButton(ns("agg_dl"), "Download Aggrescan CSV")
                   )
          )
        ))
      }
      
      if ("FoldDisco" %in% methods) {
        tabs <- c(tabs, list(
          tabPanel("FoldDisco",
                   br(),
                   fluidRow(
                     column(
                       12,
                       h4("Solubility ↑ (solubility_up.csv)"),
                       DTOutput(ns("fold_up_tbl"))
                     )
                   ),
                   br(),
                   fluidRow(
                     column(
                       12,
                       h4("Solubility ↓ (solubility_down.csv)"),
                       DTOutput(ns("fold_dn_tbl"))
                     )
                   ),
                   tags$hr(),
                   # Downloads at the bottom
                   fluidRow(
                     column(
                       6,
                       div(style="text-align:left;",
                           downloadButton(ns("fold_up_dl"), "Download solubility_up.csv")
                       )
                     ),
                     column(
                       6,
                       div(style="text-align:right;",
                           downloadButton(ns("fold_dn_dl"), "Download solubility_down.csv")
                       )
                     )
                   )
          )
        ))
      }
      
      do.call(tabsetPanel, tabs)
    })
  })
}

# =============================================================================
# APP SHELL — swap between A (inputs) and B (results)
# =============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML(page_css))),
  uiOutput("stage")
)

server <- function(input, output, session) {
  stage <- reactiveVal("A")                 # "A" → "B"
  results_rv <- reactiveValues(data = NULL)
  
  go_to_results <- function() stage("B")
  back_to_input <- function()  stage("A")
  
  output$stage <- renderUI({
    if (stage() == "A") {
      appA_ui("A")
    } else {
      resultsB_ui("B")
    }
  })
  
  appA_server("A", results_rv = results_rv, go_to_results = go_to_results)
  serverB("B",    results_rv = results_rv, back_to_input = back_to_input)
}

shinyApp(ui, server)
