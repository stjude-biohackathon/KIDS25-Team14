# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a biohackathon project (KIDS25-Team14) focused on analyzing protein solubility changes in aging mouse skeletal muscle. The project provides a comprehensive protein enrichment and aggregation analysis toolkit through an interactive Shiny web application.

**Biological Context**: Aging leads to reduced protein solubility in mammalian tissues due to oxidative damage and post-translational modifications. This project analyzes proteins with reduced solubility in aged mouse skeletal muscle to identify common functions, pathways, structural motifs, and aggregation propensities.

## Architecture

### Two-Stage Shiny Application (UI_main.R)

The main application is structured as a modular two-stage Shiny app:

1. **Stage A (Input Module)**:
   - Welcome page with start button
   - Protein ID input (text area or file upload)
   - Analysis method selection (GO, KEGG, Aggrescan, FoldDisco)
   - Uses modular design with `appA_ui()` and `appA_server()`

2. **Stage B (Results Module)**:
   - Dynamic tab generation based on selected methods
   - Results visualization with interactive plots and data tables
   - Download handlers for each analysis type
   - Uses modular design with `resultsB_ui()` and `serverB()`

**Key Design Patterns**:
- Module-based architecture with namespaced IDs using `NS()`
- Reactive values (`results_rv`) to pass data between modules
- Callback functions (`go_to_results`, `back_to_input`) for stage navigation
- Fallback data loading from pre-computed CSV files in `data/` directory
- Optional filtering: blank input shows all rows, no matches falls back to all rows

### Analysis Pipeline Components

1. **Enrichment Analysis (EnrichmentAnalysis/Enrichment.Rmd)**:
   - R Markdown workflow using the `mulea` package
   - GO enrichment analysis using GMT files
   - KEGG pathway analysis (Legacy and Medicus) via `msigdbr`
   - GSEA implementation with customizable parameters
   - Outputs CSV results and visualization figures

2. **Protein Data Acquisition (notebooks/fetch_protein_data.ipynb)**:
   - Python-based data extraction from Excel sheets
   - UniProt API integration for sequence data (FASTA format)
   - AlphaFold structure file retrieval from reference proteome
   - Aggrescan3D REST API integration for aggregation scoring
   - Batch processing with error logging

3. **Structural Motif Analysis (notebooks/Folddisco_Search.ipynb)**:
   - FoldDisco index querying for structural similarity
   - Processes multiple query sets (solubility_up, solubility_down)
   - Parallel processing with configurable thread count
   - UniProt ID extraction and mapping
   - Top-N result filtering and summarization

4. **Data Processing (notebooks/process_data_and_viz.ipynb)**:
   - Data cleaning and transformation
   - Preparation of datasets for Shiny app consumption

## Data Directory Structure

```
data/
├── reference_cleaned.csv           # Background/reference protein set
├── up_cleaned.csv                  # Proteins with increased insolubility
├── down_cleaned.csv                # Proteins with decreased insolubility
├── GO_results.csv                  # Gene Ontology enrichment
├── KEGG_Legacy_results.csv         # KEGG Legacy pathway analysis
├── KEGG_Medicus_results.csv        # KEGG Medicus pathway analysis
├── demo_proteins_waggrescan.csv    # Aggrescan3D aggregation scores
├── folddisco_cleaned.csv           # Processed FoldDisco results
├── solubility_up.csv               # FoldDisco up-regulated subset
└── solubility_down.csv             # FoldDisco down-regulated subset
```

## Running the Application

### Start the Shiny App
```bash
# From project root
R -e "shiny::runApp('UI_main.R')"
```

The app will load fallback data from `data/*.csv` files automatically. User input filtering is optional - leaving the input blank shows all data.

### Run Enrichment Analysis
```bash
cd EnrichmentAnalysis
# Open in RStudio and knit, or run via command line:
R -e "rmarkdown::render('Enrichment.Rmd')"
```

**Required R packages**: `mulea`, `tidyverse`, `readxl`, `msigdbr`

**Input Requirements**:
- Excel file with 'up' and 'down' sheets in `EnrichmentAnalysis/Inputs/data.xlsx`
- GMT files in `EnrichmentAnalysis/Genesets/`

**Outputs**: CSV files in `EnrichmentAnalysis/Outputs/` and figures in `EnrichmentAnalysis/Figs/`

### Run Python Notebooks

All Python notebooks are in the `notebooks/` directory:

```bash
# For Jupyter notebook execution
jupyter notebook notebooks/fetch_protein_data.ipynb
jupyter notebook notebooks/Folddisco_Search.ipynb
```

**Required Python packages**: `pandas`, `requests`, `numpy`, `openpyxl`

## Key Implementation Details

### Shiny App Data Flow

1. **Data Loading**:
   - CSV files loaded at startup as fallback via `load_csv()` helper
   - Last column treated as Protein_ID via `last_col_as_id()` helper
   - Reactive values prefer user-filtered data, fall back to CSV data

2. **Filtering Logic** (in `appA_server`):
   - Extract protein IDs from text input or uploaded file
   - Apply `filter_by_ids()`: empty input returns all, zero matches returns all
   - Store filtered results in `results_rv$data`

3. **Dynamic UI Generation** (in `serverB`):
   - Tabs created based on `results_rv$data$methods` vector
   - Each analysis tab includes: data table (DT), visualizations (ggplot2), download button
   - FoldDisco tab displays two separate tables (up/down) with corresponding downloads

4. **Visualization Approach**:
   - Use `ragg::agg_png()` for high-quality rasterized plots
   - Render plots to temporary files, return via `renderImage()`
   - All plots generated with `ggplot2` and `theme_minimal()`

### Aggrescan3D Integration

The Aggrescan3D API requires:
- Uncompressed PDB files
- JSON options payload with structure file
- Polling for job completion (status endpoint)
- Result parsing for A3D scores (avg, sum, min, max)

**Important**: API response key changed from `aggrescan3Dscore` to `A3Dscore` - code handles both.

### FoldDisco Workflow

```bash
folddisco query \
  -i <index_path> \
  -p <query.pdb> \
  -t <threads> \
  --per-structure \
  --sort-by-score \
  --header \
  --skip-match
```

Output is TSV format with columns: `id`, `idf_score`, `min_rmsd`, `max_node_cov`. Extract UniProt IDs using regex patterns to match against reference lists.

### GSEA with mulea

The `mulea` package workflow:
1. Format gene sets as three-column data frame: `ontology_id`, `ontology_name`, `list_of_values`
2. Create GSEA model: `gsea(gmt, element_names, element_scores, element_score_type, number_of_permutations)`
3. Run test: `run_test(gsea_model)`
4. Filter by p-value threshold
5. Annotate results with custom function to map genes/proteins to pathways
6. Reshape for visualization: `reshape_results()`
7. Visualize with `plot_graph()`, `plot_lollipop()`, `plot_heatmap()`

## Common Modifications

### Adding New Analysis Methods

1. Add CSV data file to `data/` directory
2. Load fallback data in `UI_main.R` global section
3. Add checkbox in `appA_ui()` method selection panel
4. Add method to `results$methods` vector in `appA_server()`
5. Create new `tabPanel()` in `output$tabs_dyn` in `serverB()`
6. Implement reactive data, DT output, plots, and download handler

### Updating Protein Lists

Replace or add files to `data/` directory. The app automatically loads CSVs on startup. Ensure last column contains protein IDs for filtering to work correctly.

### Modifying Visualizations

All plots use `ggplot2`. Locate the `renderImage()` block for the target plot and modify the `ggplot()` call. Use `ragg::agg_png()` for rendering and set appropriate dimensions (width, height, res).

## Development Notes

- **Shiny reactivity**: Use `reactive()` for data transformations, `eventReactive()` for user actions, `observe()/observeEvent()` for side effects
- **Module communication**: Pass reactive values objects and callback functions as parameters to module server functions
- **Error handling**: Use `validate(need())` in render functions for graceful degradation
- **ID extraction**: Helper function `extract_uniprot_from_id()` handles multiple formats: AF-UNIPROT-*, UNIPROT_*, or standalone IDs
- **Batch processing**: When calling external APIs, implement batching and rate limiting to avoid overwhelming servers

## File Organization

- `UI_main.R` - Main Shiny application (708 lines, modular structure)
- `EnrichmentAnalysis/` - R Markdown workflow for GO/KEGG analysis
- `notebooks/` - Python Jupyter notebooks for data acquisition and processing
- `data/` - Pre-computed results and reference data (CSV format)
- `assets/` - Images for README
- `www/` - Static assets for Shiny app (if needed)

## Dependencies

**R packages**: `shiny`, `bslib`, `shinyjs`, `DT`, `readr`, `dplyr`, `tidyr`, `ggplot2`, `ggrepel`, `ragg`, `tibble`, `mulea`, `msigdbr`, `readxl`

**Python packages**: `pandas`, `numpy`, `requests`, `openpyxl`

**External tools**: FoldDisco (command-line tool with pre-built index)

**APIs**: UniProt REST API, Aggrescan3D REST API
