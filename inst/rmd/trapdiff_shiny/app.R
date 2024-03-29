library(shiny)
library(shinydashboard)
library(plotly)
# Temprarily point biomart to a mirror
# options(
#   rmyknife.verbose = TRUE,
#   rmyknife.use_biomart_mirror = TRUE,
#   rmyknife.biomart_mirror_host = "useast.ensembl.org"
# )
library(rmyknife)
# Define ggplot2 options
ggplot2::theme_set(
  ggplot2::theme_minimal(
    base_size = 12
  )
)
ggplot2::theme_update(
  axis.ticks = ggplot2::element_line(color = "grey92"),
  axis.ticks.length = ggplot2::unit(.5, "lines"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.title = ggplot2::element_text(size = 12),
  legend.text = ggplot2::element_text(color = "grey30"),
  plot.title = ggplot2::element_text(size = 18, face = "bold"),
  plot.subtitle = ggplot2::element_text(size = 12, color = "grey30"),
  plot.caption = ggplot2::element_text(size = 9, margin = ggplot2::margin(t = 15)),
  plot.title.position = 'plot'
)

fc_scatterplot <- function(
  ggplot_dat,
  colour_lab,
  x_lab = "fc_ip",
  y_lab = "fc_input",
  point_alpha = 0.3,
  min_range = -5,
  max_range = 5) {
  ggplot_dat +
    ggplot2::geom_point(
      alpha = point_alpha
    ) +
    ggplot2::scale_x_continuous(limits = c(min_range, max_range)) +
    ggplot2::scale_y_continuous(limits = c(min_range, max_range)) +
    ggplot2::labs(
      # x = x_lab,
      # y = y_lab,
      colour = colour_lab
    ) +
    # ggplot2::scale_size(trans = 'reverse') +
    ggplot2::scale_color_distiller(palette = "RdBu", type = "div", direction = -1, limits = c(min_range, max_range))
}

ui <- dashboardPage(
    skin = "black",
     # options = list(sidebarExpandOnHover = TRUE),
     header = dashboardHeader(title = "⚖️ TrapDiff", disable = FALSE),
     sidebar = dashboardSidebar(
        # Have an overflowing sidebar
        tags$style(
          "#sidebarItemExpanded {
                overflow: auto;
                max-height: 100vh;
            }"
        ),
       # minified = TRUE,
       # collapsed = FALSE,
       shinybusy::add_busy_spinner(spin = "dots", position = "bottom-right"),
       helpText("Path to Trapdiff data. Typically you do not need to change this."),
       textInput(
          inputId = "trappath",
          label = "Trapdiff output path",
          value = ""
        ),
        helpText("Table of the settings of the trapdiff run to be sure we have the right labels assigned."),
        tableOutput("settings"),
        helpText("Select a gene to highlight in the plots here. Alternatively you can select one in the table."),
        selectizeInput(
          "gene_id",
          "Highlight gene",
          choices = c()
        ),
        helpText("Filter genes that are significant in the following conditions. Typically you want to look at the genes significant with the 'interaction term'"),
        shiny::checkboxGroupInput(
          "select_significant",
          label = "Filter significant genes",
          choices = c()
        ),
        helpText("This affects the filter range of the scatter plot."),
        shiny::sliderInput(
          inputId = "scatterplot_range",
          label = "Scatterplot range:",
          min = 0,
          max = 50,
          value = 10
        ),
        shiny::textOutput("selection")
     ),
     body = dashboardBody(
       fluidRow(
        box(
          width = 4,
          collapsible = TRUE,
          plotOutput("plot_group_counts")
        ),
        box(
          width = 4,
          collapsible = TRUE,
          plotOutput("stats_plot")
        ),
        box(
          width = 4,
          collapsible = TRUE,
          plotOutput("ratio_plot")
        )
        # box(
        #   width = 3,
        #   collapsible = TRUE,
        #   plotly::plotlyOutput("main_scatterplot")
        # )
      ),
      fluidRow(
        box(
          width = 12,
          collapsible = TRUE,
          DT::DTOutput("table")
        )
      )
     ),
     # controlbar = dashboardControlbar(),
     title = "TrapDiff"
   )

server <- function(input, output, session) {
  # Look for URL search changes
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['trappath']])) {
      updateTextInput(session, "trappath", value = query[['trappath']])
    }
  })
  apply_filter <- function(de_wide) {
    if (main_effect_a_comparison_name() %in% input$select_significant) {
      de_wide <- de_wide %>%
        dplyr::filter(get(glue::glue("padj_{main_effect_a_comparison_name()}")) <= 0.05)
    }
    if (main_effect_b_comparison_name() %in% input$select_significant) {
      de_wide <- de_wide %>%
        dplyr::filter(get(glue::glue("padj_{main_effect_b_comparison_name()}")) <= 0.05)
    }
    if (source_effect_a_comparison_name() %in% input$select_significant) {
      de_wide <- de_wide %>%
        dplyr::filter(get(glue::glue("padj_{source_effect_a_comparison_name()}")) <= 0.05)
    }
    if (source_effect_b_comparison_name() %in% input$select_significant) {
      de_wide <- de_wide %>%
        dplyr::filter(get(glue::glue("padj_{source_effect_b_comparison_name()}")) <= 0.05)
    }
    if ("interaction_effect" %in% input$select_significant) {
      de_wide <- de_wide %>%
        dplyr::filter(get(glue::glue("padj_interaction_effect")) <= 0.05)
    }
    return(de_wide)
  }

  # Read data
  de <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    # Attach ratio t_test to the data frame to be able to compare its results
    de_out <-
      glue::glue("{input$trappath}/de.rds") %>%
      readRDS()
    
    path_ratio_t_test <- paste0(input$trappath, "/ratio_t_test.csv.gz")
    # Attach t-test if it exists
    if (file.exists(path_ratio_t_test)) {
      ratio_t_test <- readr::read_csv(path_ratio_t_test)
      de_out <-
        de_out %>%
        dplyr::bind_rows(
          ratio_t_test %>%
          dplyr::select(ensembl_gene_id, external_gene_name, log2FoldChange, p_value) %>%
          dplyr::mutate(comparison = "ratio") %>%
          dplyr::rename(padj = p_value)
        )
    }
    # Return result
    de_out %>%
      dplyr::mutate(gene_id = glue::glue("{ensembl_gene_id}_{external_gene_name}"))
  })
  col_data <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    readRDS(glue::glue("{input$trappath}/col_data.rds"))
  })
  tpms <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    readRDS(glue::glue("{input$trappath}/tpms.rds"))
  })
  cpms <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    readRDS(glue::glue("{input$trappath}/cpms.rds"))
  })
  counts <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    readr::read_csv(glue::glue("{input$trappath}/deseq_normalized_counts.csv.gz"))
  })
  ratio_counts <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    glue::glue("{input$trappath}/ratio_counts.csv.gz") %>%
      readr::read_csv() %>%
      dplyr::mutate(gene_id = glue::glue("{ensembl_gene_id}_{external_gene_name}"))

  })

  de_wide <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    de() %>%
      tidyr::pivot_wider(
        names_from = comparison,
        values_from = c(padj, log2FoldChange),
        # For some reason we'll get problems with duplicate entries if we don't do this
        values_fn = mean
      ) %>%
      rmyknife::attach_biomart(attributes = "description")
  })
    treatment_a <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    dat <- tryCatch(
      col_data()$treatment %>%
        levels() %>%
        .[1],
      error = function(cond){ NULL }
    )
    return(dat)
  })
  treatment_b <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    dat <- tryCatch(
      col_data()$treatment %>%
        levels() %>%
        .[2],
      error = function(cond){ NULL }
      )
    return(dat)
  })
  source_a <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    dat <- tryCatch(
      col_data()$source %>%
        levels() %>%
        .[1],
      error = function(cond){ NULL }
      )
    return(dat)
  })
  source_b <- shiny::reactive({
    shiny::validate(
      need(input$trappath != "", "Please provide valid data set")
    )
    dat <- tryCatch(
      col_data()$source %>%
        levels() %>%
        .[2],
      error = function(cond){ NULL }
      )
    return(dat)
  })
# main_effect_a_comparison_name <- glue::glue("msa_{treatment_a()}_{treatment_b()}_{source_a()}")
#   main_effect_b_comparison_name <- glue::glue("msb_{treatment_a()}_{treatment_b()}_{source_b()}")
#   source_effect_a_comparison_name <- glue::glue("sma_{source_a()}_{source_b()}_{treatment_a()}")
#   source_effect_b_comparison_name <- glue::glue("smb_{source_a()}_{source_b()}_{treatment_b()}")
  main_effect_a_comparison_name <- shiny::reactive({ glue::glue("{treatment_a()}_{treatment_b()}_{source_a()}") })
  main_effect_b_comparison_name <- shiny::reactive({ glue::glue("{treatment_a()}_{treatment_b()}_{source_b()}") })
  source_effect_a_comparison_name <- shiny::reactive({ glue::glue("{source_a()}_{source_b()}_{treatment_a()}") })
  source_effect_b_comparison_name <- shiny::reactive({ glue::glue("{source_a()}_{source_b()}_{treatment_b()}") })

  gene_select <- shiny::reactive({
    selected_id <- get_first_selected_in_dt(input)
    if (!is.null(selected_id)) {
      print(selected_id)
      selected_gene <- de_wide() %>%
        apply_filter() %>%
        dplyr::slice(selected_id) %>%
        .$gene_id
      print(selected_gene)
      return(selected_gene)
    }
  })

  # Update filter box
  observe({
    # We need to to the length-thing to catch a really weird behavior
    # in R which outputs NULL == "bla" to logical(0). For details check
    # https://stackoverflow.com/questions/27350636/r-argument-is-of-length-zero-in-if-statement
    if (length(treatment_a()) > 0) {
      if (!is.na(treatment_a())) {
        updateCheckboxGroupInput(
          session,
          "select_significant",
          choices = 
            tibble::tibble(
              !!(main_effect_a_comparison_name()) := main_effect_a_comparison_name(),
              !!(main_effect_b_comparison_name()) := main_effect_b_comparison_name(),
              !!(source_effect_a_comparison_name()) := source_effect_a_comparison_name(),
              !!(source_effect_b_comparison_name()) := source_effect_b_comparison_name(),
              interaction_effect = "interaction_effect"
            ) %>% as.list()
        )
      }
    }
  })

  # Update Gene select input
  observe({
    selected_gene <- tryCatch(
      de_wide() %>%
        apply_filter() %>%
        .$gene_id %>%
        .[1],
      error = function(cond){ NULL }
    )
    gene_choices <- tryCatch(
      de_wide() %>%
        dplyr::select(gene_id) %>%
        dplyr::distinct() %>%
        dplyr::pull(),
      error = function(cond){ NULL }
    )
    if (!is.null(input$table_rows_selected)) {
      selected_gene <- de_wide() %>%
        apply_filter() %>%
        dplyr::slice(input$table_rows_selected) %>%
        .$gene_id
    }
    updateSelectizeInput(
      session,
      "gene_id",
      choices = gene_choices,
      selected = selected_gene,
      server = TRUE
    )
  })

  output$settings <- renderTable({
    tibble::tibble(
      condition = c("a", "b"),
      treatment = c(treatment_a(), treatment_b()),
      source = c(source_a(), source_b())
    )
  })

  output$table <- DT::renderDataTable({
    de_wide() %>%
      dplyr::select(ensembl_gene_id, external_gene_name, description, dplyr::everything()) %>%
      # dplyr::select(ensembl_gene_id, external_gene_name, dplyr::everything()) %>%
      apply_filter() %>%
      # Replace long names to make list smaller
      (function(dat) {
        colnames(dat) <- stringr::str_replace_all(colnames(dat), "log2FoldChange", "fc")
        colnames(dat) <- stringr::str_replace_all(colnames(dat), "padj", "p")
        colnames(dat) <- stringr::str_replace_all(colnames(dat), "interaction_effect", "interaction")
        return(dat)
      }) %>%
      # Round all values to 4 decimals
      dplyr::mutate_if(is.numeric, round, 4) %>%
      DT::datatable(
        extensions = c("Scroller", "Buttons"),
        selection = "single",
        filter = list(position = "top"),
        options = list(
          # Remove search bar but leave filter
          # https://stackoverflow.com/a/35627085
          # sDom  = '<"top">lrt<"bottom">ip',
          # Enable colvis button
          dom = "Bfrtip",
          # Define hidden columns
          columnDefs = list(list(
            targets = 0:3, visible = FALSE
          )),
          buttons = c("colvis", "copy", "csv", "excel", "pdf", "print"),
          scrollX = TRUE
          # Makes the table more responsive when it's really big
          # scrollY = 250,
          # scroller = TRUE
        )
      )
  })
  # Quick and dirty scatterplot selection.
  # See https://stackoverflow.com/questions/48939382/how-can-i-grab-the-row-of-data-from-a-ggplotly-in-shiny
  output$selection <- renderPrint({
    s <- event_data("plotly_click")
    cat("You selected: \n\n")
    data.frame(s)
    updateSelectInput(
      session,
      "gene_id",
      selected = s$key
    )
    return(s)
  })
  # scatterplot_selection <- reactive({
  #   s <- event_data("plotly_click")
  #   cat("You selected: \n\n")
  #   df <- data.frame(s)
  # })
  output$main_scatterplot <- plotly::renderPlotly({
    plot <- de_wide() %>%
      apply_filter() %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes_string(
          x = glue::glue("log2FoldChange_{main_effect_a_comparison_name()}"),
          y = glue::glue("log2FoldChange_{main_effect_b_comparison_name()}"),
          color = "log2FoldChange_interaction_effect",
          label = "external_gene_name",
          key = "gene_id"
        )
      ) %>%
      fc_scatterplot(
        colour_lab = "FC inter.",
        point_alpha = 1,
        min_range = -input$scatterplot_range,
        max_range = input$scatterplot_range
      ) +
      ggplot2::geom_hline(yintercept = 0, size = 0.05) +
      ggplot2::geom_vline(xintercept = 0, size = 0.05) +
      ggplot2::geom_point(
        data = de_wide() %>% dplyr::filter(gene_id == input$gene_id),
        color = "red",
        size = 2
      )
    plotly::ggplotly(plot)
  })

  output$plot_group_counts <- shiny::renderPlot({
    counts() %>%
      dplyr::mutate(
        gene_id = paste0(ensembl_gene_id, "_", external_gene_name),
        sample_id = glue::glue("{source}_{sample_id}"),
        group = glue::glue("{source}_{treatment}")
      ) %>%
      dplyr::filter(gene_id %in% input$gene_id) %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = group,
          y = count
        )
      ) +
      # Remove outliers from boxplot to not confuse them with the jitter data points
      ggplot2::geom_boxplot(outlier.shape = NA, color = "grey") +
      ggplot2::geom_jitter(
        mapping = ggplot2::aes(color = sample_id),
        width = 0.1,
        height = 0
      ) +
      ggplot2::ggtitle(
        "DESeq2 Counts per condition",
        input$gene_id
      )
  })

  output$plot_each_group <- shiny::renderPlot({
    # Create a tidy long format of tpms
    plot_dat <- dplyr::left_join(
      x = col_data() %>%
        dplyr::mutate(group = glue::glue("{source}_{treatment}")) %>%
        dplyr::select(id, group),
      y = cpms(),
      by = c("id" = "sample")
    ) %>%
      # Filter all genes where we do not get a gene name
      dplyr::filter(!is.na(external_gene_name)) %>%
      # Create a unique gene_id variable for filtering
      dplyr::mutate(
        gene_id = paste0(ensembl_gene_id, "_", external_gene_name)
      ) %>%
      # Attach p-value from interaction effect
      dplyr::left_join(
        de() %>%
          tidyr::pivot_wider(
            names_from = comparison,
            values_from = c(padj, log2FoldChange)
          ) %>%
          dplyr::select(ensembl_gene_id, padj_interaction_effect),
        by = "ensembl_gene_id"
      ) %>%
      # Attach sample ID
      dplyr::mutate(sample_id = stringr::str_extract(id, "\\d*$")) %>%
      dplyr::mutate(sample_id = glue::glue("{source}_{sample_id}"))

    plot_dat %>%
      dplyr::filter(gene_id %in% input$gene_id) %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = group,
          y = cpm
        )
      ) +
      # Remove outliers from boxplot to not confuse them with the jitter data points
      ggplot2::geom_boxplot(outlier.shape = NA, color = "grey") +
      ggplot2::geom_jitter(
        mapping = ggplot2::aes(color = sample_id),
        width = 0.1,
        height = 0
      ) +
      ggplot2::ggtitle(
        "CPMs per condition",
        input$gene_id
      )
  })

  output$ratio_plot <- shiny::renderPlot({
    ratio_counts() %>%
      dplyr::mutate(
        gene_id = paste0(ensembl_gene_id, "_", external_gene_name),
        sample = glue::glue("{source}_{sample_id}")
      ) %>%
      dplyr::filter(gene_id %in% input$gene_id) %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = source,
          y = log2(count_norm)
        )
      ) +
      # Remove outliers from boxplot to not confuse them with the jitter data points
      ggplot2::geom_boxplot(outlier.shape = NA, color = "grey") +
      ggplot2::geom_jitter(
        width = 0.1,
        height = 0
      ) +
      ggplot2::ggtitle(
        "IP/Input ratios per condition",
        input$gene_id
      )
  })

  output$stats_plot <- shiny::renderPlot({
    p1 <- de() %>%
      dplyr::filter(gene_id %in% input$gene_id) %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = padj,
          y = comparison
        )
      ) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed")
    
    p2 <- de() %>%
      dplyr::filter(gene_id %in% input$gene_id) %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = log2FoldChange,
          y = comparison
        )
      ) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = 0)
      # Create title
      title <- cowplot::ggdraw() +
        cowplot::draw_label(
          input$gene_id,
          x = 0,
          hjust = 0
        )
      cowplot::plot_grid(
        title, p1, p2,
        ncol = 1,
        # rel_heights values control vertical title margins
        rel_heights = c(0.1, 1, 1)
      )
  })
}

shinyApp(ui, server)
