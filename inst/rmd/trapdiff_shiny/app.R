library(shiny)
library(shinydashboard)

fc_scatterplot <- function(
  ggplot_dat,
  colour_lab,
  title,
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
    ggplot2::scale_color_distiller(palette = "RdBu", type = "div", direction = -1, limits = c(min_range, max_range)) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(title)
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
       textInput(
          inputId = "trappath",
          label = "Trapdiff output path",
          value = "/beegfs/scratch/bruening_scratch/pklemm/2020-05-yiyi-bactrap/release/trapdiff"
        ),
        tableOutput("settings"),
        selectInput(
          "gene_id",
          "Highlight gene",
          choices = c()
        ),
        shiny::checkboxGroupInput(
          "select_significant",
          label = "Filter significant genes",
          choices = c()
        ),
        shiny::sliderInput(
          inputId = "scatterplot_range",
          label = "Scatterplot range:",
          min = 0,
          max = 50,
          value = 5
        )
     ),
     body = dashboardBody(
       fluidRow(
        box(
          width = 4,
          collapsible = TRUE,
          plotOutput("plot_each_group")
        ),
        box(
          width = 4,
          collapsible = TRUE,
          plotOutput("stats_plot")
        ),
        box(
          width = 4,
          collapsible = TRUE,
          plotOutput("main_scatterplot")
        ),
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
    glue::glue("{input$trappath}/de.rds") %>%
      readRDS() %>%
      dplyr::mutate(gene_id = glue::glue("{ensembl_gene_id}_{external_gene_name}"))
  })
  col_data <- shiny::reactive({
    readRDS(glue::glue("{input$trappath}/col_data.rds"))
  })
  counts <- shiny::reactive({
    readRDS(glue::glue("{input$trappath}/counts.rds"))
  })
  de_wide <- shiny::reactive({
    de() %>%
      tidyr::pivot_wider(
        names_from = comparison,
        values_from = c(padj, log2FoldChange)
      )
  })
  treatment_a <- shiny::reactive({
    col_data()$treatment %>%
      levels() %>%
      .[1]
  })
  treatment_b <- shiny::reactive({
    col_data()$treatment %>%
      levels() %>%
      .[2]
  })
  source_a <- shiny::reactive({
    col_data()$source %>%
      levels() %>%
      .[1]
  })
  source_b <- shiny::reactive({
    col_data()$source %>%
      levels() %>%
      .[2]
  })
  main_effect_a_comparison_name <- shiny::reactive({ glue::glue("main_in_source_a___{treatment_a()}_vs_{treatment_b()}_in_{source_a()}") })
  main_effect_b_comparison_name <- shiny::reactive({ glue::glue("main_in_source_b___{treatment_a()}_vs_{treatment_b()}_in_{source_b()}") })
  source_effect_a_comparison_name <- shiny::reactive({ glue::glue("source_in_main_a___{source_a()}_vs_{source_b()}_in_{treatment_a()}") })
  source_effect_b_comparison_name <- shiny::reactive({ glue::glue("source_in_main_b___{source_a()}_vs_{source_b()}_in_{treatment_b()}") })

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
  })

  # Update Gene select input
  observe({
    selected_gene <-
      de_wide() %>%
      apply_filter() %>%
      .$gene_id %>%
      .[1]
    if (!is.null(input$table_rows_selected)) {
      # browser()
      selected_gene <- de_wide() %>%
        apply_filter() %>%
        dplyr::slice(input$table_rows_selected) %>%
        .$gene_id
    }
    updateSelectInput(
      session,
      "gene_id",
      choices = de_wide() %>%
        dplyr::select(gene_id) %>%
        dplyr::distinct() %>%
        dplyr::pull(),
      # selected = gene_select()
      selected = selected_gene
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
      # rmyknife::attach_biomart(attributes = "description", ensembl_version = 100) %>%
      # dplyr::select(ensembl_gene_id, external_gene_name, description, dplyr::everything()) %>%
      dplyr::select(ensembl_gene_id, external_gene_name, dplyr::everything()) %>%
      apply_filter() %>%
      DT::datatable(
        # extensions = c("Scroller", "Buttons"),
        extensions = c("Scroller"),
        selection = "single",
        filter = list(position = "top"),
        options = list(
          # Remove search bar but leave filter
          # https://stackoverflow.com/a/35627085
          sDom  = '<"top">lrt<"bottom">ip',
          # dom = "Bfrtip",
          # buttons = c("copy", "csv", "excel", "pdf", "print"),
          scrollX = TRUE,
          # Makes the table more responsive when it's really big
          scrollY = 250,
          scroller = TRUE
        )
      )
  })

  output$main_scatterplot <- renderPlot({
    de_wide() %>%
      apply_filter() %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes_string(
          x = glue::glue("log2FoldChange_{main_effect_a_comparison_name()}"),
          y = glue::glue("log2FoldChange_{main_effect_b_comparison_name()}"),
          color = "log2FoldChange_interaction_effect",
          label = "external_gene_name"
        )
      ) %>%
      fc_scatterplot(
        colour_lab = "FC interaction effect",
        title = "FC main effect coloured by FC interaction effect",
        point_alpha = 1,
        min_range = -input$scatterplot_range,
        max_range = input$scatterplot_range
      ) +
      ggplot2::geom_point(
        data = de_wide() %>% dplyr::filter(gene_id == input$gene_id),
        color = "red",
        size = 10
      )
  })

  output$plot_each_group <- shiny::renderPlot({
    # Create a tidy long format of counts
    plot_dat <- dplyr::left_join(
      x = col_data() %>%
        dplyr::mutate(group = glue::glue("{source}_{treatment}")) %>%
        dplyr::select(id, group),
      y = counts(),
      by = c("id" = "sample")
    ) %>%
      # Filtr all genes where we do not get a gene name
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
      )

    plot_dat %>%
      dplyr::filter(gene_id %in% input$gene_id) %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = group,
          y = value
        )
      ) +
      ggplot2::geom_violin() +
      ggplot2::geom_jitter(width = 0.1) +
      ggplot2::ggtitle(input$gene_id)
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
