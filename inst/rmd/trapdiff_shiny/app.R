library(shiny)

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

ui <- fluidPage(
  titlePanel("censusVis"),
  sidebarLayout(
    sidebarPanel(
      helpText("Trapdiff Help Text"),
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
      shiny::sliderInput(
        inputId = "scatterplot_range",
        label = "Scatterplot range:",
        min = 0,
        max = 50,
        value = 5
      ),
      shiny::checkboxGroupInput(
        "select_significant",
        label = "Filter significant genes",
        choices = c()
      )
    ),

    mainPanel(
      tableOutput("show_path"),
      tableOutput("show_de_wide"),
      textOutput("debug_text"),
      DT::DTOutput("table"),
      plotOutput("main_scatterplot")
    )
  )
)

server <- function(input, output, session) {
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

  # Debug output
  output$show_path <- renderTable({
    head(de())
  })
  output$show_de_wide <- renderTable({
    head(col_data())
  })
  output$debug_text <- renderText({
    print(glue::glue("main_effect_a_comparison_name: {main_effect_a_comparison_name()}"))
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
        extensions = c("Scroller", "Buttons"),
        selection = "single",
        filter = list(position = "top"),
        options = list(
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel", "pdf", "print"),
          scrollX = TRUE,
          # Makes the table more responsive when it's really big
          scrollY = 300,
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
}

shinyApp(ui, server)
