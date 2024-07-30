library(shiny)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(DT)
library(bslib)
library(gridExtra)


options(shiny.maxRequestSize = 10 * 1024^3)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("scRNA-seq Quality Control"),
  htmlOutput("main_header"),
  sidebarLayout(
    sidebarPanel(
      fileInput("RDS_object", "Load your pre-filtered scRNA-seq RDS", placeholder = "Load RDS object", accept = ".rds", buttonLabel = "Select RDS File..."),
      sliderInput("UMIs_per_cell_threshold", "UMIs per cell threshold:",
                  min = 0, max = 5000,
                  value = 500, step = 50),
      sliderInput("nGenes_per_cell_threshold", "Genes per cell threshold:",
                  min = 0, max = 2000,
                  value = 250, step = 25),
      sliderInput("mito_threshold", "Percent Mitochondria threshold:",
                  min = 0, max = 1,
                  value = 0.2),
      sliderInput("ribo_threshold", "Percent Ribosomal threshold:",
                  min = 0, max = 1,
                  value = 0),
      sliderInput("Log10GenesPerUMI_threshold", "Log10 Genes Per UMI threshold:",
                  min = 0, max = 1,
                  value = 0.8),
      actionButton("evaluate", "Evaluate!", class = "btn-primary"),
      conditionalPanel(
        condition = "input.evaluate > 0",
        downloadButton("download_RDS", "Download Filtered RDS", class = "btn-success")
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("UMIs per cell", br(),
                 htmlOutput("UMIs_per_cell_header"),
                 navset_card_tab(
                   nav_panel("Density Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("density_UMI_per_cell_raw"), plotOutput("density_UMI_per_cell_post_filter"),
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                conditionalPanel(
                                  condition = "input.evaluate > 0",
                                  downloadButton("download_density_UMIs_raw", "Download raw PNG", class = "btn-success", align = "center")
                               )),
                               column(1),
                               column(5,
                               conditionalPanel(
                                 condition = "input.evaluate > 0",
                                 downloadButton("download_density_UMIs_filtered", "Download filtered PNG", class = "btn-success")
                               ))
                             )
                   ),
                   nav_panel("Violin Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("violin_UMI_per_cell_raw"), plotOutput("violin_UMI_per_cell_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_UMIs_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_UMIs_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   )
                 )
        ),
        tabPanel("Genes per cell", br(),
                 htmlOutput("genes_per_cell_header"),
                 navset_card_tab(
                   nav_panel("Density Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("density_nGene_raw"), plotOutput("density_nGene_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_density_nGene_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_density_nGene_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   ),
                   nav_panel("Violin Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("violin_nGene_raw"), plotOutput("violin_nGene_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_nGene_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_nGene_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   )
                 )
        ),
        tabPanel("Percent Mitochondria", br(),
                 htmlOutput("mitoRatio_header"),
                 navset_card_tab(
                   nav_panel("Density Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("density_mito_raw"), plotOutput("density_mito_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_density_mito_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_density_mito_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   ),
                   nav_panel("Violin Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("violin_mito_raw"), plotOutput("violin_mito_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_mito_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_mito_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   )
                 )
        ),
        tabPanel("Percent Ribosomal", br(),
                 htmlOutput("riboRatio_header"),
                 navset_card_tab(
                   nav_panel("Density Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("density_ribo_raw"), plotOutput("density_ribo_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_density_ribo_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_density_ribo_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   ),
                   nav_panel("Violin Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("violin_ribo_raw"), plotOutput("violin_ribo_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_ribo_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_ribo_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   )
                 )
        ),
        tabPanel("Log10 Genes per UMI", br(),
                 htmlOutput("Log10GenesPerUMI_header"),
                 navset_card_tab(
                   nav_panel("Density Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("density_Log10GenesPerUMI_raw"), plotOutput("density_Log10GenesPerUMI_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_density_Log10GenesPerUMI_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_density_Log10GenesPerUMI_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   ),
                   nav_panel("Violin Plots",
                             fluidRow(
                               splitLayout(cellWidths = c("50%", "50%"), plotOutput("violin_Log10GenesPerUMI_raw"), plotOutput("violin_Log10GenesPerUMI_post_filter")
                               )
                             ),
                             fluidRow(
                               column(1),
                               column(5, 
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_Log10GenesPerUMI_raw", "Download raw PNG", class = "btn-success", align = "center")
                                      )),
                               column(1),
                               column(5,
                                      conditionalPanel(
                                        condition = "input.evaluate > 0",
                                        downloadButton("download_violin_Log10GenesPerUMI_filtered", "Download filtered PNG", class = "btn-success")
                                      ))
                             )
                   )
                 )
        ),
        tabPanel("UMIs vs. Genes", br(),
                 htmlOutput("UMIs_vs_genes_header"),
                 splitLayout(cellWidths = c("50%", "50%"), plotOutput("UMI_vs_Genes_raw"), plotOutput("UMI_vs_Genes_post_filter")
                 ),
                 fluidRow(
                   column(1),
                   column(5, 
                          conditionalPanel(
                            condition = "input.evaluate > 0",
                            downloadButton("download_UMI_vs_Genes_raw", "Download raw PNG", class = "btn-success", align = "center")
                          )),
                   column(1),
                   column(5,
                          conditionalPanel(
                            condition = "input.evaluate > 0",
                            downloadButton("download_UMI_vs_Genes_filtered", "Download filtered PNG", class = "btn-success")
                          ))
                 )
        )
      ),
      DTOutput("filtered_table"),
      conditionalPanel(
        condition = "input.evaluate > 0",
        downloadButton("download_raw_filtered_table", "Download Table as CSV", class = "btn-success")
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  seurat_merge <- reactive({
    req(input$RDS_object)
    readRDS(input$RDS_object$datapath)
  })
  output$main_header <- renderText ({
    paste0("In this RShiny app you will be able to examine 6 different quality metrics in your dataset. The goals of scRNA QC are to <b>(1) filter the data to only include true cells that are of high quality</b> and remove doublets, empty droplets and stressed/dying cells and <b>(2) to identify any failed samples</b> and either try to salvage the data or remove from analysis. We also to try to understand why the sample failed. High data quality will make downstream steps such as identifying distinct cell type populations much easier.")
  })
  output$UMIs_per_cell_header <- renderText ({
    paste0("Here, we look at the distribution of UMIs (unique molecular identifiers, or sequenced reads) per cell (droplet) in the dataset. Before QC, we expect a biomodal distribution with a first <i>small</i> peak at low numbers of UMIs (<250) corresponding to droplets that encapsulated background/dying cells, and a second higher peak centered at >1000. After QC we expect a unimodal  distribution. If UMI counts are between 500-1000 counts, the data is usable  but the cells probably should have been sequenced more deeply.")
  })
  output$genes_per_cell_header <- renderText({
    paste0("Here, we look at the number of different genes that were detected in each cell. By \"detected\", we mean genes with a non-zero read count  measurement. Gene detection in the range of 500 to 5000 is normal for most single-cell experiments. For high quality data, the proportional histogram  should contain <b>a single large peak that represents cells that were  encapsulated</b>. If we see a <b>small shoulder</b> to the left of the major peak (not present in our data), or a bimodal distribution of the cells, that can indicate a couple of things. It might be that there are a set of <b>cells that failed</b> for  some reason. It could also be that there are <b>biologically different types of cells</b> (i.e. quiescent cell populations, less complex cells of interest), and/or one type is much smaller than the other (i.e. cells with high counts may be cells that are larger in size). Therefore, this threshold should be assessed with other metrics.")
  })
  output$mitoRatio_header <- renderText({
    paste0("We evaluate overall mitochondrial gene expression as a biomarker of cellular stress during sample preparation. Typically, we expect mitochondrial genes  to account for <20% of overall transcripts in each cell. Note that if you are  working with single-nuclei dataset you should expect values almost at zero  for this mitochondrial ratio as the cytoplasm (where mitochondria DNA lives) is removed in the sample preparation stage.")
  })
  output$riboRatio_header <- renderText({
    paste0("The expression of ribosomal genes give us insight into the cellular activity of a cell. Different cells types are expected to have different levels of ribosomal expression. Due to this there is no suggested cutoff for Ribosomal ratio. We merely expect it to be similar among samples with similar cellular composition and note that extremely high levels can indicate low quality reads.")
  })
  output$Log10GenesPerUMI_header <- renderText({
    paste0("We can evaluate each cell in terms of how complex the RNA species are byusing a measure called the novelty score. The novelty score is computed by taking the ratio of nGenes over nUMI. If there are many captured transcripts (high nUMI) and a low number of genes detected in a cell, this likely means that you only captured a low number of genes and simply sequenced transcripts from those lower number of genes repeatedly. Typical values for  this metric are >0.8 for most cells. Cells with lower diversity in the genes  they express may be low-complexity cell types such as red blood cells.")
  })
  output$UMIs_vs_genes_header <- renderText({
    paste0("Considering any of these QC metrics in isolation can lead to misinterpretationof cellular signals. For example, cells with a comparatively high fraction of mitochondrial counts may be involved in respiratory processes and may be cells that you would like to keep. Likewise, other metrics can have other biological interpretations. A general rule of thumb when performing QC is to <b>set thresholds for individual metrics to be as permissive as possible, and always consider the joint effects</b> of these metrics. In this way, you reduce the risk of filtering out any viable cell populations. Two metrics that are often evaluated together are the number of UMIs and the number of genes detected per cell.<br>By plotting the number of UMIs per cell (x-axis) vs. the number of genes per cell (y-axis), we can visually assess whether there is a large proportion of lowquality cells with low read counts and/or gene detection (bottom left quadrant of the plot). In the following representation, cells are further color-coded based on the percentage of mitochondrial genes found among total detected genes.")
  })
  bindEvent(observe({
    filter_seurat <- reactive({
      req(seurat_merge())
      subset(x = seurat_merge(), 
             subset = (nCount_RNA >= isolate(input$UMIs_per_cell_threshold))
             & (nFeature_RNA >= isolate(input$nGenes_per_cell_threshold))
             & (mitoRatio < isolate(input$mito_threshold))
             & (riboRatio > isolate(input$ribo_threshold))
             & (Log10GenesPerUMI > isolate(input$Log10GenesPerUMI_threshold))
      )
    })
    raw_filtered_table <- reactive({
      cbind("Raw" = table(seurat_merge()$orig.ident), 
            "Filtered" = table(filter_seurat()$orig.ident), 
            "Percent retained" = paste0(round((table(filter_seurat()$orig.ident) / table(seurat_merge()$orig.ident) * 100), digits = 2), "%"))  %>% 
        datatable(options = list(dom = 't', pageLength = length(unique(seurat_merge()@meta.data$orig.ident))), caption = "Table with an overview of our cells passing various filtering standards")
    })
    raw_filtered_table_download <- reactive({
      cbind("Raw" = table(seurat_merge()$orig.ident), 
            "Filtered" = table(filter_seurat()$orig.ident), 
            "Percent retained" = paste0(round((table(filter_seurat()$orig.ident) / table(seurat_merge()$orig.ident) * 100), digits = 2), "%"))
    })
    output$filtered_table <- renderDT({
      raw_filtered_table()
    })
    density_UMI_per_cell_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(x = nCount_RNA, color = orig.ident, fill = orig.ident)) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() + 
        geom_vline(xintercept = isolate(input$UMIs_per_cell_threshold)) +
        facet_wrap(. ~ orig.ident) +
        ggtitle("Detected genes per cell in raw dataset") +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_UMI_per_cell_raw <- renderPlot({
      density_UMI_per_cell_raw_plot()
    })
    density_UMI_per_cell_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x = nCount_RNA, color = orig.ident, fill = orig.ident)) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() + 
        facet_wrap(. ~ orig.ident) +
        ggtitle("Detected genes per cell in filtered dataset") +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_UMI_per_cell_post_filter <- renderPlot({
      density_UMI_per_cell_post_filter_plot()
    })
    violin_UMI_per_cell_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=log10(nCount_RNA), fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        geom_hline(yintercept = isolate(log10(input$UMIs_per_cell_threshold))) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("Detected genes per cell in raw dataset")
    })
    output$violin_UMI_per_cell_raw <- renderPlot({
      violin_UMI_per_cell_raw_plot()
    })
    violin_UMI_per_cell_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=log10(nCount_RNA), fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("Detected genes per cell in filtered dataset")
    })
    output$violin_UMI_per_cell_post_filter <- renderPlot({
      violin_UMI_per_cell_post_filter_plot()
    })
    density_nGene_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(x = nFeature_RNA, color = orig.ident, fill = orig.ident)) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() + 
        geom_vline(xintercept = isolate(input$nGenes_per_cell_threshold)) +
        facet_wrap(. ~ orig.ident) +
        ggtitle("Detected genes per cell in raw dataset") +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_nGene_raw <- renderPlot({
      density_nGene_raw_plot()
    })
    density_nGene_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x = nFeature_RNA, color = orig.ident, fill = orig.ident)) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() + 
        facet_wrap(. ~ orig.ident) +
        ggtitle("Detected genes per cell in filtered dataset") +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_nGene_post_filter <- renderPlot({
      density_nGene_post_filter_plot()
    })
    violin_nGene_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        geom_hline(yintercept = isolate(log10(input$nGenes_per_cell_threshold))) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("Detected genes per cell in raw dataset")
    })
    output$violin_nGene_raw <- renderPlot({
      violin_nGene_raw_plot()
    })
    violin_nGene_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("Detected genes per cell in filtered dataset")
    })
    output$violin_nGene_post_filter <- renderPlot({
      violin_nGene_post_filter_plot()
    })
    density_mito_raw_plot <-reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(color = orig.ident, x = mitoRatio, fill = orig.ident)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        geom_vline(xintercept = isolate(input$mito_threshold)) + 
        facet_wrap(. ~ orig.ident) +
        ggtitle("Percentage of mitochondrial gene expression per cell in raw dataset") +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_mito_raw <- renderPlot({
      density_mito_raw_plot()
    })
    density_mito_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(color = orig.ident, x = mitoRatio, fill = orig.ident)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        facet_wrap(. ~ orig.ident) +
        ggtitle("Percentage of mitochondrial gene expression per cell in filtered dataset") +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_mito_post_filter <- renderPlot({
      density_mito_post_filter_plot()
    })
    violin_mito_raw_plot <- reactive({
      seurat_merge()@meta.data %>%
        ggplot(aes(x=orig.ident, y=mitoRatio, fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        geom_hline(yintercept = isolate(input$mito_threshold)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("Percentage of mitochondrial gene expression per cell in raw dataset")
    })
    output$violin_mito_raw <- renderPlot({
      violin_mito_raw_plot()
    })
    violin_mito_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=mitoRatio, fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("Percentage of mitochondrial gene expression per cell in filtered dataset")
    })
    output$violin_mito_post_filter <- renderPlot({
      violin_mito_post_filter_plot()
    })
    density_ribo_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(color = orig.ident, x = riboRatio, fill = orig.ident)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        geom_vline(xintercept = isolate(input$ribo_threshold)) + 
        facet_wrap(. ~ orig.ident) +
        ggtitle("Percentage of ribosomal gene expression per cell in raw dataset") +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_ribo_raw <- renderPlot({
      density_ribo_raw_plot()
    })
    density_ribo_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(color = orig.ident, x = riboRatio, fill = orig.ident)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        facet_wrap(. ~ orig.ident) +
        ggtitle("Percentage of ribosomal gene expression per cell in filtered dataset") +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_ribo_post_filter <- renderPlot({
      density_ribo_post_filter_plot()
    })
    violin_ribo_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=riboRatio, fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        geom_hline(yintercept = isolate(input$ribo_threshold)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("Percentage of ribosomal gene expression per cell in raw dataset")
    })
    output$violin_ribo_raw <- renderPlot({
      violin_ribo_raw_plot()
    })
    violin_ribo_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=riboRatio, fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ggtitle("Percentage of ribosomal gene expression per cell in filtered dataset")
    })
    output$violin_ribo_post_filter <- renderPlot({
      violin_ribo_post_filter_plot()
    })
    density_Log10GenesPerUMI_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(x = Log10GenesPerUMI, color = orig.ident, fill = orig.ident)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        geom_vline(xintercept = isolate(input$Log10GenesPerUMI_threshold)) +
        facet_wrap(. ~ orig.ident) +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_Log10GenesPerUMI_raw <- renderPlot({
      density_Log10GenesPerUMI_raw_plot()
    })
    density_Log10GenesPerUMI_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x = Log10GenesPerUMI, color = orig.ident, fill = orig.ident)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        facet_wrap(. ~ orig.ident) +
        theme(axis.text.x = element_text(angle = 90))
    })
    output$density_Log10GenesPerUMI_post_filter <- renderPlot({
      density_Log10GenesPerUMI_post_filter_plot()
    })
    violin_Log10GenesPerUMI_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=Log10GenesPerUMI, fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        geom_hline(yintercept = isolate(input$Log10GenesPerUMI_threshold)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    })
    output$violin_Log10GenesPerUMI_raw <- renderPlot({
      violin_Log10GenesPerUMI_raw_plot()
    })
    violin_Log10GenesPerUMI_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x=orig.ident, y=Log10GenesPerUMI, fill=orig.ident)) + 
        geom_violin() + geom_boxplot(width = 0.1, fill = alpha("white", 0.7), outliers = FALSE) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    })
    output$violin_Log10GenesPerUMI_post_filter <- renderPlot({
      violin_Log10GenesPerUMI_post_filter_plot()
    })
    UMI_vs_Genes_raw_plot <- reactive({
      seurat_merge()@meta.data %>% 
        ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) + 
        geom_point() + 
        stat_smooth(method=lm) +
        scale_x_log10() + 
        scale_y_log10() + 
        theme_classic() +
        geom_vline(xintercept = isolate(input$UMIs_per_cell_threshold)) +
        geom_hline(yintercept = isolate(input$nGenes_per_cell_threshold)) + 
        ggtitle("Genes vs. nUMIs in raw dataset") + 
        xlab("nUMI") + ylab("nGene") + 
        facet_wrap(. ~ orig.ident)
    })
    output$UMI_vs_Genes_raw <- renderPlot({
      UMI_vs_Genes_raw_plot()
    })
    UMI_vs_Genes_post_filter_plot <- reactive({
      filter_seurat()@meta.data %>% 
        ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) + 
        geom_point() + 
        stat_smooth(method=lm) +
        scale_x_log10() + 
        scale_y_log10() + 
        theme_classic() +
        ggtitle("Genes vs. nUMIs in filtered dataset") + 
        xlab("nUMI") + ylab("nGene") + 
        facet_wrap(. ~ orig.ident)
    })
    output$UMI_vs_Genes_post_filter <- renderPlot({
      UMI_vs_Genes_post_filter_plot()
    })
    output$download_RDS <- downloadHandler(
      filename = function(){
        "filtered_seurat.RDS"
      },
      content = function(file){
            saveRDS(filter_seurat, file)
      }
    )
    output$download_density_UMIs_raw <- downloadHandler(
      filename = function(){
        "density_UMIs_per_cell_raw.png"
      },
      content = function(file){
        png(file)
        print(density_UMI_per_cell_raw_plot())
        dev.off()
      }
    )
    output$download_density_UMIs_filtered <- downloadHandler(
      filename = function(){
        "density_UMIs_per_cell_filtered.png"
      },
      content = function(file){
        png(file)
        print(density_UMI_per_cell_post_filter_plot())
        dev.off()
      }
    )
    output$download_violin_UMIs_raw <- downloadHandler(
      filename = function(){
        "violin_UMIs_per_cell_raw.png"
      },
      content = function(file){
        png(file)
        print(violin_UMI_per_cell_raw_plot())
        dev.off()
      }
    )
    output$download_violin_UMIs_filtered <- downloadHandler(
      filename = function(){
        "violin_UMIs_per_cell_filtered.png"
      },
      content = function(file){
        png(file)
        print(violin_UMI_per_cell_post_filter_plot())
        dev.off()
      }
    )
    output$download_density_nGene_raw <- downloadHandler(
      filename = function(){
        "density_nGene_raw.png"
      },
      content = function(file){
        png(file)
        print(density_nGene_raw_plot())
        dev.off()
      }
    )
    output$download_density_nGene_filtered <- downloadHandler(
      filename = function(){
        "density_nGene_filtered.png"
      },
      content = function(file){
        png(file)
        print(density_nGene_post_filter_plot())
        dev.off()
      }
    )
    output$download_violin_nGene_raw <- downloadHandler(
      filename = function(){
        "violin_nGene_raw.png"
      },
      content = function(file){
        png(file)
        print(violin_nGene_raw_plot())
        dev.off()
      }
    )
    output$download_violin_nGene_filtered <- downloadHandler(
      filename = function(){
        "violin_nGene_filtered.png"
      },
      content = function(file){
        png(file)
        print(violin_nGene_post_filter_plot())
        dev.off()
      }
    )
    output$download_density_mito_raw <- downloadHandler(
      filename = function(){
        "density_mito_raw.png"
      },
      content = function(file){
        png(file)
        print(density_mito_raw_plot())
        dev.off()
      }
    )
    output$download_density_mito_filtered <- downloadHandler(
      filename = function(){
        "density_mito_filtered.png"
      },
      content = function(file){
        png(file)
        print(density_mito_post_filter_plot())
        dev.off()
      }
    )
    output$download_violin_mito_raw <- downloadHandler(
      filename = function(){
        "violin_mito_raw.png"
      },
      content = function(file){
        png(file)
        print(violin_mito_raw_plot())
        dev.off()
      }
    )
    output$download_violin_mito_filtered <- downloadHandler(
      filename = function(){
        "violin_mito_filtered.png"
      },
      content = function(file){
        png(file)
        print(violin_mito_post_filter_plot())
        dev.off()
      }
    )
    output$download_density_ribo_raw <- downloadHandler(
      filename = function(){
        "density_ribo_raw.png"
      },
      content = function(file){
        png(file)
        print(density_ribo_raw_plot())
        dev.off()
      }
    )
    output$download_density_ribo_filtered <- downloadHandler(
      filename = function(){
        "density_ribo_filtered.png"
      },
      content = function(file){
        png(file)
        print(density_ribo_post_filter_plot())
        dev.off()
      }
    )
    output$download_violin_ribo_raw <- downloadHandler(
      filename = function(){
        "violin_ribo_raw.png"
      },
      content = function(file){
        png(file)
        print(violin_ribo_raw_plot())
        dev.off()
      }
    )
    output$download_violin_ribo_filtered <- downloadHandler(
      filename = function(){
        "violin_ribo_filtered.png"
      },
      content = function(file){
        png(file)
        print(violin_ribo_post_filter_plot())
        dev.off()
      }
    )
    output$download_density_Log10GenesPerUMI_raw <- downloadHandler(
      filename = function(){
        "density_Log10GenesPerUMI_raw.png"
      },
      content = function(file){
        png(file)
        print(density_Log10GenesPerUMI_raw_plot())
        dev.off()
      }
    )
    output$download_density_Log10GenesPerUMI_filtered <- downloadHandler(
      filename = function(){
        "density_Log10GenesPerUMI_filtered.png"
      },
      content = function(file){
        png(file)
        print(density_Log10GenesPerUMI_post_filter_plot())
        dev.off()
      }
    )
    output$download_violin_Log10GenesPerUMI_raw <- downloadHandler(
      filename = function(){
        "violin_Log10GenesPerUMI_raw.png"
      },
      content = function(file){
        png(file)
        print(violin_Log10GenesPerUMI_raw_plot())
        dev.off()
      }
    )
    output$download_violin_Log10GenesPerUMI_filtered <- downloadHandler(
      filename = function(){
        "violin_Log10GenesPerUMI_filtered.png"
      },
      content = function(file){
        png(file)
        print(violin_Log10GenesPerUMI_post_filter_plot())
        dev.off()
      }
    )
    output$download_UMI_vs_Genes_raw <- downloadHandler(
      filename = function(){
        "UMI_vs_Genes_raw.png"
      },
      content = function(file){
        png(file)
        print(UMI_vs_Genes_raw_plot())
        dev.off()
      }
    )
    output$download_UMI_vs_Genes_filtered <- downloadHandler(
      filename = function(){
        "UMI_vs_Genes_filtered.png"
      },
      content = function(file){
        png(file)
        print(UMI_vs_Genes_post_filter_plot())
        dev.off()
      }
    )
    output$download_raw_filtered_table <- downloadHandler(
      filename = function(){
        "raw_filtered_table.csv"
      },
      content = function(file){
        write.csv(raw_filtered_table_download(), file, quote = FALSE)
      }
    )
  }), input$evaluate)
}

# Run the application 
shinyApp(ui = ui, server = server)
