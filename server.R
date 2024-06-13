library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)

options(shiny.maxRequestSize = 3 * 1024^3)

server <- function(input, output, session) {
  # Inicializa shinyjs
  shinyjs::useShinyjs()
  
  observe({
    inFile1 <- input$file1
    inFile2 <- input$file2
    
    if (is.null(inFile1) || is.null(inFile2)) {
      return(NULL)
    }
  })
  
  # Observa el evento del botón reset y actualiza los fileInputs
  observeEvent(input$reset, {
    shinyjs::reset("file1")
    shinyjs::reset("file2")
  })
  
  output$uploadedFiles <- renderUI({
    if (!is.null(input$file1) || !is.null(input$file2)) {
      fileList1 <- if (!is.null(input$file1)) input$file1$name else NULL
      fileList2 <- if (!is.null(input$file2)) input$file2$name else NULL
      fileList <- c(fileList1, fileList2)
      HTML(paste("Files uploaded:", paste(fileList, collapse = ", ")))
    } else {
      HTML("No files uploaded.")
    }
  })
  
  observeEvent(input$qc, {
    inFile1 <- input$file1
    inFile2 <- input$file2
    
    if (is.null(inFile1) || is.null(inFile2)) {
      showNotification("Please upload both forward and reverse fastq files.", type = "error")
      return(NULL)
    }
    
    # Directorios de subida y salida
    upload_directory <- "/home/miquela/RNASeqAnalysis/Pipelines/arxius"
    output_directory <- "/home/miquela/RNASeqAnalysis/Pipelines/fastqc"
    
    # Verifica si los directorios existen y tienen permisos de escritura
    if (!dir.exists(upload_directory)) {
      dir.create(upload_directory, recursive = TRUE)
    }
    if (!dir.exists(output_directory)) {
      dir.create(output_directory, recursive = TRUE)
    }
    
    # Procesa los archivos subidos
    for (i in seq_along(inFile1$datapath)) {
      # Copia el archivo subido al directorio de subida
      file.copy(from = inFile1$datapath[i], to = file.path(upload_directory, inFile1$name[i]))
      # Ejecuta FastQC en el archivo copiado
      system(paste("fastqc", file.path(upload_directory, inFile1$name[i]), "-o", output_directory))
    }
    
    for (i in seq_along(inFile2$datapath)) {
      file.copy(from = inFile2$datapath[i], to = file.path(upload_directory, inFile2$name[i]))
      system(paste("fastqc", file.path(upload_directory, inFile2$name[i]), "-o", output_directory))
    }
    
    # Muestra una notificación de éxito
    showNotification("Quality control completed successfully.", type = "message")
    
    # Actualiza el reporte de QC
    output$qc_report <- renderUI({
      report_files <- list.files(output_directory, pattern = "*.html", full.names = TRUE)
      if (length(report_files) > 0) {
        tags$ul(
          lapply(report_files, function(file) {
            tags$li(tags$a(href = file, target = "_blank", basename(file)))
          })
        )
      } else {
        HTML("No QC reports generated.")
      }
    })
    
    
    
    # Ejecuta MultiQC para generar un reporte consolidado
    system(paste("python3 -m multiqc", output_directory, "-o", output_directory, "-f", sep = " "))
    
    # Muestra el reporte de MultiQC en la aplicación
    output$qc_report_initial <- renderUI({
      includeHTML(file.path(output_directory, "multiqc_report.html"))
    })
    
    showNotification("QC analysis completed.", type = "message")
  })
  
  output$fileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/arxius/", pattern = "*.fastq", full.names = FALSE)
    checkboxGroupInput("selectedFiles", "Select files:", choices = files)
  })
  
  observeEvent(input$run_sortmerna, {
    selectedFiles <- input$selectedFiles
    if (is.null(selectedFiles)) {
      showNotification("No files selected for SortMeRNA.", type = "error")
      return(NULL)
    }
    
    # Ejecuta sortmerna en los archivos seleccionados
    for (file in selectedFiles) {
      # Separa los archivos de lecturas 1 y 2
      if (grepl("_1.fastq$", file)) {
        read1 <- file
        read2 <- gsub("_1.fastq$", "_2.fastq", file)
        
        system(paste(
          "/home/miquela/bin/sortmerna",
          "--ref /media/disc_sda4/rRNA_databases/rfam-5.8s-database-id98.fasta",
          "--ref /media/disc_sda4/rRNA_databases/rfam-5s-database-id98.fasta",
          "--ref /media/disc_sda4/rRNA_databases/silva-arc-16s-id95.fasta",
          "--ref /media/disc_sda4/rRNA_databases/silva-arc-23s-id98.fasta",
          "--ref /media/disc_sda4/rRNA_databases/silva-bac-16s-id90.fasta",
          "--ref /media/disc_sda4/rRNA_databases/silva-bac-23s-id98.fasta",
          "--ref /media/disc_sda4/rRNA_databases/silva-euk-18s-id95.fasta",
          "--ref /media/disc_sda4/rRNA_databases/silva-euk-28s-id98.fasta",
          "--reads", paste0("/home/miquela/RNASeqAnalysis/Pipelines/arxius/", read1),
          "--workdir", paste0("/home/miquela/RNASeqAnalysis/Pipelines/remove_rrna/", read1, "/"),
          "--fastx",
          "--a 12",
          "--other",
          "--out2",
          "--paired_in"
        ))
      }
    }
  })
  
  output$trimmingFileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/trimming/", pattern = "*_[12].fastq", full.names = FALSE)
    checkboxGroupInput("selectedTrimmingFiles", "Select files:", choices = files)
  })
  
  observeEvent(input$run_trimming, {
    selectedFiles <- input$selectedTrimmingFiles
    if (is.null(selectedFiles)) {
      return(NULL)
    }
    
    # Ejecuta fastp en los archivos seleccionados
    for (file in selectedFiles) {
      # Separa los archivos de lecturas 1 y 2
      if (grepl("_1.fastq$", file)) {
        read1 <- file
        read2 <- gsub("_1.fastq$", "_2.fastq", file)
        
        system(paste(
          "fastp -i", paste0("/home/miquela/RNASeqAnalysis/Pipelines/trimming/", read1),
          "-I", paste0("/home/miquela/RNASeqAnalysis/Pipelines/trimming/", read2),
          "-o", paste0("/home/miquela/RNASeqAnalysis/Pipelines/trimming/trimmed_", read1),
          "-O", paste0("/home/miquela/RNASeqAnalysis/Pipelines/trimming/trimmed_", read2)
        ))
      }
    }
  })
  
  output$alignmentFileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/arxius/", pattern = "*_[12].fastq", full.names = FALSE)
    checkboxGroupInput("selectedAlignmentFiles", "Select files:", choices = files)
  })
  
  output$rrnaRemovedFileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/remove_rrna/", pattern = "*_[12].fastq", full.names = FALSE)
    checkboxGroupInput("selectedRrnaRemovedFiles", "Select files:", choices = files)
  })
  
  output$download_rrna <- downloadHandler(
    filename = function() {
      return(paste("selected_files.zip"))
    },
    content = function(file) {
      # Obtenemos los archivos seleccionados
      selectedFiles <- input$selectedRrnaRemovedFiles
      if (is.null(selectedFiles)) {
        return(NULL)
      }
      
      # Definimos la ruta completa de los archivos seleccionados
      filePaths <- file.path("/home/miquela/RNASeqAnalysis/Pipelines/remove_rrna/", selectedFiles)
      
      # Creamos un archivo zip temporal
      zip(file, files = filePaths)
    }
  )
  
  
  output$trimmedFileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/trimming/", pattern = "*_[12].fastq", full.names = FALSE)
    checkboxGroupInput("selectedTrimmedFiles", "Select files:", choices = files)
  })
  
  observeEvent(input$run_alignment, {
    selectedFiles <- c(input$selectedAlignmentFiles, input$selectedRrnaRemovedFiles, input$selectedTrimmedFiles)
    if (is.null(selectedFiles)) {
      return(NULL)
    }
    
    # Ejecuta STAR en los archivos seleccionados
    for (file in selectedFiles) {
      # Separa los archivos de lecturas 1 y 2
      if (grepl("_1.fastq$", file)) {
        read1 <- file
        read2 <- gsub("_1.fastq$", "_2.fastq", file)
        
        # Determina la carpeta de origen del archivo
        if (file %in% input$selectedAlignmentFiles) {
          folder <- "/home/miquela/RNASeqAnalysis/Pipelines/arxius/"
        } else if (file %in% input$selectedRrnaRemovedFiles) {
          folder <- "/home/miquela/RNASeqAnalysis/Pipelines/remove_rrna/"
        } else if (file %in% input$selectedTrimmedFiles) {
          folder <- "/home/miquela/RNASeqAnalysis/Pipelines/trimming/"
        }
        
        system(paste(
          "STAR --runThreadN 10 --genomeDir /home/miquela/RNASeqAnalysis/genome/STAR",
          "--sjdbGTFfile /home/miquela/RNASeqAnalysis/genome/annotation.gtf",
          "--readFilesIn", paste0(folder, read1),
          paste0(folder, read2),
          "--outFileNamePrefix", paste0("/home/miquela/RNASeqAnalysis/Pipelines/alignment/aligned_"),
          "--outSAMtype BAM SortedByCoordinate",
          "--quantMode GeneCounts"
        ))
      }
      
    }
  })
  
  output$dupFileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/alignment", pattern = "*.bam", full.names = FALSE)
    checkboxGroupInput("selectedDupFiles", "Select files:", choices = files)
  })
  
  observeEvent(input$run_remove_dups, {
    selectedFiles <- input$selectedDupFiles
    if (is.null(selectedFiles)) {
      return(NULL)
    }
    
    # Ejecuta samtools rmdup en los archivos seleccionados
    for (file in selectedFiles) {
      # Define la ruta completa del archivo
      filepath <- paste0("/home/miquela/RNASeqAnalysis/Pipelines/alignment/", file)
      
      # Ordenar por nombre
      system(paste("samtools sort -n", filepath, "-o", paste0(filepath, ".sorted_name.bam")))
      
      # Aplicar samtools fixmate
      system(paste("samtools fixmate -m", paste0(filepath, ".sorted_name.bam"), "temp.bam"))
      
      # Ordenar por coordenadas
      system(paste("samtools sort temp.bam -o sorted_temp.bam"))
      
      # Eliminar duplicados con samtools markdup
      system(paste("samtools markdup -r sorted_temp.bam", paste0(filepath, ".dedup.bam")))
    }
  })
  
  output$indexFileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/alignment", pattern = "*.bam", full.names = FALSE)
    checkboxGroupInput("selectedDupFiles", "Select files:", choices = files)
  })
  
  observeEvent(input$run_index, {
    # Aquí iría el código para ejecutar samtools index
  })
  
  output$quantFileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/alignment", pattern = "*.bam", full.names = FALSE)
    checkboxGroupInput("selectedDupFiles", "Select files:", choices = files)
  })
  
  output$resultsFileCheckboxes <- renderUI({
    files <- list.files(path = "/home/miquela/RNASeqAnalysis/Pipelines/featurecounts", pattern = "*.txt", full.names = FALSE)
    checkboxGroupInput("selectedResultsFiles", "Select files:", choices = files)
  })
  
  output$download_results <- downloadHandler(
    filename = function() {
      return(paste("selected_results.zip"))
    },
    content = function(file) {
      selectedFiles <- input$selectedResultsFiles
      if (is.null(selectedFiles)) {
        return(NULL)
      }
      
      filePaths <- file.path("/home/miquela/RNASeqAnalysis/Pipelines/featurecounts", selectedFiles)
      
      zip(file, files = filePaths)
    }
  )
  
  
  observeEvent(input$run_quant, {
    selectedFiles <- input$selectedDupFiles
    if (is.null(selectedFiles)) {
      return(NULL)
    }
    
    # Define la ruta completa de los archivos seleccionados
    filePaths <- file.path("/home/miquela/RNASeqAnalysis/Pipelines/alignment", selectedFiles)
    
    # Ejecuta featureCounts en los archivos seleccionados
    for (file in filePaths) {
      system(paste(
        "featureCounts -T 10",
        "-a /home/miquela/RNASeqAnalysis/genome/annotation.gtf",
        "-o /home/miquela/RNASeqAnalysis/Pipelines/featurecounts/featurecounts.txt",
        file,
        "-s 1 -p -M -t exon -g gene_id"
      ))
    }
  })
  
  
  observeEvent(input$run_qc_rrna, {
    selectedFiles <- input$selectedFiles
    if (is.null(selectedFiles)) {
      return(NULL)
    }
    
    # Define la ruta de salida para los archivos FastQC
    output_directory <- "/home/miquela/RNASeqAnalysis/Pipelines/fastqc_rrna"
    
    # Ejecuta FastQC en los archivos seleccionados
    for (file in selectedFiles) {
      system(paste("fastqc", paste0("/home/miquela/RNASeqAnalysis/Pipelines/remove_rrna/", file), "-o", output_directory))
    }
    
    # Ejecuta MultiQC en los resultados de FastQC
    system(paste("python3 -m multiqc", output_directory, "-o", output_directory, "-f", sep=" "))
    
    # Muestra el informe de MultiQC
    output$qc_report_rrna <- renderUI({
      includeHTML(paste0(output_directory, "/multiqc_report.html"))
    })
  })
  
  observeEvent(input$run_qc_trimming, {
    selectedFiles <- input$selectedTrimmingFiles
    if (is.null(selectedFiles)) {
      return(NULL)
    }
    
    # Define la ruta de salida para los archivos FastQC
    output_directory <- "/home/miquela/RNASeqAnalysis/Pipelines/fastqc_trimming"
    
    # Ejecuta FastQC en los archivos seleccionados
    for (file in selectedFiles) {
      system(paste("fastqc", paste0("/home/miquela/RNASeqAnalysis/Pipelines/trimming/", file), "-o", output_directory))
    }
    
    # Ejecuta MultiQC en los resultados de FastQC
    system(paste("python3 -m multiqc", output_directory, "-o", output_directory, "-f", sep=" "))
    
    # Muestra el informe de MultiQC
    output$qc_report_trimming <- renderUI({
      includeHTML(paste0(output_directory, "/multiqc_report.html"))
    })
  })
  
}

