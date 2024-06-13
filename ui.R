library(shiny)
library(shinydashboardPlus)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  skin="purple",
  
  dashboardHeader(title = "GenoFlow"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("File upload", tabName = "upload", icon = icon("upload")),
      menuItem("Preprocessing", tabName = "preprocessing", icon = icon("edit")),
      menuItem("RNA-Seq", tabName = "rna_seq", icon = icon("dna")),
      menuItem("SNV and indels", tabName = "snv_indels", icon = icon("edit")),
      menuItem("Methylation", tabName = "methylation", icon = icon("dna"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),  # Asegurarse de que shinyjs esté cargado en la UI
    tabItems(
      tabItem(tabName = "home",
              h2("Welcome to GenoFlow App for Sequences Data Analysis"),
              p("Here you can upload fastq files and do sequence data analysis with different pipelines. 
                Start uploading your files at the File Upload tab."),
              tags$img(src = "geno_flow.png", height = "auto", width = "50%")
      ),
      
      tabItem(tabName = "upload",
              h2("Upload raw files"),
              p("Here you can upload your fastq files to start the analysis. We accept .fq, .fastq, .fq.gz and .fastq.gz"),
              fileInput("file1", "Upload your forward fastq files here:",
                        accept = c(
                          ".txt",
                          ".fastq",
                          ".fastq.gz",
                          ".fq",
                          ".fq.gz"
                        ), multiple = TRUE),
              fileInput("file2", "Upload your reverse fastq files here:",
                        accept = c(
                          ".txt",
                          ".fastq",
                          ".fastq.gz",
                          ".fq",
                          ".fq.gz"
                        ), multiple = TRUE),
              
              actionButton("goToPreprocesssing", "Go to Preprocessing"),
              actionButton("reset", "Reset Files"),
              actionButton("qc", "Quality Control"),
              tabsetPanel(
                tabPanel("Quality Control Results",
                         htmlOutput("qc_report_initial")
                )
              )
      ),
      tabItem(tabName = "preprocessing",
       h2("Preprocessing"),
       p("In the preprocessing of RNA-Seq, DNA variant, or methylation data, quality control (QC) and trimming are pivotal steps. QC assesses the data for issues like sequencing errors, which can significantly distort analytical results. Trimming removes low-quality bases and sequencing adapters, improving the reliability of downstream analyses. These steps are crucial as they ensure the accuracy of variant calling, gene expression analysis, and methylation pattern identification, ultimately leading to more valid biological conclusions.
                Integrating quality control and trimming into a single tab in a data analysis pipeline offers a streamlined and efficient workflow. It allows for immediate identification and correction of sequencing issues, ensuring only high-quality data progresses to downstream analysis. This consolidation simplifies the process, reduces the risk of errors in data handling, and accelerates the transition from raw data to actionable insights. Moreover, it enhances user experience by providing a cohesive interface for these critical initial steps."),
       actionButton("go_RNAseq", "Go to RNA-Seq pipeline"),
       actionButton("go_SNV", "Go to SNV and indels pipeline"),
       actionButton("go_methylation", "Go to Methylation pipeline"),
       tabsetPanel(
         tabPanel("Quality control",
                  h2("Quality Control of Raw Data"),
                  p("The quality control step using FastQC involves analyzing raw sequencing data to identify potential issues such as poor quality scores and adapter contamination. MultiQC aggregates these results across multiple samples, providing a comprehensive report that highlights areas needing attention. Together, they form an essential part of the preprocessing in sequencing workflows, ensuring data integrity before further analysis."),
                  actionButton("qc", "Quality Control"),
                  tabsetPanel(
                    tabPanel("Quality Control Results",
                             htmlOutput("qc_report_prepro")
                             )
                  )),
         tabPanel("Trimming",
                  h2("Trimming"),
                  p("The trimming step with fastp is an essential process that trims low-quality bases and removes adapters from sequencing reads. It enhances data quality by correcting mismatches and filtering out sequences that are too short or of low complexity. Fastp’s efficiency and speed streamline the preprocessing of large datasets, making it a valuable tool in next-generation sequencing workflows."),
                  actionButton("trimming", "Trimming"),
                  downloadButton("download_trimming1", "Download trimmed files")),
         
         tabPanel("Alignment",
                  sidebarPanel(
                    radioButtons("genomeChoice", "Select reference genome for Alignment",
                                 choices = list("Homo sapiens" = "Homo",
                                                "Mus musculus" = "Mus",
                                                "Rattus norvegicus" = "Rattus"))
                  ),
                  mainPanel(
                    h2("Alignment"),
                    p("In the context of a Next-Generation Sequencing (NGS) data analysis pipeline, aligning to a reference genome refers to the process of mapping the short DNA sequence reads obtained from NGS to a known reference genome sequence. This step is crucial because it determines the origin of each read by comparing it to the reference, allowing for the identification of variations from the reference sequence."),
                    p("STAR and BWA are tailored for RNA-Seq and SNV detection because they efficiently handle splicing and short DNA sequences, respectively. Bismark is preferred for methylation analysis as it specifically aligns bisulfite-converted reads, which is not a strength of general aligners like STAR or BWA. For this reason, if you are working with methylation data, please go to Methylation tab and use Bismark aligner."),
                    actionButton("align_STAR", "Align with STAR"),
                    actionButton("align_BWA", "Align with BWA"),
                    actionButton("align_bismark", "Align with Bismark (for Methylation)")
                  )
         ),
         tabPanel(
           "Preprocessing Image",
           tags$img(src = "Preprocessing.png", height = "auto", width = "75%")
         )
       )),
      
      tabItem(tabName="rna_seq",
              h2("RNA-Seq Data Analysis Pipeline"),
              p("An RNA-Seq Data Analysis Pipeline is a series of automated steps that cleans, organizes, and analyzes the massive amount of data generated by RNA-Sequencing. This pipeline is crucial for researchers to efficiently unlock the secrets hidden within gene expression levels."),
              tabsetPanel(
                tabPanel("Remove rRNA",
                         h2("Remove rRNA"),
                         p("Ribosomal RNA (rRNA) depletion is a critical step in RNASeq data analysis pipelines. It involves the removal of highly abundant rRNA species from the sample, which allows for the efficient detection of functionally relevant coding as well as non-coding transcripts. The importance of rRNA depletion lies in the fact that in eukaryotic cells, there are about 10 times more molecules of rRNA than all other types of RNA combined. If total RNA is sequenced without rRNA depletion, about 90% of RNA-seq reads will be ribosomal, leaving only ~10% of reads relevant to your experiment. Here you can select your uploaded files and remove rRNA from selected files"),
                         uiOutput("fileCheckboxes"),
                         actionButton("run_sortmerna", "Run SortMeRNA"),
                         actionButton("run_qc_rrna", "Run Quality Control"),
                         downloadButton("download_rrna", "Download selected files"),
                         tabsetPanel(
                           tabPanel("Quality Control Results",
                                    htmlOutput("qc_report_rrna"))
                         )),
                
                tabPanel("Trimming",
                         h2("Trimming"),
                         p("Trimming is another crucial step in RNASeq data analysis pipelines. It involves the removal of low-quality bases and adapter sequences from the sequence reads. This process is often referred to as quality control. The importance of trimming lies in its ability to improve the accuracy of downstream analyses. Low-quality bases and adapter sequences can introduce errors in read alignment and gene expression estimation. By removing these, trimming increases the quality of the reads, thereby improving the accuracy of the alignment and the subsequent analyses"),
                         uiOutput("trimmingFileCheckboxes"),
                         actionButton("run_trimming", "Run Trimming with fastp"),
                         actionButton("run_qc_trimming", "Run Quality Control"),  
                         downloadButton("download_trimming", "Download selected files"),
                         tabsetPanel(
                           tabPanel("Quality Control Results",
                                    htmlOutput("qc_report_trimming"))
                           )
                         ),
                tabPanel("Alignment",
                         sidebarPanel(
                           radioButtons("genomeChoice2", "Select reference genome for Alignment:",
                                        choices = list("GRCh38 (Human)" = "Homo",
                                                       "GRCm38 (Mus musculus)" = "Mus",
                                                       "Rattus norvegicus" = "Rattus"))
                         ),
                         mainPanel(
                           h2("Alignment"),
                           p("Alignment in RNASeq data analysis pipelines is the process of mapping short RNA sequencing reads to a reference genome or transcriptome. This step is crucial as it determines the origin of each read in the genome, which is essential for subsequent analyses such as quantification of gene and transcript levels, differential gene expression, alternative splicing, and others. Therefore, accurate alignment is key to obtaining reliable results from RNASeq data"),
                           uiOutput("alignmentFileCheckboxes"),
                           uiOutput("rrnaRemovedFileCheckboxes"),
                           uiOutput("trimmedFileCheckboxes"),
                           actionButton("run_alignment", "Run Alignment with STAR"),
                           actionButton("run_alignment2", "Run Alignment with BWA"),
                           downloadButton("download_alignment", "Download selected files")
                         )
                         
                ),
                tabPanel("Remove dups",
                         h2("Remove duplicates"),
                         p("Post-alignment duplicate removal, also known as deduplication, is a step in RNASeq data analysis pipelines where identical sequence reads that map to the same location in the genome are identified and removed. These duplicates often arise from PCR amplification during library preparation, which can lead to overrepresentation of certain reads. The importance of this step lies in its ability to reduce bias in downstream analyses. By removing these duplicates, we can obtain a more accurate representation of the original RNA population. However, it’s important to note that this process should be applied carefully, as overzealous removal of duplicates can also introduce bias and affect the results"),
                         uiOutput("dupFileCheckboxes"),
                         actionButton("run_remove_dups", "Run Remove duplicates with samtools"),
                         downloadButton("download_removedups", "Download selected files")
                ),
                tabPanel("Index",
                         h2("Index data"),
                         p("Indexing, in the context of RNASeq data analysis pipelines, is a step that follows alignment. It involves creating an index or a searchable database from the alignment results, which typically include a large file of aligned sequencing reads. The importance of indexing lies in its ability to facilitate efficient access to the alignment results. This is particularly useful for downstream analyses such as variant calling, differential expression analysis, and others, which often require random access to the alignment data. By creating an index, these analyses can quickly retrieve the alignment information for specific genomic regions without having to scan the entire alignment file"),
                         uiOutput("indexFileCheckboxes"),
                         actionButton("run_index", "Run Index with samtools"),
                         downloadButton("download_index", "Download selected files")
                ),
                tabPanel("Quantification",
                         h2("Quantification - count matrix"),
                         p("Quantification using featureCounts is a step in RNASeq data analysis pipelines where the number of reads that align to each gene or other genomic features (such as exons, introns, or UTRs) are counted. This is often referred to as “counting features” and the output is a count matrix, with genes as rows and samples as columns. The importance of this step lies in its ability to measure gene expression. The number of reads that align to a gene is proportional to the gene’s expression level in the sample. Therefore, by counting the number of reads that align to each gene, we can obtain a measure of gene expression that can be used for downstream analyses such as differential expression analysis"),
                         uiOutput("quantFileCheckboxes"),
                         actionButton("run_quant", "Run FeatureCounts"),
                         downloadButton("download_quant", "Download selected files"),
                         tabsetPanel(
                           tabPanel("Results",
                                    uiOutput("resultsFileCheckboxes"),
                                    downloadButton("download_results", "Download selected files"))
                         )
                )
                
                
                
                )
                ),
      tabItem(tabName= "methylation",
              h2("Methylation Data Analysis Pipeline"),
              p("A Methylation Data Analysis Pipeline automates processing DNA methylation data, like finding differentially methylated regions. This helps researchers understand how gene activity is regulated beyond just gene expression."),
              tabsetPanel(
                tabPanel("Alignment",
                         h2("Alignment"),
                         p("Alignment with Bismark is a crucial step in a Methylation Data Analysis Pipeline, where bisulfite-treated sequencing reads are aligned to a reference genome. This process identifies methylated and unmethylated cytosines, providing detailed insights into DNA methylation patterns at single-base resolution. It’s essential for understanding the epigenetic roles of DNA methylation in cell differentiation, development, and disease"),
                         actionButton("run_alignment_bismark", "Run Bismark for Alignment"),
                         downloadButton("download_alignment_bismark", "Download selected files")),
                tabPanel("Methylation Calling",
                         h2("Methylation extraction or methylation calls with Bismark"),
                         p("Methylation extraction with Bismark is a step where methylation information is extracted from aligned sequencing reads. This process discerns the methylation status of cytosines in different contexts (CpG, CHG, CHH), which is pivotal for assessing DNA methylation patterns. It’s a key component in understanding epigenetic modifications and their implications in gene regulation and disease states within a Methylation Data Analysis Pipeline"),
                         actionButton("run_extraction", "Run Methylation Calling"),
                         downloadButton("download_extraction", "Download selected files")),
                
                tabPanel(
                  "Methylation Pipeline Image",
                  tags$img(src = "Pipeline_metilacio.drawio.png", height = "auto", width = "75%")
                )
              )),
      
      
      tabItem(tabName = "snv_indels",
              h2("SNV and indels Detection Pipeline"),
              p("An SNV and Indels Data Analysis Pipeline processes DNA sequencing data to identify Single Nucleotide Variations (SNVs) and insertions/deletions (Indels) in the genome. This pinpoints genetic variations that may be linked to diseases, traits, or individual differences."),
              tabsetPanel(
                tabPanel("Mapping",
                         h2("Mapping"),
                         p("The mapping step aligns short DNA sequencing reads to a reference genome. It's crucial because accurate alignment allows for pinpointing true SNVs and Indels from sequencing errors. BWA-MEM is a popular mapping tool known for its speed and accuracy in handling complex regions like repetitive DNA."),
                         actionButton("run_mapping_snv", "Run BWA-MEM Mapping for SNV"),
                         downloadButton("download_mapped_snv", "Download selected files")),
                
                tabPanel("Convert to BAM, sort, remove duplicates and index",
                         h2("Samtools for converting to BAM file, sort BAM file and remove duplicates"),
                         p("Samtools acts as a multi-tool in the SNV/Indel pipeline, streamlining data preparation for accurate variant calling. It efficiently converts raw sequencing data (SAM) to the compact BAM format, saving storage space. Additionally, Samtools sorts the BAM file by read name or position, which is crucial for downstream variant callers to operate correctly."),
                         actionButton("run_samtools_snv", "Run Samtools"),
                         downloadButton("download_bam_snv", "Download selected files")),
                
                tabPanel("Call Variants",
                         h2("Call Variants with samtools"),
                         p("The 'Call Variants' step with Samtools identifies genetic variations such as SNVs and Indels in the aligned sequencing data. This is essential for pinpointing specific genetic differences that may be associated with various traits or diseases."),
                         actionButton("run_variant_calling", "Call Variants"),
                         downloadButton("download_variants", "Download selected files")),
                
                tabPanel(
                  "SNV and indels Pipeline Image",
                  tags$img(src = "SNV_indel_pipeline.drawio.png", height = "auto", width = "75%")
                )
              ))
      
      
      
      
              ))
      
  

      
      
      
 
)

