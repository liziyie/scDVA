# Load packages----
getPackage <- function(pkg, check = TRUE, load = TRUE, silent = FALSE, github = NULL) {
  if(check) {
    if(!suppressMessages(suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE)))) {
      if(is.null(github)){
        try(install.packages(pkg), silent = TRUE)
      }
      else{
        try(remotes::install_github(github))
      }
    }
  }
  if(load) suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE))
  if(load & !silent) message("Loaded ", pkg)
}
packages <- c("dplyr","psych","Matrix", "limma")
lapply(packages, getPackage)

# Read the dataset map file----
dataset_map <<- read.csv("./data/dataset_map.csv", row.names = 1, head = T, stringsAsFactors = F)

# Transfer the dataset map to a list tree----
datasetmap2list <- function(dataset_map){
  dataset_map <- dataset_map[dataset_map$DatasetName != "Initialize",]
  all_datasource <- unique(dataset_map$DatasetSource)
  dt <- NULL
  for(i in 1:length(all_datasource)){
    dt[all_datasource[i]] <- list(NULL)
    all_datasetname <- unique(dataset_map[dataset_map$DatasetSource == all_datasource[i], "DatasetName"])
    for(j in 1:length(all_datasetname)){
      dt[[i]][all_datasetname[j]] <- list(all_datasetname[j])
    }
    dt[[i]] <- structure(dt[[i]], stid=i, stopened=TRUE, stclass='project')
  }
  dt <- structure(dt, stopend = TRUE)
  return(dt)
}
dataset_tree <<- datasetmap2list(dataset_map)

# Main functions----
loadDataSet <- function(LoadData_selected) {
  ## Load the required dataset from the LoadData_selected parameter.
  ## 
  ## Args:
  #' @LoadData_selected: Record the dataset to load into the memory.
  ## 
  ## Returns:
  ## A dataframe with two columns. The first column is the dataset name. The second column 
  ## is the loading status, which could be "Failed", "Existed" or "Succeeded". 
  LoadData.status <- data.frame(Dataset = "Dataset", Status = "Status", stringsAsFactors = F)
  for(dataset in LoadData_selected){
    if(!exists(dataset_map[dataset,"Expression"])){
      expression.filepath <- sprintf("./data/%s.rda", dataset_map[dataset,"Expression"])
      metadata.filepath <- sprintf("./data/%s.rda", dataset_map[dataset,"Metadata"])
      if(file.exists(expression.filepath) & file.exists(metadata.filepath)){
        load(expression.filepath, envir = .GlobalEnv)
        load(metadata.filepath, envir = .GlobalEnv)
        LoadData.status <- rbind(LoadData.status, c(dataset, "Succeeded"))
      }else{
        LoadData.status <- rbind(LoadData.status, c(dataset, "Failed"))
      }
    }else{
      LoadData.status <- rbind(LoadData.status, c(dataset, "Existed"))
    }
  }
  LoadData.status <- LoadData.status[-1,]
  return(LoadData.status)
}

getPanelList <- function(SelectData_dataset) {
  ## Get the panel list from the dataset selected.
  ##
  ## Args:
  #'  @SelectData_dataset: Dataset used for visulization.
  ##
  ## Returns:
  ##  A list including Global_Cluster, Sub_Cluster, Tissue, Sample, Treatment and
  ## Day list. Sub_Cluster is also a list while the others are vectors.
  metadata_used <- get(dataset_map[SelectData_dataset, "Metadata"])
  SubCluster_list <- list()
  for (global_cluster in levels(as.factor(metadata_used[,"Global_Cluster"]))) {
    SubCluster_list[[global_cluster]] <-
      metadata_used %>% dplyr::filter(Global_Cluster == global_cluster) %>% pull(Sub_Cluster) %>% as.factor() %>% levels()
  }
  Tissue_list <- metadata_used %>% pull(Tissue) %>% as.factor() %>% levels()
  Sample_list <- metadata_used %>% pull(Sample) %>% as.factor() %>% levels()
  Treatment_list <- metadata_used %>% pull(Treatment) %>% as.factor() %>% levels()
  Day_list <- metadata_used %>% pull(Day) %>% as.factor() %>% levels()
  SubsetData_panel_list <- list(
    GlobalCluster = names(SubCluster_list),
    SubCluster_list = SubCluster_list,
    Tissue = Tissue_list,
    Sample = Sample_list,
    Treatment = Treatment_list,
    Day = Day_list
  )
  return(SubsetData_panel_list)
}

getGeneList <- function(SelectData_dataset, SelectData_normalization) {
  ## Get all the genes saved in the expression data from the platform input.
  ## 
  ## Args:
  #'  @SelectData_dataset: Data is used for visulization.
  #'  @SelectData_normalization: "tpm" or "counts".
  ## 
  ## Returns:
  ##  A vector of gene names in the expression data.
  GeneList <- row.names(get(dataset_map[SelectData_dataset,"Expression"])[[SelectData_normalization]])
  return(GeneList)
}

checkGeneList <- function(Genes, SelectData_dataset, SelectData_normalization) {
  ## Check whether all the genes input are in the expression data from the 
  ## platform input.
  ##
  ## Args:
  #'  @Genes: A long character including the genes input, should be separated by ",".
  #'  @SelectData_dataset: Data is used for visulization.
  #'  @SelectData_normalization: "tpm" or "counts".
  ##   
  ## Returns:
  ##  A list including the gene list info: genes input after removing the comma,
  ##  genes number, whethere all the genes are avaliable the expression data 
  ##  right genes name and the wrong genes name.
  Genes <- toupper(Genes)  # All the genes input are converted to uppercase.
  genes_input <- unlist(strsplit(Genes, ",| |\n|\t"))
  genes_input <- genes_input[genes_input != ""]
  all_genes_avaliable_flag = 1
  right_gene <- c()
  wrong_gene <- c()
  for (gene in genes_input) {
    if (!(gene %in% getGeneList(SelectData_dataset, SelectData_normalization))) {
      all_genes_avaliable_flag = 0
      wrong_gene <- c(wrong_gene, gene)
    }else{
      right_gene <- c(right_gene, gene)
    }
  }
  GeneListInfo <- list(
    genes_input = toupper(genes_input),
    all_genes_avaliable_flag = all_genes_avaliable_flag,
    gene_number = length(right_gene),
    right_gene = right_gene,
    wrong_gene = wrong_gene
  )
  return(GeneListInfo)
}

getPlotData <- function(Genes, 
                        SelectData_dataset, SelectData_normalization, 
                        GlobalCluster.selected, SubCluster.selected, 
                        Sample.selected, Tissue.selected, 
                        Treatment.selected, Day.selected) {
  ## Get the data frame used or plot.
  ## 
  ## Args:
  #'  @Genes: A long character including the genes input, should be separated 
  ##  by ",".
  #'  @SelectData_dataset: Data is used for visulization.
  #'  @SelectData_normalization: "tpm" or "counts"
  #'  @GlobalCluster.selected: A vector including the Global_Cluster to use.
  #'  @SubCluster.selected: A vector including the Sub_Cluster to use, the data used
  #'  was the intersection of GlobalCluster.selected and SubCluster.selected.
  #'  @Sample.selected: A vector including the patients to use.
  #'  @Tissue.selected: A vector including the tissues to use.
  #'  @Treatment.selected: A vector including the tretments to use.
  #'  @Day.selected: A vector including the days to use.
  ##   
  ## Returns:
  ##   A data frame including the gene expression data and patients' metadata. 
  ##   When multiple genes as input, the expression data was geometric mean of 
  ##   the tpm + 1. The expression for each gene was also in the return data.
  ## 
  Metadata_PlotData <- get(dataset_map[SelectData_dataset,"Metadata"]) %>% 
    filter(Sample %in% Sample.selected, Tissue %in% Tissue.selected, Global_Cluster %in% GlobalCluster.selected,Sub_Cluster %in% SubCluster.selected, Treatment %in% Treatment.selected, Day %in% Day.selected)
  row.names(Metadata_PlotData) <- Metadata_PlotData$CellName
  GenesInfo <- checkGeneList(Genes, SelectData_dataset, SelectData_normalization)
  Expression_PlotData <- as.matrix(get(dataset_map[SelectData_dataset,"Expression"])[[SelectData_normalization]][GenesInfo$right_gene, Metadata_PlotData$CellName])
  if (GenesInfo$gene_number > 1) {
    Mean_expression <- apply(Expression_PlotData + 1, 2,
                              psych::geometric.mean) - 1
    PlotData <- data.frame(
      Expression = Mean_expression, Metadata_PlotData,
      t(Expression_PlotData), check.names = F
    )
  }else if(GenesInfo$gene_number == 1){
    PlotData <- data.frame(
      Expression = Expression_PlotData, Metadata_PlotData, 
      input_gene = Expression_PlotData)
    names(PlotData)[names(PlotData) == "input_gene"] = GenesInfo$right_gene
  }
  return(PlotData)
}

getSigData <- function(Plot.data, Group.by, Per.cutoff) {
  ## Get the data frame used or gene significance calculation.
  ## 
  ## Args:
  #'  @Plot.data: The plot data created from getPlotData function.
  #'  @Group.by: To calculated the signifance between which groups, should be one
  #'  of a categorical variable in Plot.data, such as "Global_Cluster", 
  #'  "Sub_Cluster", "Patient" or "Tissue", "Treatment", "Day".
  #'  @Per.cutoff: The cutoff of a gene seen as expressed in the cells.
  ##   
  ## Returns:
  ##   A list including the result of gene significance calculation. Group_def is
  ##   the group name matched with the group number. Fvalue, coeff and pvalue is 
  ##   the ANOVA and tukeyHSD results. Percentage is the percentage of cells with
  ##   the gene expression in each group.
  ##   
 
  Expression.Per <- function(x) {
    return(sum(x > Per.cutoff) / length(x))
  }
  
  group.cha.id <- Plot.data[, Group.by]
  group.cha <- levels(as.factor(group.cha.id))
  ngroup <- length(group.cha)
  group_mapping <- data.frame(
    Group.cha = group.cha,
    Group.num = paste0("Grp", formatC(
      1:ngroup, width = 2, flag = "0"
    )),
    row.names = group.cha,
    stringsAsFactors = F
  )  # The group number id - group character id
  group.num.id <- group_mapping[as.character(group.cha.id), "Group.num"]
  temp <- data.frame(
    exprs = Plot.data$Expression,
    group = group.num.id,
    stringsAsFactors = F
  )  # Use the group number id
  exp.percent <- round(c(aggregate(temp$exprs, by = list(temp$group),  
                                   Expression.Per)$x), 3)
  exp.mean <- c(aggregate(temp$exprs, by = list(temp$group),  
                          mean)$x)
  exp.sd <- c(aggregate(temp$exprs, by = list(temp$group),  
                        sd)$x)
  row.names(group_mapping) <- group_mapping[, "Group.num"]
  group_mapping[, "Exp.percent"] <- exp.percent  # Expressed cells percentage
  group_mapping[, "Exp.mean"] <- exp.mean # Mean expression in each cluster
  group_mapping[, "Exp.sd"] <- exp.sd # Sd of the expression in each cluster
  
  if(length(unique(temp$group)) > 1){
    fml <- aov(exprs ~ group, data = temp)
    out.tukeyHSD <- TukeyHSD(fml)$group  # Do ANOVA
    Fvalue <- c(summary(fml)[[1]][1, 'F value'], summary(fml)[[1]][1, 'Pr(>F)'])
    names(Fvalue) <- c("F value", "Pr(>F)")
    
    coeff <- matrix(1, nrow = ngroup, ncol = ngroup)
    row.names(coeff) <- levels(as.factor(group.num.id))
    colnames(coeff) <- levels(as.factor(group.num.id))
    pvalue <- matrix(0, nrow = ngroup, ncol = ngroup)
    row.names(pvalue) <- levels(as.factor(group.num.id))
    colnames(pvalue) <- levels(as.factor(group.num.id))
    percentage <- matrix(diag(exp.percent), nrow = ngroup, ncol = ngroup)
    row.names(percentage) <- levels(as.factor(group.num.id))
    colnames(percentage) <- levels(as.factor(group.num.id))  # Initalize the matrix
    for (i in 1:nrow(out.tukeyHSD)) {
      lab1_name <- strsplit(row.names(out.tukeyHSD)[i], "-")[[1]][1]
      lab2_name <- strsplit(row.names(out.tukeyHSD)[i], "-")[[1]][2]
      coeff[lab1_name, lab2_name] <- out.tukeyHSD[i, "diff"]
      coeff[lab2_name, lab1_name] <- 0
      pvalue[lab1_name, lab2_name] <- out.tukeyHSD[i, "p adj"]
      pvalue[lab2_name, lab1_name] <- 1
    }  # Fill the matrix
    
    result <- list(
      group_def = group_mapping,
      Fvalue = Fvalue,
      coeff = coeff,
      pvalue = pvalue,
      percentage = percentage
    )
    return(result)
  }
}

getLIMMAData <- function(Plot.data, SelectData_dataset, SelectData_normalization, genes, x.cutoff, y.cutoff, group1, group2) {
  ## Get the data frame used or gene significance calculation.
  ## 
  ## Args:
  #'  @Plot.data: The plot data created from getPlotData function.
  #'  @genes: The genes used to do in-silico FACS.
  #'  @x.cutoff: The cutoff of x axis.
  #'  @y.cutoff: The cutoff of y axis.
  #'  @group1: The first group used to calculate DEGenes.
  #'  @group2: The second group used to calculate DEGenes.
  ##   
  ## Returns:
  ##   A table with the result of gene significance calculation. The differentially
  ##   expressed genes are calculated by LIMMA.
  ##   
  de.genes.all <- data.frame(Warning = "There should be two genes as input!", logFC = 0, adj.P.Val = 0)
  if(length(genes) == 2){
    gene.A <- genes[1]
    gene.B <- genes[2]
    Plot.data$FACS_groups <- "Group1"
    Plot.data[Plot.data[,gene.A] <= x.cutoff & Plot.data[,gene.B] <= y.cutoff, "FACS_groups"] <- "Group2"
    Plot.data[Plot.data[,gene.A] > x.cutoff & Plot.data[,gene.B] > y.cutoff, "FACS_groups"] <- "Group3"
    Plot.data[Plot.data[,gene.A] > x.cutoff & Plot.data[,gene.B] <= y.cutoff, "FACS_groups"] <- "Group4"
    de.genes.all <- data.frame(Warning = "Two groups should not be the same one!", logFC = 0, adj.P.Val = 0)
    if(group1 != group2){
      group1.used <- Plot.data %>% filter(FACS_groups == group1) %>% pull(CellName)
      if(sum(Plot.data$FACS_groups == group1) > 1000){
        group1.used <- sample(group1.used, 1000)
      }
      group2.used <- Plot.data %>% filter(FACS_groups == group2) %>% pull(CellName)
      if(sum(Plot.data$FACS_groups == group2) > 1000){
        group2.used <- sample(group2.used, 1000)
      }
      Plot.data <- Plot.data[Plot.data$CellName %in% c(group1.used, group2.used),]
      de.genes.all <- data.frame(Warning = "There should be at least one cell in each group!", logFC = 0, adj.P.Val = 0)
      if(length(unique(Plot.data$FACS_groups)) == 2){
        Expression_Data <- as.matrix(get(dataset_map[SelectData_dataset,"Expression"])[[SelectData_normalization]][, Plot.data$CellName])
        expression_matrix <- as.matrix(Expression_Data)
        expression_matrix <- expression_matrix[!duplicated(row.names(expression_matrix)),]
        groupid <- as.character(Plot.data$FACS_group)
        contrast <<- paste0(levels(factor(groupid)), collapse = "-")
        design <- model.matrix( ~ 0 + factor(groupid))
        colnames(design) <- levels(factor(groupid))
        rownames(design) <-
          colnames(expression_matrix)  # design data used in limma
        contrast.matrix <- makeContrasts(contrast, levels = design)
        fit <- lmFit(expression_matrix, design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        tempOutput <-  topTable(fit2, coef = 1, n = Inf)
        nrDEG <-  na.omit(tempOutput)
        nrDEG$Symbol <- row.names(nrDEG)
        positive_group <-
          row.names(fit2$contrasts)[fit2$contrasts == 1]  # high expression when logFC > 0
        negative_group <-
          row.names(fit2$contrasts)[fit2$contrasts == -1]  # low expression when logFC < 0
        nrDEG$Grp <-
          c(negative_group, positive_group)[as.numeric(nrDEG$logFC > 0) + 1]
        cell.Grp1 <- which(groupid == levels(as.factor(groupid))[1])
        cell.Grp2 <- which(groupid == levels(as.factor(groupid))[2])
        Exp.Mean.Grp1 <-
          rowMeans(expression_matrix[nrDEG$Symbol, cell.Grp1])  # mean expression in the group
        Exp.Mean.Grp2 <-
          rowMeans(expression_matrix[nrDEG$Symbol, cell.Grp2])  # mean expression out the group
        Exp.Per.Grp1 <-
          apply(expression_matrix[nrDEG$Symbol, cell.Grp1], 1, Expression.Per)  # expression percentage in the group
        Exp.Per.Grp2 <-
          apply(expression_matrix[nrDEG$Symbol, cell.Grp2], 1, Expression.Per)  # expression percentage out the group
        de.genes.all <- cbind(nrDEG, Exp.Mean.Grp1, Exp.Per.Grp1, Exp.Mean.Grp2, Exp.Per.Grp2)
        colnames(de.genes.all)[9:12] <-
          c(paste0(c("Exp.Mean.", "Exp.Per."), levels(as.factor(groupid))[1]), paste0(c("Exp.Mean.", "Exp.Per."), levels(as.factor(groupid))[2]))
        if (!is.na(de.genes.all[1, 1])) {
          de.genes.all <- de.genes.all %>% arrange(desc(logFC))
          row.names(de.genes.all) <- de.genes.all$Symbol
        }
      }
    }
  }
  return(de.genes.all)
}

Expression.Per <- function(x, cutoff = 0.1) {
  # percent of gene-expressed cell, cutoff = 0.1
  return(sum(x > cutoff) / length(x))
}