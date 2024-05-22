#### Result 5 code for validation ####

library(ggpubr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(ggrepel)
library(randomcoloR)
library(dplyr)
library(edgeR)
library(DESeq2)
library(data.table)
library(ggalt)
library(stringr)
library(tidyr)
library(ggradar)

#### Analysis for finding sex-biased genes in embryonic stage E4.5, E6.5, E8.5, E10.5 ####
meta_data <- fread("./Data/Figure5_data/Metadata.txt")
variable <- unique(meta_data$Stage[which(meta_data$Seq_type=="Single")])


for(i in 1:length(variable)){
  
  ## Load Pseudo bulk data ##
  Bulk_data <- fread("./Data/Quantification_chicken_embryo/Pseudo_bulk_gene_expression_data_galgal6.tsv")
  
  ## Loading featureCount data
  meta_data <- fread("./Data/Figure5_data/Metadata.txt")
  
  meta_data <- meta_data[which(meta_data$Age=="Embryo" & meta_data$Species=="Chicken" & meta_data$Stage==variable[i]),]
  
  
  #### Extracting orthologue gene list ####
  
  orthologue_list <- fread("./Data/Figure5_data/Orthologue_filtered_galgal6.txt")
  
  Bulk_data <- Bulk_data[which(Bulk_data$geneSymbol %in% orthologue_list$Chicken_geneSymbol),]
  
  
  Annotation <- fread("./Data/Figure5_data/mart_export_galgal6.txt")
  
  
  total_table <- merge(Annotation,
                       Bulk_data,
                       by.x="Gene stable ID",
                       by.y="geneSymbol")
  
  Annotation <- total_table %>% dplyr::select(c(colnames(Annotation)))
  
  geneSymbol <- Annotation$`Gene stable ID`
  
  
  count <- total_table %>% dplyr::select(meta_data$Run)
  
  
  unique(Annotation$`Chromosome/scaffold name`)
  ## Removing unnecessary chromosome
  indexRemoval <- c(grep("KZ",Annotation$`Chromosome/scaffold name`),grep("MT",Annotation$`Chromosome/scaffold name`))
  
  
  ## Removing all zero counted genes    
  count <- count[-indexRemoval, ]
  Annotation <- Annotation[-indexRemoval, ]
  geneSymbol <- geneSymbol[-indexRemoval]
  
  
  
  ## Defining experimental variables
  
  meta_data$sex <- str_to_sentence(meta_data$sex)
  meta_data$Tissue <- str_to_sentence(meta_data$Tissue)
  Tissue <- factor(meta_data$Tissue)
  Stage <- factor(meta_data$Stage)
  Sex <- factor(meta_data$Sex)
  Seq_type <- factor(meta_data$Seq_type)
  
  
  design <- model.matrix( ~ Sex)
  
  
  ## TMM normalization using edgeR package
  Y <- DGEList(counts = count, genes = geneSymbol)
  
  Y <- calcNormFactors(Y, method = "TMM") ## TMM Normalization
  
  Y <- estimateDisp(Y, design)
  
  fit <- glmFit(Y, design)## GLM
  
  logCPM <- edgeR::cpm(Y, normalized.lib.sizes = TRUE, log = T)
  
  
  Result_sex <- glmLRT(fit, 2)
  Result_sex_table <-  topTags(Result_sex , n = dim(logCPM)[1], sort.by = "none")$table
  
}


#### Analysis for finding sex-biased genes in embryonic stage E5.5, E6, E10, E18.5, E19 ####

meta_data <- fread("./Data/Figure5_data/Metadata.txt")
variable <- unique(meta_data$Stage[which(meta_data$Seq_type=="Bulk" & meta_data$Age=="Embryo" & meta_data$Species=="Chicken")])

variable

for(k in 1:length(variable)){
  
  meta_data <- fread("./Data/Figure5_data/Metadata.txt")
  
  meta_data <- meta_data[which(meta_data$Age=="Embryo" & meta_data$Species=="Chicken"),]
  
  meta_data <- meta_data[which(meta_data$Stage==variable[k]),]
  
  ### Extracting orthologue gene list ###
  
  orthologue_list <- fread("./Data/Figure5_data/Orthologue_filtered_galgal6.txt")
  
  fileList <- paste0("./Data/Quantification_chicken_embryo/",meta_data$Run,'.txt')
  
  
  fileList
  temp <- fread(fileList[1])
  dim(temp)
  str(temp)
  sampleID <- c(fileList[1])
  
  
  geneSymbol <- temp$Geneid
  
  orthologue_idx <- which(geneSymbol %in% orthologue_list$Chicken_geneSymbol)
  
  geneSymbol <- geneSymbol[orthologue_idx]
  
  
  
  count <- data.frame(temp[, 7])
  
  
  for (i in 2:length(fileList)) {
    temp <- data.frame(fread(fileList[i]))
    count <- data.frame(count, temp[, 7])
    sampleID <- c(sampleID,fileList[i])
  }
  
  
  count <- count[orthologue_idx,]
  
  
  ## Loading gene annotation
  temp <- data.frame(geneSymbol)
  
  
  Annotation <- fread("./Data/Figure5_data/mart_export_galgal6.txt")
  temp <-  merge(
    x = temp,
    y = Annotation,
    by.x = "geneSymbol",
    by.y = "Gene stable ID",
    sort = F
  )
  
  Annotation <- temp
  
  all.equal(as.character(Annotation$geneSymbol), geneSymbol)
  
  data.frame(as.character(Annotation$geneSymbol), geneSymbol)
  
  
  
  
  ## Removing unnecessary chromosome
  
  indexRemoval <- c(grep("KZ",Annotation$`Chromosome/scaffold name`),grep("MT",Annotation$`Chromosome/scaffold name`),grep("LGE",Annotation$`Chromosome/scaffold name`),grep("MT",Annotation$`Chromosome/scaffold name`),grep("MU",Annotation$`Chromosome/scaffold name`),grep("JAENSK",Annotation$`Chromosome/scaffold name`))
  # 
  
  unique(Annotation$`Chromosome/scaffold name`)
  
  ## Removing all zero counted genes    
  
  indexRemoval <- c(indexRemoval,(which(rowSums(count) == 0)))
  count <- count[-indexRemoval, ]
  Annotation <- Annotation[-indexRemoval, ]
  geneSymbol <- geneSymbol[-indexRemoval]
  
  
  
  ## Defining experimental variables
  
  meta_data$sex <- str_to_sentence(meta_data$sex)
  meta_data$Tissue <- str_to_sentence(meta_data$Tissue)
  Tissue <- factor(meta_data$Tissue)
  Stage <- factor(meta_data$Stage)
  Sex <- factor(meta_data$Sex)
  Seq_type <- factor(meta_data$Seq_type)
  
  
  design <- model.matrix( ~ Sex)
  
  
  
  
  
  ## TMM normalization using edgeR package
  Y <- DGEList(counts = count, genes = geneSymbol)
  
  Y <- calcNormFactors(Y, method = "TMM") ## TMM Normalizatio
  
  Y <- estimateDisp(Y, design)
  
  fit <- glmFit(Y, design)## GLM
  
  logCPM <- edgeR::cpm(Y, normalized.lib.sizes = TRUE, log = T)
  
  
  #### DEG analysis in sex
  
  colnames(fit)
  
  Result_sex <- glmLRT(fit, 2)
  Result_sex_table <-  topTags(Result_sex , n = dim(logCPM)[1], sort.by = "none")$table
  
  
  
  
}



#### Figure 5A ####


load("./Data/Figure5_data/Figure5A.RData")
MDS_data <- plotMDS(Y)


plot_data <- data.frame(ID=meta_data$Run,Tissue:Sex,Tissue,Sex,Stage,instrument,layout,Tissue:Stage,Tissue:Stage:Seq_type,Seq_type:Stage:Tissue,Stage:Tissue, X = MDS_data$x, Y = MDS_data$y)


plot_data$Stage_factor <- recode_factor(plot_data$Stage,
                                        "E4.5"="E4.5-E6",
                                        "E5.5"="E4.5-E6",
                                        "E6"="E4.5-E6",
                                        "E8.5"="E8.5-E19",
                                        "E10"="E8.5-E19",
                                        "E10.5"="E8.5-E19",
                                        "E18.5"="E8.5-E19",
                                        "E19"="E8.5-E19")


MDS_PLOT <- ggplot(plot_data, aes(x = X, y = Y)) +
  geom_encircle(aes(group = Stage_factor, color=Stage_factor), s_shape = 2, expand = 0.1,size=3, linetype="dashed") +
  geom_jitter(size = 8.5,aes(fill=Tissue),alpha=0.7, shape=21, stroke=NA) +
  theme_classic() +
    coord_cartesian(xlim = c(min(plot_data$X)-0.035,max(plot_data$X)+0.035),
                   ylim = c(min(plot_data$Y)-0.05,max(plot_data$Y)+0.05)) +
  scale_fill_manual(values=c("#F8BF28","#38629A")) +
  scale_color_manual(values=c("#34A853","#EA4335","#f5a002")) +
  theme(legend.title = element_blank(),
        legend.position = "None",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16,color="black"),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = NA)) 

MDS_PLOT

ggsave("Figure5A.pdf",
       plot = MDS_PLOT,
       height = 5,
       width = 5)


#### Figure 5B ####
load("./Data/Figure5_data/Figure5B.RData")




MDS_PLOT <- ggplot(plot_data, aes(x = X, y = Y)) +
  geom_encircle(aes(group = Stage, color=Stage), s_shape = 2, expand = 0.25,size=3, linetype="dashed") +
  geom_jitter(size = 8.5,aes(fill=Tissue),alpha=0.7, shape=21, stroke=NA) +
  theme_classic() +
  coord_cartesian(xlim = c(min(plot_data$X)-0.5,max(plot_data$X)+1.3),
                  ylim = c(min(plot_data$Y)-0.1,max(plot_data$Y)+0.18)) +
  scale_fill_manual(values=c("#F8BF28","#38629A")) +
  scale_color_manual(values=c("#EA4335","#34A853","#f5a002")) +
  theme(legend.title = element_blank(), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 16,color="black"),
        legend.position = "None",
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = NA)) 


MDS_PLOT

ggsave("Figure2D.pdf",
       plot = MDS_PLOT,
       height = 5,
       width = 5)

#### Figure 5C ####


### Load chicken ###
load("./Data/Figure5_data/Figure5C_chicken.RData")

gene_list <- c("DMRT1")


gene_expression_chicken <- as.data.frame(logCPM[which(Annotation$`Gene name` %in% gene_list),])


colnames(gene_expression_chicken) <- "Expression"
gene_expression_chicken$Stage <-  meta_data$Stage

gene_expression_chicken$Sex <- meta_data$sex
gene_expression_chicken$tissue <- meta_data$Tissue

gene_expression_chicken$species <- "Chicken"

# Calculate standard deviation
gene_expression_chicken <- gene_expression_chicken %>%
  group_by(Sex, Stage) %>%
  mutate(
    StdDev = sd(Expression),   
    Mean = mean(Expression)    
  ) %>%
  ungroup() 



### Load quail ###

load("./Data/Figure5_data/Figure5C_quail.RData")
logCPM <- cpm(Y, normalized.lib.sizes = TRUE, log = T)

gene_expression_quail <- as.data.frame(logCPM[which(Annotation$`Gene name` %in% gene_list),])

colnames(gene_expression_quail) <- "Expression"
gene_expression_quail$Stage <-  meta_data$Stage

gene_expression_quail$Sex <- meta_data$sex
gene_expression_quail$tissue <- meta_data$Tissue

gene_expression_quail$species <- "quail"


gene_expression_quail <- gene_expression_quail %>%
  group_by(Sex, Stage) %>%
  mutate(
    StdDev = sd(Expression),   
    Mean = mean(Expression)    
  ) %>%
  ungroup()  


df <- rbind(gene_expression_chicken,gene_expression_quail)

df$label <- interaction(df$species,df$tissue)

df$Stage <- factor(df$Stage,levels = c("E4.5","E5", "E5.5", "E6", "E6.5", "E8.5", "E10", "E10.5", "E18.5", "E19"))

plot <- ggplot(df, aes(x=Stage, y=Mean)) +
  geom_line(aes(group=label, color=label),size=1.5,alpha=0.8) +
  geom_errorbar(aes(ymin=Mean-StdDev,ymax=Mean+StdDev),size=1,width=0.2, alpha=0.7) +
  geom_point(aes(group=label,color=label),size=3, shape=15) +
  labs(x="Stage", y="TMM normalized logCPM of\n Z chromosome genes") +
  scale_colour_manual(values = c("#d81159","#9A161F","#00239C","#0d3b66")) +
  #0d3b66
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12), 
        axis.title.x = element_blank(),
        legend.position = "None",
        panel.grid.major.y = element_line(color = "#d4d0cf", linetype = "dashed"),
        panel.grid.minor.y = element_line(color = "#d4d0cf", linetype = "dashed"))


plot

ggsave("DMRT1.pdf",width=5,height=3)


#### Figure 5D ####


data <- fread("./Data/Figure5_data/Figure5D_data.txt")


geneName <- data$geneName

data <- data[,-75]
gender_info <- ifelse(grepl("Female", colnames(data)), "Female", "Male")

stage_info <- str_extract(colnames(data), "E[0-9]+\\.?[0-9]*")


col_df <- data.frame(Stage = stage_info,
                     Gender = gender_info,
                     OriginalCol = colnames(data))


data_long <- pivot_longer(data, cols = everything(), names_to = "OriginalCol", values_to = "Expression")
data_long$geneName <- rep(geneName, each = ncol(data))

data_long <- merge(data_long, col_df, by = "OriginalCol",sort=F)


data_summary <- data_long %>%
  group_by(geneName, Stage, Gender) %>%
  summarise(AvgExpression = mean(Expression, na.rm = TRUE)) %>%
  pivot_wider(names_from = c("Stage", "Gender"), values_from = AvgExpression)


stages_order <- c("E4.5", "E5.5", "E6", "E6.5", "E8.5", "E10", "E10.5", "E18.5", "E19")


gender_order <- c("Female", "Male")


final_col_names <- c("geneName")  
for (stage in stages_order) {
  for (gender in gender_order) {
    final_col_names <- c(final_col_names, paste(stage, gender, sep = "_"))
  }
}


final_col_names <- final_col_names[final_col_names %in% names(data_summary)]
data_summary_ordered <- data_summary[, final_col_names, drop = FALSE]

geneName <- data_summary$geneName
gene_expression <- data_summary_ordered[,-1]

rownames(gene_expression) <- geneName

data <- gene_expression[, c(which(grepl("_Female", names(gene_expression))), which(grepl("_Male", names(gene_expression))))]
rownames(data) <- rownames(gene_expression)
data <- t(apply(data, 1, function(x) (x - mean(x)) / sd(x)))

# Make annotation column 
anno_col <- data.frame(Sex=factor(gsub(".*(Male|Female)", "\\1", colnames(data))))
rownames(anno_col) <- colnames(data)

pheatmap::pheatmap(
  data,
  border_color = "black",
  angle_col = 45,
  color= rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(10)),
  fontsize = 10,
  cluster_cols = F,
  annotation_col = anno_col,
  cellwidth = 10,
  cellheight=15,
  annotation_colors = list(`Sex` = c(Female="#d81159",Male="#0d3b66"))
)









#### Figure 5E ####

data <- fread("./Data/Figure5_data/Figure5E_data.txt")

geneName <- data$geneName

data <- data[,-7]
gender_info <- ifelse(grepl("Female", colnames(data)), "Female", "Male")

stage_info <- str_extract(colnames(data), "E[0-9]+\\.?[0-9]*")

col_df <- data.frame(Stage = stage_info,
                     Gender = gender_info,
                     OriginalCol = colnames(data))



data_long <- pivot_longer(data, cols = everything(), names_to = "OriginalCol", values_to = "Expression")
data_long$geneName <- rep(geneName, each = ncol(data))

data_long <- merge(data_long, col_df, by = "OriginalCol",sort=F)

data_summary <- data_long %>%
  group_by(geneName, Stage, Gender) %>%
  summarise(AvgExpression = mean(Expression, na.rm = TRUE)) %>%
  pivot_wider(names_from = c("Stage", "Gender"), values_from = AvgExpression)


stages_order <- c("E5","E6","E10")

gender_order <- c("Female", "Male")


final_col_names <- c("geneName")  
for (stage in stages_order) {
  for (gender in gender_order) {
    final_col_names <- c(final_col_names, paste(stage, gender, sep = "_"))
  }
}


final_col_names <- final_col_names[final_col_names %in% names(data_summary)]
data_summary_ordered <- data_summary[, final_col_names, drop = FALSE]


geneName <- data_summary$geneName
gene_expression <- data_summary_ordered[,-1]

rownames(gene_expression) <- geneName


data <- gene_expression[, c(which(grepl("_Female", names(gene_expression))), which(grepl("_Male", names(gene_expression))))]
rownames(data) <- rownames(gene_expression)
data <- t(apply(data, 1, function(x) (x - mean(x)) / sd(x)))


# Make annotation column 
anno_col <- data.frame(Sex=factor(gsub(".*(Male|Female)", "\\1", colnames(data))))
rownames(anno_col) <- colnames(data)


pheatmap::pheatmap(
  data,
  border_color = "black",
  angle_col = 45,
  color= rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(50)),
  fontsize = 13,
  cluster_cols = F,
  annotation_col = anno_col,
  annotation_colors = list(`Sex` = c(Female="#d81159",Male="#0d3b66"))
)









#### Figure 5F ####

options(scipen=30)
color_label <- function(data){
  label <- c()
  
  color <- c()
  #ABF200
  for(i in 1:nrow(data)){
    if(abs(data$logFC[i]) == 0){
      label <- c(label,"logFC=0")
    }else if(abs(data$logFC[i]) > 0  && abs(data$logFC[i]) <= 0.5){
      label <- c(label,"0<logFC<0.5")
    } else if(abs(data$logFC[i]) > 0.5  && abs(data$logFC[i]) <= 1){
      label <- c(label,"0.5<logFC<1")
    } else if(abs(data$logFC[i]) > 1 && abs(data$logFC[i]) <= 2 ){
      label <- c(label,"1<logFC<2")
    }else if(abs(data$logFC[i]) > 2 && abs(data$logFC[i]) <= 3){
      label <- c(label,"2<logFC<3")
    }else if(abs(data$logFC[i]) > 3 && abs(data$logFC[i]) <= 4){
      label <- c(label,"3<logFC<4")
    }else if(abs(data$logFC[i]) > 4 && abs(data$logFC[i]) <= 5){
      label <- c(label,"4<logFC<5")
    }else if(abs(data$logFC[i]) > 5 && abs(data$logFC[i]) <= 6){
      label <- c(label,"5<logFC<6")
    }else if(abs(data$logFC[i]) > 6){
      label <- c(label,"logFC>6")
    }
    
  }
  
  return(label)
}

data <- fread("./Data/Figure5_data/Figure5F_data.txt")



Chromosome_index <- which(data$`Chromosome/scaffold name`=="Z")


data <- data[Chromosome_index]


label_list <- color_label(data)


data$label <- factor(label_list,levels = c("0<logFC<0.5","0.5<logFC<1","1<logFC<2","2<logFC<3","3<logFC<4","4<logFC<5","5<logFC<6"))


data <- data %>% arrange(desc(label))

data$`Gene name`[which(data$`Gene name` == "")] <- data$geneSymbol[which(data$`Gene name` == "")]

num_colors <- length(unique(data$label))
colors <- rev(brewer.pal(num_colors, "RdYlGn"))


Dosage_brain <- ggplot(data=data,aes(x=`Gene start (bp)`,y=logFC,label=`Gene name`,color=label))+
  geom_jitter(alpha=0.8,size=5)+
  scale_color_manual(values=colors)+
  geom_smooth(method = "loess",aes(group=1),color="black",size=2)+
  geom_hline(yintercept = mean(data$logFC),size=1,alpha=0.6,linetype='longdash',color="red")+
  geom_text_repel(fontface="italic",color="black")+
  theme_bw()+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=12,color="black"),
    strip.text=element_text(size=12,color="black"),
    legend.position = "None"
  )

Dosage_brain

ggsave("Fig6F.pdf",width=10,height=4)

#### Figure 5G ####

options(scipen=30)
color_label <- function(data){
  label <- c()
  
  color <- c()
  #ABF200
  for(i in 1:nrow(data)){
    if(abs(data$logFC[i]) == 0){
      label <- c(label,"logFC=0")
    }else if(abs(data$logFC[i]) > 0  && abs(data$logFC[i]) <= 0.5){
      label <- c(label,"0<logFC<0.5")
    } else if(abs(data$logFC[i]) > 0.5  && abs(data$logFC[i]) <= 1){
      label <- c(label,"0.5<logFC<1")
    } else if(abs(data$logFC[i]) > 1 && abs(data$logFC[i]) <= 2 ){
      label <- c(label,"1<logFC<2")
    }else if(abs(data$logFC[i]) > 2 && abs(data$logFC[i]) <= 3){
      label <- c(label,"2<logFC<3")
    }else if(abs(data$logFC[i]) > 3 && abs(data$logFC[i]) <= 4){
      label <- c(label,"3<logFC<4")
    }else if(abs(data$logFC[i]) > 4 && abs(data$logFC[i]) <= 5){
      label <- c(label,"4<logFC<5")
    }else if(abs(data$logFC[i]) > 5 && abs(data$logFC[i]) <= 6){
      label <- c(label,"5<logFC<6")
    }else if(abs(data$logFC[i]) > 6){
      label <- c(label,"logFC>6")
    }
    
  }
  
  return(label)
}

data <- fread("./Data/Figure5_data/Figure5G_data.txt")


Chromosome_index <- which(data$`Chromosome/scaffold name`=="Z")


data <- data[Chromosome_index]

label_list <- color_label(data)

unique(label_list)

data$label <- factor(label_list,levels = c("0<logFC<0.5","0.5<logFC<1","1<logFC<2","2<logFC<3","3<logFC<4","4<logFC<5","5<logFC<6","logFC>6"))


data <- data %>% arrange(desc(label))

# Replace geneName
data$`Gene name`[which(data$`Gene name` == "")] <- data$geneSymbol[which(data$`Gene name` == "")]

# Palette
num_colors <- length(unique(data$label))
colors <- rev(brewer.pal(num_colors, "RdYlGn"))

colors
unique(data$label)



Dosage_brain <- ggplot(data=data,aes(x=`Gene start (bp)`,y=logFC,label=`Gene name`,color=label))+
  geom_jitter(alpha=0.8,size=5)+
  scale_color_manual(values=colors)+
  geom_smooth(method = "loess",aes(group=1),color="black",size=2)+
  geom_hline(yintercept = mean(data$logFC),size=1,alpha=0.6,linetype='longdash',color="red")+
  geom_text_repel(fontface="italic",color="black")+
  theme_bw()+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=12,color="black"),
    strip.text=element_text(size=12,color="black"),
    legend.position = "None"
  )

Dosage_brain

ggsave("Figure5G.pdf",width=10,height=4)




#### Figure 5H ####
stage <- c("E4.5","E5.5-E6","E6.5", "E8.5","E10-E10.5", "E18.5-E19")

filelist <- c("galgal6_4_5.txt","galgal6_5_5_6.txt","galgal6_6_5.txt","galgal6_8_5.txt","galgal6_10_10_5.txt","galgal6_18_5_E19.txt")


filelist <- paste0("./Data/Figure5_data/Sex_biased_genes/",filelist)

logFC <- c()
Zchromosome <- c()


for(i in 1:length(filelist)){
  data <- fread(filelist[i])
  data <- data[data$`Chromosome/scaffold name`=="Z",]
  Zchromosome <- c(Zchromosome,nrow(data))
  logFC <- c(logFC,sum(data$logFC>0 & data$FDR<=0.05))
}



df <- data.frame(Stage=stage,
                 All=Zchromosome,
                 NUM=logFC)



df$percentage <- (df$NUM/df$All)*100

df$Stage <- factor(df$Stage,levels=c("E4.5","E5.5-E6","E6.5", "E8.5","E10-E10.5", "E18.5-E19"))


ggplot(df, aes(x = Stage)) +
  geom_bar(aes(y = All), stat = "identity", fill = "grey",width = 1,color="black") +
  scale_y_continuous(expand=c(0,0))+
  geom_bar(aes(y = logFC, fill = Stage), stat = "identity",width = 1,color="black",size=1) +  
  scale_fill_viridis_d(begin = 0, end = 1, option = "D") + 
  theme_classic()+
  theme(legend.position = "None",
        axis.title=element_blank(),
        axis.text=element_text(size=10,color="black")
  )



ggsave("Barplot_galgal6.pdf",width=4,height=3)



#### Figure 5I ####
load("./Data/Figure5_data/Figure5I.RData")
### Extraction of Z chromosome genes ###


logCPM_Z <- logCPM[which(Annotation$`Chromosome/scaffold name`=="Z"),]

Annotation_Z <- Annotation[which(Annotation$`Chromosome/scaffold name`=="Z"),]



temp_df <- data.frame(MEAN=colMeans(logCPM_Z),
                      SE=apply(logCPM_Z, 2, function(x) sd(x) / sqrt(length(x))),
                      Sex=meta_data$sex,
                      Age=meta_data$Stage)

plot <- ggplot(temp_df, aes(x=Age, y=MEAN)) +
  geom_line(aes(group=Sex, color=Sex),size=1.5,alpha=0.8) +
  geom_errorbar(aes(ymin=MEAN-SE,ymax=MEAN+SE),size=1,width=0.1, alpha=0.7) +
  geom_point(aes(group=Sex,color=Sex),size=5, shape=15) +
  scale_colour_manual(values = c("#F5D042","#0A174E")) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        legend.position="None",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14,color="Black"),
        legend.text = element_text(size = 12), 
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = "#d4d0cf", linetype = "dashed"),
        panel.grid.minor.y = element_line(color = "#d4d0cf", linetype = "dashed"))


plot

ggsave("Quail_Zchromosome.pdf",
       plot = plot,
       height = 3,
       width = 3)


#### Figure 5J ####


### Load chicken ###
load("./Data/Figure5_data/Figure5C_chicken.RData")

gene_list <- c("ENSGALG00000051419", "ENSCJPG00005003198")


gene_expression_chicken <- as.data.frame(logCPM[which(Annotation$geneSymbol %in% gene_list),])


colnames(gene_expression_chicken) <- "Expression"
gene_expression_chicken$Stage <-  meta_data$Stage

gene_expression_chicken$Sex <- meta_data$sex
gene_expression_chicken$tissue <- meta_data$Tissue

gene_expression_chicken$species <- "Chicken"

gene_expression_chicken <- gene_expression_chicken %>%
  group_by(Sex, Stage) %>%
  mutate(
    StdDev = sd(Expression),   
    Mean = mean(Expression)   
  ) %>%
  ungroup()  



### Load quail ###

load("./Data/Figure5_data/Figure5C_quail.RData")
logCPM <- cpm(Y, normalized.lib.sizes = TRUE, log = T)

gene_expression_quail <- as.data.frame(logCPM[which(Annotation$geneSymbol %in% gene_list),])






colnames(gene_expression_quail) <- "Expression"
gene_expression_quail$Stage <-  meta_data$Stage

gene_expression_quail$Sex <- meta_data$sex
gene_expression_quail$tissue <- meta_data$Tissue

gene_expression_quail$species <- "quail"


gene_expression_quail <- gene_expression_quail %>%
  group_by(Sex, Stage) %>%
  mutate(
    StdDev = sd(Expression),   
    Mean = mean(Expression)    
  ) %>%
  ungroup()  


df <- rbind(gene_expression_chicken,gene_expression_quail)




df$label <- interaction(df$species,df$tissue)


df$Stage <- factor(df$Stage,levels = c("E4.5","E5", "E5.5", "E6", "E6.5", "E8.5", "E10", "E10.5", "E18.5", "E19"))

plot <- ggplot(df, aes(x=Stage, y=Mean)) +
  geom_line(aes(group=label, color=label),size=1.5,alpha=0.8) +
  geom_errorbar(aes(ymin=Mean-StdDev,ymax=Mean+StdDev),size=1,width=0.2, alpha=0.7) +
  geom_point(aes(group=label,color=label),size=3, shape=15) +
  labs(x="Stage", y="TMM normalized logCPM of\n Z chromosome genes") +
  scale_colour_manual(values = c("#d81159","#9A161F","#00239C","#0d3b66")) +
  #0d3b66
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12), 
        axis.title.x = element_blank(),
        legend.position = "None",
        panel.grid.major.y = element_line(color = "#d4d0cf", linetype = "dashed"),
        panel.grid.minor.y = element_line(color = "#d4d0cf", linetype = "dashed"))


plot

ggsave("lncRNA.pdf",width=5,height=3)


#### Figure 5K ####
data <- read.csv("./Data/Figure5_data/Figurek5K_L_data.csv",stringsAsFactors = TRUE)
colnames(data) <- gsub("X1.","1-",colnames(data))
data <- data[,-c(2,3,14)]

color <- c("#998dcb","#8da0cb","#66c2a5","#bfa37a","#ff91fd","#fadf8c","#a6d854","#e27e80")

df <- data %>% dplyr::filter(Species=="Coturnix japonica"|Species=="Mus musculus"|Species=="Homo sapiens" | Species=="Arabidopsis thaliana") 

colnames(df)


plt <- ggradar(df,base.size = 15,
               values.radar = c(0,0.5,1.0),
               axis.labels = colnames(df)[-1],
               axis.label.size = 7,
               group.colours = c("#180270","#F08B9A","#013B03","#590018"),
               gridline.min.colour="grey0",
               gridline.mid.colour="grey0",
               gridline.max.colour="grey0",
               grid.label.size = 10,
               group.line.width = 3,
               axis.line.colour = "gray0",
               background.circle.colour = "grey0",
               legend.position = "none",
               fill=TRUE,
               fill.alpha = 0.3) +
  theme(axis.text = element_text(hjust=0.5, vjust = 1),
        plot.margin = margin(0, 5, 0, 5, 'cm')) +
  coord_cartesian(clip = "off")


plt

ggsave(filename = "Figure6K.pdf",plot = plt,width = 17,height = 10)



#### Figure 5L ####

data <- read.csv("./Data/Figure5_data/Figurek5K_L_data.csv",stringsAsFactors = TRUE)
colnames(data) <- gsub("X1.","1-",colnames(data))
data <- data[,-c(2,3,14)]


color <- c("#998dcb","#8da0cb","#66c2a5","#bfa37a","#ff91fd","#fadf8c","#a6d854","#e27e80")

df <- data %>% dplyr::filter(Species=="Gallus gallus"|Species=="Mus musculus"|Species=="Homo sapiens" | Species=="Arabidopsis thaliana") 

colnames(df)


plt <- ggradar(df,base.size = 15,
               values.radar = c(0,0.5,1.0),
               axis.labels = colnames(df)[-1],
               axis.label.size = 7,
               group.colours = c("#180270","#F08B9A","#013B03","#590018"),
               gridline.min.colour="grey0",
               gridline.mid.colour="grey0",
               gridline.max.colour="grey0",
               grid.label.size = 10,
               group.line.width = 3,
               axis.line.colour = "gray0",
               background.circle.colour = "grey0",
               legend.position = "none",
               fill=TRUE,
               fill.alpha = 0.3) +
  theme(axis.text = element_text(hjust=0.5, vjust = 1),
        plot.margin = margin(0, 5, 0, 5, 'cm')) +
  coord_cartesian(clip = "off")


plt

ggsave(filename = "Figure6L.pdf",plot = plt,width = 17,height = 10)

