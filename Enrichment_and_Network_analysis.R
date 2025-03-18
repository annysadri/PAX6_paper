
# ----------       Fig 1C      -----------------
data_cathomas = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/Cathomas.txt", header = TRUE, sep = "\t")
data_oh = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/OH, et al.txt", header = TRUE, sep = "\t") #ok
data_ramaker = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/remarker_GSE80655.txt", header = TRUE, sep = "\t") #ok
data_trang = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/trang et al.txt", header = TRUE, sep = "\t")

#Nomeando os grupos
data_oh$Group <- "Oh, H. et al., 2022"
data_cathomas$Group <- "Cathomas, F. et al., 2022"
data_ramaker$Group <- "Ramaker, R. C. et al., 2017"
data_trang$Group <- "Trang T. L. et al., 2018"

#Vendo BP em comum e filtrando
# Lista dos dataframes
dataframes <- list(data_cathomas,
                   data_oh,
                   data_ramaker,
                   data_trang)

# Função para fazer o merge baseado na coluna "Term"
merge_data <- function(df1, df2) {
  merge(df1, df2, by = "Term", all = FALSE)
}
# Reduzir a lista de dataframes usando a função merge_data
all_enrich <- Reduce(merge_data, dataframes)
#Vetor de Terms em comum
commum_vector = as.character(all_enrich$Term)

#Filtrando
{
  #cathomas
  commum_cathomas <- data_cathomas %>% 
    filter(Term %in% commum_vector)
  #oh
  commum_oh <- data_oh %>% 
    filter(Term %in% commum_vector)
  #ramaker
  commum_ramaker <- data_ramaker %>% 
    filter(Term %in% commum_vector)
  #trang
  commum_trang <- data_trang %>% 
    filter(Term %in% commum_vector)
}

all_commum <- rbind(commum_cathomas,
                    commum_oh,
                    commum_ramaker,
                    commum_trang)

#Plotando
data_all <- all_commum %>%
  separate_rows(Genes, sep = ";")

data_all <- data_all %>%
  distinct(Term,
           Genes,
           Group, .keep_all = TRUE)
#Extraindo ID de GO
data_all <- data_all %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_all$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_all$Term)

#Completando vales faltantes
# data_all_sig <- data_all[data_all$Adjusted.P.value < 0.05, ]
data_all_sig = data_all
# all_combinations <- expand.grid(Genes = unique(data_all_sig$Genes),
#                                 Term = unique(data_all_sig$Term),
#                                 System = unique(data_all_sig$System))
# all_combinations <- expand.grid(Genes = unique(data_all_sig$Genes),
#                                 Term = unique(data_all_sig$Term))
all_combinations <- expand.grid(Genes = unique(data_all_sig$Genes),
                                Term = unique(data_all_sig$Term),
                                Group = unique(data_all_sig$Group))
# Combinar com os dados existentes
# merged_data <- merge(all_combinations,
#                      data_all_sig, by = c("Genes",
#                                           "Term",
#                                           "System"),
#                      all.x = TRUE)
# merged_data <- merge(all_combinations,
#                      data_all_sig, by = c("Genes",
#                                           "Term"),
#                      all.x = TRUE)
merged_data <- merge(all_combinations,
                     data_all_sig, by = c("Genes",
                                          "Term",
                                          "Group"),
                     all.x = TRUE)
# Ordenar os dados de acordo com o "Adjusted P-value"
merged_data <- merged_data %>% arrange(Odds.Ratio)

# Substituir NA na coluna Group
# merged_data$Group <- ifelse(is.na(merged_data$Group),
#                             "Oh, H. et al., 2022", merged_data$Group)
# Plot do heatmap
# Ordenando a variável group de acordo
merged_data <- merged_data %>%
  mutate(Group = factor(Group, levels = c("Oh, H. et al., 2022",
                                          "Ramaker, R. C. et al., 2017",
                                          "Cathomas, F. et al., 2022", #Boa
                                          "Trang T. L. et al., 2018")))

ordered_terms <- merged_data %>%
  dplyr::select(Term, Odds.Ratio, Group) %>%
  filter(Group == "Oh, H. et al., 2022") %>%
  arrange(Odds.Ratio)

merged_data$Term <- factor(merged_data$Term,
                           levels = unique(ordered_terms$Term))

# Criar o dotplot usando ggplot2
tiff("/Users/adrielnobile/Parcerias/Haroldo/Enrich/Datasets Separados/enrich_dotplot.tiff",
     width = 21,
     height = 8,
     res = 300, units = 'in')
ggplot(merged_data, aes(x = Term, y = Odds.Ratio)) +
  geom_point(aes(size = Odds.Ratio, color = Combined.Score)) +
  scale_size_continuous(range = c(2, 8), breaks = c(0.001, 0.01, 0.1),
                        name = "Combined Score") +  # Nome para a escala de tamanho
  labs(
    x = "Biological Process",
    y = "Odds Ratio",
    color = "Combined Score"  # Atualizado para "Combined Score"
  ) +
  scale_color_gradient(low = "#21918c", high = "#fde725") +  # Escala de cores contínuas
  theme_bw(base_size = 14) +  # Tamanho base da fonte
  facet_wrap(~Group, nrow = 1) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 18)
  ) +
  coord_flip()
dev.off()

# ----------       Fig 1D      -----------------
#Networking
# Combine os dataframes bp_oh, bp_cathomas, bp_ramaker e bp_trang
dfnet <- rbind(bp_oh, bp_cathomas, bp_ramaker, bp_trang)
# dfnet <- rbind(data_oh, data_cathomas, data_ramaker, data_trang)

# Remover duplicatas com base em Term, Genes, System e Group
net_data <- dfnet %>% distinct(Term, Genes, System, Group, .keep_all = TRUE)

# Separar as linhas das colunas de Genes
net_data <- net_data %>%
  separate_rows(Genes, sep = ";")

# Criar o objeto de rede
net <- network(net_data[, c("Term", "Genes")], directed = FALSE)

# Adicionar atributo System ao objeto network (CORES)
# Criar uma lista de atributos
vertex_color <- unique(net_data[, c("Genes", "System")])
vertex_color <- setNames(vertex_color$System, vertex_color$Genes)

# Adicionar atributos ao objeto network
set.vertex.attribute(net, "System", vertex_color[network.vertex.names(net)])

# Definir as cores manualmente
unique_systems <- unique(dfnet$System)  # Analise quais são os nomes do system
# Crie um data.frame com esses nomes definindo as cores para eles
df_color_dicio <- data.frame("System" = c("Immune", "Nervous", "BP"), 
                             "color" = c("#21918c", "#fde725", "gray"))

# Defina os nomes em um data.frame
color_palette <- setNames(df_color_dicio$color, df_color_dicio$System)  # O valor da cor

# Obtenha os System correspondentes aos nós
df_color <- data.frame("System" = get.vertex.attribute(net, "System"))
print(df_color)
#Substitua os valores NA por "Other"/ou qualquer outra coisa (nesse caso Processos Biológicos/BP)
df_color %>% 
  mutate(System = ifelse(is.na(System), "BP", System)) -> df_color
set.vertex.attribute(net, "System", df_color$System)

# Adicionar atributo Group ao objeto network (SHAPES)
# Criar uma lista de atributos
vertex_shape <- unique(net_data[, c("Genes", "Group")])
vertex_shape <- setNames(vertex_shape$Group, vertex_shape$Genes)

# Adicionar atributos ao objeto network
set.vertex.attribute(net, "Group", vertex_shape[network.vertex.names(net)])

# Definir as cores manualmente
unique(dfnet$Group)  # Analise quais são os nomes do system
# Crie um data.frame com esses nomes definindo as cores para eles
df_shape_dicio <- data.frame("Group" = c("Oh, H. et al., 2022", 
                                         "Cathomas, F. et al., 2022", 
                                         "Ramaker, R. C. et al., 2017", 
                                         "Trang T. L. et al., 2018",
                                         "BP"), 
                             "shape" = c(19, 18, 17, 15, 10))

# Defina os nomes em um data.frame
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$Group)  # O valor da cor

# Obtenha os Group correspondentes aos nós
df_shape <- data.frame("Group" = get.vertex.attribute(net, "Group"))
print(df_shape)
#Substitua os valores NA por "Other"/ou qualquer outra coisa (nesse caso Processos Biológicos/BP)
df_shape %>% 
  mutate(Group = ifelse(is.na(Group), "BP", Group)) -> df_shape
set.vertex.attribute(net, "Group", df_shape$Group)

# Plotar a rede usando ggnet2
tiff("/Users/adrielnobile/Parcerias/Haroldo/Haroldo_networking_labels.tiff",
     width = 16,
     height = 16,
     res = 300, units = 'in')
ggnet2(net, 
       size = "degree",
       # label = TRUE,
       # label.size = 4,
       repel = TRUE,
       alpha = 0.9,
       shape = "Group",
       shape.palette = shape_palette,
       size.legend.title = "Degree",  
       size.legend.position = "bottom",
       shape.mapping = aes(shape = degree),  # Use degree para shape diretamente
       shape.legend = TRUE,
       legend.size = 18,
       color = "System", # Adicione a cor baseada no grupo
       palette = color_palette
) +
  theme(legend.position = "right")
dev.off()

#SVG
ggsave("Haroldo_networking_labels.svg",  # Nome do arquivo e formato
       ggnet2(net, 
              size = "degree",
              label = TRUE,
              label.size = 4,
              repel = TRUE,
              alpha = 0.9,
              shape = "Group",
              shape.palette = shape_palette,
              size.legend.title = "Degree",  
              size.legend.position = "bottom",
              shape.mapping = aes(shape = degree),  # Use degree para shape diretamente
              shape.legend = TRUE,
              legend.size = 18,
              color = "System", # Adicione a cor baseada no grupo
              palette = color_palette
       ) +
         theme(legend.position = "right"),
       device = "svg",  # Define o formato do arquivo
       width = 12,  # Largura do gráfico em polegadas
       height = 9,  # Altura do gráfico em polegadas
       units = "in")  # Unidade de medida

# ----------       Fig 2A.     -----------------
#ScatterPlot
data.cluster = read.xlsx("~/Documents/DataSets/MetaVolcano - ArrayxBulk/3 - Metanalise HBCDV - PBMC & Liver/Enrichment/HBV_Liver/ScatterPlot/HBV_liver_pathways.xlsx")
data_cathomas = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/Cathomas.txt", header = TRUE, sep = "\t")
data_oh = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/OH, et al.txt", header = TRUE, sep = "\t") #ok
data_ramaker = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/remarker_GSE80655.txt", header = TRUE, sep = "\t") #ok
data_trang = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/trang et al.txt", header = TRUE, sep = "\t")

#Editando tabelas
#Criando tabela com identificacao de GO e deletando de "Terms"
data.cluster <- data.cluster %>%
  mutate(
    GO = str_extract(term, "(?<=\\().*(?=\\))"))
data.cluster$term <- sub("\\s*\\(.*?\\)\\s*", "", data.cluster$term)
data_cathomas <- data_cathomas %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_cathomas$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_cathomas$Term)
data_oh <- data_oh %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_oh$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_oh$Term)
data_trang <- data_trang %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_trang$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_trang$Term)
data_ramaker <- data_ramaker %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_ramaker$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_ramaker$Term)

#Unindo as tabelas por GO ID
data_cathomas_join <- left_join(data.cluster, data_cathomas,
                                by = "GO", suffix = c(".cluster", "_cathomas"),
                                keep = TRUE)
data_oh_join <- left_join(data.cluster, data_oh,
                          by = "GO", suffix = c(".cluster", "_oh"),
                          keep = TRUE)
data_trang_join <- left_join(data.cluster, data_trang,
                             by = "GO", suffix = c(".cluster", "_trang"),
                             keep = TRUE)
data_ramaker_join <- left_join(data.cluster, data_ramaker,
                               by = "GO", suffix = c(".cluster", "_ramaker"),
                               keep = TRUE)
View(data_cathomas)
#Transferindo outros Terms
data_cathomas_join <- data_cathomas_join %>%
  mutate(Term = if_else(is.na(Term), term, Term))
data_oh_join <- data_oh_join %>%
  mutate(Term = if_else(is.na(Term), term, Term))
data_trang_join <- data_trang_join %>%
  mutate(Term = if_else(is.na(Term), term, Term))
data_ramaker_join <- data_ramaker_join %>%
  mutate(Term = if_else(is.na(Term), term, Term))

#Inserindo P-valores nao significantes
data_cathomas_join$P.value[is.na(data_cathomas_join$P.value)] <- 1
data_oh_join$P.value[is.na(data_oh_join$P.value)] <- 1
data_trang_join$P.value[is.na(data_trang_join$P.value)] <- 1
data_ramaker_join$P.value[is.na(data_ramaker_join$P.value)] <- 1

#Ajeitando tabela pro Plot
data_1 = data_cathomas_join[,c(2:4,14,16)] 
data_2 = data_oh_join[,c(2:4,14,16)] 
data_3 = data_trang_join[,c(2:4,14,16)] 
data_4 = data_ramaker_join[,c(2:4,14,16)] 

#Agrupando coluna clusters por numero total
data_1 = data_1 |>
  group_by(cluster) |>
  mutate(total_number = n())
data_2 = data_2 |>
  group_by(cluster) |>
  mutate(total_number = n())
data_3 = data_3 |>
  group_by(cluster) |>
  mutate(total_number = n())
data_4 = data_4 |>
  group_by(cluster) |>
  mutate(total_number = n())
# Definindo "group" no data
data_1$group = "Cathomas, F. et al., 2022"
data_2$group = "Oh, H. et al., 2022"
data_3$group = "Trang T. L. et al., 2018"
data_4$group = "Ramaker, R. C. et al., 2017"

# Selecionando "pvalues" significantes
data_1$sig = ifelse(data_1$P.value < 0.05, 1, 0)
data_2$sig = ifelse(data_2$P.value < 0.05, 1, 0)
data_3$sig = ifelse(data_3$P.value < 0.05, 1, 0)
data_4$sig = ifelse(data_4$P.value < 0.05, 1, 0)

#Contabilizando Significantes
data_1 = data_1 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))
data_2 = data_2 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))
data_3 = data_3 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))
data_4 = data_4 |>
  group_by(cluster) |>
  mutate(filtered = sum(sig))
#Unindo dados
data = rbind(data_1,
             data_2,
             data_3,
             data_4)
#Separando "Cluster"
# Remover o prefixo "Cluster " e converter para fatores
data$cluster <- as.factor(gsub("Cluster ", "", data$cluster))
# clusters = paste("Cluster", 0:22, sep = " ")
clusters = paste(0:22)
data$cluster = factor(data$cluster, levels = clusters)
# Scatterplot
table(data$filtered)
mid = mean(data$filtered)

data_all = subset(data, filtered >= 5)

data_all = data_all |>
  group_by(cluster, group) |>
  summarise(x_mean = mean(x),
            y_mean = mean(y))

data$alpha = ifelse(data$filtered < 5, 0.1, 1)
#Separando Neuro e Imune
neuro_keywords <- c("Synaptic", "neuro", "Synapse",
                    "Perisynaptic", "Postsynaptic",
                    "Neurotransmitter", "Neuroinflammatory",
                    "Neuron", "Postsynapse", "Presynapse",
                    "Presynaptic", "Nervous", "Dopamine", "Axon",
                    "Neurogenesis", "Adrenergic", "Dendritic", "Neural",
                    "Neuromuscular", "Neuroblast", "Sensory", "Gliogenesis",
                    "Learning", "Brain")
immune_keywords <- c("leukocyte", "T cell", "mononuclear",
                     "lymphocyte", "immune response", "B cell",
                     "interferon", "Inflammatory", "Immune", "Defense",
                     "Cytokine", "Kinase", "ERBB2", "G Protein", 
                     "Cysteine", "ERK1 And ERK2", "Interleukin",
                     "Macrophage", "Necrosis", "Tumor", "MAPK",
                     "Healing", "Humoral", "Neutrophil", "Leukocyte",
                     "Myeloid", "Response", "Negative Regulation Of Type II Interferon-Mediated Signaling Pathway",
                     "Positive Regulation Of Growth", "MHC", "Cell Growth",
                     "Toll-Like Receptor", "Monocyte", "Interferon", "STAT",
                     "Prostaglandin", "NF-kappaB","GTPase", "Apoptotic", "Viral",
                     "Modulation", "Immunoglobulin", "NIK/NF-kappaB", "B Cell",
                     "T Cell", "I-kappaB kinase/NF-kappaB", "kinase", "Death",
                     "Ubiquitin", "Ubiquitin-Dependent","Ubiquitination",
                     "Proliferation", "Growth", "Mast", "Regulation")
data <- data %>%
  mutate(pathway = case_when(
    str_detect(Term, paste(neuro_keywords,
                           collapse = "|")) ~ "Nervous",
    str_detect(Term, paste(immune_keywords,
                           collapse = "|")) ~ "Immune",
    TRUE ~ "Others"
  ))
data <- data %>%
  mutate(Cores = case_when(
    pathway == "Immune" ~ "#21918c",
    pathway == "Nervous" ~ "#fde725",
    pathway == "Others" ~ "#440154",
    TRUE ~ NA_character_  # Se nenhum dos casos acima for atendido, definir como NA
  ))
# Ordenando a variável group de acordo com Brain e PBMC
data <- data %>%
  mutate(group = factor(group, levels = c("Oh, H. et al., 2022",
                                          "Ramaker, R. C. et al., 2017",
                                          "Cathomas, F. et al., 2022",
                                          "Trang T. L. et al., 2018")))
# Plotando
write.xlsx(data,
           "ScatterPlot_leiden.xlsx")
# Com Labels
tiff("/Users/adrielnobile/Parcerias/Haroldo/Enrich/Datasets Separados/ScatterStudyLabels.tiff",
     width = 10,
     height = 8,
     res = 300, units = 'in')
ggplot(data, aes(x, y)) +
  geom_point(colour = data$Cores, 
             shape = 16,
             aes(alpha = alpha)) +
  geom_label_repel(data = data_all,
                   mapping = aes(x = x_mean, y = y_mean,
                                 label = paste(cluster)),
                   colour = 'black',
                   size = 8,
                   # nudge_x = -1.5,  # Ajuste manual de deslocamento horizontal
                   # nudge_y = 1
  ) +  # Ajuste manual de deslocamento vertical +
  facet_wrap(~ group, dir = "h") +  
  theme_bw(base_size = 25) +  # Ajusta o tamanho das fontes
  theme(axis.text = element_text(size = 20),  # Ajusta o tamanho do texto dos eixos
        strip.text = element_text(size = 20))  # Ajusta o tamanho do texto do facet_wrap
dev.off()
# ----------       Fig 2B. -------------
#rrvgo Heatmap
#Carregando biblioteca
{
  library(rrvgo)
  library(stringr)
  library(shinydashboard)
  library(heatmaply)
  library(seriation)
  library(wordcloud)
  library(grid)
}

#Unindo
heat_data = rbind(data_cathomas,
                  data_oh,
                  data_ramaker,
                  data_trang)
#Atribuir os IDs para simMatrix
simMatrix <- calculateSimMatrix(heat_data$GO, #ele usa os IDs dos terms para fazer a análise
                                orgdb="org.Hs.eg.db",
                                ont="BP", #Selecione a Ontology
                                method="Rel")
scores <- setNames(-log10(heat_data$P.value), heat_data$GO) #Clusterizando
reducedTerms <- reduceSimMatrix(simMatrix, #Terms reduzidos
                                scores,
                                threshold=0.1,
                                orgdb="org.Hs.eg.db")
#Vetor pra Reduzir "ParentTerm"
parent_vector = c("regulation of short-term neuronal synaptic plasticity",
                  "neuron fate commitment", "positive regulation of neuron projection development",
                  "positive regulation of neuron projection development",
                  "regulation of nervous system process",
                  "regulation of neuron projection development",
                  "neuron projection extension", "positive regulation of neuron differentiation",
                  "regulation of AMPA receptor activity", "regulation of neurogenesis",
                  "neuropeptide signaling pathway", "positive regulation of neuron differentiation",
                  "neuron projection morphogenesis", "neuron differentiation",
                  "cellular response to interleukin-1", "positive regulation of interleukin-12 production",
                  "regulation of leukocyte degranulation", "regulation of myeloid cell differentiation",
                  "regulation of myeloid leukocyte mediated immunity", "inflammatory response",
                  "immune response-regulating cell surface receptor signaling pathway", "negative regulation of immune response",
                  "negative regulation of innate immune response", "central nervous system development",
                  "learning", "long-term synaptic depression", "positive regulation of synaptic plasticity",
                  "positive regulation of B cell proliferation", "positive regulation of MAPK cascade",
                  "response to cytokine", "lymphocyte proliferation", "calcineurin-NFAT signaling cascade",
                  "positive regulation of I-kappaB kinase/NF-kappaB signaling", "synaptic membrane adhesion", "positive regulation of synaptic transmission, glutamatergic")
heat_filtered = subset(reducedTerms,
                       parentTerm %in% parent_vector)
go_vector = as.data.frame(heat_filtered$go)
names(go_vector)[1] = "go"
#Filtrando matriz
# Extrair os elementos únicos do vetor
go_vector_unique <- unique(go_vector$go)
# Identificar os elementos que estão no data.frame e não no vetor
rows_to_remove <- setdiff(rownames(simMatrix), go_vector_unique)
cols_to_remove <- setdiff(colnames(simMatrix), go_vector_unique)
# Remover as linhas e colunas do data.frame
simMatrix_filtered <- simMatrix[!rownames(simMatrix) %in%
                                  rows_to_remove, !colnames(simMatrix) %in%
                                  cols_to_remove]
# Defina as cores da sua escala
cores_escala <- c("#440154","#21918c", "#fde725")
# Crie uma paleta de cores interpoladas entre as cores desejadas
escala_cores <- colorRampPalette(cores_escala)
# Defina o número de cores na sua escala (ajuste conforme necessário)
num_cores <- 100
heatmap_colors <- escala_cores(num_cores)
#Plotando e salvando
p <- heatmapPlot(simMatrix_filtered,
                 heat_filtered,
                 annotateParent=TRUE,
                 annotationLabel="parentTerm",
                 fontsize=16, col=heatmap_colors)
p
ggsave(filename = "/Users/adrielnobile/Parcerias/Haroldo/Enrich/Datasets Separados/Heatmap_study_pathway_size.tiff",
       plot = p,
       device = "tiff",  # Define o formato do arquivo
       width = 22,  # Largura do gráfico em polegadas
       height = 12,  # Altura do gráfico em polegadas
       units = "in",  # Unidade de medida
       dpi = 300)  # Resolução da imagem em dpi

#OutPut
write.xlsx(simMatrix_filtered,
           "SimMatrix_heatmap.xlsx", rowNames = TRUE)
#Input
write.xlsx(heat_data,
           "SimMatrix_heatmap_Input.xlsx")

# ----------       Fig 2C. -------------
# Enriquecimento separado
data_cathomas = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/Cathomas.txt", header = TRUE, sep = "\t")
data_oh = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/OH, et al.txt", header = TRUE, sep = "\t") #ok
data_ramaker = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/remarker_GSE80655.txt", header = TRUE, sep = "\t") #ok
data_trang = read.delim("~/Parcerias/Haroldo/Enrich/Datasets Separados/trang et al.txt", header = TRUE, sep = "\t")

#Criando vetores de filtragem
{
  #Cathomas
  vetor_cathoma_immune = c(
    "Negative Regulation Of Complement Activation (GO:0045916)",
    "Regulation Of Immune Effector Process (GO:0002697)",
    "Lymphocyte Mediated Immunity (GO:0002449)",
    "T Cell Mediated Immunity (GO:0002456)",
    "CD4-positive, Alpha-Beta T Cell Activation (GO:0035710)",
    "Regulation Of Toll-Like Receptor 7 Signaling Pathway (GO:0034155)",
    "Oxygen Transport (GO:0015671)",
    "Regulation Of T-helper 2 Cell Cytokine Production (GO:2000551)",
    "Alpha-Beta T Cell Activation (GO:0046631)",
    "Cellular Response To Cytokine Stimulus (GO:0071345)"
  )
  vetor_cathoma_nervous = c(
    "Regulation Of Excitatory Synapse Assembly (GO:1904889)",
    "Regulation Of Axon Guidance (GO:1902667)",
    "Positive Regulation Of T-helper 2 Cell Cytokine Production (GO:2000553)",
    "Regulation Of Presynapse Assembly (GO:1905606)",
    "Oligodendrocyte Differentiation (GO:0048709)",
    "Positive Regulation Of Nervous System Development (GO:0051962)",
    "Positive Regulation Of Immune Response (GO:0050778)",
    "Synapse Organization (GO:0050808)",
    "Nervous System Development (GO:0007399)",
    "Central Nervous System Development (GO:0007417)"
  )
  #oh
  vetor_oh_immune = c(
    "Regulation Of Dendritic Cell Chemotaxis (GO:2000508)",
    "Natural Killer Cell Activation (GO:0030101)",
    "Regulation Of Type II Interferon Production (GO:0032649)",
    "Positive Regulation Of Mononuclear Cell Migration (GO:0071677)",
    "Regulation Of T Cell Differentiation (GO:0045580)",
    "Positive Regulation Of Leukocyte Chemotaxis (GO:0002690)",
    "Inflammatory Response (GO:0006954)",
    "Regulation Of Inflammatory Response (GO:0050727)",
    "Interleukin-2-Mediated Signaling Pathway (GO:0038110)",
    "Dendritic Cell Differentiation (GO:0097028)",
    "T Cell Activation (GO:0042110)",
    "B Cell Proliferation (GO:0042100)"
  )
  vetor_oh_nervous = c(
    "Negative Regulation Of Glial Cell Apoptotic Process (GO:0034351)",
    "B Cell Chemotaxis (GO:0035754)",
    "Regulation Of Nervous System Development (GO:0051960)",
    "Regulation Of Neurogenesis (GO:0050767)",
    "Regulation Of Interleukin-8 Production (GO:0032677)",
    "B Cell Activation (GO:0042113)",
    "Generation Of Neurons (GO:0048699)",
    "Neuron Differentiation (GO:0030182)"
  )
  #ramaker
  vetor_ramaker_immune = c(
    "Response To Cytokine (GO:0034097)",
    "Regulation Of Inflammatory Response (GO:0050727)",
    "Response To Interleukin-1 (GO:0070555)",
    "Regulation Of I-kappaB kinase/NF-kappaB Signaling (GO:0043122)",
    "Response To Interferon-Beta (GO:0035456)",
    "Phagocytosis (GO:0006909)",
    "Leukocyte Cell-Cell Adhesion (GO:0007159)",
    "Interleukin-15-Mediated Signaling Pathway (GO:0035723)",
    "T Cell Proliferation (GO:0042098)",
    "Interleukin-4-Mediated Signaling Pathway (GO:0035771)"
  )
  vetor_ramaker_nervous = c(
    "Regulation Of Nervous System Process (GO:0031644)",
    "Peripheral Nervous System Neuron Development (GO:0048935)",
    "Excitatory Chemical Synaptic Transmission (GO:0098976)",
    "Immunological Synapse Formation (GO:0001771)",
    "Postsynapse Assembly (GO:0099068)",
    "Presynapse Organization (GO:0099172)",
    "Regulation Of Microglial Cell Activation (GO:1903978)",
    "Glial Cell Differentiation (GO:0010001)",
    "Synapse Assembly (GO:0007416)",
    "Nervous System Development (GO:0007399)"
  )
  
  #trang
  vetor_trang_immune = c(
    "Regulation Of T Cell Migration (GO:2000404)",
    "Regulation Of Cytokine Production Involved In Immune Response (GO:0002718)",
    "Regulation Of Type II Interferon Production (GO:0032649)",
    "Positive Regulation Of Cytokine Production (GO:0001819)",
    "Positive Regulation Of Interleukin-6 Production (GO:0032755)",
    "Positive Regulation Of Lymphocyte Migration (GO:2000403)",
    "Regulation Of Alpha-Beta T Cell Proliferation (GO:0046640)",
    "Regulation Of T-helper 1 Type Immune Response (GO:0002825)",
    "Regulation Of Interleukin-6 Production (GO:0032675)",
    "Positive Regulation Of T-helper 17 Type Immune Response (GO:2000318)"
  )
  vetor_trang_nervous = c(
    "Regulation Of Postsynaptic Neurotransmitter Receptor Activity (GO:0098962)",
    "Regulation Of Short-Term Neuronal Synaptic Plasticity (GO:0048172)",
    "Neurogenesis (GO:0022008)",
    "Regulation Of Presynapse Assembly (GO:1905606)",
    "Regulation Of Nervous System Process (GO:0031644)",
    "Immunological Synapse Formation (GO:0001771)",
    "Central Nervous System Development (GO:0007417)",
    "Nervous System Development (GO:0007399)",
    "Sensory Perception Of Smell (GO:0007608)",
    "Synapse Assembly (GO:0007416)"
  )
}
#Filtrando
{
  #cathomas
  immune_cathomas <- data_cathomas %>% 
    filter(Term %in% vetor_cathoma_immune)
  nervous_cathomas <- data_cathomas %>% 
    filter(Term %in% vetor_cathoma_nervous)
  #oh
  immune_oh <- data_oh %>% 
    filter(Term %in% vetor_oh_immune)
  nervous_oh <- data_oh %>% 
    filter(Term %in% vetor_oh_nervous)
  #ramaker
  immune_ramaker <- data_ramaker %>% 
    filter(Term %in% vetor_ramaker_immune)
  nervous_ramaker <- data_ramaker %>% 
    filter(Term %in% vetor_ramaker_nervous)
  #trang
  immune_trang <- data_trang %>% 
    filter(Term %in% vetor_trang_immune)
  nervous_trang <- data_trang %>% 
    filter(Term %in% vetor_trang_nervous)
}
#Iserindo Coluna de Nervous e Immune
{
  #Cathomas
  immune_cathomas$System <- "Immune"
  nervous_cathomas$System <- "Nervous"
  #Oh
  immune_oh$System <- "Immune"
  nervous_oh$System <- "Nervous"
  
  #Ramaker
  immune_ramaker$System <- "Immune"
  nervous_ramaker$System <- "Nervous"
  #Trang
  immune_trang$System <- "Immune"
  nervous_trang$System <- "Nervous"
}
#Unindo
{
  bp_oh <- rbind(immune_oh, nervous_oh)
  bp_cathomas <- rbind(immune_cathomas, nervous_cathomas)
  bp_ramaker <- rbind(immune_ramaker, nervous_ramaker)
  bp_trang <- rbind(immune_trang, nervous_trang)
}
#Definindo Grupos
bp_oh$Group <- "Oh, H. et al., 2022"
bp_cathomas$Group <- "Cathomas, F. et al., 2022"
bp_ramaker$Group <- "Ramaker, R. C. et al., 2017"
bp_trang$Group <- "Trang T. L. et al., 2018"

#Todos juntos
data_all <- rbind(bp_oh, bp_cathomas, bp_ramaker, bp_trang) #INPUT
#Salvando
write.xlsx(data_all,
           "BiologicalProcess_input.xlsx")

#Plotando separado
data_all <- bp_ramaker %>%
  separate_rows(Genes, sep = ";")

data_all <- data_all %>%
  distinct(Term,
           Genes,
           System, .keep_all = TRUE)
#Extraindo ID de GO
data_all <- data_all %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_all$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_all$Term)

#Completando vales faltantes
# data_all_sig <- data_all[data_all$Adjusted.P.value < 0.05, ]
data_all_sig = data_all
# all_combinations <- expand.grid(Genes = unique(data_all_sig$Genes),
#                                 Term = unique(data_all_sig$Term),
#                                 System = unique(data_all_sig$System))
all_combinations <- expand.grid(Genes = unique(data_all_sig$Genes),
                                Term = unique(data_all_sig$Term))

# Combinar com os dados existentes
# merged_data <- merge(all_combinations,
#                      data_all_sig, by = c("Genes",
#                                           "Term",
#                                           "System"),
#                      all.x = TRUE)
merged_data <- merge(all_combinations,
                     data_all_sig, by = c("Genes",
                                          "Term"),
                     all.x = TRUE)

# Ordenar os dados de acordo com o "Adjusted P-value"
merged_data <- merged_data %>% arrange(Odds.Ratio)

# Substituir NA na coluna Group
merged_data$Group <- ifelse(is.na(merged_data$Group),
                            "Ramaker, R. C. et al., 2017", merged_data$Group)
# Plot do heatmap
tiff("/Users/adrielnobile/Parcerias/Haroldo/Enrich/Datasets Separados/ram_2022_enrich_heatmap.tiff",
     width = 18, #18 #10
     height = 5, #5 #5
     res = 300, units = 'in')
ggplot(merged_data, aes(x = reorder(Genes, Odds.Ratio),
                        y = Term,
                        fill = Odds.Ratio)) +
  geom_tile(color = "black",
            na.rm = FALSE) +  # Adiciona contornos pretos, incluindo valores NA
  # coord_flip() +
  scale_fill_gradient(low = "#21918c",
                      high = "#fde725",
                      na.value = "white") +  # Define a cor para valores NA como branco
  theme_bw(14) +
  facet_wrap(~Group) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1,
                                   size = 9),  # Defina o tamanho da fonte para o eixo X
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16)) +  # Defina o tamanho da fonte para o eixo Y
  labs(x = "Genes", y = "Enrichments", fill = "Odds.Ratio")
dev.off()

# ----------       Fig 3E. -------------

#ggnet - Genes
data.blood = read.delim("~/Parcerias/Haroldo/enrich_blood_haroldo.txt", header = TRUE, sep = "\t")
data.blood2 = read.delim("~/Parcerias/Haroldo/data.blood_2.txt", header = TRUE, sep = "\t")

# data.brain = read.delim("~/Parcerias/Haroldo/enrich_brain_haroldo.txt", header = TRUE, sep = "\t")
data.brain = read.delim("~/Parcerias/Haroldo/GO_Biological_Process_2023_table.txt", header = TRUE, sep = "\t")
data_syn = read.delim("~/Parcerias/Haroldo/Enrich/SynGO_2024_table.txt", header = TRUE, sep = "\t")
#Nomeando
data.blood$group = c("PBMC")
data.blood2$group = c("PBMC")
data.brain$group = c("Cingulate")
data_syn$group = c("PBMC")

#Unindo
data_all = rbind(data.blood,
                 data.brain,
                 data_syn,
                 data.blood2)
net_data <- data_all %>%
  separate_rows(Genes, sep = ";")

# net_data <- net_data[net_data$P.value < 0.1, ]
net_data <- net_data[,c(1,9,10)]

#Criando data.frame de genes nao enriquecidos
# Exemplos de vetores de dados
Term <- c("Not mapped", 
          "Not mapped", 
          "Not mapped", 
          "Not mapped",
          "Not mapped")
Genes <- c("POGZ", "NRD1", "HIST1H3I", "RERE", "SAMD5")
group <- c("PBMC", "PBMC", "PBMC", "PBMC", "Cingulate")

# Criando o data.frame customizado
uncharacterized <- data.frame(Term, Genes, group, stringsAsFactors = FALSE)

net_total <- rbind(net_data,
                   uncharacterized)
net_data <- net_total
#Jeito antigo
dfnet_total <- net_data %>% distinct(Term,
                                     Genes,
                                     .keep_all = TRUE)

#Criando Vetores de separacao
gen_trang <- c(
  "SORCS3", "NEGR1", "SHISA9", "TCF4", "PAX6", "NRD1", 
  "MST1", "ACVR1B", "STK32A", "PPP6C", "MLEC", "USP19", 
  "HARS2", "FLRT1", "RBM4", "ARIH2", "POGZ"
)
gen_cat <- c("HIST1H3I", "KCNJ13")
gen_wit <- c("CTNNA3", "RERE", "SHANK2", "NUP43", "BTG3", "PPP3CC")
gen_oh <- c("OR2B2")
gen_ram <- c("CD40", "SERPING1", "SAMD5", "FES", "CHML")

#Separando
filt_trang <- dfnet_total[dfnet_total$Genes %in% gen_trang, ]
filt_cat <- dfnet_total[dfnet_total$Genes %in% gen_cat, ]
filt_wit <- dfnet_total[dfnet_total$Genes %in% gen_wit, ]
filt_oh <- dfnet_total[dfnet_total$Genes %in% gen_oh, ]
filt_ram <- dfnet_total[dfnet_total$Genes %in% gen_ram, ]

#Nomeando para shapes
unique(filt_trang$Genes)
filt_trang$Group <- "(B) Als x Trang (2018)"

unique(filt_cat$Genes)
filt_cat$Group <- "(B) Als x Cathomas (2022)"

unique(filt_wit$Genes)
filt_wit$Group <- "(C) Als x Witteenberg (2020)"

unique(filt_ram$Genes)
filt_ram$Group <- "(A) Als x Ramaker (2017)"

unique(filt_trang$Genes)
filt_oh$Group <- "(A) Als x Oh (2022)"

#reunindo
dfnet_total <- rbind(filt_cat,
                     filt_oh,
                     filt_ram,
                     filt_trang,
                     filt_wit)
#Salvando Input
write.xlsx(dfnet_total,
           "Network_bystudy.xlsx")

#Definindo genes Immune e nervoso
net <- network(dfnet_total[, c("Term", "Genes")], directed = FALSE)

# Adicionar atributo System ao objeto network (CORES)
# Criar uma lista de atributos
vertex_color <- unique(dfnet_total[, c("Genes", "Group")])
vertex_color <- setNames(vertex_color$Group, vertex_color$Genes)

# Adicionar atributos ao objeto network
set.vertex.attribute(net, "Group", vertex_color[network.vertex.names(net)])

# Definir as cores manualmente
unique_systems <- unique(dfnet_total$Group)  # Analise quais são os nomes do system
# Crie um data.frame com esses nomes definindo as cores para eles
df_color_dicio <- data.frame("Group" = c("(B) Als x Cathomas (2022)",
                                         "(A) Als x Oh (2022)",
                                         "(A) Als x Ramaker (2017)",
                                         "(B) Als x Trang (2018)",
                                         "(C) Als x Witteenberg (2020)",
                                         "BP"), 
                             "color" = c("#fde725",
                                         "#21918c",
                                         "#21918c",
                                         "#fde725",
                                         "#440154",
                                         "gray"))

# Defina os nomes em um data.frame
color_palette <- setNames(df_color_dicio$color, df_color_dicio$Group)  # O valor da cor

# Obtenha os Group correspondentes aos nós
df_color <- data.frame("Group" = get.vertex.attribute(net, "Group"))

#Substitua os valores NA por "Other"/ou qualquer outra coisa (nesse caso Processos Biológicos/BP)
df_color %>% 
  mutate(Group = ifelse(is.na(Group), "BP", Group)) -> df_color
set.vertex.attribute(net, "Group", df_color$Group)

# Adicionar atributo Group ao objeto network (SHAPES)
# Criar uma lista de atributos
vertex_shape <- unique(dfnet_total[, c("Genes", "Group")])
vertex_shape <- setNames(vertex_shape$Group, vertex_shape$Genes)

# Adicionar atributos ao objeto network
set.vertex.attribute(net, "Group", vertex_shape[network.vertex.names(net)])

# Definir os shapes manualmente
unique(dfnet_total$Group)  # Analise quais são os nomes do system
# Crie um data.frame com esses nomes definindo as cores para eles
df_shape_dicio <- data.frame("Group" = c("(B) Als x Cathomas (2022)",
                                         "(A) Als x Oh (2022)",
                                         "(A) Als x Ramaker (2017)",
                                         "(B) Als x Trang (2018)",
                                         "(C) Als x Witteenberg (2020)",
                                         "BP"), 
                             "shape" = c(17, 17, 19, 19, 19, 10))

# Defina os nomes em um data.frame
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$Group)  # O valor da cor

# Obtenha os Group correspondentes aos nós
df_shape <- data.frame("Group" = get.vertex.attribute(net, "Group"))
print(df_shape)
#Substitua os valores NA por "Other"/ou qualquer outra coisa (nesse caso Processos Biológicos/BP)
df_shape %>% 
  mutate(Group = ifelse(is.na(Group), "BP", Group)) -> df_shape
set.vertex.attribute(net, "Group", df_shape$Group)

# Visualize a rede com ggnet2
ggsave("NetworkHaroldo_Venn.svg",  # Nome do arquivo e formato
       ggnet2(net, 
              size = "degree",
              label = TRUE,
              label.size = 4,
              repel = TRUE,
              alpha = 0.9,
              shape = "Group",
              shape.palette = shape_palette,
              size.legend.title = "Degree",  
              size.legend.position = "bottom",
              shape.mapping = aes(shape = degree),  # Use degree para shape diretamente
              shape.legend = TRUE,
              legend.size = 18,
              color = "Group", # Adicione a cor baseada no grupo
              palette = color_palette
       ) +
         theme(legend.position = "right"),
       device = "svg",  # Define o formato do arquivo
       width = 9,  # Largura do gráfico em polegadas
       height = 6,  # Altura do gráfico em polegadas
       units = "in")  # Unidade de medida
ggsave("Haroldo_networking_labels.svg",  # Nome do arquivo e formato
       ggnet2(net, 
              size = "degree",
              label = TRUE,
              label.size = 5,
              repel = TRUE,
              alpha = 0.9,
              size.legend.title = "Degree",  
              size.legend.position = "bottom",
              shape.mapping = aes(shape = degree),  # Use degree para shape diretamente
              shape.legend = TRUE,
              legend.size = 18,
              color = "group", # Adicione a cor baseada no grupo
              palette = color_palette) +
         theme(legend.position = "right"),
       device = "svg",  # Define o formato do arquivo
       width = 8,  # Largura do gráfico em polegadas
       height = 8,  # Altura do gráfico em polegadas
       units = "in")  # Unidade de medida

# ----------       Fig 4C. -------------
#DandelionPlot
library(DandEFA)
library(writexl)
library(openxlsx)
library(psych)
library(viridis)
# Defina a primeira coluna como rownames
rownames(df) <- df[, 1]
# Remova a primeira coluna do dataframe
df <- df[, -1]
df <- na.omit(df)
#To Omit NAs
df_test = df
#Carregando Tabelas
metadado = read.xlsx("~/Parcerias/Haroldo/Dandelion/metadata.xlsx")
df_blood <- read.xlsx("~/Parcerias/Haroldo/Dandelion/count_blood_genes_comum.xlsx", rowNames = TRUE)
df_brain <- read.xlsx("~/Parcerias/Haroldo/Dandelion/count_brain_genes_comum.xlsx", rowNames = TRUE)
df_brain <- df_brain[, -1]
str(df_brain)
df_brain[, 1:48] <- lapply(df_brain[, 1:48], function(x) as.numeric(gsub("[^0-9.-]", "", x)))

#Para refazer
df = df_brain
df = df_blood
View(df)
#Transpor para numerico
df <- as.data.frame(t(df))

#combinando
df_new <- cbind(df,
                metadado)
df_new <- df_new[,c(1:17,22)]
df_new <- df_new[,c(1:17,22)]
# Filtrando apenas os elementos "Control" da coluna "Diagnostico"
df_control <- df_new %>%
  filter(Diagnostico == "Control")
df_MDD <- df_new %>%
  filter(Diagnostico == "MDD")
View(df_control)
View(df_MDD)
df <- df_control[,c(1:17)]
df <- df_MDD[,c(1:17)]
#Vetor de cor
colfunc <- colorRampPalette(c("#21918c","#fde725"))
dandpal <- colfunc(2)
# color names at: www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/
# Extrair fatores usando o número máximo de fatores
facl <- factload(df,nfac=4,method="mle",cormeth="pearson")

#To obtain the table of loading factor
load <- as.matrix.data.frame(facl)
load_df <- as.data.frame(load)
gene <- as.data.frame(row.names(facl))
load_fac <- cbind(gene, load_df)
# setwd("/Users/otavio.cmarques/Documents/Otavio/DISCIPLINA POS/R - Atualizado/Dandelion plot")
# setwd("/Users/adrielnobile/Parcerias/Haroldo/Dandelion")
# write_xlsx(load_fac, "load_fac.xlsx")
#Build the graphic: Define the directory where the abbreviation script is located
# setwd("/Users/otavio.cmarques/Documents/Otavio/DISCIPLINA POS/R - Atualizado/Dandelion plot")
# source("abbreviation.R")
# tiff("/Users/adrielnobile/Parcerias/Haroldo/dandelion_blood_infected.tiff",
#      width = 15, #15
#      height = 12, #12
#      res = 150, units = 'in') 
# dandelion(facl,bound=0,mcex=c(2.12,2.12),palet=dandpal)
# dev.off()

# for additional information see: https://rdrr.io/cran/DandEFA/man/factload.html
#FOr instance: mcex=c(1,1.2) it represent the size of letters in the factors and uniqueness/communalities respectively
#Make a heatmap of loading factors
#Access packages
library(BiocManager)
library(ComplexHeatmap)
library(circlize)

#PLOT RETANGULAR HEATMAP
col_fun1 = colorRamp2(c(-1, 0, 1), c("#440154","#21918c", "#fde725"))
num_colors <- 100
col_fun1 = viridis(num_colors)

row.names(load_df) <- row.names(facl)
load_df <- as.matrix(load_df)
View(load_df)
tiff("/Users/adrielnobile/Parcerias/Haroldo/Heatmap_dandelion_blood_control.tiff",
     width = 5,
     height = 8,
     res = 200, units = 'in') 
Heatmap(load_df,
        name="Factor Loading",
        row_title = "",
        row_title_side = c("right"),
        row_title_gp = gpar(fontsize = 22, fontface = "bold"),
        row_names_gp = gpar(fontsize= 16),
        row_title_rot = 270,
        row_dend_side = "left",
        row_dend_width = unit(1, "cm"),
        column_dend_side = c("top"),
        column_title = "Factors",
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        column_names_gp = gpar(fontsize= 12),
        column_dend_height = unit(1, "cm"), col = col_fun1)
dev.off()








# ----------       Fig 5C. ----------------
#Criando o Heatmap de Enriquecimento
#Carregando tabela
data_bp = read.delim("~/Parcerias/Haroldo/Enrich/GO_Biological_Process_2023_table.txt", header = TRUE, sep = "\t")
data_cc = read.delim("~/Parcerias/Haroldo/Enrich/GO_Cellular_Component_2023_table.txt", header = TRUE, sep = "\t")
data_mf = read.delim("~/Parcerias/Haroldo/Enrich/GO_Molecular_Function_2023_table.txt", header = TRUE, sep = "\t")
# data_reac = read.delim("~/Parcerias/Haroldo/Enrich/Reactome_2022_table.txt", header = TRUE, sep = "\t")
data_syn = read.delim("~/Parcerias/Haroldo/Enrich/SynGO_2024_table.txt", header = TRUE, sep = "\t")

#Definindo grupos
data_bp$enrich = c("BP")
data_cc$enrich = c("CC")
data_mf$enrich = c("MF")
# data_reac$enrich = c("Reactome")
data_syn$enrich = c("SynGO 2024")

#Unindo
data_all = rbind(data_bp,
                 data_mf,
                 data_syn)
data_all <- data_all %>%
  separate_rows(Genes, sep = ";")

data_all <- data_all %>%
  distinct(Genes, Term,
           enrich, .keep_all = TRUE)

#Extraindo ID de GO
data_all <- data_all %>%
  mutate(
    GO = str_extract(Term, "(?<=\\().*(?=\\))"))
data_all$Term <- sub("\\s*\\(.*?\\)\\s*", "", data_all$Term)

#Completando vales faltantes
data_all_sig <- data_all[data_all$Adjusted.P.value < 0.05, ]
all_combinations <- expand.grid(Genes = unique(data_all_sig$Genes),
                                Term = unique(data_all_sig$Term),
                                enrich = unique(data_all_sig$enrich))
# Combinar com os dados existentes
merged_data <- merge(all_combinations,
                     data_all_sig, by = c("Genes",
                                          "Term",
                                          "enrich"),
                     all.x = TRUE)

# Ordenar os dados de acordo com o "Adjusted P-value"
merged_data <- merged_data %>% arrange(Adjusted.P.value)
merged_data <- merged_data %>% filter(Term != "Negative Regulation Of Nervous System Development")
merged_data <- merged_data %>% filter(Term != "Negative Regulation Of Neurogenesis")

#Salvando input
write.xlsx(data_bp,
           "enrich_bp.xlsx")
write.xlsx(data_cc,
           "enrich_cc.xlsx")
write.xlsx(data_mf,
           "enrich_mf.xlsx")
write.xlsx(data_syn,
           "enrich_syngo.xlsx")

#Salvando output
write.xlsx(merged_data,
           "Heatmap_enrich_characterization.xlsx")

# Plot do heatmap
tiff("/Users/adrielnobile/Parcerias/Haroldo/Enrich/enrich_heatmap.tiff",
     width = 16,
     height = 10,
     res = 300, units = 'in') 
ggplot(merged_data, aes(x = reorder(Genes, Adjusted.P.value),
                        y = Term,
                        fill = Adjusted.P.value)) +
  geom_tile(color = "black",
            na.rm = FALSE) +  # Adiciona contornos pretos, incluindo valores NA
  # coord_flip() +
  scale_fill_gradient(low = "#21918c",
                      high = "#fde725",
                      na.value = "white") +  # Define a cor para valores NA como branco
  theme_bw(14) +
  facet_wrap(~enrich) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1,
                                   size = 12),  # Defina o tamanho da fonte para o eixo X
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16)) +  # Defina o tamanho da fonte para o eixo Y
  labs(x = "Genes",
       y = "Enrichments",
       fill = "Adjusted.P.value")
dev.off()
# ----------       Fig 7B.  --------------
#Enrichment by Brain Parts
setwd(dir = "~/Parcerias/Anny/Enrich/")

#Carregando dados
data_aINS = read.delim("aINS.txt", header = TRUE, sep = "\t")
data_amy = read.delim("enrich_amygdala.txt", header = TRUE, sep = "\t")
data_Cg25 = read.delim("Cg25_enriq.txt", header = TRUE, sep = "\t")
data_DPLFC = read.delim("DPLFC.txt", header = TRUE, sep = "\t")
data_Nac = read.delim("Nac.txt", header = TRUE, sep = "\t")
data_Orb = read.delim("orbitofrontal.txt", header = TRUE, sep = "\t")
data_Sub = read.delim("Sub.txt", header = TRUE, sep = "\t")

#Tranderindo GO ID e isolando Terms
data_aINS <- data_aINS %>%
  mutate(GO = str_extract(Term, "\\((.*?)\\)"), 
         Term = str_remove(Term, "\\s*\\(.*?\\)"))
data_amy <- data_amy %>%
  mutate(GO = str_extract(Term, "\\((.*?)\\)"), 
         Term = str_remove(Term, "\\s*\\(.*?\\)"))
data_Cg25 <- data_Cg25 %>%
  mutate(GO = str_extract(Term, "\\((.*?)\\)"), 
         Term = str_remove(Term, "\\s*\\(.*?\\)"))
data_DPLFC <- data_DPLFC %>%
  mutate(GO = str_extract(Term, "\\((.*?)\\)"), 
         Term = str_remove(Term, "\\s*\\(.*?\\)"))
data_Nac <- data_Nac %>%
  mutate(GO = str_extract(Term, "\\((.*?)\\)"), 
         Term = str_remove(Term, "\\s*\\(.*?\\)"))
data_Orb <- data_Orb %>%
  mutate(GO = str_extract(Term, "\\((.*?)\\)"), 
         Term = str_remove(Term, "\\s*\\(.*?\\)"))
data_Sub <- data_Sub %>%
  mutate(GO = str_extract(Term, "\\((.*?)\\)"), 
         Term = str_remove(Term, "\\s*\\(.*?\\)"))

#Vendo o que tem em comum
x = list(aINS = data_aINS$Term,
         Amygdala = data_amy$Term,
         Cg25 = data_Cg25$Term,
         DPLFC = data_DPLFC$Term,
         Nac = data_Nac$Term,
         Orb = data_Orb$Term,
         Sub = data_Sub$Term
)
x1 = make_comb_mat(x, mode = "intersect")

UpSet(x1)
Reduce(intersect, x)
#Definindo Grupos
data_aINS$Group <- "aINS"
data_amy$Group <- "Amygdala"
data_Cg25$Group <- "Cg25"
data_DPLFC$Group <- "DPLFC"
data_Nac$Group <- "Nac"
data_Orb$Group <- "Orb"
data_Sub$Group <- "Sub"

#Combinando tabelas
data_enrich <- rbind(data_aINS, data_amy, data_Cg25,
                     data_DPLFC, data_Nac, data_Orb, data_Sub)

neuro_keywords <- c("Synaptic", "neuro", "Synapse",
                    "Perisynaptic", "Postsynaptic",
                    "Neurotransmitter", "Neuroinflammatory",
                    "Neuron", "Postsynapse", "Presynapse",
                    "Presynaptic", "Nervous", "Dopamine", "Axon",
                    "Neurogenesis", "Adrenergic", "Dendritic", "Neural",
                    "Neuromuscular", "Neuroblast", "Sensory", "Gliogenesis",
                    "Learning", "Brain")
immune_keywords <- c("leukocyte", "T cell", "mononuclear",
                     "lymphocyte", "immune response", "B cell",
                     "interferon", "Inflammatory", "Immune", "Defense",
                     "Cytokine", "Kinase", "ERBB2", "G Protein", 
                     "Cysteine", "ERK1 And ERK2", "Interleukin",
                     "Macrophage", "Necrosis", "Tumor", "MAPK",
                     "Healing", "Humoral", "Neutrophil", "Leukocyte",
                     "Myeloid", "Response", "Negative Regulation Of Type II Interferon-Mediated Signaling Pathway",
                     "Positive Regulation Of Growth", "MHC", "Cell Growth",
                     "Toll-Like Receptor", "Monocyte", "Interferon", "STAT",
                     "Prostaglandin", "NF-kappaB","GTPase", "Apoptotic", "Viral",
                     "Modulation", "Immunoglobulin", "NIK/NF-kappaB", "B Cell",
                     "T Cell", "I-kappaB kinase/NF-kappaB", "kinase", "Death",
                     "Ubiquitin", "Ubiquitin-Dependent","Ubiquitination",
                     "Proliferation", "Growth", "Mast", "Regulation")

data_enrich <- data_enrich %>%
  mutate(System = case_when(
    str_detect(Term, paste(neuro_keywords,
                           collapse = "|")) ~ "Nervous",
    str_detect(Term, paste(immune_keywords,
                           collapse = "|")) ~ "Immune",
    TRUE ~ "Others"
  ))
#Removendo "Others"
data_enrich <- data_enrich %>% filter(System != "Others")

#Pegando top 5 Immune e Nervous
str(data_enrich)
#Plotando Network
# Remover duplicatas com base em Term, Genes, System e Group
names(data_enrich)[12] <- "System"
colnames(data_enrich)

dfnet <- data_enrich

genes <- c("CTNNA3", "FLRT1", "POGZ",
           "KCNJ13", "PAX6", "ACVR1B")

# Separar as linhas das colunas de Genes
dfnet <- dfnet %>%
  separate_rows(Genes, sep = ";")

net_data <- dfnet %>% distinct(Term, Genes,
                               .keep_all = TRUE)

write.xlsx(net_data,
           "Output_Network_brain_parts.xlsx")

# Criar o objeto de rede
net <- network(net_data[, c("Term", "Genes")], directed = FALSE)

# Criar uma lista de atributos
vertex_color <- unique(dfnet[, c("Genes", "Group")])
vertex_color <- setNames(vertex_color$Group, vertex_color$Genes)

# Adicionar atributos ao objeto network
# set.vertex.attribute(net, "Group", vertex_color[network.vertex.names(net)])

set.vertex.attribute(net, "Genes", vertex_color[network.vertex.names(net)])

# target <- c("Term", "Genes")
# Columns <- colnames(dfnet)[!colnames(dfnet)%in%target]

# dfnet %>% tidyr::gather("Term", "Column", -Columns) -> dfnet2
# Definir as cores manualmente
unique_systems <- unique(dfnet$Group)  # Analise quais são os nomes do system
# Crie um data.frame com esses nomes definindo as cores para eles
# df_color_dicio <- data.frame("Group" = c("aINS", "Amygdala", "Cg25",
#                                          "DPLFC", "Nac", "Orb", "Sub", "BP"), 
#                              "color" = c("#440154FF", "#fde725", "#3E4A89FF",
#                                          "#31688EFF", "#26828EFF", "#1F9D89FF",
#                                          "#35B779FF", "gray"))

df_color_dicio <- data.frame("Group" = c("aINS", "Amygdala", "Cg25",
                                         "DPLFC", "Nac", "Orb", "Sub", "BP"), 
                             "color" = c("#01204E", "#028391", "#F6DCAC",
                                         "#FEAE6F", "#799351", "#AF8F6F",
                                         "#A1DD70", "gray"))

# Defina os nomes em um data.frame
color_palette <- setNames(df_color_dicio$color, df_color_dicio$Group)  # O valor da cor

# Obtenha os Group correspondentes aos nós
df_color <- data.frame("Group" = get.vertex.attribute(net, "Group"))

#Substitua os valores NA por "Other"/ou qualquer outra coisa (nesse caso Processos Biológicos/BP)
df_color %>% 
  mutate(Group = ifelse(is.na(Group), "BP", Group)) -> df_color
set.vertex.attribute(net, "Group", df_color$Group)
print(df_color)
# Adicionar atributo System ao objeto network (SHAPES)
# Criar uma lista de atributos
vertex_shape <- unique(dfnet[, c("Genes", "System")])
vertex_shape <- setNames(vertex_shape$System, vertex_shape$Genes)

# Adicionar atributos ao objeto network
set.vertex.attribute(net, "System", vertex_shape[network.vertex.names(net)])

# Definir os shapes manualmente
unique(dfnet$System)  # Analise quais são os nomes do system
# Crie um data.frame com esses nomes definindo as cores para eles
df_shape_dicio <- data.frame("System" = c("Immune", "Nervous", "BP"), 
                             "shape" = c(19, 18, 17))

# Defina os nomes em um data.frame
shape_palette <- setNames(df_shape_dicio$shape, df_shape_dicio$System)  # O valor da cor

# Obtenha os System correspondentes aos nós
df_shape <- data.frame("System" = get.vertex.attribute(net, "System"))
print(df_shape)
#Substitua os valores NA por "Other"/ou qualquer outra coisa (nesse caso Processos Biológicos/Genes)
df_shape %>% 
  mutate(System = ifelse(is.na(System), "BP", System)) -> df_shape
set.vertex.attribute(net, "System", df_shape$System)

node_names = network.vertex.names(net)
genes_palette = sapply(node_names, function(gene) {
  if (gene %in% genes) {
    return('red')
  } else {
    return('black')
  }
})
bps <- network.vertex.names(net)[1:768]
genes_size = sapply(node_names, function(gene) {
  if (gene %in% bps) {
    return(0)
  } else {
    return(5)
  }
})

p <- ggnet2(net, 
            size = "degree",
            label = TRUE,
            label.size = genes_size,
            label.color = genes_palette,
            # mode = layout,
            repel = TRUE,
            alpha = 0.9,
            shape = "System",
            shape.palette = shape_palette,
            size.legend.title = "Degree",  
            size.legend.position = "bottom",
            edge.color = "gray90",
            shape.mapping = aes(shape = degree),  # Use degree para shape diretamente
            shape.legend = TRUE,
            legend.size = 18,
            color = "Group", # Adicione a cor baseada no grupo
            palette = color_palette
) +
  theme(legend.position = "right")
p
ggsave(filename = "NetworkBrainParts.svg",
       plot = p,
       device = "svg",  # Define o formato do arquivo
       width = 15,  # Largura do gráfico em polegadas
       height = 12,  # Altura do gráfico em polegadas
       units = "in")  # Unidade de medida
