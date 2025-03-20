#CORRELOGRAM OF CIRCULES
#FROM: https://rpubs.com/bigcat/258548
install.packages("qgraph")
library(qgraph)

PCAf_Control <- read.delim("~/Parcerias/Haroldo/Haroldo/PCA/PCAf_Control.txt", header = TRUE, row.names = 1, sep = "\t")
PCAf_MDD = read.delim("~/Parcerias/Haroldo/Haroldo/PCA/PCAf_MDD.txt", header = TRUE, row.names = 1, sep = "\t")

#Control
PCAf_Control <- as.matrix(sapply(PCAf_Control, as.numeric))
cormat_c=cor(PCAf_Control, method="spearman") # or cormat=cor(dat[,1:12],dat[,1:12])

pdf("qgraph_control.pdf", width = 12, height = 12)
qgraph(cormat_c, shape = "circle",
       minimum = 0.4,
       posCol = "#58a391",
       negCol = "darkred",
       layout = "groups",
       vsize = 5,
       labels = colnames(cormat_c),
       label.cex = 2,
       edge.width = 0.4)  # Define uma espessura menor para as linhas
dev.off()

#MDD
PCAf_MDD <- as.matrix(sapply(PCAf_MDD, as.numeric))
cormat_MDD=cor(PCAf_MDD, method="spearman")

pdf("qgraph_MDD.pdf", width = 12, height = 12)
qgraph(cormat_MDD, shape = "circle",
       minimum = 0.4,
       posCol = "#58a391",
       negCol = "darkred",
       layout = "groups",
       vsize = 5,
       labels = colnames(cormat_MDD),
       label.cex = 2,
       edge.width = 0.4)  # Define uma espessura menor para as linhas
dev.off()

#Regressao
tabelaregress <- rbind(PCAf_Control, PCAf_MDD)

tabelaregress$group <- rownames(tabelaregress)

tabelaregress$group <- ifelse(grepl("Control", tabelaregress$group) == TRUE, 0, 1)

#modelo de regressao 
result <- glm(group ~ PAX6 + NEGR1 + PPP6C + SORCS3 , 
              data = tabelaregress, 
              family = binomial(link = "logit"))

data_values <- glm(group ~ PAX6 + NEGR1 + PPP6C + SORCS3 , 
                   data = tabelaregress, 
                   family = binomial(link = "logit"))

summary(data_values)

#plotar resultado
# Carregar os pacotes necessários
library(ggplot2)

genes = c("PAX6", "NEGR1", "PPP6C", "SORCS3")
x <- tabelaregress[colnames(tabelaregress) %in%genes]
x$group <- tabelaregress$group

result <- x %>% tidyr::gather("genes", "values", -group)

# Criando um data frame com os p-values correspondentes aos genes
p_values <- data.frame(
  genes = c("NEGR1", "PAX6", "PPP6C", "SORCS3"),
  pvalue = c(0.1354261, 0.0008477, 0.1919619, 0.0430961) # Valores extraídos do modelo
)

# Formatando os valores de p-value com 7 casas decimais (ou conforme necessário)
p_values$pvalue_text <- formatC(p_values$pvalue, format = "f", digits = 7)

# Modificando os rótulos dos genes no dataset result para incluir os valores de p-value
result$genes <- factor(result$genes, levels = p_values$genes,
                       labels = paste0(p_values$genes, "\npvalue = ", p_values$pvalue_text))

pdf("regression.pdf", width = 10, height = 3)
ggplot(result, aes(x = values, y = group)) +  
  geom_smooth(method = "glm", method.args = list(family = "binomial"),  
              color = "darkgreen", fill = "grey60", alpha = 0.2, size = 1.2, linetype = "solid") +  
  facet_wrap(~genes, nrow = 1, scales = "free_x") +  
  geom_point(aes(fill = factor(group, labels = c("no", "yes"))), shape = 21, size = 1, stroke = 1, color = "black") +  
  scale_fill_manual(name = "MDD", values = c("no" = "darkgreen", "yes" = "#693382")) +  
  labs(x = "Gene expression", y = "Probability") +  
  theme_minimal() +  
  theme(axis.title.x = element_text(size = 10, colour = "black", face = "bold"),  
        axis.title.y = element_text(size = 10, colour = "black", face = "bold"),
        strip.text =  element_text(size = 14, colour = "black"))
dev.off()