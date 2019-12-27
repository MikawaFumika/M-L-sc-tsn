library(Rtsne)
library(ggplot2)

gse103322_df <- read.csv('GSE103322_high_expression_clean.csv', sep=',', head=T)

tsne.unc <- Rtsne(gse103322_df, distance=FALSE, perplexity=30)
unc.tsne <- data.frame(tsne.unc$Y)
colnames(unc.tsne) <- c("Dim1", "Dim2")

ggplot(unc.tsne, aes(x=Dim1, y=Dim2)) +
  geom_point(size=2, shape=21, color='red') + theme_classic() +
  scale_fill_manual(values=interact.cols) +
  scale_y_continuous(limits=c(-25, 25)) +
  scale_x_continuous(limits=c(-25, 25)) +
  labs(x="tSNE 1", y="tSNE 2") +
  #guides(fill=FALSE) +
  theme(axis.title=element_text(size=24), axis.text=element_text(size=16))
