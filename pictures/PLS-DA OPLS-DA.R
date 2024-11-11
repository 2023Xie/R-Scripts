if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")


library(readxl)
library(ropls)
library(ggplot2)
library(ggsci)
library(Cairo)
library(tidyverse)
library(extrafont)
loadfonts() #加载系统字体



#示例数据
# load data
data(sacurine)
names(sacurine)
# view data information
attach(sacurine)
strF(dataMatrix)
strF(sampleMetadata)
strF(variableMetadata)


# PCA analysis
pca = opls(dataMatrix)
genderFc = sampleMetadata[, "gender"]

pdf(file = 'figures/PCA.pdf', width = 5, height = 5)
plot(pca, typeVc = "x-score",
     parAsColFcVn = genderFc, parEllipsesL = TRUE)
dev.off()

#导入数据
serum <- read_xlsx("serum.xlsx", col_names = T)
serum_group <- read_xlsx("serum_group.xlsx", col_names = T)

# PLSDA analysis
plsda = opls(dataMatrix,genderFc)
plsda = opls(serum[,2:ncol(serum)],serum_group$groups)


# sample scores plot
sample.score = plsda@scoreMN %>% 
  as.data.frame() %>%
  mutate(group = serum_group$groups)

p1 = ggplot(sample.score, aes(p1, p2, color = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() +
  geom_point(aes(-4,-3), color = 'white') +
  labs(x = 'P1',y = 'P2') +
  stat_ellipse(level = 0.95, linetype = 'solid', 
               size = 1, show.legend = FALSE) +
  scale_color_manual(values = c("#1597A5","#FFC24B","#FEB3AE")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(color = 'black',size = 12, face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, face = 'plain'),
        axis.title = element_text(color = 'black',size = 15, face = 'plain'),
        axis.ticks = element_line(color = 'black'))
ggsave(p1, filename = 'pls.pdf', 
       width = 5, height = 5, device = cairo_pdf)

# VIP scores plot
vip.score = as.data.frame(plsda@vipVn)
colnames(vip.score) = 'vip'
vip.score$metabolites = rownames(vip.score)
vip.score = vip.score[order(-vip.score$vip),]
vip.score$metabolites = factor(vip.score$metabolites,
                               levels = vip.score$metabolites)

loading.score = plsda@loadingMN %>% as.data.frame()
loading.score$metabolites = rownames(loading.score)

all.score = merge(vip.score, loading.score, by = 'metabolites')

all.score$cat = paste('A',1:nrow(all.score), sep = '')

p2 = ggplot(all.score[all.score$vip >= 1,], aes(metabolites, vip)) +
  geom_segment(aes(x = metabolites, xend = metabolites,
                   y = 0, yend = vip)) +
  geom_point(shape = 21, size = 5, color = '#008000' ,fill = '#008000') +
  geom_point(aes(1,2.5), color = 'white') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'VIP value') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(color = 'black',size = 12, face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, face = 'plain'),
        axis.text.x = element_text(angle = 0),
        axis.title = element_text(color = 'black',size = 15, face = 'plain'),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank())
ggsave(p2, filename = 'pls_VIP.pdf', 
       width = 8, height = 5, device = cairo_pdf)




# OPLS-DA analysis 二分类
oplsda = opls(serum[,2:ncol(serum)],serum_group$groups, predI = 1, orthoI = NA)

# sample scores plot
sample.score = plsda@scoreMN %>% 
  as.data.frame() %>%
  mutate(group = serum_group$groups)

# sample scores plot
sample.score = oplsda@scoreMN %>% 
  as.data.frame() %>%
  mutate(group = serum_group$groups,
         o1 = oplsda@orthoScoreMN[,1])

p3 = ggplot(sample.score, aes(p1, o1, color = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() +
  #geom_point(aes(-10,-10), color = 'white') +
  labs(x = 'P1',y = 'to1') +
  stat_ellipse(level = 0.95, linetype = 'solid', 
               size = 1, show.legend = FALSE) +
  scale_color_manual(values = c('#008000','#FFA74F')) +
  theme_bw() +
  theme(legend.position = c(0.1,0.85),
        legend.title = element_blank(),
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'))
ggsave(p3, filename = 'figures/opls.pdf', 
       width = 5, height = 5, device = cairo_pdf)


# VIP scores plot
vip.score = as.data.frame(oplsda@vipVn)
colnames(vip.score) = 'vip'
vip.score$metabolites = rownames(vip.score)
vip.score = vip.score[order(-vip.score$vip),]
vip.score$metabolites = factor(vip.score$metabolites,
                               levels = vip.score$metabolites)

loading.score = oplsda@loadingMN %>% as.data.frame()
loading.score$metabolites = rownames(loading.score)

all.score = merge(vip.score, loading.score, by = 'metabolites')

all.score$cat = paste('A',1:nrow(all.score), sep = '')

p4 = ggplot(all.score[all.score$vip >= 1,], aes(cat, vip)) +
  geom_segment(aes(x = cat, xend = cat,
                   y = 0, yend = vip)) +
  geom_point(shape = 21, size = 5, color = '#008000' ,fill = '#008000') +
  geom_point(aes(1,2.5), color = 'white') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'VIP value') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank())
p4
ggsave(p4, filename = 'figures/opls_VIP.pdf', 
       width = 8, height = 5, device = cairo_pdf)


# model training
oplsda.2 = opls(dataMatrix, genderFc, predI = 1, orthoI = NA,subset = "odd") 


OPLS-DA
92 samples x 109 variables and 1 response
standard scaling of predictors and response(s)
R2X(cum) R2Y(cum) Q2(cum)
Total     0.26    0.825   0.608
RMSEE RMSEP pre ort
Total 0.213 0.341   1   2
Warning message:
  'permI' set to 0 because train/test partition is selected.

# 模型在训练集上的准确率
trainVi = getSubsetVi(oplsda.2)
tab = table(genderFc[trainVi], fitted(oplsda.2))
print(paste('模型准确率：',round(sum(diag(tab))/sum(tab)*100, 2),'%', sep = ''))

M  F
M 50  0
F  0 42
[1] "模型准确率：100%"

# model on test data
tab2 = table(genderFc[-trainVi],predict(oplsda.2, dataMatrix[-trainVi, ]))
print(paste('模型准确率：',round(sum(diag(tab2))/sum(tab2)*100, 2),'%', sep = ''))

M  F
M 43  7
F  7 34
[1] "模型准确率：84.62%"

# volcano plot
df = dataMatrix %>% as.data.frame()
df$gender = sacurine[["sampleMetadata"]][["gender"]]
df = df[order(df$gender),]
df = df[,-110]

M.mean = apply(df[1:100,],2,FUN = mean)
F.mean = apply(df[101:183,],2,FUN = mean)

FC = M.mean / F.mean
log2FC = log(FC,2)

pvalue = apply(df, 2, function(x)
{t.test(x[1:100],x[101:183])$p.value})

p.adj = p.adjust(pvalue, method = 'BH')
p.adj.log = -log10(p.adj)

colcano.df = data.frame(log2FC,p.adj, p.adj.log)
colcano.df$cat = ifelse(colcano.df$log2FC >= 1 & colcano.df$p.adj < 0.05,'Up',
                        ifelse(colcano.df$log2FC <= -1 & colcano.df$p.adj < 0.05,'Down','NS'))

p5 = ggplot(colcano.df, aes(log2FC, p.adj.log)) +
  geom_point() +
  labs(y = '-log10(p-value.adj)') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank())
p5
ggsave(p5, filename = '20201214PLSDA分析/figures/volcano.pdf', 
       width = 5, height = 5, device = cairo_pdf)
