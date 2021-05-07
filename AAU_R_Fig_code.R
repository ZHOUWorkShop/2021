## All codes are written on R-4.0.3
## Fig 1 A
setwd("F:\\AAU")
data <- read.csv("AUU_lianxu.csv",stringsAsFactors = F,row.names = 1)
data <- data[order(data$relapse_event),]
data <- data[,c(-1,-4,-22)]
heat_cli <- data[,1:5]
heat_exp <- data[,6:43]

annotation_col = data.frame( gender = heat_cli[,2],
                             #predict_label = test_heat_cli[,2],
                             relapse_event = as.factor(heat_cli[,1]),
                             Ankylosing_spondylitis = as.factor(heat_cli[,3]),
                             hypertension = as.factor(heat_cli[,5]),
                             diabetes = as.factor(heat_cli[,4])
                             
)
rownames(annotation_col) <- rownames(heat_cli)
ann_colors = list(
  #predict_label = c( "low"= "#00676B", "high" = "#BC3C28"),
  gender = c( "male"= "#D95F02", "female" = "#7570B3"),
  relapse_event = c("0" = "#E7298A","1" = "#66A61E"),
  Ankylosing_spondylitis = c("0" = "#F781BF","1" = "#A65628"),
  hypertension = c("0" = "#5175A4","1" = "#D99748"),
  diabetes = c("0" = "#D94E48","1" = "#800080")
) 
red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
heat_exp <- t(heat_exp)
linshi <- apply(heat_exp, 1, scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(heat_exp)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
pheatmap(linshi,fontsize=6,cutree_col = 4,cellheight = 5,cellwidth = 1 ,
         color  = colorRampPalette(c(blue,white,red))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F
)  

## Fig 1 B
setwd("F:\\AAU")
data <- read.csv("AUU_lianxu.csv",stringsAsFactors = F,row.names = 1)
cox <- matrix(,1,7)
colnames(cox) <- c('lower .95','upper .95','coef', 'exp(coef)','se(coef)','z','Pr(>|z|)' )
for (i in 3:46){
  a <- summary(coxph(Surv(data$relapse_time,data$relapse_event) ~ data[,i]))
  colnames(data)[i]
  cox <- rbind(cox,c(a$conf.int[,3:4],a$coefficients[,1:5]))
  #cox <- rbind(cox,cbind(a$conf.int[,3:4],a$coefficients[,1:5]))
}
cox <- cox[-1,]
rownames(cox) <- colnames(data)[3:46]
cox[,7] <- round(cox[,7],3)
cox <- cox[order(cox[,7]),]
cox <- cbind(cox,paste(round(cox[,4],3)," (",round(cox[,1],3)," - ",round(cox[,2],3),")",sep = ""))
#write.csv(cox,"cox_forest.csv",quote = F)
data <- read.csv("cox_forest.csv",stringsAsFactors = F)
forestplot(as.matrix(data[,c(1,4,9)]),data$exp.coef.,data$lower..95,
           data$upper..95,zero = 1,xlog = F,
           clip = c(0,2),
           colgap = unit(6,"mm"),graphwidth=unit(70,"mm"),
           lineheight = unit(0.8,"cm"),graph.pos = 3,
           col = fpColors(box="black", lines="black", zero = "black"),boxsize = 0.3,
           ci.vertices = T,ci.vertices.height = 0.2,
           lty.ci = 7,lwd.zero=0.4,lwd.ci = 3,
           txt_gp=fpTxtGp(label = gpar(cex=0.8),
                          ticks = gpar(cex=0.8)) ,
           is.summary=c(TRUE,rep(FALSE,100))
)

## Fig 2,3 A
score = train$Ankylosing.spondylitis*0.09230 + train$HLA.B27*0.19863 + 
  train$MO*(-0.59456) + train$HDL*0.36348 + train$LDL*(-0.12934) +  0.3287
train <- read.csv("train_clinical.csv",stringsAsFactors = F,row.names = 1)
test <- read.csv("test_clinical.csv",stringsAsFactors = F,row.names = 1)
train <- train[,c(-8)]
mod <- glm(relapse_event ~., data=train[,-1])
result_train <- predict(mod,train,type = "response")
aa <- roc(train$relapse_event,result_train)

result_test <- predict(mod,test,type = "response")
bb <- roc(test$relapse_event,result_test)

roc(train$relapse_event,result_train, plot=TRUE, print.thres=TRUE, ci=TRUE,
    print.auc=TRUE,legacy.axes = TRUE,col = "red")

roc(test$relapse_event,result_test, plot=TRUE, print.thres=TRUE, ci=TRUE,
    print.auc=TRUE,legacy.axes = TRUE,col = "red")

## Fig 2 B,C,E,F,D
train_label <- c()
result_train <- result_train[rownames(train)]
cutoff <- 0.328
for (i in 1:length(result_train)){
  if (result_train[i] > cutoff){
    train_label <- c(train_label,"high")
  }
  else {
    train_label <- c(train_label,"low")
  }
}
blue <- "#00676B"
red <- "#BC3C28"
train_surv <- survfit(Surv(train$relapse_time,train$relapse_event)~train_label,data = train)
train_surv
survdiff(Surv(train$relapse_time,train$relapse_event)~train_label,data = train)
summary(train_surv)
ggsurvplot(train_surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#
           legend.labs=c("high","low")
)
summary(coxph(Surv(train$relapse_time,train$relapse_event) ~ result_train))  ## HR
data2 <- data.frame(score = result_train, 
                    cluster = as.factor(train$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  ylab("Predict score") # 
###   heatmap 
data <- read.csv("AUU_连续.csv",stringsAsFactors = F,row.names = 1)
train_heat <- cbind(data[rownames(train),3],train_label,train[,-1])
train_heat <- train_heat[names(result_train[order(result_train)]),]
result_train <- result_train[rownames(train_heat)]
train_heat_cli <- train_heat[,1:5]
train_heat_exp <- cbind(train_heat[,6:8],result_train)
train_heat_exp[,1] <- scale(train_heat_exp[,1])
train_heat_exp[,2] <- scale(train_heat_exp[,2])
train_heat_exp[,3] <- scale(train_heat_exp[,3])
#train_heat_exp[,4] <- scale(train_heat_exp[,4])
aa <- train_heat_exp[rownames(train_heat_cli)[which(train_heat_cli$train_label == "high")],]
bb <- train_heat_exp[rownames(train_heat_cli)[which(train_heat_cli$train_label == "low")],]
aa <- aa[rev(order(aa$result_train)),]
bb <- bb[rev(order(bb$result_train)),]
train_heat_exp <- as.matrix(rbind(aa,bb))
annotation_col = data.frame( gender = train_heat_cli[,1],
                             predict_label = train_heat_cli[,2],
                             relapse_event = as.factor(train_heat_cli[,3]),
                             Ankylosing_spondylitis = as.factor(train_heat_cli[,4]),
                             HLA.B27 = as.factor(train_heat_cli[,5])
)
rownames(annotation_col) <- rownames(train_heat_cli)
ann_colors = list(
  predict_label = c( "low"= "#00676B", "high" = "#BC3C28"),
  gender = c( "male"= "#7570B3", "female" = "#D2A6C7"),
  relapse_event = c("0" = "#5F98C6","1" = "#FF9933"),
  Ankylosing_spondylitis = c("0" = "#FBB4AE","1" = "#BC805E"),
  HLA.B27 = c("0" = "#B1C7BA","1" = "#E6BD9F")
) 
red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
hist(train_heat_exp)
train_heat_exp[train_heat_exp>2] <- 2
train_heat_exp[train_heat_exp<(-2)] <- c(-2)
result_train <- result_train[rownames(train_heat_exp)]
pheatmap(t(train_heat_exp[,1:3]),fontsize=6,cutree_col = 4,cellheight = 25,cellwidth = 5 ,
         color  = colorRampPalette(c(blue,white,red))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T
)  #
pheatmap(t(train_heat_exp[,4]),fontsize=6,cutree_col = 4,cellheight = 25,cellwidth = 5 ,
         color  = colorRampPalette(c(white,"#0B0B0B"))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T
)  #
train_fufa <- train_heat[which(train_heat$relapse_event == 1),]
train_bufufa <- train_heat[which(train_heat$relapse_event == 0),]
table(train_fufa$HLA.B27)
a <- matrix(c(13,2,28,13),2,2)
fisher.test(a)
### 相关性分析
train_heat_exp <- train_heat[,6:8]
result_train <- result_train[rownames(train_heat_exp)]
data22 <- as.data.frame(cbind(result_train,train_heat_exp))
a1 <- ggplot(data = data22 , aes(x = MO, y = result_train)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + 
  ylab("Score") + theme_bw()+ #背景变为白色 
  stat_cor(method = "pearson", label.x = 0, label.y = 0.6) +
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #
        axis.text.y=element_text(size=16,colour="black"), #
        axis.title.y=element_text(size = 20,colour="black"), #
        legend.text=element_text( colour="black",  #
                                  size=16),
        legend.title=element_text( colour="black", #
                                   size=18),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  theme_bw()
a2 <- ggplot(data = data22 , aes(x = HDL, y = result_train)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + 
  ylab("Score") + theme_bw()+ # 
  stat_cor(method = "pearson", label.x = 0, label.y = 0.6) +
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #
        axis.text.y=element_text(size=16,colour="black"), #
        axis.title.y=element_text(size = 20,colour="black"), 
        legend.text=element_text( colour="black",  #
                                  size=16),
        legend.title=element_text( colour="black", #
                                   size=18),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  theme_bw()
a3 <- ggplot(data = data22 , aes(x = LDL, y = result_train)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + ## 
  ylab("Score") + theme_bw()+ # 
  stat_cor(method = "pearson", label.x = 1, label.y = 0.6) + ## 
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #
        axis.text.y=element_text(size=16,colour="black"), 
        axis.title.y=element_text(size = 20,colour="black"), #
        legend.text=element_text( colour="black",  #
                                  size=16),
        legend.title=element_text( colour="black", #
                                   size=18),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  theme_bw()
ggarrange(a1,a2,a3,nrow = 2,ncol = 2)
## 
setwd("F:\\AAU")
train <- read.csv("train_clinical.csv",stringsAsFactors = F,row.names = 1)
a1 <- ggplot(data = train, aes(x = train_label , y = MO,fill=train_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #
  theme_bw()+ #
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())  #
a2 <- ggplot(data = train, aes(x = train_label , y = HDL,fill=train_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())  #
a3 <- ggplot(data = train, aes(x = train_label , y = LDL,fill=train_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())  #
a= "#00AE72"
b= "#5175A4"
c= "#D94E48"
d= "#8C63A4"
data <- data.frame(type = train$train_label,
                   Ank = train$Ankylosing.spondylitis,
                   clu1 = paste(train$train_label,train$Ankylosing.spondylitis,sep = ""),
                   clu2 = paste(train$train_label,train$HLA.B27,sep = "") )
data$number <- 1
data <- ddply(data,'type',transform,percent = 1/sum(number)*100)
a4 <- ggplot(data,aes(type,percent,fill=clu1))+
  scale_fill_manual(values = c(a,b,c,d))+ #
  geom_bar(stat="identity",position="stack")+
  theme_bw()+ ylab("Ankylosing.spondylitis") +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())
a5 <- ggplot(data,aes(type,percent,fill=clu2))+
  scale_fill_manual(values = c(a,b,c,d))+ #
  geom_bar(stat="identity",position="stack")+
  theme_bw()+ ylab("HLA.B27") +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())
ggarrange(a1,a2,a3,a4,a5,nrow = 2,ncol = 3)

## Fig 3 B,C,E,F,D
test_label <- c()
result_test <- result_test[rownames(test)]
cutoff = 0.328 ## training cutoff
for (i in 1:length(result_test)){
  if (result_test[i] > cutoff){
    test_label <- c(test_label,"high")
  }
  else {
    test_label <- c(test_label,"low")
  }
}
blue <- "#00676B"
red <- "#BC3C28"
test_surv <- survfit(Surv(test$relapse_time,test$relapse_event)~test_label,data = test)
test_surv
survdiff(Surv(test$relapse_time,test$relapse_event)~test_label,data = test)
summary(test_surv)
ggsurvplot(test_surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#添加置信区间
           legend.labs=c("high","low")
)
# HR
summary(coxph(Surv(test$relapse_time,test$relapse_event) ~ result_test))  ## HR
# 
data2 <- data.frame(score = result_test, 
                    cluster = as.factor(test$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  ylab("Predict score") # 
###   heatmap 
test_heat <- cbind(data[rownames(test),3],test_label,test[,-1])
test_heat <- test_heat[names(result_test[order(result_test)]),]
result_test <- result_test[rownames(test_heat)]
test_heat_cli <- test_heat[,1:5]
test_heat_exp <- cbind(test_heat[,6:8],result_test)
test_heat_exp[,1] <- scale(test_heat_exp[,1])
test_heat_exp[,2] <- scale(test_heat_exp[,2])
test_heat_exp[,3] <- scale(test_heat_exp[,3])
#test_heat_exp[,4] <- scale(test_heat_exp[,4])
aa <- test_heat_exp[rownames(test_heat_cli)[which(test_heat_cli$test_label == "high")],]
bb <- test_heat_exp[rownames(test_heat_cli)[which(test_heat_cli$test_label == "low")],]
aa <- aa[rev(order(aa[,4])),]
bb <- bb[rev(order(bb[,4])),]
test_heat_exp <- as.matrix(rbind(aa,bb))
annotation_col = data.frame( gender = test_heat_cli[,1],
                             predict_label = test_heat_cli[,2],
                             relapse_event = as.factor(test_heat_cli[,3]),
                             Ankylosing_spondylitis = as.factor(test_heat_cli[,4]),
                             HLA.B27 = as.factor(test_heat_cli[,5])
)
rownames(annotation_col) <- rownames(test_heat_cli)
ann_colors = list(
  predict_label = c( "low"= "#00676B", "high" = "#BC3C28"),
  gender = c( "male"= "#7570B3", "female" = "#D2A6C7"),
  relapse_event = c("0" = "#5F98C6","1" = "#FF9933"),
  Ankylosing_spondylitis = c("0" = "#FBB4AE","1" = "#BC805E"),
  HLA.B27 = c("0" = "#B1C7BA","1" = "#E6BD9F")
)  
red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
hist(test_heat_exp)
test_heat_exp[test_heat_exp>2] <- 2
test_heat_exp[test_heat_exp<(-2)] <- c(-2)
pheatmap(t(test_heat_exp[,1:3]),fontsize=6,cutree_col = 4,cellheight = 25,cellwidth = 5 ,
         color  = colorRampPalette(c(blue,white,red))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T
)  #
pheatmap(t(test_heat_exp[,4]),fontsize=6,cutree_col = 4,cellheight = 25,cellwidth = 5 ,
         color  = colorRampPalette(c(white,"#0B0B0B"))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T
)  #
test_fufa <- test_heat[which(test_heat$relapse_event == 1),]
test_bufufa <- test_heat[which(test_heat$relapse_event == 0),]
table(test_fufa$HLA.B27)
a <- matrix(c(16,3,24,12),2,2)
chisq.test(a)
### 
test_heat_exp <- test_heat[,6:8]
result_test <- result_test[rownames(test_heat_exp)]
data22 <- as.data.frame(cbind(result_test,test_heat_exp))
a1 <- ggplot(data = data22 , aes(x = MO, y = result_test)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + 
  ylab("Score") + theme_bw()+ # 
  stat_cor(method = "pearson", label.x = 0, label.y = 0.6) +
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #
        axis.text.y=element_text(size=16,colour="black"), #
        axis.title.y=element_text(size = 20,colour="black"), #
        legend.text=element_text( colour="black",  #
                                  size=16),
        legend.title=element_text( colour="black", #
                                   size=18),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  theme_bw()
a2 <- ggplot(data = data22 , aes(x = HDL, y = result_test)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + 
  ylab("Score") + theme_bw()+ # 
  stat_cor(method = "pearson", label.x = 0, label.y = 0.6) +
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #
        axis.text.y=element_text(size=16,colour="black"), #
        axis.title.y=element_text(size = 20,colour="black"), #
        legend.text=element_text( colour="black",  #
                                  size=16),
        legend.title=element_text( colour="black", #
                                   size=18),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  theme_bw()
a3 <- ggplot(data = data22 , aes(x = LDL, y = result_test)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + ## 
  ylab("Score") + theme_bw()+ # 
  stat_cor(method = "pearson", label.x = 1, label.y = 0.6) + ## 
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #
        axis.text.y=element_text(size=16,colour="black"), 
        axis.title.y=element_text(size = 20,colour="black"), #
        legend.text=element_text( colour="black",  #
                                  size=16),
        legend.title=element_text( colour="black", #
                                   size=18),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  theme_bw()
ggarrange(a1,a2,a3,nrow = 2,ncol = 2)
## 
setwd("F:\\AAU回顾性研究")
test <- read.csv("test_clinical.csv",stringsAsFactors = F,row.names = 1)
a1 <- ggplot(data = test, aes(x = test_label , y = MO,fill=test_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())  #
a2 <- ggplot(data = test, aes(x = test_label , y = HDL,fill=test_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())  #
a3 <- ggplot(data = test, aes(x = test_label , y = LDL,fill=test_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())  #
a= "#00AE72"
b= "#5175A4"
c= "#D94E48"
d= "#8C63A4"
data <- data.frame(type = test$test_label,
                   Ank = test$Ankylosing.spondylitis,
                   clu1 = paste(test$test_label,test$Ankylosing.spondylitis,sep = ""),
                   clu2 = paste(test$test_label,test$HLA.B27,sep = "") )
data$number <- 1
data <- ddply(data,'type',transform,percent = 1/sum(number)*100)
a4 <- ggplot(data,aes(type,percent,fill=clu1))+
  scale_fill_manual(values = c(a,b,c,d))+ #
  geom_bar(stat="identity",position="stack")+
  theme_bw()+ ylab("Ankylosing.spondylitis") +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())
a5 <- ggplot(data,aes(type,percent,fill=clu2))+
  scale_fill_manual(values = c(a,b,c,d))+ #
  geom_bar(stat="identity",position="stack")+
  theme_bw()+ ylab("HLA.B27") +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())
ggarrange(a1,a2,a3,a4,a5,nrow = 2,ncol = 3)

## Fig 4 
setwd("F:\\AAU")
data1 <- read.csv("data_new.csv",stringsAsFactors = F,row.names = 1)
data2 <- read.csv("AUU_分组.csv",stringsAsFactors = F,row.names = 1)
data2 <- data2[rownames(data1),]
data <- data.frame(relapse_time = c(data2$relapse_time,data2$relapse_time,data2$relapse_time),
                   relapse_event = c(data2$relapse_event,data2$relapse_event,data2$relapse_event),
                   label = c(paste("Ank",data2$Ankylosing.spondylitis),paste("HLA",data2$HLA.B27),data1$label))
data_surv <- survfit(Surv(data$relapse_time,data$relapse_event)~label,data = data)
data_surv
survdiff(Surv(data$relapse_time,data$relapse_event)~label,data = data)
summary(data_surv)
ggsurvplot(data_surv,data=data,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E6BD9F","#8C63A4","#D94E48","#00AE72","#FB9900","#5175A4"),
           conf.int = F#
           #legend.labs=c("high","low")
)
data_new1 <- data[c(which(data$label == "HLA 0"),which(data$label == "HLA 1")),]
data_new2 <- data[c(which(data$label == "Ank 0"),which(data$label == "Ank 1")),]
data_new3 <- data[c(which(data$label == "low"),which(data$label == "high")),]
survdiff(Surv(data_new1$relapse_time,data_new1$relapse_event)~label,data = data_new1)
survdiff(Surv(data_new2$relapse_time,data_new2$relapse_event)~label,data = data_new2)
survdiff(Surv(data_new3$relapse_time,data_new3$relapse_event)~label,data = data_new3)

score = data1$Ankylosing.spondylitis*0.09230 + data1$HLA.B27*0.19863 + 
  data1$MO*(-0.59456) + data1$HDL*0.36348 + data1$LDL*(-0.12934) +  0.3287

roc1 <- roc(data1$relapse_event,score, plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
roc2 <- roc(data2$relapse_event,data2$Ankylosing.spondylitis, plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
roc3 <- roc(data2$relapse_event,data2$HLA.B27, plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
plot(roc1,print.thres=TRUE, ci=TRUE,print.auc=TRUE,legacy.axes = TRUE,col = "#BC3C28")
plot(roc2,add=T,print.thres=TRUE, ci=TRUE,print.auc=TRUE,legacy.axes = TRUE,col = "#8C63A4")
plot(roc3,add=T,print.thres=TRUE, ci=TRUE,print.auc=TRUE,legacy.axes = TRUE,col = "#FB9900")

roc.test(roc1,roc2)
roc.test(roc1,roc3)
roc.test(roc2,roc3)

####  
setwd("F:\\AAU")
data1 <- read.csv("data_new.csv",stringsAsFactors = F,row.names = 1)
data2 <- read.csv("AUU_分组.csv",stringsAsFactors = F,row.names = 1)
data2 <- data2[rownames(data1),]
data1 <- cbind(data1,data2$Ankylosing.spondylitis,data2$HLA.B27)
data1 <- data1[,-10]
data_Ank0 <- data1[which(data1$`data2$Ankylosing.spondylitis` == 0),]
data_Ank1 <- data1[which(data1$`data2$Ankylosing.spondylitis` == 1),]
data_HLA0 <- data1[which(data1$`data2$HLA.B27` == 0),]
data_HLA1 <- data1[which(data1$`data2$HLA.B27` == 1),]
## Ank0 
score = data_Ank0$Ankylosing.spondylitis*0.09230 + data_Ank0$HLA.B27*0.19863 + 
  data_Ank0$MO*(-0.59456) + data_Ank0$HDL*0.36348 + data_Ank0$LDL*(-0.12934) +  0.3287
data_Ank0_label <- c()
cutoff <- 0.328
for (i in 1:length(score)){
  if (score[i] > cutoff){
    data_Ank0_label <- c(data_Ank0_label,"high")
  }
  else {
    data_Ank0_label <- c(data_Ank0_label,"low")
  }
}
blue <- "#00676B"
red <- "#BC3C28"
data_Ank0_label_surv <- survfit(Surv(data_Ank0$relapse_time,data_Ank0$relapse_event)~data_Ank0_label,data = data_Ank0)
data_Ank0_label_surv
survdiff(Surv(data_Ank0$relapse_time,data_Ank0$relapse_event)~data_Ank0_label,data = data_Ank0)
summary(data_Ank0_label_surv)
ggsurvplot(data_Ank0_label_surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#
           legend.labs=c("high","low")
)
# HR
summary(coxph(Surv(data_Ank0$relapse_time,data_Ank0$relapse_event) ~ score))  ## HR
# 
data2 <- data.frame(score = data_Ank0, 
                    cluster = as.factor(data_Ank0$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  ylab("Predict score") # 
## Ank1 
score = data_Ank1$Ankylosing.spondylitis*0.09230 + data_Ank1$HLA.B27*0.19863 + 
  data_Ank1$MO*(-0.59456) + data_Ank1$HDL*0.36348 + data_Ank1$LDL*(-0.12934) +  0.3287
data_Ank1_label <- c()
cutoff <- 0.328
for (i in 1:length(score)){
  if (score[i] > cutoff){
    data_Ank1_label <- c(data_Ank1_label,"high")
  }
  else {
    data_Ank1_label <- c(data_Ank1_label,"low")
  }
}
blue <- "#00676B"
red <- "#BC3C28"
data_Ank1_label_surv <- survfit(Surv(data_Ank1$relapse_time,data_Ank1$relapse_event)~data_Ank1_label,data = data_Ank1)
data_Ank1_label_surv
survdiff(Surv(data_Ank1$relapse_time,data_Ank1$relapse_event)~data_Ank1_label,data = data_Ank1)
summary(data_Ank1_label_surv)
ggsurvplot(data_Ank1_label_surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#
           legend.labs=c("high","low")
)
# HR
summary(coxph(Surv(data_Ank1$relapse_time,data_Ank1$relapse_event) ~ score))  ## HR
# 
data2 <- data.frame(score = data_Ank1, 
                    cluster = as.factor(data_Ank1$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  ylab("Predict score") # 
## HLA0 的
score = data_HLA0$Ankylosing.spondylitis*0.09230 + data_HLA0$HLA.B27*0.19863 + 
  data_HLA0$MO*(-0.59456) + data_HLA0$HDL*0.36348 + data_HLA0$LDL*(-0.12934) +  0.3287
data_HLA0_label <- c()
cutoff <- 0.328
for (i in 1:length(score)){
  if (score[i] > cutoff){
    data_HLA0_label <- c(data_HLA0_label,"high")
  }
  else {
    data_HLA0_label <- c(data_HLA0_label,"low")
  }
}
blue <- "#00676B"
red <- "#BC3C28"
data_HLA0_label_surv <- survfit(Surv(data_HLA0$relapse_time,data_HLA0$relapse_event)~data_HLA0_label,data = data_HLA0)
data_HLA0_label_surv
survdiff(Surv(data_HLA0$relapse_time,data_HLA0$relapse_event)~data_HLA0_label,data = data_HLA0)
summary(data_HLA0_label_surv)
ggsurvplot(data_HLA0_label_surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#
           legend.labs=c("high","low")
)
# HR
summary(coxph(Surv(data_HLA0$relapse_time,data_HLA0$relapse_event) ~ score))  ## HR
# 
data2 <- data.frame(score = data_HLA0, 
                    cluster = as.factor(data_HLA0$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  ylab("Predict score") # 
## HLA1
score = data_HLA1$Ankylosing.spondylitis*0.09230 + data_HLA1$HLA.B27*0.19863 + 
  data_HLA1$MO*(-0.59456) + data_HLA1$HDL*0.36348 + data_HLA1$LDL*(-0.12934) +  0.3287
data_HLA1_label <- c()
cutoff <- 0.328
for (i in 1:length(score)){
  if (score[i] > cutoff){
    data_HLA1_label <- c(data_HLA1_label,"high")
  }
  else {
    data_HLA1_label <- c(data_HLA1_label,"low")
  }
}
blue <- "#00676B"
red <- "#BC3C28"
data_HLA1_label_surv <- survfit(Surv(data_HLA1$relapse_time,data_HLA1$relapse_event)~data_HLA1_label,data = data_HLA1)
data_HLA0_label_surv
survdiff(Surv(data_HLA1$relapse_time,data_HLA1$relapse_event)~data_HLA1_label,data = data_HLA1)
summary(data_HLA1_label_surv)
ggsurvplot(data_HLA1_label_surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#
           legend.labs=c("high","low")
)
# HR
summary(coxph(Surv(data_HLA1$relapse_time,data_HLA1$relapse_event) ~ score))  ## HR
# 
data2 <- data.frame(score = data_HLA1, 
                    cluster = as.factor(data_HLA1$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #
  theme_bw()+ #
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #
        axis.text.y=element_text(size=10,colour="black"), #
        axis.title.y=element_text(size = 10,colour="black"), #
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  ylab("Predict score") # 



























