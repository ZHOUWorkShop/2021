#######   AAU  ########
## 总体数据的热图
setwd("F:\\AAU回顾性研究")
data <- read.csv("AUU_连续.csv",stringsAsFactors = F,row.names = 1)
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
)  #聚类热图
all_fufa <- data[which(data$relapse_event == 1),]
all_bufufa <- data[which(data$relapse_event == 0),]
a <- matrix(c(55,45,105,28),2,2)
chisq.test(a)

## 连续值
setwd("F:\\AAU回顾性研究")
data <- read.csv("AUU_连续.csv",stringsAsFactors = F,row.names = 1)
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
#write.csv(cox,"连续单cox森林图.csv",quote = F)
data <- read.csv("连续单cox森林图.csv",stringsAsFactors = F)
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



## 离散值
setwd("F:\\AAU回顾性研究")
data <- read.csv("AUU_分组.csv",stringsAsFactors = F,row.names = 1)
cox <- matrix(,1,7)
colnames(cox) <- c('lower .95','upper .95','coef', 'exp(coef)','se(coef)','z','Pr(>|z|)' )
name <- c()
for (i in 3:46){
  a <- summary(coxph(Surv(data$relapse_time,data$relapse_event) ~ data[,i]))
  
  if (length(a$conf.int[,3:4]) == 2){
    cox <- rbind(cox,c(a$conf.int[,3:4],a$coefficients[,1:5]))
    name <- c(name,colnames(data)[i])
  }
  else {
    cox <- rbind(cox,cbind(a$conf.int[,3:4],a$coefficients[,1:5]))
    name <- c(name,paste(colnames(data)[i],substr(rownames(a$conf.int),10,11),sep = ": "))
  }
}
cox <- cox[-1,]
rownames(cox) <- name
cox[,7] <- round(cox[,7],3)
cox <- cox[order(cox[,7]),]
cox <- cbind(cox,paste(round(cox[,4],3)," (",round(cox[,1],3)," - ",round(cox[,2],3),")",sep = ""))
#write.csv(cox,"分组单cox森林图.csv",quote = F)
data <- read.csv("分组单cox森林图.csv",stringsAsFactors = F)
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

### 复发和不复发组之间的 t-test
setwd("F:\\AAU回顾性研究")
data <- read.csv("AUU_连续.csv",stringsAsFactors = F,row.names = 1)
fufa <- data[which(data$relapse_event == 1),4:46]
bufufa <- data[which(data$relapse_event == 0),4:46]
fufa <- fufa[,-19]
bufufa <- bufufa[,-19]
jieguo <- matrix(,,4)
for (i in 1:dim(fufa)[2]){
  a <- t.test(as.numeric(fufa[,i]),as.numeric(bufufa[,i]))
  jieguo <- rbind(jieguo,c(colnames(fufa)[i],
                           length(which(is.na(fufa[,i]) == F)),
                           length(which(is.na(bufufa[,i]) == F)),
                           a$p.value
                           )
                  )
}
jieguo <- jieguo[-1,]
jieguo <- jieguo[order(jieguo[,4]),]
colnames(jieguo) <- c("Infomation","relapse=1","relapse=0","p-value")
#write.csv(jieguo,"t_test_result.csv",quote = F)


##############################################################
#####             分 训练 + 测试进行模型构建             #####
##############################################################
setwd("F:\\AAU回顾性研究")
data <- read.csv("data_new.csv",stringsAsFactors = F,row.names = 1)
data <- data[,-4]
data_model <- data[,c(2:8)]
data_model <- data_model[,-4]
#n <- dim(data_model)[1]-1    #计算数据集中自变量个数
#rate=1     #设置模型误判率向量初始值
##for(i in 1:n){
#  set.seed(1234)
#  rf_train<-randomForest(as.factor(data_model$relapse_event)~.,data=data_model,mtry=i,ntree=1000)
#  rate[i]<-mean(rf_train$err.rate)   #计算基于OOB数据的模型误判率均值
#  print(rf_train)    
#}
#rate     #展示所有模型误判率的均值
#plot(rate)
#set.seed(1234)
#rf_train<-randomForest(as.factor(data_model$relapse_event)~.,data=data_model,mtry=2,ntree=300)
#plot(rf_train)    #绘制模型误差与决策树数量关系图  

set.seed(1234)
otu_forest <- randomForest(as.factor(data_model$relapse_event)~., data = data_model,
                           ,mtry=2,ntree=200,importance=TRUE,proximity=TRUE) 
importance<-importance(otu_forest) 

varImpPlot(otu_forest, n.var = min(10, nrow(otu_forest$importance)),
           main = 'Importance')
#write.csv(importance,"RF_Importance.csv")

###  根据重要性 和 cox结果 选择构建模型特征  ###
setwd("F:\\AAU回顾性研究")
data <- read.csv("data_new.csv",stringsAsFactors = F,row.names = 1)
data_model <- data[,c(-4,-6,-10)]
fufa <- data_model[which(data_model$relapse_event == 1),]
bufufa <- data_model[which(data_model$relapse_event == 0),]

set.seed(396)
a <- c(1:30)
b <- c(1:81)
train_a <- sample(1:30,15)
train_b <- sample(1:81,41)
test_a <- a[-train_a]
test_b <- b[-train_b]
train <- rbind(fufa[train_a,],bufufa[train_b,])
test <- rbind(fufa[test_a,],bufufa[test_b,])

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

###  KM曲线  ###
##################### 训练集  ###################
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
           pval = TRUE,#添加P值
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#添加置信区间
           legend.labs=c("high","low")
           )
# 计算HR
summary(coxph(Surv(train$relapse_time,train$relapse_event) ~ result_train))  ## HR
# 绘制score的boxplot
data2 <- data.frame(score = result_train, 
                    cluster = as.factor(train$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Predict score") #设置y轴的标题 
###  绘制 heatmap 
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
)  #聚类热图
pheatmap(t(train_heat_exp[,4]),fontsize=6,cutree_col = 4,cellheight = 25,cellwidth = 5 ,
         color  = colorRampPalette(c(white,"#0B0B0B"))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T
)  #聚类热图
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
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=16,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 20,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text( colour="black",  #设置图例的子标题的字体属性
                                  size=16),
        legend.title=element_text( colour="black", #设置图例的总标题的字体属性
                                   size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  theme_bw()
a2 <- ggplot(data = data22 , aes(x = HDL, y = result_train)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + 
  ylab("Score") + theme_bw()+ #背景变为白色 
  stat_cor(method = "pearson", label.x = 0, label.y = 0.6) +
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=16,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 20,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text( colour="black",  #设置图例的子标题的字体属性
                                  size=16),
        legend.title=element_text( colour="black", #设置图例的总标题的字体属性
                                   size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  theme_bw()
a3 <- ggplot(data = data22 , aes(x = LDL, y = result_train)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + ## colour控制线的颜色，fill控制置信区间的颜色
  ylab("Score") + theme_bw()+ #背景变为白色 
  stat_cor(method = "pearson", label.x = 1, label.y = 0.6) + ## 添加拟合结果
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=16,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 20,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text( colour="black",  #设置图例的子标题的字体属性
                                  size=16),
        legend.title=element_text( colour="black", #设置图例的总标题的字体属性
                                   size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  theme_bw()
ggarrange(a1,a2,a3,nrow = 2,ncol = 2)
## 两组间5个因素的比较
setwd("F:\\AAU回顾性研究")
train <- read.csv("train_clinical.csv",stringsAsFactors = F,row.names = 1)
a1 <- ggplot(data = train, aes(x = train_label , y = MO,fill=train_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
a2 <- ggplot(data = train, aes(x = train_label , y = HDL,fill=train_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
a3 <- ggplot(data = train, aes(x = train_label , y = LDL,fill=train_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
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
  scale_fill_manual(values = c(a,b,c,d))+ #设置填充的颜色
  geom_bar(stat="identity",position="stack")+
  theme_bw()+ ylab("Ankylosing.spondylitis") +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())
a5 <- ggplot(data,aes(type,percent,fill=clu2))+
  scale_fill_manual(values = c(a,b,c,d))+ #设置填充的颜色
  geom_bar(stat="identity",position="stack")+
  theme_bw()+ ylab("HLA.B27") +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())
ggarrange(a1,a2,a3,a4,a5,nrow = 2,ncol = 3)

###########################  测试集  ##################
test_label <- c()
result_test <- result_test[rownames(test)]
cutoff = 0.328 ## 训练集的cutoff
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
           pval = TRUE,#添加P值
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#添加置信区间
           legend.labs=c("high","low")
)
# 计算HR
summary(coxph(Surv(test$relapse_time,test$relapse_event) ~ result_test))  ## HR
# 绘制score的boxplot
data2 <- data.frame(score = result_test, 
                    cluster = as.factor(test$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Predict score") #设置y轴的标题 
###  绘制 heatmap 
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
)  #聚类热图
pheatmap(t(test_heat_exp[,4]),fontsize=6,cutree_col = 4,cellheight = 25,cellwidth = 5 ,
         color  = colorRampPalette(c(white,"#0B0B0B"))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T
)  #聚类热图
test_fufa <- test_heat[which(test_heat$relapse_event == 1),]
test_bufufa <- test_heat[which(test_heat$relapse_event == 0),]
table(test_fufa$HLA.B27)
a <- matrix(c(16,3,24,12),2,2)
chisq.test(a)
### 相关性分析
test_heat_exp <- test_heat[,6:8]
result_test <- result_test[rownames(test_heat_exp)]
data22 <- as.data.frame(cbind(result_test,test_heat_exp))
a1 <- ggplot(data = data22 , aes(x = MO, y = result_test)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + 
  ylab("Score") + theme_bw()+ #背景变为白色 
  stat_cor(method = "pearson", label.x = 0, label.y = 0.6) +
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=16,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 20,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text( colour="black",  #设置图例的子标题的字体属性
                                  size=16),
        legend.title=element_text( colour="black", #设置图例的总标题的字体属性
                                   size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  theme_bw()
a2 <- ggplot(data = data22 , aes(x = HDL, y = result_test)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + 
  ylab("Score") + theme_bw()+ #背景变为白色 
  stat_cor(method = "pearson", label.x = 0, label.y = 0.6) +
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=16,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 20,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text( colour="black",  #设置图例的子标题的字体属性
                                  size=16),
        legend.title=element_text( colour="black", #设置图例的总标题的字体属性
                                   size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  theme_bw()
a3 <- ggplot(data = data22 , aes(x = LDL, y = result_test)) + 
  geom_point(colour = "#426671", size = 1) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7') + ## colour控制线的颜色，fill控制置信区间的颜色
  ylab("Score") + theme_bw()+ #背景变为白色 
  stat_cor(method = "pearson", label.x = 1, label.y = 0.6) + ## 添加拟合结果
  theme(axis.text.x=element_text(angle = 45,hjust = 1,colour="black",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=16,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 20,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text( colour="black",  #设置图例的子标题的字体属性
                                  size=16),
        legend.title=element_text( colour="black", #设置图例的总标题的字体属性
                                   size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  theme_bw()
ggarrange(a1,a2,a3,nrow = 2,ncol = 2)
## 两组间5个因素的比较
setwd("F:\\AAU回顾性研究")
test <- read.csv("test_clinical.csv",stringsAsFactors = F,row.names = 1)
a1 <- ggplot(data = test, aes(x = test_label , y = MO,fill=test_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
a2 <- ggplot(data = test, aes(x = test_label , y = HDL,fill=test_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
a3 <- ggplot(data = test, aes(x = test_label , y = LDL,fill=test_label)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#BC3C28","#00676B"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
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
  scale_fill_manual(values = c(a,b,c,d))+ #设置填充的颜色
  geom_bar(stat="identity",position="stack")+
  theme_bw()+ ylab("Ankylosing.spondylitis") +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())
a5 <- ggplot(data,aes(type,percent,fill=clu2))+
  scale_fill_manual(values = c(a,b,c,d))+ #设置填充的颜色
  geom_bar(stat="identity",position="stack")+
  theme_bw()+ ylab("HLA.B27") +
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())
ggarrange(a1,a2,a3,a4,a5,nrow = 2,ncol = 3)


####  分层分析, 在全体数据上
setwd("F:\\AAU回顾性研究")
data1 <- read.csv("data_new.csv",stringsAsFactors = F,row.names = 1)
data2 <- read.csv("AUU_分组.csv",stringsAsFactors = F,row.names = 1)
data2 <- data2[rownames(data1),]
data1 <- cbind(data1,data2$Ankylosing.spondylitis,data2$HLA.B27)
data1 <- data1[,-10]
data_Ank0 <- data1[which(data1$`data2$Ankylosing.spondylitis` == 0),]
data_Ank1 <- data1[which(data1$`data2$Ankylosing.spondylitis` == 1),]
data_HLA0 <- data1[which(data1$`data2$HLA.B27` == 0),]
data_HLA1 <- data1[which(data1$`data2$HLA.B27` == 1),]
## Ank0 的分层
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
           pval = TRUE,#添加P值
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#添加置信区间
           legend.labs=c("high","low")
)
# 计算HR
summary(coxph(Surv(data_Ank0$relapse_time,data_Ank0$relapse_event) ~ score))  ## HR
# 绘制score的boxplot
data2 <- data.frame(score = data_Ank0, 
                    cluster = as.factor(data_Ank0$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Predict score") #设置y轴的标题 
## Ank1 的分层
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
           pval = TRUE,#添加P值
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#添加置信区间
           legend.labs=c("high","low")
)
# 计算HR
summary(coxph(Surv(data_Ank1$relapse_time,data_Ank1$relapse_event) ~ score))  ## HR
# 绘制score的boxplot
data2 <- data.frame(score = data_Ank1, 
                    cluster = as.factor(data_Ank1$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Predict score") #设置y轴的标题 
## HLA0 的分层
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
           pval = TRUE,#添加P值
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#添加置信区间
           legend.labs=c("high","low")
)
# 计算HR
summary(coxph(Surv(data_HLA0$relapse_time,data_HLA0$relapse_event) ~ score))  ## HR
# 绘制score的boxplot
data2 <- data.frame(score = data_HLA0, 
                    cluster = as.factor(data_HLA0$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Predict score") #设置y轴的标题 
## HLA1 的分层
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
           pval = TRUE,#添加P值
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(red,blue),
           conf.int = F,#添加置信区间
           legend.labs=c("high","low")
)
# 计算HR
summary(coxph(Surv(data_HLA1$relapse_time,data_HLA1$relapse_event) ~ score))  ## HR
# 绘制score的boxplot
data2 <- data.frame(score = data_HLA1, 
                    cluster = as.factor(data_HLA1$relapse_event)) 
ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ #绘制箱线图
  scale_fill_manual(values = c("#00676B","#BC3C28"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Predict score") #设置y轴的标题 

















### 配对boxplot
setwd("F:\\AAU回顾性研究")
data <- read.csv("连续血样检验指标配对.csv",stringsAsFactors = F)
# Mon
name <- data[which(is.na(data$Mon) == T),1]
weizhi <- c()
for (i in 1:length(name)){
  weizhi <- c(weizhi,which(name[i] == data[,1]))
}
data1 <- data[-weizhi,]
p1 <- ggpaired(data1, x = "label", y = "Mon",xlab = "Mon",
               color = "label", palette = "jama", 
               line.color = "gray", line.size = 0.4,
                short.panel.labs = FALSE) + stat_compare_means(label = "p.format", paired = TRUE)
# TG
p2 <- ggpaired(data, x = "label", y = "TG",xlab = "TG",
               color = "label", palette = "jama", 
               line.color = "gray", line.size = 0.4,
               short.panel.labs = FALSE) + stat_compare_means(label = "p.format", paired = TRUE)
# TC
p3 <- ggpaired(data, x = "label", y = "TC",xlab = "TC",
               color = "label", palette = "jama", 
               line.color = "gray", line.size = 0.4,
               short.panel.labs = FALSE) + stat_compare_means(label = "p.format", paired = TRUE)
# HDL
p4 <- ggpaired(data, x = "label", y = "HDL",xlab = "HDL",
               color = "label", palette = "jama", 
               line.color = "gray", line.size = 0.4,
               short.panel.labs = FALSE) + stat_compare_means(label = "p.format", paired = TRUE)
# LDL
p5 <- ggpaired(data, x = "label", y = "LDL",xlab = "LDL",
               color = "label", palette = "jama", 
               line.color = "gray", line.size = 0.4,
               short.panel.labs = FALSE) + stat_compare_means(label = "p.format", paired = TRUE)
ggarrange(p1,p2,p3,p4,p5,nrow = 3,ncol = 2)



















