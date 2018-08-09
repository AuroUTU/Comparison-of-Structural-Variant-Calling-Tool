D_precision <-  read.table("30X_D_Precision2.txt",header = FALSE)

D_data_Precision <- vector()
for (i in 2:6)
{
  for (j in 1:10)
  {
    item <- vector()
    item <- c(as.character(D_precision[i,1]),as.numeric(as.character(D_precision[i,j+1])),-j*200+100)
    D_data_Precision <- rbind(D_data_Precision,item)
  }
}

D_data_Precision <- data.frame(D_data_Precision)
#library(ggplot2)
#sp <- ggplot(data=D_data_Precision, aes(x=as.numeric(as.character(X3)), y=as.numeric(as.character(X2)), group=X1)) +
#  geom_line(aes(color=X1),size=1.2)+
#  geom_point(aes(color=X1))
#sp+scale_color_manual(values=c("#0d61f2","#5d5d35","#00f6dc","#ca2116","#da2491","#f8c30e","#31d3ae","#fa8072"))
#sp + ylim(0, 1)

I_precision <-  read.table("30X_I_Precision2.txt",header = FALSE)

I_data_Precision <- vector()
for (i in 2:6)
{
  for (j in 1:10)
  {
    item <- vector()
    item <- c(as.character(I_precision[i,1]),as.numeric(as.character(I_precision[i,j+1])),(j-1)*200+100)
    I_data_Precision <- rbind(I_data_Precision,item)
  }
}

I_data_Precision <- data.frame(I_data_Precision)
#library(ggplot2)
#sp <- ggplot(data=I_data_Precision, aes(x=as.numeric(as.character(X3)), y=as.numeric(as.character(X2)), group=X1)) +
#  geom_line(aes(color=X1),size=1.2)+
#  geom_point(aes(color=X1))
#sp+scale_color_manual(values=c("#0d61f2","#5d5d35","#00f6dc","#ca2116","#da2491","#f8c30e","#31d3ae","#fa8072"))
#sp + ylim(0, 1)

Data_Precision <- rbind(D_data_Precision,I_data_Precision)
Data_Precision <- Data_Precision[order(Data_Precision$X1),]
colnames(Data_Precision) <- c("Callers","precision","size")
library(ggplot2)
sp_precision <- ggplot(data=Data_Precision, 
                       aes(x=as.numeric(as.character(Data_Precision$size)), 
                           y=as.numeric(as.character(Data_Precision$precision)), group=Data_Precision$Caller)) +
  geom_line(aes(color=Data_Precision$Caller,linetype=Data_Precision$Caller),size=1.2)+
  geom_point(aes(color=Data_Precision$Caller))

sp_precision <- sp_precision+ggtitle("Indel Callers Precision")+
  labs(y="Precision", x = "Indel Size (bp)",size=18)#+labs(colour = "Callers")

sp_precision <- sp_precision+
  theme(
    legend.key = element_rect(colour = "transparent", fill = NA),
    axis.text.x = element_text(color = "black",size=18),
    axis.text.y = element_text(color = "black",size=18),
    axis.title=element_text(size=18),
    title =element_text(size=15),
    legend.text=element_text(size=14),
    legend.title=element_blank(),
    legend.position = c(0.92,0.3),
    legend.background = element_rect(fill=alpha('transparent', 0)
    ))

sp_precision <- sp_precision+scale_linetype_manual(values=c("solid","twodash","dashed","solid",
                                                            "dotdash","longdash","twodash","solid","solid","dashed"))

sp_precision <- sp_precision+scale_color_manual(values=c("#0d61f2","#ca2116",
                                                         "#f8c30e","#00f6dc","#7CFC00"))

sp_precision <- sp_precision+ylim(0, 1)

sp_precision
