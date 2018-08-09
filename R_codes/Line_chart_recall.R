D_recall <-  read.table("30X_D_Recall2.txt")

D_data_recall <- vector()
for (i in 2:6)
{
  for (j in 1:10)
  {
    item <- vector()
    item <- c(as.character(D_recall[i,1]),as.numeric(as.character(D_recall[i,j+1])),-j*200+100)
    D_data_recall <- rbind(D_data_recall,item)
  }
}

D_data_recall <- data.frame(D_data_recall)
#library(ggplot2)
#sp <- ggplot(data=D_data_recall, aes(x=as.numeric(as.character(X3)), y=as.numeric(as.character(X2)), group=X1)) +
#  geom_line(aes(color=X1),size=1.2)+
#  geom_point(aes(color=X1))
#sp + ylim(0, 1)
#sp+scale_color_manual(values=c("#0d61f2","#5d5d35","#00f6dc","#ca2116","#da2491","#f8c30e","#31d3ae","#fa8072"))


I_recall <-  read.table("30X_I_Recall2.txt",header = FALSE)

I_data_recall <- vector()
for (i in 2:6)
{
  for (j in 1:10)
  {
    item <- vector()
    item <- c(as.character(I_recall[i,1]),as.numeric(as.character(I_recall[i,j+1])),(j-1)*200+100)
    I_data_recall <- rbind(I_data_recall,item)
  }
}

I_data_recall <- data.frame(I_data_recall)
#library(ggplot2)
#sp <- ggplot(data=I_data_recall, aes(x=as.numeric(as.character(X3)), y=as.numeric(as.character(X2)), group=X1)) +
#  geom_line(aes(color=X1),size=1.2)+
#  geom_point(aes(color=X1))
#sp+scale_color_manual(values=c("#0d61f2","#5d5d35","#00f6dc","#ca2116","#da2491","#f8c30e","#31d3ae","#fa8072"))
#sp + ylim(0, 1)


Data_Recall <- rbind(D_data_recall,I_data_recall)
Data_Recall <- Data_Recall[order(Data_Recall$X1),]
colnames(Data_Recall) <- c("Callers","recall","size")
library(ggplot2)
sp_recall <- ggplot(data=Data_Recall, aes(x=as.numeric(as.character(Data_Recall$size)), 
                                          y=as.numeric(as.character(Data_Recall$recall)), group=Data_Recall$Caller)) +
  geom_line(aes(color=Data_Recall$Caller,linetype=Data_Recall$Caller),size=1.2)+
  geom_point(aes(color=Data_Recall$Caller))

sp_recall <- sp_recall+ggtitle("Indel Callers Recall")+
  labs(y="Recall", x = "Indel Size (bp)",size=18)#+labs(colour = "Callers")

sp_recall <-sp_recall+
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

#sp_recall <-sp_recall+theme(legend.title=element_blank())+
#  theme(legend.position = c(0.95,0.3)),legend.background = element_rect(fill=alpha('transparent', 0)))

sp_recall <-sp_recall+scale_linetype_manual(values=c("solid","twodash","dashed","solid",
                                                     "dotdash","longdash","twodash","solid","solid","dashed"))

sp_recall <- sp_recall+scale_color_manual(values=c("#0d61f2","#ca2116",
                                                   "#f8c30e","#00f6dc","#7CFC00"))

#sp_recall <-sp_recall + guides(color=guide_legend(override.aes=list(fill=NA)))
sp_recall <- sp_recall+ylim(0,1)
sp_recall



















