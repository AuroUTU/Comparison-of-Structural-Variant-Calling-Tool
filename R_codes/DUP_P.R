TD_precision <-  read.table("30X_TD_Precision2.txt",header = FALSE)

TD_data_Precision <- vector()
for (i in 2:6)
{
  for (j in 1:10)
  {
    item <- vector()
    item <- c(as.character(TD_precision[i,1]),as.numeric(as.character(TD_precision[i,j+1])),(j-1)*200+100)
    TD_data_Precision <- rbind(TD_data_Precision,item)
  }
}

TD_data_Precision <- data.frame(TD_data_Precision)
colnames(TD_data_Precision) <- c("Callers","precision","size")
library(ggplot2)
sp_precision <- ggplot(data=TD_data_Precision, 
                       aes(x=as.numeric(as.character(TD_data_Precision$size)),
                           y=as.numeric(as.character(TD_data_Precision$precision)), group=TD_data_Precision$Caller)) +
  geom_line(aes(color=TD_data_Precision$Caller,linetype=TD_data_Precision$Caller),size=1.2)+
  geom_point(aes(color=TD_data_Precision$Caller))

sp_precision <- sp_precision+ggtitle("Indel Callers TD Precision")+
  labs(y="Precision", x = "Indel Size (bp)",size=18)#+labs(colour = "Callers")

sp_precision <- sp_precision +
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

sp_precision <- sp_precision + scale_linetype_manual(values=c("solid","twodash","dashed","solid",
                                                              "dotdash","longdash","twodash","solid","solid","dashed"))

sp_precision <- sp_precision + scale_color_manual(values=c("#0d61f2","#ca2116",
                                                           "#f8c30e","#00f6dc","#7CFC00"))

sp_precision
