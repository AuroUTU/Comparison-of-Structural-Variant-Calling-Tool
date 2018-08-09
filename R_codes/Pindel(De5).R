# Get D from simulated data
CAL <- read.table("CAL.txt")
D_truth <- CAL[CAL[,1]=="deletion",]

De10 <- read.table("Pindel_simulated_D.vcf")


# make a function
parse_length <- function(P_input,truthset)
{
  position_in <- unique(as.numeric(P_input[,2]))
  position_truth <- truthset[,6]
  
  TP <- vector()
  #TP_sub <- intersect(position_in,position_truth)
  #position_truth <- position_truth[!(position_truth %in% TP_sub)]
  #TP <- c(TP,TP_sub)
  for (i in -5:5)
  {
    TP_sub <- intersect(position_in+i,position_truth)
    for (P in TP_sub)
    {
      P_input_info <- strsplit(as.character(P_input[P_input[,2]==P-i,8]),";")[[1]]
      P_input_len <- P_input_info[grep("SVLEN=",P_input_info)]
      P_input_LEN <- substr(as.character(P_input_len),7,nchar(P_input_len)) 
      if (abs(as.numeric(truthset[truthset[,6]==P,3])-abs(as.numeric(P_input_LEN))) < 10)
      {
        position_truth <- position_truth[!(position_truth %in% P)]
        TP <- c(TP,P)
      }
    }
  }
  print ("Total TP:")
  print (length(TP)) # how many TP
  print ("Total FP:")
  print (length(position_in)-length(TP))
  print ("")
  print("Precision")
  print(length(TP)/dim(P_input)[1])
  print("Recall")
  print((length(TP)/dim(truthset)[1]))
}

parse_length(De2,D_truth)
parse_length(De5,D_truth)
parse_length(De10,D_truth)
parse_length(De15,D_truth)
parse_length(De20,D_truth)


dat <- read.table(text = "2   5  10   15   20 
TP 790 790 791 787 784
FP 3224  102  25  12  8", header = TRUE)

barplot(as.matrix(dat),legend = rownames(dat),col=c("darkblue","red"),
        main="Deletion FP rate ",xlab="min_read_support",ylab="Total variants",
        cex.lab=1.5,cex.main=1.5)
text(5.8,3000, "Truth=800",cex = 1.2)
text(5.7,2800, "Precision/Recall",cex = 1)

#add Precision and Recall rate
text(1.7,3800, "19.6%/98.7%",cex = 1.2)
text(1.9,1100, "88.5%/98.7%",cex = 1.2)
text(3.1,1000, "96.9%/98.8%",cex = 1.2)
text(4.3,1000, "98.4%/98.3%",cex = 1.2)
text(5.5,1000, "98.9%/98.0%",cex = 1.2)

#add the total numner of truth
abline(h=800,col="green",lwd=3)






# split into subsets
Pindel_del_data <- vector()
Pindel_del <- De10
for (i in 1:dim(Pindel_del)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- Pindel_del[i,2]
  Pindel_del_info<- strsplit(as.character(Pindel_del[i,8]),";")[[1]]
  Pindel_del_len <- Pindel_del_info[grep("SVLEN=",Pindel_del_info)]
  Pindel_del_LEN <- substr(as.character(Pindel_del_len),8,nchar(Pindel_del_len)) # move "-" out
  item <- c(position,Pindel_del_LEN)
  Pindel_del_data <- rbind(Pindel_del_data,item)
}

# split dataset into 10 sbusets
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  for (i in 0:9)
  {
    truthset <- truthset_get[as.numeric(truthset_get[,3]) > i*200 & as.numeric(truthset_get[,3]) <= (i+1)*200,]
    P_input <- P_input_get[as.numeric(P_input_get[,2]) > i*200 & as.numeric(P_input_get[,2]) <= (i+1)*200,]
    position_in <- as.numeric(P_input[,1])
    position_truth <- truthset[,6]
    TP <- vector()
    #TP_sub <- intersect(position_in,position_truth)
    #position_truth <- position_truth[!(position_truth %in% TP_sub)]
    #TP <- c(TP,TP_sub)
    for (j in -5:5) # position tolerate
    {
      TP_sub <- intersect(as.numeric(position_in)+j,position_truth)
      for (P in TP_sub)
      {
        P_input_LEN <- P_input[P_input[,1]==P-j,2]
        if (abs(as.numeric(truthset[truthset[,6]==P,3])-abs(as.numeric(P_input_LEN))) < 10) # length tolerate
        {
          # remove match position, avoid getting repeat results, due to one indel maybe report twice near by the break point
          position_truth <- position_truth[!(position_truth %in% P)] 
          TP <- c(TP,P)
        }
      }
    }
    print (c(i*200,(i+1)*200))
    print ("Total TP:")
    print (length(TP)) # how many TP
    print ("Total FP:")
    print (length(position_in)-length(TP))
    print ("")
    print("Precision")
    print(length(TP)/length(position_in))
    precision <- c(precision,length(TP)/length(position_in))
    print("Recall")
    print((length(TP)/dim(truthset)[1]))
    recall <- c(recall,length(TP)/dim(truthset)[1])
  }
  print (precision)
  print (recall)
}

parse_length(Pindel_del_data,D_truth)




# Only parse 50bp--150bp
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  truthset <- truthset_get[as.numeric(truthset_get[,3]) > 50 & as.numeric(truthset_get[,3]) <= 150,]
  P_input <- P_input_get[as.numeric(P_input_get[,2]) > 50 & as.numeric(P_input_get[,2]) <= 150,]
  position_in <- unique(as.numeric(P_input[,1]))
  position_truth <- truthset[,6]
  TP <- vector()
  #TP_sub <- intersect(position_in,position_truth)
  #position_truth <- position_truth[!(position_truth %in% TP_sub)]
  #TP <- c(TP,TP_sub)
  for (j in -5:5) # position tolerate
  {
    TP_sub <- intersect(as.numeric(position_in)+j,position_truth)
    for (P in TP_sub)
    {
      P_input_LEN <- P_input[P_input[,1]==P-j,2]
      if (abs(as.numeric(truthset[truthset[,6]==P,3])-abs(as.numeric(P_input_LEN))) < 10) # length tolerate
      {
        # remove match position, avoid getting repeat results, due to one indel maybe report twice near by the break point
        position_truth <- position_truth[!(position_truth %in% P)] 
        TP <- c(TP,P)
      }
    }
  }
  print (c(50,150))
  print ("Total TP:")
  print (length(TP)) # how many TP
  print ("Total FP:")
  print (length(position_in)-length(TP))
  print ("")
  print("Precision")
  print(length(TP)/length(position_in))
  precision <- c(precision,length(TP)/length(position_in))
  print("Recall")
  print((length(TP)/dim(truthset)[1]))
  recall <- c(recall,length(TP)/dim(truthset)[1])
  print (precision)
  print (recall)
}
parse_length(Pindel_del_data,D_truth)






























