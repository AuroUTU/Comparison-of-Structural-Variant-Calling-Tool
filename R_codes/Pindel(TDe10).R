# Get TD from simulated data
CAL <- read.table("CAL.txt")
TD_truth <- CAL[CAL[,1]=="duplication",]


TDe2 <- read.table("Pindel_simulated_TDe2.vcf")
TDe5 <- read.table("Pindel_simulated_TDe5.vcf")
TDe10 <- read.table("Pindel_simulated_TD.vcf")
TDe15 <- read.table("Pindel_simulated_TDe15.vcf")
TDe20 <- read.table("Pindel_simulated_TDe20.vcf")

# make a function
parse_length <- function(P_input,truthset)
{
  position_in <- P_input[,2]
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

parse_length(TDe2,TD_truth)
parse_length(TDe5,TD_truth)
parse_length(TDe10,TD_truth)
parse_length(TDe15,TD_truth)
parse_length(TDe20,TD_truth)


dat <- read.table(text = "2   5  10   15   20 
TP 385 385 384 383 378
FP 202  40  10  2  1", header = TRUE)
barplot(as.matrix(dat),legend = rownames(dat),col=c("darkblue","red"),
        main="Tandem Duplication FP rate ",xlab="min_read_support",ylab="Total variants",
        cex.lab=1.5,cex.main=1.5)s
text(5.8,460, "Truth=400",cex = 1.2)
text(5.7,440, "Precision/Recall",cex = 1)

#add Precision and Recall rate
text(1.7,560, "65.5%/96.2%",cex = 1.2)
text(1.9,440, "90.5%/96.2%",cex = 1.2)
text(3.1,410, "97.4%/96.0%",cex = 1.2)
text(4.3,400, "99.4%/95.7%",cex = 1.2)
text(5.5,400, "100%/94.5%",cex = 1.2)

abline(h=400,col="green",lwd=3)




# split into sunsets
Pindel_dup_data <- vector()
Pindel_dup <- TDe10
for (i in 1:dim(Pindel_dup)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- Pindel_dup[i,2]
  Pindel_dup_info<- strsplit(as.character(Pindel_dup[i,8]),";")[[1]]
  Pindel_dup_len <- Pindel_dup_info[grep("SVLEN=",Pindel_dup_info)]
  Pindel_dup_LEN <- substr(as.character(Pindel_dup_len),7,nchar(Pindel_dup_len)) 
  item <- c(position,Pindel_dup_LEN)
  Pindel_dup_data <- rbind(Pindel_dup_data,item)
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
    print (c(i*200,(i+1)*200))
    #print ("Total TP:")
    #print (length(TP)) # how many TP
    #print ("Total FP:")
    #print (length(position_in)-length(TP))
    #print ("")
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

parse_length(Pindel_dup_data,TD_truth)




parse_length(TDe10,TD_truth)




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
parse_length(Pindel_dup_data,TD_truth)









