# in Pindel, the short insertion is less than 82

# Get SI from simulated data
CAL <- read.table("CAL.txt")
insertion <- CAL[CAL[,1]=="insertion",]

LIe10 <- read.table("Pindel_simulated_LI.vcf")

Pindel_ins_data <- vector()
for (i in 1:dim(LIe10)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- LIe10[i,2]
  Pindel_ins_info<- strsplit(as.character(LIe10[i,8]),";")[[1]]
  Pindel_ins_len <- Pindel_ins_info[grep("SVLEN=",Pindel_ins_info)]
  Pindel_ins_LEN <- substr(as.character(Pindel_ins_len),7,nchar(Pindel_ins_len)) 
  item <- c(position,Pindel_ins_LEN)
  Pindel_ins_data <- rbind(Pindel_ins_data,item)
}


# can not do, because don't know the length of input, no way to split by length, only can have a gernal precision and recall
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  for (i in 0:9)
  {
    truthset <- truthset_get[as.numeric(truthset_get[,3]) > i*200 & as.numeric(truthset_get[,3]) <= (i+1)*200,]
    position_in <- P_input_get[,1]
    position_truth <- truthset[,6]
    TP <- vector()
    #TP_sub <- intersect(position_in,position_truth)
    #position_truth <- position_truth[!(position_truth %in% TP_sub)]
    #TP <- c(TP,TP_sub)
    for (i in -20:20)
    {
      TP_sub <- intersect(position_in+i,position_truth)
      position_truth <- position_truth[!(position_truth %in% TP_sub)]
      TP <- c(TP,TP_sub)
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

parse_length(Pindel_ins_data,insertion)


TP_sub <- intersect(Pindel_ins_data[,1],insertion[insertion[,3]>200,6])















