VarScan <- read.table("varscan_indel.vcf")


VarScan_indel <- VarScan[nchar(as.character(VarScan[,4])) > 1 | nchar(as.character(VarScan[,5])) > 1,]
rownames(VarScan_indel) <- seq(length=nrow(VarScan_indel))# reset row name
VarScan_indel_del <- VarScan_indel[nchar(as.character(VarScan_indel[,4])) > 1,]
VarScan_indel_ins <- VarScan_indel[nchar(as.character(VarScan_indel[,5])) > 1,]
# make all indels together
VarScan_indel_ins1 <- vector()
VarScan_indel_del1 <- vector()
for (j in 1:dim(VarScan_indel_ins)[1])
{
  item <- vector()
  input_len <- nchar(as.character(VarScan_indel_ins[j,5]))
  item <- c(VarScan_indel_ins[j,2],input_len)
  VarScan_indel_ins1 <- rbind(VarScan_indel_ins1,item)
}
for (j in 1:dim(VarScan_indel_del)[1])
{
  item <- vector()
  input_len <- nchar(as.character(VarScan_indel_del[j,4]))
  item <- c(VarScan_indel_del[j,2],input_len)
  VarScan_indel_del1 <- rbind(VarScan_indel_del1,item)
}




# truthset
CAL <- read.table("CAL.txt")
insertion <- CAL[CAL[,1]=="insertion",]
deletion <- CAL[CAL[,1]=="deletion",]
duplication <- CAL[CAL[,1]=="duplication",]


# split dataset into 10 sbusets
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  all_TP <- vector()
  for (i in 0:9)
  {
    truthset <- truthset_get[as.numeric(truthset_get[,3]) > i*200 & as.numeric(truthset_get[,3]) <= (i+1)*200,]
    P_input <- P_input_get[as.numeric(as.character(P_input_get[,2])) > i*200 & as.numeric(as.character(P_input_get[,2])) <= (i+1)*200,]
    position_in <- as.numeric(as.character(P_input[,1]))
    position_truth <- as.numeric(truthset[,6])
    TP <- vector()
    #TP_sub <- intersect(position_in,position_truth)
    #position_truth <- position_truth[!(position_truth %in% TP_sub)]
    #TP <- c(TP,TP_sub)
    for (j in -5:5) # position tolerate
    {
      TP_sub <- intersect(as.numeric(position_in)+j,position_truth)
      for (P in TP_sub)
      {
        P_input_LEN <- as.numeric(as.character(P_input[P_input[,1]==P-j,2]))
        if (abs(as.numeric(truthset[truthset[,6]==P,3])-abs(as.numeric(P_input_LEN))) < 10) # length tolerate
        {
          # remove match position, avoid getting repeat results, due to one indel maybe report twice near by the break point
          position_truth <- position_truth[!(position_truth %in% P)] 
          TP <- c(TP,P)
          all_TP <- c(all_TP,P-j)
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
  return(all_TP)
}

#parse_length(VarScan_indel_del1,deletion)
#parse_length(VarScan_indel_ins1,insertion)
#parse_length(VarScan_indel_ins1,duplication)

parse_length(VarScan_indel_del1,deletion)
ins_correct <- parse_length(VarScan_indel_ins1,insertion)
dup_correct <-parse_length(VarScan_indel_ins1,duplication)

ins_post <- VarScan_indel_ins1[!(VarScan_indel_ins1[,1] %in% dup_correct),]
parse_length(ins_post,insertion)
dup_post <- VarScan_indel_ins1[!(VarScan_indel_ins1[,1] %in% ins_correct),]
parse_length(dup_post,duplication)













