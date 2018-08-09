# in Pindel, the short insertion is less than 82

# Get SI from simulated data
CAL <- read.table("CAL.txt")
insertion <- CAL[CAL[,1]=="insertion",]
SI_truth <- insertion[as.numeric(insertion[,3]) <= 200,]

SIe5 <- read.table("Pindel_simulated_SI.vcf")

Pindel_ins_data <- vector()
for (i in 1:dim(SIe5)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- SIe5[i,2]
  Pindel_ins_info<- strsplit(as.character(SIe5[i,8]),";")[[1]]
  Pindel_ins_len <- Pindel_ins_info[grep("SVLEN=",Pindel_ins_info)]
  Pindel_ins_LEN <- substr(as.character(Pindel_ins_len),7,nchar(Pindel_ins_len)) 
  item <- c(position,Pindel_ins_LEN)
  Pindel_ins_data <- rbind(Pindel_ins_data,item)
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

parse_length(Pindel_ins_data,insertion)


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
parse_length(Pindel_ins_data,insertion)














