VarDict_input <- read.table("vardict.simulated.vcf")
position_del <- vector()
position_ins <- vector()
position_SNV <- vector()
#type <- vector()
for (i in 1:dim(VarDict_input)[1])
{
  VarDict_input_info <- strsplit(as.character(VarDict_input[i,8]),";")[[1]]
  VarDict_input_type <- VarDict_input_info[grep("TYPE=",VarDict_input_info)]
  VarDict_input_TYPE <- substr(as.character(VarDict_input_type),6,nchar(VarDict_input_type))
  #type <- c(type,VarDict_input_TYPE) "Deletion"  "Insertion" "SNV"       "Complex"
  if (VarDict_input_TYPE == "Deletion")
  {
    position_del <- c(position_del,VarDict_input[i,2])
  }
  else if (VarDict_input_TYPE == "Insertion")
  {
    position_ins <- c(position_ins,VarDict_input[i,2])
  }
  else if (VarDict_input_TYPE == "SNV")
  {
    position_SNV <- c(position_SNV,VarDict_input[i,2])
  }
}

VarDict_del <- VarDict_input[VarDict_input[,2] %in% position_del,]
VarDict_ins <- VarDict_input[VarDict_input[,2] %in% position_ins,]

CAL <- read.table("CAL.txt")
insertion <- CAL[CAL[,1]=="insertion",]
deletion <- CAL[CAL[,1]=="deletion",]
duplication <- CAL[CAL[,1]=="duplication",]

# make a function
parse_length <- function(P_input,truthset,Indel_type)
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
      P_input_LEN <- nchar(as.character(P_input[P_input[,2]==P-i,as.numeric(Indel_type)]))
      measure <- abs(as.numeric(truthset[truthset[,6]==P,3])-as.numeric(P_input_LEN))
      if (measure < 10)
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

parse_length(VarDict_del,deletion,"4")
parse_length(VarDict_ins,insertion,"5")
parse_length(VarDict_ins,duplication,"5")


VarDict_del_data <- vector()
for (i in 1:dim(VarDict_del)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- VarDict_del[i,2]
  VarDict_del_LEN <- nchar(as.character(VarDict_del[i,4])) 
  item <- c(position,VarDict_del_LEN)
  VarDict_del_data <- rbind(VarDict_del_data,item)
}

VarDict_ins_data <- vector()
for (i in 1:dim(VarDict_ins)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- VarDict_ins[i,2]
  VarDict_ins_LEN <- nchar(as.character(VarDict_ins[i,5])) 
  item <- c(position,VarDict_ins_LEN)
  VarDict_ins_data <- rbind(VarDict_ins_data,item)
}


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
    P_input <- matrix(unlist(P_input), ncol = 2) #as.numeric(as.character(P_input[,1]))
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

#parse_length(VarDict_del_data,deletion)
#parse_length(VarDict_ins_data,insertion)
#parse_length(VarDict_ins_data,duplication)

parse_length(VarDict_del_data,deletion)
ins_correct <- parse_length(VarDict_ins_data,insertion)
dup_correct <-parse_length(VarDict_ins_data,duplication)

ins_post <- VarDict_ins_data[!(VarDict_ins_data[,1] %in% dup_correct),]
parse_length(ins_post,insertion)
dup_post <- VarDict_ins_data[!(VarDict_ins_data[,1] %in% ins_correct),]
parse_length(dup_post,duplication)

# Only parse 50bp--150bp
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  all_TP <- vector()
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
        all_TP <- c(all_TP,P-j)
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
  return(all_TP)
}
#parse_length(VarDict_del_data,deletion)
#parse_length(VarDict_ins_data,insertion)
#parse_length(VarDict_ins_data,duplication)

parse_length(VarDict_del_data,deletion)
ins_correct <- parse_length(VarDict_ins_data,insertion)
dup_correct <-parse_length(VarDict_ins_data,duplication)

ins_post <- VarDict_ins_data[!(VarDict_ins_data[,1] %in% dup_correct),]
parse_length(ins_post,insertion)
dup_post <- VarDict_ins_data[!(VarDict_ins_data[,1] %in% ins_correct),]
parse_length(dup_post,duplication)


























