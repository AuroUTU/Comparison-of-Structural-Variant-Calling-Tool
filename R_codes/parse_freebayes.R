freebayes_input <- read.table("freebayes.indel.vcf")

CAL <- read.table("CAL.txt")
insertion <- CAL[CAL[,1]=="insertion",]
deletion <- CAL[CAL[,1]=="deletion",]
duplication <- CAL[CAL[,1]=="duplication",]

SD_truth <- deletion[as.numeric(deletion[,3]) <= 40,]
SI_truth <- insertion[as.numeric(insertion[,3]) <= 40,]
SDU_truth <- duplication[as.numeric(insertion[,3]) <= 40,]
Short_indel <- rbind(SD_truth,SI_truth)
Short_indel <- Short_indel[order(Short_indel$V2),]

found <- intersect(freebayes_input[,2]+2,Short_indel[,6])

# only guess type from REF and ALT
freebayes_ins <- vector()
freebayes_del <- vector()
for (i in 1: dim(freebayes_input)[1])
{
  item <- vector()
  input_info <- strsplit(as.character(freebayes_input[i,8]),";")[[1]]
  input_type <- input_info[grep("TYPE=",input_info)]
  input_len <- input_info[grep("LEN=",input_info)]
  input_LEN <- substr(as.character(input_len),5,nchar(input_len))
  if (input_type == "TYPE=ins")
  {
    item <- c(freebayes_input[i,2],input_LEN)
    freebayes_ins <- rbind(freebayes_ins,item)  # first column is position, second column is SV length
  }
  else if (input_type == "TYPE=del")
  {
    item <- c(freebayes_input[i,2],input_LEN)
    freebayes_del <- rbind(freebayes_del,item)  # first column is position, second column is SV length
  }
}

#make a function
parse_length <- function(P_input,truthset)
{
  position_in <- as.numeric(as.character(P_input[,1]))
  position_truth <- truthset[,6]
  
  TP <- vector()
  #TP_sub <- intersect(position_in,position_truth)
  #position_truth <- position_truth[!(position_truth %in% TP_sub)]
  #TP <- c(TP,TP_sub)
  for (i in -5:5) # position tolerate
  {
    TP_sub <- intersect(position_in+i,position_truth)
    for (P in TP_sub)
    {
      P_input_LEN <- as.character(P_input[P_input[,1]==P-i,2])
      if (abs(as.numeric(truthset[truthset[,6]==P,3])-abs(as.numeric(P_input_LEN))) < 10) # length tolerate
      {
        # remove match position, avoid getting repeat results, due to one indel maybe report twice near by the break point
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

parse_length(freebayes_del,SD_truth)
parse_length(freebayes_ins,SI_truth)
parse_length(freebayes_ins,SDU_truth)






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

#parse_length(freebayes_del,deletion)
#parse_length(freebayes_ins,insertion)
#parse_length(freebayes_ins,duplication)

parse_length(freebayes_del,deletion)
ins_correct <- parse_length(freebayes_ins,insertion)
dup_correct <-parse_length(freebayes_ins,duplication)

ins_post <- freebayes_ins[!(freebayes_ins[,1] %in% dup_correct),]
parse_length(ins_post,insertion)
dup_post <- freebayes_ins[!(freebayes_ins[,1] %in% ins_correct),]
parse_length(dup_post,duplication)




