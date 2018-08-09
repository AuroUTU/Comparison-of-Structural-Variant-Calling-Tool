ScanIndel_input <- read.table("simulated_30X.merged.indel.vcf")
position_del <- vector()
position_ins <- vector()
# according INFO column, split results into different SV_type 
for (i in 1:dim(ScanIndel_input)[1])
{
  ScanIndel_input_info <- strsplit(as.character(ScanIndel_input[i,8]),";")[[1]]
  ScanIndel_input_type <- ScanIndel_input_info[grep("TYPE=",ScanIndel_input_info)]
  ScanIndel_input_TYPE <- substr(as.character(ScanIndel_input_type),6,nchar(ScanIndel_input_type))
  if (ScanIndel_input_TYPE == "del")
  {
    position_del <- c(position_del,ScanIndel_input[i,2])
  }
  else if (ScanIndel_input_TYPE == "ins")
  {
    position_ins <- c(position_ins,ScanIndel_input[i,2])
  }
}
ScanIndel_del <- ScanIndel_input[ScanIndel_input[,2] %in% position_del,]
ScanIndel_ins <- ScanIndel_input[ScanIndel_input[,2] %in% position_ins,]

# truthset split into different SV_type
# column 1: SV_ytpe
# column 2: position in indel_genome.fa, the one remove all the gap
# column 3: the length of SV
# column 4: accelerate length change in that position(useless)
# column 5: new position after calcuate column 4 (stupid)
# column 6: simulated indels position in original genome, consider removed gaps(can be used to compare with VCF results)
CAL <- read.table("CAL.txt")
insertion <- CAL[CAL[,1]=="insertion",]
deletion <- CAL[CAL[,1]=="deletion",]
duplication <- CAL[CAL[,1]=="duplication",]

# make a function
# set a criteria:
# if position deviation is -5 to 5
# and indel length deviation is -10 to 10
parse_length <- function(P_input,truthset)
{
  position_in <- P_input[,2]
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
      # get length information from INFO column
      P_input_info <- strsplit(as.character(P_input[P_input[,2]==P-i,8]),";")[[1]]
      P_input_len <- P_input_info[grep("LEN=",P_input_info)]
      P_input_LEN <- substr(as.character(P_input_len),5,nchar(P_input_len)) 
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

parse_length(ScanIndel_del,deletion)
parse_length(ScanIndel_ins,insertion)
parse_length(ScanIndel_ins,duplication)


ScanIndel_del_data <- vector()
for (i in 1:dim(ScanIndel_del)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- ScanIndel_del[i,2]
  ScanIndel_del_info<- strsplit(as.character(ScanIndel_del[i,8]),";")[[1]]
  ScanIndel_del_len <- ScanIndel_del_info[grep("LEN=",ScanIndel_del_info)]
  ScanIndel_del_LEN <- substr(as.character(ScanIndel_del_len),5,nchar(ScanIndel_del_len)) 
  item <- c(position,ScanIndel_del_LEN)
  ScanIndel_del_data <- rbind(ScanIndel_del_data,item)
}
ScanIndel_del_data <- unique(ScanIndel_del_data)

ScanIndel_ins_data <- vector()
for (i in 1:dim(ScanIndel_ins)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- ScanIndel_ins[i,2]
  ScanIndel_ins_info<- strsplit(as.character(ScanIndel_ins[i,8]),";")[[1]]
  ScanIndel_ins_len <- ScanIndel_ins_info[grep("LEN=",ScanIndel_ins_info)]
  ScanIndel_ins_LEN <- substr(as.character(ScanIndel_ins_len),5,nchar(ScanIndel_ins_len)) 
  item <- c(position,ScanIndel_ins_LEN)
  ScanIndel_ins_data <- rbind(ScanIndel_ins_data,item)
}
ScanIndel_ins_data <- unique(ScanIndel_ins_data)

# split into length interval
# make a function
# set a criteria:
# if position deviation is -5 to 5
# and indel length deviation is -10 to 10
# split dataset into 10 sbusets
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  all_TP <- vector()
  for (i in 0:9)
  {
    truthset <- truthset_get[as.numeric(truthset_get[,3]) > i*200 & as.numeric(truthset_get[,3]) <= (i+1)*200,]
    P_input <- P_input_get[as.numeric(P_input_get[,2]) > i*200 & as.numeric(P_input_get[,2]) <= (i+1)*200,]
    #position_in <- unique(as.numeric(P_input[,1]))
    P_input <- matrix(unlist(P_input), ncol = 2) #as.numeric(as.character(P_input[,1]))
    position_in <- unique(as.numeric(as.character(P_input[,1])))
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
        for (thing in P_input_LEN)
        {
          if(is.na(thing)==FALSE)
          {
            P_input_LEN_real <- thing
            #print (P_input_LEN_real)
            break
          }
        }
        if (abs(as.numeric(truthset[truthset[,6]==P,3])-abs(as.numeric(P_input_LEN_real))) < 10) # length tolerate
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
    a <- length(sort(position_in))
    b <- length(sort(TP))
    print (c(b,a))
    print (a-b)  # don't know what reason, "position_in" always has a "ghost" NA data
    print ("")
    print("Precision")
    print(b/a)
    precision <- c(precision,b/a)
    print("Recall")
    print((length(TP)/dim(truthset)[1]))
    recall <- c(recall,length(TP)/dim(truthset)[1])
  }
  print (precision)
  print (recall)
  return(all_TP)
}

#parse_length(ScanIndel_del_data,deletion)
#parse_length(ScanIndel_ins_data,insertion)
#parse_length(ScanIndel_ins_data,duplication)

parse_length(ScanIndel_del_data,deletion)
ins_correct <- parse_length(ScanIndel_ins_data,insertion)
dup_correct <-parse_length(ScanIndel_ins_data,duplication)

ins_post <- ScanIndel_ins_data[!(ScanIndel_ins_data[,1] %in% dup_correct),]
parse_length(ins_post,insertion)
dup_post <- ScanIndel_ins_data[!(ScanIndel_ins_data[,1] %in% ins_correct),]
parse_length(dup_post,duplication)

# Only parse 50bp--150bp
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  all_TP <- vector()
  truthset <- truthset_get[as.numeric(truthset_get[,3]) > 50 & as.numeric(truthset_get[,3]) <= 150,]
  P_input <- P_input_get[as.numeric(P_input_get[,2]) > 50 & as.numeric(P_input_get[,2]) <= 150,]
  P_input <- matrix(unlist(P_input), ncol = 2) #as.numeric(as.character(P_input[,1]))
  position_in <- unique(as.numeric(as.character(P_input[,1])))
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
      for (thing in P_input_LEN)
      {
        if(is.na(thing)==FALSE)
        {
          P_input_LEN_real <- thing
          # print (P_input_LEN_real)
          break
        }
      }
      if (abs(as.numeric(truthset[truthset[,6]==P,3])-abs(as.numeric(P_input_LEN_real))) < 10) # length tolerate
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
#parse_length(ScanIndel_del_data,deletion)
#parse_length(ScanIndel_ins_data,insertion)
#parse_length(ScanIndel_ins_data,duplication)

parse_length(ScanIndel_del_data,deletion)
ins_correct <- parse_length(ScanIndel_ins_data,insertion)
dup_correct <-parse_length(ScanIndel_ins_data,duplication)

ins_post <- ScanIndel_ins_data[!(ScanIndel_ins_data[,1] %in% dup_correct),]
parse_length(ins_post,insertion)
dup_post <- ScanIndel_ins_data[!(ScanIndel_ins_data[,1] %in% ins_correct),]
parse_length(dup_post,duplication)





























