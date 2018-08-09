fermikit_sv <- read.table("fermikit.sv.vcf")
# FermiKit only has 3 SVTYPE: INS,DEL,COMPLEX
fermikit_sv_ins <- vector()
fermikit_sv_del <- vector()
fermikit_sv_cpx <- vector()
# split indels from SV result
# summary(as.numeric(fermikit_sv_ins[,2])), min_length is 104
# summary(as.numeric(fermikit_sv_del[,2])), min_length is 88
for (i in 1: dim(fermikit_sv)[1])
{
  item <- vector()
  input_info <- strsplit(as.character(fermikit_sv[i,8]),";")[[1]]
  input_type <- input_info[grep("SVTYPE=",input_info)]
  input_len <- input_info[grep("SVLEN=",input_info)]
  input_LEN <- substr(as.character(input_len),7,nchar(input_len))
  if (input_type == "SVTYPE=INS")
  {
    item <- c(fermikit_sv[i,2],input_LEN)
    fermikit_sv_ins <- rbind(fermikit_sv_ins,item)  # first column is position, second column is SV length
  }
  else if (input_type == "SVTYPE=DEL")
  {
    item <- c(fermikit_sv[i,2],input_LEN)
    fermikit_sv_del <- rbind(fermikit_sv_del,item)  # first column is position, second column is SV length
  }
  else if (input_type == "SVTYPE=COMPLEX")
  {
    item <- c(fermikit_sv[i,2],input_LEN)
    fermikit_sv_cpx <- rbind(fermikit_sv_cpx,item)  # first column is position, second column is SV length
  }
}

# get insertion and deletion from flt result
fermikit_flt <- read.table("fermikit.flt.vcf")
fermikit_flt_indel <- fermikit_flt[nchar(as.character(fermikit_flt[,4])) > 1 | nchar(as.character(fermikit_flt[,5])) > 1,]
rownames(fermikit_flt_indel) <- seq(length=nrow(fermikit_flt_indel))# reset row name
fermikit_flt_ins <- fermikit_flt_indel[nchar(as.character(fermikit_flt_indel[,5])) > 1,]
fermikit_flt_del <- fermikit_flt_indel[nchar(as.character(fermikit_flt_indel[,4])) > 1,]
# make all indels together
fermikit_flt_ins1 <- vector()
fermikit_flt_del1 <- vector()
for (j in 1:dim(fermikit_flt_ins)[1])
{
  item <- vector()
  input_len <- nchar(as.character(fermikit_flt_ins[j,5]))-1 # insertion length is column 5 
  item <- c(fermikit_flt_ins[j,2],input_len)
  fermikit_flt_ins1 <- rbind(fermikit_flt_ins1,item)
  fermikit_sv_ins <- rbind(fermikit_sv_ins,item)  
}
for (j in 1:dim(fermikit_flt_del)[1])
{
  item <- vector()
  input_len <- nchar(as.character(fermikit_flt_del[j,4]))-1 # deletion length is column 4 
  item <- c(fermikit_flt_del[j,2],input_len)
  fermikit_flt_del1 <- rbind(fermikit_flt_del1,item)
  fermikit_sv_del <- rbind(fermikit_sv_del,item)  
}

# summary(fermikit_flt_del1[,2]) max_length:197
# summary(fermikit_flt_ins1[,2]) max_length:200

# make matrix into dataframe and sort it by position
rownames(fermikit_sv_del) <- seq(length=nrow(fermikit_sv_del))
fermikit_sv_del <- data.frame(fermikit_sv_del)
fermikit_sv_del <- fermikit_sv_del[order(fermikit_sv_del$X1),]
rownames(fermikit_sv_del) <- seq(length=nrow(fermikit_sv_del))

rownames(fermikit_sv_ins) <- seq(length=nrow(fermikit_sv_ins))
fermikit_sv_ins <- data.frame(fermikit_sv_ins)
fermikit_sv_ins <- fermikit_sv_ins[order(fermikit_sv_ins$X1),]
rownames(fermikit_sv_ins) <- seq(length=nrow(fermikit_sv_ins))

# truthset
CAL <- read.table("CAL.txt")
insertion <- CAL[CAL[,1]=="insertion",]
deletion <- CAL[CAL[,1]=="deletion",]
duplication <- CAL[CAL[,1]=="duplication",]

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

parse_length(fermikit_sv_del,deletion)
parse_length(fermikit_sv_ins,insertion)
parse_length(fermikit_sv_ins,duplication)
parse_length(fermikit_sv_cpx,duplication) # complex has no SVLEN 


# split dataset into 10 sbusets
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  all_TP <- vector() # get all the TP insertion or duplication, and remove them from other category, otherwise, the TP will treat as FP 
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
  return (all_TP)
}

parse_length(fermikit_sv_del,deletion)
ins_correct <- parse_length(fermikit_sv_ins,insertion)
dup_correct <-parse_length(fermikit_sv_ins,duplication)

ins_post <- fermikit_sv_ins[!(fermikit_sv_ins[,1] %in% dup_correct),]
parse_length(ins_post,insertion)
dup_post <- fermikit_sv_ins[!(fermikit_sv_ins[,1] %in% ins_correct),]
parse_length(dup_post,duplication)

# Only parse 50bp--150bp
parse_length <- function(P_input_get,truthset_get)
{
  precision <- vector()
  recall <- vector()
  all_TP <- vector()
  truthset <- truthset_get[as.numeric(truthset_get[,3]) > 50 & as.numeric(truthset_get[,3]) <= 150,]
  P_input <- P_input_get[as.numeric(as.character(P_input_get[,2])) > 50 & as.numeric(as.character(P_input_get[,2])) <= 150,]
  position_in <- as.numeric(as.character(P_input[,1]))
  position_truth <- as.numeric(truthset[,6])
  TP <- vector()
  for (j in -5:5) # position tolerate
  {
    TP_sub <- intersect(as.numeric(position_in)+j,position_truth)
    #print (TP_sub)
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
  print (c(50,150))
  #print ("Total TP:")
  #print (length(TP)) # how many TP
  #print ("Total FP:")
  #print (length(position_in)-length(TP))
  #print ("")
  #print("Precision")
  #print(length(TP)/length(position_in))
  precision <- c(precision,length(TP)/length(position_in))
  #print("Recall")
  #print((length(TP)/dim(truthset)[1]))
  recall <- c(recall,length(TP)/dim(truthset)[1])
  print (precision)
  print (recall)
  return (all_TP)
}
#parse_length(fermikit_sv_del,deletion)
#parse_length(fermikit_sv_ins,insertion)
#parse_length(fermikit_sv_ins,duplication)

parse_length(fermikit_sv_del,deletion)
ins_correct <- parse_length(fermikit_sv_ins,insertion)
dup_correct <-parse_length(fermikit_sv_ins,duplication)

ins_post <- fermikit_sv_ins[!(fermikit_sv_ins[,1] %in% dup_correct),]
parse_length(ins_post,insertion)
dup_post <- fermikit_sv_ins[!(fermikit_sv_ins[,1] %in% ins_correct),]
parse_length(dup_post,duplication)
