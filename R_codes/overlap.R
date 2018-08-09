# return TP 
parse_length <- function(P_input,truthset)
{
  position_in <- as.numeric(as.character(P_input[,1]))
  position_truth <- truthset[,6]
  TP <- vector()
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
  return(TP)
}

CAL <- read.table("CAL.txt")
insertion <- CAL[CAL[,1]=="insertion",]
deletion <- CAL[CAL[,1]=="deletion",]
duplication <- CAL[CAL[,1]=="duplication",]

# Delly
setwd("/home/kning/Bio/Indels_calling/mfranberg_simulated/Results_30/delly")
Delly_input <- read.table("delly_default.vcf")
position_del <- vector()
position_ins <- vector()
position_dup <- vector()
# type <- vector()
for (i in 1:dim(Delly_input)[1])
{
  Delly_input_info <- strsplit(as.character(Delly_input[i,8]),";")[[1]]
  Delly_input_type <- Delly_input_info[grep("SVTYPE=",Delly_input_info)]
  Delly_input_TYPE <- substr(as.character(Delly_input_type),8,nchar(Delly_input_type))
  # type <- c(type,Delly_input_TYPE) contain "DEL" "DUP" "INS" "INV"
  if (Delly_input_TYPE == "DEL")
  {
    position_del <- c(position_del,Delly_input[i,2])
  }
  else if (Delly_input_TYPE == "INS")
  {
    position_ins <- c(position_ins,Delly_input[i,2])
  }
  else if (Delly_input_TYPE == "DUP")
  {
    position_dup <- c(position_dup,Delly_input[i,2])
  }
}
Delly_del <- Delly_input[Delly_input[,2] %in% position_del,]
Delly_ins <- Delly_input[Delly_input[,2] %in% position_ins,]
Delly_dup <- Delly_input[Delly_input[,2] %in% position_dup,]
Delly_del_data <- vector()
for (i in 1:dim(Delly_del)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- Delly_del[i,2]
  Delly_del_SV_start <- Delly_del[i,2]
  Delly_del_info <- strsplit(as.character(Delly_del[i,8]),";")[[1]]
  Delly_del_input_end <- Delly_del_info[grep("END=",Delly_del_info)][1]
  Delly_del_SV_END <- as.numeric(substr(as.character(Delly_del_input_end),5,nchar(Delly_del_input_end)))
  Delly_del_input_LEN <- Delly_del_SV_END-Delly_del_SV_start
  item <- c(position,Delly_del_input_LEN)
  Delly_del_data <- rbind(Delly_del_data,item)
}

# Delly insertion have INSLEN in VCF format column
# While deletion and duplication INSLEN are all 0, length can only calculate by END-position
Delly_ins_data <- vector()
for (i in 1:dim(Delly_ins)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- Delly_ins[i,2]
  Delly_ins_SV_start <- Delly_ins[i,2]
  Delly_ins_info <- strsplit(as.character(Delly_ins[i,8]),";")[[1]]
  Delly_ins_input_len <- Delly_ins_info[grep("INSLEN=",Delly_ins_info)][1]
  Delly_ins_input_LEN <- as.numeric(substr(as.character(Delly_ins_input_len),8,nchar(Delly_ins_input_len)))
  item <- c(position,Delly_ins_input_LEN)
  Delly_ins_data <- rbind(Delly_ins_data,item)
}

Delly_dup_data <- vector()
for (i in 1:dim(Delly_dup)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- Delly_dup[i,2]
  Delly_dup_SV_start <- Delly_dup[i,2]
  Delly_dup_info <- strsplit(as.character(Delly_dup[i,8]),";")[[1]]
  Delly_dup_input_end <- Delly_dup_info[grep("END=",Delly_dup_info)][1]
  Delly_dup_SV_END <- as.numeric(substr(as.character(Delly_dup_input_end),5,nchar(Delly_dup_input_end)))
  Delly_dup_input_LEN <- Delly_dup_SV_END-Delly_dup_SV_start
  item <- c(position,Delly_dup_input_LEN)
  Delly_dup_data <- rbind(Delly_dup_data,item)
}
Delly_del <- parse_length(Delly_del_data,deletion)
Delly_ins <- parse_length(Delly_ins_data,insertion)
Delly_dup <- parse_length(Delly_dup_data,duplication)
Delly_TP <- c(Delly_del,Delly_ins,Delly_dup)

# fermikit
setwd("/home/kning/Bio/Indels_calling/mfranberg_simulated/Results_30/fermikit")
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
  input_len <- nchar(as.character(fermikit_flt_ins[j,5]))-1
  item <- c(fermikit_flt_ins[j,2],input_len)
  fermikit_flt_ins1 <- rbind(fermikit_flt_ins1,item)
  fermikit_sv_ins <- rbind(fermikit_sv_ins,item)  
}
for (j in 1:dim(fermikit_flt_del)[1])
{
  item <- vector()
  input_len <- nchar(as.character(fermikit_flt_del[j,4]))-1
  item <- c(fermikit_flt_del[j,2],input_len)
  fermikit_flt_del1 <- rbind(fermikit_flt_del1,item)
  fermikit_sv_del <- rbind(fermikit_sv_del,item)  
}

# make matrix into dataframe and sort it by position
rownames(fermikit_sv_del) <- seq(length=nrow(fermikit_sv_del))
fermikit_sv_del <- data.frame(fermikit_sv_del)
fermikit_sv_del <- fermikit_sv_del[order(fermikit_sv_del$X1),]
rownames(fermikit_sv_del) <- seq(length=nrow(fermikit_sv_del))

rownames(fermikit_sv_ins) <- seq(length=nrow(fermikit_sv_ins))
fermikit_sv_ins <- data.frame(fermikit_sv_ins)
fermikit_sv_ins <- fermikit_sv_ins[order(fermikit_sv_ins$X1),]
rownames(fermikit_sv_ins) <- seq(length=nrow(fermikit_sv_ins))

fermikit_sv_del <- as.matrix(fermikit_sv_del)
fermikit_sv_ins <- as.matrix(fermikit_sv_ins)

fermikit_del <- parse_length(fermikit_sv_del,deletion)
fermikit_ins <- parse_length(fermikit_sv_ins,insertion)
fermikit_dup <- parse_length(fermikit_sv_ins,duplication)
fermikit_TP <- c(fermikit_del,fermikit_ins,fermikit_dup)

# freebayes
setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/freebayes")
freebayes_input <- read.table("freebayes.indel.vcf")
# only guess type from REF and ALT
freebayes_ins1 <- vector()
freebayes_del1 <- vector()
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
    freebayes_ins1 <- rbind(freebayes_ins1,item)  # first column is position, second column is SV length
  }
  else if (input_type == "TYPE=del")
  {
    item <- c(freebayes_input[i,2],input_LEN)
    freebayes_del1 <- rbind(freebayes_del1,item)  # first column is position, second column is SV length
  }
}
freebayes_del <- parse_length(freebayes_del1,deletion)
freebayes_ins <- parse_length(freebayes_ins1,insertion)
freebayes_dup <- parse_length(freebayes_ins1,duplication)
freebayes_TP <- c(freebayes_del,freebayes_ins,freebayes_dup)

# lumpy
setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/LUMPY")
lumpy_input <- read.table("sample.vcf")
lumpy_del <- vector()
lumpy_dup <- vector()
lumpy_bnd <- vector()
for (i in 1:dim(lumpy_input)[1])
{
  lumpy_input_info <- strsplit(as.character(lumpy_input[i,8]),";")[[1]]
  lumpy_input_type <- lumpy_input_info[grep("SVTYPE=",lumpy_input_info)]
  lumpy_input_TYPE <- substr(as.character(lumpy_input_type),8,nchar(lumpy_input_type))
  lumpy_input_len <- lumpy_input_info[grep("SVLEN=",lumpy_input_info)]
  lumpy_input_LEN <- substr(as.character(lumpy_input_len),7,nchar(lumpy_input_len)) 
  if (lumpy_input_TYPE == "DEL")
  {
    item <- vector()
    item <- c(lumpy_input[i,2],abs(as.numeric(lumpy_input_LEN)))
    lumpy_del <- rbind(lumpy_del,item)
  }
  else if (lumpy_input_TYPE == "DUP")
  {
    item <- vector()
    item <- c(lumpy_input[i,2],as.numeric(lumpy_input_LEN))
    lumpy_dup <- rbind(lumpy_dup,item)
  }
  else if (lumpy_input_TYPE == "BND")
  {
    item <- vector()
    item <- c(lumpy_input[i,2],as.numeric(lumpy_input_LEN))
    lumpy_bnd <- rbind(lumpy_bnd,item)
  }
}
lumpy_del <- parse_length(lumpy_del,deletion)
lumpy_dup <- parse_length(lumpy_dup,duplication)
lumpy_TP <- c(lumpy_del,lumpy_dup)

# Platypus
setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/Platypus")
Platypus_input <- read.table("variants.vcf")
Platypus_ins <- Platypus_input[nchar(as.character(Platypus_input[,4])) < nchar(as.character(Platypus_input[,5])),]
Platypus_del <- Platypus_input[nchar(as.character(Platypus_input[,4])) > nchar(as.character(Platypus_input[,5])),]
# get some SNV, which can be false positive
Platypus_ins_data <- vector()
for (i in 1:dim(Platypus_ins)[1])
{
  item <- vector()
  position <- as.numeric(Platypus_ins[i,2])#+nchar(as.character(Platypus_ins[i,5]))
  length <-  abs(nchar(as.character(Platypus_ins[i,4])) - nchar(as.character(Platypus_ins[i,5])))
  item <- c(position,length)
  Platypus_ins_data <- rbind(Platypus_ins_data,item)
}
Platypus_del_data <- vector()
for (i in 1:dim(Platypus_del)[1])
{
  item <- vector()
  position <- as.numeric(Platypus_del[i,2])#+nchar(as.character(Platypus_del[i,5]))
  length <-  abs(nchar(as.character(Platypus_del[i,4])) - nchar(as.character(Platypus_del[i,5])))
  item <- c(position,length)
  Platypus_del_data <- rbind(Platypus_del_data,item)
}
Platypus_del <- parse_length(Platypus_del_data,deletion)
Platypus_ins <- parse_length(Platypus_ins_data,insertion)
Platypus_dup <- parse_length(Platypus_ins_data,duplication)
Platypus_TP <- c(Platypus_del,Platypus_ins,Platypus_dup)

# ScanIndel
setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/ScanIndel")
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
ScanIndel_del <- parse_length(ScanIndel_del_data,deletion)
ScanIndel_ins <- parse_length(ScanIndel_ins_data,insertion)
ScanIndel_dup <- parse_length(ScanIndel_ins_data,duplication)
ScanIndel_TP <- c(ScanIndel_del,ScanIndel_ins,ScanIndel_dup)

# VarDict
setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/VarDict")
VarDict_input <- read.table("vardict.simulated.vcf")
position_del <- vector()
position_ins <- vector()
position_SNV <- vector()
for (i in 1:dim(VarDict_input)[1])
{
  VarDict_input_info <- strsplit(as.character(VarDict_input[i,8]),";")[[1]]
  VarDict_input_type <- VarDict_input_info[grep("TYPE=",VarDict_input_info)]
  VarDict_input_TYPE <- substr(as.character(VarDict_input_type),6,nchar(VarDict_input_type))
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

VarDict_del <- parse_length(VarDict_del_data,deletion)
VarDict_ins <- parse_length(VarDict_ins_data,insertion)
VarDict_dup <- parse_length(VarDict_ins_data,duplication)
VarDict_TP <- c(VarDict_del,VarDict_ins,VarDict_dup)

# VarScan
setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/VarScan")
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
VarScan_del <- parse_length(VarScan_indel_del1,deletion)
VarScan_ins <- parse_length(VarScan_indel_ins1,insertion)
VarScan_dup <- parse_length(VarScan_indel_ins1,duplication)
VarScan_TP <- c(VarScan_del,VarScan_ins,VarScan_dup)

#GATK
setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/GATK")
GATK <- read.table("qual_indels.vcf")
#GATK2 <- read.table("filtered_indels.vcf")
#GATK <- GATK[as.numeric(GATK[,6]) > 280 ,] GATK hard filter apply this
GATK <- GATK[GATK[,7]=="qual_test" ,]
GATK_indel <- GATK[nchar(as.character(GATK[,4])) > 1 | nchar(as.character(GATK[,5])) > 1,]
rownames(GATK_indel) <- seq(length=nrow(GATK_indel))# reset row name
GATK_indel_del <- GATK_indel[nchar(as.character(GATK_indel[,4])) > 1,]
GATK_indel_ins <- GATK_indel[nchar(as.character(GATK_indel[,5])) > 1,]
# make all indels together
GATK_indel_ins1 <- vector()
GATK_indel_del1 <- vector()
for (j in 1:dim(GATK_indel_ins)[1])
{
  item <- vector()
  input_len <- nchar(as.character(GATK_indel_ins[j,5]))
  item <- c(GATK_indel_ins[j,2],input_len)
  GATK_indel_ins1 <- rbind(GATK_indel_ins1,item)
}
for (j in 1:dim(GATK_indel_del)[1])
{
  item <- vector()
  input_len <- nchar(as.character(GATK_indel_del[j,4]))
  item <- c(GATK_indel_del[j,2],input_len)
  GATK_indel_del1 <- rbind(GATK_indel_del1,item)
}
GATK_del <- parse_length(GATK_indel_del1,deletion)
GATK_ins <- parse_length(GATK_indel_ins1,insertion)
GATK_dup <- parse_length(GATK_indel_ins1,duplication)
GATK_TP <- c(GATK_del,GATK_ins,GATK_dup)

# Pindel
setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/Pindel")
D_truth <- CAL[CAL[,1]=="deletion",]
De10 <- read.table("Pindel_simulated_D.vcf")
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
Pindel_del <- parse_length(Pindel_del_data,D_truth)

insertion <- CAL[CAL[,1]=="insertion",]
SI_truth <- insertion[as.numeric(insertion[,3]) <= 200,]
SIe5 <- read.table("Pindel_simulated_SI.vcf")
Pindel_ins_data_S <- vector()
for (i in 1:dim(SIe5)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- SIe5[i,2]
  Pindel_ins_info<- strsplit(as.character(SIe5[i,8]),";")[[1]]
  Pindel_ins_len <- Pindel_ins_info[grep("SVLEN=",Pindel_ins_info)]
  Pindel_ins_LEN <- substr(as.character(Pindel_ins_len),7,nchar(Pindel_ins_len)) 
  item <- c(position,Pindel_ins_LEN)
  Pindel_ins_data_S <- rbind(Pindel_ins_data_S,item)
}
Pindel_ins_S <- parse_length(Pindel_ins_data_S,insertion)

LIe10 <- read.table("Pindel_simulated_LI.vcf")
Pindel_ins_data_L <- vector()
for (i in 1:dim(LIe10)[1])
{
  item <- vector()
  # get length information from INFO column
  position <- LIe10[i,2]
  Pindel_ins_info<- strsplit(as.character(LIe10[i,8]),";")[[1]]
  Pindel_ins_len <- Pindel_ins_info[grep("SVLEN=",Pindel_ins_info)]
  Pindel_ins_LEN <- substr(as.character(Pindel_ins_len),7,nchar(Pindel_ins_len)) 
  item <- c(position,Pindel_ins_LEN)
  Pindel_ins_data_L <- rbind(Pindel_ins_data_L,item)
}
Pindel_ins_L <- intersect(Pindel_ins_data_L[,1],insertion[insertion[,3]>200,6])

CAL <- read.table("CAL.txt")
TD_truth <- CAL[CAL[,1]=="duplication",]
TDe10 <- read.table("Pindel_simulated_TD.vcf")
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
Pindel_dup <- parse_length(Pindel_dup_data,TD_truth)
Pindel_ins <- c(Pindel_ins_S,Pindel_ins_L)
Pindel_TP <- c(Pindel_del,Pindel_ins_S,Pindel_ins_L,Pindel_dup)

length(intersect(as.numeric(Delly_TP),as.numeric(lumpy_TP)))

setwd("~/Bio/Indels_calling/mfranberg_simulated/Results_30/Similarity")

similarity <- read.table("TD_test.txt")
total_TP <- read.table("Total_DUP.txt")
for (i in 1:9)
{
  for (j in 1:9)
  {
    similarity[i,j] <- (similarity[i,j]*2)/(total_TP[i+1,2]+ total_TP[j,2])
  }
}
similarity <- as.matrix(similarity)
#for (i in 1:9)
#{
#  for (j in 1:9)
#  {
#    if (i!=j)
#    {
#      similarity[i,j]=similarity[i,j]-0.5
#    }
#  }
#}
library(corrplot)
#corrplot(similarity,method = "circle",col = c("black", "white"), bg = "lightblue",
#        type = "upper",p.mat = similarity,insig = "p-value",sig.level = -1)

corrplot(similarity,method = "circle",col = colorRampPalette(brewer.pal(n=11, name="Spectral"))(100), bg = "white",
         type = "lower",tl.col="black",cl.lim = c(0,1),mar=c(0,0,2,0),main="DUP Similarity")


corrplot(similarity, is.corr = FALSE,
         cl.lim=c(0,100),
         p.mat = perc_matrix,
         insig = "p-value",
         sig.level = -1,
         tl.col="black",
         col = colorRampPalette(brewer.pal(n=11, name="Spectral"))(100),
         mar=c(0,0,1,0))





