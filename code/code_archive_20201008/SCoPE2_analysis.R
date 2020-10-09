source("code/functions_parameters.R")

# Import ------------------------------------------------------------------

# Load raw data and annotations

  # ev<-read_csv("dat/ev_tmt11_tmt16_2.csv")
  #
  # parse_row<-grep("|",ev$Leading.razor.protein, fixed=T)
  #
  # split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
  # split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  #
  # ev$Leading.razor.protein[parse_row]<-split_prot2
  #
  # #grep("|",ev$Leading.razor.protein, fixed=T)
  #
  # design<-read.csv("dat/annotation.csv")
  #
  # batch<-read.csv("dat/batch.csv")
  #
  # # Attach batch data to protein data
  # ev[,colnames(batch)[-1]]<-NA
  # for(X in batch$set){
  #
  #   ev$lcbatch[ev$Raw.file==X] <- as.character(batch$lcbatch[batch$set%in%X])
  #   ev$sortday[ev$Raw.file==X] <- as.character(batch$sortday[batch$set%in%X])
  #   ev$digest[ev$Raw.file==X] <- as.character(batch$digest[batch$set%in%X])
  #
  # }
  #
  #
  #  save(ev,batch,design,file="dat/raw.RData")

load("dat/raw.RData")

write.csv(paste0(unique(ev$Raw.file),".raw | "), "all_raw_files_used_besides_ladder.csv")

xx<-unique(ev$Raw.file)

length(grep("SQC", xx))
length(grep("SQC", xx))


# Create unique peptide+charge column:
ev$modseq<-paste0(ev$Modified.sequence,ev$Charge)

# If somehow we have NA for raw files?! Remove 'em
ev<-ev[!is.na(ev$Raw.file),]

# Add X in front of experiment names because R doesn't like column names starting with numbers
ev$Raw.file<-paste0("X",ev$Raw.file)
design$Set<-paste0("X",design$Set)

# Group the experimental runs by type: single cell, 10 cell, 100 cell, 1000 cell, or ladder
all.runs<-unique(ev$Raw.file); length(all.runs)
c10.runs<-c( all.runs[grep("col19", all.runs)] , all.runs[grep("col20", all.runs)] ); length(unique(c10.runs))
c100.runs<-c( all.runs[grep("col21", all.runs)] , all.runs[grep("col22", all.runs)] ); length(unique(c100.runs))
c1000.runs<-c( all.runs[grep("col23", all.runs)] ); length(unique(c1000.runs))
ladder.runs<-c( all.runs[grep("col24", all.runs)] ); length(unique(ladder.runs))
carrier.runs<-c( all.runs[grep("arrier", all.runs)] ); length(unique(carrier.runs))
other.runs<-c( all.runs[c(grep("Ref", all.runs), grep("Master", all.runs), grep("SQC", all.runs), grep("blank", all.runs) )] ); length(unique(other.runs))
sc.runs<-all.runs[!all.runs%in%c(c10.runs,c1000.runs,c100.runs,ladder.runs, carrier.runs, other.runs)]; length(unique(sc.runs))

# To demonstrate raison d'etre for reference channel, only use unbalanced single cell sets and normalize to mean:
if(ref_demo){

  m0_u_ratio<-c()
  rawt<-c()

  design2<-design[design$Set%in%sc.runs, ]

  for(i in 1:nrow(design2)){

    xt<-design2$Set[i]
    m0count<-length(which(design2[i,-1]=="sc_m0"))
    ucount<-length(which(design2[i,-1]=="sc_u"))

    rawt<-c(rawt, xt)
    m0_u_ratio<-c(m0_u_ratio, m0count / ucount)

  }

  unbalanced_sets_m0<-rawt[m0_u_ratio < quantile(m0_u_ratio, prob=0.1) | m0_u_ratio > quantile(m0_u_ratio, prob=0.9) ]

  sc.runs<-unbalanced_sets_m0

}

# Remove blank runs, if any
if(length(grep("blank",sc.runs))>0){
  sc.runs<-sc.runs[-grep("blank",sc.runs)]
}

# Remove experimental sets concurrent with low mass spec performance
i1<-grep("FP103", sc.runs)
i2<-grep("FP95", sc.runs)
i3<-grep("FP96", sc.runs)

if(length(c(i1,i2,i3))>0){

  remove.sc.runs<-sc.runs[sc.runs%in%sc.runs[c(i1,i2,i3)]]

  ev<-ev[!ev$Raw.file%in%remove.sc.runs, ]

}


# Which columns hold the TMT Reporter ion (RI) data
ri.index<-which(colnames(ev)%in%paste0("Reporter.intensity.",1:16))

# Make sure all runs are described in design, if not, print and remove them:
not.described<-unique(ev$Raw.file)[ !unique(ev$Raw.file) %in% paste0(design$Set) ]
ev<-ev[!ev$Raw.file%in%not.described,]

# Filter out reverse hits, contaminants, and contaminated spectra...
ev<-ev[-which(ev$Reverse=="+"),]
if(length(grep("REV", ev$Leading.razor.protein))>0){ ev<-ev[-grep("REV", ev$Leading.razor.protein),] }
if(length(grep("CON", ev$Leading.razor.protein))>0){ ev<-ev[-grep("CON", ev$Leading.razor.protein),] }
if(length(which(ev$Potential.contaminant=="+"))>0){ ev<-ev[-which(ev$Potential.contaminant=="+"),] }
ev_standard_filtering<-ev[,c("Raw.file","dart_PEP","Leading.razor.protein","PEP")]
ev<-ev[!is.na(ev$PIF),]
ev<-ev[ev$PIF>0.8,]


# Remove sets with less than X peptides from 10, 100 cell runs:
ev<-ev[ !( (ev$Raw.file%in%names(table(ev$Raw.file))[table(ev$Raw.file)<thres_ids_c10c100])&(ev$Raw.file%in%c(c10.runs, c100.runs)) ) , ]

# Remove sets with less than X peptides from single cell runs:
ev<-ev[ !( (ev$Raw.file%in%names(table(ev$Raw.file))[table(ev$Raw.file)<thres_ids_sc])&(ev$Raw.file%in%sc.runs) ), ]

# CHECK
evv<-remove.duplicates(ev, c("Raw.file","Leading.razor.protein"))
table(evv$Raw.file)

# Remove peptides that are more the 10% the intensity of the carrier in the single cell runs (only)
ev<-as.data.frame(ev)
ev$mrri<-0
ev$mrri[ev$Raw.file%in%sc.runs] <- rowMeans(ev[ev$Raw.file%in%sc.runs, ri.index[4:16]] / ev[ev$Raw.file%in%sc.runs, ri.index[1]], na.rm = T)
ev<-ev[ev$mrri < 0.1, ]


# Take DART-ID or spectra peptide identifications from proteins that have been detected with < 1% FDR
ev$runtype<-"other"
ev$runtype[ev$Raw.file%in%sc.runs]<-"sc"
ev$runtype[ev$Raw.file%in%c10.runs]<-"c10"
ev$runtype[ev$Raw.file%in%c100.runs]<-"c100"


if(dart_or_spectra=="dart"){

  ev <- ev %>% group_by(runtype, modseq) %>% mutate(pep_fdr = calc_fdr(dart_PEP))

}

if(dart_or_spectra=="spectra"){

  ev <- ev %>% group_by(runtype, modseq) %>% mutate(pep_fdr = calc_fdr(PEP))

}

ev<-ev[ev$pep_fdr < 0.01,]


# Normalize single cell runs to normalization channel
ev<-as.data.frame(ev)
ev[ev$Raw.file%in%sc.runs, ri.index] <- ev[ev$Raw.file%in%sc.runs, ri.index] / ev[ev$Raw.file%in%sc.runs, ri.index[2]]

if(ref_demo){

  ev[ev$Raw.file%in%sc.runs, ri.index] <- ev[ev$Raw.file%in%sc.runs, ri.index] / rowMeans(ev[ev$Raw.file%in%sc.runs, ri.index], na.rm=T)


}

# Organize data into a more convenient data structure:


# Create empty data frame
ev.melt<-melt(ev[0, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev)[ri.index]) ],
              id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"))

colnames(ev.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")


# Record mapping of cell type to Channel:
ct.v<-c()
qt.v<-c()

# Create a unique ID string
unique.id.numeric<-1:16
unique.id<-paste0("i",unique.id.numeric)

RI_keep<-ri.index

# Give each sample a unique identifier
for(X in unique(ev$Raw.file)){

  # Subset data by X'th experiment
  ev.t<-ev[ev$Raw.file%in%X, ]

  if(is.na(X)){next}

  # Name the RI columns by what sample type they are: carrier, single cell, unused, etc...
  colnames(ev.t)[ri.index]<-paste0(as.character(unlist(design[design$Set==X,-1])),"-", unique.id)

  if(length(RI_keep)>0){

    if( X%in%c(c10.runs, c100.runs) ){

      ev.t[,RI_keep]<-ev.t[,RI_keep] / apply(ev.t[, RI_keep], 1, median.na)

    }

    # Melt it! and combine with other experimental sets
    ev.t.melt<-melt(ev.t[, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev.t)[RI_keep]) ],
                    id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"));

    # Record mapping of cell type to Channel:
    ct.v<-c(ct.v, unique.id[which(ri.index%in%RI_keep)] )
    qt.v<-c(qt.v, colnames(ev)[RI_keep] )

    colnames(ev.t.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")

    ev.melt<-rbind(ev.melt, ev.t.melt)

  }

  # Update unique ID string
  unique.id.numeric<-unique.id.numeric + 16
  unique.id<-paste0("i", unique.id.numeric)

}

c2q<-data.frame(ct.v, qt.v); colnames(c2q)<-c("celltype","channel")

# Grab the unique number associate to each and every cell, carrier channel, and empty channel
ev.melt$id<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(2,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$celltype<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(1,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$id<-as.factor(ev.melt$id)

# Remove duplicate observations of peptides from a single experiment
ev.melt<-remove.duplicates(ev.melt,c("sequence","id") )
ev.melt<-ev.melt[!is.na(ev.melt$protein), ]

# Create additional meta data matrices
ev.melt.uniqueID<-remove.duplicates(ev.melt,"id")
ev.melt.pep<-remove.duplicates(ev.melt, c("sequence","protein") )

# Create data frame of peptides x cells, populated by quantitation
ev.unmelt<-dcast(ev.melt, sequence ~ id, value.var = "quantitation", fill=NA)

# Also create matrix of same shape
ev.matrix<-as.matrix(ev.unmelt[,-1]); row.names(ev.matrix)<-ev.unmelt$sequence

# Replace all 0s with NA
ev.matrix[ev.matrix==0]<-NA
ev.matrix[ev.matrix==Inf]<-NA
ev.matrix[ev.matrix==-Inf]<-NA

# Divide matrix into single cells (including intentional blanks) and carriers
sc_cols<-unique(ev.melt$id[(ev.melt$celltype%in%c("sc_u","sc_m0", "sc_0"))&(ev.melt$Raw.file%in%sc.runs)])
ev.matrix.sc<-ev.matrix[, sc_cols]

# 10 cells
b_cols<-unique(ev.melt$id[ev.melt$celltype%in%c("u_10","m0_10")&(ev.melt$Raw.file%in%c10.runs)])
ev.matrix.10<-ev.matrix[, b_cols]

# 100 cells
d_cols<-unique(ev.melt$id[ev.melt$celltype%in%c("u_100","m0_100")&(ev.melt$Raw.file%in%c100.runs)])
ev.matrix.100<-ev.matrix[, d_cols]


# Save data

if(save_tmp == T){
save(c2q,
     ev.matrix.sc,
     ev.melt,
     ev.melt.pep,
     ev.melt.uniqueID,
     sc.runs,
     c10.runs,
     c100.runs,
     ev.matrix.10,
     ev.matrix.100,
     file="dat/imported.RData" )
}






# Filter single cells ----------------------------------------------------------------------

load("dat/imported.Rdata")

sc.melt<-ev.melt[ev.melt$Raw.file%in%sc.runs,]

xd<-as_tibble( sc.melt )

xd <- xd %>% group_by(id) %>% mutate(med_per_c = median(quantitation, na.rm=T)); length(unique(xd$id))

if(ref_demo==F){
  xd <- xd %>% filter(med_per_c > 1/50); length(unique(xd$id))
}

length(unique(xd$id))

xd$quantitation[(xd$quantitation)==Inf]<-NA
xd$quantitation[(xd$quantitation)==0]<-NA

xd <- xd %>% mutate_if(is.factor, as.character)

xd1 <- xd %>%
  group_by(id) %>%
  mutate(norm_q1 = quantitation / median(quantitation, na.rm=T))

xd2 <- xd1 %>%
  group_by(sequence, Raw.file) %>%
  mutate(norm_q = quantitation / mean(norm_q1, na.rm=T))

xd3<- xd2 %>%
  filter(celltype%in%c("sc_m0", "sc_u","sc_0"))

xd4<- xd3 %>%
  group_by(protein, id) %>%
  mutate(cvq = cv(norm_q))

xd5<- xd4 %>%
  group_by(protein, id) %>%
  mutate(cvn = cvna(norm_q))

xd6<- xd5 %>%
  filter(cvn > 5)

#length(unique(ev.melt$Raw.file[ev.melt$id%in%colnames(ev.matrix.sc.f)]))

xd7<-xd6 %>% group_by(id) %>% mutate(cvm=median(cvq, na.rm=T))

xdf<-xd7

hist(unique(xdf$cvm[xdf$celltype!="sc_0"]), col=rgb(0,1,0,1/4), prob=T, breaks=50, main = "X single cells ", xlab="CV")
hist(unique(xdf$cvm[xdf$celltype=="sc_0"]), col=rgb(1,0,0,1/4), prob=T, add=T, breaks=40)

length(unique(xdf$id[!is.na(xdf$cvm)]))
remove_sets<-unique( ev.melt$Raw.file[ev.melt$id%in%unique(xdf$id[xdf$celltype=="sc_0" & xdf$cvm < 0.35])] )

xdf<-xdf[!xdf$Raw.file%in%remove_sets, ]


#length(unique(xdf$id[xdf$celltype!="sc_0" & xdf$cvm < 0.35]))
kid<-(unique(xdf$id[xdf$celltype!="sc_0" & xdf$cvm < 0.4]))
#kid<-(unique(xdf$id[xdf$cvm < 0.35]))
length(kid)

hist(unique(xdf$cvm[xdf$celltype!="sc_0"]), col=rgb(0,1,0,1/4), prob=T, breaks=50, main = "X single cells ", xlab="CV")
hist(unique(xdf$cvm[xdf$celltype=="sc_0"]), col=rgb(1,0,0,1/4), prob=T, add=T, breaks=40)

sc0_kept<-unique( xdf$id[xdf$celltype=="sc_0" & xdf$cvm < 0.4])
sc0_total<-unique( xdf$id[xdf$celltype=="sc_0"])
sc0rate<-round(length(sc0_kept) / length(sc0_total),2)*100

sc_kept<-unique( xdf$id[xdf$celltype!="sc_0" & xdf$cvm < 0.4])
sc_total<-unique( xdf$id[xdf$celltype!="sc_0"])
scrate<-round(length(sc_kept) / length(sc_total),2)*100

ev.matrix.sc.f<-ev.matrix.sc[,colnames(ev.matrix.sc)%in%kid]; dim(ev.matrix.sc.f)
ev.matrix.sc.f[ev.matrix.sc.f==Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==-Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==0]<-NA

xdf$control<-"sc"
xdf$control[xdf$celltype=="sc_0"]<-"ctl"

my_col3<-c( "black", "purple2")

# Plot!
px<-ggplot(data=xdf, aes(x=cvm)) + geom_density(aes(fill=control, alpha=0.5), adjust=4) + theme_pubr() +
  scale_fill_manual(values=my_col3[c(1,2)]) +
  xlab("Quantification variability") + ylab("Density") + rremove("y.ticks") + rremove("y.text") +
  font("xylab", size=35) +
  font("x.text", size=30) +
  #xlim(c(-0.15, 0.35)) +
  # annotate("text", x=0.27, y= 14, label=paste0(scrate,"% single cells passed"), size=8, color=my_col3[c(2)])+
  # annotate("text", x=0.27, y= 12.5, label=paste0(sc0rate,"% control wells passed"), size=8, color=my_col3[c(1)])+
  annotate("text", x=0.272, y= 14, label=paste0(length(sc_kept)," single cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.265, y= 12, label=paste0(length(sc0_kept)," control wells"), size=10, color=my_col3[c(1)])+
  annotate("text", x=0.5, y= 14, label=paste0(length(sc_total) -length(sc_kept)," single cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.5, y= 12, label=paste0(length(sc0_total) - length(sc0_kept)," control wells"), size=10, color=my_col3[c(1)])+
  #annotate("text", x=0.25, y= 3, label="Macrophage-like", size=6) +
  rremove("legend") + geom_vline(xintercept=0.4, lty=2, size=2, color="gray50")


ggsave("figs/cv_2.pdf", plot=px, device="pdf", width=9, height=5)

if(save_tmp == T){

# Save data
save(c2q,
     ev.matrix.sc,
     ev.melt,
     ev.melt.pep,
     ev.melt.uniqueID,
     sc.runs,
     c10.runs,
     c100.runs,
     ev.matrix.10,
     ev.matrix.100,
     ev.matrix.sc.f,
     xdf,
     file="dat/imported.RData" )

}



# Data transformations ----------------------------------------------------

if(save_tmp == T){
  rm(list = ls())
  load("dat/imported.RData")
  source("code/functions_parameters.R")
}


# Perform normalizations / transformations in multiple steps with visual sanity checks:
b.t<-"FD"
xlim.t<-c(-2,2)
par(mfrow=c(3,3))

# Original data, normalized to reference channel, filtered for failed wells:
t0<-ev.matrix.sc.f


hist(c(t0), breaks=b.t, xlim=xlim.t)

# Column then row normalize by median or mean (see source functions):
t1<-cr_norm(t0)
hist(c(t1), breaks=b.t, xlim=xlim.t)



# Filter for missing data:
t2<-filt.mat.rc(t1, na.row, na.col)
hist(c(t2), breaks=b.t, xlim=xlim.t)



# Log2 transform:
t3<-log2(t2)
t3[t3==Inf]<-NA
t3[t3==-Inf]<-NA
t3[t3==0]<-NA
hist(c(t3), breaks=b.t, xlim=xlim.t)

mean(ncol(t3)-na.count(t3))

# # Collapse to protein level by median:
t3m<-data.frame(t3)
t3m$pep<-rownames(t3)
t3m$prot <- ev.melt.pep$protein[match(t3m$pep, ev.melt.pep$sequence)]
t3m<-melt(t3m, variable.names = c("pep", "prot"))
colnames(t3m) <-c("pep","prot","id","quantitation")
t3m2<- t3m %>% group_by(prot, id) %>% summarize(qp = median(quantitation, na.rm=T))
t4m<-dcast(t3m2, prot ~ id, value.var = "qp", fill=NA)
t4<-as.matrix(t4m[,-1]); row.names(t4)<-t4m[,1]
hist(c(t4), breaks=b.t, xlim=xlim.t)

# Re-column and row normalize:
t4b<-cr_norm_log(t4)
hist(c(t4b), breaks=b.t, xlim=xlim.t)

# Assign to a final variable name:
ev.matrix.sc.f.n<-t4b

mean(ncol(ev.matrix.sc.f.n)-na.count(ev.matrix.sc.f.n))


# Perform similar operations for the 10 and 100-cell data (10-cell data not used in publication)
ev.matrix.10.n<-( cr_norm(ev.matrix.10) )
ev.matrix.100.n<-( cr_norm(ev.matrix.100) )

## Impute single celldata
imp.input<-ev.matrix.sc.f.n
sc.imp <- hknn(imp.input, k.t)
t5<-sc.imp
sum(is.na(sc.imp))
dim(sc.imp)

sc.imp[(is.na(sc.imp))]<-0

# Batch correction with ComBat
batch.N<-table(ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id%in%colnames(sc.imp)])
sc.imp<-sc.imp[,!colnames(sc.imp)%in%ev.melt.uniqueID$id[ev.melt.uniqueID$Raw.file%in%names(batch.N)[batch.N==1]]]

if(ref_demo){
  batch.covs<-ev.melt.uniqueID$Raw.file[ev.melt.uniqueID$id %in% colnames(sc.imp)]
  rr<-c()
  rawx<-c()
  for(Y in unique(batch.covs)){

    xt<-sc.imp[,batch.covs%in%Y]

    vt<-rowVars(xt)

    rawx<-c(rawx, rep(Y, length(which(vt==0))))

    rr<-c(rr, which(vt==0) )

  }

  table(rawx)

  sc.imp<-sc.imp[-rr,]; dim(sc.imp)
}

# Single cells
# Define the batches and model:
batch.covs <- ev.melt.uniqueID$Raw.file[match(colnames(sc.imp), ev.melt.uniqueID$id)]
mod<-data.frame(ev.melt.uniqueID$celltype[match(colnames(sc.imp), ev.melt.uniqueID$id)]); colnames(mod)<-"celltype"
mod<-model.matrix(~as.factor(celltype), data=mod)

matrix.sc.batch <- ComBat(sc.imp, batch=batch.covs, mod=mod, par.prior=T)
t6<-matrix.sc.batch

# visual sanity checks post-imputation:
hist(c(t5), breaks=b.t, xlim=xlim.t)
hist(c(t6), breaks=b.t, xlim=xlim.t)

par(mfrow=c(1,1))

# Save data
save(c2q,
     ev.matrix.sc,
     ev.melt,
     ev.melt.pep,
     ev.melt.uniqueID,
     sc.runs,
     c10.runs,
     c100.runs,
     ev.matrix.10,
     ev.matrix.100,
     ev.matrix.sc.f,
     ev.matrix.sc.f.n,
     sc.imp,
     t3,
     matrix.sc.batch,
     ev.matrix.10.n,
     ev.matrix.100.n,
     file="dat/imported.RData" )





# Determine monocyte and macrophage markers from bulk data ----------------

# Which genes up or down regulated in mono vs. mac in bulk proteomic data?
# Take the 100 cell TMT11-plex data, put on log2 scale and renormalize:

# p100<-collapse_to_protein(log2(ev.matrix.100.n), ev.melt, T)
# p100<-cr_norm_log(p100)
# save(p100, "dat/p100.RData")
load("dat/p100.RData")

# Calculate fold-change between cell types and significance of that fold-change using the distributions of the protein values in
# each of those cell types:

# Initialize variables to store protein, p-value, fold-change:
pc<-c()
pval<-c()
fc<-c()

# Iterate over every protein:
for(X in rownames(p100)){

  # Separate data according to cell type
  m0.id<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="m0_100"]
  u.id<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="u_100"]

  # As long as there is data for both cell types:
  if( !all(is.na(p100[rownames(p100)==X, colnames(p100)%in%m0.id])) & !all(is.na(p100[rownames(p100)==X, colnames(p100)%in%u.id])) ){

    fc.t<-NA
    pval.t<-NA

    pval.t<-t.test( p100[rownames(p100)==X, colnames(p100)%in%m0.id],
                    p100[rownames(p100)==X, colnames(p100)%in%u.id],
                    alternative = "two.sided" )$p.value

    fc.t<- mean(p100[rownames(p100)==X, colnames(p100)%in%m0.id], na.rm = T) -
      mean(p100[rownames(p100)==X, colnames(p100)%in%u.id], na.rm = T)

    pc<-c(pc, X); pval<-c(pval, pval.t); fc<-c(fc, fc.t)

  }

}

# Organize data
df100<-data.frame(pc, pval, fc)

# Take significantly changing proteins
df100<-df100[df100$pval<0.01, ]
df100<-df100[order(df100$fc, decreasing=T),]

# Take top 30 most differentiatial proteins, up and down:
mono_genes<-df100[df100$fc<0, ]
mono_genes<-mono_genes[order(mono_genes$fc, decreasing=F), ]
mono_genes30<-mono_genes$pc[1:30]
mac_genes<-df100[df100$fc>0, ]
mac_genes30<-mac_genes$pc[1:30]





# Bulk vs. single cell macrophage / monocyte ratios -----------------------

mat.input.sc<-cr_norm(filt.mat.rc(ev.matrix.sc.f, 0.6, 0.8))
mat.input.10<-cr_norm(filt.mat.rc(ev.matrix.10, 0.6, 0.8))
mat.input.100<-cr_norm(filt.mat.rc(ev.matrix.100, 0.6, 0.8))


inset.ratios<-function(mat.input, ev.melt.f, m0i, ui){

  rat.p<-c()
  rat.raws<-c()
  for(Y in unique(ev.melt.f$Raw.file)){

    #print(Y)

    ### Get col ids ###
    m0.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == m0i )&(ev.melt.f$Raw.file==Y), "id"]))
    u.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == ui )&(ev.melt.f$Raw.file==Y), "id"]))

    #print(u.ids)

    rat.t<-NA

    if(length(m0.ids)>1 & length(u.ids)>1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / rowMeans(mat.input[,u.ids], na.rm=T)

      # print(rat.t[1])

    }

    if(length(m0.ids)>1 & length(u.ids)==1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / (mat.input[,u.ids])

    }

    rat.raws<-cbind(rat.raws, rat.t)


  }

  rat.p<-rowMeans(rat.raws, na.rm=T)


  names(rat.p)<-rownames(mat.input)
  return(rat.p)
}

c1<-inset.ratios(mat.input.sc, ev.melt[ev.melt$id%in%colnames(mat.input.sc),], "sc_m0","sc_u")
c10<-inset.ratios(mat.input.10, ev.melt[ev.melt$id%in%colnames(mat.input.10),], "m0_10","u_10")
c100<-inset.ratios(mat.input.100, ev.melt[ev.melt$id%in%colnames(mat.input.100),], "m0_100","u_100")


# Compile into a data frame
df.c<-data.frame(c1v<-c1)
df.c<-cbind(df.c, c10v<-c10[ match(names(c1),names(c10)) ] )
df.c<-cbind(df.c, c100v<-c100[ match(names(c1),names(c100)) ] )
colnames(df.c) <- c("c1","c10", "c100")

# Log2 transform
dflc <- log2(df.c[, 1:3])
#dflc <- (df.c[, 1:3])


# Replace any problematic values
dflc[dflc == -Inf] <- NA
dflc[dflc == Inf] <- NA


library(viridisLite)
library(MASS)

k <- 60

x <- dflc[,3]
y <- dflc[,1]

temp.df<-data.frame(x, y)

# Record the positions of the NA values in either x or y
na.v<-c()
for(i in 1:nrow(temp.df)){

  na.v<-c( na.v, any(is.na(temp.df[i,])) )

}

# Remove those points
temp.df.na<-temp.df[!na.v, ]

x<-temp.df.na$x
y<-temp.df.na$y

contour_cols <- viridis(k, alpha = 0.5)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])

}

dens <- get_density(x, y, k)

def.par<-par()

pdf(file="figs/bulk_scatter.pdf", width = 5, height =5)

par(mar=c(6,8,2,2))
plot(x, y, col = contour_cols[findInterval(dens, seq(0, max(dens), length.out = k))], pch = 16, xlab="", ylab="",
     xlim=c(-2.3,3.5), ylim=c(-1.3,2.5),
     #xlim=c(-3,3.5), ylim=c(-3,3.5),
     #xlim=c(-0.2,1), ylim=c(-0.2,1),
     cex.axis=2)

text(0, 2, paste0("  = ", round(cor(x,y),2)), cex=3)
mtext(expression("SCoPE2, log"[2]), side=2, padj = -1.4, cex=3)
mtext(expression("Bulk, log"[2]), side=1, padj = 1.8, cex=3)
#mtext(expression("SCoPE2, m0/u"), side=2, padj = -1.4, cex=3)
#mtext(expression("100-cell bulk, m0/u"), side=1, padj = 1.8, cex=3)
# abline(a=0, b=1)
# abline(lm(y~x), col="black", lty=2)
# print(lm(y~x))


dev.off()


# Bulk vs single cell: signaling proteins ---------------------------------

load("dat/imported.RData")
gn<-read.csv("C:/R-dir/Github/mpop_final/signaling_gset.csv")
gns<-unique(gn$Entry)

TF<-read.delim("C:/R-dir/TFs.txt")
gnstf<-unique(TF$Entry)
xxxt<-unique(ev.melt$sequence[ev.melt$protein%in%gnstf])

gns<-c(as.character(gns), as.character(gnstf))



xx<-gns
#xx<-plow
xxx<-unique(ev.melt$sequence[ev.melt$protein%in%xx])

xxk<-unique( gn[grep("kinase",gn$Protein.names),"Entry"] )
xxxk<-unique(ev.melt$sequence[ev.melt$protein%in%xxk])

xxr<-unique( gn[grep("receptor",gn$Protein.names),"Entry"] )
xxxr<-unique(ev.melt$sequence[ev.melt$protein%in%xxr])


mat.input.sc<-cr_norm(filt.mat.rc(ev.matrix.sc.f[rownames(ev.matrix.sc.f)%in%xxx,], 0.99, 0.99))


#mat.input.sc<-cr_norm(filt.mat.rc(ev.matrix.sc.f, 0.5, 0.8))
mat.input.10<-cr_norm(filt.mat.rc(ev.matrix.10, 0.99, 0.99))
mat.input.100<-cr_norm(filt.mat.rc(ev.matrix.100, 0.99, 0.99))
inset.ratios<-function(mat.input, ev.melt.f, m0i, ui){

  rat.p<-c()
  rat.raws<-c()
  for(Y in unique(ev.melt.f$Raw.file)){

    #print(Y)

    ### Get col ids ###
    m0.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == m0i )&(ev.melt.f$Raw.file==Y), "id"]))
    u.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == ui )&(ev.melt.f$Raw.file==Y), "id"]))

    #print(u.ids)

    rat.t<-rep(NA, nrow(mat.input))


    if(length(m0.ids)>1 & length(u.ids)>1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / rowMeans(mat.input[,u.ids], na.rm=T)

      # print(rat.t[1])

    }

    if(length(m0.ids)>1 & length(u.ids)==1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / (mat.input[,u.ids])

    }

    rat.raws<-cbind(rat.raws, rat.t)


  }

  rownames(rat.raws)<-rownames(mat.input)
  rat.raws2<-filt.mat.rc(rat.raws, (1-(10/ncol(rat.raws))), 1 )

  rat.p<-rowMeans(rat.raws2, na.rm=T)

  names(rat.p)<-rownames(mat.input)[rownames(mat.input)%in%rownames(rat.raws2)]
  return(rat.p)
}

inset.ratios100<-function(mat.input, ev.melt.f, m0i, ui){

  rat.p<-c()
  rat.raws<-c()
  for(Y in unique(ev.melt.f$Raw.file)){

    #print(Y)

    ### Get col ids ###
    m0.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == m0i )&(ev.melt.f$Raw.file==Y), "id"]))
    u.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == ui )&(ev.melt.f$Raw.file==Y), "id"]))

    #print(u.ids)

    rat.t<-NA

    if(length(m0.ids)>1 & length(u.ids)>1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / rowMeans(mat.input[,u.ids], na.rm=T)

      # print(rat.t[1])

    }

    if(length(m0.ids)>1 & length(u.ids)==1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / (mat.input[,u.ids])

    }

    rat.raws<-cbind(rat.raws, rat.t)


  }

  rownames(rat.raws)<-rownames(mat.input)
  rat.raws2<-rat.raws

  rat.p<-rowMeans(rat.raws2, na.rm=T)

  names(rat.p)<-rownames(mat.input)[rownames(mat.input)%in%rownames(rat.raws2)]
  return(rat.p)
}


c1<-inset.ratios(mat.input.sc, ev.melt[ev.melt$id%in%colnames(mat.input.sc),], "sc_m0","sc_u")
#c10<-inset.ratios(mat.input.10, ev.melt[ev.melt$id%in%colnames(mat.input.10),], "m0_10","u_10")
c100<-inset.ratios100(mat.input.100, ev.melt[ev.melt$id%in%colnames(mat.input.100),], "m0_100","u_100")






# Compile into a data frame
df.c<-data.frame(c1v<-c1)
#df.c<-cbind(df.c, c10v<-c10[ match(names(c1),names(c10)) ] )
df.c<-cbind(df.c, c100v<-c100[ match(names(c1),names(c100)) ] )
colnames(df.c) <- c("c1", "c100")

# Log2 transform
dflc <- log2(df.c[, 1:2])
#dflc <- (df.c[, 1:3])


# Replace any problematic values
dflc[dflc == -Inf] <- NA
dflc[dflc == Inf] <- NA


library(viridisLite)
library(MASS)

k <- 60



x <- dflc[,2]
y <- dflc[,1]

temp.df<-data.frame(x, y)

# Record the positions of the NA values in either x or y
na.v<-c()
for(i in 1:nrow(temp.df)){

  na.v<-c( na.v, any(is.na(temp.df[i,])) )

}

# Remove those points
temp.df.na<-temp.df[!na.v, ]

x<-temp.df.na$x
y<-temp.df.na$y

contour_cols <- viridis(k, alpha = 0.5)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])

}

dens <- get_density(x, y, k)

def.par<-par()



cols<- rep(rgb(0.1,0.1,0.1,1/4),length(x))

cols[rownames(df.c)%in%xxxk]<-"red"

cols[rownames(df.c)%in%xxxr]<-"orange"

cols[rownames(df.c)%in%xxxt]<-"green2"

# Take only proteins that have quantification in both bulk and sc to count number of each type:
dfxx<-rowSums(df.c)
dfxxx<-dfxx[!is.na(dfxx)]

nkinase<-length(unique(ev.melt$protein[ev.melt$sequence%in%(names(dfxxx)[names(dfxxx)%in%xxxk])]))
nreceptor<-length(unique(ev.melt$protein[ev.melt$sequence%in%(names(dfxxx)[names(dfxxx)%in%xxxr])]))
ntf<-length(unique(ev.melt$protein[ev.melt$sequence%in%(names(dfxxx)[names(dfxxx)%in%xxxt])]))

pdf(file="figs/signaling_scatter.pdf", width = 6, height =5)


par(xpd=TRUE)
par(mar=c(6,8,2,2))
plot(x, y, col = rgb(0.1,0.1,0.1,1/4), pch = 16, cex=1, xlab="", ylab="",
     xlim=c(-1.3,3), ylim=c(-1.6,3),
     #    xlim=c(-0,4), ylim=c(0,4),
     cex.axis=2)

cols<- rep(NA,length(x))

cols[rownames(df.c)%in%xxxk]<-"red"

cols[rownames(df.c)%in%xxxr]<-"orange"

cols[rownames(df.c)%in%xxxt]<-"green2"

points(x, y, cex=1.5, pch=16, col=cols)

text(0, 2.5, paste0("  = ", round(cor(x,y),2)), cex=3)

text(2.32-0.5-0.1, -1.3, paste0(nreceptor," receptors"), cex=2.5,col="orange")
text(1.9, -0.75, paste0(nkinase, " kinases"), cex=2.5,col="red")
text(2.3, -0.3, paste0(ntf, " TFs"), cex=2.5,col="green2")


mtext(expression("SCoPE2, log"[2]), side=2, padj = -1.4, cex=3)
mtext(expression("Bulk, log"[2]), side=1, padj = 1.8, cex=3)

dev.off()

# Bulk vs. single cell: low abundance proteins ----------------------------

# # Calculating low, medium, and high abundance proteins based on AUC of precursor ion
# evx<-ev[ev$Raw.file%in%sc.runs, ]
# protx<-c()
# aucx<-c()
# for(X in unique(evx$Leading.razor.protein)){
#
#   protx<-c(protx, X)
#   aucx<-c(aucx, median(evx$Intensity[evx$Leading.razor.protein%in%X], na.rm=T))
#
# }
#
# pframe<-data.frame(protx, aucx)
#
# plow<-pframe$protx[pframe$aucx < quantile(pframe$aucx, probs = 0.33,na.rm = T)]
# pmed<-pframe$protx[(pframe$aucx >= quantile(pframe$aucx, probs = 0.33,na.rm = T) ) & (pframe$aucx <= quantile(pframe$aucx, probs = 0.66,na.rm = T) )]
# phigh<-pframe$protx[pframe$aucx > quantile(pframe$aucx, probs = 0.66,na.rm = T)]
#
# save(plow,pframe, file="dat/plow.RData")

load("dat/plow.RData")


xx<-plow
xxx<-unique(ev.melt$sequence[ev.melt$protein%in%xx])

mat.input.sc<-cr_norm(filt.mat.rc(ev.matrix.sc.f[rownames(ev.matrix.sc.f)%in%xxx,], 0.99, 0.99))
mat.input.10<-cr_norm(filt.mat.rc(ev.matrix.10, 0.99, 0.99))
mat.input.100<-cr_norm(filt.mat.rc(ev.matrix.100, 0.99, 0.99))

inset.ratios<-function(mat.input, ev.melt.f, m0i, ui){

  rat.p<-c()
  rat.raws<-c()
  for(Y in unique(ev.melt.f$Raw.file)){

    #print(Y)

    ### Get col ids ###
    m0.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == m0i )&(ev.melt.f$Raw.file==Y), "id"]))
    u.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == ui )&(ev.melt.f$Raw.file==Y), "id"]))

    #print(u.ids)

    rat.t<-rep(NA, nrow(mat.input))

    if(length(m0.ids)>1 & length(u.ids)>1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / rowMeans(mat.input[,u.ids], na.rm=T)

      # print(rat.t[1])

    }

    if(length(m0.ids)>1 & length(u.ids)==1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / (mat.input[,u.ids])

    }

    rat.raws<-cbind(rat.raws, rat.t)


  }

  rownames(rat.raws)<-rownames(mat.input)
  rat.raws2<-filt.mat.rc(rat.raws, (1-(10/ncol(rat.raws))), 1 )

  rat.p<-rowMeans(rat.raws2, na.rm=T)

  names(rat.p)<-rownames(mat.input)[rownames(mat.input)%in%rownames(rat.raws2)]
  return(rat.p)
}



inset.ratios100<-function(mat.input, ev.melt.f, m0i, ui){

  rat.p<-c()
  rat.raws<-c()
  for(Y in unique(ev.melt.f$Raw.file)){

    #print(Y)

    ### Get col ids ###
    m0.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == m0i )&(ev.melt.f$Raw.file==Y), "id"]))
    u.ids <- as.character(unique(ev.melt.f[ (ev.melt.f$celltype == ui )&(ev.melt.f$Raw.file==Y), "id"]))

    #print(u.ids)

    rat.t<-NA

    if(length(m0.ids)>1 & length(u.ids)>1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / rowMeans(mat.input[,u.ids], na.rm=T)

      # print(rat.t[1])

    }

    if(length(m0.ids)>1 & length(u.ids)==1){

      rat.t<-rowMeans(mat.input[,m0.ids], na.rm=T) / (mat.input[,u.ids])

    }

    rat.raws<-cbind(rat.raws, rat.t)


  }

  rownames(rat.raws)<-rownames(mat.input)
  rat.raws2<-rat.raws

  rat.p<-rowMeans(rat.raws2, na.rm=T)

  names(rat.p)<-rownames(mat.input)[rownames(mat.input)%in%rownames(rat.raws2)]
  return(rat.p)
}


c1<-inset.ratios(mat.input.sc, ev.melt[ev.melt$id%in%colnames(mat.input.sc),], "sc_m0","sc_u")
#c10<-inset.ratios(mat.input.10, ev.melt[ev.melt$id%in%colnames(mat.input.10),], "m0_10","u_10")
c100<-inset.ratios100(mat.input.100, ev.melt[ev.melt$id%in%colnames(mat.input.100),], "m0_100","u_100")


# Compile into a data frame
df.c<-data.frame(c1v<-c1)
#df.c<-cbind(df.c, c10v<-c10[ match(names(c1),names(c10)) ] )
df.c<-cbind(df.c, c100v<-c100[ match(names(c1),names(c100)) ] )
colnames(df.c) <- c("c1", "c100")

# Log2 transform
dflc <- log2(df.c[, 1:2])
#dflc <- (df.c[, 1:3])


# Replace any problematic values
dflc[dflc == -Inf] <- NA
dflc[dflc == Inf] <- NA


library(viridisLite)
library(MASS)

k <- 60



x <- dflc[,2]
y <- dflc[,1]

temp.df<-data.frame(x, y)

# Record the positions of the NA values in either x or y
na.v<-c()
for(i in 1:nrow(temp.df)){

  na.v<-c( na.v, any(is.na(temp.df[i,])) )

}

# Remove those points
temp.df.na<-temp.df[!na.v, ]

x<-temp.df.na$x
y<-temp.df.na$y

contour_cols <- viridis(k, alpha = 0.5)

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])

}

dens <- get_density(x, y, k)

def.par<-par()


pframeL<- pframe[pframe$aucx < quantile(pframe$aucx, probs = 0.33),]

par(mar=c(6,8,2,2))
plot(x, y, col = contour_cols[findInterval(dens, seq(0, max(dens), length.out = k))], pch = 16, xlab="", ylab="",
     xlim=c(-1.3,3), ylim=c(-1.3,3),
     #    xlim=c(-0,4), ylim=c(0,4),
     cex.axis=2)




text(0, 2.5, paste0("  = ", round(cor(x,y),2)), cex=3)
mtext(expression("SCoPE2, log"[2]), side=2, padj = -1.4, cex=3)
mtext(expression("Bulk, log"[2]), side=1, padj = 1.8, cex=3)
#mtext(expression("SCoPE2, m0/u"), side=2, padj = -1.4, cex=3)
#mtext(expression("100-cell bulk, m0/u"), side=1, padj = 1.8, cex=3)
# abline(a=0, b=1)
# abline(lm(y~x), col="black", lty=2)
# print(lm(y~x))


cols<- pframeL$aucx[match(unique(ev.melt$protein[ev.melt$sequence%in%rownames(df.c)]), pframeL$protx)]


par(xpd=TRUE)
par(mar=c(6,8,2,8))
plot(x, y, col = colorRampPalette(colors = c('black','yellow'))(findInterval(cols, seq(0, max(cols,na.rm = T), length.out = k))), pch = 16, xlab="", ylab="",
     xlim=c(-2.1,3), ylim=c(-2.1,3),
     #    xlim=c(-0,4), ylim=c(0,4),
     cex.axis=2)

lgd_<-rep(NA,11)
lgd_[c(1,6,11)] = names(quantile(pframeL$aucx[match(unique(ev.melt$protein[ev.melt$sequence%in%rownames(df.c)]), pframeL$protx)],probs=seq(0,1,0.1), na.rm=T))[c(1,6,11)]

legend(x = 3.5, y = 2,
       legend = lgd_,
       fill = colorRampPalette(colors = c('black','yellow'))(11),
       border = NA,
       y.intersp = 0.5,
       cex = 1.5, text.font = 2)



text(0, 2.5, paste0("  = ", round(cor(x,y),2)), cex=3)
mtext(expression("SCoPE2, log"[2]), side=2, padj = -1.4, cex=3)
mtext(expression("Bulk, log"[2]), side=1, padj = 1.8, cex=3)


cols<- pframeL$aucx[match(ev.melt.pep$protein[match(rownames(df.c),ev.melt.pep$sequence)], pframeL$protx)]
pdf(file="figs/lowabundance_scatter.pdf", width = 7, height =5)


par(xpd=TRUE)
par(mar=c(6,8,2,8))
plot(x, y, col = colorRampPalette(colors = c('black','yellow'))(findInterval(cols, seq(0, max(cols,na.rm = T), length.out = k))), pch = 16, xlab="", ylab="",
     xlim=c(-1.5,3.2), ylim=c(-1.5,2.5),
     #    xlim=c(-0,4), ylim=c(0,4),
     cex.axis=2)

lgd_<-rep(NA,11)
lgd_[c(1,6,11)] = names(quantile(pframeL$aucx[match(unique(ev.melt$protein[ev.melt$sequence%in%rownames(df.c)]), pframeL$protx)],probs=seq(0,1,0.1), na.rm=T))[c(1,6,11)]

legend(x = 3.5, y = 2,
       legend = lgd_,
       fill = colorRampPalette(colors = c('black','yellow'))(11),
       border = NA,
       y.intersp = 0.5,
       cex = 1.5, text.font = 2)



text(0, 2.3, paste0("  = ", round(cor(x,y),2)), cex=3)
mtext(expression("SCoPE2, log"[2]), side=2, padj = -1.4, cex=3)
mtext(expression("Bulk, log"[2]), side=1, padj = 1.8, cex=3)
dev.off()

# PCA ------------------------------------------------------------------------

# PC's to display:
PCx<-"PC1"
PCy<-"PC2"

# Data to use:
matrix.sc.batch<-matrix.sc.batch[!is.na(rownames(matrix.sc.batch)), ]
mat.sc.imp<-cr_norm_log(matrix.sc.batch)

# Dot product of each protein correlation vector with itself
r1<-cor(t(matrix.sc.batch))
rsum<-rowSums(r1^2)

# Calculate the weighted data matrix:

X.m <- mat.sc.imp
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- cor(X.m)

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data
pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
add.cols<-colnames(ev.melt)[4:7]
pca.melt[,add.cols]<-NA

for(X in unique(pca.melt$id)){

  pca.melt[pca.melt$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]

}


# Re map ...
pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)

pca.display[,add.cols]<-NA

for(X in unique(pca.display$id)){

  pca.display[pca.display$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]

}


# Map to TMT channel
#pca.display$channel<-c2q[c2q$celltype%in%pca.display$id, "channel"]

# Load in bulk proteomic data to color the PCA plot:
m.m<-mat.sc.imp[rownames(mat.sc.imp)%in%mono_genes30, match(pca.display$id, colnames(mat.sc.imp))]; dim(m.m)
mac.m<-mat.sc.imp[rownames(mat.sc.imp)%in%mac_genes30,match(pca.display$id, colnames(mat.sc.imp)) ]; dim(mac.m)

# Color by median value across replicates
int.m<-rowMedians((t(m.m)), na.rm=T)
int.mac<-rowMedians((t(mac.m)), na.rm=T)
pca.display$mono<-int.m#[length(int.m):1]
pca.display$mac<-int.mac#[length(int.mac):1]


# Display celltype
pg1<-ggscatter(pca.display, x =PCx, y = PCy , color="celltype", size = 2, alpha=0.5) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  rremove("legend") +
  scale_color_manual(values = my_colors[2:3]) +
  annotate("text", x=-0.025, y=0.21,label="Macrophage", color=my_colors[2], size=10)  +
  annotate("text", x=0.03, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
  annotate("text", x=0.05-0.02, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  annotate("text", x=0.062-0.03, y=-0.11, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8)

# Adjust color scale:
varx<-int.m
qs.varx<-quantile(rescale(varx), probs=c(0,0.3,0.5,0.7,1))

# Display celltype
pg2<-ggscatter(pca.display, x =PCx, y = PCy , color="mono", size = 1, alpha=0.7) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  scale_color_gradientn(colors=c("purple","yellow2"), values=qs.varx) +
  #annotate("text", x=-0.04, y=0.21,label="Macrophage", color=my_colors[2], size=10)  +
  # annotate("text", x=0.04, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
  #annotate("text", x=0.05-0.01, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  #annotate("text", x=0.062-0.014, y=-0.13, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) +
  rremove("xylab")+
  rremove("xy.text") +
  rremove("ticks") +
  rremove("legend") +
  ggtitle("Monocyte genes")

# Adjust color scale:
varx<-int.mac
qs.varx<-quantile(rescale(varx), probs=c(0,0.3,0.5,0.7,1))

# Display celltype
pg3<-ggscatter(pca.display, x =PCx, y = PCy , color="mac", size = 1, alpha=0.7) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  scale_color_gradientn(colors=c("purple","yellow2"), values=qs.varx) +
  #annotate("text", x=-0.04, y=0.21,label="Macrophage", color=my_colors[2], size=10)  +
  #annotate("text", x=0.04, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
  #annotate("text", x=0.05-0.01, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  #annotate("text", x=0.062-0.014, y=-0.13, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) +
  rremove("xylab")+
  rremove("xy.text") +
  rremove("ticks") +
  ylab("Macrophage-genes") +
  rremove("legend")+
  ggtitle("Macrophage genes")

# ggscatter(pca.display, x =PCx, y = PCy , color="mac", size = 3, alpha=0.7) +
#   xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
#   ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
#   font("ylab",size=30) +
#   font("xlab",size=30) +
#   font("xy.text", size=20) +
#   scale_color_gradientn(colors=c("purple","yellow2"), values=qs.varx) +
#   #annotate("text", x=-0.04, y=0.21,label="Macrophage", color=my_colors[2], size=10)  +
#   #annotate("text", x=0.04, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
#   #annotate("text", x=0.05-0.01, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
#   #annotate("text", x=0.062-0.014, y=-0.13, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) +
#   rremove("xylab")+
#   rremove("xy.text") +
#   rremove("ticks") +
#   ylab("Macrophage-genes") +
#   #rremove("legend")+
#   ggtitle("Macrophage genes")


px <- pg1 + (pg2 + pg3 + plot_layout(ncol=1, nrow=2, heights=c(1,1))) + plot_layout(ncol=2, nrow=1, widths=c(2,1))

ggsave(px, filename = "figs/pca.pdf", device="pdf", width = 8, height = 5)







# PCA on subsets of proteins: signaling proteins ---------

# PC's to display:
PCx<-"PC1"
PCy<-"PC2"

# Data to use:
mat.sc.imp<-cr_norm_log(matrix.sc.batch[rownames(matrix.sc.batch)%in%gns, ])
dim(mat.sc.imp)

# Dot product of each protein correlation vector with itself
r1<-cor(t(mat.sc.imp))
rsum<-rowSums(r1^2)

# Calculate the weighted data matrix:

X.m <- mat.sc.imp
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- cor(X.m)

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data
pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
add.cols<-colnames(ev.melt)[4:7]
pca.melt[,add.cols]<-NA

for(X in unique(pca.melt$id)){

  pca.melt[pca.melt$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]

}


# Re map ...
pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)

pca.display[,add.cols]<-NA

for(X in unique(pca.display$id)){

  pca.display[pca.display$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]

}


# Display celltype
px<-ggscatter(pca.display, x =PCx, y = PCy , color="celltype", size = 2, alpha=0.5) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  rremove("legend") +
  scale_color_manual(values = my_colors[2:3]) +
  annotate("text", x=-0.025, y=0.21,label="Macrophage", color=my_colors[2], size=10)  +
  annotate("text", x=0.03, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
  annotate("text", x=0.05-0.02, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  annotate("text", x=0.062-0.03, y=-0.11, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8)

ggsave(px, filename = "figs/pca_signaling.pdf", device="pdf", width = 6, height = 5)


# PCA on subsets of proteins: low abundance -------------------------------

# PC's to display:
PCx<-"PC1"
PCy<-"PC2"

# Data to use:
mat.sc.imp<-cr_norm_log(matrix.sc.batch[rownames(matrix.sc.batch)%in%plow, ])
dim(mat.sc.imp)

# Dot product of each protein correlation vector with itself
r1<-cor(t(mat.sc.imp))
rsum<-rowSums(r1^2)

# Calculate the weighted data matrix:

X.m <- mat.sc.imp
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- cor(X.m)

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data
pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
add.cols<-colnames(ev.melt)[4:7]
pca.melt[,add.cols]<-NA

for(X in unique(pca.melt$id)){

  pca.melt[pca.melt$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]

}


# Re map ...
pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)

pca.display[,add.cols]<-NA

for(X in unique(pca.display$id)){

  pca.display[pca.display$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]

}


# Display celltype
px<-ggscatter(pca.display, x =PCx, y = PCy , color="celltype", size = 2, alpha=0.5) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  rremove("legend") +
  scale_color_manual(values = my_colors[2:3]) +
  annotate("text", x=-0.025, y=0.21,label="Macrophage", color=my_colors[2], size=10)  +
  annotate("text", x=0.03, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
  annotate("text", x=0.05-0.02, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  annotate("text", x=0.062-0.03, y=-0.11, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8)

ggsave(px, filename = "figs/pca_lowabundance.pdf", device="pdf", width = 6, height = 5)



# Spectral ordering: macrophages and monocytes ------------------------------------------

# Calculate the Laplacian of the correlation matrix (+1 to all values, so that no values are negative)
W <- 1 + pca.imp.cor
D <- diag(rowSums(W))
L <- D - W

# Get the eigenvalues and eigenvectors of the Laplacian
eigs<-eigen(L)
vec<-eigs$vectors
val<-eigs$values

# Sanity check
plot(val)

# Order the eigenvectors by eigenvalue
vec<-vec[,order(val, decreasing=F)]
vec.m<-vec[,1:2]

# Reformat eigenvector into data frame (taking second eigen vector, first should have eigenvalue 0 and all same values in eigenvector)
vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")

# Plot the eigenvector
pl1<-ggline(vdf, x="num", y="eigen", size=0.001, color="gray60") +
  theme_pubr() +
  rremove("xy.text") +
  rremove("xylab") +
  rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0))

# Record re-ordered eigenvector
mac_mono_vec<-vec.m[,2]

# Reorder the data matrix
mat.c<-mat.sc.imp[, order(vec.m[,2]) ]

# Record rowMedians of the outer 40 cells on either side of matrix
rm1<-rowMedians(mat.c[,1:40], na.rm=T)
rm2<-rowMedians(mat.c[, (ncol(mat.c)-40):ncol(mat.c) ], na.rm=T)

# Create a score to reorder rows for visualization
mfc<-(rm1-rm2)
names(mfc)<-rownames(mat.sc.imp)
mfc2<-mfc; names(mfc2)<-rownames(matrix.sc.batch)
names(mfc)<-1:length(mfc)

mfc.g<-mfc[mfc>0]
mfc.l<-mfc[mfc<0]

mfc.g<-mfc.g[order(mfc.g, decreasing = F)]
mfc.l<-mfc.l[order(mfc.l, decreasing = T)]

mfc.reorder<-c(mfc.g, mfc.l)

names(rm1)<-1:length(rm1)

rm1g<-rm1[names(mfc.g)]
rm1l<-rm1[names(mfc.l)]

rm1g<-rm1g[order(rm1g, decreasing = F)]
rm1l<-rm1l[order(rm1l, decreasing = T)]

rm1.reorder<-c(rm1g,rm1l)

rows.keep<-as.numeric(names(mfc.reorder[abs(mfc.reorder)>quantile(abs(mfc.reorder),probs=0.8)]))

rm.keep<-rm1.reorder[as.numeric(names(rm1.reorder))%in%rows.keep]


# Create a new matrix with re-ordered rows and columns
mat.p<-mat.c[as.numeric(names(rm.keep)),]

# Normalize
mat.p<-cr_norm_log(mat.p)

# mat.p[mat.p>2] <- 2
# mat.p[mat.p< (-2)] <- (-2)

# Reverse order to have monocytes on the left
mat.p<-mat.p[, ncol(mat.p):1]

nrow(mat.p)

# Create the color scale
varx<-melt(mat.p)$value
qs.varx<-quantile(rescale(varx), probs=c(0,0.05,0.5,0.95,1))
qs.varx<-quantile(rescale(varx), probs=c(0,0.10,0.5,0.9,1))

my_col2<-c(rgb(0,0,1,0.5),rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),rgb(1,0,0,0.5))
my_col2<-c("blue","blue","white","red","red")
my_col2<-c("blue",rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),"red")

# Plot!
pl2<-ggplot(melt(mat.p), aes(x=Var2, y=Var1, fill=value))+
  geom_tile() +
  scale_fill_gradientn(colors=my_col2, values=qs.varx) +
  #scale_fill_gradient2(low="blue", mid="white",high="red", midpoint=0) +
  rremove("xy.text") +
  rremove("ticks") +
  ylab("Proteins") +
  xlab("Cells") +
  font("xylab", size=20)


dfx<-data.frame(var1<-ev.melt.uniqueID$celltype[match(colnames(mat.c), ev.melt.uniqueID$id)]); colnames(dfx)<-c("var")
dfx$x<-nrow(dfx):1
dfx$y<-1
vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")
vdf$celltype<-dfx$var[order(dfx$x, decreasing=F)]

# Plot!

pl1<-ggbarplot(vdf, x="num", y="eigen", fill="celltype", color="celltype") + theme_pubr() +
  rremove("xy.text") + rremove("xylab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0)) +
  scale_fill_manual(values=my_colors[c(2,3)])  +
  scale_color_manual(values=my_colors[c(2,3)]) +
  rremove("legend") + rremove("axis")

plot_cells<-ggplot(dfx, aes(x=x,y=y, fill=var))+geom_tile() +
  theme_pubr() +
  rremove("legend") +
  rremove("xylab") +
  rremove("axis") +
  rremove("ticks") +
  rremove("xy.text") +
  scale_fill_manual(values=my_colors[c(2,3)])

# Final figure

px<-pl1 + pl2 + plot_layout(ncol=1, nrow=2, heights = c(2,7))

#ggsave("figs/monomac.pdf", device="pdf", plot=px, width=5, height=7)
#ggsave("figs/monomac.png", device="png", plot=px, width=5, height=7)


save(vec.m, file="dat/vecm.RData")




# Spectral ordering: macrophages only  ------------------------------------------

m0.ind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_m0"]
#m0.ind<-ev.melt.uniqueID$id[intersect(which(ev.melt.uniqueID$celltype=="sc_m0"), grep("FP94",ev.melt.uniqueID$Raw.file))]
mat.sc.imp.m0<-matrix.sc.batch[, colnames(matrix.sc.batch)%in%m0.ind]
mat.sc.imp.m0<-cr_norm_log(mat.sc.imp.m0)

dim(mat.sc.imp.m0)

# Dot product of each protein correlation vector with itself
r1<-cor(t(mat.sc.imp.m0))
rsum<-rowSums(r1^2)

# Calculate the weighted data matrix:
X.m <- mat.sc.imp.m0
X.m <- diag(rsum) %*%  X.m
pca.imp.mac <- cor(X.m)

# Calculate the Laplacian of the correlation matrix (+1 to all values, so that no values are negative)
W <-  1+ (pca.imp.mac)
D <- diag( rowSums(W))
L <- D - W

# Get the eigenvalues and eigenvectors of the Laplacian
eigs<-eigen(L)
vec<-eigs$vectors
val<-eigs$values
plot(val)

# Order the eigenvectors by eigenvalue
vec<-vec[,order(val, decreasing=F)]
vec.m<-vec[,1:2]

# Reformat eigenvector into data frame (taking second eigen vector, first should have eigenvalue 0 and all same values in eigenvector)
vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")

# Plot the eigenvector
pl1<-ggline(vdf, x="num", y="eigen", size=0.001, color="gray60") + theme_pubr() +
  rremove("xy.text") + rremove("xylab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0))


mac_vec<-vec.m[,2]

# Reorder the data matrix
mat.c<-mat.sc.imp.m0[, order(vec.m[,2]) ]

dim(mat.c)

# Record rowMedians of the outer 40 cells on either side of matrix
rm1<-rowMedians(mat.c[,1:40], na.rm=T)
rm2<-rowMedians(mat.c[, (ncol(mat.c)-40):ncol(mat.c) ], na.rm=T)

# Create a score to reorder rows for visualization
mfc<-(rm1-rm2)
names(mfc)<-rownames(mat.sc.imp.m0)
mfc2<-mfc; names(mfc2)<-rownames(matrix.sc.batch); write.csv(mfc2,"fc_2019.csv")
names(mfc)<-1:length(mfc)

mfc.g<-mfc[mfc>0]
mfc.l<-mfc[mfc<0]

mfc.g<-mfc.g[order(mfc.g, decreasing = F)]
mfc.l<-mfc.l[order(mfc.l, decreasing = T)]

mfc.reorder<-c(mfc.g, mfc.l)

names(rm1)<-1:length(rm1)

rm1g<-rm1[names(mfc.g)]
rm1l<-rm1[names(mfc.l)]

rm1g<-rm1g[order(rm1g, decreasing = F)]
rm1l<-rm1l[order(rm1l, decreasing = T)]

rm1.reorder<-c(rm1g,rm1l)

rows.keep<-as.numeric(names(mfc.reorder[abs(mfc.reorder)>quantile(abs(mfc.reorder),probs=0.75)]))

rm.keep<-rm1.reorder[as.numeric(names(rm1.reorder))%in%rows.keep]

# Create a new matrix with re-ordered rows and columns
mat.p<-mat.c[as.numeric(names(rm.keep)),]

# Normalize
mat.p<-cr_norm_log(mat.p)
nrow(mat.p)
# mat.p[mat.p>2] <- 2
# mat.p[mat.p< (-2)] <- (-2)

# Reverse order
mat.p<-mat.p[, ncol(mat.p):1]
mac.spectrum.ps<-rownames(mat.p)

# Create the color scale

varx<-melt(mat.p)$value
qs.varx<-quantile(rescale(varx), probs=c(0,0.05,0.5,0.95,1))
qs.varx<-quantile(rescale(varx), probs=c(0,0.10,0.5,0.9,1))

my_col2<-c("blue",rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),"red")

# Plot!
dim(mat.p)
pl2<-ggplot(melt(mat.p), aes(x=Var2, y=Var1, fill=value))+
  geom_tile() +
  scale_fill_gradientn(colors=my_col2, values=qs.varx) +
  #scale_fill_gradient2(low="blue", mid="white",high="red", midpoint=0) +
  rremove("xy.text") +
  rremove("ticks") +
  ylab("Proteins") +
  xlab("Cells") +
  font("xylab", size=20)

dfx<-data.frame(var1<-ev.melt.uniqueID$celltype[match(colnames(mat.c), ev.melt.uniqueID$id)]); colnames(dfx)<-c("var")
dfx$x<-nrow(dfx):1
dfx$y<-1

vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")

vdf$celltype<-dfx$var[order(dfx$x, decreasing=F)]

pl1<-ggbarplot(vdf, x="num", y="eigen", fill="celltype", color="celltype") + theme_pubr() +
  rremove("xy.text") + rremove("xylab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0)) +
  scale_fill_manual(values=my_colors[c(2,3)])  +
  scale_color_manual(values=my_colors[c(2,3)]) +
  rremove("legend") + rremove("axis")

plot_cells<-ggplot(dfx, aes(x=x,y=y, fill=var))+geom_tile() +
  theme_pubr() +
  rremove("legend") +
  rremove("xylab") +
  rremove("axis") +
  rremove("ticks") +
  rremove("xy.text") +
  scale_fill_manual(values=my_colors[c(2,3)]) +
  theme(plot.margin = margin(0, 1, 0, 5, "cm"))

# Extract markers of M1 and M2 polarization and look at their quantitation across spectrum
M2<-read.csv("dat/M2overM1_Geneset.csv")
M1<-read.csv("dat/M1overM2_Geneset.csv")
mark_map<-read.csv("dat/mark_map.csv")
head(mark_map)

m2marks<-mark_map$uniprot[mark_map$rna%in%M2$M2overM12]
m1marks<-mark_map$uniprot[mark_map$rna%in%M1$M1overM22]

M2.m<-( mat.c[which(rownames(mat.c)%in%as.character(m2marks)), ] )
M1.m<-( mat.c[which(rownames(mat.c)%in%as.character(m1marks)), ] )
m2.v<-rowMedians(t(M2.m), na.rm = T)
m1.v<-rowMedians(t(M1.m), na.rm = T)

dim(M2.m)
dim(M1.m)

m2d<-data.frame(m2.v, 1:length(m2.v))
colnames(m2d)<-c("value", "index")

m1d<-data.frame(m1.v, 1:length(m1.v))
colnames(m1d)<-c("value", "index")

colnames(M2.m)<-1:ncol(M2.m)

colnames(M2.m)<-1:ncol(M2.m)
m2d<-melt(M2.m)

colnames(M1.m)<-1:ncol(M1.m)
m1d<-melt(M1.m)

se<-function(x){ sd(x, na.rm=T)/sqrt(length(x))}

qs<-quantile(m1d$Var2, probs=seq(0,1,0.1))[-1]

m1d$quantile<-NA

for(i in round(qs[order(qs, decreasing = T)],0)){

  print(i)
  m1d$quantile[m1d$Var2%in%1:i]<-i

}

qs<-quantile(m2d$Var2, probs=seq(0,1,0.1))[-1]

m2d$quantile<-NA

for(i in round(qs[order(qs, decreasing = T)],0)){

  print(i)
  m2d$quantile[m2d$Var2%in%1:i]<-i

}

sdf<-data.frame( aggregate(value~quantile, data=m1d, FUN = median), aggregate(value~quantile, data=m1d, FUN = se) )

colnames(sdf)<-c("quantile","median","quantile1","se")

sdf$quantile1<-sdf$quantile1-13


rhom1<-cor(sdf$quantile1, sdf$median)

sp1<-ggscatter(sdf, x="quantile1",y="median") + geom_smooth(method = "lm") +
  geom_errorbar(data = sdf, aes(x = quantile1, y = median, ymin = median - se, ymax = median + se)) +
  rremove("xy.text") +
  rremove("xlab") +
  rremove("ticks")+
  xlim(c(0, ncol(mat.p)))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(-0.1,0,0.15), labels=c("-100%", "0%", "100%") ) +
  rremove("legend")  +  geom_text(x=200, y=0.1, label=paste0("r = ",round(rhom1,2)), size=7) +
  ylab("M1-genes") + font("ylab", size=20) #+ geom_hline(yintercept=0, color="black", size=1)


sdf<-data.frame( aggregate(value~quantile, data=m2d, FUN = median), aggregate(value~quantile, data=m1d, FUN = se) )


colnames(sdf)<-c("quantile","median","quantile1","se")
rhom2<-cor(sdf$quantile1, sdf$median)

sdf$quantile1<-sdf$quantile1-13

sp2<-ggscatter(sdf, x="quantile1",y="median") + geom_smooth(method = "lm") +
  geom_errorbar(data = sdf, aes(x = quantile1, y = median, ymin = median - se, ymax = median + se)) +
  rremove("xy.text") +
  rremove("xlab") +
  rremove("ticks")+
  xlim(c(0, ncol(mat.p)))+
  scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(breaks=c(-0.1,0,0.15), labels=c("-100%", "0%", "100%") ) +
  rremove("legend")  +  #geom_text(x=50, y=0.1, label=paste0("r = ",round(rhom2,2)), size=7) +
  ylab("M2-genes") + font("ylab", size=20) #+ geom_hline(yintercept=0, color="black", size=1)

print(rhom1); print(rhom2)

# As a summary point and standard error over bins of 26 cells:
px<- pl1 + pl2 + sp1 + sp2 + plot_layout(ncol=1, nrow=4, heights = c(1,5,1.5,1.5))

ggsave("figs/mac.pdf", device="pdf", plot=px, width=5, height=7)
ggsave("figs/mac.png", device="png", plot=px, width=5, height=7)



# Fiedler vector values distributions  ------------------------------------------------------------------------


# # Eigenvectors used to perform spectral clustering on both macrophage-like + monocyte, and macrophage-like only data matrices:
# mac_vec<-c(mac_vec, rep(NA, length(mac_mono_vec) - length(mac_vec)))
# dfxx<-data.frame(mac_mono_vec, mac_vec)
# dfxm<-melt(dfxx)
#
# # Plot!
# px<-ggplot(data=dfxm, aes(x=value, y=variable)) + geom_density_ridges(aes(fill=variable)) + theme_pubr() +
#   scale_fill_manual(values=my_colors[c(1,2)]) +
#   xlab("Eigenvector value") + ylab("Density") + rremove("y.ticks") + rremove("y.text") +
#   font("xylab", size=20) +
#   font("x.text", size=20) +
#   #xlim(c(-0.15, 0.35)) +
#   #annotate("text", x=0.25, y= 1.5, label="Monocyte and\n macrophage-like", size=6)+
#   #annotate("text", x=0.25, y= 3, label="Macrophage-like", size=6) +
#   rremove("legend")
#
# ggsave("figs/fiedler.png", plot=px, width=5, height=5)



# Coverage 1/2 ----------------------------------------------------------------

# Spectra only:
ev_standard_filtering_so<-ev_standard_filtering[ev_standard_filtering$Raw.file %in% ev.melt$Raw.file[ev.melt$id%in%colnames(ev.matrix.sc.f)], ]
ev_standard_filtering_so$fdr<-calc_fdr(ev_standard_filtering_so$PEP)
ev_standard_filtering_so<-ev_standard_filtering_so[ev_standard_filtering_so$fdr < 0.01, ]
ev_std_prot_so<-remove.duplicates(ev_standard_filtering_so, c("Raw.file","Leading.razor.protein"))
peptides_per_cell<-as.numeric(table(ev_standard_filtering_so$Raw.file))
proteins_per_cell<-as.numeric(table(ev_std_prot_so$Raw.file))
coverage_df<-data.frame(peptides_per_cell, proteins_per_cell)
colnames(coverage_df)<-c("Peptides","Proteins")
coverage_df_melt<-melt(coverage_df)
head(coverage_df_melt)

#write.csv(coverage_df_melt, "dat/coverage_spectra_only.csv")


# Dart ID:
ev_standard_filtering<-ev_standard_filtering[ev_standard_filtering$Raw.file %in% ev.melt$Raw.file[ev.melt$id%in%colnames(ev.matrix.sc.f)], ]
ev_standard_filtering$fdr<-calc_fdr(ev_standard_filtering$dart_PEP)
ev_standard_filtering<-ev_standard_filtering[ev_standard_filtering$fdr < 0.01, ]
ev_std_prot<-remove.duplicates(ev_standard_filtering, c("Raw.file","Leading.razor.protein"))

peptides_per_cell<-as.numeric(table(ev_standard_filtering$Raw.file))
proteins_per_cell<-as.numeric(table(ev_std_prot$Raw.file))

ev_standard_filtering2<-ev[ev$Raw.file %in% ev.melt$Raw.file[ev.melt$id%in%colnames(ev.matrix.sc.f)], ]
ev_std_prot2<-remove.duplicates(ev_standard_filtering2, c("Raw.file","Leading.razor.protein"))
peptides_per_cell_f<-as.numeric(table(ev_standard_filtering2$Raw.file))
proteins_per_cell_f<-as.numeric(table(ev_std_prot2$Raw.file))

coverage_df<-data.frame(peptides_per_cell, peptides_per_cell_f, proteins_per_cell, proteins_per_cell_f)
colnames(coverage_df)<-c("Peptides", "Peptides,\n filtered","Proteins", "Proteins,\n filtered")
coverage_df_melt<-melt(coverage_df)
coverage_df_melt$type<-"1% FDR"; coverage_df_melt$type[grep("filtered", coverage_df_melt$variable)]<-"strict filtering"

px<-ggboxplot(coverage_df_melt, x="variable", y="value", fill="type") +
  theme(text=element_text(size=20)) +
  xlab("")+
  ylab("# quantified / run\n") +
  font("ylab",size=30)+
  font("xlab",size=30)+
  scale_fill_manual(values = c("white","lightpink") ) +
  rremove("x.ticks") +
  font("x.text", size=20) +
  theme(axis.text.x  = element_text(angle=30, vjust=0.5)) +
  theme(legend.position = c(0.7, 0.9)) +
  theme(legend.background = element_rect(fill="white",
                                         size=1, linetype="solid",
                                         colour ="white")) +
  font("legend.title", size= 20) +
  font("legend.text", size= 20) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.key.size = unit(1.2,"line")) +
  rremove("legend.title") +
  rremove("x.ticks") +
  rremove("legend.title") +
  rremove("xlab") +
  ylim(0,3000)

#write.csv(coverage_df_melt, "dat/coverage.csv")

ggsave(px, filename = "figs/coverage.pdf", device="pdf", width = 5, height = 10)

# Individual protein distributions from enriched GO terms -----------------

library(GSA)

mat.sc.imp<-cr_norm_log(ev.matrix.sc.f.n)

entrez_2_unip<-read.delim("dat/FUNS.tab")
row.names(mat.sc.imp)<-entrez_2_unip$To[match(row.names(mat.sc.imp), entrez_2_unip$From)]

gs<-GSA.read.gmt("dat/c5.all.v6.2.entrez.gmt")

table(entrez_2_unip$To%in%unlist(gs$genesets))

m0.ids <- ev.melt.uniqueID[ev.melt.uniqueID$celltype == "sc_m0", "id"]
u.ids <- ev.melt.uniqueID[ev.melt.uniqueID$celltype == "sc_u", "id"]

u.enriched<-c("GO_CHROMATOID_BODY", "GO_HISTONE_KINASE_ACTIVITY", "GO_NUCLEOLAR_PART")
m0.enriched<-c("GO_INTERMEDIATE_FILAMENT_BINDING", "GO_REGULATION_OF_RECEPTOR_BINDING","GO_REGULATION_OF_GLUCONEOGENESIS" )

quant<-c()
celltypes<-c()
GO<-c()
prot_topx<-c()
for(X in c(u.enriched, m0.enriched)){

  i <- which(gs$geneset.names==X)

  prots.t<-sapply(gs$genesets[i], paste0)

  mat.t<-mat.sc.imp[rownames(mat.sc.imp)%in%prots.t, ]
  print(nrow(mat.t))
  fct_top<-0
  jk<-NA
  for(j in 1:nrow(mat.t)){



    # c1<-c(mat.t[j,colnames(mat.sc.imp)%in%u.ids])
    # c2<-c(mat.t[j,colnames(mat.sc.imp)%in%m0.ids])
    c1<-c(rank(mat.t[j,],na.last=NA)[colnames(mat.sc.imp)%in%u.ids])
    c2<-c(rank(mat.t[j,],na.last=NA)[colnames(mat.sc.imp)%in%m0.ids])

    fct<-abs(median(c2, na.rm = T) - median(c1, na.rm = T))

    if(fct > fct_top){

      fct_top<-fct
      jk<-j
      prot_top<-which(rownames(mat.sc.imp)==rownames(mat.t)[j])

    }
  }

  c1<-NA; c2<-NA


  c1<-c(mat.t[jk,colnames(mat.sc.imp)%in%u.ids])
  c2<-c(mat.t[jk,colnames(mat.sc.imp)%in%m0.ids])

  quant<-c(quant, c1, c2)
  celltypes<-c(celltypes, rep("u",length(c1)), rep("m0",length(c2)))
  GO<-c(GO, rep(X, length(c(c1,c2))))
  prot_topx<-c(prot_topx, rep(prot_top, length(c(c1,c2))))


}

selectdf<-data.frame(GO,celltypes,quant, prot_topx)
selectdf$prot_topx<- rownames(ev.matrix.sc.f.n)[selectdf$prot_topx]
selectdf$GO <- factor(selectdf$GO, levels = c(u.enriched, m0.enriched ))
selectdf$prot_topx <- factor(selectdf$prot_topx, levels = c(unique(selectdf$prot_topx)))

# Plot!
px<-ggplot(data=selectdf, aes(x=quant, y=prot_topx)) +
  geom_density_ridges(aes(fill=celltypes, alpha=0.5), scale=1) +
  theme_pubr() +
  scale_fill_manual(values=my_colors[c(2,3)]) +
  scale_x_continuous(limits=c(-2.5,2.5)) +
  coord_flip() +
  xlab("Protein level, log2") +
  ylab("Gene ontology") +
  rremove("y.ticks") +
  #rremove("x.text") +
  rremove("xlab") +
  font("ylab", size=30) +
  font("xy.text", size=30) +
  rremove("legend") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot=px, "figs/go.pdf", width = 10, height = 5)

# Loading RNA data ----------------------------------------------------------

# ldl<-read.delim("dat/ldl.conlnames.txt", header = F)
#
# library(Matrix)
# r1<-readMM("dat/matrix1.mtx")
# r2<-readMM("dat/matrix2.mtx")
#
# r1r<-read.delim("dat/s1_features.tsv", header=F)
# r2r<-read.delim("dat/s2_features.tsv", header=F)
#
# r1c<-read.delim("dat/s1_barcodes.tsv", header=F)
# r2c<-read.delim("dat/s2_barcodes.tsv", header=F)
#
#
# r1<-as.matrix(r1)
# r2<-as.matrix(r2)
# colnames(r1)<-r1c$V1
# colnames(r2)<-r2c$V1
#
# rownames(r1)<-r1r$V2
# rownames(r2)<-r2r$V2
#
# all(rownames(r1)==rownames(r2))
#
# rmx<-cbind(r1[,paste0("RNA1_",colnames(r1))%in%ldl$V1],r2[,paste0("RNA2_",colnames(r2))%in%ldl$V1])
# sum(paste0("RNA1_",colnames(r1))%in%ldl$V1); sum(paste0("RNA2_",colnames(r2))%in%ldl$V1)
# #502 cells from r1, 516 cells from r2 (for batch correction purposes)
# dim(rmx)
#
#save(rmx, file="dat/rna.RData")

load("dat/rna.RData")
load("dat/imported.RData")

gn_uni<-read.csv("dat/gn_uni.csv")

kg<-as.character(gn_uni$gn[gn_uni$uniprot%in%rownames(matrix.sc.batch)])

rm.f2<-rmx[rownames(rmx)%in%kg,]

rm.f2<-rm.f2[!duplicated(rownames(rm.f2)), ]

rm.f2[rm.f2==0]<-NA

rm.f3<-cr_norm_log(log2(rm.f2))

rm.i<-hknn(rm.f3,3)

# Define the batches and model:
batch.covs <- c(rep("A", 502), rep("B",516))

rr<-c()
rawx<-c()
for(Y in unique(batch.covs)){

  xt<-rm.i[,batch.covs%in%Y]

  vt<-rowVars(xt)

  rawx<-c(rawx, rep(Y, length(which(vt==0))))

  rr<-c(rr, which(vt==0) )

}

table(rawx)

rm.i<-rm.i[-rr,]

rm.batch <- ComBat(rm.i, batch=batch.covs, par.prior=T)

rm.f2<-rm.f2[rownames(rm.f2)%in%rownames(rm.batch), ]
rm.f3<-rm.f3[rownames(rm.f3)%in%rownames(rm.batch), ]

save(gn_uni, rm.f2, rm.f3,rm.i,rm.batch,file="dat/rna_dat.RData")



# PCA (RNA) ---------------------------------------------------------------

# PC's to display:
PCx<-"PC1"
PCy<-"PC2"

# Data to use:
mat.sc.imp<-cr_norm_log(rm.batch)
mat.sc.imp<-rmat


# Dot product of each protein correlation vector with itself
r1<-cor(t(mat.sc.imp))
rsum<-rowSums(r1^2)

X.m <- mat.sc.imp
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- cor(X.m)

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data
pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")

# Re map ...
pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)


# Display celltype
ggscatter(pca.display, x =PCx, y = PCy, size = 2, alpha=0.3) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  rremove("legend") +
  annotate("text", x=0.05-0.03, y=-0.155, label=paste0(dim(mat.sc.imp)[1], " RNA"), size=8) +
  annotate("text", x=0.062-0.034, y=-0.13, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8)




# Merge imputed and batch-corrected Protein and RNA data ----------------------------------------------

load("dat/rna_dat.RData")

pmat<-cr_norm_log(matrix.sc.batch)
rmat<-cr_norm_log(rm.batch)

gnk<-gn_uni$gn[gn_uni$uniprot%in%rownames(pmat)]
length(which(rownames(rmat)%in%gnk))

length(unique(gn_uni$uniprot)) / length(unique(gn_uni$gn)) # ugh

gn_uni.f<-gn_uni[gn_uni$uniprot%in%rownames(pmat),]

row.names(rmat)<-gn_uni.f$uniprot[match(rownames(rmat), gn_uni.f$gn)]

pmat<-pmat[rownames(pmat)%in%rownames(rmat), ]
rmat<-rmat[rownames(rmat)%in%rownames(pmat), ]

dim(rmat)
dim(pmat)

rmat<-rmat[match(rownames(pmat), rownames(rmat) ), ]

# Compute fiedler vectors:
mat.sc.imp<-cr_norm_log(pmat)
r1<-cor(t(mat.sc.imp))
rsum<-rowSums(r1^2)
X.m <- mat.sc.imp
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- cor(X.m)
W <- 1 + pca.imp.cor
D <- diag(rowSums(W))
L <- D - W
eigs<-eigen(L)
vec<-eigs$vectors
val<-eigs$values
vec<-vec[,order(val, decreasing=F)]
vec.m<-vec[,1:2]
mac_mono_vec_p<-vec.m[,2]

# Compute fiedler vectors:
mat.sc.imp<-cr_norm_log(rmat)
r1<-cor(t(mat.sc.imp))
rsum<-rowSums(r1^2)
X.m <- mat.sc.imp
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- cor(X.m)
W <- 1 + pca.imp.cor
D <- diag(rowSums(W))
L <- D - W
eigs<-eigen(L)
vec<-eigs$vectors
val<-eigs$values
vec<-vec[,order(val, decreasing=F)]
vec.m<-vec[,1:2]
mac_mono_vec_r<-vec.m[,2]

mac_mono_vec_p

pmat<-cr_norm_log(pmat[,order(mac_mono_vec_p)])
rmat<-cr_norm_log(rmat[,order(mac_mono_vec_r)])

all(rownames(pmat)==rownames(rmat))

prmat<-cbind(pmat,rmat)
dim(prmat)
write.csv(prmat, "dat/prot_rna_mat_ordered_2.csv")

# Macrophage polarization in RNA data  ------------------------------------

con3h.e<-read.delim("dat/con3h.embedding.txt")

plot(con3h.e$V1, con3h.e$V2)

maclike<-unlist(str_split(rownames(con3h.e)[con3h.e$V2>1], "_"))[seq(2, 2*length(unlist(str_split(rownames(con3h.e)[con3h.e$V2>1], "_"))),2)]
monolike<-unlist(str_split(rownames(con3h.e)[con3h.e$V2< (-1)], "_"))[seq(2, 2*length(unlist(str_split(rownames(con3h.e)[con3h.e$V2>1], "_"))),2)]



rmat.mac<-rmat[,colnames(rmat)%in%maclike]
dim(rmat.mac)

mat.sc.imp.m0<-cr_norm_log(rmat.mac)

# Dot product of each protein correlation vector with itself
r1<-cor(t(mat.sc.imp.m0))
rsum<-rowSums(r1^2)

# Calculate the weighted data matrix:
X.m <- mat.sc.imp.m0
X.m <- diag(rsum) %*%  X.m
pca.imp.mac <- cor(X.m)

# Calculate the Laplacian of the correlation matrix (+1 to all values, so that no values are negative)
W <-  1+ (pca.imp.mac)
D <- diag( rowSums(W))
L <- D - W

# Get the eigenvalues and eigenvectors of the Laplacian
eigs<-eigen(L)
vec<-eigs$vectors
val<-eigs$values
plot(val)

# Order the eigenvectors by eigenvalue
vec<-vec[,order(val, decreasing=F)]
vec.m<-vec[,1:2]

# Reformat eigenvector into data frame (taking second eigen vector, first should have eigenvalue 0 and all same values in eigenvector)
vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")

# Plot the eigenvector
pl1<-ggline(vdf, x="num", y="eigen", size=0.001, color="gray60") + theme_pubr() +
  rremove("xy.text") + rremove("xylab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0))


mac_vec<-vec.m[,2]

# Reorder the data matrix
mat.c<-mat.sc.imp.m0[, order(vec.m[,2]) ]

dim(mat.c)

# Record rowMedians of the outer 40 cells on either side of matrix
rm1<-rowMedians(mat.c[,1:40], na.rm=T)
rm2<-rowMedians(mat.c[, (ncol(mat.c)-40):ncol(mat.c) ], na.rm=T)

# Create a score to reorder rows for visualization
mfc<-(rm1-rm2)
names(mfc)<-rownames(mat.sc.imp.m0)
mfc2<-mfc; names(mfc2)<-rownames(mat.sc.imp.m0)
names(mfc)<-1:length(mfc)

mfc.g<-mfc[mfc>0]
mfc.l<-mfc[mfc<0]

mfc.g<-mfc.g[order(mfc.g, decreasing = F)]
mfc.l<-mfc.l[order(mfc.l, decreasing = T)]

mfc.reorder<-c(mfc.g, mfc.l)

names(rm1)<-1:length(rm1)

rm1g<-rm1[names(mfc.g)]
rm1l<-rm1[names(mfc.l)]

rm1g<-rm1g[order(rm1g, decreasing = F)]
rm1l<-rm1l[order(rm1l, decreasing = T)]

rm1.reorder<-c(rm1g,rm1l)

rows.keep<-as.numeric(names(mfc.reorder[abs(mfc.reorder)>quantile(abs(mfc.reorder),probs=0.75)]))

rm.keep<-rm1.reorder[as.numeric(names(rm1.reorder))%in%rows.keep]

# Create a new matrix with re-ordered rows and columns
mat.p<-mat.c[as.numeric(names(rm.keep)),]

# Normalize
mat.p<-cr_norm_log(mat.p)
dim(mat.p)
# mat.p[mat.p>2] <- 2
# mat.p[mat.p< (-2)] <- (-2)

# Reverse order
mat.p<-mat.p[, ncol(mat.p):1]
mac.spectrum.ps<-rownames(mat.p)

# Create the color scale

varx<-melt(mat.p)$value
qs.varx<-quantile(rescale(varx), probs=c(0,0.05,0.5,0.95,1))
qs.varx<-quantile(rescale(varx), probs=c(0,0.10,0.5,0.9,1))

my_col2<-c("blue",rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),"red")

# Plot!
dim(mat.p)
pl2<-ggplot(melt(mat.p), aes(x=Var2, y=Var1, fill=value))+
  geom_tile() +
  scale_fill_gradientn(colors=my_col2, values=qs.varx) +
  #scale_fill_gradient2(low="blue", mid="white",high="red", midpoint=0) +
  rremove("xy.text") +
  rremove("ticks") +
  ylab("Proteins") +
  xlab("Cells") +
  font("xylab", size=20)

dfx<-data.frame(var1<-ev.melt.uniqueID$celltype[match(colnames(mat.c), ev.melt.uniqueID$id)]); colnames(dfx)<-c("var")
dfx$x<-nrow(dfx):1
dfx$y<-1

vdf<-data.frame(vec.m[order(vec.m[,2], decreasing=T),2]); vdf$num<-1:nrow(vdf); colnames(vdf)<-c("eigen","num")

vdf$celltype<-dfx$var[order(dfx$x, decreasing=F)]

pl1<-ggbarplot(vdf, x="num", y="eigen", fill="#048ABF", color="#048ABF") + theme_pubr() +
  rremove("xy.text") + rremove("xylab") + rremove("ticks")+
  scale_x_continuous( limits = c(min(vdf$num)-1, max(vdf$num))+1, expand = c(0, 0)) +
  scale_fill_manual(values=my_colors[c(2,3)])  +
  scale_color_manual(values=my_colors[c(2,3)]) +
  rremove("legend") + rremove("axis")

plot_cells<-ggplot(dfx, aes(x=x,y=y, fill=var))+geom_tile() +
  theme_pubr() +
  rremove("legend") +
  rremove("xylab") +
  rremove("axis") +
  rremove("ticks") +
  rremove("xy.text") +
  scale_fill_manual(values=my_colors[c(2,3)]) +
  theme(plot.margin = margin(0, 1, 0, 5, "cm"))

# Extract markers of M1 and M2 polarization and look at their quantitation across spectrum
M2<-read.csv("dat/M2overM1_Geneset.csv")
M1<-read.csv("dat/M1overM2_Geneset.csv")
mark_map<-read.csv("dat/mark_map.csv")
head(mark_map)

m2marks<-mark_map$uniprot[mark_map$rna%in%M2$M2overM12]
m1marks<-mark_map$uniprot[mark_map$rna%in%M1$M1overM22]

M2.m<-( mat.c[which(rownames(mat.c)%in%as.character(m2marks)), ] )
M1.m<-( mat.c[which(rownames(mat.c)%in%as.character(m1marks)), ] )
m2.v<-rowMedians(t(M2.m), na.rm = T)
m1.v<-rowMedians(t(M1.m), na.rm = T)

dim(M2.m)
dim(M1.m)

m2d<-data.frame(m2.v, 1:length(m2.v))
colnames(m2d)<-c("value", "index")

m1d<-data.frame(m1.v, 1:length(m1.v))
colnames(m1d)<-c("value", "index")

colnames(M2.m)<-1:ncol(M2.m)

colnames(M2.m)<-1:ncol(M2.m)
m2d<-melt(M2.m)

colnames(M1.m)<-1:ncol(M1.m)
m1d<-melt(M1.m)

se<-function(x){ sd(x, na.rm=T)/sqrt(length(x))}

qs<-quantile(m1d$Var2, probs=seq(0,1,0.1))[-1]

m1d$quantile<-NA

for(i in round(qs[order(qs, decreasing = T)],0)){

  print(i)
  m1d$quantile[m1d$Var2%in%1:i]<-i

}

qs<-quantile(m2d$Var2, probs=seq(0,1,0.1))[-1]

m2d$quantile<-NA

for(i in round(qs[order(qs, decreasing = T)],0)){

  print(i)
  m2d$quantile[m2d$Var2%in%1:i]<-i

}

sdf<-data.frame( aggregate(value~quantile, data=m1d, FUN = median), aggregate(value~quantile, data=m1d, FUN = se) )

colnames(sdf)<-c("quantile","median","quantile1","se")

sdf$quantile1<-sdf$quantile1-13


rhom1<-cor(sdf$quantile1, sdf$median)

sp1<-ggscatter(sdf, x="quantile1",y="median") + geom_smooth(method = "lm") +
  geom_errorbar(data = sdf, aes(x = quantile1, y = median, ymin = median - se, ymax = median + se)) +
  rremove("xy.text") +
  rremove("xlab") +
  rremove("ticks")+
  xlim(c(0, ncol(mat.p)))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks=c(-0.1,0,0.15), labels=c("-100%", "0%", "100%") ) +
  rremove("legend")  +  geom_text(x=200, y=0.1, label=paste0("r = ",round(rhom1,2)), size=7) +
  ylab("M1-genes") + font("ylab", size=20) #+ geom_hline(yintercept=0, color="black", size=1)


sdf<-data.frame( aggregate(value~quantile, data=m2d, FUN = median), aggregate(value~quantile, data=m1d, FUN = se) )


colnames(sdf)<-c("quantile","median","quantile1","se")
rhom2<-cor(sdf$quantile1, sdf$median)

sdf$quantile1<-sdf$quantile1-13

sp2<-ggscatter(sdf, x="quantile1",y="median") + geom_smooth(method = "lm") +
  geom_errorbar(data = sdf, aes(x = quantile1, y = median, ymin = median - se, ymax = median + se)) +
  rremove("xy.text") +
  rremove("xlab") +
  rremove("ticks")+
  xlim(c(0, ncol(mat.p)))+
  scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(breaks=c(-0.1,0,0.15), labels=c("-100%", "0%", "100%") ) +
  rremove("legend")  +  #geom_text(x=50, y=0.1, label=paste0("r = ",round(rhom2,2)), size=7) +
  ylab("M2-genes") + font("ylab", size=20) #+ geom_hline(yintercept=0, color="black", size=1)

print(rhom1); print(rhom2)

# As a summary point and standard error over bins of 26 cells:
px<- pl1 + pl2 + sp1 + sp2 + plot_layout(ncol=1, nrow=4, heights = c(1,5,1.5,1.5))

#ggsave("figs/mac_rna.pdf", device="pdf", plot=px, width=5, height=7)
ggsave("figs/mac_rna.png", device="png", plot=px, width=5, height=7)

dim(mat.p)

# Ion counting ------------------------------------------------------------

load("dat/imported.RData")

# Import data (Proteome Discoverer (demo) search, confident PSM-level output, 1%FDR):
tev<-read.delim("dat/pd_ioncounts_PSMs.txt")
rmat<-rmx
gn_uni.f<-gn_uni[gn_uni$uniprot%in%tev$Master.Protein.Accessions,]
row.names(rmat)<-gn_uni.f$uniprot[match(rownames(rmat), gn_uni.f$gn)]
rmat<-rmat[rownames(rmat)%in%(tev$Master.Protein.Accessions), ]
dim(rmat)
tev$Raw.file<-paste0("X", unlist(str_split(tev$Spectrum.File, fixed(".")))[seq(1, length(unlist(str_split(tev$Spectrum.File, fixed(".")))), 2)])

# Take only single cell runs used in analysis
tev<-tev[tev$Raw.file%in%ev.melt$Raw.file[ev.melt$id%in%colnames(matrix.sc.batch)],]

table(tev$Raw.file)

# Take only proteins seen in mRNA data and final protein data
tev<-tev[tev$Master.Protein.Accessions%in%rownames(rmat), ]


ev.meltq<-ev.melt[ev.melt$id%in%colnames(matrix.sc.batch), ]
ev.meltq$ch<-(c2q$channel[match(ev.meltq$id, c2q$celltype)])
head(ev.meltq)

colnames(tev)[30:40]<-paste0("Reporter.intensity.",1:11)


tevm<-melt(tev, id.vars = colnames(tev)[-30:-40])

tevm<-tevm[tevm$variable%in%paste0("Reporter.intensity.",1:11)[4:11], ]

length(unique(tevm$Master.Protein.Accessions))


tevm$rq<-paste0(tevm$Raw.file,tevm$variable)
ev.meltq$rq<-paste0(ev.meltq$Raw.file,ev.meltq$ch)


tevm<-tevm[tevm$rq%in%ev.meltq$rq, ]

length(unique(tevm$Master.Protein.Accessions))


rmat<-rmat[rownames(rmat)%in%(tev$Master.Protein.Accessions), ]
rmat<-rmat[!duplicated(rownames(rmat)), ]

ppm<-tevm$value
pm<-tevm %>% group_by(Master.Protein.Accessions, rq) %>% summarise(pm = sum(value))
pm<-pm$pm
rm<-c(rmat)
rm[rm%in%c(0)]<-NA

ppm<-c(ppm, rep(NA, length(rmat)-length(ppm))); length(ppm)
pm<-c(pm, rep(NA, length(rmat)-length(pm))); length(pm)

df<-data.frame(rm, pm, ppm)

# Convert S/N to number of ions on QE:
noise<- 3.5*sqrt(240000 / 70000)
df[,2:3]<-df[,2:3]*noise

# Rearrange data for convenience
dfm<-melt(df)

# Calculate Poisson error:
pep.error<-round( 1 / sqrt( median ( df$ppm , na.rm=T)),2)*100
prot.error<-round( 1 / sqrt( median ( df$pm , na.rm=T)),2)*100
rna.error<-round( 1 / sqrt( median ( df$rm , na.rm=T)),2)*100

# Log10 transform data
dfm$valuelog10<-log10(dfm$value)


library(scales)

pois_error<-function(x){ 1/sqrt(x) }

x1<-log10(seq(1,10000,1))
y1<-pois_error(10^(x1))

dat<-data.frame(x1,y1); colnames(dat)<-c("x","y")



# Visualize
p1<-ggplot(dfm, aes(x = valuelog10, y = variable, fill=variable)) +
  geom_density_ridges() + theme_ridges(center_axis_labels = TRUE) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(-0.1, 0)) +
  rremove("legend") +
  xlab("\ Copy number measured\n per single cell") +
  scale_fill_manual(values = c("blue", "red", "orange"))+
  scale_x_continuous(limits = c(0, 4.6), breaks = c(0,1,2,3,4),
                     labels = c(bquote(10^0),bquote(10^1),bquote(10^2),bquote(10^3),bquote(10^4) )) +
  annotate("text", x=3.5, y =2.5, label=paste0(length(unique(tevm$Master.Protein.Accessions)), " Proteins"), color="red", size=9) +
  annotate("text", x=3.5, y =1.5, label=paste0(length(rownames(rmat)), " mRNAs"), color="blue", size=9)+
  annotate("text", x=3.5, y =3.6, label=paste0(length(unique(tevm$Annotated.Sequence)), " Peptides"), color="orange3", size=9)+
  ylab("Density") +
  theme(text=element_text(size=15)) +
  theme(text = element_text(size = 25)) + rremove("y.text") + rremove("y.ticks")+
  #ggtitle("Sampling molecules in single cells\n") +
  font("title", size=20) +
  font("x.text", size= 30)  +
  font("ylab",size= 30)


p2<-ggplot(dat, aes(x,y)) + geom_line(size=1, color="black") +
  scale_x_continuous(limits = c(0, 4.6), breaks = c(0,1,2,3,4)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0,50,100)) +
  ylab("CV,  %") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix="")) +
  theme_ridges(center_axis_labels = TRUE) +
  #theme(panel.grid.major.y = element_blank())+
  rremove("x.text") +
  rremove("xlab") +
  geom_point(aes(x=log10( median ( df$rm , na.rm=T) ), y=rna.error/100), colour="blue", size=7) +
  geom_point(aes(x=log10( median ( df$pm , na.rm=T) ), y=prot.error/100), colour="red", size=7) +
  geom_point(aes(x=log10( median ( df$ppm , na.rm=T) ), y=pep.error/100), colour="orange", size=7) +
  annotate("text",x=3.5, y=0.45, label = expression(paste(CV==frac(1, sqrt(n)))), size=8) +
  annotate("text",x=3.5, y=0.9, label = "Counting error:", size=9) +
  font("title", size=20) +
  font("ylab",size=28) +
font("y.text",size=25)
#theme(axis.line.x = element_line(color="black"))

px<-plot_spacer() + p2 +plot_spacer()+ p1 + plot_layout(ncol = 1, heights = c(0.1,2,0.2, 5))

ggsave("figs/ion_count.pdf", plot=px, height = 8, width=6.5)




# Density correlation: Protein and mRNA -----------------------------------
ct1<-(rm.batch[,colnames(rm.batch)%in%monolike])
ct2<-(rm.batch[,colnames(rm.batch)%in%maclike])

# ct1<-(rm.f3[,!colnames(rm.f3)%in%maclike])
# ct2<-(rm.f3[,colnames(rm.f3)%in%maclike])

c1<-c(cor(ct1, use="pairwise.complete.obs", method="spearman"))
c2<-c(cor(ct2, use="pairwise.complete.obs", method="spearman"))

# Combine into data frame:
c2<-c(c2, rep(NA, length(c1)-length(c2)))
dfxxx<-data.frame(c1, c2)
dxm<-melt(dfxxx)

# Grab cells
m0.ind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_m0"]
u.ind<-ev.melt.uniqueID$id[ev.melt.uniqueID$celltype=="sc_u"]

ct1<-(matrix.sc.batch[,colnames(matrix.sc.batch)%in%u.ind])
ct2<-(matrix.sc.batch[,colnames(matrix.sc.batch)%in%m0.ind])

# ct1<-(ev.matrix.sc.f.n[,colnames(ev.matrix.sc.f.n)%in%u.ind])
# ct2<-(ev.matrix.sc.f.n[,colnames(ev.matrix.sc.f.n)%in%m0.ind])

c1<-c(cor(ct1, use="pairwise.complete.obs", method="spearman"))
c2<-c(cor(ct2, use="pairwise.complete.obs", method="spearman"))

mean(c1); mean(c2)

# Combine into data frame:
c1<-c(c1, rep(NA, length(c2)-length(c1)))
dfxxx<-data.frame(c1, c2)
dxm2<-melt(dfxxx)

dxm$dataset<-"RNA"
dxm2$dataset<-"Protein"

dxm3<-rbind(dxm, dxm2)


px<-ggplot(data = dxm3, aes(y=dataset, x= value, alpha=0.5)) +
  geom_density_ridges(aes(fill=variable, scale=0.75)) +
  scale_fill_manual(values=my_colors[c(3,2)]) +
  xlim(c(-0.1,0.65)) +
  theme_pubr() +
  rremove("legend") +
  rremove("y.ticks") +
  rremove("y.text") +
  ylab("Density") +
  xlab("Spearman correlation between cells") +
  ggtitle("Similarity within cell type") +
  font("title", size=40) +
  font("xylab", size=30) +
  font("x.text", size=30) +
  annotate("text", x=-0.1, y=1.5, label="Protein", size=15) +
  annotate("text", x=-0.1, y=2.5, label="mRNA", size=15) +
  annotate("text", x=0.5, y=2.45, label="Macrophage", size=10, color=my_colors[2]) +
  annotate("text", x=0.5, y=2.25, label="Monocyte", size=10, color=my_colors[3]) +
  theme(panel.grid.minor = element_line(colour="gray", size=0.5))

ggsave("figs/density_correlation.pdf", plot=px, width=6, height=5)


# Export raw data------------------------------------------------------------------

load("dat/raw.RData")

# Group the experimental runs by type: single cell, 10 cell, 100 cell, 1000 cell, or ladder
all.runs<-unique(ev$Raw.file)
c10.runs<-c( all.runs[grep("col19", all.runs)] , all.runs[grep("col20", all.runs)] ); length(unique(c10.runs))
c100.runs<-c( all.runs[grep("col21", all.runs)] , all.runs[grep("col22", all.runs)] ); length(unique(c100.runs))
c1000.runs<-c( all.runs[grep("col23", all.runs)] ); length(unique(c1000.runs))
ladder.runs<-c( all.runs[grep("col24", all.runs)] ); length(unique(ladder.runs))
carrier.runs<-c( all.runs[grep("arrier", all.runs)] ); length(unique(carrier.runs))
sc.runs<-all.runs[!all.runs%in%c(c10.runs,c1000.runs,c100.runs,ladder.runs, carrier.runs)]; length(unique(sc.runs))

write.csv(ev[ev$Raw.file%in%sc.runs, ], "dat/ev_SCoPE2_runs.csv")
write.csv(ev[!ev$Raw.file%in%sc.runs, ], "dat/ev_additional_runs.csv")


# Export processed data ---------------------------------------------------

load("dat/imported.RData")

# Peptides-raw.csv
peptides.raw<-data.frame(t3)
peptides.raw$peptide<-rownames(t3)
peptides.raw$protein<-ev.melt.pep$protein[match(rownames(t3), ev.melt.pep$sequence)]

write.csv(peptides.raw[,c(ncol(peptides.raw),ncol(peptides.raw)-1 , 1:(ncol(peptides.raw)-2))], "dat/Peptides-raw.csv", row.names = F)

# Proteins-processed.csv
prot.pro<-data.frame(matrix.sc.batch)
prot.pro$protein<-rownames(matrix.sc.batch)
write.csv(prot.pro, "dat/Proteins-processed.csv")

# Cells.csv
cells<-ev.melt.uniqueID[ev.melt.uniqueID$id%in%colnames(ev.matrix.sc.f.n), c("id","celltype","digest","sortday","lcbatch","Raw.file")]
row.names(cells) <- NULL
colnames(cells)<-c("id","celltype","batch_digest","batch_sort","batch_chromatography","raw.file")
cells.t<-t(cells)
colnames(cells.t)<-cells$id
cells.t<-cells.t[-1,]

write.csv(cells.t, "dat/Cells.csv", row.names = T)
