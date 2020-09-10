source("code/functions_parameters.R")

# Import ------------------------------------------------------------------

# Import S/N for all scans

sn_files<-list.files("/Users/hs/GoogleDrive/My\ Drive/MS/SCoPE/public_uploads/sn_scope2/", pattern="sn.txt_")
sn<-read.delim(paste0("/Users/hs/GoogleDrive/My\ Drive/MS/SCoPE/public_uploads/sn_scope2/", sn_files[1]) )

for (X in sn_files[-1]){
  
  snt<-read.delim(paste0("/Users/hs/GoogleDrive/My\ Drive/MS/SCoPE/public_uploads/sn_scope2/", X) )
  sn<-rbind(sn, snt)
  
}

sn$raw_file<-paste0("X", sn$raw_file)

# Load raw data and annotations

#ev<-read_csv("dat/ev_tmt11_tmt16_2.csv")
ev<-read_csv("dat/evidence_unfiltered_v2.csv")

ev$scan<-paste0(ev$Raw.file,"_",ev$MS.MS.scan.number)
sn$scan<-paste0(sn$raw_file,"_",sn$scan_number)

ev<-merge(ev, sn, by="scan")

ev<-ev[ev$Raw.file%in%sn$raw_file, ]

parse_row<-grep("|",ev$Leading.razor.protein, fixed=T)

if(length(parse_row)>0){
  split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
  split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  
  ev$Leading.razor.protein[parse_row]<-split_prot2
}
#grep("|",ev$Leading.razor.protein, fixed=T)

design<-read.csv("dat/annotation.csv")

batch<-read.csv("dat/batch.csv")

# Attach batch data to protein data
ev[,colnames(batch)[-1]]<-NA
for(X in batch$set){

  ev$lcbatch[ev$Raw.file==X] <- as.character(batch$lcbatch[batch$set%in%X])
  ev$sortday[ev$Raw.file==X] <- as.character(batch$sortday[batch$set%in%X])
  ev$digest[ev$Raw.file==X] <- as.character(batch$digest[batch$set%in%X])

}


#save(ev,batch,design,file="dat/raw.RData")
#load("dat/raw.RData")
#write.csv(paste0(unique(ev$Raw.file),".raw | "), "all_raw_files_used_besides_ladder.csv")

# Create unique peptide+charge column:
ev$modseq<-paste0(ev$Modified.sequence,ev$Charge)

# If somehow we have NA for raw files?! Remove 'em
ev<-ev[!is.na(ev$Raw.file),]

# Add X in front of experiment names because R doesn't like column names starting with numbers
#ev$Raw.file<-paste0("X",ev$Raw.file)
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
#ri.index<-which(colnames(ev)%in%paste0("Reporter.intensity.",1:16))
sn_names<-c("X126_sn",
            "X127N_sn","X127C_sn",
            "X128N_sn","X128C_sn",
            "X129N_sn","X129C_sn",
            "X130N_sn","X130C_sn",
            "X131_sn","X131C_sn",
            "X132N_sn","X132C_sn",
            "X133N_sn","X133C_sn",
            "X134_sn")
ri.index<-which(colnames(ev)%in%sn_names)

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
# ev[ev$Raw.file%in%sc.runs, ri.index] <- ev[ev$Raw.file%in%sc.runs, ri.index] / ev[ev$Raw.file%in%sc.runs, ri.index[2]]

# if(ref_demo){
#   
#   ev[ev$Raw.file%in%sc.runs, ri.index] <- ev[ev$Raw.file%in%sc.runs, ri.index] / rowMeans(ev[ev$Raw.file%in%sc.runs, ri.index], na.rm=T)
#   
#   
# }

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

load("dat/rna.RData")
tev<-ev.melt[ev.melt$celltype%in%c("sc_m0","sc_u"),]

rmat<-rmx
gn_uni<-read.csv("dat/gn_uni.csv")
gn_uni.f<-gn_uni[gn_uni$uniprot%in%tev$protein,]
row.names(rmat)<-gn_uni.f$uniprot[match(rownames(rmat), gn_uni.f$gn)]
rmat<-rmat[rownames(rmat)%in%(tev$protein), ]
dim(rmat)


# Take only proteins seen in mRNA data and final protein data
tev<-tev[tev$protein%in%rownames(rmat), ]

rmat<-rmat[rownames(rmat)%in%(tev$protein), ]
rmat<-rmat[!duplicated(rownames(rmat)), ]

tev$quantitation<-as.numeric(tev$quantitation)

ppm<-tev$quantitation

#library(dplyr)
pm<-tev %>% group_by(protein, Raw.file) %>% summarise(pm = sum(quantitation))
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
  annotate("text", x=3.5, y =2.6, label=paste0(length(unique(tev$protein)), " Proteins"), color="red", size=7) +
  annotate("text", x=3.5, y =1.5, label=paste0(length(rownames(rmat)), " mRNAs"), color="blue", size=7)+
  annotate("text", x=3.5, y =3.3, label=paste0(length(unique(tev$sequence)), " Peptides"), color="orange3", size=7)+
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
px

ggsave("ion_count.pdf", plot=px, height = 8, width=6.5)

