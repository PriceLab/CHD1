library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(GeneHancerDB)
library(trena)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tmse"))
   tpe <- TrenaProjectErythropoiesis()

if(!exists("igv"))
    igv <- start.igv("CHD1", "hg38")

tms <- TrenaMultiScore(tpe, "CHD1")

tbl.ghFullRegion <- getGeneHancerRegion(tms)    # 240k

ghdb <- GeneHancerDB()
tbl.gh <- retrieveEnhancersFromDatabase(ghdb, "CHD1", tissues="all")
dim(tbl.gh)
tbl.gh$score <- asinh(tbl.gh$combinedscore)
track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "score")], color="random", autoscale=TRUE)
displayTrack(igv, track)

tbl.atac.merged <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
tbl.atac.chd1 <- subset(tbl.atac.merged, chrom==tbl.ghFullRegion$chrom &
                                         start >= tbl.ghFullRegion$start &
                                         end <= tbl.ghFullRegion$end)
dim(tbl.atac.chd1)
track <- DataFrameAnnotationTrack("atac", tbl.atac.chd1, color="gray")
displayTrack(igv, track)

with(tbl.gh, findOpenChromatin(tmse, chrom, start, end))
tbl.oc <- getOpenChromatin(tmse)

    #------------------------------------------------
    # build a tf-only model in oc in gh-ish promoter
    #------------------------------------------------

mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
dim(mtx)

    #--------------------------------------------------
    # identify the genes with high correlation to CHD1
    #--------------------------------------------------

cors <- lapply(rownames(mtx), function(gene) cor(mtx["CHD1",], mtx[gene,], method="spearman"))
names(cors) <- rownames(mtx)
high.cors <- which(abs(as.numeric(cors)) > 0.7)
high.cor.genes <- names(cors)[high.cors]
length(high.cor.genes)  # 5300
head(high.cor.genes)
mdb.indices <- as.data.frame(findMatches(high.cor.genes,  mcols(MotifDb)$geneSymbol))$subjectHits
length(mdb.indices)     # 1035
tf.candidates <- mcols(MotifDb)$geneSymbol[mdb.indices]
length(tf.candidates)
length(unique(tf.candidates))  # 288

mdb.filtered <- MotifDb[mdb.indices]

tbl.big.atac.promoter <- data.frame(chrom="chr5", start=98928153, end=98929885, stringsAsFactors=FALSE)
with(tbl.big.atac.promoter, 1 + end - start)  # 1733 bases

source("~/github/fimoService/batchMode/fimoBatchTools.R")
meme.filename <- "chd1.atac.promoter.meme"
motifs <- query(mdb.filtered, c("sapiens"), c("jaspar2018", "hocomocov11-core"))
length(motifs) # 168
export(motifs, con=meme.filename, format="meme")
tbl.fimo <- fimoBatch(tbl.big.atac.promoter, matchThreshold=1e-6, genomeName="hg38", pwmFile=meme.filename)
dim(tbl.fimo)  # 37 9
tbl.freq <- as.data.frame(sort(table(tbl.fimo$tf), decreasing=TRUE))
       #--------------------------------------------------
       # fimo threshold of 1e-6
       #     Var1 Freq
       # 1  PATZ1   12
       # 2    SP3    9
       # 3    MAZ    8
       # 4 ZNF263    4
       # 5  KLF12    1
       # 6  NFKB1    1
       # 7  STAT3    1
       # 8  ZNF76    1
       #--------------------------------------------------
tfs <- as.character(tbl.freq$Var1)
for(tf.name in tfs){
   tbl.tf <- subset(tbl.fimo, tf==tf.name)
   track <- DataFrameAnnotationTrack(tf.name, tbl.tf, color="blue", trackHeight=25)
   displayTrack(igv, track)
   }

solver <- EnsembleSolver(mtx,
                         targetGene="CHD1",
                         candidateRegulators=tfs,
                         geneCutoff=1.0,
                         solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"))
tbl.out <- run(solver)

    #----------------------------------------------------------------------------------
    #  tbl.out[order(tbl.out$rfScore, decreasing=TRUE),]
    #     gene  betaLasso  betaRidge spearmanCoeff pearsonCoeff   rfScore     xgboost
    # 3  NFKB1 0.00000000 0.04542964     0.9458128    0.7823592 0.7729057 0.063811469
    # 8  ZNF76 0.08004888 0.12219364     0.9161079    0.8008851 0.6645926 0.153028981
    # 7 ZNF263 0.10864989 0.08761415     0.9293924    0.8407644 0.6084442 0.104255861
    # 5    SP3 1.02151166 0.44297837     0.8845101    0.8838637 0.5822912 0.272322825
    # 1  KLF12 0.03519295 0.02486689     0.7137773    0.7398377 0.5817719 0.376467555
    # 6  STAT3 0.00000000 0.05885936     0.7914614    0.7911963 0.5618033 0.001947006
    # 4  PATZ1 0.00000000 0.01763782     0.7733990    0.6811361 0.4789502 0.001945737
    # 2    MAZ 0.00000000 0.07077797     0.7733990    0.6904603 0.3217621 0.026220566
    #----------------------------------------------------------------------------------

summary(lm(CHD1 ~ 1 + ., data=as.data.frame(t(mtx[c(tfs, "CHD1"),]))))

    #-------------------------------------------------------------------
    # Coefficients:
    #             Estimate Std. Error t value Pr(>|t|)
    # (Intercept) -4.62038    1.42172  -3.250 0.004215 **
    # PATZ1       -0.74506    0.18091  -4.118 0.000585 ***
    # SP3          0.83512    0.19331   4.320 0.000369 ***
    # MAZ          0.60299    0.28803   2.093 0.049955 *
    # ZNF263       0.36739    0.20787   1.767 0.093211 .
    # KLF12        0.09188    0.01905   4.823 0.000118 ***
    # NFKB1        0.55814    0.20937   2.666 0.015274 *
    # STAT3       -0.15436    0.12837  -1.202 0.243976
    # ZNF76       -0.15133    0.21010  -0.720 0.480124
    # ---
    # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    #
    # Residual standard error: 0.09977 on 19 degrees of freedom
    # Multiple R-squared:  0.9614,	Adjusted R-squared:  0.9451
    # F-statistic: 59.13 on 8 and 19 DF,  p-value: 8.481e-12
    #-------------------------------------------------------------------


    #----------------------------------------
    # now add in the rna binding proteins
    #----------------------------------------

library(RPostgreSQL)
db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")

query <- "select * from rbp where target='CHD1'"   # AND celltype='K562'  # multipotent hematopoietic cells
tbl <- dbGetQuery(db, query)
dim(tbl)  # 7562 12
length(unique(tbl$gene)) # 120 unique RBP

query <- "select * from rbp where target='CHD1' AND celltype='K562'"  # multipotent hematopoietic cells
tbl <- dbGetQuery(db, query)
dim(tbl)  # 343 12
rbps <- unique(tbl$gene)
length(rbps) # 54

length(intersect(rbps, rownames(mtx))) # 54

solver <- EnsembleSolver(mtx,
                         targetGene="CHD1",
                         candidateRegulators=c(tfs, rbps),
                         geneCutoff=0.5,
                         solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"))
tbl.out <- run(solver)
dim(tbl.out)
tbl.out <- tbl.out[order(tbl.out$rfScore, decreasing=TRUE),]
head(tbl.out, n=20)

     #--------------------------------------------------------------------------------------
     #      gene  betaLasso    betaRidge spearmanCoeff pearsonCoeff    rfScore      xgboost
     # 15  NFKB1 0.00000000  0.018089433     0.9458128   0.78235916 0.26418898 5.464711e-02
     # 27  TRA2A 0.39526691  0.054860516     0.9107827   0.93226637 0.25642550 0.000000e+00
     # 29  ZNF76 0.00000000  0.032139407     0.9161079   0.80088513 0.25448943 1.347558e-01
     # 28 ZNF263 0.00000000  0.028079517     0.9293924   0.84076438 0.22538513 1.051584e-01
     # 25    SP3 0.45078447  0.149073486     0.8845101   0.88386373 0.17172273 2.065055e-01
     # 22  SF3B4 0.00000000  0.016982726     0.7794198   0.69898888 0.16395890 0.000000e+00
     # 24 SMNDC1 0.00000000  0.034883192     0.8571429   0.80365231 0.15603043 6.353706e-04
     # 21  SF3B1 0.01786503  0.074917718     0.8795840   0.88989950 0.15460884 3.901526e-05
     # 19   RPS3 0.04320629  0.042567141     0.8817734   0.88537142 0.14081685 6.048630e-02
     # 23   SLTM 0.00000000  0.037362655     0.8883415   0.82007208 0.12778471 0.000000e+00
     # 18  RBM22 0.00000000  0.042865847     0.8861522   0.87747312 0.12035724 2.263045e-02
     # 26  STAT3 0.00000000  0.018942204     0.7914614   0.79119632 0.11335231 0.000000e+00
     # 7  GTF2F1 0.00000000  0.026038723     0.8664477   0.80196640 0.11334644 0.000000e+00
     # 20   SBDS 0.11711038  0.058302482     0.7903667   0.88729716 0.11209869 7.167678e-06
     # 17   PUM2 0.00000000 -0.010960838     0.5747126   0.28297437 0.10286898 7.162598e-03
     # 10   ILF3 0.00000000  0.027283590     0.8385331   0.79302823 0.09768723 0.000000e+00
     # 1   DDX3X 0.00000000  0.094080218     0.7958402   0.78812121 0.09383062 1.144630e-02
     # 12  LARP4 0.00000000 -0.029454532     0.5298303  -0.09790866 0.08194445 6.098436e-03
     # 13    MAZ 0.00000000  0.024664978     0.7733990   0.69046026 0.06101670 1.223065e-02
     # 11  KLF12 0.00000000  0.006305911     0.7137773   0.73983770 0.06052956 3.253615e-01
     #--------------------------------------------------------------------------------------


summary(lm(CHD1 ~ 1 + ., data=as.data.frame(t(mtx[unique(c(tfs, head(tbl.out$gene, n=20), "CHD1")),]))))

     #-------------------------------------------------------------
     # Residuals:
     #       Min        1Q    Median        3Q       Max
     # -0.055377 -0.008679  0.003791  0.013276  0.033927
     #
     # Coefficients:
     #              Estimate Std. Error t value Pr(>|t|)
     # (Intercept)  4.542851   2.056374   2.209  0.06922 .
     # PATZ1       -0.554474   0.210069  -2.639  0.03857 *
     # SP3          1.332820   0.371862   3.584  0.01158 *
     # MAZ          0.855970   0.352949   2.425  0.05150 .
     # ZNF263       1.886234   0.555905   3.393  0.01462 *
     # KLF12        0.018279   0.021306   0.858  0.42387
     # NFKB1        0.177954   0.220481   0.807  0.45040
     # STAT3        0.238749   0.194260   1.229  0.26508
     # ZNF76        0.769852   0.322760   2.385  0.05438 .
     # TRA2A        0.168729   0.267305   0.631  0.55117
     # SF3B4        0.103942   0.387834   0.268  0.79767
     # SMNDC1      -0.394985   0.266998  -1.479  0.18953
     # SF3B1        0.007952   0.379437   0.021  0.98396
     # RPS3        -1.040891   0.204971  -5.078  0.00227 **  # ribosome protein, and also endonuclease
     # SLTM         0.221743   0.257923   0.860  0.42296
     # RBM22       -0.011041   0.209952  -0.053  0.95977
     # GTF2F1      -0.776679   0.368270  -2.109  0.07948 .
     # SBDS         0.066047   0.290934   0.227  0.82795
     # PUM2         0.228439   0.444470   0.514  0.62565
     # ILF3        -0.872259   0.468243  -1.863  0.11178
     # DDX3X       -0.746716   0.388909  -1.920  0.10326
     # LARP4       -0.503798   0.223484  -2.254  0.06506 .
     # ---
     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
     #
     # Residual standard error: 0.03891 on 6 degrees of freedom
     # Multiple R-squared:  0.9981,	Adjusted R-squared:  0.9917
     # F-statistic: 153.8 on 21 and 6 DF,  p-value: 1.583e-06
     #-------------------------------------------------------------

print(load("~/github/TrenaProjectErythropoiesis/viz/srm.vs.mrna/shinyapps.io/srm.rna.averaged.clean.RData"))
# mtx.rna, mtx.srm
dim(mtx.rna) #
dim(mtx.srm)
lines(mtx.srm["CHD1",]/8000, type="b", col="blue", ylim=c(0, 1))
plot(mtx["RBM22", seq(1, 24, by=2)]-8.6, type="b", col="red")
legend(2, 0.1, c("srm CHD1", "mRNA RBM22"), c("blue", "red"))

system("open 'chd1-proteinLevels-RBM22-mrnaLevels.png'")
system("open'CHD1-promoter-igv-with-GH-and-8-tf-binding-sites.png')"
