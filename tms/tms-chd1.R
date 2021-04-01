library(R6)
library(TrenaProjectErythropoiesis)
library(TrenaMultiScore)
library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
TMS = R6Class("TMS",

    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   tp=NULL,
                   tms=NULL,
                   tbl.gh=NULL,
                   ghRegion=NULL,
                   tbl.fimo=NULL,
                   tbl.rbp=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(targetGene, trenaProject){
           suppressWarnings(db.access.test <-
                                try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
           if(length(db.access.test) == 0)
              stop("khaleesi database server unavailable")
           printf("initializing TMS('%s')", targetGene)
           private$targetGene <- targetGene
           private$tp <- trenaProject
           private$tms <- TrenaMultiScore(private$tp, targetGene)
           },
        addGeneHancer = function(){
           private$tbl.gh <- getEnhancers(private$tp, tissues="all", maxSize=20000)
           private$ghRegion <- getGeneHancerRegion(private$tms)
           },
        getGeneHancer = function(){
           private$tbl.gh
           },
        getGeneHancerRegion = function(){
           private$ghRegion
           },
        addOpenChromatin = function(){
           invisible(findOpenChromatin(private$tms))
           },
        getOpenChromatin = function(){
           getOpenChromatin(private$tms)
           },
        addAndScoreFimoTFBS = function(fimo.threshold=1e-4){
           findFimoTFBS(private$tms, motifs=list(), fimo.threshold=fimo.threshold)
           tbl.tmp <- getMultiScoreTable(private$tms)
           if(nrow(tbl.tmp) == 0){
               s <- sprintf("no motif hits at threshold %f", fimo.threshold)
               stop(s)
               }
           addChIP(private$tms)
           scoreMotifHitsForConservation(private$tms)
           scoreMotifHitsForGeneHancer(private$tms)
           addGenicAnnotations(private$tms)
           addDistanceToTSS(private$tms)
           },
        addChIP = function(){
           addChIP(private$tms)
           },
        getTfTable = function(){
           getMultiScoreTable(private$tms)
           },

        addRBP = function(){
           db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
             #  AND celltype='K562'
           query <- sprintf("select * from rbp where target='%s'", private$targetGene)
           tbl.rbp <- dbGetQuery(db, query)
           colnames(tbl.rbp)[grep("endpos", colnames(tbl.rbp))] <- "end"
           private$tbl.rbp <- tbl.rbp
           dbDisconnect(db)
           },

        getRbpTable = function(){
            private$tbl.rbp
            },

        add.tf.mrna.correlations = function(mtx.mrna, featureName){
           addGeneExpressionCorrelations(private$tms, mtx.mrna, featureName)
           },

        add.rbp.mrna.correlations = function(mtx.mrna, featureName){
           f <- function(rbp){
              if(rbp %in% rownames(mtx))
                 return(cor(mtx.mrna[private$targetGene,], mtx.mrna[rbp,], method="spearman"))
              else return(NA)
              }
           suppressWarnings(cor <- unlist(lapply(private$tbl.rbp$gene, f)))
           private$tbl.rbp[, featureName] <- cor
           }

        #------------------------------------------------------------
       ) # public
    ) # class

#--------------------------------------------------------------------------------
