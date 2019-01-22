

my.dmp.finder <- function(betas, chip = c("450k", "EPIC"),
                          group.pos, group.neg,
                          rm.na = TRUE, offset = 100) {
  require(qvalue)
  chip <- match.arg(chip)
  if (chip == "450k") {
    print("Using 450k chip annotation")
    data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19
  }
  if (chip == "EPIC") {
    print("Using EPIC chip annotation")
    data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    anno <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19
  } 
  locs <- getLocations(anno, mergeManifest=FALSE, orderByLocation=TRUE)    
return(locs)
  locs$str <- anno$strand
  locs$dist.to.next <- distance(locs, c(locs[-1], locs[1]))
  locs$dist.to.prev <- distance(
                         locs,
                         c(locs[length(locs)], locs[-length(locs)])
                       )
  locs$singleton <- (locs$dist.to.next > 500 | is.na(locs$dist.to.next)) &
                    (locs$dist.to.prev > 500 | is.na(locs$dist.to.prev))
  #mcols(locs) <- cbind(mcols(locs), anno[,c(19,18,24,26,25,32)])
  t.list <- lapply(1:nrow(betas), function(r) {
                                    my.t <- t.test(
                                              x = betas[r,group.neg],
                                              y = betas[r,group.pos],
                                              alternative = "two.sided",
                                              mu = 0,
                                              paired = FALSE,
                                              var.equal = FALSE,
                                              conf.level = 0.95
                                            )
                                    means <- unname(my.t$estimate)
                                    c(
                                      group.pos.mean=means[2],
                                      group.neg.means=means[1],
                                      abs.mean.diff=abs(means[2]-means[1]),
                                      mean.diff=means[2]-means[1],
                                      p.value=my.t$p.value
                                    )
                                    })
  t.df<- as.data.frame(do.call("rbind", t.list))
  t.df$q.value <- qvalue(t.df[,"p.value"])$qvalues
  mcols(locs) <- cbind(mcols(locs), t.df)
  locs
}


