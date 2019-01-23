
#p.q.mean.spans <- function(p.q.vals, abs.mean.diff, p.q.thresh, mean.diff.thresh, span, promote.singles){
#
#    below <- p.q.vals <= p.q.thresh & abs.mean.diff >= mean.diff.thresh
#	if(promote.singles & length(below) > 2) below <- below | (c(1,below[1:(length(below)-1)]) & c(below[2:length(below)],1))
#	below <- c(FALSE, below, FALSE)
#	d <- diff(below)
#	rising <- which(d == 1)
#	falling <- which(d == -1)
#	width <- falling - rising
#	sat.span <- width >= span
#	start.pos <- rising[sat.span]
#	end.pos <- falling[sat.span]-1
#	cbind(start.pos, end.pos, width=end.pos-start.pos+1)
#}
#
#p.q.mean.spans2 <- function(p.q.vals, abs.mean.diff, singleton, p.q.thresh, mean.diff.thresh, span, promote.singles){
#    below <- p.q.vals <= p.q.thresh & abs.mean.diff >= mean.diff.thresh
#    #if (length(below) == 0){
#    #    return(NULL)
#    #}
#	if(promote.singles & length(below) > 2) below <- below | (c(1,below[1:(length(below)-1)]) & c(below[2:length(below)],1))
#    below <- below & !singleton
#	below <- c(FALSE, below, FALSE)
#	d <- diff(below)
#	rising <- which(d == 1)
#	falling <- which(d == -1)
#	width <- falling - rising
#	sat.span <- width >= span
#	start.pos <- rising[sat.span]
#	end.pos <- falling[sat.span]-1
#	cbind(start.pos, end.pos, width=end.pos-start.pos+1)
#}
#
#p.spans <- function(x, threshold = 0.05, span = 2, allow.singles = TRUE){
#
#	below <- x <= threshold
#	if(!allow.singles) below <- below | (c(1,below[1:(length(below)-1)]) & c(below[2:length(below)],1))
#	below <- c(FALSE, below, FALSE)
#	d <- diff(below)
#	rising <- which(d == 1)
#	falling <- which(d == -1)
#	width <- falling - rising
#	sat.span <- width >= span
#	start.pos <- rising[sat.span]
#	end.pos <- falling[sat.span]-1
#	cbind(start.pos, end.pos, width=end.pos-start.pos+1)
#}
#
#p.mean.spans <- function(t.table, p.thresh=0.001, mean.diff.thresh=0.3, span = 3, promote.singles = TRUE){
#
#    pvals <- t.table$p.value
#    mean.diff <- t.table$abs.mean.diff
#	below <- pvals <= p.thresh & mean.diff >= mean.diff.thresh
#	if(promote.singles) below <- below | (c(1,below[1:(length(below)-1)]) & c(below[2:length(below)],1))
#	below <- c(FALSE, below, FALSE)
#	d <- diff(below)
#	rising <- which(d == 1)
#	falling <- which(d == -1)
#	width <- falling - rising
#	sat.span <- width >= span
#	start.pos <- rising[sat.span]
#	end.pos <- falling[sat.span]-1
#	cbind(start.pos, end.pos, width=end.pos-start.pos+1)
#}
#
#q.mean.spans <- function(t.table, q.thresh=0.01, mean.diff.thresh=0.3, span = 2, promote.singles = TRUE){
#    qvals <- t.table$q.value
#    mean.diff <- t.table$abs.mean.diff
#	below <- qvals <= q.thresh & mean.diff >= mean.diff.thresh
#	if(promote.singles) below <- below | (c(1,below[1:(length(below)-1)]) & c(below[2:length(below)],1))
#	below <- c(FALSE, below, FALSE)
#	d <- diff(below)
#	rising <- which(d == 1)
#	falling <- which(d == -1)
#	width <- falling - rising
#	sat.span <- width >= span
#	start.pos <- rising[sat.span]
#	end.pos <- falling[sat.span]-1
#	cbind(start.pos, end.pos, width=end.pos-start.pos+1)
#}
#
#t.test.list2table <- function(tList){
#    require(qvalue)
#	m <- sapply(tList, function(test) {
#	                                   means <- unname(test$estimate)
#									   c(x.mean=means[1], y.means=means[2],
#									     abs.mean.diff=abs(means[2]-means[1]),
#										 mean.diff=means[2]-means[1],
#										 p.value=test$p.value)
#										})
#	m <- as.data.frame(t(m))
#	m$q.value <- qvalue(m$p.value)$qvalues
#	m
#}



dmr_finder <- function(dmp.table,
                       use.value.dmr=c("p", "q"),
                       p.q.thresh.dmr=0.01,
                       mean.diff.thresh.dmr=0.3,
                       span=3,
                       promote.singles=TRUE,
                       find.singletons=TRUE,
                       use.value.sing=NULL, p.q.thresh.sing=NULL,
                       mean.diff.thresh.sing=NULL, merge.tables=FALSE) {
    require(GenomicRanges)
    use.value.dmr <- match.arg(use.value.dmr)
    if (is.null(use.value.sing)) {
      use.value.sing <- use.value.dmr
    } else {
      use.value.sing <- match.arg(use.value.dmr)
    }
    if (is.null(p.q.thresh.sing)) {
      p.q.thresh.sing <- p.q.thresh.dmr
    }
    if (is.null(mean.diff.thresh.sing)) {
      mean.diff.thresh.sing <- mean.diff.thresh.dmr
    }
    if (!find.singletons) {
      merge.tables <- FALSE
    }
    #mcols.slice <- c(20:21,1:4,8:13,17:19,22:23)
    #mcols.slice <- c(20:21,1:4,8:13,17:19,22:23,7)
    #dmp.glist <- split(dmp.table, factor(as.vector(seqnames(dmp.table))))
    dmp.glist <- split(dmp.table, seqnames(dmp.table))
    #return(dmp.glist)
    dmr.list <- lapply(dmp.glist,
                  function(x) {
                    if(use.value.dmr == "p") {
                      spans <- p.q.mean.spans2(
                                 x$p.value,
                                 x$abs.mean.diff,
                                 x$singleton,
                                 p.q.thresh.dmr,
                                 mean.diff.thresh.dmr,
                                 span, promote.singles
                               )
                    } else {
                      spans <- p.q.mean.spans2(
                                 x$q.value,
                                 x$abs.mean.diff,
                                 x$singleton,
                                 p.q.thresh.dmr,
                                 mean.diff.thresh.dmr,
                                 span, promote.singles
                               )
                    }
                    if(nrow(spans) > 0) {
                      do.call("c",
                        lapply(1:nrow(spans),
                          function(y) {
                            r <- spans[y,]
                            g.slice <- x[r[1]:r[2],]
                            range.w.prev <- c(0, diff(g.slice$pos)-1)
                            dmr.range <- g.slice$pos[length(g.slice$pos)] - g.slice$pos[1] - 1
                            g.slice$window.num <- y
                            g.slice$window.size <- r[3]
                            g.slice$width.to.prev <- range.w.prev
                            g.slice$dmr.width <- dmr.range
                            g.slice[, mcols.slice]
                          }))
                    }
                  })

    if(find.singletons) { 
        singleton.list <- lapply(dmp.glist, function(xx) {
                                            if(use.value.sing == "p") {
                                                sig.singletons <- xx$singleton & xx$p.value <= p.q.thresh.sing & xx$abs.mean.diff >= mean.diff.thresh.sing
                                            } else {
                                                sig.singletons <- xx$singleton & xx$q.value <= p.q.thresh.sing & xx$abs.mean.diff >= mean.diff.thresh.sing
                                            }
                                            sig.singletons.spans <- cbind(start=which(sig.singletons), end=which(sig.singletons), width=rep(1,sum(sig.singletons)))
                                            if(nrow(sig.singletons.spans) > 0) {
                                                do.call("c", lapply(1:nrow(sig.singletons.spans), function(yy) {
                                                                    r <- sig.singletons.spans[yy,1]
                                                                    g.slice <- xx[r,]
                                                                    g.slice$window.num <- yy
                                                                    g.slice$window.size <- 1
                                                                    g.slice$width.to.prev <- 0
                                                                    g.slice$dmr.width <- 0
                                                                    g.slice[, mcols.slice]
                                                                    
                                                                  }))
                                            }
                                            })
    }
#return(dmr.list)
    dmr.granges <- unlist(GRangesList(dmr.list[!sapply(dmr.list, is.null)]), use.names=FALSE)
    if(find.singletons) {
        singleton.granges <- unlist(GRangesList(singleton.list[!sapply(singleton.list, is.null)]), use.names=FALSE)
    } else {
        singleton.granges <- NULL
    }
    if(!merge.tables) {
        if(!is.null(dmr.granges)) {
            win.num <- sum(dmr.granges$width.to.prev == 0)
            dmr.granges$window.num <- rep(1:win.num, diff(c(which(dmr.granges$width.to.prev == 0), length(dmr.granges$width.to.prev)+1)))
        }
        if(!is.null(singleton.granges)) {
            singleton.granges$window.num <- 1:length(singleton.granges)
        }
        if(find.singletons) {
            return(list(dmr.granges=dmr.granges, singleton.granges=singleton.granges))
        } else {
            return(dmr.granges)
        }
    } else {
        if(!is.null(dmr.granges) && !is.null(singleton.granges)) {
            all.granges <- sort(c(dmr.granges, singleton.granges))
            #all.granges <- all.granges[order(all.granges$chr, all.granges$pos)]
            win.num <- sum(all.granges$width.to.prev == 0)
            all.granges$window.num <- rep(1:win.num, diff(c(which(all.granges$width.to.prev == 0), length(all.granges$width.to.prev)+1)))
        } else if (!is.null(singleton.granges)) {
            all.granges <- singleton.granges
        } else {
            all.granges <- dmr.granges
        }
        return(all.granges)
        if(!is.null(all.granges)) {
            all.granges <- all.granges[order(all.granges$chr, all.granges$pos)]
            win.num <- sum(all.granges$width.to.prev == 0)
            all.granges$window.num <- rep(1:win.num, diff(c(which(all.granges$width.to.prev == 0), length(all.granges$width.to.prev)+1)))
        }
        all.granges
    }
}

