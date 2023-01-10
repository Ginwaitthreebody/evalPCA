## copied from /home/genis/recent_admix/new2022/analyses/simulations/scripts/simulateRecAdmix.R
### and /home/genis/recent_admix/new2022/analyses/simulations/scripts/simFreqs.R


library(bnpsd)
## this is a more diverged ghost
treestring <- "(((pop1:0.1,popGhost:0.2):0.05,pop2:0.3):0.1,pop3:0.5)T;"

## this is a less diverged ghost
#treestring <- "(((pop1:0.1,popGhost:0.05):0.2,pop2:0.3):0.1,pop3:0.5)T;"

M <- 10e6

models <- list("ghost1"="(((pop1:0.1,popGhost:0.2):0.05,pop2:0.3):0.1,pop3:0.5)T;",
               "ghost2"="(((pop1:0.1,popGhost:0.05):0.2,pop2:0.3):0.1,pop3:0.5)T;",
               "ghost1b"="(((pop1:0.01,popGhost:0.02):0.005,pop2:0.03):0.01,pop3:0.05)T;")


## draw ancestral frequencies enriched for rare variants to mimic human sfs (based on bnpsd documentation)
f_anc<- draw_p_anc(M, beta = 0.03)
#k <- f_anc > 0.05 & f_anc < 0.95


for(scenario in names(models)){


    tree <- ape::read.tree(text=models[[scenario]])
    outdir <- paste0("/home/genis/other/evalPCA/sims_tosong/ghosts/", scenario)
    outped <- paste0(outdir, "/data/", scenario)
    outpng <- paste0(outdir,"/",  scenario, "_tree.png")

    bitmap(outpng, h=4,w=4,res=300)
    plot(tree)
    dev.off()


    #f_anc <- draw_p_anc(M, p_min = 0.05, p_max = 0.95)

    f <- draw_p_subpops_tree(f_anc, tree)

    q <- matrix(
        c(rep(c(1, 0, 0, 0), 50),
          rep(c(0,0,1,0), 50),
          rep(c(0,0,0,1), 50),
          rep(c(0, 0.3, 0.7,0), 50)),
        ncol=4, byrow=T)


    pi <- q %*% t(f)

    # sample genotypes, and do maf filte to keep only 0.05
    g <- matrix(rbinom(length(pi), 2, pi), ncol=ncol(pi))
    f_sample <- colMeans(g) / 2
    k <- f_sample > 0.05 & f_sample < 0.95

    ### now separate to convert go ped
    a1 <- g[,k] > 0
    a2 <- g[,k] > 1

    # genotypes in ped format
    gped <- matrix(paste(as.integer(a1) + 1, as.integer(a2) + 1), ncol=ncol(a1), nrow=nrow(a1))


    iid <- paste0("ind",1:nrow(q))
    famid <- rep(c("pop1", "pop2", "pop3", "pop4"), each=50)

    # ped file with info on indivdiual and their genotypes
    ped <- cbind(iid, famid,
                 matrix(0, ncol=3,nrow=nrow(q)),
                 rep("-9", nrow(q)),
                 gped)
    write.table(file=paste0(outped, ".ped"), x = ped, col.names=F, row.names=F, quote=F, sep=" ")

    # .map file with info on snps
    map <- cbind(rep(1, sum(k)),
                 paste0("snp", 1:sum(k)),
                 rep(0, sum(k)),
                 1:sum(k))
    write.table(file=paste0(outped, ".map"), x = map, col.names=F, row.names=F, quote=F, sep=" ")


    ## file to force plink to preserve major/minor allele so individual freqs are consistent
    refinfo <- cbind(map[,2],
                     rep(1,sum(k)))
    write.table(file=paste0(outped, ".refallele"), x = refinfo, col.names=F, row.names=F, quote=F, sep=" ")

}
