source('intervention.R')
source('bayes_net_test.R')
source('seedselect.R')
source('generalization_aware_pruning.R')
source('greedy_wracc_baseline.R')
source('ssgm.R')
## require(xlsx)
## require(rJava)

evaluate <- function( topk, name, D, targetname, relevant ) {
    ccov    <- cumulative_diversity( topk, coverage_divfun, D )
    cent    <- cumulative_diversity( topk, entropy_divfun, D )
    cent2   <- cumulative_diversity( topk, entropy_divfun,
                                    D[D[,targetname]>mean(D[,targetname]),])
    ok <- sum( sapply( topk$description, function(desc) descr_quality_causal( desc ) ) )
    prec <- ok/nrow(topk)
    reca <- ok/relevant
    qfnc    <- get_quality_effsize( D, targetname )
    C       <- cover_matrix( topk$description, D )
    qualsum <- sum( apply( C, 2, function(S) {
        qfnc( S==1 )
    } ) )
    pcd <- pairwise_cover_dist( topk, D )
    n <- nrow(topk)
    as <- mean( apply( C, 2, sum )/nrow(D) )
    data.frame( algname=name, ccov=ccov[n], cent=cent[n], cent2=cent2[n],
               qual=qualsum, avgqual=qualsum/n, pcd=pcd, n=n, avgsize=as,
               prec=prec, recall=reca, f1=2*(prec*reca)/(prec+reca) )
}

run_script <- function( N, antecedents, noise, blocks, k ) {
    require(rJava)
    .jinit( classpath="tools.jar", parameters=c("-Xmx4G") )
    D <- generate_combined_path_data(num_antecedents=antecedents,
                                     num_noise_attrs=noise,
                                     num_blocks=blocks, N=N,
                                     prx1=0.1, a=0.1, b=0.7)
    D <- to_numeric_target( D, 2, 1 )
    targetname <- 'T'
    thr <- 0
    mr  <- run_miner( D, targetname, 0.01, 3, thr )
    x   <- list( mr=mr, D=D, targetname=targetname )
    x$C <- cover_matrix( x$mr$description, x$D )
    relevant <- sum( sapply( x$mr$description, function(desc) descr_quality_causal( desc ) ) )

    x$mr$wracc <- sqrt(x$mr$size)*x$mr$quality

    ## If x$mr has more than MAXGROUPS rows, keep only those having high enough wracc
    MAXGROUPS <- 5000
    n <- nrow(x$mr)
    if ( n > MAXGROUPS ) {
        wracc_sorted <- sort(x$mr$wracc, decreasing=TRUE)
        wracc_lb     <- wracc_sorted[ MAXGROUPS ]
        idx  <- x$mr$wracc >= wracc_lb
        x$mr <- x$mr[ idx, ]
        x$C  <- x$C[,idx]
    }

    ## x$mr <- subset( x$mr, x$mr$quality >= thr )
    included   <- x$mr$quality >= thr

    E <- compute_E_matrix( x$C, x$D[, x$targetname],
                          included=included, idx_out=x$mr$idx_out, sketchSize=0 )
    E <- E[included,included]
    E <- to_markov( E )

    x$mr <- x$mr[included,]
    x$C  <- x$C[,included]

    cat( nrow(x$mr), 'subgroups remain.\n' )

    results <- c()
    
    ## add some baseline rankers and pagerank based score
    ## compute pagerank using quality as teleportation weight
    x$mr$pagerank <- pagerank( E, weighted_pr_update(0.7, x$mr$quality ) )
    pagerank <- head( order( x$mr$pagerank, decreasing=T ), k )
    results <- rbind( results, evaluate( x$mr[ pagerank, ], 'pagerank', x$D, x$targetname, relevant ) )

    Esparse  <- summary(Matrix(E, sparse=TRUE))

    ## use the entire matrix E to compute seeds
    print( system.time( seeds <- seedselect( Esparse$j, Esparse$i, Esparse$x, k ) ) )
    results <- rbind( results, evaluate( x$mr[ seeds, ], 'seeds', x$D, x$targetname, relevant ) )

    ## run generalization aware filtering
    x$mr$genscore  <- generalization_qual( x$mr$description, x$D, x$D[,x$targetname] )
    x$mr$genscore2 <- sqrt( x$mr$size ) * x$mr$genscore
    gap2    <- head( order( x$mr$genscore2, decreasing=T ), k )
    results <- rbind( results, evaluate( x$mr[ gap2, ], 'gap2', x$D, x$targetname, relevant ) )

    ## add baselines
    wracc   <- head( order( x$mr$wracc, decreasing=T ), k )
    quality <- head( order( x$mr$quality, decreasing=T ), k )
    results <- rbind( results, evaluate( x$mr[ wracc, ], 'wracc', x$D, x$targetname, relevant ) )
    results <- rbind( results, evaluate( x$mr[ quality, ], 'quality', x$D, x$targetname, relevant ) )

    ## run greedy wracc baseline
    greedywracc <- greedy_wracc( x$mr$description, x$C, x$D, x$targetname, x$mr$wracc, k )
    results <- rbind( results, evaluate( x$mr[greedywracc,], 'greedy-wracc', x$D, x$targetname, relevant ) )

    ## normalized version of measures to compute final score
    ## all scores will be normalized by their maximum value in results
    results$cent2_norm   <- results$cent2/max(results$cent2)
    results$avgqual_norm <- results$avgqual/max(results$avgqual)
    results$avgsize_norm <- results$avgsize/max(results$avgsize)
    results$ccov_norm    <- results$ccov/max(results$ccov)

    results$f1[ is.nan(results$f1) ] <- 0
    results <- results[ order( results$f1 ), ]

    ## compute final score and save table
    results$score <- (results$avgqual_norm * results$avgsize_norm * results$cent2_norm * results$ccov_norm)^0.25

    print( results, digits=2 )
    ## write.xlsx( results, file=sprintf('../exp/tables/bayes/N%d_ant%d_noise%d_bl%d.xlsx', N, antecedents, noise, blocks) )
}


### SCRIPT BEGINS HERE ####
k <- 20
N           <- 2000
antecedents <- 5
noise       <- 5
blocks      <- 10

run_script( N, antecedents, noise, blocks, k )
