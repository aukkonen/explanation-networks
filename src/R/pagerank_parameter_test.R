require(AUC)
source('intervention.R')
source('bayes_net_test.R')
source('exp_data_loader.R')
source('generalization_aware_pruning.R')
require( rJava )
.jinit( classpath="tools.jar", parameters=c("-Xmx4G") )

run_bayes_net_test <- function( ROUNDS=20, E_fnc=compute_E_matrix ) {
    res <- data.frame( alg=character(0), algparam=double(0), auc=double(0) )
    for ( r in 1:ROUNDS ) {
        D <- generate_combined_path_data(num_antecedents=3, num_noise_attrs=6,
                                         num_blocks=3, N=5000,
                                         prx1=0.1, a=0.1, b=0.7)
        D <- to_numeric_target( D, 2, 1 )
        targetname <- 'T'
        thr <- 0

        x <- run_miner( D, targetname, 0.01, 3, thr )
        included <- x$quality >= thr
        x$wracc <- sqrt( x$size ) * x$quality

        C <- cover_matrix( x$description, D )

        E <- E_fnc( C, D[, targetname], included=included, idx_out=x$idx_out, sketchSize=0 )
        E <- E[included,included]

        E <- to_markov( E )

        x <- x[included, ]

        ## parent <- parent_node( x$idx_in, x$idx_out )

        ## Es <- summary(Matrix(E, sparse=TRUE))
        ## EE <- partition_E_by_search_tree_fast( Es, parent )
        ## EE$inter <- as.matrix(sparseMatrix( i=EE$inter$i, j=EE$inter$j, x=EE$inter$x ))
        ## EE$intra <- as.matrix(sparseMatrix( i=EE$intra$i, j=EE$intra$j, x=EE$intra$x ))

        causal <- as.factor( sapply( x$description, function(d) descr_quality_causal( d ) ) )

        weight <- x$quality

        cat( 'using all explanations\n' )
        pr_all <- alpha_vs_auc( E, causal, weight )
        pr_all <- data.frame( alg=rep( 'pr_all', nrow(pr_all) ), algparam=pr_all[,1], auc=pr_all[,2] )

        ## cat( 'using explanations across disjoint search tree branches\n' )
        ## tmp1 <- alpha_vs_auc( EE$inter, causal, weight )
        ## tmp1 <- data.frame( alg=rep( 'pr_inter', nrow(tmp1) ), algparam=tmp1[,1], auc=tmp1[,2] )
        
        ## cat( 'using explanations from within same search tree branch\n' )
        ## tmp2 <- alpha_vs_auc( EE$intra, causal, weight )
        ## tmp2 <- data.frame( alg=rep( 'pr_intra', nrow(tmp2) ), algparam=tmp2[,1], auc=tmp2[,2] )

        ## compute gap scores
        x$genscore <- generalization_qual( x$description, D, D[,targetname] )
        x$genscore2 <- sqrt( x$size ) * x$genscore
        
        tmp <- pr_all
        tmp <- rbind( tmp, data.frame( alg='wracc',
                                      algparam=NA,
                                      auc=auc(roc(x$wracc, causal)) ) )
        tmp <- rbind( tmp, data.frame( alg='quality',
                                      algparam=NA,
                                      auc=auc(roc(x$quality, causal)) ) )
        tmp <- rbind( tmp, data.frame( alg='genscore2',
                                      algparam=NA,
                                      auc=auc(roc(x$genscore2, causal)) ) )

        res <- rbind( res, tmp )
    }
    res
}

run_real_data_test <- function( dataname, targetclass, qualityThr ) {
    x <- load_data_and_subgroups( dataname, targetclass )

    mr <- x$mr
    mr <- subset( mr, quality >= qualityThr )
    cat( nrow(mr), 'subgroups remain after quality filtering\n' )
    mr$wracc <- sqrt( mr$size ) * mr$quality
    
    E <- compute_E_matrix( x$D, mr$description, x$D[, x$targetname] )
    E <- to_markov( E )

    alpha_vs_div( E, mr$quality, mr$description, x$D )
}

alpha_vs_div <- function( E, weight, description, data, k=10 ) {
    alphas <- seq( 0.05, 0.95, length.out=20 )
    t(sapply( alphas, function( alpha ) {
        score <- pagerank( E, weighted_pr_update( alpha, weight ) )
        tmp <- description[ order( score, decreasing=TRUE ) ]
        tmp <- tmp[1:k]
        C   <- cover_matrix( tmp, data )
        c( alpha, entropy_divfun( C ), coverage_divfun( C ) )
    } ))
}

alpha_vs_auc <- function( E, causal, weight ) {
    alphas <- seq( 0.05, 0.95, length.out=20 )
    t(sapply( alphas, function( alpha ) {
        cat( 'using alpha = ', alpha, '\n' )
        score <- pagerank( E, weighted_pr_update( alpha, weight ) )
        c( alpha, auc(roc(score, causal)) )
    } ))
}

## testoutput is the data frame built by run_bayes_net_test
make_all_algs_plot <- function( testoutput, dolegend=TRUE ) {
    require(plyr)
    foo <- ddply( testoutput, c('alg','algparam'),
                 function(x) data.frame(mu=mean(x$auc), dev=sd(x$auc)) )

    draw <- function( x, mu, dev, color ) {
        lines( x, mu, col=color, lwd=2 )
        lines( x, mu+3*dev, col=color, lwd=1, lty=2 )
        lines( x, mu-3*dev, col=color, lwd=1, lty=2 )
    }
    
    dev.new( width=5, height=4)
    par( mai=c(0.9, 0.9, 0.5, 0.1) )
    plot(NA, NA, xlim=c(0,1), ylim=c(0.5,1),
         main=sprintf( 'alpha vs AUC, n = %d', sum(testoutput$alg=='wracc')),
         bty='n', xlab=expression(alpha), ylab='AUC')

    pr0 <- subset( foo, alg=='pr_all' )
    pr0 <- pr0[ order( pr0$algparam ), ]
    ## pr1 <- subset( foo, alg=='pr_inter' )
    ## pr1 <- pr1[ order( pr1$algparam ), ]
    ## pr2 <- subset( foo, alg=='pr_intra' )
    ## pr2 <- pr2[ order( pr2$algparam ), ]
    wracc   <- subset( foo, alg=='wracc' )
    quality <- subset( foo, alg=='quality' )
    gap     <- subset( foo, alg=='genscore2' )

    draw( pr0$algparam, pr0$mu, pr0$dev, 'black' )
    ## draw( pr1$algparam, pr1$mu, pr1$dev, 'red' )
    ## draw( pr2$algparam, pr2$mu, pr2$dev, 'blue' )
    draw( c(0,1), rep( wracc$mu, 2 ), rep( wracc$dev, 2 ), 'red' )
    draw( c(0,1), rep( quality$mu, 2 ), rep( 0, 2 ), 'blue' )
    draw( c(0,1), rep( gap$mu, 2 ), rep( 0, 2 ), 'darkgreen' )

    if ( dolegend ) {
        legend( 'topleft',
               legend=c('pagerank', 'wracc', 'quality', 'gap'),
               lwd=2,
               bty='n',
               col=c('black', 'red','blue','darkgreen') )
    }
}

## testoutput is the data frame built by run_bayes_net_test
make_all_divs_plot <- function( testoutput ) {
    dev.new( width=10, height=3.5 )
    par( mfrow=c(1,3) )
    plot( testoutput[,1], testoutput[,2], main='ediv', xlab='alpha' )
    plot( testoutput[,1], testoutput[,3], main='cdiv', xlab='alpha' )
    plot( testoutput[,2], testoutput[,3], xlab='ediv', ylab='cdiv' )
}
