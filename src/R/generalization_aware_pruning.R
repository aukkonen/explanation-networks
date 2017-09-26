generalization_pruning <- function( description, C, target, thr=0.01 ) {
    pval_qual <- generalization_pval_qual( description, C, target )
    res <- pval_qual[,1] <= thr
    print( table(res) )
    hist( pval_qual[,1] )
    cbind( res, pval_qual[,2] )
}

subsets <- function(k) {
    L <- lapply( 1:(k-1), function(i) {
        t(apply( combn( k, i ), 2, function(x) {
            q <- rep(0,k)
            q[x] <- 1
            q
        } ))
    } )
    res <- L[[1]]
    if ( length(L) > 1 ) {
        for ( i in 2:length(L) ) {
            res <- rbind( res, L[[i]] )
        }
    }
    res
}

generalization_qual <- function( description, data, target ) {
    cat( 'running gap...\n' )
    
    is.subset <- function(X, Y) all( is.element(X,Y) )

    description_mean <- sapply( description, function(descr) {
        mean( target[ get_subgroup_cover( data, descr ) == 1 ] )
    } )
    
    description_sets <- lapply( description, function(descr) {
        strsplit( descr, split=' & ' )[[1]]
    } )

    cat( 'got', length(description_sets), 'description sets.\n' )

    ss <- list()
    ss[[16]] <- 1 ## this is a hack to make the ss list contain something

    global_mean <- mean( target )

    sapply( 1:length(description_sets), function(i) {
        dset <- description_sets[[i]]
        n <- length( dset )
        if ( is.null( ss[[n]] ) ) {
            ss[[n]] <- subsets( n )
        }
        genmean <- apply( ss[[n]], 1, function(idx) {
            descr <- paste( dset[idx==1], collapse=' & ' )
            mean( target[ get_subgroup_cover( data, descr )==1 ] )
        } )
        ref <- max(c(genmean, global_mean))
        description_mean[i] - ref
    } )
}

## this might be a bit stupid
generalization_pval_qual <- function( description, C, target ) {
    cat( 'running gap...\n' )
    
    is.subset <- function(X, Y) all( is.element(X,Y) )

    description_sets <- sapply( description, function(descr) {
        strsplit( descr, split=' & ' )[[1]]
    }, USE.NAMES=FALSE )

    cat( 'got', length(description_sets), 'description sets.\n' )

    unique_conditions <- unique( unlist( description_sets ) )

    totlen <- sapply( description_sets, length )
    
    targetmean <- sapply( 1:length(description_sets), function(i) mean( target[ C[,i]==1 ] ) )
    globalmean <- mean( target )

    X <- sapply( unique_conditions, function( cond ) {
        sapply( description_sets, function(dset) {
            is.subset( cond, dset )
        } )
    } )

    print( length( unique_conditions ) )
    print( length( description_sets ) )
    print( dim(X) )

    ## maxidx will have two columns, 1st with the maxidx, the 2nd with genscore
    maxidx <- t(sapply( 1:length(description_sets), function(i) {
        dset       <- description_sets[[i]]
        subsets    <- apply( as.matrix(X[,dset]), 1, sum ) == totlen
        subsets[i] <- FALSE ## remove dset itself
        subsets    <- which( subsets )
        if ( length( subsets ) > 0 ) {
            idx <- subsets[which.max( targetmean[subsets] )]
            out <- c(idx, targetmean[i]-targetmean[idx])
        } else {
            ## this means the only subset of dset is dset itself
            out <- c(i, targetmean[i]-globalmean)
        }
    } ))

    res <- sapply( 1:length(description_sets), function(i) {
        if ( maxidx[i] == i ) {
            ## we get here if dset was the only subset of dset,
            ## such cases will have a pvalue of 0, as we would like
            ## to always include them.
            out <- 0
        } else {
            out <- t.test( target[ C[,i]==1 ], target[ C[,maxidx[i,1]]==1 ], alternative='g' )$p.value
        }
        out
    } )
    
    cat( 'gap done.\n' )
    cbind(res, maxidx[,2])
}

test_prune <- function() {
    descr <- c( 'a', 'a & b', 'c', 'b', 'b & c', 'a & c' )
    qual  <- c( 0.7,     0.5, 0.6, 0.1,    0.3,      0.9 )

    ## data.frame( descr=descr, qual=qual, genqual=prune( descr, qual, NA, NA, NA ) )
    prune( descr, qual )
}
