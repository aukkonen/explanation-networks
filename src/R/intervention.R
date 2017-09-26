require(Matrix)
source('small_utils.R')

##
fSminusT_fast <- function( C, targetvar, fS, included, idx_out, sketchSize=16 ) {
    CC <- matrix( 1, nrow=nrow(C), ncol=ncol(C) ) - C

    btime = proc.time()[1]    
    scorer <- .jnew( 'PruningScorer',
                    as.integer(CC),
                    as.integer(nrow(CC)),
                    as.integer(ncol(CC)),
                    as.integer(sketchSize) )
    fSminusT <- .jcall( scorer, '[D', 'scores',
                       as.integer(C), as.integer(ncol(C)), included,
                       as.integer(idx_out), targetvar, fS )
    fSminusT <- matrix( fSminusT, nrow=ncol(C), ncol=ncol(CC), byrow=TRUE )
    cat( 'fSminusT matrix done in', proc.time()[1]-btime, 'seconds.\n' )
    cat( 'PruningScorer accessed a data element', .jcall(scorer, 'J', 'getAndResetCost'),
        'times.\n' )

    fSminusT
}

fSminusT_R <- function( C, targetvar ) {
    CC <- matrix( 1, nrow=nrow(C), ncol=ncol(C) ) - C

    tvCC <- targetvar * CC

    cat( 'Computing fSminusT...\n' )
    print( system.time( fSminusT <- crossprod( C, tvCC )  ) )
    cat( 'fSminusT done!\n' )

    fSminusT
}

compute_E_matrix <- function( C, targetvar, use.java=TRUE, ... ) {
    n     <- nrow(C)
    tvC   <- targetvar * C
    fS    <- apply( tvC, 2, sum )
    Tsize <- apply( C, 2, sum )

    Exp_fSminusT <- as.matrix(fS) %*% t(as.matrix((1-Tsize/n)))

    if ( use.java ) {
        fSminusT <- fSminusT_fast( C, targetvar, fS, ... )
    } else {
        fSminusT <- fSminusT_R( C, targetvar )
    }
    
    E <- Exp_fSminusT - fSminusT
    E[ E<0 ] <- 0

    E
}

compute_E_matrix_old <- function( C, targetvar, use.java=TRUE, ... ) {
    n     <- nrow(C)
    tvC   <- targetvar * C
    fS    <- apply( tvC, 2, sum )
    Tsize <- apply( C, 2, sum )

    fSminusT <- fSminusT_fast( C, targetvar, fS, ... )

    E <- 1 - (fSminusT / ( as.matrix(fS) %*% rep(1, length(fS)) ))
    E <- n * E / ( as.matrix( rep(1, length(Tsize) ) ) %*% Tsize )

    E[ E<0 ] <- 0

    E
}

to_markov <- function( E ) {
    rs <- apply( E, 1, sum )
    for ( idx in which( rs==0 ) ) {
        E[idx,idx] <- 1
        rs[idx]    <- 1
    }
    E/rs
}

pagerank_update <- function( E, s, alpha, p ) {
    (1.0-alpha) * p + alpha * as.numeric(s %*% E)
}

basic_pagerank_update <- function( alpha=1.0, p=c(), n=0 ) {
    if ( length(p)==0 ) {
        p <- rep( 1/n, n )
    }
    function( E, s ) {
        pagerank_update( E, s, alpha, p )
    }
}

weighted_pr_update <- function( alpha, weight ) {
    p <- weight/sum(weight)
    basic_pagerank_update( alpha, p )
}

## This is the standard score that uses a time-invariant transition matrix.
## E MUST be row-normalised, i.e. all rows must sum to 1.
pagerank <- function( E, pr_update_fnc, tol=sqrt(.Machine$double.eps) ) {
    s <- t(runif( nrow(E) ))
    s <- s/sum(s)
    while ( TRUE ) {
        s_new <- pr_update_fnc( E, s )
        delta <- max(abs(s - s_new))
        cat( delta, '\n' )
        if ( delta < tol ) {
            break
        }
        s <- s_new
    }
    s_new
}

## for numeric targets, no not set targetclass (or set it to NA)
compute_score <- function(data,
                          subgroup_descriptions,
                          subgroup_sizes,
                          targetname,
                          targetclass=NA) {
    if ( is.na( targetclass ) ) {
        targetvar <- data[, targetname]
    } else {
        targetvar <- as.numeric( data[, targetname]==targetclass )
    }
    E <- compute_E_matrix( data, subgroup_descriptions, targetvar )
    E <- to_markov( E )

    prupfnc <- weighted_pr_update( 0.99, subgroup_sizes )

    pagerank( E, prupfnc )
}
