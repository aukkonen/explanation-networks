require( compiler )

cover_matrix <- function( description_list, data ) {
    cat( 'collecting cover matrix...' )
    Cmat <- sapply( 1:length(description_list), function(i) {
        x <- as.numeric( get_subgroup_cover( data, description_list[i] ) )
        x[ is.na(x) ] <- 0
        x
    } )
    cat( 'done! Got matrix of size', dim(Cmat), '\n' )
    Cmat
}

test_search_tree_edges <- function() {
    name   <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O' )
    idxin  <- 1:15
    idxout <- c( 5 ,  3 ,  1 ,  2 ,  4 , 10 ,  8 ,  6 ,  7 ,  9 , 15 , 11 , 14 , 12 , 13 )
    tmp <- data.frame( name=name, idxin=idxin, idxout=idxout )
    tmp$name <- as.character( tmp$name )
    tmp
}

parent_node <- function( idx_in, idx_out ) {

    cat( 'computing parent nodes...\n' )

    ancestors <- function( Y ) which( (idx_in < idx_in[Y]) & (idx_out > idx_out[Y]) )

    parent    <- function( Y ) {
        ancY <- ancestors( Y )
        ancY[which.max( idx_in[ ancY ] )]
    }

    res <- sapply( 1:length( idx_in ), function( Y ) {
        p <- parent( Y )
        if ( length( p ) == 0 ) {
            return(-1)
        }
        else {
            return(p)
        }
    } )

    cat( 'done!\n' )
    res
}

make_parent_edges <- function( parent, name ) {
    t(sapply( 1:length(parent), function( Y ) {
        p <- parent[Y]
        if ( p == -1 ) {
            return(c('root', name[Y]))
        }
        else {
            return(c(name[p], name[Y]))
        }
    } ))
}

output_dot <- function( edges ) {
    cat( 'digraph T {\n' )
    for ( i in 1:nrow(edges) ) {
        edgestr <- sprintf( '"%s" -> "%s"\n', edges[i,1], edges[i,2] )
        cat( edgestr )
    }
    cat( '}\n' )
}

output_dot_E <- function( E, names, thr ) {
    maxweight <- max(E[])
    for ( T in 1:ncol(E) ) {
        for ( S in 1:nrow(E) ) {
            if ( (S != T) & (E[S,T] >= thr) & (E[T,S] < thr) ) {
                cat( sprintf('"%s" -> "%s" [color="red", penwidth=%.3f]\n',
                             names[T], names[S], ((E[S,T]/maxweight)^2)*100) )
            }
        }
    }
}

## returns TRUE when u and v are in the same subtree
same_subtree <- cmpfun(function( u, v, parent ) {
    
    path_to_root <- cmpfun(function( vertex, parent ) {
        par <- c(vertex)
        while ( parent[vertex] > -1 ) {
            par <- c( par, parent[vertex] )
            vertex <- parent[vertex]
        }
        par
    })

    pu <- path_to_root(u, parent)
    pv <- path_to_root(v, parent)

    length( intersect(pu, pv) ) > 0
})

inter_subtree_explanations <- cmpfun( function( E, parent ) {
    res <- c()
    for ( i in 1:(nrow(E)-1) ) {
        for ( j in (i+1):ncol(E) ) {
            ss <- same_subtree( i, j, parent )
            res <- rbind( res, c( i, j, E[i,j], ss ) )
            res <- rbind( res, c( j, i, E[j,i], ss ) )
        }
    }
    res <- data.frame( res )
    names(res) <- c('i', 'j', 'E', 'ss' )
    res$ss <- (res$ss==1)
    res
} )

## creates two copies of E, one where only within subtree branch edges
## are present, and another that is the complement of this.
partition_E_by_search_tree <- cmpfun( function( E, parent ) {
    cat( 'building inter/intra edge split...\n' )
    E_inter <- matrix( 0, nrow=nrow(E), ncol=ncol(E) )
    E_intra <- matrix( 0, nrow=nrow(E), ncol=ncol(E) )

    for ( i in 1:(nrow(E)-1) ) {
        for ( j in (i+1):ncol(E) ) {
            if ( same_subtree( i, j, parent ) ) {
                E_intra[i,j] <- E[i,j]
                E_intra[j,i] <- E[j,i]
            }
            else {
                E_inter[i,j] <- E[i,j]
                E_inter[j,i] <- E[j,i]
            }
        }
    }

    diag( E_inter ) <- diag( E )
    diag( E_intra ) <- diag( E )

    ## renormalise rows!
    E_inter <- E_inter / apply( E_inter, 1, sum )
    E_intra <- E_intra / apply( E_intra, 1, sum )

    cat( 'done!\n' )
    
    list( inter=E_inter, intra=E_intra )
} )

partition_E_by_search_tree_fast <- function( E, parent ) {
    cat( 'building inter/intra edge split...\n' )

    fu <- .jnew( 'FastUtils' )
    isIntraEdge <- .jcall( fu, '[Z', 'treePartitioner',
                          as.integer(E$i-1), as.integer(E$j-1), as.integer(parent-1) )

    Einter <- list( i=E$i[ !isIntraEdge ], j=E$j[ !isIntraEdge ], x=E$x[ !isIntraEdge ] )
    Eintra <- list( i=E$i[ isIntraEdge ], j=E$j[ isIntraEdge ], x=E$x[ isIntraEdge ] )
    
    cat( 'done!\n' )
    list( inter=Einter, intra=Eintra )
}


## X is a matrix with subgroups
## sorted in decreasing order of some quality measure or score
##
## divfun is a function that takes
## a 0-1 matrix with the cover vectors as columns
##
## D is the original dataset
cumulative_diversity <- function( X, divfun, D ) {
    C <- c()
    d <- c()
    for ( i in 1:nrow(X ) ) {
        x <- get_subgroup_cover( D, X$description[i] )
        x[ is.na(x) ] <- 0
        C <- cbind( C, x )
        d <- c( d, divfun( C ) )
    }
    d
}

## returns TRUE if descriptor doesn't contain the letter 'R',
## and returns FALSE if it does.
descr_quality_causal <- function( descriptor ) {
    !("R" %in% strsplit( descriptor, split='' )[[1]])
}

## returs fraction of rows that have at least one 1 in C
coverage_divfun <- function( C ) {
    sum(as.numeric( apply( C, 1, sum ) >= 1 ))/nrow(C)
}

## returns entropy of joint distribution
entropy_divfun <- function( C ) {
    p <- table( apply( C, 1, function(x) paste( x, collapse='' ) ) )
    p <- p/sum(p)
    -1 * sum(log2(p)*p)
}

prediction_RMSE <- function( X, D, targetname ) {
    mu <- mean( D[, targetname] )
    C  <- sapply( X$description, function(desc) {
        x <- get_subgroup_cover( D, desc )
        x[ is.na(x) ] <- 0
        x
    } )
    means <- apply( C, 2, function(x) mean( D[ x==1, targetname ] ) )
    C <- C * as.matrix(rep(1,nrow(C))) %*% means
    pred <- apply( cbind(C,mu), 1, max )
    sqrt(mean((pred - D[,targetname])^2))
}

pairwise_cover_dist <- function( X, D, targetname ) {
    C <- cover_matrix( X$description, D )
    dmat <- as.matrix( dist( t(C), method='binary' ) )
    sum( triu( dmat, k=1 ) )
}

split_train_test <- function( D, C, testfraction=0.25 ) {
    n <- nrow(D)
    testidx <- sample( n, floor(0.25*n) )
    list( Test=D[ testidx, ], Train=D[ -testidx, ], CTrain=C[ -testidx, ] )
}
