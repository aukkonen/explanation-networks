joint_distr_vectors <- function( k ) {
    ## build joint distribution probabilities:
    tmp <- sapply( 1:k, function(i) {
        t(apply( combn( k, i ), 2, function(y) {
            x <- rep(1,k)
            x[y] <- 2
            x
        }))
    } )
    foo <- tmp[[1]]
    for ( i in 2:length(tmp) ) {
        foo <- rbind( foo, tmp[[i]] )
    }
    
    rbind( rep(1,k), foo )
}

sample_from_jointp <- function( k, p, N ) {
    vectors <- joint_distr_vectors( k )
    jointp  <- apply( vectors, 1, p )
    jointp  <- jointp/sum(jointp)

    vectors <- vectors - 1

    ## sample N instances from jointp
    data.frame( vectors[ sample( nrow(vectors), N, replace=TRUE, prob=jointp ), ] )
}

to_factors <- function( data ) {
    for ( i in 1:ncol(data) ) {
        data[,i] <- as.factor( data[,i] )
    }
    data
}

to_logical <- function( data ) {
    for ( i in 1:ncol(data) ) {
        data[,i] <- data[,i]==1
    }
    data
}

generate_test_data <- function( N=10000 ) {
    ## first build conditional probability tables:
    prB <- c(0.95,0.05)
    prE <- c(0.90,0.10)
    prAgivenBE <- matrix( c(0.10,0.90,0.35,0.95), nrow=2 )
    prCgivenA <- c(0.05,0.9)

    p <- function( x ) {
        pr <- prB[ x[1] ] * prE[ x[2] ]
        if ( x[3] == 2 ) {
            pr <- pr * prAgivenBE[ x[1], x[2] ]
        }
        else {
            pr <- pr * (1 - prAgivenBE[ x[1], x[2] ])
        }
        if ( x[4] == 2 ) {
            pr <- pr * prCgivenA[ x[3] ]
        }
        else {
            pr <- pr * (1 - prCgivenA[ x[3] ])
        }
        pr
    }

    data <- sample_from_jointp( 4, p, N )

    names(data) <- c( 'B', 'E', 'A', 'C' )
    data <- to_factors( data )
    data$R1 <- as.factor( runif(N)>0.7 )
    data$R2 <- as.factor( runif(N)>0.7 )
    data
}

generate_test_data_2 <- function( num_noise_attrs=10, N=10000 ) {
    ## first build conditional probability tables:
    prB <- c(0.95,0.15)
    prE <- c(0.90,0.10)
    prA1givenBE  <- matrix( c(0.01,0.05,0.05,0.95), nrow=2 )
    prA2givenB   <- c( 0.05, 0.9 )
    prCgivenA1A2 <- matrix( c(0.01, 0.05, 0.05, 0.95), nrow=2 )

    p <- function( x ) {
        B  <- x[1]
        E  <- x[2]
        A1 <- x[3]
        A2 <- x[4]
        C  <- x[5]
        pr <- prB[ B ] * prE[ E ]
        if ( A1 == 2 ) {
            pr <- pr * prA1givenBE[ B, E ]
        }
        else {
            pr <- pr * (1 - prA1givenBE[ B, E ])
        }
        if ( A2 == 2 ) {
            pr <- pr * prA2givenB[ B ]
        }
        else {
            pr <- pr * (1 - prA2givenB[ B ])
        }
        if ( C == 2 ) {
            pr <- pr * prCgivenA1A2[ A1, A2 ]
        }
        else {
            pr <- pr * (1 - prCgivenA1A2[ A1, A2 ])
        }
        pr
    }

    data <- sample_from_jointp( 5, p, N )
    names( data ) <- c( 'B', 'E', 'A1', 'A2', 'C' )
    data <- to_logical( data )
    for ( i in 1:num_noise_attrs ) {
        oldnames <- names( data )
        data <- cbind( data, runif(N)>0.7 )
        names( data ) <- c( oldnames, sprintf('R%d', i) )
    }
    data  
}

## Generates data from a bayesian network with structure X1 -> X2 -> ... -> T,
## where Xi are antecedents and T is the target attribute.
## In addition to this, there are some random attributes that are uncorrelated
## with each other as well as the Xi and T.
##
## Setting num_antecedents=2 creates the model X1 -> X2 -> T.
generate_path_data <- function(num_antecedents=2,
                               num_noise_attrs=2,
                               N=10000,
                               prx1=0.2, a=0.05, b=0.90) {
    ## prior probabilitis for X1
    prior <- c( 1-prx1, prx1 )

    ## conditional prob table for Xi+1 (as well as T) given Xi:
    ## [ Pr( Xi+1 = 0 | Xi = 0 ), Pr( Xi+1 = 0 | Xi = 1 ) ]
    ## [ Pr( Xi+1 = 1 | Xi = 0 ), Pr( Xi+1 = 1 | Xi = 1 ) ]
    ##
    ## or
    ##
    ## [ 1-a, 1-b ]
    ## [   a,   b ]
    cp <- matrix( c(1-a, a, 1-b, b), nrow=2 )

    p <- function( x ) {
        pr <- prior[ x[1] ]
        for ( i in 2:length(x) ) {
            pr <- pr * cp[ x[i], x[i-1] ]
        }
        pr
    }

    data <- sample_from_jointp( num_antecedents + 1, p, N )
    names( data ) <- c( paste( 'X', 1:num_antecedents, sep='' ), 'T' )
    data <- to_logical( data )
    ## every noise attr is a permutation of a randomly chosen antecendet.
    for ( i in 1:num_noise_attrs ) {
        oldnames <- names( data )
        data <- cbind( data, data[sample(nrow(data)), sample(num_antecedents,1)] )
        names( data ) <- c( oldnames, sprintf('R%d', i) )
    }
    data
}

generate_combined_path_data <- function(num_antecedents=2,
                                        num_noise_attrs=2,
                                        num_blocks=2,
                                        N=10000,
                                        prx1=0.05, a=0, b=0.8) {
    replace_names <- function( current_names, i ) {
        sapply( current_names, function(name) {
            name <- strsplit( name, split='' )[[1]]
            if ( name[1] == 'T' ) {
                sprintf( 'T%d', i )
            }
            else {
                sprintf( '%s%d%s', name[1], i, name[2] )
            }
        } )
    }

    D <- generate_path_data(num_antecedents, num_noise_attrs, N, prx1, a, b)
    names(D) <- replace_names( names(D), 1 )
    for ( i in 2:num_blocks ) {
        d <- generate_path_data(num_antecedents, num_noise_attrs, N, prx1, a, b)
        names(d) <- replace_names( names(d), i )
        D <- cbind( D, d )
    }

    varnames <- names(D)
    targets <- startsWith( varnames, 'T' )
    D$T <- apply( D[, targets], 1, function(x) max(as.numeric(x)) )
    
    for ( tname in varnames[targets] ) {
        D[,tname] <- NULL
    }
    
    D
}

## replaces all 1s in D$T with random numbers from N(mu, sigma),
## and all zeros with a sample from N(0, sigma).
to_numeric_target <- function( D, mu, sigma ) {
    target <- D$T
    D$T[ target==FALSE ] <- rnorm( sum(target==FALSE), mean=0, sd=sigma )
    D$T[ target==TRUE ]  <- rnorm( sum(target==TRUE), mean=mu, sd=sigma )
    D$T <- D$T - min( D$T ) + 1 ## make sure all target values are >= 1.
    D
}
