greedy_wracc <- function( description, C, data, targetname, wracc, k ) {
    solution <- c()
    size <- rep(0, ncol(C))
    while( (nrow(data) > 0) & (length(solution) < k) ) {
        idx <- which.max( wracc )
        cat( 'max(wracc) =', wracc[idx], 'size =', size[idx], '\n' )
        if ( wracc[idx] < 0 ) {
            ## only continue for positive wracc
            break
        }
        solution <- c( solution, idx )
        cat( 'added', description[idx], '\n' )

        ## update data by removing everything that matches the selected description
        keep <- C[,idx]==0
        data <- data[ keep, ]
        C    <- C[ keep, ]
        size <- apply( C, 2, sum )

        cat( 'nrow(data) =', nrow(data), '\n' )

        ## recompute wracc
        qfnc <- get_quality_effsize( data, targetname )
        wracc <- sqrt( size ) * apply( C, 2, function(S) qfnc(S==1) )

        if ( sum( is.na( wracc ) ) == length( wracc ) ) {
            ## this means we have covered all target instances
            ## and there are no more interesting subgroups
            break
        }
    }
    
    solution
}
