## require( rJava );
## .jinit(classpath="seedselect.jar", parameters="-Xmx4G");

## graph: n x 3 matrix with i, j, prob triplets
## k:     number of seeds to pick
seedselect <- function( i, j, p, k, num_samples=250 ) {
    ss <- .jnew( 'SeedSelector' )
    .jcall( ss, '[I', 'selectSeeds',
           as.integer( k ),
           as.integer( i ),
           as.integer( j ),
           as.numeric( p ),
           as.integer( max( c(i,j) ) ),
           as.integer( num_samples ) )
}
