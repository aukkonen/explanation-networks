## The MIT License (MIT)

## Copyright (c) 2015 Antti Ukkonen (antti.ukkonen@gmail.com)

## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:

## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.

require(compiler)

run_miner <- function( data, target_name, minsupp, maxdepth, quality_thr,
                      small_target=FALSE, qualityfnc=get_quality_effsize,
                      targetclass=NA, output.data.frame=TRUE ) {
    ## clean up data st. target has no NAs
    data    <- subset( data, !is.na(data[,target_name]) )
    minsupp <- floor( minsupp*nrow(data) )
    cat( 'N = ', nrow(data), ' after removing target NAs, minsupp = ', minsupp, '\n' )
    quality <- qualityfnc( data, target_name, small_target, targetclass )
    selections <- prepare( data, target_name, minsupp )
    cat( 'Running miner...' )
    search_tree <- recurse(quality, 1:length(selections), minsupp,
                           selections, parent_cond=NULL, active_rows=rep(TRUE, nrow(data)),
                           depth=1, maxdepth=maxdepth)
    cat( 'done!\n' )
    result <- NULL
    if ( output.data.frame ) {
        results <- collect_results(quality_thr,
                                   list( cond='', size=0, quality=NA, subtree=search_tree ),
                                   1, 1)
        resdf <- data.frame(description=results$conds,
                            quality=results$qualities,
                            size=results$sizes,
                            idx_in=results$in_indices,
                            idx_out=results$out_indices)
        resdf <- resdf[ order( resdf$idx_in ), ]
        resdf <- resdf[2:nrow(resdf),]
        resdf$description <- as.character( resdf$description )
        ## resdf <- resdf[ order( resdf$quality, decreasing=TRUE ), ]
        ## this only works if quality is tstat
        ## resdf$effect_size <- resdf$quality/sqrt(resdf$size)
        row.names(resdf) <- NULL
        cat( 'Got', nrow(resdf), 'subgroups.\n' )
        result <- resdf
    }
    else {
        result <- search_tree
    }
    result
}

## Collects the results from a given search tree (as produced by recurse)
## into a list of lists. These can be combined into a data frame later.
collect_results <- cmpfun( function( quality_thr, node, in_idx, out_idx ) {
    conds       <- c()
    sizes       <- c()
    qualities   <- c()
    in_indices  <- c()
    out_indices <- c()

    my_in_idx <- in_idx

    in_idx <- in_idx + 1
    for ( childnode in node$subtree ) {
        subtree_results <- collect_results( quality_thr, childnode, in_idx, out_idx )

        conds       <- c( conds, subtree_results$conds )
        sizes       <- c( sizes, subtree_results$sizes )
        qualities   <- c( qualities, subtree_results$qualities )
        in_indices  <- c( in_indices, subtree_results$in_indices )
        out_indices <- c( out_indices, subtree_results$out_indices )

        out_idx <- subtree_results$out_idx
        in_idx  <- subtree_results$in_idx
    }
    if ( node$quality >= quality_thr | length(conds) > 0 ) {
        conds       <- c( conds, node$cond )
        sizes       <- c( sizes, node$size )
        qualities   <- c( qualities, node$quality )
        in_indices  <- c( in_indices, my_in_idx )
        out_indices <- c( out_indices, out_idx )
    }

    list( conds=conds, sizes=sizes, qualities=qualities,
         in_indices=in_indices, out_indices=out_indices, in_idx=in_idx, out_idx=out_idx+1 )
} )

plot_density_comparison <- function( data, descriptor_string, target_attr ) {
    sg <- select_subgroup( data, descriptor_string )
    global_target <- data[,target_attr]
    sg_target     <- sg[,target_attr]
    d_global <- density( global_target, na.rm=T )
    d_sg     <- density( sg_target, na.rm=T )
    quartz()
    plot( d_global$x, d_global$y, lwd=2, type='l', xlab=NA, ylab=NA )
    lines( d_sg$x, d_sg$y, col='red', lwd=2 )
    legend( 'topright', legend=c('global', 'subgroup'), fill=c('black', 'red') )
}

plot_comparison <- function( data, descriptor_string, target_attr ) {
    sg <- select_subgroup( data, descriptor_string )
    global_target <- factor( data[,target_attr] )
    sg_target     <- factor( sg[,target_attr], levels=levels(global_target) )
    t_global <- table( global_target )/nrow(data)
    t_sg     <- table( sg_target )/nrow(sg)
    quartz()
    barplot( cbind( t_global, t_sg ), beside=F, width=0.7, xlim=c(0,2),
            names.arg=c('global','subgroup'), legend.text=levels(global_target) )
}

get_subgroup_cover <- function( data, descriptor_string ) {
    with( data=data, eval( parse( text=descriptor_string ) ) )
}

get_all_covers <- function( data, subgroup_list ) {
    t(sapply( 1:length(subgroup_list), function(i) {
        as.numeric( get_subgroup_cover( data, subgroup_list[i] ) )
    } ))
}

## returns that subset of data that satisfies the given descriptor string
select_subgroup <- function( data, descriptor_string ) {
    data[ get_subgroup_cover( data, descriptor_string ), ]
}

cluster_topk_subgroups <- function( data, subgroup_list ) {
    covers <- get_all_covers( data, subgroup_list )
    row.names(covers) <- subgroup_list$description
    d <- dist( covers, method='binary' )
    hc <- hclust( d, method='average' )
    quartz(width=10, height=7);
    par(mai=c(1, 1, 1, 5));
    plot(as.dendrogram(hc), horiz=T)
    hc
}

plot_topk_subgroups <- function( data, subgroup_list ) {
    covers <- get_all_covers( data, subgroup_list )
    pca <- princomp( t(covers), scores=TRUE )
    plot( pca$scores[,1], pca$scores[,2] )
}

## Implements a simple depth-first exhaustive mining algorithm.
## Returns a nested list of lists that represent the search tree.
##
## The structure of the list is
## [ [node, node_children], [node, node_children], ... ],
## where node is a tuple that encodes the condition at that
## node of the search tree and node_children is
## a list of nodes of the same type.
recurse <- cmpfun( function(quality, attrs, minsupp,
                    selections, parent_cond, active_rows,
                    depth, maxdepth ) {
    tree <- list()
    activesize <- sum( active_rows, na.rm=TRUE )
    for ( a in attrs ) {
        selectionsA <- selections[[a]]
        for ( sel in selectionsA ) {
            subgroup <- active_rows & sel$rows
            sgsize   <- sum( subgroup, na.rm=TRUE )
            ## we only insert this as a new node if it
            ## satisfies the minsupp constraint, as well as
            ## produces a nontrivial specialisation, i.e.
            ## does not match all current active rows.
            if ( (sgsize >= minsupp) && (sgsize < activesize) ) {
                node_cond <- sel$cond
                if ( !is.null( parent_cond ) ) {
                    node_cond <- sprintf( '%s & %s', parent_cond, sel$cond )
                }
                node_quality <- quality(subgroup)
                if ( is.na( node_quality ) ) {
                    print( table( subgroup ) )
                    cat( 'quality is NA!!\n' )
                    stop( node_cond )
                }
                node_size    <- sum(subgroup)
                node_subtree <- list()
                if ( depth < maxdepth ) {
                    node_subtree <- recurse(quality,
                                            attrs[ attrs > a ],
                                            minsupp,
                                            selections,
                                            node_cond,
                                            subgroup,
                                            depth+1, maxdepth)
                }
                tree <- c( tree, list(list(cond=node_cond,
                                           quality=node_quality,
                                           size=node_size,
                                           subtree=node_subtree)) )
            }
        }
    }
    tree
} )

get_quality_tstat <- function( data, target_name, small_target=FALSE ) {
    target      <- as.numeric(data[,target_name])
    global_mean <- mean(target)
    ## global_sdev <- sd(target)
    if ( small_target ) {
        cmpfun( function( S ) {
            targetS <- target[S]
            length(targetS)^0.5 * (global_mean - mean(targetS)) / sd(targetS)
        } )
    }
    else {
        cmpfun( function( S ) {
            targetS <- target[S]
            length(targetS)^0.5 * (mean( targetS ) - global_mean) / sd(targetS)
        } )
    }
}

get_quality_effsize <- function( data, target_name, small_target=FALSE, targetclass=NA ) {
    if ( is.na( targetclass ) ) {
        cat( 'get_quality_effsize: WARNING: targetclass not specified! using numeric target.\n' )
        target <- as.numeric( data[, target_name ] )
    } else {
        target <- as.numeric( data[, target_name ]==targetclass )
    }
    global_mean <- mean(target)
    global_sdev <- sd(target)
    if ( small_target ) {
        cmpfun( function( S ) {
            (global_mean - mean( target[ S ] ))/global_sdev
        } )
    }
    else {
        cmpfun( function( S ) {
            (mean( target[ S ] ) - global_mean)/global_sdev
        } )
    }
}

get_quality_f1measure <- function( data, target_name, small_target=FALSE, targetclass=NA ) {
    ## this only applies to binary targets! make sure this is the case:
    if ( length(unique(data[,target_name])) != 2 ) {
        stop( 'must have binary target!' )
    }
    target <- as.numeric( data[,target_name] )
    print( table( target ) )
    N      <- sum(target)
    cmpfun( function(S) {
        x    <- sum(target*S)
        if ( x == 0 ) {
            return(0)
        }
        prec <- x/sum(S)
        reca <- x/N
        2*prec*reca/(prec+reca)
    } )
}

get_quality_jaccard <- function( data, target_name, small_target=FALSE, targetclass=NA ) {
    ## this only applies to binary targets! make sure this is the case:
    if ( length(unique(data[,target_name])) != 2 ) {
        stop( 'must have binary target!' )
    }
    target <- as.numeric( data[,target_name] )
    print( table( target ) )
    N      <- sum(target)
    cmpfun( function(S) {
        x <- sum(target*S)
        x / (N + sum(S) - x)
    } )
}

llratio_fnc <- function( target, condition=c() ) {
    target <- factor( target )
    ## this assumes target is a factor
    target <- as.numeric( target )

    if ( length( condition ) == 0 ) {
        targetZero <- target==1
        targetOne  <- target==2
    }
    else {
        targetZero <- target==1 & condition
        targetOne  <- target==2 & condition
    }

    ## +2 for Laplace smoothing of probability estimates...
    targetZeroSum <- sum(targetZero) + 2
    targetOneSum  <- sum(targetOne) + 2

    cat('targetZeroSum', targetZeroSum, ', targetOneSum', targetOneSum, '\n' )

    gamma <- targetZeroSum/targetOneSum

    cmpfun( function(S) {
        gamma * (sum(targetOne & S)+1)/(sum(targetZero & S)+1)
    } )
}

## this only works with binary targets
get_quality_likelihoodratio <- function( data, target_name, small_target=FALSE ) {
    if ( length(unique(data[,target_name])) != 2 ) {
        stop( "Only binary targets allowed for likelihood ratio quality.")
    }
    llratio_fnc( data[,target_name] )
}

## Runs a simple preprocessing that precomputes the
## selections of individual conditions.
prepare <- function( data, target_name, minsupp ) {
    cat( 'Preprocessing...' )
    attr_names <- names(data)
    selections <- unlist(lapply( which( attr_names != target_name ), function(attr_idx) {
        compute_attr_selections( data, attr_idx, minsupp )
    } ), recursive=FALSE, use.names=FALSE)
    cat( 'done!\n' )
    selections
}

## data: a data frame
## attribute: an index
compute_attr_selections <- function( data, attribute, minsupp, max.unique.numeric=9 ) {
    attr_name = names(data)[attribute]
    output     <- NULL
    output_leq <- NULL
    output_gt  <- NULL
    attribute_class <- class( data[,attribute] )
    if ( attribute_class == 'numeric' || attribute_class == 'integer' ) {
        uniqvals <- unique(data[,attribute])
        if ( length(uniqvals) > max.unique.numeric ) {
            splitpoints <- quantile( data[,attribute], probs=seq(0.1,0.9,by=0.1), na.rm=TRUE )
            splitpoints <- unique( splitpoints )
        }
        else {
            splitpoints <- uniqvals[1:(length(uniqvals)-1)]
        }

        ## create selections that use <=
        output_leq <- lapply( splitpoints, function(sp) {
            rows <- data[,attribute] <= sp
            rows[ is.na(rows) ] <- FALSE
            list(
                cond=sprintf('%s <= %.2f', attr_name, sp),
                rows=rows
            )
        } )
        output_leq <- Filter(function(sel) sum(sel$rows, na.rm=TRUE) >= minsupp &
                                 sum(sel$rows, na.rm=TRUE) < nrow(data),
                             output_leq)

        ## create selections that use >
        output_gt <- lapply( splitpoints, function(sp) {
            rows <- data[,attribute] > sp
            rows[ is.na(rows) ] <- FALSE
            list(
                cond=sprintf('%s > %.2f', attr_name, sp),
                rows=rows
            )
        } )
        output_gt <- Filter(function(sel) sum(sel$rows, na.rm=TRUE) >= minsupp &
                                sum(sel$rows, na.rm=TRUE) < nrow(data),
                            output_gt)

        output <- list( output_leq, output_gt )
    }
    else if ( attribute_class == 'factor' ) {
        splitpoints <- levels( data[,attribute] )
        output <- lapply( splitpoints, function(sp) {
            rows <- data[,attribute] == sp
            rows[ is.na(rows) ] <- FALSE
            list(
                cond=sprintf("%s == '%s'", attr_name, sp),
                rows=rows
            )
        } )
        output <- Filter(function(sel) sum(sel$rows, na.rm=TRUE) >= minsupp &
                             sum(sel$rows, na.rm=TRUE) < nrow(data),
                         output)
        output <- list( output )
    }
    ## NOTICE: For attributes of type 'logical' we only consider TRUE values,
    ##         ie, no negative conditions will occurr in the resulting descriptions.
    ## To include negative conditions, convert all logical attributes to factors
    ## before running the miner.
    else if ( attribute_class == 'logical' ) {
        rows <- data[,attribute]
        rows[ is.na(rows) ] <- FALSE
        output <- list( list(cond=sprintf("%s", attr_name),
                             rows=rows) )
        output <- Filter(function(sel) sum(sel$rows, na.rm=TRUE) >= minsupp &
                             sum(sel$rows, na.rm=TRUE) < nrow(data),
                         output)
        output <- list( output )
    }
    else {
        print(sprintf('Unknown attribute type: %s. Only numeric and factors supported.\n',
              attribute_class))
    }
    output
}







## collect_results2 <- cmpfun( function( search_tree ) {
##     conds     <- c()
##     sizes     <- c()
##     qualities <- c()
##     for ( node in search_tree ) {
##         conds           <- c( conds, node$cond )
##         sizes           <- c( sizes, node$size )
##         qualities       <- c( qualities, node$quality )

##         subtree_results <- collect_results2( node$subtree )

##         conds           <- c( conds, subtree_results$conds )
##         sizes           <- c( sizes, subtree_results$sizes )
##         qualities       <- c( qualities, subtree_results$qualities )
##     }
##     list( conds=conds, sizes=sizes, qualities=qualities )
## } )

