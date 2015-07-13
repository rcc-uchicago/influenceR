# pick k numbers less than n. return S, T, where |S|=k and |T|=n-k, S union T = {1,2,..,n} 
# no longer used, switched to using a mask (see get_random_mask)
get_random_set <- function(n, k) {
    range <- 1:n
    samp <- sample(n, k)
    s <- range[range %in% samp]
    t <- range[! range %in% samp]
    return (list(s,t));
}



# KPP-Neg metric:
# Borghatti (9)
# maximize D_F = 1 - 2 * (sum(i>j) (1 / d_ij) / (n * n-1)
# equivalent to minimizing just the sum
metric_9 <- function(g) {
  n <- igraph::vcount(g)
  D <- igraph::shortest.paths(g)
  Dr <- 1/D
  s <- sum(Dr[lower.tri(D)]);
  df <- 1-2*s/(n*(n-1))
  
  return(df);
}

# return something like [TRUE, FALSE, FALSE, TRUE, ...] where there are n total values and k TRUEs
get_random_mask = function(n, k) {
    s <- sample(n, k);
    t <- 1:n %in% s;
    return (t);
}

kpp_neg = function(g, k, tol) {
    n <- igraph::vcount(g)

    s <- get_random_mask(n, k); # random index matrix to start with.
                                # our "target" nodes are {i | t[i] is TRUE}

    
    g_ <- igraph::delete.vertices(g, which(s))
    fit <- metric_9(g_)

    print(paste("Fit: ", fit))
    
    while (TRUE) {
        Dfit = 0
        pair = NULL
        for (u in which(s)) {
            for (v in which(!s)) {
                s_ <- s # clone
                s_[u] = FALSE;
                s_[v] = TRUE;
                
                g_ <- igraph::delete.vertices(g, which(s_))
                fit_ = metric_9(g_) # maximize fit
                d = fit_ - fit
                if ((d >= 0) && (d > Dfit)) {
                    Dfit = d
                    pair = list(u,v)
                }
                #print(paste("test",u,v,fit_))
            }
        }
        if (Dfit < tol)
            break
        u <- pair[[1]]
        v <- pair[[2]]
        s[u] = FALSE
        s[v] = TRUE
        fit = fit + Dfit
        print(paste("Fit: ", fit, " Time: ", format(Sys.time(), "%H:%M:%S")))
    }

    print(paste("New fit: ", fit));
    return(which(s));
}


