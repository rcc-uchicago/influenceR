library(sna)


al2mat = function(al) {
    n <- max(unlist(al))
    m <- matrix(0, n, n)
    for (i in 1:length(al)) {
        x = al[i];
        m[x[0],x[1]] = 1
    }
    return(m);
}

# pick k numbers less than n. return S, T, where |S|=k and |T|=n-k, S union T = {1,2,..,n} 
# no longer used, switched to using a mask (see get_random_mask)
get_random_set <- function(n, k) {
    range <- 1:n
    samp <- sample(n, k)
    s <- range[range %in% samp]
    t <- range[! range %in% samp]
    return (list(s,t));
}


graph_distance <- function(x) { 
    return(geodist(x)$gdist);
}

reachability = function(A) {
    D = graph_distance(A)
    return (D > 0);
}


# KPP-Neg metric:

# number of disconnected pairs (Borgatti (3))
# minimize sum(i>j) R_ij
# Equation (4) should be more tractable in the optimization context (although lose out on elegance of equation)
metric_3 <- function(A) {
    R <- reachability(A);
    s <- sum(R[lower.tri(R)]);
    return(s);
}

# Borghatti (9)
# maximize D_F = 1 - 2 * (sum(i>j) (1 / d_ij) / (n * n-1)
# equivalent to minimizing just the sum
metric_9 <- function(A) {
    D <- graph_distance(A);
    Dr <- 1/D
    s <- sum(Dr[lower.tri(D)]);
    return(s);
}

# return something like [TRUE, FALSE, FALSE, TRUE, ...] where there are n total values and k FALSEs
get_random_mask = function(n, k) {
    s <- sample(n, k);
    t <- !1:n %in% s;
    return (t);
}

# trim array based on True, False mask (see get_random_mask above)
trimmed_array = function(A, m) {
    i <- which(m)
    B <- A[i,][,i]
    return(B);
}


greedy_optimize = function(A, k, metric, tolerance) {
    n <- length(A[1,])

    t <- get_random_mask(n, k); # random index matrix to start with.
                                # our "target" nodes are {i | t[i] is FALSE}

    B <- trimmed_array(A, t);
    fit <- metric(B)

    print(paste("Fit: ", fit))
    while (TRUE) {
        Dfit = Inf
        pair = NULL
        for (u in which(!t)) {
            for (v in which(t)) {
                t_ <- t # clone
                t_[v] = FALSE;
                t_[u] = TRUE;
                
                B_ = trimmed_array(A, t_)
                fit_ = metric(B_)
                d = fit - fit_
                if ((d >= 0) && (d < Dfit)) {
                    Dfit = d
                    pair = list(u,v)
                }
            }
        }
        if (Dfit < tolerance)
            break
        u <- pair[1]
        v <- pair[2]
        t[v] = FALSE
        t[u] = TRUE
        fit = fit - d
    }

    print(paste("New fit: ", fit));
    return(which(!t));
}


main = function(argv) {
    n = as.numeric(argv[1]) # number of nodes in graph
    p = as.numeric(argv[2]) # probability two nodes are connected
    k = as.numeric(argv[3]) # number of nodes to find
    tol = as.numeric(argv[4]) # tolerance to stop optimize algorithm at
    t = as.numeric(argv[5]) # times to repeat the optimize algorithm
    m = as.numeric(argv[6]) # metric (3 or 9)
   
    if (m == 3)
        metric <- metric_3
    else if (m == 9)
        metric <- metric_9
    else {
        print("Invalid metric! Use 3 or 9 only.")
        q()
    }

    G <- rgraph(n, p)

    for (i in 1:t) {
        S <- greedy_optimize(G, k, metric, tol)
        print(sprintf("nodes: %s", paste(S, collapse=",")))
        cat('\n')
    }
}


argv <- commandArgs(trailingOnly = TRUE);
if (length(argv) > 1)
    main(argv)

