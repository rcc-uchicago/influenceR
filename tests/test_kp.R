library(influenceR)

g <- dimacs.to.graph("/project/jschnei1/sdjacobs-workspace/R0_net_results/R0-undir.dim")


keyplayer(g, 50, maxsec=15, roundsec=5)

