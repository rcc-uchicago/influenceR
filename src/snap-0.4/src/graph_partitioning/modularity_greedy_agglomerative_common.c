#include "graph_defs.h"
#include "modularity_greedy_agglomerative.h"

int keytype; /* 0:dq This is the index to use for the value used in the heap */
double MaxKeyHistory[MAXHISTORY];
int MaxKeyHistoryIdx;
int MaxKeyHistoryFilled;
double sumMKH, sumMKH2;
