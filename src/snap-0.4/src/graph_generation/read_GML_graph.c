#include "graph_defs.h"
#include "graph_gen.h"

#define LINELENGTH 1000

/* Types */

typedef struct {
    int target;        /* Index in the vertex[] array of 
                          neighboring vertex.
                          (Note that this is not necessarily equal to the GML
                          ID of the neighbor if IDs are nonconsecutive or do
                          not start at zero.*/
    double weight;     /* Weight of edge. 1 if no weight is specified. */
    int id;
} EDGE;

typedef struct {
    int id;            /* GML ID number of vertex */
    int degree;        /* Degree of vertex (out-degree for directed nets) */
    char *label;       /* GML label of vertex.  NULL if no label specified */
    EDGE *edge;        /* Array of EDGE structs, one for each neighbor */
} VERTEX;

typedef struct {
    int nvertices;     /* Number of vertices in network */
    int directed;      /* 1 = directed network, 0 = undirected */
    VERTEX *vertex;    /* Array of VERTEX structs, one for each vertex */
} NETWORK;


typedef struct line {
    char *str;
    struct line *ptr;
} LINE;

/* Globals */

LINE *first;
LINE *current;

/* Function to read one line from a specified stream.  Return value is */
/* 1 if an EOF was encountered.  Otherwise 0. */
int read_line(FILE *stream, char line[LINELENGTH]) {
    if (fgets(line,LINELENGTH,stream)==NULL) return 1;
    line[strlen(line)-1] = '\0';   /* Erase the terminating NEWLINE */
    return 0;
}


/* Function to read in the whole file into a linked-list buffer, so that we
   can do several passes on it, as required to read the GML format
   efficiently
 */
int fill_buffer(FILE *stream) {
    int length;
    char line[LINELENGTH];
    LINE *previous;

    if (read_line(stream,line)!=0) {
        first = NULL;                /* Indicates empty buffer */
        return 1;
    }
    length = strlen(line) + 1;
    first = malloc(sizeof(LINE));
    first->str = malloc(length*sizeof(char));
    strcpy(first->str,line);

    previous = first;
    while (read_line(stream,line)==0) {
        length = strlen(line) + 1;
        previous->ptr = malloc(sizeof(LINE));
        previous = previous->ptr;
        previous->str = malloc(length*sizeof(char));
        strcpy(previous->str,line);
    }
    previous->ptr = NULL;          /* Indicates last line */

    return 0;
}


/* Function to free up the buffer again */
void free_buffer(void) {
    LINE *thisptr;
    LINE *nextptr;

    thisptr = first;
    while (thisptr!=NULL) {
        nextptr = thisptr->ptr;
        free(thisptr->str);
        free(thisptr);
        thisptr = nextptr;
    }
}


/* Function to reset to the start of the buffer again */
void reset_buffer(void) {
    current = first;
}


/* Function to get the next line in the buffer.  Returns 0 if there was
   a line or 1 if we've reached the end of the buffer.
 */
int next_line(char line[LINELENGTH]) {
    if (current==NULL) return 1;
    strcpy(line,current->str);
    current = current->ptr;
    return 0;
}



/* Function to establish whether the network read from a given stream is
   directed or not.  Returns 1 for a directed network, and 0 otherwise.  If
   the GML file contains no "directed" line then the graph is assumed to be
   undirected, which is the GML default behavior.
 */
int is_directed(void) {
    int result=0;
    char *ptr;
    char line[LINELENGTH];

    reset_buffer();

    while (next_line(line)==0) {
        ptr = strstr(line,"directed");
        if (ptr==NULL) continue;
        sscanf(ptr,"directed %i",&result);
        break;
    }
    result=0;		/* Making it undirected for all cases; */
    return result;
}



/* Function to count the vertices in a GML file.  Returns number of vertices.
 */
int count_vertices(void) {
    int result=0;
    char *ptr;
    char line[LINELENGTH];

    reset_buffer();

    while (next_line(line)==0) {
        ptr = strstr(line,"node");
        if (ptr!=NULL) result++;
    }
    return result;
}


/* Function to compare the IDs of two vertices */

int cmpid(const void *v1p, const void *v2p) {
    const VERTEX *v1, *v2;
    v1 = (const VERTEX *) v1p;
    v2 = (const VERTEX *) v2p;
    if (v1->id>v2->id) return 1;
    if (v1->id<v2->id) return -1;
    return 0;
}


/* Function to allocate space for a network structure stored in a GML file
   and determine the parameters (id, label) of each of the vertices.
 */

void create_network(NETWORK *network) {
    int i;
    int length;
    char *ptr;
    char *start,*stop;
    char line[LINELENGTH];
    char *label;
    int temp;

    /* Determine whether the network is directed */

    network->directed = is_directed();

    /* Count the vertices */

    network->nvertices = count_vertices();

    /* Make space for the vertices */

    network->vertex = calloc(network->nvertices,sizeof(VERTEX));

    label = (char *) malloc(LINELENGTH * sizeof(char));

    /* Go through the file reading the details of each vertex one by one */

    reset_buffer();
    temp=0; 
    for (i=0; i<network->nvertices; i++) {
        temp++;
        /* 
           printf("Reading node %d\n",temp);
           Skip to next "node" entry
         */
        do {
            next_line(line);
        } while (strstr(line,"node")==NULL);

        /* Read in the details of this vertex */

        do {

            /* Look for ID */

            ptr = strstr(line,"id");
            if (ptr!=NULL) sscanf(ptr,"id %i",&network->vertex[i].id);

            /* Look for label */

            ptr = (strstr(line,"label"));
            if (ptr!=NULL) {
                start = strchr(line,'"');
                if (start==NULL) {
                    sscanf(ptr,"label %s",label);
                } else {
                    stop = strchr(++start,'"');
                    if (stop==NULL) length = strlen(line) - (start-line);
                    else length = stop - start;
                    strncpy(label,start,length);
                    label[length] = '\0';
                    network->vertex[i].label = malloc((length+1)*sizeof(char));
                    strcpy(network->vertex[i].label,label);
                }
            }

            /* If we see a closing square bracket we are done */

            if (strstr(line,"]")!=NULL) break;

        } while (next_line(line)==0);

    }

    /* Sort the vertices in increasing order of their IDs so we can find them
       quickly later
     */
    qsort(network->vertex,network->nvertices,sizeof(VERTEX), cmpid);
    free(label);
}

/*
   Function to find a vertex with a specified ID using binary search.
   Returns the element in the vertex[] array holding the vertex in question,
   or -1 if no vertex was found.
 */
int find_vertex(int id, NETWORK *network) {
    int top,bottom,split;
    int idsplit;

    top = network->nvertices;
    if (top<1) return -1;
    bottom = 0;
    split = top/2;

    do {
        idsplit = network->vertex[split].id;
        if (id>idsplit) {
            bottom = split + 1;
            split = (top+bottom)/2;
        } else if (id<idsplit) {
            top = split;
            split = (top+bottom)/2;
        } else return split;
    } while (top>bottom);

    return -1;
}

/*
   Function to determine the degrees of all the vertices by going through
   the edge data
 */
void get_degrees(NETWORK *network) {
    int s,t;
    int vs,vt;
    char *ptr;
    char line[LINELENGTH];

    reset_buffer();

    while (next_line(line)==0) {

        /* Find the next edge entry */

        ptr = strstr(line,"edge");
        if (ptr==NULL) continue;

        /* Read the source and target of the edge */

        s = t = -1;

        do {

            ptr = strstr(line,"source");
            if (ptr!=NULL) sscanf(ptr,"source %i",&s);
            ptr = strstr(line,"target");
            if (ptr!=NULL) sscanf(ptr,"target %i",&t);

            /* If we see a closing square bracket we are done */

            if (strstr(line,"]")!=NULL) break;

        } while (next_line(line)==0);

        /* Increment the degrees of the appropriate vertex or vertices */

        if ((s>=0)&&(t>=0)) {
            vs = find_vertex(s,network);
            network->vertex[vs].degree++;
            if (network->directed==0) {
                vt = find_vertex(t,network);
                network->vertex[vt].degree++;
            }
        }

    }

    return;
}


/* Function to read in the edges */
void read_edges(NETWORK *network) {
    int i;
    int s,t;
    int vs,vt;
    int *count;
    double w;
    char *ptr;
    char line[LINELENGTH];
    int temp;
    /* Malloc space for the edges and temporary space for the edge counts
       at each vertex
     */
    for (i=0; i<network->nvertices; i++) {
        network->vertex[i].edge = malloc(network->vertex[i].degree*sizeof(EDGE));
    }
    count = calloc(network->nvertices,sizeof(int));

    /* Read in the data */

    reset_buffer();
    temp = 0;
    while (next_line(line)==0) {
        temp++;
        /* Find the next edge entry */

        ptr = strstr(line,"edge");
        if (ptr==NULL) continue;

        /* Read the source and target of the edge and the edge weight */

        s = t = -1;
        w = 1.0;

        do {

            ptr = strstr(line,"source");
            if (ptr!=NULL) sscanf(ptr,"source %i",&s);
            ptr = strstr(line,"target");
            if (ptr!=NULL) sscanf(ptr,"target %i",&t);
            ptr = strstr(line,"value");
            if (ptr!=NULL) sscanf(ptr,"value %lf",&w);

            /* If we see a closing square bracket we are done */

            if (strstr(line,"]")!=NULL) break;

        } while (next_line(line)==0);

        /* Add these edges to the appropriate vertices */

        if ((s>=0)&&(t>=0)) {
            vs = find_vertex(s,network);
            vt = find_vertex(t,network);
            network->vertex[vs].edge[count[vs]].target = vt;
            network->vertex[vs].edge[count[vs]].weight = w;
            network->vertex[vs].edge[count[vs]].id = temp;
            count[vs]++;
            if (network->directed==0) {
                network->vertex[vt].edge[count[vt]].target = vs;
                network->vertex[vt].edge[count[vt]].weight = w;
                network->vertex[vt].edge[count[vt]].id = temp;
                count[vt]++;
            }
        }

    }

    free(count);
    return;
}


/* Function to read a complete network */
int read_network(NETWORK *network, FILE *stream) {
    fill_buffer(stream);
    create_network(network);

    get_degrees(network);

    read_edges(network);
    free_buffer();

    return 0;
}


/* Function to free the memory used by a network again */
void free_network(NETWORK *network) {
    int i;
    for (i=0; i<network->nvertices; i++) {
        free(network->vertex[i].edge);
        free(network->vertex[i].label);
    }
    free(network->vertex);
}

void network_to_graph(graph_t *G, NETWORK *N) {
    int i,numEdges;
    VERTEX vertex;
    int j,degree,target,sumDegree;
    double weight;
    attr_id_t start;
    int count;

    G->n = N->nvertices;
    numEdges = 0;
    sumDegree = 0;
    for(i=0; i<N->nvertices; i++)
    {
        vertex = N->vertex[i];
        numEdges+=vertex.degree;
    }
    assert(numEdges%2==0);
    G->m = numEdges;
    assert(G->n > 0);
    assert(G->m > 0);
    G->undirected=1;	/* Hard-coded undirected */
    G->zero_indexed=0;	/* Hard-coded false */
    G->weight_type=4;	/* Hard-coded weight_type is double */
    /* Allocating memory */
    G->numEdges = (attr_id_t*) calloc(G->n+1, sizeof(attr_id_t));
    G->endV = (attr_id_t*) calloc(G->m, sizeof(attr_id_t));
    G->edge_id = (attr_id_t *) calloc(G->m, sizeof(attr_id_t));
    G->dbl_weight_e = (double*) malloc(sizeof(double)* G->m);

    assert(G->numEdges != NULL);
    assert(G->endV != NULL);
    assert(G->dbl_weight_e != NULL);

    G->numEdges[0]=0;
    count=0;
    for(i=0; i<N->nvertices; i++)
    {
        vertex = N->vertex[i];
        degree = vertex.degree;
        G->numEdges[i+1] = G->numEdges[i] + degree;
        start = G->numEdges[i];

        for(j=0; j<degree;j++)
        {
            target = vertex.edge[j].target;
            weight = vertex.edge[j].weight;
            G->edge_id[start+j] = vertex.edge[j].id;
            G->endV[start+j] = target;
            G->dbl_weight_e[start+j] = weight;
            count++;
        }
    }

    free_network(N);
    free(N);

}

void read_GML_graph(graph_t* G, char* filename) {
    NETWORK *N = (NETWORK*)malloc(sizeof(NETWORK));
    FILE *ifile = fopen(filename, "r");
    read_network(N,ifile);	/* Note that readgml is hardcoded to read 
                               undirected graphs,even if the input graph 
                               is directed.
                             */
    network_to_graph(G,N);
    fclose(ifile);
}	

