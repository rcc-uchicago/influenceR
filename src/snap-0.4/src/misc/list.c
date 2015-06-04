#include "utils.h"

list_t* makeList()
{
    list_t *l  = (list_t*)malloc(sizeof(list_t));
    l->head=NULL;
    l->tail=NULL;
    l->size=0;
    return l;
}

node_t* makeNode(int id)
{
    node_t *newnode = (node_t*)malloc(sizeof(node_t));
    newnode->id = id;
    newnode->next = NULL;
    return newnode;
}

/* note append will append at the last. */
void append(list_t *L, node_t *n)
{
    if(L->size==0)
    {
        L->head=n;
        L->tail=n;
    }
    else
    {
        L->tail->next = n;
        L->tail = n;
    }
    L->size++;

}

node_t* getFirst(list_t *L)
{
    return L->head;
}

void deleteFirst(list_t *L)
{
    node_t *n = L->head;
    L->head = n->next;
    free(n);
    L->size--;
}

void printList(list_t *L)
{
    node_t *n;
    printf("Printing list of size:%d\n",L->size);
    if(L->size==0) return;
    for(n=L->head; n!=L->tail; n=n->next)
    {
        printf("%d,",n->id);
    }
    printf("%d,",n->id);
    printf("\n\n\n\n");
}
