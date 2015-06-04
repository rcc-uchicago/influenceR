#include "utils.h"

void printDoubleVector(double *arr,int start, int end)
{
    int i;
    for (i=start; i<end; i++)
    {
        printf("%g", arr[i]);
        if(i != end-1)
            printf(" ");
    }
    printf("\n");
}

void printIntVector(int *arr,int start, int end)
{
    int i;
    for (i=start; i<end; i++)
    {
        printf("%d", arr[i]);
        if(i != end-1)
            printf(" ");
    }
    printf("\n");
}

void print_attr_id_t_Vector(attr_id_t *arr,attr_id_t start, attr_id_t end)
{
    attr_id_t i;
    for (i=start; i<end; i++)
    {
        printf("%d", arr[i]);
        if(i != end-1)
            printf(" ");
    }
    printf("\n");
}


