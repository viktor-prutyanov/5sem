/**
 *  Task 0
 *
 *  @file main.c
 *  
 *  @date 11.2016
 * 
 *  @copyright GNU GPL v2.0
 * 
 *  @author Viktor Prutyanov mailto:viktor.prutyanov@phystech.edu 
 *  
 */

#include <stdio.h>
#include <pthread.h>

#define NT 3

pthread_mutex_t mutex;
pthread_cond_t cond;

void *threadfn(void *data)
{
    int x = *(int *)data;
    int y = x * x;
    return (void *)y;
}

int main()
{
    pthread_t tid[NT];
    int x[NT] = {1, 2, 3}, x2;

    for (int i = 0; i < NT; i++) 
    {
        x2 = x[i] * 2;
        pthread_create(&tid[i], NULL, threadfn, (void*)&x2);
    }

    for (int i = 0; i < NT; i++)
    {
        int result;
        pthread_join(tid[i], (void **)&result);
        printf("Result of thread %d: %d\n", i, result);
    }

    return 0;
}
