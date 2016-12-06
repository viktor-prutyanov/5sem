#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

long int get_pos_long_opt(char *optarg, const char *name)
{
    char *endptr;
    errno = 0;
    
    long int arg = strtol(optarg, &endptr, 10);

    if (errno == ERANGE)
    {
        fprintf(stderr, "Args error: %s is out of range.\n", name);
    }
    else if (*endptr != '\0')
    {
        fprintf(stderr, "Args error: invalid %s.\n", name);
        errno = EINVAL;
    }
    else if (arg <= 0)
    {
        fprintf(stderr, "Args error: nonpositive %s.\n", name);
        errno = EINVAL;
    }
    
    return arg;
}

double get_pos_double_opt(char *optarg, const char *name)
{
    char *endptr;
    errno = 0;
    
    double arg = strtod(optarg, &endptr);
    
    if (errno == ERANGE)
    {
        fprintf(stderr, "Args error: %s is out of range.\n", name);
    }
    else if (*endptr != '\0')
    {
        fprintf(stderr, "Args error: invalid %s.\n", name);
        errno = EINVAL;
    }
    else if (arg <= 0)
    {
        fprintf(stderr, "Args error: nonpositive %s.\n", name);
        errno = EINVAL;
    }

    return arg;    
}
