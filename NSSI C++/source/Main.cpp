#include <unistd.h>
#include <chrono>
#include <omp.h>

#include "Scheme.h"
#include "Containers.h"
#include "Files.h"

int main(int argc, char* argv[])
{
    // Obtain a number of the threads;
    // The default one is 1
    char opt; int Threads=1;
    while ( (opt = getopt(argc, argv, "t:")) != -1)
    {
        if (opt == 't')
        {
            Threads = atoi(optarg);
        }
    }
    // Set obtained number of threads
    omp_set_num_threads(Threads);

    return 0;
}