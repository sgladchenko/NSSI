#include <getopt.h>
#include <chrono>
#include <omp.h>

#include "Scheme.h"
#include "Files.h"

int main(int argc, char* argv[])
{
    // Obtain a number of the threads;
    // The default one is 1
    struct option longoptions[] = 
    {
        {"numthreads", required_argument, NULL, 't'},
        {"periodN_x",  required_argument, NULL, 'x'},
        {"periodN_z",  required_argument, NULL, 'z'},
        {"root",       required_argument, NULL, 'r'},
        {"parameters", required_argument, NULL, 'p'},
        {"noise",      required_argument, NULL, 'n'}
    };

    char opt;
    // The parameters given as arguements in the command line
    int Threads=1;
    int periodN_x=1;
    int periodN_z=1;
    String root = "./Data/";
    String parameters = "./Parameters.json";
    String noise = "./Noise.json";

    while ((opt = getopt_long(argc, argv, "", longoptions, NULL)) != -1)
    {
        switch (opt)
        {
            case 't': Threads = atoi(optarg); break;
            case 'x': periodN_x = atoi(optarg); break;
            case 'z': periodN_z = atoi(optarg); break;
            case 'r': root = optarg; break;
            case 'p': parameters = optarg; break;
            case 'n': noise = optarg; break;
        }
    }

    // Check if the root ends with /
    if (root[root.length()-1] != '/')
    {
        root += "/";
    }

    // Set obtained number of threads
    omp_set_num_threads(Threads);
    
    std::cout << "Set number of threads: " << Threads << std::endl;
    std::cout << "Period of the displayed x grid: " << periodN_x << std::endl;
    std::cout << "Period of the displayed z grid: " << periodN_z << std::endl;
    std::cout << "Set root directory: " << root << std::endl;

    // Time marks
    auto start = std::chrono::high_resolution_clock::now();

    Scheme Setup(parameters, noise, root, periodN_x, periodN_z);
    Setup.Solve();

    // Time marks
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    std::cout << "Calculation completed" << std::endl;
    std::cout << "Elapsed time: " << elapsed.count() << " sec" << std::endl;
    std::cout << "Location: " << std::endl << Setup.location() << std::endl;

    return 0;
}