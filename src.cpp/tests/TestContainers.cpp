#include <unistd.h>
#include <chrono>
#include <omp.h>
#include <fstream>

#include "Containers.h"

/*
    Tested 10th March
    Seems eveything's OK. Acceleration in the test here is about 3.2 on 4 cores
*/

void Test_Constructors()
{
    TwoLines tl1;
    TwoLines tl2(4);

    std::vector<Matrix> left3(4);
    std::vector<Matrix> right3(4);
    for (int i = 0; i < 4; ++i)
    {
        left3[i]  = Matrix::Identity() * (i+1);
        right3[i] = Matrix::Identity() * (5-i);
    }

    TwoLines tl3(left3, right3, 4);
    TwoLines tl4 = tl3;
    TwoLines tl5; tl5 = tl3;

    std::cout << "tl1 >>>" << std::endl << tl1 << std::endl << std::endl
              << "tl2 >>> "<< std::endl << tl2 << std::endl << std::endl
              << "tl3 >>> "<< std::endl << tl3 << std::endl << std::endl
              << "tl4 >>> "<< std::endl << tl4 << std::endl << std::endl
              << "tl5 >>> "<< std::endl << tl5 << std::endl << std::endl;

    tl1.init(left3, right3, 4);
    std::cout << "tl1.init >>> "<< std::endl << tl1 << std::endl << std::endl;
}

void Test_Mod()
{
    int N = 20;

    for (int i = -3*N; i < 3*N; ++i)
    {
        std::cout << i << " % "  << N  << " = "<< ( (i % N) + N) % N << std::endl;
    }
}

void Test_Operations()
{
    int N = 20;

    std::vector<Matrix> left6(N);
    std::vector<Matrix> right6(N);
    for (int i = 0; i < N; ++i)
    {
        left6[i]  = Matrix::Identity() * (i+2);
        right6[i] = Matrix::Identity() * (7-i);
    }

    std::vector<Matrix> left7(N);
    std::vector<Matrix> right7(N);
    for (int i = 0; i < N; ++i)
    {
        left7[i]  = Matrix::Identity() * (i+1);
        right7[i] = Matrix::Identity() * (5-i);
    }

    TwoLines tl6(left6, right6, N);
    TwoLines tl7(left7, right7, N);
    
    std::cout << "tl6 >>>" << std::endl << tl6 << std::endl << std::endl
              << "tl7 >>> "<< std::endl << tl7 << std::endl << std::endl
              << "tl6 + tl7 >>> "<< std::endl << tl6 + tl7 << std::endl << std::endl
              << "tl6 - tl7 >>> "<< std::endl << tl6 - tl7 << std::endl << std::endl
              << "tl6 * 2.0 >>> "<< std::endl << tl6 * 2.0 << std::endl << std::endl
              << "2.0 * tl6 >>> "<< std::endl << 2.0 * tl6 << std::endl << std::endl
              << "tl6 / 2.0 >>> "<< std::endl << tl6 / 2.0 << std::endl << std::endl;

    std::cout << "t6.left >>>" << std::endl;
    for (int i = -N; i < 2*N; ++i)
    {
        std::cout << "[" << i << "]:" << std::endl << tl6.left(i) << std::endl << std::endl;
    }
    std::cout << "t6.right >>>" << std::endl;
    for (int i = -N; i < 2*N; ++i)
    {
        std::cout << "[" << i << "]:" << std::endl << tl6.right(i) << std::endl << std::endl;
    }

    std::cout << "tl6.derivative(order=2) >>>" << std::endl << tl6.centralDerivative(1.0, 2) << std::endl << std::endl;
    std::cout << "tl6.derivative(order=4) >>>" << std::endl << tl6.centralDerivative(1.0, 4) << std::endl << std::endl;
}

double Test_Acc(int N)
{
    // Let's make a real grid of coordinates
    Real X  = 2.0*3.1415;
    Real dx = X / N; Real x;
    std::vector<Matrix> left(N);
    std::vector<Matrix> right(N);
    // Initialise the matrices
    for (int k = 0; k < N; ++k)
    {
        x = k * dx;
        left[k] =  Matrix::Identity() * std::sin(x);
        right[k] = Matrix::Identity() * std::cos(x);
    }
    // Container
    TwoLines tl(left, right, N);

    // Calculations
    auto start = std::chrono::high_resolution_clock::now();

    TwoLines der = tl.centralDerivative(dx, 2);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    // Just a formal output
    std::ofstream output("testSins.txt");
    for (int i = 0; i < N; ++i)
    {
        output << der.right(i).trace().real() / 4.0 << " ";
    }
    output.close();

    return elapsed.count();
}

struct Kek
{
    std::vector<Real> vec1;
    std::vector<Real> vec2;
};

int main(int argc, char *argv[])
{
    // Finally, the test of the acceleration

    char opt; int threads=1;
    while ( (opt = getopt(argc, argv, "t:")) != -1)
    {
        if (opt == 't')
        {
            threads = atoi(optarg);
        }
    }

    omp_set_num_threads(threads);
    std::cout << "threads=" << threads << std::endl;
    std::cout << "Elapsed time: " << Test_Acc(10000) << " s" << std::endl;

    Kek kek = {{1, 2, 3}, {4, 5, 6}};
    Kek heh  = kek;

    for (auto it = kek.vec1.begin(); it != kek.vec1.end(); ++it)
    {
        std::cout << *it << " ";
    }

    return 0;
}