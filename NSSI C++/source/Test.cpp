#include "su4.h"
#include <iostream>
#include <algorithm>

#include "Constants.h"
#include "Noise.h"

class Class
{
    public:
        Class() : v {}
        {
            std::cout << this << ": default constructor"  << std::endl;
        }
        Class(std::vector<int> vv) : v {std::move(vv)}
        {
            std::cout << this << ": arg constructor"  << std::endl;
        }
        Class(const Class& obj) : v {obj.v}
        {
            std::cout << this << ": copy constructor"  << std::endl;
        }
        Class(Class&& obj) : v {std::move(obj.v)}
        {
            std::cout << this << ": move constructor"  << std::endl;
        }
        Class operator+(const Class& obj) // Elementwise summation
        {   
            if (obj.v.size() == v.size())
            {
                Class res;
                for (int i = 0; i < v.size(); ++i) res.v.push_back(v[i] + obj.v[i]);
                return res;
            }
            else
            {
                std::cout << "The vectors have different sizes!";
                return Class();
            }
        }
        void print() // Output of the vector v
        {
            for (auto it = v.begin(); it != v.end(); ++it) std::cout << *it << " ";
            std::cout << std::endl;
        }

    private:
        std::vector<int> v;
};

int main()
{
    /*
    // Cook test #1
    std::vector<int> arr(10);
    for (auto it = arr.begin(); it != arr.end(); ++it) *it = 5;

    Class o1(arr);
    Class o2(o1);
    Class o3 = o1 + o2;

    o1.print();
    o2.print();
    o3.print();
    */

    /*
    //Test: Constants class
    Constants c;
    c.gather("./ParametersFancy.json");
    c.output();
    c.dump("./Parameters.json");
    */

    /*
    su4::Vector vec1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
    su4::Vector vec2({5, 5, 5, 5, 5, 5, 5, 5, 5, 6,  6,  6,  6,  6,  6});

    su4::Vector vecP   = vec1 + vec2; vecP.output();
    su4::Vector vecM   = vec1 - vec2; vecM.output();
    su4::Vector vecPP  = vec1 * 2.0;  vecPP.output();
    su4::Vector vecPPR = 2.0 * vec1;  vecPPR.output();
    su4::Vector vecD   = vec1 / 2.0;  vecD.output();
    */

   /*
   su4::Vector vec;
   vec = su4::Vector({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
   vec.output();
   */

    Noise n;
    n.gather("./Noise.json");

    return 0;
}