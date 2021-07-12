#include "su4.h"

// The default one that allocates the memory and actually sets all the entries to zeros
su4::Vector::Vector() : arr {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} {} // Omg...

// The constructor that grabs the data from a std::array<Real,15>
su4::Vector::Vector(Array15 stdArray) : arr {std::move(stdArray)} {}

// Copy constructor
su4::Vector::Vector(const su4::Vector& vec) : arr {vec.arr} {}

// Move constructor
su4::Vector::Vector(su4::Vector&& vec) : arr {std::move(vec.arr)} {}

// Copy assignment operator
su4::Vector& su4::Vector::operator=(const su4::Vector& vec)
{
    arr = vec.arr;
    return *this;
}

// Move assignment operator
su4::Vector& su4::Vector::operator=(su4::Vector&& vec)
{
    arr = std::move(vec.arr);
    return *this;
}

// Access method
const Real& su4::Vector::operator[](int k) const
{
    return arr[k];
}

// Basic arithmetic operations

// Elementwise summation
su4::Vector su4::Vector::operator+(const su4::Vector& vec) const
{
    Vector tmp;
    for (int i = 0; i < 15; ++i)
    {
        tmp.arr[i] = arr[i] + vec.arr[i];
    }
    return tmp;
}

// Elementwise subtraction
su4::Vector su4::Vector::operator-(const su4::Vector& vec) const
{
    Vector tmp;
    for (int i = 0; i < 15; ++i)
    {
        tmp.arr[i] = arr[i] - vec.arr[i];
    }
    return tmp;
}

// Elementwise multiplication by a constant
su4::Vector su4::Vector::operator*(Real z) const
{
    Vector tmp;
    for (int i = 0; i < 15; ++i)
    {
        tmp.arr[i] = arr[i]*z;
    }
    return tmp;
}

// Elementwise division by a constant
su4::Vector su4::Vector::operator/(Real z) const
{
    Vector tmp;
    for (int i = 0; i < 15; ++i)
    {
        tmp.arr[i] = arr[i]/z;
    }
    return tmp;
}

// Add elementwise another vector to this one
su4::Vector& su4::Vector::operator+=(const su4::Vector& vec)
{
    for (int i = 0; i < 15; ++i)
    {
        arr[i] += vec.arr[i];
    }
    return *this;
}

// Multiplication on a real constant, while vector is in the RHS
su4::Vector operator*(Real z, const su4::Vector& vec)
{
    return vec * z;
}

// An output in the STL-style
std::ostream& operator<<(std::ostream& os, const su4::Vector& vec)
{
    os << "(";
    for (int i = 0; i < 14; ++i)
    {
        os << vec[i] << ", ";
    }
    os << vec[14] << ")";
    return os;
}