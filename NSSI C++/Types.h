#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <array>

#include "Dense"

using Real    = double;
using Complex = std::complex<Real>;
using Matrix  = Eigen::Matrix<Complex, 4, 4, Eigen::RowMajor>;
using String  = std::string;
using Array15 = std::array<Real,15>;
using namespace std::complex_literals;

// long double complex zero
const Complex ComplexZero = 0.0;