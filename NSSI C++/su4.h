#pragma once

#include "Types.h"
#include <iostream>

namespace su4
{
	// An additional structure needed to save the 15D real-valued vectors
	class Vector
	{
		public:
			// The default one; all the entries are set to zeros
			Vector();
			// Copy from the vector
			Vector(Array15 stdvec);
			// Copy constructor
			Vector(const Vector& vec);
			// Move constructor
			Vector(Vector&& vec);

			// Copy assignment operator
			Vector& operator=(const Vector& vec);
			// Move assignment operator
			Vector& operator=(Vector&& vec);

			// Access method
			const Real& operator[](int k) const;

			// Basic arithmetic operations
			Vector operator+(const Vector& vec) const;
			Vector operator-(const Vector& vec) const;
			Vector operator*(Real z) const;
			Vector operator/(Real z) const;

			Vector& operator+=(const Vector& vec);

			// Some output
			void output()
			{
				for (auto it = arr.begin(); it != arr.end(); ++it) std::cout << *it << " ";
				std::cout << std::endl;
			}
		private:
			Array15 arr;
	};

	// The operations that I'm going to use in all the parts of my code
	Vector Decompose(const Matrix& m);
	Matrix Compose(const Vector& v);
	Matrix Normalise(const Matrix& m, Real initialNorm);
	Real   Norm(const Matrix& m);

	Matrix Diag(const Matrix& m);
	Matrix Offdiag(const Matrix& m);
	Matrix Round(const Matrix& m);

	Matrix MakeDiagonal(const std::vector<Real>& v);
	Matrix MakeDiagonal(Real a, Real b, Real c, Real d);

	// My lovely su(4) basis in the Lie algebra
	// used in different operations over the density matrices
	// First, I need to define the following constants matrices

	inline const Complex lambda1[16] = {0.0, 1.0, 0.0, 0.0,
								        1.0, 0.0, 0.0, 0.0,
								        0.0, 0.0, 0.0, 0.0,
								        0.0, 0.0, 0.0, 0.0};

	inline const Complex lambda2[16] = {0.0, -1.0i, 0.0, 0.0,
								        1.0i, 0.0,  0.0, 0.0,
								        0.0,  0.0,  0.0, 0.0,
								        0.0,  0.0,  0.0, 0.0};

	inline const Complex lambda3[16] = {1.0, 0.0, 0.0,  0.0,
								        0.0, -1.0, 0.0, 0.0,
								        0.0,  0.0, 0.0, 0.0,
								        0.0,  0.0, 0.0, 0.0};

	inline const Complex lambda4[16] = {0.0, 0.0, 1.0, 0.0,
								        0.0, 0.0, 0.0, 0.0,
								        1.0, 0.0, 0.0, 0.0,
								        0.0, 0.0, 0.0, 0.0};

	inline const Complex lambda5[16] = {0.0,  0.0, -1.0i, 0.0,
								        0.0,  0.0,   0.0, 0.0,
								        1.0i, 0.0,   0.0, 0.0,
								        0.0,  0.0,   0.0, 0.0};

	inline const Complex lambda6[16] = {0.0, 0.0, 0.0, 0.0,
								        0.0, 0.0, 1.0, 0.0,
								        0.0, 1.0, 0.0, 0.0,
								        0.0, 0.0, 0.0, 0.0};

	inline const Complex lambda7[16] = {0.0,  0.0,   0.0, 0.0,
								        0.0,  0.0, -1.0i, 0.0,
								        0.0, 1.0i,   0.0, 0.0,
								        0.0,  0.0,   0.0, 0.0};

	inline const Complex lambda8[16] = {1.0/std::sqrt(3),  0.0, 0.0, 0.0,
								        0.0, 1.0/std::sqrt(3),  0.0, 0.0,
								        0.0, 0.0, -2.0/std::sqrt(3), 0.0,
								        0.0, 0.0,               0.0, 0.0};

	inline const Complex lambda9[16] = {0.0, 0.0, 0.0, 1.0,
								        0.0, 0.0, 0.0, 0.0,
								        0.0, 0.0, 0.0, 0.0,
								        1.0, 0.0, 0.0, 0.0};

	inline const Complex lambda10[16] = {0.0,  0.0, 0.0, -1.0i,
								         0.0,  0.0, 0.0,   0.0,
								         0.0,  0.0, 0.0,   0.0,
								         1.0i, 0.0, 0.0,   0.0};

	inline const Complex lambda11[16] = {0.0, 0.0, 0.0, 0.0,
								         0.0, 0.0, 0.0, 1.0,
								         0.0, 0.0, 0.0, 0.0,
								        0.0, 1.0, 0.0, 0.0};

	inline const Complex lambda12[16] = {0.0,  0.0, 0.0,   0.0,
								         0.0,  0.0, 0.0, -1.0i,
								         0.0,  0.0, 0.0,   0.0,
								         0.0, 1.0i, 0.0,   0.0};

	inline const Complex lambda13[16] = {0.0, 0.0, 0.0, 0.0,
								         0.0, 0.0, 0.0, 0.0,
								         0.0, 0.0, 0.0, 1.0,
								         0.0, 0.0, 1.0, 0.0};

	inline const Complex lambda14[16] = {0.0, 0.0,  0.0,   0.0,
								         0.0, 0.0,  0.0,   0.0,
								         0.0, 0.0,  0.0, -1.0i,
								         0.0, 0.0, 1.0i,   0.0};

	inline const Complex lambda15[16] = {1.0/std::sqrt(6),  0.0, 0.0, 0.0,
								         0.0, 1.0/std::sqrt(6),  0.0, 0.0,
								         0.0, 0.0, 1.0/std::sqrt(6),  0.0,
								         0.0, 0.0, 0.0, -3.0/std::sqrt(6)};

	// Then the Eigen::Matrix objects themselves
	// Here I also normalise them to ||.||=1.0
	inline const Matrix Basis[15] = {Matrix(lambda1)/std::sqrt(2),  Matrix(lambda2)/std::sqrt(2),
						             Matrix(lambda3)/std::sqrt(2),  Matrix(lambda4)/std::sqrt(2),
						             Matrix(lambda5)/std::sqrt(2),  Matrix(lambda6)/std::sqrt(2),
						             Matrix(lambda7)/std::sqrt(2),  Matrix(lambda8)/std::sqrt(2),
						             Matrix(lambda9)/std::sqrt(2),  Matrix(lambda10)/std::sqrt(2),
						             Matrix(lambda11)/std::sqrt(2), Matrix(lambda12)/std::sqrt(2),
						             Matrix(lambda13)/std::sqrt(2), Matrix(lambda14)/std::sqrt(2),
						             Matrix(lambda15)/std::sqrt(2)};

	// And some useful matrices engaged in the noncollective terms
	inline const Complex AMMLike_arr[16] = {0.0,    0.0,   0.0,  -1.0i,
					  					    0.0,    0.0,  1.0i,    0.0,
					  					    0.0,  -1.0i,   0.0,    0.0,
					 					    1.0i,   0.0,   0.0,    0.0};

	inline const Complex flavourSigma1_arr[16] = {0.0, 1.0, 0.0, 0.0,
						                          1.0, 0.0, 0.0, 0.0,
						    				      0.0, 0.0, 0.0, 1.0,
						    				      0.0, 0.0, 1.0, 0.0};

	inline const Complex flavourSigma3_arr[16] = {1.0,  0.0, 0.0,  0.0,
											      0.0, -1.0, 0.0,  0.0,
											      0.0,  0.0, 1.0,  0.0,
											      0.0,  0.0, 0.0, -1.0};

	// And the Matrix instances
	inline const Matrix AMMLike(AMMLike_arr);
	inline const Matrix flavourSigma1(flavourSigma1_arr);
	inline const Matrix flavourSigma3(flavourSigma3_arr);
}

// A multiplication of a vector on a Complex, but the vector is on the RHS
su4::Vector operator*(Real z, const su4::Vector& vec);