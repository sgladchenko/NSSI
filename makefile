Compiler    = g++

Sources =  ./NSSI C++/source/
Tests   =  ./NSSI C++/tests/
Headers =  -I "./NSSI C++/" -I "./NSSI C++/Eigen/" -I "./NSSI C++/JSON/single_include/"

# The main part

Main:
	$(Compiler) $(Headers) -o nssi "$(Sources)Main.cpp" \
								   "$(Sources)su4.cpp" \
								   "$(Sources)Vector.cpp" \
								   "$(Sources)Containers.cpp" \
								   "$(Sources)Constants.cpp" \
								   "$(Sources)Files.cpp" \
								   "$(Sources)Functions.cpp" \
								   "$(Sources)Hamiltonians.cpp" \
								   "$(Sources)Noise.cpp" \
								   "$(Sources)Scheme.cpp" \
						   		   -std=c++17 \
								   -fopenmp

# Test section

TestConstants:
	$(Compiler) $(Headers) -o test "$(Tests)TestConstants.cpp" \
								   "$(Sources)Constants.cpp" \
								   "$(Sources)su4.cpp" \
								   "$(Sources)Vector.cpp" \
						   		   -std=c++17

TestVector:
	$(Compiler) $(Headers) -o test "$(Tests)TestVector.cpp" \
								   "$(Sources)su4.cpp" \
								   "$(Sources)Vector.cpp" \
						   		   -std=c++17

Testsu4:
	$(Compiler) $(Headers) -o test "$(Tests)Testsu4.cpp" \
								   "$(Sources)su4.cpp" \
								   "$(Sources)Vector.cpp" \
						   		   -std=c++17

TestContainers:
	$(Compiler) $(Headers) -o test "$(Tests)TestContainers.cpp" \
								   "$(Sources)su4.cpp" \
								   "$(Sources)Vector.cpp" \
								   "$(Sources)Containers.cpp" \
						   		   -std=c++17 \
								   -fopenmp

TestFiles:
	$(Compiler) $(Headers) -o test "$(Tests)TestFiles.cpp" \
								   "$(Sources)Files.cpp" \
						   		   -std=c++17 \
								   -fopenmp
