Compiler    = g++

Sources =  ./NSSI C++/source/
Tests   =  ./NSSI C++/tests/
Headers =  -I "./NSSI C++/" -I "./NSSI C++/Eigen/" -I "./NSSI C++/JSON/single_include/"

test:
	$(Compiler) $(Headers) -o test "$(Tests)TestConstants.cpp" \
								   "$(Sources)Constants.cpp" \
								   "$(Sources)su4.cpp" \
								   "$(Sources)Vector.cpp" \
						   		   -std=c++17