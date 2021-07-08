// Undefine to use C++17 std::filesystem
// #define __DONT_USE_FILESYSTEM__

#include "Files.h"

#include <chrono>
#include <ctime> 

#include <iostream>
#include <fstream>

#ifndef __DONT_USE_STD_FILESYSTEM__

	#include <filesystem>
	namespace fs = std::filesystem;
	
#else // Add support for older gcc compilers/linkers [added by O.K., 20Jun2021]

	#include <sstream>
	
	#if !defined(__LINUX__) && !defined(__linux__) && !defined(__gnu_linux__)
		#include <direct.h>
	#else
		#include <sys/stat.h>
		inline int _mkdir(const char* dirName)
		{
			return mkdir(dirName, S_IRUSR | S_IWUSR | S_IXUSR);
		}
	#endif
	
	template <typename T>
	std::string ToString(const T& x)
	{
		std::ostringstream os;
		os << x;
		return os.str();
	}
	
	inline bool DirectoryExists(const std::string& fileName)
	{
		#if defined(__LINUX__) || defined(__linux__) || defined(__gnu_linux__)
			struct stat sb;
			return stat(fileName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode);
		#else  // Windows
			std::string nulFile = fileName + dirSep + "nul";
			std::ofstream os(nulFile.c_str(), std::ios::app);
			return (bool)os;
		#endif
	}
	
	inline void CreateDirectory(const std::string& dir)
	{
		if(DirectoryExists(dir))
			return;
		auto error = _mkdir(dir.c_str());
		if(error != 0)
			throw "Cannot create directory \"" + dir + "\", error code " + ToString(error);
	}
	  
#endif


namespace ph = physicalConstants;
 
// Strange way to write something in the binary mode
template<typename T>
void bin_write(std::ofstream& stream, T* buf, int count)
{
	stream.write(reinterpret_cast<char*>(buf), count*sizeof(T));
}

// Return a ctime time stamp
String timeStamp()
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
    String ctime_stamp = std::ctime(&now_time_t);
    return ctime_stamp.substr(0, ctime_stamp.length()-1); // Reduce an idiotic newline
}

// Make a directory where the binaries and JSONs will be put
// (and make an inner directory for the binaries)
String setDir(String root)
{
    String time = timeStamp();
    String dir = root + "NSSI NLM " + time + "/";
	
	#ifndef __DONT_USE_STD_FILESYSTEM__
		fs::create_directories(dir + "bin/"); //  This creates both dir and bin inside the dir
		fs::create_directory(dir + "plots/");
	#else
		CreateDirectory(root);
		CreateDirectory(dir);
		CreateDirectory(dir + "bin/");
		CreateDirectory(dir + "plots/");
	#endif
	
    return dir;
}

// Dump the x-grid
void dumpXGrid(const Constants& c, int periodN_x, String dir)
{
    // Note! Here I aim to also ssave the periodical 
    // continuation of the numerical solution, so I need
    // the (N_x+1)th point! This point is same with the first point
    // in the solution, that this is better for displaying
    // (similar to addFirst in NSSI Py)
    std::ofstream out(dir + "XGrid.txt");
    for (int i=0; i<c.N_x+1; ++i)
    {
        if (i % periodN_x == 0)
        {
            out << i*c.dx/ph::km << " ";
        }
    }
    out.close();
}

// Dump TwoLines and add position in z axis
void dumpTwoLines(const Constants& c, Real z, int periodN_x, const TwoLines& tl, String dir)
{
    // Add z position
    std::ofstream outZ(dir + "ZGrid.txt", std::ios::app);
    outZ << z / ph::km << " ";
    outZ.close();
    // Add probabilities
    std::ofstream outLeft(dir  + "bin/left.bin",  std::ios::app | std::ios::binary);
    std::ofstream outRight(dir + "bin/right.bin", std::ios::app | std::ios::binary); // std::ios::app | std::ios::binary
    for (int i=0;i<c.N_x+1; ++i)
    {
        if (i % periodN_x == 0)
        {
            // left beam probabilities
            Real ProbsL[4] = {tl.left(i)(0,0).real(),
                              tl.left(i)(1,1).real(),
                              tl.left(i)(2,2).real(),
                              tl.left(i)(3,3).real()};
            // right beam probabilities
            Real ProbsR[4] = {tl.right(i)(0,0).real(),
                              tl.right(i)(1,1).real(),
                              tl.right(i)(2,2).real(),
                              tl.right(i)(3,3).real()};
            
            bin_write(outLeft,  ProbsL, 4);
            bin_write(outRight, ProbsR, 4);
        }
    }
    outLeft.close();
    outRight.close();
}

// Writes the message in the stdout
void stdoutLog(String message)
{
    std::cout << "[" << timeStamp() << "] " << message << std::endl;
}