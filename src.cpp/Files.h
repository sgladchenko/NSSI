#pragma once

#include "Types.h"
#include "Constants.h"
#include "Containers.h"

// Return a time stampt in the format of ctime
String timeStamp();

// Make a directory where the binaries and JSONs will be put
// (and make a inner directory for the binaries)
String setDir(String root);

// Dump the x-grid
void dumpXGrid(const Constants& c, int periodN_x, String dir);

// Dump TwoLines and add position in z axis
void dumpTwoLines(const Constants& c, Real z, int periodN_x, const TwoLines& tl, String dir);

// Stdout logger
void stdoutLog(String message);