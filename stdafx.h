// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <random>
#include <stdexcept>

#define WXUSINGDLL 1

//#include <wx/image.h>
//#include <wx/quantize.h>

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

#include <bcrypt.h>

#define WXUNUSED(identifier) /* identifier */
typedef int16_t wxInt16;
typedef uint16_t wxUint16;
typedef int32_t wxInt32;

typedef unsigned char JSAMPLE;
typedef JSAMPLE *JSAMPROW;

void DoQuantize(unsigned w, unsigned h, unsigned char **in_rows, unsigned char **out_rows,
	unsigned char *palette, int desiredNoColours);
