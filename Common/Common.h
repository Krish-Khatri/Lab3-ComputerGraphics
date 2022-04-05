#ifndef COMMON_H
#define COMMON_H

#include <Fl/Fl.H>
#include <string>
#include <vector>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#ifndef FINF32
#define FINF32 1e9f
#endif

#ifndef DINF
#define DINF 1e9
#endif

using namespace std;

const Fl_Color WIN_COLOR = fl_rgb_color(244, 247, 251);

namespace Str {
	double parseDouble(const string& str);
	int parseInt(const string& str);
	bool isLetter(char c);
	bool isNumber(char c);
	vector<string> split(const string& s, char c);
}

#endif