#include <Fl/Fl.H>
#include <string>
#include <sstream>
#include <iostream>
#include "Common/Common.h"

using namespace std;

double Str::parseDouble(const string& str) {
	double ret;
	stringstream sst(str);
	sst>>ret;
	return ret;
}

int Str::parseInt(const string& str) {
	int ret;
	stringstream sst(str);
	sst>>ret;
	return ret;
}

bool Str::isLetter(char c) {
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

bool Str::isNumber(char c) {
	return (c >= '0' && c <= '9');
}

vector<string> Str::split(const string& s, char c) {
	vector<string> ret;
	string cur = "";
	for(int j = 0;j<s.length();j++) {
		if(s[j] == c) {
			ret.push_back(cur);
			cur = "";
		}
		else {
			cur += s[j];
		}
	}
	if(cur.length() > 0)
		ret.push_back(cur);

	vector<string> ret2;
	for(int j = 0; j < ret.size(); j++)
		if(ret[j].length() > 0)
			ret2.push_back(ret[j]);

	return move(ret2);
}
