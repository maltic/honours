#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

vector<string> split(const string& str)
{
	stringstream ss(str);
	string item;
	vector<string> elems;
	while(getline(ss,item, ' '))
		elems.push_back(item);
	return elems;
}

int main()
{
	string in;
	vector<string> line;
	int inint;
	string delim = "\t";
	while(cin.good())
	{
		getline(cin, in);
		line = split(in);
		cout << line[3] << delim;
		getline(cin, in);
		line = split(in);
		cout << line[2] << delim << line[4] << delim;
		getline(cin, in);
		getline(cin, in);
		line = split(in);
		if (line[2] == "NOT")
		{
			cout << "0" << delim << "0";
		}
		else
		{
			cout << line[5] << delim << line[7];
		}
		cout << endl;
		getline(cin, in);
		getline(cin, in);
		getline(cin, in);
		cin >> ws;
	}
}