#include <iostream>
#include <string>

using namespace std;


int main()
{
	string line;
	while(cin.good())
	{
		getline(cin, line);
		if(line[0] == '#')
			continue;
		
	}
	cout << "The end" << endl;
}