#include <iostream>
#include <cctype>

using namespace std;

int main()
{
	string in;
	while(cin >> in)
	{
		for(int i = 0; i < in.size(); ++i)
			in[i] = toupper(in[i]);
		if(in == ",")
			cout << endl;
		else
			cout << in;
	}
	return 0;
}
