#include <iostream>
#include <string>

using namespace std;

int main()
{
	string tmp;
	cin >> tmp;
	while(cin.good())
	{
		cout << tmp.substr(2, tmp.size()-2) << endl;
		getline(cin, tmp);
		getline(cin, tmp);
		getline(cin, tmp);
		getline(cin, tmp);
		getline(cin, tmp);
		getline(cin, tmp);
		while(true)
		{
			cin >> tmp;
			if (tmp.find('.') != string::npos 
				|| tmp.find('(') != string::npos 
				|| tmp.find(')') != string::npos)
			{
				break;
			}
			cout << tmp;
		}
		cout << endl;
		cout << tmp;
		while(cin.good())
		{
			cin >> tmp;
			if(tmp[0] == '%')
				break;
			cout << tmp;
		}
		cout << endl;

	}

}