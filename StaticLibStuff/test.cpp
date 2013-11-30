#include "ViennaLib.h"

int main()
{
	SStruct s = RNAfold("UAGCUACUGUACGUAUCGUAUCGAC");
	cout << s.dbNotation << endl;
	return 0;
}