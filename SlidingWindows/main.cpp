//Compile from project root
//g++ SlidingWindows/main.cpp vienna/libVienna.a
#include <iostream>
#include <string>
#include "prediction.h"

int main()
{
	// testSelection();
	std::string rna, name;
	std::string targetStructure;
	while(std::cin.good())
	{
		
		std::cin >> name >> rna >> targetStructure;
		std::cout << "SSTRAND ID = " << name << std::endl;
		//normal RNAfold
		RNAInterval vanilla = zuker_fold(rna, 0, rna.size()-1);
		int vanillaErrors = count_errors(targetStructure, vanilla.sstruct);
		std::cout << "RNA size: " << rna.size() << " with " << vanillaErrors << " errors." << std::endl;

		two_window_prediction(rna, targetStructure);


		
		std::cin >> std::ws;

	}
	
	return 0;
}
