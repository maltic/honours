//Compile from project root
//g++ SlidingWindows/main.cpp vienna/libVienna.a
#include <iostream>
#include <string>
#include "prediction.h"
#include <chrono>

int test_splat()
{
	std::string rna, name;
	std::string targetStructure;
	std::vector<std::string> rnas;
	std::vector<std::string> targets;
	while(std::cin.good())
	{
		std::cin >> name >> rna >> targetStructure;
		rnas.push_back(rna);
		targets.push_back(targetStructure);
		std::cin >> std::ws;
	}




	for(int sz = 20; sz <= 50; ++sz)
	{
		int total_errors = 0;
		for(int i = 0; i < rnas.size(); ++i)
		{
			total_errors += splat_prediction(sz, rnas[i], targets[i]);
		}
		std::cout << sz << " " << total_errors << std::endl;
	}

}

int main()
{
	std::string rna, name;
	std::string targetStructure;

	int tot = 0;
	while(std::cin.good())
	{
		
		std::cin >> name >> rna >> targetStructure;
		std::cout << "SSTRAND ID = " << name << std::endl;


		// typedef std::chrono::high_resolution_clock Clock;
	 //    typedef std::chrono::milliseconds milliseconds;
	 //    Clock::time_point t0 = Clock::now();
	    RNAInterval vanilla = zuker_fold(rna, 0, rna.size()-1);
	    // Clock::time_point t1 = Clock::now();
	    // milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	    // std::cout << ms.count() << "ms\n";

		//normal RNAfold
		int vanillaErrors = count_errors(targetStructure, vanilla.sstruct);
		tot += vanillaErrors;
		//std::cout << "RNA size: " << rna.size() << " with " << vanillaErrors << " errors." << std::endl;

		// splat_prediction(rna, targetStructure);


		
		std::cin >> std::ws;

	}

	std::cout << tot << std::endl;
	
	return 0;
}
