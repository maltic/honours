//Compile from project root
//g++ src/slidingwindows/main.cpp src/common/vienna/libRNA.a -std=c++11 -O2
#include <iostream>
#include <string>
#include "prediction.h"
#include "magic_seq.h"
#include "precomputed_windows.h"
#include "../common/rna_util.h"
#include <chrono>


void magic_seq_train()
{
	std::string rna, name;
	std::string targetStructure;
	std::vector<std::string> rnas;
	std::vector<std::string> targets;
	while(std::cin.good())
	{
		std::cin >> name >> rna >> targetStructure;
		std::cin >> std::ws;
		// if (rna.size() < 300)
		// 	continue;
		rnas.push_back(rna);
		targets.push_back(targetStructure);
	}
	MagicSequenceOptimizer mso (32, 100);

	mso.optimize (rnas, targets);
}

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


	for (float mult = 2.1; mult <= 2.1; mult += 0.1)
	{
		for (float cutoff = 3; cutoff <= 5; cutoff += 1.0)
		{
			for (int start = 10; start < 30; start += 4)
			{
				int total_errors = 0;
				for (int i = 0; i < rnas.size(); ++i)
					total_errors += splat_tester(mult, cutoff, start, rnas[i], targets[i]);

				std::cout << mult << " " << cutoff << " " << start << " " << total_errors << std::endl;

			}
		}
	}

}

// This function is used to test the various different selection strategies
// Note to self, link to thesis section
void run_selection_tests ()
{
	std::cout << "Testing selection algorithns..." << std::endl;

	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);

	std::sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);

	std::cout << "Finished loading precomputed windows!" << std::endl;

	std::cout << "Weighted Activity Selection: " << std::endl;
	selection_algorithm_test (weighted_activity_selection, precomp);
	std::cout << "Top Down Selection: " << std::endl;
	selection_algorithm_test (top_down_selection, precomp);
	std::cout << "Bottom Up Selection: " << std::endl;
	selection_algorithm_test (bottom_up_selection, precomp);
	std::cout << "MFE Selection: " << std::endl;
	selection_algorithm_test (greedy_MFE_selection, precomp);

	std::cout << "Done!" << std::endl;
}

int main()
{
	run_selection_tests();
	return 0;


	//all bullshit bellow this point
	magic_seq_train();
	return 0;

	std::string rna, name;
	std::string targetStructure;

	//for splat prediction
	int magic[] = {49, 95, 233};
	std::vector<int> magicVec(magic, magic+3);

	int tot = 0, zukTot = 0;
	while(std::cin.good())
	{
		
		std::cin >> name >> rna >> targetStructure;

		if (rna.size() < 300)
			continue;

		std::cout << "SSTRAND ID = " << name << std::endl;

		typedef std::chrono::high_resolution_clock Clock;
	    typedef std::chrono::milliseconds milliseconds;
	    Clock::time_point t0 = Clock::now();
		RNAInterval vanilla = zuker_fold(rna, 0, rna.size()-1);
		Clock::time_point t1 = Clock::now();
		milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	    std::cout << "RNAfold took " << ms.count() << "ms\n";

	    std::cout << "RNAfold Sensitivity: " << calc_sensitivity (targetStructure, vanilla.sstruct) 
	    	<< " PPV: " << calc_ppv(targetStructure, vanilla.sstruct) << std::endl;

		int vanillaErrors = count_errors(targetStructure, vanilla.sstruct);
		std::cout << "RNA size: " << rna.size() << " with " << vanillaErrors << " Zuker errors." << std::endl;

		int err = splat_prediction(magicVec, rna, targetStructure);

		tot += err;
		zukTot += vanillaErrors;


		// typedef std::chrono::high_resolution_clock Clock;
	 //    typedef std::chrono::milliseconds milliseconds;
	 //    Clock::time_point t0 = Clock::now();
	    //RNAInterval vanilla = zuker_fold(rna, 0, rna.size()-1);
	    // Clock::time_point t1 = Clock::now();
	    // milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	    // std::cout << ms.count() << "ms\n";

		//normal RNAfold
		//int vanillaErrors = count_errors(targetStructure, vanilla.sstruct);
		//tot += vanillaErrors;
		//std::cout << "RNA size: " << rna.size() << " with " << vanillaErrors << " errors." << std::endl;

		// int fib[] = {21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181};
		// std::vector<int> fibVec(fib, fib+12);

		

		// tot += splat_prediction(magicVec, rna, targetStructure);


		
		std::cin >> std::ws;

	}

	std::cout << "Total splat errors = " << tot << ", total Zuker errors = " << zukTot << std::endl;
	
	return 0;
}
