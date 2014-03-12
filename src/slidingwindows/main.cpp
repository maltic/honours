//Compile from project root
//g++ src/slidingwindows/main.cpp src/common/vienna/libRNA.a -std=c++11 -O2
#include <iostream>
#include <string>
#include "prediction.h"
#include "magic_seq.h"
#include "precomputed_windows.h"
#include "../common/rna_util.h"
#include <chrono>

void test_ab_splat (int a, float b)
{


	// this is extremely suboptimal, as we could just read in the scraped RNA data
	// however this saves me having to write more code, and wirte a new way to sort stuff
	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);
	std::sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);
	std::cout << "Finished loading precomputed windows!" << std::endl;

	std::cout << "Starting AB-splat test with a = " << a << " and b = " << b << std::endl;

	for (auto & pc : precomp)
	{
		std::string ss = ab_splat (a, b, pc);
		float f1score = calc_f1score (pc.actual_sstruct, ss);
		std::cout << f1score << "\t" << pc.rna.size() << std::endl;
	}

	std::cout << "Done!" << std::endl;
}


void brute_force_ab()
{
	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);

	// shuffle randomly so odd and even indexes are random sets of even* size
  	// * as even as possible!
  	std::random_shuffle ( precomp.begin(), precomp.end() );

	std::cout << "Finished loading precomputed windows!" << std::endl;


  	

  	std::cout << "Starting brute force..." << std::endl;

	for (int a = 10; a <= 30; ++a)
	{
		for (float b = 1.5; b <= 4.0; b += 0.1)
		{
			float sumf1_even = 0.0, sumf1_odd = 0.0;
			int ne = 0, no = 0;
			for (int i = 0; i < precomp.size(); ++i)
			{
				std::string ss = ab_splat (a, b, precomp[i]);
				float f1score = calc_f1score (precomp[i].actual_sstruct, ss);
				if (i % 2 == 0)
				{
					ne++;
					sumf1_even += f1score;
				}
				else
				{
					no++;
					sumf1_odd += f1score;
				}
				//std::cout << f1score << "\t" << i.rna.size() << std::endl;
			}
			std::cout << a << "\t" << b << "\t" << (sumf1_even / ne) << "\t" << (sumf1_odd / no) << std::endl;
		}
	}


	std::cout << "Done!" << std::endl;


}


// Runs Zukers algorithm (RNAfold) on the same testing data as used by run_selection_tests
// Used for comparison
void run_zuker_test()
{
	// this is extremely suboptimal, as we could just read in the scraped RNA data
	// however this saves me having to write more code, and wirte a new way to sort stuff
	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);
	std::sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);
	std::cout << "Finished loading precomputed windows!" << std::endl;


	std::cout << "Testing RNAfold..." << std::endl;
	for (auto & pc : precomp)
	{
		RNAInterval folded = zuker_fold(pc.rna, 0, pc.rna.size()-1);
		float f1score = calc_f1score (pc.actual_sstruct, folded.sstruct);
		std::cout << f1score << "\t" << pc.rna.size() << std::endl;
	}
	std::cout << "Done!" << std::endl;


}

// This function is used to test the various different selection strategies
// Note to self, link to thesis section
void run_selection_tests (const int n_windows)
{
	SelectionTester st;
	std::cout << "Testing selection algorithns..." << std::endl;

	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);

	std::sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);

	std::cout << "Finished loading precomputed windows!" << std::endl;

	std::cout << "Weighted Activity Selection: " << std::endl;
	st.selection_algorithm_test (n_windows, weighted_activity_selection, precomp);
	std::cout << "Top Down Selection: " << std::endl;
	st.selection_algorithm_test (n_windows, top_down_selection, precomp);
	std::cout << "Bottom Up Selection: " << std::endl;
	st.selection_algorithm_test (n_windows, bottom_up_selection, precomp);
	std::cout << "MFE Selection: " << std::endl;
	st.selection_algorithm_test (n_windows, greedy_MFE_selection, precomp);

	std::cout << "Done!" << std::endl;
}


void run_magic_seq_training()
{
	std::cout << "Training a magic sequence!" << std::endl;
	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin, 300);
	std::cout << "Finished loading precomputed windows!" << std::endl;

	MagicSequenceOptimizer mso (2048, 256);

	mso.optimize (precomp);
	std::cout << "Done!" << std::endl;
}

int main()
{

	test_ab_splat (24, 1.8);
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
