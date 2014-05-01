#ifndef TESTING_H
#define TESTING_H

/*
	This file contains functions and classes for testing various algorithms
	related to the honours dissertation of Max Ward.
*/

#include <iostream>
#include <string>
#include <chrono>
#include <random>


#include "prediction.h"
#include "magic_seq.h"
#include "precomputed_windows.h"
#include "../common/rna_util.h"



// A SelectionTester object can be used to compare various selection algorithms
// Because I wanted to be able to test using an arbitrary number of best windows
// I needed to make a class to 'hide' the state variables for the recursion
class SelectionTester
{
private:

	// hidden recursion variables
	float best_f1score;
	std::vector<int> best_windows;
	std::vector<int> chosen;

	void selection_test_recursion (const int n, const t_selection_algorithm algo, const PrecomputedWindows& precomp)
	{
		if (n == 0)
		{
			std::vector<RNAInterval> windows;

			for (int i = 0; i < chosen.size(); ++i)
			{
				if (i == 0 || chosen[i] > chosen[i-1])
					windows.insert ( windows.end(), precomp.windows[chosen[i]].begin(), precomp.windows[chosen[i]].end() );
			}

			std::vector<int> selected = algo (windows);
			std::string sstruct = get_dotbracket (precomp.rna.size(), windows, selected);
			float f1score = calc_f1score (precomp.actual_sstruct, sstruct);

			if (f1score > best_f1score)
			{
				best_f1score = f1score;
				best_windows = chosen;
			}


		}
		else
		{
			// start the loop at the index of the previously selected window
			// this optimization is possible since index order doesnt matter
			// windows 1,2,3 are the same as 2,1,3
			int start = chosen.size() == 0 ? 0 : chosen.back();

			for (int i = start; i < precomp.windows.size(); ++i)
			{
				chosen.push_back (i);
				selection_test_recursion (n-1, algo, precomp);
				chosen.pop_back();
			}
		}
	}


public:

	// Given a n (the number of windows to use), a selection algorithm and a collection of precomputed windows
	// This function will perform a test on the accuracy of the selection function
	void test_selection_algorithm (const int n, const t_selection_algorithm algo, std::vector<PrecomputedWindows>& precomp)
	{
		float avg_f1 = 0.0;
		for (int i = 0; i < precomp.size(); ++i)
		{
			// init recursion variables
			best_windows.clear();
			best_f1score = -1.0;
			chosen.clear();
			//do recursion
			selection_test_recursion (n, algo, precomp[i]);

			// std::cout << precomp[i].name << ": Best F1 score = " << bestf1 << " at window size = " 
			// 	<< (best_window_index + MIN_WINDOW_SIZE) << " for RNA size = " << precomp[i].rna.size() << std::endl;

			// this commented line was used to produce spreadsheet compatible output
			std::cout << best_f1score << "\t";
			for (int & w : best_windows)
				std::cout << (w + MIN_WINDOW_SIZE) << " ";
			std::cout << "\t" << precomp[i].rna.size() << std::endl;

			avg_f1 += best_f1score;
		}
		avg_f1 = avg_f1 / precomp.size();
		std::cout << "Average best F1 score = " << avg_f1 << std::endl;
	}

	void test_accuracy_landscapes (const t_selection_algorithm algo, std::vector<PrecomputedWindows>& precomp) {
		for (int i = 0; i < precomp.size(); ++i)
		{
			std::cout << "Exploring accuracy landscape for " << precomp[i].name << std::endl;
			for (int sz = 0; sz < precomp[i].windows.size(); ++sz)
			{
				std::vector<int> selected = algo (precomp[i].windows[sz]);
				std::string sstruct = get_dotbracket (precomp[i].rna.size(), precomp[i].windows[sz], selected);
				float f1score = calc_f1score (precomp[i].actual_sstruct, sstruct);
				std::cout << "Window Size = " << sz + MIN_WINDOW_SIZE << ", f1score = " << f1score << std::endl;
			}
		}


	}
};

// Creates a new random number generator with the current system time as the seed
std::default_random_engine make_prng()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	return std::default_random_engine (seed);
}

// Runs tests for zukers algorithm (RNAfold)
// Loads precomputed windows from standard in, prints to standard out
void test_zukers()
{
	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;

	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);
	std::sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);
	std::cout << "Finished loading precomputed windows!" << std::endl;
	for (auto & pc : precomp)
	{

		Clock::time_point t0 = Clock::now();
		RNAInterval folded = zuker_fold(pc.rna, 0, pc.rna.size()-1);
		Clock::time_point t1 = Clock::now();
		milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
		float f1score = calc_f1score (pc.actual_sstruct, folded.sstruct);
		std::cout << f1score << "\t" << ms.count() << "\t" << pc.rna.size() << std::endl;
	}
	std::cout << "Done!" << std::endl;


}


// Runs tests for ab-splat algorithm
// Loads precomputed windows from standard in (if precompute=true), prints to standard out
void test_ab_splat (int a, float b, bool precompute = true)
{
	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;

	std::cout << "Loading precomputed windows..." << std::endl;

	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);
	std::sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);

	std::cout << "Finished loading precomputed windows!" << std::endl;

	std::cout << "Starting AB-splat test with a = " << a << " and b = " << b << std::endl;
		
	for (auto & pc : precomp)
	{

		Clock::time_point t0 = Clock::now();
		std::string ss = precompute ? ab_splat (a, b, pc) : ab_splat (a, b, pc.rna);
		Clock::time_point t1 = Clock::now();
		milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
		float f1score = calc_f1score (pc.actual_sstruct, ss);
		std::cout << f1score << "\t" << ms.count() << "\t" << pc.rna.size() << std::endl;
	}
	std::cout << "Done!" << std::endl;


}


// Brute force tests values for ab in ab splat
// Loads precomputed windows from standard in
void brute_force_ab()
{
	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);

	// shuffle randomly so odd and even indexes are random sets of even* size
  	// * as even as possible!
  	std::shuffle ( precomp.begin(), precomp.end(), make_prng() );

	std::cout << "Finished loading precomputed windows!" << std::endl;


  	std::cout << "Starting brute force..." << std::endl;


  	const int mina = 10, maxa = 30;
  	const float minb = 1.5, maxb = 4.0, bincrement = 0.1;

	for (int a = mina; a <= maxa; ++a)
	{
		for (float b = minb; b <= maxb; b += bincrement)
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


// This function is used to test the various different selection strategies
// Loads precomputed windows from standard in as per usual
void test_selection_algorithms (const int n_windows)
{
	SelectionTester st;
	std::cout << "Testing selection algorithns..." << std::endl;

	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);

	std::stable_sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);

	std::cout << "Finished loading precomputed windows!" << std::endl;

	std::cout << "Weighted Activity Selection: " << std::endl;
	st.test_selection_algorithm (n_windows, weighted_activity_selection, precomp);
	std::cout << "Top Down Selection: " << std::endl;
	st.test_selection_algorithm (n_windows, top_down_selection, precomp);
	std::cout << "Bottom Up Selection: " << std::endl;
	st.test_selection_algorithm (n_windows, bottom_up_selection, precomp);
	std::cout << "MFE Selection: " << std::endl;
	st.test_selection_algorithm (n_windows, greedy_MFE_selection, precomp);

	std::cout << "Done!" << std::endl;
}

// Takes precomputed windows from standard in
// Reports the f-scores for various window sizes
void test_accuracy_landscapes()
{
	SelectionTester st;
	std::cout << "Testing accuracy landscapes.." << std::endl;

	std::cout << "Loading precomputed windows..." << std::endl;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);

	std::stable_sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);

	std::cout << "Finished loading precomputed windows!" << std::endl;

	st.test_accuracy_landscapes (weighted_activity_selection, precomp);

	std::cout << "Done!" << std::endl;
}


// Takes precomputed windows from standard in, and uses them to compute
// a magic sequence using a GA
// for more info, see magic_seq.h
void train_magic_seq()
{
	std::cout << "Training a magic sequence!" << std::endl;
	std::cout << "Loading precomputed windows..." << std::endl;

	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);

	// random shuffle
	std::shuffle ( precomp.begin(), precomp.end(), make_prng() );

	std::cout << "Finished loading precomputed windows!" << std::endl;

	std::cout << "Running GA to find a 'magic sequence'..." << std::endl;

	MagicSequenceOptimizer mso (5000, 256);



	// first half is the training set
	std::vector<MagicSequence> pop = mso.optimize ( precomp.begin(), precomp.begin() + precomp.size() / 2 );

	std::vector<int> best = pop[0].sequence;

	std::cout << "Done! Best sequence was... " << std::endl;

	for (int & i : best)
		std::cout << i << " ";
	std::cout << std::endl;

	std::cout << "Using best found sequence on testing set..." << std::endl;


	// now compare to test set, which is second half
	for (std::vector<PrecomputedWindows>::iterator it = (precomp.begin() + precomp.size() / 2 + 1);
			it != precomp.end();
			++it)
	{
		std::string ss = splat_prediction_ga (best, *it);
		float f1score = calc_f1score (it->actual_sstruct, ss);
		std::cout << f1score << "\t" << it->rna.size() << std::endl;
	}

	std::cout << "Done!" << std::endl;
}


#endif