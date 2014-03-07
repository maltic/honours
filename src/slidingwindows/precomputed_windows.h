#ifndef PRECOMPUTED_WINDOWS_H
#define PRECOMPUTED_WINDOWS_H

#include <vector>
#include <string>
#include <iostream>
#include "../common/rnainterval.h"

const unsigned MIN_WINDOW_SIZE = 5;
const unsigned MAX_WINDOW_SIZE = 500;

// This header contains functions to precompute windows for RNA and print it to std out
// In addition, it allows one to load this precomputed RNA info from std in

// Represents all the precomputed windows for a given RNA record
struct PrecomputedWindows
{
	std::string name, rna, actual_sstruct;
	std::vector<std::vector<RNAInterval> > windows;
	PrecomputedWindows (std::istream& strm)
	{
		int num, sz;
		strm >> this->name >> this->rna >> this->actual_sstruct;
		std::string in;
		for (int i = MIN_WINDOW_SIZE; i <= MAX_WINDOW_SIZE; ++i)
		{
			strm >> sz >> num;
			if (sz != i)
				std::cerr << "Something went wrong, window size was wrong ( i = " 
					<< i << " vs sz = " << sz << ")" << std::endl;
			this->windows.push_back ( std::vector<RNAInterval>() );
			for (int j = 0; j < num; ++j)
				this->windows.back().push_back ( RNAInterval (strm) );
			if (this->windows.back().size() != num)
				std::cerr << "Something went wrong, number of windows was wrong" << std::endl;

		}

	}
};

bool precomputed_windows_size_cmp(const PrecomputedWindows& a, const PrecomputedWindows& b)
{
	return a.rna.size() < b.rna.size();
}

// Load precomputed windows from an input stream
std::vector<PrecomputedWindows> load_precomputed_windows (std::istream& strm)
{
	std::vector<PrecomputedWindows> precomputed;
	while ( std::cin.good() )
	{
		precomputed.push_back ( PrecomputedWindows (strm) );
		std::cin >> std::ws;
	}
	return precomputed;
}

// This function computes all possible sliding windows from the input stream
void precompute_windows()
{
	std::string rna, name, targetStructure;
	while(std::cin.good())
	{
		std::cin >> name >> rna >> targetStructure >> std::ws;

		std::cout << name << std::endl << rna << std::endl << targetStructure << std::endl;

		for (int sz = MIN_WINDOW_SIZE; sz <= MAX_WINDOW_SIZE; ++sz)
		{
			std::vector<RNAInterval> windows = rnal_fold(rna, sz);
			std::cout << sz << " " << windows.size() << std::endl;
			for (int i = 0; i < windows.size(); ++i)
				std::cout << windows[i].to_string() << std::endl;
			std::cerr << "Done sz: " << sz << std::endl;
		}

		std::cerr << "Done rna: " << name << std::endl;
	}
}

#endif