
//g++ rnadesign/main.cpp common/vienna/libVienna.a -std=c++11 -O2

#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <vector>
#include <stack>
#include "../common/vienna.h"

char bases [4] = { 'A', 'U', 'G', 'C' };


bool valid_bond (const char a, const char b) 
{
    return  (a == 'A' && b == 'U')
         || (a == 'U' && b == 'A')
         || (a == 'G' && b == 'C')
         || (a == 'C' && b == 'G')
         || (a == 'U' && b == 'G')
         || (a == 'G' && b == 'U');
}


std::vector<int> get_matching_bonds (const std::string& structure)
{
	std::stack<int> s;
	std::vector<int> matches (structure.size(), -1);

	for ( int i = 0; i < structure.size(); ++i )
	{
		if (structure[i] == ')')
		{

			if ( s.empty() )
			{
				std::cerr << "Empty stack, therefore invalid secondary structure" << std::endl;
				throw 1;
			}

			matches[i] = s.top();
			s.pop();
		}
		else if (structure[i] == '(')
			s.push (i);

	}

	return matches;
}



std::string design_rna (const int max_iterations, const std::string& target_structure)
{

	// init the random number generator
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 prng (seed);

	//get all the bond pairs in the target structure
	std::vector<int> matching = get_matching_bonds(target_structure);


	// init the proband
	std::string proband (target_structure.size(), 'X');

	for ( int i = 0; i < matching.size(); ++i )
	{

		// in already defined, skip
		if (proband[i] != 'X')
			continue;

		// choose a random base
		char base = bases [prng() % 4];
		proband[i] = (base);

		// select a random but valid bonding base, if this base has a bond
		// as defined by the target_structure
		if (matching[i] != -1)
		{
			std::vector<char> valid;
			for ( int j = 0; j < 4; ++j )
				if ( valid_bond (base, bases[j]) )
					valid.push_back (bases[j]);

			proband[matching[i]] = valid[ prng() % valid.size() ];
		}
	}

	// iteratively improve the proband
	for (int i = 0; i < max_iterations; ++i)
	{
		// fold using RNAfold
		std::string folded = zuker_fold (proband, 0, proband.size() - 1).sstruct;

		// find all the errors
		std::vector<int> errors;

		for (int j = 0; j < folded.size(); ++j)
			if (folded[j] != target_structure[j])
			{
				errors.push_back (j);
				proband[j] = 'X';
			}

		// replace the errors with something different
		for (auto &error : errors)
		{
			// in already defined, skip
			if (proband[error] != 'X')
				continue;

			// choose a random base
			char base = bases [prng() % 4];
			proband[error] = (base);

			// select a random but valid bonding base, if this base has a bond
			// as defined by the target_structure
			if (matching[error] != -1)
			{
				std::vector<char> valid;
				for ( int j = 0; j < 4; ++j )
					if ( valid_bond (base, bases[j]) )
						valid.push_back (bases[j]);

				proband[matching[error]] = valid[ prng() % valid.size() ];
			}
		}

		std::cout << "Iteration #" << (i + 1) << ": Total Errors = " 
			<< errors.size() << " Structure = '" << folded << "'" << std::endl;

	}

	return proband;
}


int main()
{
	std::string target_structure;
	std::cin >> target_structure;
	std::cout << design_rna (10000, target_structure) << std::endl;
	return 0;
}