
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



std::string design_rna (const std::string& target_structure)
{


	std::vector<int> matching = get_matching_bonds (target_structure);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 prng (seed);

	// init the proband
	std::string proband (target_structure.size(), 'X');

	for ( int i = 0; i < matching.size(); ++i )
	{

		if (proband[i] != 'X')
			continue;

		char base = bases [prng() % 4];
		proband[i] = (base);

		if (matching[i] != -1)
		{
			std::vector<char> valid;
			for ( int j = 0; j < 4; ++j )
				if ( valid_bond (base, bases[j]) )
					valid.push_back (bases[j]);

			proband[ matching[i] ] = valid[ prng() % valid.size() ];
		}
	}


	// iteratively improve the proband

	return proband;
}


int main()
{
	std::string target_structure;
	std::cin >> target_structure;
	std::cout << design_rna (target_structure) << std::endl;
	return 0;
}