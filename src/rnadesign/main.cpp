
//g++ rnadesign/main.cpp common/vienna/libVienna.a -std=c++11 -O2

#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <vector>
#include <stack>
#include <algorithm>
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

class GARNADesigner
{
private:


	class Genome
	{
	private:
		GARNADesigner *parent = NULL;
		std::vector<int> errors;
	public:
		int fitness = 0;
		std::string rna;

		bool operator< (const Genome& other) const
		{
			return fitness < other.fitness;
		}

		Genome (GARNADesigner *p, const std::string& rna)
		{
			this->parent = p;
			this->rna = rna;
		}

		Genome (GARNADesigner *p)
		{
			this->parent = p;
			this->rna = std::string ( (*parent).target.size(), 'X' );

			for (int i = 0; i < (*parent).target.size(); ++i)
			{

				// in already defined, skip
				if (this->rna[i] != 'X')
					continue;

				// choose a random base
				char base = bases [(*parent).prng() % 4];
				this->rna[i] = base;

				// select a random but valid bonding base, if this base has a bond
				// as defined by the target_structure
				if ( (*parent).matching[i] != -1 )
				{
					std::vector<char> valid;
					for (int j = 0; j < 4; ++j)
						if ( valid_bond (base, bases[j]) )
							valid.push_back (bases[j]);

					this->rna[(*parent).matching[i]] = valid[ (*parent).prng() % valid.size() ];
				}
			}

			calc_fitness();

		}

		void calc_fitness()
		{
			// eventually this will change to a boltzmann based fold
			// or better yet a fully custom CFG approach
			std::string sstruct = zuker_fold ( rna, 0, rna.size() - 1 ).sstruct;
			fitness = 0;
			for (int i = 0; i < sstruct.size(); ++i)
				if (sstruct[i] != (*parent).target[i])
				{
					fitness += 1;
					this->errors.push_back (i);
				}
		}

		Genome mutate()
		{
			Genome next (parent, rna);

			for (int & error : errors)
				next.rna[error] = 'X';

			// replace the errors with something different
			for (int & error : errors)
			{
				// in already defined, skip
				if (next.rna[error] != 'X')
					continue;

				// choose a random base
				char base = bases [parent->prng() % 4];
				next.rna[error] = (base);

				// select a random but valid bonding base, if this base has a bond
				// as defined by the target_structure
				if (parent->matching[error] != -1)
				{
					std::vector<char> valid;
					for ( int j = 0; j < 4; ++j )
						if ( valid_bond (base, bases[j]) )
							valid.push_back (bases[j]);

					next.rna[parent->matching[error]] = valid[ parent->prng() % valid.size() ];
				}
			}

			next.calc_fitness();

			return next;
		}

		Genome crossover (const Genome& other)
		{
			Genome next (parent, rna);

			// replace any bases that caused an error with the equivalent bases in other
			for (int & error : errors)
				next.rna[error] = other.rna[error];

			next.calc_fitness();

			return next;

		}

	};



	void design_rna_step ()
	{
		int keepers = this->population_size * this->trunc_pcnt;

		std::vector<Genome> next_pop;

		int i = 0;

		while (next_pop.size() < this->population_size)
		{
			int odds = prng() % (elitism_odds + mutation_odds + crossover_odds + new_odds);

			if (odds < elitism_odds)
				next_pop.push_back (this->population[i]);

			else if (odds < elitism_odds + mutation_odds)
				next_pop.push_back ( this->population[i].mutate() );

			else if (odds < elitism_odds + mutation_odds + crossover_odds)
				next_pop.push_back ( this->population[i].crossover (population[prng() % keepers]) );

			else
				next_pop.push_back ( Genome (this) );

			i = (i + 1) % keepers;
		}

		population = next_pop;

		std::sort ( this->population.begin(), this->population.end() );
	}

	unsigned seed;
	std::minstd_rand0 prng;
	std::vector<Genome> population;
	std::vector<int> matching;
	std::string target;


public:

	// Number of evolutionary steps taken in design_rna,
	// basically the number of calls to design_rna_step
	int evolutionary_iterations = 1000;

	int population_size = 1000;

	float trunc_pcnt = 0.4;

	// odds for various breeding events
	int elitism_odds = 1;
	int mutation_odds = 10;
	int crossover_odds = 3;
	int new_odds = 1;

	GARNADesigner()
	{
		this->seed = std::chrono::system_clock::now().time_since_epoch().count();
		prng = std::minstd_rand0 (seed);
	}


	// Takes a dot bracket version of the target structure
	// Returns a string version of RNA the algorithm deduced
	// This RNA is the most likely to fold into the target_structure which could be found
	std::string design_rna (const std::string& target_structure)
	{
		target = target_structure;
		matching = get_matching_bonds(target_structure);

		for (int i = 0; i < population_size; ++i)
			population.push_back ( Genome (this) );

		std::sort ( population.begin(), population.end() );

		for (int i = 0; i < evolutionary_iterations; ++i)
		{
			std::cout << "RNA Design step " << i << ". Best fitness = " << population[0].fitness << std::endl;
			design_rna_step();
		}

		return population[0].rna;
	}


};



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
	GARNADesigner designer;
	std::cout << designer.design_rna (target_structure) << std::endl;
	//std::cout << design_rna (10000, target_structure) << std::endl;
	return 0;
}