#ifndef MAGICSEQ_H
#define MAGICSEQ_H

// This file contains all the classes related to 'Magic Sequence Prediction'.
// A magic sequence is a simply a sequence of window sizes used for generating good sliding windows.
// I have called it 'magic' because I am either guessing good windows, or using a GA to generate them.
// There is no underlying methodology used, or deep understanding of RNA behind these numbers.

#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <iostream>


// The follow two classes (MagicSequence, and MagicSequenceOptimizer) implement a Generic Algorithm (GA).
// This GA is used to find good magical sequence to use during sliding window prediction.

struct MagicSequence
{
	float fitness;
	std::vector<int> sequence;

	MagicSequence(const std::vector<int>& seq)
	{
		this->set_seq (seq);
		this->fitness = 0.0;
	}

	void set_seq(const std::vector<int>& seq)
	{
		this->sequence = seq;
		std::sort (this->sequence.begin(), this->sequence.end());
		std::vector<int>::iterator it = std::unique ( this->sequence.begin(), this->sequence.end() );
		this->sequence.resize ( std::distance (this->sequence.begin(), it) );
	}

	void calc_fitness( const std::vector<std::string>& rnas, const std::vector<std::string>& targets )
	{
		this->fitness = 0.0;
	}

	bool operator< (const MagicSequence& other) const
	{
		// want to sort in ascending order
		return fitness > other.fitness;
	}

};


class MagicSequenceOptimizer
{

protected:
	// protected class members
	std::vector<MagicSequence> genomes;
	std::vector<std::string> rnas;
	std::vector<std::string> targets;
	std::minstd_rand0 generator;

	int valid_random_num()
	{
		return this->generator() % (max_num - min_num) + min_num;
	}

	MagicSequence gen_magic_seq()
	{
		// get a sequence random but valid numbers
		int num = this->generator() % max_nums + 1;
		std::vector<int> seq (num);
		for (int j = 0; j < num; ++j)
			seq.push_back (this->valid_random_num());
		// turn them into a magic sequence as if by magic
		MagicSequence ms (seq);
		ms.calc_fitness (this->rnas, this->targets);
		return ms;
	}

	MagicSequence breed (const MagicSequence& a, const MagicSequence& b)
	{
		// breed two sequences together
		
		// pick a pivot point on sequence a
		int pivot = this->generator() % a.sequence.size();

		std::vector<int> seq;

		// get itmes from a
		for (int i = 0; i <= pivot; ++i)
			seq.push_back (a.sequence[i]);

		// add all items from b that are valid
		for (int i = 0; i < b.sequence.size() && seq.size() < max_nums; ++i)
			if (b.sequence[i] > seq.back())
				seq.push_back (b.sequence[i]);

		MagicSequence ms (seq);
		ms.calc_fitness (rnas, targets);
		return ms;
	}

	MagicSequence mutate (const MagicSequence& ms)
	{
		// mutate a magic sequence

		std::vector<int> v = ms.sequence;

		// deletions
		int deletions = this->generator() % max_deletions + 1;
		if (deletions > ms.sequence.size())
			deletions = ms.sequence.size();
		
		std::shuffle (v.begin(), v.end(), generator); //so we can delete randomly
		for (int i = 0; i < deletions; ++i)
			v.pop_back();

		// insertions
		int insertions = this->generator() % max_insertions + 1;
		if (insertions + v.size() > max_nums)
			insertions = max_nums - v.size();
		for (int i = 0; i < insertions; ++i)
			v.push_back (this->valid_random_num());

		MagicSequence ns (v);
		ns.calc_fitness (this->rnas, this->targets);

		return ns;
	}



public:

	// min and max for a magic sequnce number
	int min_num = 10;
	int max_num = 2500;

	// maximum number of items in a magic sequence
	int max_nums = 20;

	// GA settings
	int num_genomes;
	int num_generations;
	float truncation_pcnt = 0.4;
	int elitism = 1;

	// odds for selection events
	int breed_odds = 4;
	int mutation_odds = 4;
	int replacement_odds = 1;

	// mutation operators
	int max_deletions = 2;
	int max_insertions = 2;


	MagicSequenceOptimizer (int n_genomes, int gens)
	{
		// init prng
  		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  		this->generator = std::minstd_rand0 (seed);

  		this->num_generations = gens;
  		this->num_genomes = n_genomes;
	}


	void evolve()
	{
		// select by sorting into ascending order
		std::sort (genomes.begin(), genomes.end());
		
		int survivors = (int) (truncation_pcnt * genomes.size());

		int total_odds = breed_odds + mutation_odds + replacement_odds;

		// generate the next generation

		std::vector<MagicSequence> next_gen;

		for (int i = 0; i < elitism; ++i)
			next_gen.push_back (genomes[i]);

		// the current genome we're doing a selection event on
		int curr_g = 0;

		std::cout << "about to select events " << survivors << " " << genomes.size() << std::endl;

		while (next_gen.size() < num_genomes)
		{
			int odd = this->generator() % total_odds;

			//breed
			if (odd < breed_odds)
			{
				next_gen.push_back ( breed (genomes[curr_g], genomes[this->generator() % survivors]) );
			}
			//mutate
			else if (odd < breed_odds + mutation_odds)
			{
				next_gen.push_back ( mutate (genomes[curr_g]) );
			}
			//replace
			else
			{
				next_gen.push_back (this->gen_magic_seq());
			}

			// go to next selected genome with wraparound
			curr_g = (curr_g + 1) % survivors;
		}

		genomes = next_gen;
	}

	std::vector<MagicSequence> optimize ( const std::vector<std::string>& rnas, const std::vector<std::string>& targets )
	{ 
		// fill the test data
		this->rnas = rnas;
		this->targets = targets;

		// init the genomes
		this->genomes = std::vector<MagicSequence> (genomes);
		for (int i = 0; i < num_genomes; ++i)
			genomes.push_back (this->gen_magic_seq());


		// evolve through generations
		for (int i = 0; i < num_generations; ++i)
			evolve();

		return this->genomes;
	}


};

#endif