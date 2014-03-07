

#ifndef RNAINTERVAL_H
#define RNAINTERVAL_H

#include <string>
#include <sstream>

/*
	The RNA Inteval class
	Used to represent an interval of an RNA Secondary structure.
	For use with algorithms that require sub-structures, such as an RNA sliding window algorithm.

	The term RNA Interval is interchanable with 'window' in the thesis
*/
struct RNAInterval
{
	//Left and right indexes in the entire RNA molecule
	unsigned int left, right;
	//The score (usually free energy, but doesn't have to be)
	double score;
	//The dot bracket notation
	std::string sstruct;
	RNAInterval(unsigned int l, unsigned int r, double s, const std::string& ss)
	{
		this->left = l;
		this->right = r;
		this->score = s;
		this->sstruct = ss;
	}
	RNAInterval(unsigned int l, unsigned int r, double s)
	{
		this->left = l;
		this->right = r;
		this->score = s;
		this->sstruct = "";
	}

	RNAInterval (const std::string& ri)
	{
		std::stringstream ss (ri);
		ss >> this->sstruct >> this->left >> this->right >> this->score;
	}


	RNAInterval (std::istream& strm)
	{
		strm >> this->sstruct >> this->left >> this->right >> this->score;
	}

	int size() const
	{
		return this->right - this->left;
	}
	//Check if two RNA Intevals do not overlap
	bool compatible_with(const RNAInterval& other) const
	{
		return other.right < this->left || other.left > this->right;
	}
	bool overlaps(const RNAInterval& other) const
	{
		return !this->compatible_with(other);
	}
	bool operator< (const RNAInterval& other) const
	{
		if (right == other.right)
			return left < other.left;
		return right < other.right;
	}
	std::string to_string() const
	{
		std::stringstream ss;
		ss << sstruct << " " << left << " " << right << " " << score;
		return ss.str();
	}
};

bool rna_int_size_comp(const RNAInterval& a, const RNAInterval& b)
{
	if ( a.size() == b.size() )
		return a < b;
	return a.size() < b.size();
}

bool rna_int_fe_comp(const RNAInterval& a, const RNAInterval& b)
{
	if (a.score == b.score)
		return a < b;
	return a.score < b.score;
}

#endif