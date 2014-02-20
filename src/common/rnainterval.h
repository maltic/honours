

#ifndef RNAINTERVAL_H
#define RNAINTERVAL_H

#include <string>

/*
	The RNA Inteval class
	Used to represent an interval of an RNA Secondary structure.
	For use with algorithms that require sub-structures, such as an RNA sliding window algorithm.
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
		return right < other.right;
	}
};

bool size_comp(const RNAInterval& a, const RNAInterval& b)
{
	return a.size() < b.size();
}

bool fe_comp(const RNAInterval& a, const RNAInterval& b)
{
	return a.score < b.score;
}

#endif