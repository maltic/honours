#ifndef RNA_UTIL_H
#define RNA_UTIL_H

#include <vector>
#include <stack>
#include <utility>
#include "rnainterval.h"


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

// Calculates all the bonding pairs from an RNA structure
// For example: '((..))' -> {5, 4, -1, -1, 1, 0} 
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
			matches[s.top()] = i;
			s.pop();
		}
		else if (structure[i] == '(')
			s.push (i);

	}

	return matches;
}


int count_true_positives (const std::vector<int>& matching_true, const std::vector<int>& matching_proband)
{
	unsigned tps = 0;
	std::vector<bool> marked (matching_true.size(), false);

	for (int i = 0; i < matching_true.size(); ++i)
	{
		if (marked[i])
			continue;
		if (matching_true[i] != -1 && matching_true[i] == matching_proband[i])
		{
			tps += 1;
			marked[matching_true[i]] = true;
		}
	}
	return tps;
}

int count_false_negatives (const std::vector<int>& matching_true, const std::vector<int>& matching_proband)
{
	unsigned fns = 0;
	std::vector<bool> marked (matching_true.size(), false);

	for (int i = 0; i < matching_true.size(); ++i)
	{
		if (marked[i])
			continue;
		if (matching_true[i] != -1 && matching_proband[i] == -1)
		{
			fns += 1;
			marked[matching_true[i]] = true;
		}
	}
	return fns;
}

int count_false_positives (const std::vector<int>& matching_true, const std::vector<int>& matching_proband)
{
	unsigned fps = 0;
	std::vector<bool> marked (matching_true.size(), false);

	for (int i = 0; i < matching_true.size(); ++i)
	{
		if (marked[i])
			continue;
		if (matching_proband[i] != -1 && matching_true[i] != matching_proband[i])
		{
			fps += 1;
			marked[matching_proband[i]] = true;
		}
	}
	return fps;
}

float calc_sensitivity (const std::string& target, const std::string& prediction)
{
	std::vector<int> mtrue = get_matching_bonds (target);
	std::vector<int> mpred = get_matching_bonds (prediction); 
	float tp = count_true_positives (mtrue, mpred);
	float fn = count_false_negatives (mtrue, mpred);
	return tp / (tp + fn);
}

float calc_ppv (const std::string& target, const std::string& prediction)
{
	std::vector<int> mtrue = get_matching_bonds (target);
	std::vector<int> mpred = get_matching_bonds (prediction); 
	float tp = count_true_positives (mtrue, mpred);
	float fp = count_false_positives (mtrue, mpred);
	return tp / (tp + fp);
}


#endif