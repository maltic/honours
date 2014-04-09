
#ifndef SELECTION_H
#define SELECTION_H

#include <vector>
#include <algorithm>
#include "../common/vienna.h"
#include "../common/rna_util.h"
#include "precomputed_windows.h"

// The type signature of an RNA window selection algorithm 
typedef std::vector<int> (*t_selection_algorithm)(std::vector<RNAInterval>&);

// A SelectionTester object can be used to compare various selection algorithms
// Because I wanted to be able to test using an arbitrary number of combinatorically best windows
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
	void selection_algorithm_test (const int n, const t_selection_algorithm algo, std::vector<PrecomputedWindows>& precomp)
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
			for (int & i : best_windows)
				std::cout << (i + MIN_WINDOW_SIZE) << " ";
			std::cout << "\t" << precomp[i].rna.size() << std::endl;

			avg_f1 += best_f1score;
		}
		avg_f1 = avg_f1 / precomp.size();
		std::cout << "Average best F1 score = " << avg_f1 << std::endl;
	}
};



std::vector<int> bottom_up_selection(std::vector<RNAInterval>& intervals)
{
	//sort by size, asc
	std::sort(intervals.begin(), intervals.end(), rna_int_size_comp);

	std::vector<int> chosen;
	for(int i = 0; i < intervals.size(); ++i)
	{
		bool choose = true;
		for(int j = 0; j < chosen.size(); ++j)
		{
			if ( !intervals[i].compatible_with(intervals[chosen[j]]) )
			{
				choose = false;
				break;
			}
		}
		if(choose)
			chosen.push_back(i);
	}
	return chosen;
}

std::vector<int> top_down_selection(std::vector<RNAInterval>& intervals)
{


	//sort by size, desc
	std::sort(intervals.begin(), intervals.end(), rna_int_size_comp);
	std::reverse(intervals.begin(), intervals.end());

	std::vector<int> chosen;
	for(int i = 0; i < intervals.size(); ++i)
	{
		bool choose = true;
		for(int j = 0; j < chosen.size(); ++j)
		{
			if(!intervals[i].compatible_with(intervals[chosen[j]]))
			{
				choose = false;
				break;
			}
		}
		if(choose)
			chosen.push_back(i);
	}
	return chosen;
}

// also called Score Selection in the Dissertation
std::vector<int> greedy_MFE_selection(std::vector<RNAInterval>& intervals)
{

	//sort by score (free energy * -1), desc
	std::sort(intervals.begin(), intervals.end(), rna_int_fe_comp);
	std::reverse(intervals.begin(), intervals.end());

	std::vector<int> chosen;
	for(int i = 0; i < intervals.size(); ++i)
	{
		bool choose = true;
		for(int j = 0; j < chosen.size(); ++j)
		{
			if(!intervals[i].compatible_with(intervals[chosen[j]]))
			{
				choose = false;
				break;
			}
		}
		if(choose)
			chosen.push_back(i);
	}
	return chosen;
}


std::vector<int> weighted_activity_selection (std::vector<RNAInterval>& intervals)
{
	//sort by right end point
	std::sort ( intervals.begin(), intervals.end() );
	std::vector<int> q(intervals.size()+1, 0);
	//compute q values
	for(int i = 1; i < intervals.size(); ++i)
	{
		//binary search for last compatible index
		int l = 0, r = i-1, mid;
		while(l <= r)
		{
			mid = l + (r-l)/2;
			if(intervals[mid].right < intervals[i].left && !(intervals[mid+1].right < intervals[i].left))
			{
				q[i+1] = mid+1;
				break;
			}
			else if(intervals[i].left > intervals[mid].right)
				l = mid+1;
			else
				r = mid-1;
		}

	}
	std::vector<double> dp(intervals.size()+1);
	//bottom up dp fill
	dp[0] = 0;
	for(int i = 1; i < dp.size(); ++i)
	{
		double iScore = intervals[i-1].score + dp[q[i]];
		if(iScore > dp[i-1])
			dp[i] = iScore;
		else
			dp[i] = dp[i-1];
	}
	//traceback
	int curr = intervals.size();
	std::vector<int> traced;
	while(curr > 0)
	{
		if (intervals[curr-1].score + dp[q[curr]] > dp[curr-1])
		{
			traced.push_back(curr-1);
			curr = q[curr];
		}
		else
			--curr;
	}
	std::reverse(traced.begin(), traced.end());
	return traced;
}

#endif