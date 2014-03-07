
#ifndef SELECTION_H
#define SELECTION_H

#include <vector>
#include <algorithm>
#include "../common/vienna.h"
#include "../common/rna_util.h"
#include "precomputed_windows.h"

// The type signature of an RNA window selection algorithm 
typedef std::vector<int> (*t_selection_algorithm)(std::vector<RNAInterval>&);


// Given a selection algorithm and a collection of precomputed windows
// This function will perform a test on the accuracy of the function
void selection_algorithm_test (t_selection_algorithm algo, std::vector<PrecomputedWindows>& precomp)
{
	float avg_f1 = 0.0;
	for (int i = 0; i < precomp.size(); ++i)
	{
		int best_window_index;
		float bestf1 = -1.0;
		for (int j = 0; j < precomp[i].windows.size(); ++j)
		{
			std::vector<int> selected = algo (precomp[i].windows[j]);
			std::string sstruct = get_dotbracket (precomp[i].rna.size(), precomp[i].windows[j], selected);
			float f1score = calc_f1score (precomp[i].actual_sstruct, sstruct);
			if (f1score > bestf1)
			{
				bestf1 = f1score;
				best_window_index = j;
			}
		}

		std::cout << precomp[i].name << ": Best F1 score = " << bestf1 << " at window size = " 
			<< (best_window_index + MIN_WINDOW_SIZE) << " for RNA size = " << precomp[i].rna.size() << std::endl;

		// this commented line was used to produce spreadsheet compatible output
		// std::cout << bestf1 << "\t" << (best_window_index + MIN_WINDOW_SIZE) << "\t" << precomp[i].rna.size() << std::endl;

		avg_f1 += bestf1;
	}
	avg_f1 = avg_f1 / precomp.size();
	std::cout << "Average best F1 score = " << avg_f1 << std::endl;
}

std::vector<int> bottom_up_selection(std::vector<RNAInterval>& intervals)
{
	//sort by size, asc
	std::sort(intervals.begin(), intervals.end(), size_comp);

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

std::vector<int> top_down_selection(std::vector<RNAInterval>& intervals)
{
	//sort by size, desc
	std::sort(intervals.begin(), intervals.end(), size_comp);
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

std::vector<int> greedy_MFE_selection(std::vector<RNAInterval>& intervals)
{
	//sort by score (abs(free energy)), desc
	std::sort(intervals.begin(), intervals.end(), fe_comp);
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
	std::sort(intervals.begin(), intervals.end());
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