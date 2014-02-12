

#ifndef PREDICTION_H
#define PREDICTION_H

#include <string>
#include <vector>
#include <iostream>
#include <stack>
#include <utility>
#include <chrono>
#include "rnainterval.h"
#include "vienna.h"
#include "selection.h"

int count_errors(const std::string& model, const std::string& proband)
{
	int count = 0;
	for(int i = 0; i < model.size(); ++i)
		if(model[i] != proband[i])
			++count;
	return count;
}

std::string get_dotbracket(int rna_sz, const std::vector<RNAInterval>& windows, const std::vector<int>& selected)
{
	std::string windowsStruct(rna_sz, '.');
	for(int k = 0; k < selected.size(); ++k)
	{
		int push = windows[selected[k]].left;
		for(int l = 0; l < windows[selected[k]].sstruct.size(); ++l)
			windowsStruct[push+l] = windows[selected[k]].sstruct[l];
	}
	return windowsStruct;
}


std::string get_dotbracket(int rna_sz, const std::vector<RNAInterval>& windows)
{
	std::string windowsStruct(rna_sz, '.');
	for(int k = 0; k < windows.size(); ++k)
	{
		int push = windows[k].left;
		for(int l = 0; l < windows[k].sstruct.size(); ++l)
			windowsStruct[push+l] = windows[k].sstruct[l];
	}
	return windowsStruct;
}

int splat_prediction(const std::vector<int>& splat, const std::string& rna, const std::string& target_sstruct) 
{

	std::vector<RNAInterval> all_windows;

	for(int i = 0; i < splat.size() && splat[i] <= rna.size() / 4; ++i)
	{
		std::vector<RNAInterval> windows = rnal_fold(rna, splat[i]);
		all_windows.insert(all_windows.end(), windows.begin(), windows.end());
	}

	std::vector<int> selected_windows = weighted_activity_selection(all_windows);



	std::string windows_struct = get_dotbracket(rna.size(), all_windows, selected_windows);

	int windows_errors = count_errors(target_sstruct, windows_struct);

	return windows_errors;

}


void two_window_prediction(const std::string& rna, const std::string& target_sstruct)
{
	std::vector<std::vector<RNAInterval> > all_windows;
	std::vector<RNAInterval> windows;
	std::vector<int> selectedWindows;
	std::string windowsStruct;
	std::string best_struct;
	for(int sz = 10; sz < 400; sz++)
		all_windows.push_back(rnal_fold(rna, sz));
	int besti, bestj;
	int min_errors = 9999999;
	for(int i = 0; i < all_windows.size(); ++i)
		for(int j = i; j < all_windows.size(); ++j)
		{
			windows.clear();
			windows.insert(windows.end(), all_windows[i].begin(), all_windows[i].end());
			if(i != j)
				windows.insert(windows.end(), all_windows[j].begin(), all_windows[j].end());
			selectedWindows = weighted_activity_selection(windows);
			std::vector<int> selectedWindows = weighted_activity_selection(windows);
			std::string windowsStruct = get_dotbracket(rna.size(), windows, selectedWindows);
			int windowsError = count_errors(target_sstruct, windowsStruct);
			if(windowsError < min_errors)
			{
				min_errors = windowsError;
				besti = i;
				bestj = j;
				best_struct = windowsStruct;
			}
			
		}
	std::cout << "Windows errors = " << min_errors << " at i = " << besti << ", j = " << bestj << std::endl;
	std::cout << best_struct << std::endl;
	std::cout << "---------------------------------" << std::endl;
}

int min(int a, int b)
{
	return (a < b) ? a : b;
}

int max(int a, int b)
{
	return (a > b) ? a : b;
}

void elastic_prediction(const std::string& rna, const std::string& target_sstruct)
{
	std::vector<RNAInterval> windows = rnal_fold(rna, rna.size()/20);
	std::vector<RNAInterval> elastic_windows;
	for(int i = 0; i < windows.size(); ++i)
	{
		elastic_windows.push_back(windows[i]);
		RNAInterval expanded = zuker_fold(rna, max(0, windows[i].left-4), min(rna.size()-1, windows[i].right+4));
		//RNAInterval shrunk = zuker_fold(rna, windows[i].left+5, windows[i].right-5);
		//double prev_fe;
		//prev_fe = windows[i].score;
		//int size_mod = 5;
		double proportional_score_expanded = (double) expanded.score / (expanded.right - expanded.left + 1);
		double proportional_score_prev = (double) windows[i].score / (windows[i].right - windows[i].left + 1);

		for(int j = 8; 
			proportional_score_expanded > proportional_score_prev * 1.1 && j < rna.size()/4;
			j *= 2)
		{
			proportional_score_prev = proportional_score_expanded;
			elastic_windows.push_back(expanded);
			expanded = zuker_fold(rna, max(0, windows[i].left-j), min(rna.size()-1, windows[i].right+j));
			proportional_score_expanded = (double) expanded.score / (expanded.right - expanded.left + 1);
		}

		// while(expanded.score > prev_fe && size_mod < rna.size()/5)
		// {
		// 	size_mod *= 2;
		// 	prev_fe = expanded.score;
		// 	elastic_windows.push_back(expanded);
		// 	expanded = zuker_fold(rna, max(0, expanded.left-size_mod), min(rna.size()-1, expanded.right+size_mod));
		// }
		// prev_fe = windows[i].score;
		// size_mod = 5;
		// while(shrunk.score > prev_fe)
		// {
		// 	size_mod *= 2;
		// 	prev_fe = shrunk.score;
		// 	elastic_windows.push_back(shrunk);
		// 	if (shrunk.right - shrunk.left < size_mod + 5)
		// 		break;
		// 	std::cout << shrunk.left << " " << shrunk.right << std::endl;
		// 	shrunk = zuker_fold(rna, shrunk.left+size_mod, shrunk.right-size_mod);
		// }
		// std::cout << elastic_windows.size() << std::endl;
	}
	std::vector<int> selected_windows = weighted_activity_selection(elastic_windows);
	std::string windowsStruct = get_dotbracket(rna.size(), elastic_windows, selected_windows);
	int windowsError = count_errors(target_sstruct, windowsStruct);
	std::cout << "Windows errors = " << windowsError << " total elastic windows = " << elastic_windows.size() << std::endl;
	std::cout << windowsStruct << std::endl;
	std::cout << "---------------------------------" << std::endl;


}

void consultation_prediction(const std::string& rna, const std::string& target_sstruct)
{
	RNAInterval zuker_struct = zuker_fold(rna, 0, rna.size()-1);
	std::stack<int> s;
	std::vector<int> stems;
	for (int i = 0; i < zuker_struct.sstruct.size(); ++i)
	{
		if (zuker_struct.sstruct[i] == '(')
			s.push(i);
		else if(zuker_struct.sstruct[i] == ')')
		{
			int stem_size = i - s.top() + 1;
			s.pop();
			stems.push_back(stem_size);
		}
	}
	std::sort(stems.begin(), stems.end());
	int sum = 0;
	for(int i = 0; i < stems.size(); ++i)
		sum += stems[i];
	int med_stem_sz = stems[stems.size()/2 + stems.size()/4];
	std::cout << "Median Stem Size = " << med_stem_sz << std::endl;
	//std::vector<RNAInterval> min_windows = rnal_fold(rna, min_stem_sz);
	std::vector<RNAInterval> all_windows = rnal_fold(rna, med_stem_sz);
	//all_windows.insert(all_windows.end(), min_windows.begin(), min_windows.end());
	std::vector<int> selected_windows = weighted_activity_selection(all_windows);
	std::string windows_struct = get_dotbracket(rna.size(), all_windows, selected_windows);
	int windows_errors = count_errors(target_sstruct, windows_struct);
	std::cout << "Windows errors = " << windows_errors << std::endl;
	std::cout << windows_struct << std::endl;
	std::cout << "---------------------------------" << std::endl;

}



void incremental_prediction(const std::string& rna, const std::string& target_sstruct)
{

	const double score_threshold = 3.2;
	const int min_sz = 40, max_sz = 400, sz_inc = 30;

	std::vector<std::pair<RNAInterval, int> > selected;

	int broken_windows = 0;

	int fold_epoch = 0;

	for(int sz = min_sz; sz <= max_sz; sz += sz_inc)
	{
		//std::cerr << "Doing size " << sz << std::endl;
		std::vector<RNAInterval> sz_windows = rnal_fold(rna, sz);
		for(int i = 0; i < sz_windows.size(); ++i)
		{
			std::vector<int> to_remove;
			double total_score_loss = 0.0;
			for(int j = 0; j < selected.size(); ++j)
				if(selected[j].first.overlaps(sz_windows[i]))
				{
					total_score_loss += selected[j].first.score + (fold_epoch - selected[j].second) * score_threshold;
					to_remove.push_back(j);
				}
			if(total_score_loss < sz_windows[i].score)
			{
				broken_windows++;
				int swap_pos = selected.size()-1;
				for(int j = 0; j < to_remove.size(); ++j)
				{
					selected[to_remove[j]] = selected[swap_pos--];
				}
				for(int j = 0; j < to_remove.size(); ++j)
					selected.pop_back();
				selected.push_back(std::make_pair(sz_windows[i], fold_epoch));
			}
		}
		++fold_epoch;
	}

	std::vector<RNAInterval> windows;
	for(int i = 0; i < selected.size(); ++i)
		windows.push_back(selected[i].first);


	std::string windows_struct = get_dotbracket(rna.size(), windows);
	int windows_errors = count_errors(target_sstruct, windows_struct);
	std::cout << "Windows errors = " << windows_errors << " broken windows = " << broken_windows << std::endl;
	std::cout << windows_struct << std::endl;
	std::cout << "---------------------------------" << std::endl;

}

#endif