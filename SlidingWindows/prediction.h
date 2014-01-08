

#ifndef PREDICTION_H
#define PREDICTION_H

#include <string>
#include <vector>
#include <iostream>
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
			std::string windowsStruct(rna.size(), '.');
			for(int k = 0; k < selectedWindows.size(); ++k)
			{
				int push = windows[selectedWindows[k]].left;
				std::string str = windows[selectedWindows[k]].sstruct;
				for(int l = 0; l < str.size(); ++l)
					windowsStruct[push+l] = str[l];
			}
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

#endif