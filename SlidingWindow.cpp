#include <iostream>
#include <vector>
#include <algorithm>
#include <string>


//include external lib
extern "C" 
{
    #include "vienna/fold.h"
    #include "vienna/Lfold.h"
}

using namespace std;

class RNAInterval
{
public:
	unsigned int left, right;
	double score;
	string sstruct;
	RNAInterval(unsigned int l, unsigned int r, double s, const string& ss)
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
	bool compatibleWith(const RNAInterval& other) const
	{
		return other.right < this->left || other.left > this->right;
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

vector<int> bottomUpSelection(vector<RNAInterval>& intervals)
{
	//sort by size, asc
	sort(intervals.begin(), intervals.end(), size_comp);
	vector<int> chosen;
	for(int i = 0; i < intervals.size(); ++i)
	{
		bool choose = true;
		for(int j = 0; j < chosen.size(); ++j)
		{
			if(!intervals[i].compatibleWith(intervals[chosen[j]]))
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

vector<int> topDownSelection(vector<RNAInterval>& intervals)
{
	//sort by size, asc
	sort(intervals.begin(), intervals.end(), size_comp);
	reverse(intervals.begin(), intervals.end());
	vector<int> chosen;
	for(int i = 0; i < intervals.size(); ++i)
	{
		bool choose = true;
		for(int j = 0; j < chosen.size(); ++j)
		{
			if(!intervals[i].compatibleWith(intervals[chosen[j]]))
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

vector<int> greedyMFESelection(vector<RNAInterval>& intervals)
{
	//sort by score (abs(free energy)), asc
	sort(intervals.begin(), intervals.end(), fe_comp);
	reverse(intervals.begin(), intervals.end());
	vector<int> chosen;
	for(int i = 0; i < intervals.size(); ++i)
	{
		bool choose = true;
		for(int j = 0; j < chosen.size(); ++j)
		{
			if(!intervals[i].compatibleWith(intervals[chosen[j]]))
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

vector<int> weightedActivitySelection(vector<RNAInterval>& intervals)
{
	//sort by right end point
	sort(intervals.begin(), intervals.end());
	vector<int> q(intervals.size()+1, 0);
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
	vector<double> dp(intervals.size()+1);
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
	vector<int> traced;
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
	reverse(traced.begin(), traced.end());
	return traced;
}


RNAInterval zuker(const string& rna, const int l, const int r)
{
    //get the interval string and calculate the score
	int sz = r-l+1;
    char* ss = new char[sz+1];
    double score = fold(rna.substr(l, sz).c_str(), ss);
    //trim the string of danging ends
	//int ll, rr;
    //for(ll = l; ll <= r && rna[ll] == '.'; ++ll);
    //for(rr = r; r >= 0 && rna[rr] == '.'; --rr);
    //cerr << ss << endl;
    return RNAInterval(l, r, -score, ss);
}

vector<RNAInterval> zukerMultiFold(const string& rna, const int windowSize) {
	vector<RNAInterval> ints;
	for(int i = 0; i < rna.size()-windowSize; ++i) {
		//cerr << "Zuker: " << i << endl;
		ints.push_back(zuker(rna, i, i+windowSize));
	}
	return ints;
}

vector<RNAInterval> RNALfold(const string& rna, int windowSize)
{
    char* ss = new char[rna.size()+1];
    struct mnode* n = mLfold(rna.c_str(), ss, windowSize);
    vector<RNAInterval> tmp;
    string s;
    while(n->next != 0) 
    {
    	s = string(n->sstruct);
        tmp.push_back(RNAInterval(n->left-1, n->left+s.size()-1, -(n->fe), s));
        //cout << tmp[tmp.size()-1].left << " " << tmp[tmp.size()-1].right << " " << tmp[tmp.size()-1].score << endl;
        //int span = tmp[tmp.size()-1].right - tmp[tmp.size()-1].left;
        //if(span > windowSize)
		//	cerr << "PROBLEM" << endl;
        n = n->next;
    }
    return tmp;
}

int countErrors(const string& model, const string& proband)
{
	int count = 0;
	for(int i = 0; i < model.size(); ++i)
		if(model[i] != proband[i])
			++count;
	return count;
}

void testSelection()
{
	std::vector<RNAInterval> v;
	v.push_back(RNAInterval(1, 2, 12.4));
	v.push_back(RNAInterval(0, 29, 12.4));
	v.push_back(RNAInterval(6, 66, 8.4));
	v.push_back(RNAInterval(33, 66, 8.4));
	v.push_back(RNAInterval(67, 68, 8.4));
	v.push_back(RNAInterval(69, 70, 8.4));
	v.push_back(RNAInterval(6, 71, 22.0));
	vector<int> res = weightedActivitySelection(v);
	for(int i = 0; i < res.size(); ++i)
	{
		cout << v[res[i]].left << " " << v[res[i]].right << " " << v[res[i]].score << endl;
	}
}

int countUnbonded(const string& structure)
{
	int count = 0;
	for(int i = 0; i < structure.size(); ++i)
	{
		if(structure[i] == '.')
			++count;
	}
	return count;
}


void incremental_prediction(const string& rna, const string& target_sstruct)
{
	vector<RNAInterval> current_sstruct;
	for(int i = 10; i < 400; i*=2)
	{
		vector<RNAInterval> windows = RNALfold(rna, i);
		for(int j = 0; j < windows.size(); ++j)
		{
			int score_loss = 0;
			vector<int> to_remove;
			for(int k = 0; k < current_sstruct.size(); ++k)
				if (!current_sstruct[k].compatibleWith(windows[j]))
				{
					score_loss += current_sstruct[k].score;
					to_remove.push_back(k);
				}
			if (windows[j].score > score_loss)
			{
				for(int k = 0; k < to_remove.size(); ++k)
					current_sstruct.erase(current_sstruct.begin()+to_remove[k]);
				current_sstruct.push_back(windows[j]);
			}
		}
	}
	string sstruct(rna.size(), '.');
	for(int k = 0; k < current_sstruct.size(); ++k)
	{
		int push = current_sstruct[k].left;
		string str = current_sstruct[k].sstruct;
		for(int l = 0; l < str.size(); ++l)
			sstruct[push+l] = str[l];
	}
	int errors = countErrors(target_sstruct, sstruct);
	cout << "Windows errors = " << errors << endl;
	cout << sstruct << endl;
	cout << "---------------------------------" << endl;

}

void two_window_prediction(const string& rna, const string& target_sstruct)
{
	vector<vector<RNAInterval> > all_windows;
	vector<RNAInterval> windows;
	vector<int> selectedWindows;
	string windowsStruct;
	string best_struct;
	for(int sz = 10; sz < 400; sz++)
		all_windows.push_back(RNALfold(rna, sz));
	int besti, bestj;
	int min_errors = 9999999;
	for(int i = 0; i < all_windows.size(); ++i)
		for(int j = i; j < all_windows.size(); ++j)
		{
			windows.clear();
			windows.insert(windows.end(), all_windows[i].begin(), all_windows[i].end());
			if(i != j)
				windows.insert(windows.end(), all_windows[j].begin(), all_windows[j].end());
			selectedWindows = weightedActivitySelection(windows);
			vector<int> selectedWindows = weightedActivitySelection(windows);
			string windowsStruct(rna.size(), '.');
			for(int k = 0; k < selectedWindows.size(); ++k)
			{
				int push = windows[selectedWindows[k]].left;
				string str = windows[selectedWindows[k]].sstruct;
				for(int l = 0; l < str.size(); ++l)
					windowsStruct[push+l] = str[l];
			}
			int windowsError = countErrors(target_sstruct, windowsStruct);
			if(windowsError < min_errors)
			{
				min_errors = windowsError;
				besti = i;
				bestj = j;
				best_struct = windowsStruct;
			}
			
		}
	cout << "Windows errors = " << min_errors << " at i = " << besti << ", j = " << bestj << endl;
	cout << best_struct << endl;
	cout << "---------------------------------" << endl;
}


int main()
{
	
	string rna, name;
	string targetStructure;
	while(cin.good())
	{
		
		cin >> name >> rna >> targetStructure;
		cout << "SSTRAND ID = " << name << endl;
		//normal RNAfold
		RNAInterval vanilla = zuker(rna, 0, rna.size()-1);
		int vanillaErrors = countErrors(targetStructure, vanilla.sstruct);
		cout << "RNA size: " << rna.size() << " with " << vanillaErrors << " errors." << endl;
		// vector<RNAInterval> windows;
		// int minErrors = 999999999;
		// int size;
		// vector<RNAInterval> bestWindows;
		// string windowsStruct;
		/*
		for(int sz = 40; sz < 400; ++sz)
		{
			windowsStruct = string(rna.size(), '.');
			windows = RNALfold(rna, sz);
			//cout << "Finished sliding windows." << endl;
			//cout << "Starting activity selection... " << endl;
			vector<int> selectedWindows = weightedActivitySelection(windows);
			//cout << "Finished activity selection." << endl;
			for(int i = 0; i < selectedWindows.size(); ++i)
			{
				int push = windows[selectedWindows[i]].left;
				string str = windows[selectedWindows[i]].sstruct;
				for(int j = 0; j < str.size(); ++j)
					windowsStruct[push+j] = str[j];
			}
			int windowsError = countErrors(targetStructure, windowsStruct);
			//cout << "Trying window size of " << sz << " had " << windowsError << " errors." << endl;
			if(windowsError < minErrors)
			{
				minErrors = windowsError;
				bestWindows.clear();
				for(int i = 0; i < selectedWindows.size(); ++i)
					bestWindows.push_back(windows[selectedWindows[i]]);
				size = sz;
			}
			cout << "Size = " << sz << " errors = " << windowsError << " unbonded = " << countUnbonded(windowsStruct) << endl;
		}
		*/

		vector<vector<RNAInterval> > all_windows;
		vector<RNAInterval> windows;
		vector<int> selectedWindows;
		string windowsStruct;
		string best_struct;
		for(int sz = 10; sz < 400; sz++)
			all_windows.push_back(RNALfold(rna, sz));
		int besti, bestj;
		int min_errors = 9999999;
		for(int i = 0; i < all_windows.size(); ++i)
			for(int j = i; j < all_windows.size(); ++j)
			{
				windows.clear();
				windows.insert(windows.end(), all_windows[i].begin(), all_windows[i].end());
				if(i != j)
					windows.insert(windows.end(), all_windows[j].begin(), all_windows[j].end());
				selectedWindows = weightedActivitySelection(windows);
				vector<int> selectedWindows = weightedActivitySelection(windows);
				string windowsStruct(rna.size(), '.');
				for(int k = 0; k < selectedWindows.size(); ++k)
				{
					int push = windows[selectedWindows[k]].left;
					string str = windows[selectedWindows[k]].sstruct;
					for(int l = 0; l < str.size(); ++l)
						windowsStruct[push+l] = str[l];
				}
				int windowsError = countErrors(targetStructure, windowsStruct);
				if(windowsError < min_errors)
				{
					min_errors = windowsError;
					besti = i;
					bestj = j;
					best_struct = windowsStruct;
				}
				
			}
		cout << "Windows errors = " << min_errors << " at i = " << besti << ", j = " << bestj << endl;
		cout << best_struct << endl;

		// cout << "Zuker errors: " << vanillaErrors << endl;
		
		// if(minErrors < vanillaErrors)
		// {
		// 	cout << "Windows was better at size " << size << " with " << minErrors << " errors." << endl;
			/*
			cout << "Windows chosen... " << endl;
			for(int i = 0; i < bestWindows.size(); ++i)
			{
				cout << "left index: " << bestWindows[i].left << "  right index: " << bestWindows[i].right << endl;
				cout << bestWindows[i].sstruct << " " << "   fe: " << bestWindows[i].score << endl;
			}
			*/
		// }
		// else
		// 	cout << "Windows was NOT better!" << " with a minimum of " << minErrors << " errors." << endl;

		// cout << "Window Predicted Structure: " << endl << windowsStruct << endl;
		cout << "---------------------------------" << endl;
		
		cin >> ws;

	}
	
	return 0;
}
