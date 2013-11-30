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
	bool compatibleWith(const RNAInterval& other) const
	{
		return other.right < this->left;
	}
	bool operator< (const RNAInterval& other) const
	{
		return right < other.right;
	}
};

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
	int sz = r-l+1;
    char* ss = new char[sz+1];
    double score = fold(rna.substr(l, sz).c_str(), ss);
    /*
    for(int i = 0; i < sz; ++i)
    {
    	if(ss[i] == '.')
    		score += 1.5;
    	else
    		break;
    }
    for(int i = sz-1; i >= 0; ++i)
    {
    	if(ss[i] == '.')
    		score += 1.5;
    	else
    		break;
    }
    */
    return RNAInterval(l, r, -score, string(ss));
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


int main()
{
	
	string rna;
	string targetStructure;
	while(cin.good())
	{
		
		cin >> rna >> targetStructure;
		//normal RNAfold
		RNAInterval vanilla = zuker(rna, 0, rna.size()-1);
		int vanillaErrors = countErrors(targetStructure, vanilla.sstruct);
		cout << "RNA size: " << rna.size() << endl;
		vector<RNAInterval> windows;
		int minErrors = 999999999;
		int size;
		vector<RNAInterval> bestWindows;
		for(int sz = 40; sz < 400; ++sz)
		{
			windows = RNALfold(rna, sz);
			//cout << "Finished sliding windows." << endl;
			string windowsStruct(rna.size(), '.');
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
			if(windowsError < minErrors)
			{
				minErrors = windowsError;
				bestWindows.clear();
				for(int i = 0; i < selectedWindows.size(); ++i)
					bestWindows.push_back(windows[selectedWindows[i]]);
				size = sz;
			}
		}
		cout << "Zuker errors: " << vanillaErrors << endl;
		if(minErrors < vanillaErrors)
		{
			cout << "Windows was better at size " << size << " with " << minErrors << " errors." << endl;
			/*
			cout << "Windows chosen... " << endl;
			for(int i = 0; i < bestWindows.size(); ++i)
			{
				cout << "left index: " << bestWindows[i].left << "  right index: " << bestWindows[i].right << endl;
				cout << bestWindows[i].sstruct << " " << "   fe: " << bestWindows[i].score << endl;
			}
			*/
		}
		else
			cout << "Windows was NOT better!" << " with a minimum of " << minErrors << " errors." << endl;
		cout << "---------------------------------" << endl;
		


	}
	
	return 0;
}
