#include <iostream>
#include <string>

using namespace std;

#define maxSize 1000
#define POSINF 9999999

int dp[maxSize][maxSize];

string rnaSeq;

int max(int a, int b) { return a > b ? a : b; }

bool validBond(char a, char b) {
    return (a == 'A' && b == 'U')
            || (a == 'U' && b == 'A')
            || (a == 'G' && b == 'C')
            || (a == 'C' && b == 'G')
            || (a == 'U' && b == 'G')
            || (a == 'G' && b == 'U');
}

int nussWeight(int a, int b) {
    validBond(rnaSeq[a], rnaSeq[b]) ? 1 : -99999;
}

int zukWeight(int a, int b) {
    if ((a == 'A' && b == 'U')
            || (a == 'U' && b == 'A')) return -8;
    if ((a == 'G' && b == 'C')
            || (a == 'C' && b == 'G')) return -12;
    if ((a == 'G' && b == 'U')
            || (a == 'U' && b == 'G')) return -4;
    return 0;
}

int nussinovs(int i, int j) {
    if(i >= j) return 0;
    if(dp[i][j] != -1) return dp[i][j];
    int best = nussinovs(i+1, j-1) + nussWeight(i, j);
    for(int k = i; k < j; ++k) 
        best = max(nussinovs(i,k) + nussinovs(k+1,j), best);
    return dp[i][j] = best;
}



int zuker(int i, int j) {
    if(i >= j) return POSINF;
    if(dp[i][j] != -1) return dp[i][j];
    int best = zuker(i+1, j-1) + zukWeight(i, j);
    for(int k = i; k < j; ++k) 
        best = min(zuker(i,k) + zuker(k+1,j), best);
    return dp[i][j] = best;
}

int main() {
    while(cin>>rnaSeq) {
        for(int i = 0; i < maxSize; ++i)
            for(int j = 0; j < maxSize; ++j)
                dp[i][j] = -1;
        cout << nussinovs(0, rnaSeq.size()-1) << endl;
    }
    return 0;
}
