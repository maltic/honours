//Compile from project root
//g++ src/slidingwindows/main.cpp src/common/vienna/libRNA.a -std=c++11 -O2
#include <iostream>
#include "testing.h"

void MEA_and_zuk()
{
	typedef std::chrono::high_resolution_clock Clock;
	typedef std::chrono::milliseconds milliseconds;
	std::vector<PrecomputedWindows> precomp = load_precomputed_windows (std::cin);
	std::sort (precomp.begin(), precomp.end(), precomputed_windows_size_cmp);
	for (auto & pc : precomp)
	{
		if (pc.rna.size() > 800)
			continue;
		Clock::time_point t0 = Clock::now();
		RNAInterval mea = MEA_fold(pc.rna);
		Clock::time_point t1 = Clock::now();
		milliseconds ms_mea = std::chrono::duration_cast<milliseconds>(t1 - t0);
		t0 = Clock::now();
		RNAInterval zuk = zuker_fold(pc.rna);
		t1 = Clock::now();
		milliseconds ms_zuk = std::chrono::duration_cast<milliseconds>(t1 - t0);
		float mea_f1score = calc_f1score (pc.actual_sstruct, mea.sstruct);
		float zuk_f1score = calc_f1score (pc.actual_sstruct, zuk.sstruct);

		std::cout << pc.rna.size() << '\t' << ms_mea.count() << '\t' << mea_f1score << '\t'
			<< ms_zuk.count() << '\t' << zuk_f1score << std::endl;
	}

}


void print_help()
{
	std::cout << "Some helpful information!" << std::endl;
}




int main(int argc, char **argv)
{
	if (argc < 2)
		print_help();
	else if (argv[1] == "-zukers")
		std::cout << "ahhh yes" << std::endl;
	else
		print_help();
	return 0;
}
