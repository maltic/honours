#include <iostream>
#include <string>
#include <cstdlib>
#include "GAevolve.hpp"

int main(int argc, const char* argv[]) {
    if (argc != 4) {
        std::cout <<  "Guide for use: Provide this program with 3 arguments, the rna sequence, maximum number of insertions and no of generations" << std::endl;
        return 1;
    }
    srand (time(NULL));
    std::string rna = argv[1];
    std::string result = evolve_rna_sstruct(rna, atoi(argv[2]), atoi(argv[3]));
    std::cout << result << std::endl;
    std::cout << "free energy: " << energy_of_structure(rna.c_str(), result.c_str(), 0) << std::endl;
    return 0;
}
