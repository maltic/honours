#include "GAevolve.hpp"



void GAfold_error(const char* msg) {
    std::cerr << msg << std::endl;
    throw 1;
}

//picks a random number between two points
int rand_between(const int a, const int b) {
    if (b == 0) GAfold_error("Cant have a random number between something and a 0");
    int range = b-a;
    return rand() % range + a;
}

bool valid_bond(const char a, const char b) {
    return  (a == 'A' && b == 'U')
         || (a == 'U' && b == 'A')
         || (a == 'G' && b == 'C')
         || (a == 'C' && b == 'G')
         || (a == 'U' && b == 'G')
         || (a == 'G' && b == 'U');
}

//given a starting open bond point, finds its closing bond using the parenthesis property
int find_closing_bond(const std::string& sstruct, const int a) {
    if(sstruct[a] != OPEN_BOND) { 
        GAfold_error("Error: called find_closing_bond, but didn't start on a bond!");
    }
    int count = 1;
    for(int i = a+1; i < sstruct.size(); ++i) {
        if(sstruct[i] == OPEN_BOND)
            ++count;
        else if(sstruct[i] == CLOSE_BOND) {
            --count;
            if(count == 0) return i;
        }
    }
    return -1;
}

//makes a pair of locations bonded
void make_bond_pair(std::string& sstruct, const int a, const int b) {
    sstruct[a] = OPEN_BOND;
    sstruct[b] = CLOSE_BOND;
}

//checks if a pair of locations can be made into a bonded pair
bool possible_bond_pair(std::string& sstruct, const int a, const int b) {
    if (a >= b) return false;
    return sstruct[a] == NO_BOND && sstruct[b] == NO_BOND;
}

//checks if a pair of locations are already a bonded pair
bool is_bond_pair(std::string& sstruct, const int a, const int b) {
    if (a >= b) return false;
    return sstruct[a] == OPEN_BOND && sstruct[b] == CLOSE_BOND;
}

std::string apply_mod(const StructMod& mod, const std::string& sstruct) {
    return mod.apply_to(sstruct);
}

//inserts a new random bond into a secondary structure string
StructMod insert_rand_bond(const std::string& sstruct, const std::string& rna) {
    
    std::vector<StructMod> insertions;
    int count;
    for(int i = 0; i < sstruct.size(); ++i) {
        if(sstruct[i] != NO_BOND) continue;
        count = 0;
        for(int j = i+1; j < sstruct.size(); ++j) {
            if(sstruct[j] != NO_BOND) break;
            ++count;
            if(count > 2 && valid_bond(rna[i], rna[j])) {
                insertions.push_back(StructMod(INSERTION, i, j));
            }
        }
    }
    //if there are no possible insertions return a blank mod
    if(insertions.size() == 0) return StructMod();

    return insertions[rand_between(0, insertions.size())];
}

//deletes a randomly picked bond from a secondary structure string
StructMod delete_rand_bond(const std::string& sstruct) {
    //find all open bond locations
    std::vector<int> open_locs;
    for(int i = 0; i < sstruct.size(); ++i) {
        if(sstruct[i] == OPEN_BOND) open_locs.push_back(i);
    }
    
    //short circuit if there are no bonds to delete
    if (open_locs.size() == 0) return StructMod(); //return blank modification
    
    //pick the opening bond point
    int a = open_locs[rand_between(0, open_locs.size())];
    //now we need to find the closing bond
    int b = find_closing_bond(sstruct, a);
    
    if(b == -1) GAfold_error("Couldn't find closing bond in delete_rand_bond");

    return StructMod(DELETION, a, b);
}

std::string evolve_rna_sstruct(const std::string& rna, const int maxInsertions, const int gens) {

    std::string blank(rna.size(), '.');
    Genome best(blank, 0.0);

    for(int gen = 0; gen < gens; ++gen) {

        float fitness = 0.0;
        std::string sstruct = blank;
        for(int i = 0; i < maxInsertions; ++i) {
            StructMod sm = insert_rand_bond(sstruct, rna);
            if(sm.type == INSERTION) {
                fitness += energy_of_move(rna.c_str(), sstruct.c_str(), sm.left+1, sm.right+1);
            }
            else if (sm.type == DELETION) {
                fitness += energy_of_move(rna.c_str(), sstruct.c_str(), -(sm.left+1), -(sm.right+1));
            }
            sstruct = sm.apply_to(sstruct);
            
            if(fitness < best.fitness) {
                best.fitness = fitness;
                best.sstruct = sstruct;
            }
        }

        //std::cout << fitness << " " << sstruct << std::endl;

    }
    return best.sstruct;
}
