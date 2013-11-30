#pragma once
#include <iostream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>

//include the vienna headers
extern "C" {
    #include "vienna/fold.h"
}

//define standard character types for RNA secondary structure string
const char NO_BOND = '.';
const char OPEN_BOND = '(';
const char CLOSE_BOND = ')';
const int POS_INF = 999999999;

enum ModType { NIL, INSERTION, DELETION };
//represents the modification of a secondary structure
class StructMod {
public:
    unsigned int left, right;
    ModType type;

    StructMod(const ModType t, const unsigned int l, const unsigned int r) {
        if (t == INSERTION) {
            this->left = l;
            this->right = r;
            this->type = INSERTION;
        }
        else if (t == DELETION) {
            this->left = l;
            this->right = r;
            this->type = DELETION;
        }
        else {
            this->left = POS_INF;
            this->right = POS_INF;
            this->type = NIL;
        }
    }
    StructMod() {
        this->left = POS_INF;
        this->right = POS_INF;
        this->type = NIL;
    }

    std::string apply_to(const std::string& sstruct) const {
        
        std::string result(sstruct);
        if(this->type == INSERTION) {
            result[this->left] = OPEN_BOND;
            result[this->right] = CLOSE_BOND;
        }
        if(this->type == DELETION) {
            result[this->left] = result[this->right] = NO_BOND;
        }
        return result;

    }
};

void GAfold_error(const char* msg);

//picks a random number between two points
int rand_between(const int a, const int b);

//tells you if a bond is valid
bool valid_bond(const char a, const char b);

//given a starting open bond point, finds its closing bond using the parenthesis property
int find_closing_bond(const std::string& sstruct, const int a);

//makes a pair of locations bonded
void make_bond_pair(std::string& sstruct, const int a, const int b);

//checks if a pair of locations can be made into a bonded pair
bool possible_bond_pair(std::string& sstruct, const int a, const int b);

//checks if a pair of locations are already a bonded pair
bool is_bond_pair(std::string& sstruct, const int a, const int b);

//applies a structual modification to a secondary structure
std::string apply_mod(const StructMod& mod, const std::string& sstruct);

//inserts a new random bond into a secondary structure string
StructMod insert_rand_bond(const std::string& sstruct, const std::string& rna);

//deletes a randomly picked bond from a secondary structure string
StructMod delete_rand_bond(const std::string& sstruct);

struct Genome {
    std::string sstruct;
    float fitness;
    bool operator< (const Genome& i) const { return (this->fitness<i.fitness); }
    Genome(const std::string& ss, float f) {
        this->sstruct = ss;
        this->fitness = f;
    }
};

std::string evolve_rna_sstruct(const std::string& rna, const int maxInsertions, const int gens);
