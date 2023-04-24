#include <iostream>

#include "Network.hpp"
#include "Node.hpp"
#include "organism.hpp"
 // Technique officially works on April 12 Wednesday 
using namespace std;

int main(int, char const**) {

    int population_sizes[1]={100};
    int num_pop=1;
    int max_round=20000;
    int no_inputs=2;
    int no_output=2;
    float percentage_size=0.20;

    run_world(population_sizes, num_pop,  max_round,  no_inputs,  no_output,  percentage_size );

    return 0;
}
