
#ifndef Dna_hpp
#define Dna_hpp


#include <vector>
#include <unordered_map>
#include "Node.hpp"

using namespace std;

int count_o = 0;



// if the gene is -1 then it is a marker gene expressing that a certain interval has needed
struct genes {
    bool is_parameter;
    bool is_collection;
    bool is_value;
    /* These first three variables are boolean variables that check if the gene is either a parameter collection or value node*/
    
    bool is_active;                     /*This checks if the gene is active, if it is dead it will not be expressed*/
    bool collection_index_mutator;      /*This is the mutator for the index of a collection node, the index is a
                                        the type of operation it will be, for example an index of 0 is addition while 1 is the reLu activation */
    int collection_index;               // Type of collection node 
    int expression;                     //if this variable is -1 then the gene will be a mutator node, mutator node are 
                                        //special location within the genome that designate areas where input and output nodes will be 
    int synapse_time;                   // This is the time required for a specific type of expressed node to fire 
    int collective_time;                // This is the time given to a collective node, to gather inputs, after which iw will fire
    int current_time = 0;               // current synapse time, versus max synapse time
    int abs_time = 0;                   // absolute time since gene was released to the wil
    string hash;                        // the hash of the gene, when the program is multithreaded, the hash will have to be randomly generated so 
                                        //the likelihood of any two genes having the same hash is small
    int level;
    vector<string> past_forms;
    Node* ptr;                          /* The ptr to the node*/
    double mutator_connections;         //The mutator connections, 
    vector<string> node_connections;    //The string of hashes of other nodes(may include itself), that this node is connected too
    std::unordered_map<std::string, int> node_connections_found; // A map to check if a certain node has already been connected too
    float parameter_index;              //parameter value, if the gene is a parameter 
    genes& operator=(const genes& right) {
        node_connections.assign(right.node_connections.begin(), right.node_connections.end());
        past_forms.assign(right.past_forms.begin(), right.past_forms.end());
        node_connections_found = right.node_connections_found;
        ptr = nullptr;
        is_parameter = right.is_parameter;
        is_collection = right.is_collection;
        is_value = right.is_value;
        is_active = right.is_active;
        collection_index_mutator = right.collection_index_mutator;
        collection_index = right.collection_index;
        level=right.level;
        expression = right.expression;
        hash = right.hash;
        expression = right.expression;
        mutator_connections = right.mutator_connections;
        synapse_time = right.synapse_time;
        current_time = right.current_time;
        parameter_index = right.parameter_index;
        collective_time = right.collective_time;
        abs_time = right.abs_time;

        return *this;
    }
};

struct dna {  // dominant strand of the gene is the one that will be expressed when passed down
    vector<genes> dominant_strand;
    // vector<genes> silent_strand; do later
    vector<double> dominant_strand_mutators;
    std::unordered_map<std::string, Node*> map;// map from a hash of a gene to it's node ptr
    // vector<double> silent_strand_mutators; do later
    double mutator_add;// Addition of a gene
    double mutator_delete;//deletion of a gene 
    bool is_dead;//genome has noa active genes 

    dna& operator=(const dna& right) {
        dominant_strand.clear();
        dominant_strand_mutators.clear();
        if (right.dominant_strand.size() > 0) {
            dominant_strand.assign(right.dominant_strand.begin(), right.dominant_strand.end());
            dominant_strand_mutators.assign(right.dominant_strand_mutators.begin(),
                                            right.dominant_strand_mutators.end());
        }
        mutator_add = right.mutator_add;
        mutator_delete = right.mutator_delete;
        is_dead = right.is_dead;

        map.clear();

        return *this;
    }
};

string generateHash() {
    // hash<long long> hasher;  // define the hash function
    // std::random_device rd;
    // std::mt19937 mt(rd());
    // std::uniform_int_distribution<long long> dist(0, LLONG_MAX);
    // long long random_number = dist(mt);
    // long long random_number1 = dist(mt);

    // //        std::hash<long long> hasher;
    // //        std::size_t hash = hasher(random_number);
    // //        std::stringstream ss;
    // //        ss << hash;

    count_o++;// right now it is just generating the count since the start of the program
    return to_string(count_o);
}

// this function generates a random number given a probability of change, suppose you give it a number 0 and the probability of
// change is 0.001 so after 1000 tries you might have changed 0 by a step of change(the variable step of change)

void randomOutput(double& val, double probability_of_change, double step_of_change, double lower_bound,
                         double higher_bound) {
    float division = 1/probability_of_change;        // this i the division for the random number
    float random = uniformTest(0,division);  // our random number
    if (random==0) {  // if the guess is within the range then the addition will be added
        int direction = uniformTest(-1, 1);
        double temp_val = val;
        temp_val = temp_val + direction * step_of_change;
        if (temp_val < lower_bound || temp_val > higher_bound) {
        } else {
            val = temp_val;
        }
    }
}
// just gives you random gene from a choice a choice of two genes, and changes the value of a reference variable return value
genes& randomGene(genes& first_gene, genes& second_gene, double first_mutator, double second_mutator,
                         double& return_value, double probability_of_change) {
    double random = uniformTest(0, 1);                 // our random number
    if (random == 0) {  // if the guess is within the range then the addition will be added
        return_value = first_mutator;
        return first_gene;
    } else {
        return_value = second_mutator;
        return second_gene;
    }
}

void mutate_connection(double mutator, dna& genome, genes& x1) {
    double value = 0;
    value = uniformTest(0, 10);
        float erase = uniformTestRange(0, 1,1000);
        if (erase > 0.3 && x1.node_connections.size() > 0) {
            double index = uniformTest(0, x1.node_connections.size() - 1);

            x1.node_connections_found.erase(x1.node_connections.at(index));
            x1.node_connections.erase(x1.node_connections.begin() + index);
        } else if (genome.dominant_strand.size() > 0) {
            double index = uniformTest(0, genome.dominant_strand.size() - 1);
            auto it = x1.node_connections_found.find(genome.dominant_strand.at(index).hash);
            if (it == x1.node_connections_found.end() && genome.dominant_strand.at(index).expression != -1&&
            genome.dominant_strand.at(index).level<=x1.level+1&& x1.level-1<=genome.dominant_strand.at(index).level) {
                x1.node_connections.push_back(genome.dominant_strand.at(index).hash);
                x1.node_connections_found[genome.dominant_strand.at(index).hash] = 1;
            }

        }
}

genes* make_random_gene(dna& genome) {
    genes* new_gene = new genes;
    double param = uniformTest(0,100);
    // if (param < 10) {
    //     new_gene->is_parameter = true;
    //     new_gene->is_value = false;
    //     new_gene->is_collection = false;
    //     new_gene->expression = 0;
    //     new_gene->parameter_index = uniformTestRange(-1, 1, 1000);
    //     new_gene->synapse_time = uniformTest(1, 9);

    // } 
     if (param < 70) {
        new_gene->is_parameter = false;
        new_gene->is_value = true;
        new_gene->is_collection = false;
        new_gene->expression = 1;
        new_gene->synapse_time = uniformTestRange(0,9,1000);

    // } else if (param < 90) {
    //     new_gene->is_parameter = false;
    //     new_gene->is_value = false;
    //     new_gene->is_collection = true;
    //     new_gene->expression = 2;
    //     new_gene->synapse_time = uniformTest(1, 9);
    //     new_gene->collection_index = uniformTest(0, 6);
    //     new_gene->collection_index_mutator = uniformTestRange(0, 0.1, 1000);
    //     new_gene->collective_time = 2;

    // } else {
        new_gene->expression = -1;
    }
    new_gene->is_active = true;
    new_gene->hash = generateHash();
    new_gene->ptr = nullptr;
    new_gene->current_time = 0;

    if (param != -1) {
        new_gene->mutator_connections = uniformTestRange(0, 0.0001, 1000);
       // mutate_connection(new_gene->mutator_connections, genome, *new_gene);
    }

    return new_gene;
}

genes* instruction_gene(dna& genome, int param) {
    genes* new_gene = new genes;
    if (param < 30) {
        new_gene->is_parameter = true;
        new_gene->is_value = false;
        new_gene->is_collection = false;
        new_gene->expression = 0;
        new_gene->synapse_time = uniformTest(1, 9);

    } else if (param < 70) {
        new_gene->is_parameter = false;
        new_gene->is_value = true;
        new_gene->is_collection = false;
        new_gene->expression = 1;
        new_gene->synapse_time = uniformTest(1, 9);

    } else if (param < 95) {
        new_gene->is_parameter = false;
        new_gene->is_value = false;
        new_gene->is_collection = true;
        new_gene->expression = 2;

        new_gene->collection_index = uniformTest(0, 7);

        new_gene->collection_index_mutator = uniformTestRange(0, 0.1, 1000);

        new_gene->synapse_time = uniformTest(1, 10);

    } else {
        new_gene->expression = -1;
    }
    new_gene->is_active = true;
    new_gene->hash = generateHash();
    new_gene->ptr = nullptr;
    new_gene->current_time = 0;


    return new_gene;
}

void mutate_dna(dna& genome) {
    double val = 0;
    randomOutput(val, genome.mutator_add, 1, 0, 1);
    if (val) {
        int index = uniformTest(0, genome.dominant_strand.size() - 1);
        genes* ptr=make_random_gene(genome);
        if (index == -1) {
            index = 0;
        }

        if((index-1)>0){
            ptr->level=genome.dominant_strand.at(index-1).level;
        }
        else{
        ptr->level=0;
        }
        genome.dominant_strand.insert(genome.dominant_strand.begin() + index, *ptr);

        genome.dominant_strand_mutators.insert(genome.dominant_strand_mutators.begin() + index,
                                               uniformTestRange(0, 0.0001, 1000));
    }

    // int size = genome.dominant_strand_mutators.size();
    // for (int i = 0; i < size; i++) {
    //     double val = 0;
    //     randomOutput(val, genome.dominant_strand_mutators.at(i), 1, 0, 1);
    //     if (val && genome.dominant_strand.at(i).is_active == true) {
    //         genome.dominant_strand.at(i).is_active = false;  // deactivate gene if it is turned on
    //     } else {
    //         genome.dominant_strand.at(i).is_active = true;  // activate gene if it is turned off
    //     }
    // }

    val = 0;
    randomOutput(val, genome.mutator_delete, 1, 0, 1);

    if (val) {
        int index = uniformTest(0, genome.dominant_strand.size() - 1);
        if (index <= 0) {
            index = 0;
        } else {
            genome.dominant_strand.erase(genome.dominant_strand.begin() + index);
            genome.dominant_strand_mutators.erase(genome.dominant_strand_mutators.begin() + index);
        }
    }

    for (int i = 0; i < genome.dominant_strand.size(); i++) {
        if (1 ) {
            mutate_connection(genome.dominant_strand.at(i).mutator_connections, genome,
                              genome.dominant_strand.at(i));
        }
        if (genome.dominant_strand.at(i).is_parameter && uniformTest(0, 100) == 5) {
            genome.dominant_strand.at(i).parameter_index =
                genome.dominant_strand.at(i).parameter_index + uniformTestRange(-0.001, 0.001, 1000);
        }
        if (genome.dominant_strand.at(i).is_collection && uniformTest(0, 100) == 5) {
            genome.dominant_strand.at(i).collection_index = uniformTest(0, 7);
        }
    }
    if (genome.dominant_strand.size() == 0) {
        genome.is_dead = true;
        ;
    }
}
dna* make_random_dna() {//you need to fix this
    dna* new_dna = new dna;
    int number_of_nodes = 3;
    for (int i = 0; i < number_of_nodes; i++) {
        genes* ptr=make_random_gene(*new_dna);

        new_dna->dominant_strand.push_back(*ptr);
        new_dna->dominant_strand_mutators.push_back(uniformTestRange(0, 0.01, 1000));
    }
    new_dna->mutator_add = (double)uniformTestRange(0, 0.1, 1000);
    new_dna->mutator_delete = (double)uniformTestRange(0, 0.1, 1000);
    new_dna->is_dead = false;
    return new_dna;
}

dna* instruction_dna(int number_of_nodes, int array[], int size) {
    dna* new_dna = new dna;

    int count_marker=0;
    for (int i = 0; i < number_of_nodes; i++) {
        genes* ptr =instruction_gene(*new_dna, array[i]);
        if(ptr->expression==-1){
            count_marker++;
        }
        ptr->level=count_marker;
        new_dna->dominant_strand.push_back(*ptr);
        new_dna->dominant_strand_mutators.push_back(uniformTestRange(0, 0.001, 1000));
    }
    new_dna->mutator_add = (double)uniformTestRange(0, 0.01, 100);
    new_dna->mutator_delete = (double)uniformTestRange(0, 0.001, 100);
    new_dna->is_dead = false;
    return new_dna;
}

dna* reproduce(dna& mater, dna& matee) {  // come back to this tomorrow morning
    dna* new_dna = new dna;
    new_dna->mutator_add = 0.01;
    new_dna->mutator_delete = 0.0001;
    
    
    new_dna->is_dead = false;
    unordered_map<string, int> map;
    if (mater.dominant_strand.size() > matee.dominant_strand.size()) {
        int m = matee.dominant_strand.size();
        for (int i = 0; i < m; i++) {
            double return_mutator;
            genes* ptr = &randomGene(mater.dominant_strand.at(i), matee.dominant_strand.at(i),
                                     mater.dominant_strand_mutators.at(i), matee.dominant_strand_mutators.at(i),
                                     return_mutator, 0.5);
            auto it = map.find(ptr->hash);
            if (it == map.end()) {
                map[ptr->hash] = 1;
                ptr->current_time = 0;
                new_dna->dominant_strand.push_back(*ptr);
                new_dna->dominant_strand_mutators.push_back(return_mutator);
            }
        }

        int bound = mater.dominant_strand.size();
        for (int i = m; i < bound; i++) {
            auto it = map.find(mater.dominant_strand.at(i).hash);
            if (it == map.end()) {
                map[mater.dominant_strand.at(i).hash] = 1;
                mater.dominant_strand.at(i).current_time = 0;
                new_dna->dominant_strand.push_back(mater.dominant_strand.at(i));
                double return_mutator=mater.dominant_strand_mutators.at(i);
                new_dna->dominant_strand_mutators.push_back(return_mutator);
            }
        }

    } else {
        int m = mater.dominant_strand.size();

        for (int i = 0; i < m; i++) {
            double return_mutator;
            genes* ptr = &randomGene(mater.dominant_strand.at(i), matee.dominant_strand.at(i),
                                     mater.dominant_strand_mutators.at(i), matee.dominant_strand_mutators.at(i),
                                     return_mutator, 0.5);
            auto it = map.find(ptr->hash);
            if (it == map.end()) {
                map[ptr->hash] = 1;
                ptr->current_time = 0;
                new_dna->dominant_strand.push_back(*ptr);
                // if(ptr->abs_time>1000){
                //     if(return_mutator>0){
                //         return_mutator=return_mutator-0.0001;
                //             }
                //     }
                new_dna->dominant_strand_mutators.push_back(return_mutator);
            }
        }

        int bound = matee.dominant_strand.size();
        for (int i = m; i < bound; i++) {
            auto it = map.find(matee.dominant_strand.at(i).hash);
            if (it == map.end()) {
                map[matee.dominant_strand.at(i).hash] = 1;
                matee.dominant_strand.at(i).current_time = 0;
                double return_mutator=matee.dominant_strand_mutators.at(i);
                new_dna->dominant_strand.push_back(matee.dominant_strand.at(i));
                new_dna->dominant_strand_mutators.push_back(return_mutator);
            }
        }
    }

    mutate_dna(*new_dna);
    return new_dna;
}

bool output_now(float probability){
    if(uniformTestRange(0,1,1000)<=probability){
        return true;
    }
    return false;
}
dna* asexual_reproduce(dna& mater) {  // come back to this tomorrow morning
    dna* new_dna = new dna;
    float div=uniformTestRange(-0.0001,0.0001,10);
    if(0<new_dna->mutator_add+div){
    new_dna->mutator_add = mater.mutator_add+div;
    }
    else{
    new_dna->mutator_add = mater.mutator_add;
    }

    float div1=uniformTestRange(-0.0001,0.0001,10);
    if(0<new_dna->mutator_delete+div1){
    new_dna->mutator_delete = mater.mutator_delete+div1;
    }
    else{
    new_dna->mutator_delete = mater.mutator_delete;
    }
    new_dna->is_dead = false;
    unordered_map<string, int> map;
    int m = mater.dominant_strand.size();
    for (int i = 0; i < m; i++) {
        double return_mutator = mater.dominant_strand_mutators.at(i);
        genes* ptr = &mater.dominant_strand.at(i);
        ptr->current_time = 0;
        new_dna->dominant_strand.push_back(*ptr);
        new_dna->dominant_strand_mutators.push_back(return_mutator);
    }
    mutate_dna(*new_dna);
    return new_dna;
}


#endif /* Dna_hpp */