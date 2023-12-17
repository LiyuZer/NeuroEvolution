
//  Organisms.hpp
//  Evolutionary_Neural+Networks
//
//  Created by Liyu Zerihun on 1/7/23.

#ifndef Organisms_hpp
#define Organisms_hpp

#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdlib>
#include "Node.hpp"
#include "Utilities.hpp"
#include "dna.hpp"

using namespace std;


const int connection_maximum = 1000;  // maximum number of connections one node can have
const int synapses_wait_time = 10;    // the number of synapses that will be allowed before food disappears;


// A sims object class, the node structure and food objects will be a sub class of this parent class
class sim_objects {
   public:
    sim_objects(bool validate_organism) { is_organism = validate_organism; }

   private:
    bool is_organism;
};

class organism : public sim_objects {
   public:
    organism(double val, Network<Node>* m) : sim_objects(true) {
        genome = make_random_dna();
        life = 20000;
        time_since_reproduction = -10;
        num_of_synapses = 0;
        score = 0;
        main = m;
    }

    organism(double val, Network<Node>* m, int array[], int size) : sim_objects(true) {
        genome = instruction_dna(size, array, size);
        life = 20000;
        time_since_reproduction = -10;
        num_of_synapses = 0;
        score = 0;
        main = m;
    }

    organism() : sim_objects(true) {
        life = 20000;
        time_since_reproduction = -10;
        num_of_synapses = 0;
        score = 0;
    }

    organism(dna* ptr, Network<Node>* m) : sim_objects(true) {
        genome = (ptr);
        life = 20000;
        time_since_reproduction = -10;
        num_of_synapses = 0;
        score = 0;
        main = m;
    }

    void explore_gene(genes& x1){

        if(x1.node_connections.size()>0){
            string s=x1.node_connections.at(x1.node_connections.size()-1);
        for(int i=0; i<x1.past_forms.size(); i++){
            if(s==x1.past_forms.at(i)){
                return;
            }
        }
        x1.past_forms.push_back(s);
    }

    }
    void express_dna() {
        string s;
        int marker=0;
        for (int i = 0; i < genome->dominant_strand.size(); i++) {
            if (genome->dominant_strand.at(i).expression != -1 &&
                genome->dominant_strand.at(i)
                    .is_active) {  // If the gene is a marked gene add 1 if not then express the gene
                    genome->dominant_strand.at(i).abs_time++;
                    explore_gene(genome->dominant_strand.at(i));
                if (genome->dominant_strand.at(i).is_parameter) {
                    s=s+"p";
                    parameter_node* ptr = new parameter_node(main, true, 0);
                    ptr->getInputVector()->push_back(genome->dominant_strand.at(i).parameter_index);
                    create_nodes.push_back(ptr);  // parameter node is just created with it's value
                    active_genes.push_back(&genome->dominant_strand.at(i));
                    gene_node_map[ptr] = &genome->dominant_strand.at(i);
                    genome->map[genome->dominant_strand.at(i).hash] = ptr;
                    genome->dominant_strand.at(i).ptr = ptr;
                } else if (genome->dominant_strand.at(i).is_value) {
                    s=s+"v";
                    value_node* ptr = new value_node(main, true, 0);
                    ptr->getInputVector()->push_back(0);
                    create_nodes.push_back(ptr);
                    active_genes.push_back(&genome->dominant_strand.at(i));
                    genome->map[genome->dominant_strand.at(i).hash] = ptr;
                    genome->dominant_strand.at(i).ptr = ptr;
                    gene_node_map[ptr] = &genome->dominant_strand.at(i);

                } else {
                    s=s+"c";
                    collection_node* ptr =
                        new collection_node(genome->dominant_strand.at(i).collection_index, main, true, 0);
                    ptr->getInputVector()->push_back(1);
                    create_nodes.push_back(ptr);
                    active_genes.push_back(&genome->dominant_strand.at(i));
                    genome->map[genome->dominant_strand.at(i).hash] = ptr;
                    genome->dominant_strand.at(i).ptr = ptr;
                    gene_node_map[ptr] = &genome->dominant_strand.at(i);
                }
            }

            else {  // add one to marker
                s=s+"m";
                genome->map[genome->dominant_strand.at(i).hash] = nullptr;
                if (genome->dominant_strand.at(i).expression == -1) {
                    marker++;
                    active_genes.push_back(&genome->dominant_strand.at(i));
                }
            }

            // score=score-0.0001;// this decrease the number of nodes;
        }

        for (int i = 0; i < active_genes.size();i++) {
            // This connects the nodes together, if the node is not a nullptr and if the node does not connect
            // to a node that is a nullptr then make the connection
            Node* ptr = active_genes.at(i)->ptr;
            if (ptr != nullptr && active_genes.at(i)->expression!=-1) {
                for (int m = 0; m < active_genes.at(i)->node_connections.size(); m++) {
                    string hash = active_genes.at(i)->node_connections.at(m);
                    auto found = genome->map.find(hash);
                    if (found == genome->map.end()) {
                        active_genes.at(i)->node_connections.erase(
                            active_genes.at(i)->node_connections.begin() + m);
                    } else if (genome->map[hash] != nullptr) {
                        *active_genes.at(i)->ptr >> *(genome->map[hash]);
                    }
                }
            }
        }

    }

    vector<double>* synapse(vector<double>& food, vector<double>& opt, int round) {
        output_num.clear();
        optimal.clear();
        int count = 0;
        int der_number = 0;
        int max = 0;
        optimal = opt;

        if (active_genes.size() > 1) {
            der_number = active_genes.at(0)->current_time;
            max = active_genes.at(0)->current_time + 20;
        }

        vector<genes*> input;
        vector<genes*> output;

        size = 3;
        for (int i = 0; i < size; i++) {
            in.push_back(food.at(i));
        }

        int count_food = 0;  // count for the index of the food vector(which is the actual input for the organism
        bool found_input = false;
        bool found_output = false;
        bool found_hidden = false;
        int count_input = 0;
        for (int i = 0; i < active_genes.size(); i++) {
            if (active_genes.at(i)->expression != -1 && count < 1) {
                if (active_genes.at(i)->is_value) {

                    if (!found_input && active_genes.at(i)->is_active) {
                                found_input = true;
                    }
                    if (count_food < food.size()) {
                        active_genes.at(i)->ptr->getInputVector()->at(0) =
                            active_genes.at(i)->ptr->getInputVector()->at(0) + food.at(count_food);
                        count_food++;
                    }
                    input.push_back(active_genes.at(i));
                } else {
                    input.push_back(active_genes.at(i));
                }
            } else if (count>= 1 && count< genome->division &&active_genes.at(i)->expression != -1) {
                if (!found_hidden && active_genes.at(i)->is_active) {
                    found_hidden = true;
                }
            } else if (active_genes.at(i)->expression != -1 && count >= 2) {
                if (active_genes.at(i)->is_value && active_genes.at(i)->is_active && count==genome->division ) {
                    output.push_back(active_genes.at(i));
                    if (!found_output) {
                        found_output = true;
                    }
                }

            } else if (active_genes.at(i)->expression == -1) {
                count++;
            }
        }

        for (; der_number < max; der_number++) {
            // print();
            for (int i = 0; i < input.size(); i++) {
                cycle(input.at(i)->ptr, der_number);
            }
        }

        vector<double>* return_vec = new vector<double>;
        float op=0;
        // if(round<2000){
        //     op=0.7;
        // }
        // else if(round<4000){
        //     op=0.4;
        // }
        // else{
        //     op=0.3;
        // }
        op=0.2;
        for (int i = 0; i < output.size(); i++) {
            if (i < optimal.size() && output.at(i)->current_time>0) {
                if (input.size() > 0) {
                    double sub =optimal.at(i)- output.at(i)->ptr->getInputVector()->at(0) ;
                    double num = 0;

                   num = exp(-1 * pow(sub, 2) / (2 * pow(op, 2)));
                   // num=pow(sub+0.001,2);
                    //num=1/num;
                    score = (num) + score;
                    output_num.push_back(output.at(i)->ptr->getInputVector()->at(0));
                    out.push_back(output.at(i)->ptr->getInputVector()->at(0));
                }
            }
            return_vec->push_back(output.at(i)->ptr->getInputVector()->at(0));
            output.at(i)->ptr->getInputVector()->at(0) = 0;
        }
        if (!found_output || !found_input || !found_hidden) {
            score = 0;
            score = found_output + found_input + found_hidden;
        } else {
            score = score + found_output + found_input + found_hidden;
        }

        return return_vec;
    }

    void pr() {
        cout << "\033[36m ***This is the best organism, good job, you have won the matrix***" << endl;
        cout << "Size of input: " << size << endl;

        cout << score << endl;// prints out  the score of the organism
        int sz = 0;
        if (out.size() < in.size()) {
            sz = out.size();
        } else {
            sz = in.size();
        }
        for (int i = 0; i < sz; i++) {
            cout << "\033[33mInput: " << in.at(i) << " Output: " << out.at(i) << endl;
        }
        cout << "\033[30m****" << endl;
    }

    void cycle(Node* ptr, int num) {
        // breadth first search algorithm
        vector<genes*> vec;
        queue<genes*> list;
        list.push(gene_node_map[ptr]);
        while (!list.empty()) {
            genes* ptr_1 = list.front();
            list.pop();
            ptr = ptr_1->ptr;
            int size = ptr_1->ptr->getNextVector()->size();
            genes* gene_ptr = ptr_1;
            if (gene_ptr->current_time <= num) {
                if (gene_ptr->current_time>=gene_node_map[ptr]->synapse_time){
                    if (ptr->getNextVector()->size() > 0) {
    

                        ptr->special_activation(gene_node_map[ptr]->gene_param,gene_node_map[ptr]->bias );
                    }
                    gene_ptr->current_time++;
                }
                else {
                    gene_ptr->current_time++;
                }
                for (int m = 0; m < size; m++) {
                    list.push(gene_node_map[ptr_1->ptr->getNextVector()->at(m)]);
                }
            }
        }

    }

    void print_node(Node* ptr, string hash) {
        if (ptr != nullptr && ptr->get_is_parameter()) {
            cout << "This is a parameter node: " << hash << endl;

        } else if (ptr != nullptr && ptr->get_is_value_node()) {
            cout << "This is a value node: " << hash << endl;

        } else if (ptr != nullptr) {
            cout << "This is a collection node: " << hash << endl;
        }
    }
    void print_connections(genes& m) {
        cout << "\033[36m\033[1mConnection node*****************" << endl;
        for (int i = 0; i < m.node_connections.size(); i++) {
            print_node(genome->map[m.node_connections.at(i)], m.node_connections.at(i));
        }
        cout << "*************************" << endl;
        cout << endl;
    }
    void print() {
        cout << "\033[33m\033[1mThe size of the organism is:" << genome->dominant_strand.size() << endl;
        //The extra characters at the beginning are the color markers for linux and mac terminals 
        for (int i = 0; i < genome->dominant_strand.size(); i++) {
            cout << "\033[33m\033The hash of this node is: " << genome->dominant_strand.at(i).hash << endl;
            cout << "The synapse wait time is " << genome->dominant_strand.at(i).synapse_time << endl;
            cout << "The absolute time of this node is: " << genome->dominant_strand.at(i).abs_time << endl;
            cout << "The level: " << genome->dominant_strand.at(i).level << endl;
            cout << "The parameter size; " <<genome->dominant_strand.at(i).gene_param.size();

            cout << "The parameter: ";
            for(int k=0; k<genome->dominant_strand.at(i).gene_param.size(); k++){
                cout<<genome->dominant_strand.at(i).gene_param.at(k)<<" ";
            }
            cout<<endl;

            if(past_forms.size()>0){
            cout << "The nearest past_form: " << past_forms.at(0)<< endl;
            }

            if (genome->dominant_strand.at(i).expression != -1 && genome->dominant_strand.at(i).is_parameter) {
                cout << "\033[33m\033[1mThis is a parameter node" << endl;
                cout << "Active status: " << genome->dominant_strand.at(i).is_active << endl;
                if (genome->dominant_strand.at(i).is_active) {
                    cout << "The synapse number is: " << genome->dominant_strand.at(i).current_time << endl;
                    cout << "First input: " << genome->dominant_strand.at(i).ptr->getInputVector()->at(0) << endl;
                }

                print_connections(genome->dominant_strand.at(i));
            } else if (genome->dominant_strand.at(i).expression != -1 && genome->dominant_strand.at(i).is_value) {
                cout << "\033[33m\033[1mThis is a value node" << endl;
                cout << "Active status: " << genome->dominant_strand.at(i).is_active << endl;
                if (genome->dominant_strand.at(i).is_active) {
                    cout << "The synapse number is: " << genome->dominant_strand.at(i).current_time << endl;
                    cout << "First input: " << genome->dominant_strand.at(i).ptr->getInputVector()->at(0) << endl;
                }

                print_connections(genome->dominant_strand.at(i));
            } else if (genome->dominant_strand.at(i).expression != -1) {
                cout << "\033[33m\033[1mThis is a collection node" << endl;
                cout << "Active status: " << genome->dominant_strand.at(i).is_active << endl;
                cout << "Collection Index: " << genome->dominant_strand.at(i).collection_index << endl;
                cout << "The collection  time is " << genome->dominant_strand.at(i).collective_time << endl;
                if (genome->dominant_strand.at(i).is_active) {
                    cout << "The synapse number is: " << genome->dominant_strand.at(i).current_time << endl;
                    cout << "Inputs:";
                    if (genome->dominant_strand.at(i).ptr->getInputVector() != nullptr) {
                        for (int m = 0; m < genome->dominant_strand.at(i).ptr->getInputVector()->size(); m++) {
                            cout << genome->dominant_strand.at(i).ptr->getInputVector()->at(m) << ", ";
                        }
                        cout << endl;
                    }
                }

                print_connections(genome->dominant_strand.at(i));
            } else if (genome->dominant_strand.at(i).expression == -1) {
                cout << "This is a mutator node" << endl;
                cout << "Active status: " << genome->dominant_strand.at(i).is_active << endl;
                cout << endl;
            }
        }
    }
    string return_print() {
        string s;
        s = s + "The size of the organism is: " + to_string(create_nodes.size());
        for (int i = 0; i < create_nodes.size(); i++) {
            s = s + "The hash of this node is: " + genome->dominant_strand.at(i).hash;
            if (genome->dominant_strand.at(i).expression != -1 && genome->dominant_strand.at(i).is_parameter) {
                s = s + "This is a parameter node";
            } else if (genome->dominant_strand.at(i).expression != -1 && genome->dominant_strand.at(i).is_value) {
                s = s + "This is a value node";
            } else if (genome->dominant_strand.at(i).expression != -1) {
                s = s + "This is a collection node";
            }
        }

        return s;
    }

    void close() {  // deletes the dynamically appointed objects after creation
        for (int i = 0; i < create_nodes.size(); i++) {
            delete create_nodes.at(i);
        }
        genome->map.clear();
        create_nodes.clear();
    }

    bool status() { return genome->is_dead; }
    dna* getGenome() { return genome; }
    void set_genome(dna* genes) { genome = genes; }
    void decrease_life() {
        (life) = life - 1;
        time_since_reproduction++;
    }
    double get_life() { return life; }
    float get_score() { return score; }
    double get_time_since_reproduction() { return time_since_reproduction; }
    void set_time_since_reproduction(double d) { time_since_reproduction = d; }
    vector<genes*>* get_active_genes() { return &active_genes; }

    ~organism() { delete genome; }
    void clear() {
        create_nodes.clear();
        active_genes.clear();
        gene_node_map.clear();
        optimal.clear();
        output_num.clear();
        score = 0;
        genome->map.clear();
    }

   private:
    dna* genome;
    vector<Node*> create_nodes;
    vector<genes*> active_genes;
    unordered_map<Node*, genes*> gene_node_map;
    vector<double> optimal;
    vector<double> output_num;
    vector<double> in;
    vector<double> out;
    double life;
    double time_since_reproduction;
    int num_of_synapses;
    double score;
    Network<Node>* main;
    vector<string> past_forms;
    int size;
};

bool can_reproduce(organism* a, organism* b, double max_life) {
    double life_amountA = a->get_life() / max_life;
    double life_amountB = b->get_life() / max_life;
    if (life_amountA >= 0.3 && life_amountB >= 0.3 && a->get_time_since_reproduction() > 15 &&
        b->get_time_since_reproduction() > 15) {
        a->set_time_since_reproduction(-16);
        b->set_time_since_reproduction(-16);

        return true;
    } else {
        return false;
    }
}

bool compare(organism*& a, organism*& b) {
    if (a->get_score() > b->get_score()) {
        return true;
    }
    return false;
}


class food : public sim_objects {
   public:
    food(int num) : sim_objects(false) {
        int count = 0;
        for (int i = 0; i < num; i++) {
            int val = uniformTest(0, 1);
            food_values.push_back(val);
            if (val != 0) {
                count++;
            }
        }

        output_wanted.push_back(count);
    }

   private:
    vector<double> food_values;
    vector<double> output_wanted;
};



static void run_world(int population_sizes[],int num_pop, int max_round, int no_inputs, int no_output, float percentage_size ){
    vector<organism*> creatures;
    Network<Node>* main = new Network<Node>(1);

    vector<double> food;
    vector<double> optimal;

    int array[20];

    for(int i=0; i<20; i++){
        array[i]=55;
    }
    array[4]=99;
    array[15]=99;
 







    for (int i = 0; i < population_sizes[0]; i++) {
        organism* ptr = new organism(100, main, array, 20);
        creatures.push_back(ptr);
        //mutate_dna(*creatures.at(i)->getGenome());
    }

    int div=(float)max_round/num_pop;
    float array_num[40];
    for (size_t i = 0; i < 40; i++)
    {
        array_num[i]=uniformTestRange(0,1,100);
    }
    
    for(int round = 0; round < num_pop; round++) {

        for(int sub_rounds=0; sub_rounds<div; sub_rounds++){
                for (int k = 0; k < 10; k++) {
                    food.clear();
                    optimal.clear();
                    for (int i = 0; i < 1; i++) {

                        int count_num=0;
                        // for(int m=0; m<3; m++){
                        //     double n=uniformTestRange(0.1,1,1);
                        //     if(n<0.5){
                        //     }
                        //     else{
                        //         count_num++;
                        //     }
                        //     food.push_back(n);

                        // }

                        //     for(int m=0; m<3; m++){
                        //     if(m==count_num-1){
                        //         optimal.push_back(1);
                        //     }
                        //     else{
                        //         optimal.push_back(0.1);
                        //     }
                        // }
                        for(int m=0; m<3; m++){
                            double n=uniformTestRange(0,1,1000);
                            food.push_back(n);
                            optimal.push_back(n*2);


                        }


                        for (int i = 0; i < population_sizes[round]; i++) {
                        organism* ptr = creatures.at(i);
                        if (k == 0) {
                            ptr->express_dna();
                        }
                            ptr->synapse( food, optimal, sub_rounds);
                    }
                }
                }
                    sort(creatures.begin(), creatures.end(), compare);
                    if ((sub_rounds + 1) % 10 == 0) {
                        system("clear");
                        creatures.at(0)->pr();
                       //cout<< creatures.at(0)->getGenome()->dominant_strand.size()<<endl;
                        creatures.at(0)->print();
                    }

                        for (int i = 0; i < population_sizes[round]; i++) {
                                int choosen_one=uniformTest(0,population_sizes[round]*percentage_size);
                                int choosen_two=uniformTest(0,population_sizes[round]*percentage_size);

                               // organism* child = new organism(asexual_reproduce(*creatures.at(choosen_one)->getGenome()), main);
                               organism* child = new organism(reproduce(*creatures.at(choosen_one)->getGenome(), *creatures.at(choosen_one)->getGenome()), main);
                                creatures.push_back(child);
                        }
                            for (int i = 0; i < population_sizes[round]; i++) {
                                delete creatures.at(0);
                                creatures.erase(creatures.begin());
                            }
                }
                main->clear(0);
    }
}


#endif /* Organisms_hpp */
