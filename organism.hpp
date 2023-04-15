
//  Organisms.hpp
//  Evolutionary_Neural+Networks
//
//  Created by Liyu Zerihun on 1/7/23.


#ifndef Organisms_hpp
#define Organisms_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <list>
#include <cmath>
#include <unordered_map>
#include <random>
#include <sstream>
#include <climits>
#include "Node.hpp"
#include "Matrix.hpp"
#include "Utilities.hpp"
#include "Network.hpp"

#include <iostream>
#include <functional>
#include <chrono>
#include <thread>
#include <queue>





//the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define ANSI_RED     "\x1b[31m"
#define ANSI_GREEN   "\x1b[32m"
#define ANSI_YELLOW  "\x1b[33m"
#define ANSI_BLUE    "\x1b[34m"
#define ANSI_MAGENTA "\x1b[35m"
#define ANSI_CYAN    "\x1b[36m"
#define ANSI_BLACK   "\x1b[0m"
#include <chrono>

const int connection_maximum=1000;// maximum number of connections one node can have
const int synapses_wait_time=10; // the number of synapses that will be allopwed before food dissapears;
using namespace std;
int count_o=0;

//if the gene is -1 then it is a marker gene expressing that a certain interval has neded
 struct genes{
    bool is_parameter;
    bool is_collection;
    bool is_value;
    bool is_active;
    bool collection_index_mutator;
    int collection_index;
    int expression;
    int synapse_time;
    int collective_time;
    int current_time=0;/// current synapse time, versus max synapse time
     int abs_time=0;
    double threshold;
    string type;
    string hash;
    Node* ptr;
    double mutator_connections;
    vector<string> node_connections;
    std::unordered_map<std::string, int> node_connections_found;
     float parameter_index;
    genes& operator=(const genes& right){
        node_connections.assign(right.node_connections.begin(), right.node_connections.end());
        node_connections_found=right.node_connections_found;
        ptr=nullptr;
        is_parameter=right.is_parameter;
        is_collection=right.is_collection;
        is_value=right.is_value;
        is_active=right.is_active;
        collection_index_mutator=right.collection_index_mutator;
        collection_index=right.collection_index;

        expression=right.expression;
        type=right.type;
        hash=right.hash;
        expression=right.expression;
        threshold=right.threshold;
        mutator_connections=right.mutator_connections;
        synapse_time=right.synapse_time;
        current_time=right.current_time;
        parameter_index=right.parameter_index;
        collective_time=right.collective_time;
        abs_time=right.abs_time;

        return *this;
     }

     
    };

  struct dna{// dominant strand of the gene is the one that will be expressed when passed down
    vector<genes> dominant_strand;
   // vector<genes> silent_strand; do later
    vector<double> dominant_strand_mutators;
     std::unordered_map<std::string, Node*> map;
    //vector<double> silent_strand_mutators; do later
      double mutator_add;
      double mutator_delete;
      bool is_dead;
      
      
      dna& operator=(const dna& right){
          dominant_strand.clear();
          dominant_strand_mutators.clear();
          if(right.dominant_strand.size()>0){
          dominant_strand.assign(right.dominant_strand.begin(), right.dominant_strand.end());
          dominant_strand_mutators.assign(right.dominant_strand_mutators.begin(), right.dominant_strand_mutators.end());
          }
          mutator_add=right.mutator_add;
          mutator_delete=right.mutator_delete;
          is_dead=right.is_dead;

          map.clear();
          
          return *this;

      }
      



};


static string generateHash(){
    
    hash<long long> hasher;  // define the hash function
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<long long> dist(0, LLONG_MAX);
        long long random_number = dist(mt);
        long long random_number1 = dist(mt);

//        std::hash<long long> hasher;
//        std::size_t hash = hasher(random_number);
//        std::stringstream ss;
//        ss << hash;

    count_o++;
    return to_string(count_o);
}

static void randomOutput(double& val, double probablility_of_change, double step_of_change, double lower_bound, double higher_bound){

    double division=log10(probablility_of_change)-2;// this i the division for the random number
    double random=uniformTest(0, 1, pow(10,-1*division));// our radnom number
    double random_num=0.5;// just some arbitrary value that will be checked
    if(random_num-probablility_of_change<=random &&  random<=random_num+probablility_of_change){// if the guess is within the range then the addition will be added

        int direction=uniformTest(-1, 1,1);
        double temp_val=val;
        temp_val=temp_val+direction*step_of_change;

        if(temp_val<lower_bound|| temp_val>higher_bound){

        }
        else{
            val=temp_val;
        }
    }
}


static genes& randomGene(genes& first_gene, genes& second_gene, double first_mutator, double second_mutator, double& return_value, double probablility_of_change){
    genes* returnGene;
    double division=log10(probablility_of_change)-2;// this i the division for the random number
    double random=uniformTest(0, 1, 1);// our radnom number
    double random_num=0.5;// just some arbitrary value that will be checked
    if(random==0){// if the guess is within the range then the addition will be added
        return_value=first_mutator;
        return first_gene;
    }
    else{
        return_value=second_mutator;
        return second_gene;
    }
}


static void mutate_connection(double mutator, vector<genes>& strand_of_genes, genes& x1){
    double value=0;
    value=uniformTest(0,10,1);
    if(value<5){
        int erase=uniformTest(0, 10, 1);
        if(erase>9 &&  x1.node_connections.size()>0 ){
            double index=uniformTest(0, x1.node_connections.size()-1, 1);
            x1.node_connections_found.erase(x1.node_connections.at(index));
            x1.node_connections.erase(x1.node_connections.begin()+index);
        }
        else if(strand_of_genes.size()>0){
            double index=uniformTest(0, strand_of_genes.size()-1, 1);
            auto it = x1.node_connections_found.find(strand_of_genes.at(index).hash);
            if (it == x1.node_connections_found.end()&&  strand_of_genes.at(index).expression!=-1) {
                x1.node_connections.push_back(strand_of_genes.at(index).hash);
                x1.node_connections_found[strand_of_genes.at(index).hash]=1;
        }
        
    }
    }
}
    
static genes* make_random_gene(dna& genome){
    genes* new_gene=new genes;
    double param=uniformTest(0, 100, 1);
    if(param<10){
        new_gene->is_parameter=true;
        new_gene->is_value=false;
        new_gene->is_collection=false;
        new_gene->expression=0;
        new_gene->parameter_index=uniformTest(-1,1,1000);
        new_gene->synapse_time=uniformTest(1,3,1);

    }
    else if(param<70){
            new_gene->is_parameter=false;
            new_gene->is_value=true;
            new_gene->is_collection=false;
            new_gene->expression=1;
            new_gene->synapse_time=uniformTest(1,3,1);

        }
    else if(param<90){
            new_gene->is_parameter=false;
            new_gene->is_value=false;
            new_gene->is_collection=true;
            new_gene->expression=2;
            new_gene->synapse_time=uniformTest(1,10,1);
            new_gene->collection_index=uniformTest(0, 2, 1);
            new_gene->collection_index_mutator=uniformTest(0, 0.1, 1000);
            new_gene->collective_time=2;
            
        }
    else{
        new_gene->expression=-1;
    }
    new_gene->is_active=true;
    new_gene->hash=generateHash();
    new_gene->ptr=nullptr;
    new_gene->current_time=0;

    if(param!=-1){
        new_gene->mutator_connections=uniformTest(0, 0.1, 1000);
        mutate_connection(new_gene->mutator_connections, genome.dominant_strand, *new_gene);
    }
    
    return new_gene;
}
     static genes* instruction_gene(dna& genome, int param){
         genes* new_gene=new genes;
         if(param<30){
             new_gene->is_parameter=true;
             new_gene->is_value=false;
             new_gene->is_collection=false;
             new_gene->expression=0;
             new_gene->synapse_time=uniformTest(1,3,1);


         }
         else if(param<70){
                 new_gene->is_parameter=false;
                 new_gene->is_value=true;
                 new_gene->is_collection=false;
             new_gene->expression=1;
             new_gene->synapse_time=uniformTest(1,3,1);


             }
         else if(param<95){
                 new_gene->is_parameter=false;
                 new_gene->is_value=false;
                 new_gene->is_collection=true;
             new_gene->expression=2;

                 bool collection_index_mutator;
                 bool collection_index;
                 
                 new_gene->collection_index=uniformTest(0,7,1);

                 new_gene->collection_index_mutator=uniformTest(0, 0.1, 1000);

                 new_gene->synapse_time=uniformTest(1,10,1);

             }
         else{
             new_gene->expression=-1;
         }
         new_gene->is_active=true;
         new_gene->hash=generateHash();
         new_gene->ptr=nullptr;
         new_gene->current_time=0;

         if(param!=-1){
             new_gene->mutator_connections=uniformTest(0, 0.3, 1000);
             for(int i=0; i<10; i++){
             mutate_connection(new_gene->mutator_connections, genome.dominant_strand, *new_gene);
             }
         }
         
         return new_gene;
     }


static void mutate_dna(dna& genome){
    double val=0;
    randomOutput(val, genome.mutator_add, 1, 0, 1);

    if(val){
        int index=uniformTest(0,  genome.dominant_strand.size()-1, 1);
       // int index_first=uniformTest(0,  5, 1);
        if(index==-1){
            index=0;
        }
      //  if(index_first==5){
            
        genome.dominant_strand.insert( genome.dominant_strand.begin()+index,*make_random_gene(genome));
        
        genome.dominant_strand_mutators.insert( genome.dominant_strand_mutators.begin()+index,uniformTest(0, 0.1, 1000));
      //  }
//        else{
//            genome.dominant_strand.insert( genome.dominant_strand.begin(),*make_random_gene(genome));
//            genome.dominant_strand_mutators.insert( genome.dominant_strand_mutators.begin(),uniformTest(0, 0.1, 1000));
//        }
    }
   
        int size=genome.dominant_strand_mutators.size();
        for(int i=0; i<size; i++){
            double val=0;
            randomOutput(val, genome.dominant_strand_mutators.at(i), 1, 0, 1);
            if(val && genome.dominant_strand.at(i).is_active==true){
                genome.dominant_strand.at(i).is_active=false;//deactivate gene if it is turned on
            }
            else{
                genome.dominant_strand.at(i).is_active=true;//activate gene if it is turned off
            }
        }
    
    val=0;
    randomOutput(val, genome.mutator_delete, 1, 0, 1);

    if(val){
        int index=uniformTest(0,  genome.dominant_strand.size()-1, 1);
        if(index<=0){
            index=0;
        }
        else{
            
//            if(genome.dominant_strand.at(index).abs_time<10000){// if the gene is not old then you can erase it
                genome.dominant_strand.erase( genome.dominant_strand.begin()+index);
                genome.dominant_strand_mutators.erase( genome.dominant_strand_mutators.begin()+index);
//            }

        }
        
        
        
    }
        
    
    
    for(int i=0; i<genome.dominant_strand.size(); i++){
        
        if(uniformTest(0,100,1)==5){
        mutate_connection(genome.dominant_strand.at(i).mutator_connections, genome.dominant_strand, genome.dominant_strand.at(i));
        }
        if(genome.dominant_strand.at(i).is_parameter&& uniformTest(0,100,1)==5){
            genome.dominant_strand.at(i).parameter_index=genome.dominant_strand.at(i).parameter_index+uniformTest(-0.001,0.001,1000);
        }
        if(genome.dominant_strand.at(i).is_collection&& uniformTest(0,100,1)==5){
            genome.dominant_strand.at(i).collection_index=uniformTest(0,7,1);
        }
        
    }
    if(genome.dominant_strand.size()==0){
        genome.is_dead=true;;
    }
        
}
static dna* make_random_dna(){
    dna* new_dna=new dna;
    int number_of_nodes=3;
    for(int i=0; i<number_of_nodes; i++){
        new_dna->dominant_strand.push_back(*make_random_gene(*new_dna));
        new_dna->dominant_strand_mutators.push_back(uniformTest(0, 0.01, 1000));
    }
    new_dna->mutator_add=(double)uniformTest(0, 0.01, 1000);
    new_dna->mutator_delete=(double)uniformTest(0, 0.01, 1000);
    new_dna->is_dead=false;
    return new_dna;
}
     
     static dna* instruction_dna(int number_of_nodes, int array[], int size){
         dna* new_dna=new dna;
         
         for(int i=0; i<number_of_nodes; i++){
             new_dna->dominant_strand.push_back(*instruction_gene(*new_dna, array[i]));
             new_dna->dominant_strand_mutators.push_back(uniformTest(0, 0.01, 1000));
         }
         new_dna->mutator_add=(double)uniformTest(0, 0.01, 1000);
         new_dna->mutator_delete=(double)uniformTest(0, 0.01, 1000);
         new_dna->is_dead=false;
         return new_dna;
     }
     





dna* reproduce(dna& mater, dna& matee){// come back to this tomorrow morning
    dna* new_dna= new dna;
    new_dna->mutator_add=(double)uniformTest(0, 0.4, 1000);
    new_dna->mutator_delete=(double)uniformTest(0, 0, 1000);
    new_dna->is_dead=false;
    unordered_map<string, int> map;
    if(mater.dominant_strand.size()>matee.dominant_strand.size()){
        int m=matee.dominant_strand.size();
        for(int i=0; i<m; i++){
            double return_mutator;
            genes* ptr=&randomGene(mater.dominant_strand.at(i), matee.dominant_strand.at(i), mater.dominant_strand_mutators.at(i), matee.dominant_strand_mutators.at(i), return_mutator, 0.5);
            auto it = map.find(ptr->hash);
            if (it == map.end()){
                map[ptr->hash]=1;
                ptr->current_time=0;
                new_dna->dominant_strand.push_back(*ptr);
                new_dna->dominant_strand_mutators.push_back(return_mutator);
                
            }

            
        }

        int bound=mater.dominant_strand.size();
        for(int i=m; i<bound; i++){
            
            
            auto it = map.find(mater.dominant_strand.at(i).hash);
            if (it == map.end()){
                map[mater.dominant_strand.at(i).hash]=1;
                mater.dominant_strand.at(i).current_time=0;
                new_dna->dominant_strand.push_back(mater.dominant_strand.at(i));
                    new_dna->dominant_strand_mutators.push_back(mater.dominant_strand_mutators.at(i));
            }
                                               }

        }
    else{
        int m=mater.dominant_strand.size();
   

        for(int i=0; i<m; i++){
            double return_mutator;
            genes* ptr=&randomGene(mater.dominant_strand.at(i), matee.dominant_strand.at(i), mater.dominant_strand_mutators.at(i), matee.dominant_strand_mutators.at(i), return_mutator, 0.5);
            auto it = map.find(ptr->hash);
            if (it == map.end()){
                map[ptr->hash]=1;
                ptr->current_time=0;
                new_dna->dominant_strand.push_back(*ptr);
                new_dna->dominant_strand_mutators.push_back(return_mutator);
            }
            }

        int bound=matee.dominant_strand.size();
        for(int i=m; i<bound; i++){

            auto it = map.find(matee.dominant_strand.at(i).hash);
            if (it == map.end()){
                map[matee.dominant_strand.at(i).hash]=1;
                matee.dominant_strand.at(i).current_time=0;
                new_dna->dominant_strand.push_back(matee.dominant_strand.at(i));
                new_dna->dominant_strand_mutators.push_back(matee.dominant_strand_mutators.at(i));
            }


            }
        
        

}
    mutate_dna(*new_dna);
    return new_dna;
}


dna* asexual_reproduce(dna& mater){// come back to this tomorrow morning
    dna* new_dna= new dna;
    new_dna->mutator_add=(double)uniformTest(0, 0.1, 1000);
    new_dna->mutator_delete=(double)uniformTest(0, 0.1, 1000);
    new_dna->is_dead=false;
    unordered_map<string, int> map;
        int m=mater.dominant_strand.size();
        for(int i=0; i<m; i++){
            double return_mutator=mater.dominant_strand_mutators.at(i);
            genes* ptr=&mater.dominant_strand.at(i);
                ptr->current_time=0;
//            if(ptr->abs_time>100000){
//                if(return_mutator>0){
//                return_mutator=return_mutator-0.0001;
//                }
//            }
                new_dna->dominant_strand.push_back(*ptr);
                new_dna->dominant_strand_mutators.push_back(return_mutator);
            }
    mutate_dna(*new_dna);
    return new_dna;
}


// A sims object class, the node structure and food objects will be a sub class of this parent class
class sim_objects{
public:

    sim_objects(bool validate_organism){
        is_organism=validate_organism;
    }


private:
    bool is_organism;

};


class organism: public sim_objects{
public:
  organism(double val, Network<Node>* m):
    sim_objects(true){
        genome=make_random_dna();
        life=20000;
        time_since_reproduction=-10;
        num_of_synapses=0;
        score=0;
        main=m;
    }
    
    organism(double val, Network<Node>* m, int array[], int size):
      sim_objects(true){
          genome=instruction_dna(size,array,size);
          life=20000;
          time_since_reproduction=-10;
          num_of_synapses=0;
          score=0;
          main=m;
      }

    organism():
      sim_objects(true){
          life=20000;
          time_since_reproduction=-10;
          num_of_synapses=0;
          score=0;


      }
    
    organism(dna* ptr, Network<Node>* m):
    sim_objects(true){
        genome=(ptr);
        life=20000;
        time_since_reproduction=-10;
        num_of_synapses=0;
        score=0;
        main=m;




    }
    void express_dna(){

        int marker_number=0;// the number of times there has been a marker gene
        for(int i=0; i<genome->dominant_strand.size(); i++){
            if(genome->dominant_strand.at(i).expression!=-1 && genome->dominant_strand.at(i).is_active){//If the gene is a marked gene add 1 if not then express the gene
                if(genome->dominant_strand.at(i).is_parameter){
                    parameter_node* ptr= new parameter_node( main,true,0);
                    ptr->getInputVector()->push_back(genome->dominant_strand.at(i).parameter_index);
                        create_nodes.push_back(ptr);// paramter node is just created with it's value
                    active_genes.push_back(&genome->dominant_strand.at(i));
                    gene_node_map[ptr]=&genome->dominant_strand.at(i);
                    genome->map[genome->dominant_strand.at(i).hash]=ptr;
                    genome->dominant_strand.at(i).ptr=ptr;
                    }
                else if(genome->dominant_strand.at(i).is_value){
                        value_node* ptr= new value_node(main, true,0);
                    ptr->getInputVector()->push_back(0);
                        create_nodes.push_back(ptr);
                    active_genes.push_back(&genome->dominant_strand.at(i));
                    genome->map[genome->dominant_strand.at(i).hash]=ptr;
                    genome->dominant_strand.at(i).ptr=ptr;
                    gene_node_map[ptr]=&genome->dominant_strand.at(i);


                    }
                    else{
                        collection_node* ptr=new collection_node(genome->dominant_strand.at(i).collection_index,main,true,0);
                        ptr->getInputVector()->push_back(1);
                        create_nodes.push_back(ptr);
                        active_genes.push_back(&genome->dominant_strand.at(i));
                        genome->map[genome->dominant_strand.at(i).hash]=ptr;
                        genome->dominant_strand.at(i).ptr=ptr;
                        gene_node_map[ptr]=&genome->dominant_strand.at(i);

                    }
                }

            else{// add one to marker
                marker_number++;
                genome->map[genome->dominant_strand.at(i).hash]=nullptr;
                if(genome->dominant_strand.at(i).expression==-1){
                    active_genes.push_back(&genome->dominant_strand.at(i));
                }

            }
            
            //score=score-0.0001;// this decrease the number of nodes;

        }



        for(int i=0; i<create_nodes.size(); i++){// This connecte the nodes together, if the node is not a nullptr and if the node does not connect to a node that is a nullptr then make the connection
            Node* ptr=create_nodes.at(i);
            if(ptr!=nullptr){
                for(int m=0;m<genome->dominant_strand.at(i).node_connections.size(); m++){
                    string hash=genome->dominant_strand.at(i).node_connections.at(m);
                    auto found=genome->map.find(hash);
                    if(found==genome->map.end()){
                        genome->dominant_strand.at(i).node_connections.erase(genome->dominant_strand.at(i).node_connections.begin()+m);
                    }
                    else if( genome->map[hash]!=nullptr){
                        *create_nodes.at(i)>>*(genome->map[hash]);
                    }
                }
            }

        }

    }
    
    

    vector<double>* synapse(int index, vector<double>& food, vector<double>& opt){
        output_num.clear();
        optimal.clear();
        bool found=false;
        int count=0;
        int der_number=0;
        int max=0;
        optimal=opt;

        if(active_genes.size()>1){
            der_number=active_genes.at(0)->current_time;
            max=active_genes.at(0)->current_time+10;
        }
        
        vector<genes*> input;
        vector<genes*> output;
            
         size=2;
        for(int i=0; i<size; i++){
            in.push_back(food.at(i));

            
        }


        int count_food=0;// count for the index of the food vector(which is the actual input for the organism
        bool found_input=false;
        bool found_output=false;
        bool found_hidden=false;
        float input_num=0;
            for(int i=0; i<active_genes.size(); i++){
                active_genes.at(i)->abs_time++;
                if(active_genes.at(i)->expression!=-1 && count<1){
                     if(active_genes.at(i)->is_value){

                         if(!found_input&&active_genes.at(i)->is_active && active_genes.at(i)->node_connections.size()>0){
                             found_input=true;
                         }
                         input_num++;
                         if(count_food<food.size()){
                        active_genes.at(i)->ptr->getInputVector()->at(0)=active_genes.at(i)->ptr->getInputVector()->at(0)+food.at(count_food);
                             count_food++;
                         }
                        input.push_back(active_genes.at(i));
                    }
                     else{
                         input.push_back(active_genes.at(i));
                     }
                }
                else if(count<=1 && active_genes.at(i)->expression!=-1){
                    if(!found_hidden&&active_genes.at(i)->is_active){
                        found_hidden=true;
                    }
                }
                else if(active_genes.at(i)->expression!=-1 && count>=2){
                    if(active_genes.at(i)->is_value&&active_genes.at(i)->is_active){
                        output.push_back(active_genes.at(i));
                        if(!found_output){
                            found_output=true;
                        }

                   }

  
                }
                else if (active_genes.at(i)->expression==-1){
                    count++;
                }
                
        }


        
        for(; der_number<max; der_number++){
            //print();
            for(int i=0; i<input.size(); i++){
                cycle(input.at(i)->ptr, der_number);
            }
            
        }
        
        vector<double>* return_vec= new vector<double>;

            for(int i=0; i<output.size(); i++){
                if(i<optimal.size()){
                    if(input.size()>0){
                    double sub=output.at(i)->ptr->getInputVector()->at(0)-optimal.at(i);
///                    double num=(double)(1)/(pow(sub,2)+1);
                        double num=exp(-1*pow(sub,2)/(2*pow(0.1,2)));
                    //cout<<num<<endl;
//                        if(sub==0){
//                            num=1;
//                        }
//                        else{
//                            num=-1;
//                        }
//
                    score=num+score;
                        output_num.push_back(output.at(i)->ptr->getInputVector()->at(0));
                        out.push_back(output.at(i)->ptr->getInputVector()->at(0));
                    }
                }
                return_vec->push_back(output.at(i)->ptr->getInputVector()->at(0));
                output.at(i)->ptr->getInputVector()->at(0)=0;
            }
        if(!found_output || !found_input || !found_hidden){
            score=0;
            score=found_output+found_input+found_hidden;
        }
        else{
            score=score+found_output+found_input+found_hidden;

        }
            
        return return_vec;

    }
    
    
    void pr(){
        cout<<"***This is the best organism, good job, you have won the matrix***"<<endl;
        cout<<"Size of input: "<<size<<endl;

        cout<<score<<endl;
        
        for(int i=0; i<optimal.size(); i++){
        cout<<"Optimal single: "<<optimal.at(i)<<endl;
        }

        for(int i=0; i<output_num.size(); i++){
            cout<<"Out single: "<<output_num.at(i)<<endl;
        }
        
        int sz=0;
        if(out.size()<in.size()){
            sz=out.size();
        }
        else{
            sz=in.size();
        }
        for(int i=0; i<sz; i++){
            cout<<"Input: "<<in.at(i)<<" Output: "<<out.at(i)<<endl;
        }
        cout<<"****"<<endl;

    }
    
    void cycle(Node* ptr, int num){
        
        
        
        
        //breadth first search algorithm
        vector<genes*> vec;
        
        queue<genes*> list;
        list.push(gene_node_map[ptr]);
        int max=list.size();
        while(!list.empty()){
            genes* ptr_1=list.front();
            list.pop();
            ptr=ptr_1->ptr;
            int size=ptr_1->ptr->getNextVector()->size();
            genes* gene_ptr=ptr_1;
            if(gene_ptr->current_time<=num){
                if(gene_ptr->current_time!=0&& (gene_ptr->current_time%gene_ptr->synapse_time==0  || (gene_ptr->is_collection&& gene_ptr->collective_time<=gene_ptr->ptr->getInputVector()->size()) )){
                    if(ptr->getNextVector()->size()>0){
                        ptr->special_activation();
                    }
                    gene_ptr->current_time++;
                }
                else if(gene_ptr->current_time==0){
                    if(ptr->getNextVector()->size()>0){
                        ptr->special_activation();
                    }
                    gene_ptr->current_time++;
                }
                else{
                    gene_ptr->current_time++;

                }
                for(int m=0; m<size; m++){
                    list.push(gene_node_map[ptr_1->ptr->getNextVector()->at(m)]);
                }
            }
            
            max=list.size();

        }
        
        


        
        
        //ptr->getInputVector
    }
    
    void print_node(Node* ptr, string hash){
        if(ptr!=nullptr&&ptr->get_is_parameter()){
            cout<<"This is a parameter node: "<<hash<<endl;

        }
       else if(ptr!=nullptr&&ptr->get_is_value_node()){
            cout<<"This is a value node: "<<hash<<endl;

        }
        else if(ptr!=nullptr){
            cout<<"This is a collection node: "<<hash<<endl;

        }
        

    }
    void print_connections(genes& m){
        cout<<"Connection node*****************"<<endl;
        for(int i=0; i<m.node_connections.size(); i++){
            print_node(genome->map[m.node_connections.at(i)],m.node_connections.at(i) );
        }
        cout<<"*************************"<<endl;
        cout<<endl;


    }
    void print(){
        cout<<"The size of the organism is:"<<genome->dominant_strand.size()<<endl;
        for(int i=0; i<genome->dominant_strand.size(); i++){
            cout<<"The hash of this node is: "<< genome->dominant_strand.at(i).hash<<endl;
            cout<<"The synapse wait time is "<< genome->dominant_strand.at(i).synapse_time<<endl;
            cout<<"The absolute time of this node is: "<< genome->dominant_strand.at(i).abs_time<<endl;

            if(genome->dominant_strand.at(i).expression!=-1&&genome->dominant_strand.at(i).is_parameter  ){
                cout<<"This is a parameter node"<<endl;
                cout<<"Active status: "<<genome->dominant_strand.at(i).is_active<<endl;
                if(genome->dominant_strand.at(i).is_active){
                    cout<<"The synapse number is: "<<genome->dominant_strand.at(i).current_time<<endl;
                    cout<<"First input: "<<genome->dominant_strand.at(i).ptr->getInputVector()->at(0)<<endl;
                }

                print_connections(genome->dominant_strand.at(i));
            }
            else if(genome->dominant_strand.at(i).expression!=-1&&genome->dominant_strand.at(i).is_value){
                cout<<"This is a value node"<<endl;
                cout<<"Active status: "<<genome->dominant_strand.at(i).is_active<<endl;
                if(genome->dominant_strand.at(i).is_active){
                    cout<<"The synapse number is: "<<genome->dominant_strand.at(i).current_time<<endl;
                    cout<<"First input: "<<genome->dominant_strand.at(i).ptr->getInputVector()->at(0)<<endl;
               }

                print_connections(genome->dominant_strand.at(i));
            }
            else if(genome->dominant_strand.at(i).expression!=-1){
                cout<<"This is a collection node"<<endl;
                cout<<"Active status: "<<genome->dominant_strand.at(i).is_active<<endl;
                cout<<"Collection Index: "<<genome->dominant_strand.at(i).collection_index<<endl;
                cout<<"The collection  time is "<< genome->dominant_strand.at(i).collective_time<<endl;
                if(genome->dominant_strand.at(i).is_active){
                    cout<<"The synapse number is: "<<genome->dominant_strand.at(i).current_time<<endl;
                    cout<<"Inputs:";
                    if(genome->dominant_strand.at(i).ptr->getInputVector()!=nullptr){
                        for(int m=0; m<genome->dominant_strand.at(i).ptr->getInputVector()->size(); m++){
                        cout<<genome->dominant_strand.at(i).ptr->getInputVector()->at(m)<<", ";
                        
                    }
                    cout<<endl;
                    }
                }

                print_connections(genome->dominant_strand.at(i));
            }
            else if(genome->dominant_strand.at(i).expression==-1){
                cout<<"This is a mutator node"<<endl;
                cout<<"Active status: "<<genome->dominant_strand.at(i).is_active<<endl;
                cout<<endl;
            }
        }
                
        
    }
    string return_print( ){
        string s;
        s=s+"The size of the organism is: "+ to_string(create_nodes.size());
        for(int i=0; i<create_nodes.size(); i++){
            s=s+"The hash of this node is: "+ genome->dominant_strand.at(i).hash;
            if(genome->dominant_strand.at(i).expression!=-1&&genome->dominant_strand.at(i).is_parameter  ){
                s=s+"This is a parameter node";
            }
            else if(genome->dominant_strand.at(i).expression!=-1&&genome->dominant_strand.at(i).is_value){
                s=s+"This is a value node";
            }
            else if(genome->dominant_strand.at(i).expression!=-1){
                s=s+"This is a collection node";
            }
        }
                
        return s;
        
    }

    void close(){// delets the dynamically appointed objects after creation
        for(int i=0; i<create_nodes.size(); i++){
            delete create_nodes.at(i);
        }
        genome->map.clear();
        create_nodes.clear();
    }

    bool status(){
        return genome->is_dead;
    }
    dna* getGenome(){
        return genome;
    }
   void set_genome(dna* genes){
       genome=genes;
    }
    void decrease_life(){
        (life)=life-1;
        time_since_reproduction++;
    }
    double get_life(){
        return life;
    }
    float get_score(){
        return score;
    }
    double get_time_since_reproduction(){
        return time_since_reproduction;
    }
    void set_time_since_reproduction(double d){
         time_since_reproduction=d;
    }
    vector<genes*>* get_active_genes(){
        return &active_genes;
    }

    ~organism(){
        delete genome;
    }
    void clear(){
        create_nodes.clear();
        active_genes.clear();
        gene_node_map.clear();
        optimal.clear();
        output_num.clear();
        score=0;
        genome->map.clear();
    }
private:
    dna* genome;
    vector<Node*> create_nodes;
    vector<genes*> active_genes;
    unordered_map<Node*,genes*> gene_node_map;
    bool reproduction_type;
    vector<double> optimal;
    vector<double> output_num;
    vector<double> in;
    vector<double> out;
    double life;
    double time_since_reproduction;
    int num_of_synapses;
    float score;
    Network<Node>* main;
    int size;

};


bool can_reproduce(organism* a, organism* b, double max_life){
    double life_amountA=a->get_life()/max_life;
    double life_amoundB=b->get_life()/max_life;
    if(life_amountA>=0.3&& life_amoundB>=0.3 && a->get_time_since_reproduction()> 15&& b->get_time_since_reproduction()>15){
        a->set_time_since_reproduction(-16);
        b->set_time_since_reproduction(-16);

        return true;
    }
    else{
        return false;
    }
}


bool compare(organism*& a, organism*& b){
    if(a->get_score()>b->get_score()){
        return true;
    }
    return false;
}
class food: public sim_objects{
public:
    
    food(int num):
    sim_objects(false){
        int count=0;
        for(int i=0; i< num; i++){
            int val=uniformTest(0, 1, 1);
            food_values.push_back(val);
            if(val!=0){
                count++;
            }
            
        }
        
        output_wanted.push_back(count);
    }

private:
    vector<double> food_values;
    vector<double> output_wanted;
    
    
    
};



struct obj_organism{

   
    sf::Shape* shape;
    genes* ptr;
    obj_organism(int x, int y, genes* pt){
        sf::CircleShape* circle=new sf::CircleShape;
        
        if(pt->is_value){
            circle->setFillColor(sf::Color(255, 0, 0));

        }
        else if(pt->is_parameter){
            circle->setFillColor(sf::Color(0, 255, 0));

        }
        else{
            circle->setFillColor(sf::Color(0, 0, 255));

        }
             
        circle->setRadius(10);
        circle->setOrigin(10, 10);

        circle->setPosition(x , y );
        
        shape=circle;
        ptr=pt;
        
    }
    
    
    void set_pos(int x, int y){
        shape->setPosition(x, y);
    }
    ~obj_organism(){
        delete shape;
    }

    
};
struct coordinates{
    int x;
    int y;
    
    coordinates(int this_x, int this_y){
        x=this_x;
        y=this_y;
    }
    coordinates(){
        
    }
    
};


struct line{
    int x;
    int y;
    int end_x;
    int end_y;
    
    line(int _x,int  _y, int _end_x, int _end_y ){
        x=_x;
        y=_y;
        end_x=_end_x;
        end_y=_end_y;
        
    }
    
sf::Vertex lines[2] =
    {
        sf::Vertex(sf::Vector2f(x, y)),
        sf::Vertex(sf::Vector2f(end_x, end_y))
    };
    
    

};

 void show_organism(organism& x){
    int s=x.get_active_genes()->size();
    obj_organism* array[1000][1000];
    
    for(int i=0; i<1000; i++){
        for(int m=0; m<1000; m++){
            array[i][m]=nullptr;
        }
    }
     
    bool found=false;
    
    std::unordered_map<string, coordinates> map;
    
    int count=0;
    while(!found && count<s){
        if(x.get_active_genes()->at(count)->expression==-1){
            found=true;
        }
        
        else{
        bool place_found=false;
        while(!place_found){
            int rand_rowX=uniformTest(100, 120, 1);
            int rand_rowY=uniformTest(100, 700, 1);
            if(array[rand_rowX][rand_rowY]==nullptr){
            array[rand_rowX][rand_rowY]=new obj_organism(rand_rowX, rand_rowY, x.get_active_genes()->at(count));
                place_found=true;
                coordinates c(rand_rowX, rand_rowY);
                map[x.get_active_genes()->at(count)->hash]=c;
                map[x.get_active_genes()->at(count)->hash].x=array[rand_rowX][rand_rowY]->shape->getPosition().x;
                map[x.get_active_genes()->at(count)->hash].y=array[rand_rowX][rand_rowY]->shape->getPosition().y;

            }
        }
        
            
        }
        count++;
    }
    found=false;
    while(!found && count<s){
        if(x.get_active_genes()->at(count)->expression==-1){
            found=true;
        }
        
        else{
        bool place_found=false;
        int rand_rowX=uniformTest(300, 400, 1);
        int rand_rowY=uniformTest(100, 700, 1);

        while(!place_found){
            if(array[rand_rowX][rand_rowY]==nullptr){
            array[rand_rowX][rand_rowY]=new obj_organism(rand_rowX, rand_rowY, x.get_active_genes()->at(count));
                place_found=true;
                coordinates c(rand_rowX, rand_rowY);
                map[x.get_active_genes()->at(count)->hash]=c;
                map[x.get_active_genes()->at(count)->hash].x=array[rand_rowX][rand_rowY]->shape->getPosition().x;
                map[x.get_active_genes()->at(count)->hash].y=array[rand_rowX][rand_rowY]->shape->getPosition().y;
            }
        }
        }
        count++;


    }
    
    found=false;

    
    while(count<s){
        if(x.get_active_genes()->at(count)->expression==-1){
            found=true;
        }
        else{
        bool place_found=false;
        int rand_rowX=uniformTest(800, 820, 1);
        int rand_rowY=uniformTest(100, 700, 1);

        while(!place_found){
            if(array[rand_rowX][rand_rowY]==nullptr){
            array[rand_rowX][rand_rowY]=new obj_organism(rand_rowX, rand_rowY, x.get_active_genes()->at(count));
                place_found=true;
                coordinates c(rand_rowX, rand_rowY);
                map[x.get_active_genes()->at(count)->hash]=c;
                map[x.get_active_genes()->at(count)->hash].x=array[rand_rowX][rand_rowY]->shape->getPosition().x;
                map[x.get_active_genes()->at(count)->hash].y=array[rand_rowX][rand_rowY]->shape->getPosition().y;
            }
        }
        }
        count++;

    }
    
    
    cout<<count<<":"<<s<<endl;

    vector< line > connections;
    

  
    for(int i=0; i<1000; i++){
        for(int m=0; m<1000; m++){
            if(array[i][m]!=nullptr){
            for(int s=0; s< array[i][m]->ptr->node_connections.size(); s++){
            
                if(array[i][m]->ptr->is_active && array[i][m]->ptr->node_connections.at(s)!=array[i][m]->ptr->hash){

                    
                    
                    auto it = map.find(array[i][m]->ptr->node_connections.at(s));
                    if (it == map.end()) {
                    }
                    else{
                        line l(i,m,map[array[i][m]->ptr->node_connections.at(s)].x, map[array[i][m]->ptr->node_connections.at(s)].y);
                        connections.push_back(l);

                    }
                }
                
                }
            }
        }
    }

    

    sf::RenderWindow window(sf::VideoMode(1000, 1000), "SFML Rectangle");
    while (window.isOpen())
    {

        count=0;

        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        
        window.clear();

        
        for (int i=0; i<1000 ; i++)
            {
            for (int m=0; m<1000 ; m++){
                if(array[i][m]!=nullptr){
                window.draw(*array[i][m]->shape);
                }
                }
            }
        
        
        for(int i=0; i< connections.size(); i++){
            window.draw(connections.at(i).lines, 2, sf::Lines);
        }

        
        
        window.display();
    }
    
    
    
    for(int i=0; i<1000; i++){
        for(int m=0; m<1000; m++){
            
            if(array[i][m]!=nullptr){
                
               delete array[i][m];
            }
        }
    }
     
    window.close();
    
    
}

 
static void show_text(string a[]){
    
    
    sf::RenderWindow window(sf::VideoMode(1000, 1000), "SFML Scrolling Text Example");

       // create a view with the same size as the window
       sf::View view(sf::FloatRect(0, 0, window.getSize().x, window.getSize().y));

       // set the view to the window
       window.setView(view);

       sf::Font font;
       font.loadFromFile("/Users/liyuzerihun/Documents/Neuroevolution_project1/oswald/Oswald-Bold.ttf");

       sf::Text text1(a[0], font, 10);
       text1.setFillColor(sf::Color::Green);
       text1.setPosition(50, 0);
    
        sf::Text text2(a[1], font, 10);
        text2.setFillColor(sf::Color::Green);
        text2.setPosition(300, 0);
    
    
        sf::Text text3(a[2], font, 10);
        text3.setFillColor(sf::Color::Green);
        text3.setPosition(600, 0);

       while (window.isOpen())
       {
           sf::Event event;
           while (window.pollEvent(event))
           {
               if (event.type == sf::Event::Closed)
                   window.close();
           }

           // handle keyboard input to scroll the view
           if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
           {
               view.move(0, -10); // move the view up
           }
           else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
           {
               view.move(0, 10); // move the view down
           }
           else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
           {
               view.move(-10, 0); // move the view left
           }
           else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
           {
               view.move(10, 0); // move the view right
           }

           window.clear();
           window.setView(view); // update the view in the window
           // draw the text object using the view to transform its position
           window.draw(text1);
           window.draw(text2);
           window.draw(text3);

           window.display();
       }

       return 0;

}


#endif /* Organisms_hpp */


