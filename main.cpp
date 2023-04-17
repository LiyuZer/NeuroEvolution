#include <iostream>

#include "Network.hpp"
#include "Node.hpp"
#include "organism.hpp"
 // Technique officially works on April 12 Wednesday 
using namespace std;

int main(int, char const**) {
    vector<organism*> creatures;
    Network<Node>* main = new Network<Node>(1);

    vector<double> food;
    vector<double> optimal;
    for (int i = 0; i < 2; i++) {
        int k = uniformTest(0, 1, 10);
        int mult = uniformTest(1, 10, 1);
        // int k=23;
        food.push_back(k);
        optimal.push_back(k);
    }

    int array[6] = {55,55,98,55,98,55};
    for (int i = 0; i < 1000; i++) {
        organism* ptr = new organism(100, main, array, 19);
        creatures.push_back(ptr);
        for (int m = 0; m < 10; m++) {
            mutate_dna(*creatures.at(i)->getGenome());
        }
        ptr->express_dna();
        ptr->synapse( food, optimal);
    }
    for (int round = 0; round < 100000; round++) {
        if (round != 0) {
            for (int k = 0; k < 30; k++) {
                food.clear();
                optimal.clear();
                for (int i = 0; i < 1; i++) {
                    if (round < 2000) {
                        int ks = uniformTest(1, 20, 1);
                        int mult = uniformTest(1, 10, 1);
                        food.push_back(ks);
                        food.push_back(mult);
                        optimal.push_back(ks);
                        optimal.push_back(mult);

                    } else {
                        int ks = uniformTest(1, 3, 1);
                        int mult = uniformTest(1, 3, 1);
                        int max = 0;
                        int less = 0;
                        if (mult < ks) {
                            max = ks;
                            less = mult;
                        } else {
                            max = mult;
                            less = ks;
                        }

                        food.push_back(ks);
                        food.push_back(mult);
                        optimal.push_back(max); 
                        optimal.push_back(less); 
                    }
                }
                for (int i = 0; i < 100; i++) {
                    organism* ptr = creatures.at(i);
                    if (k == 0) {
                        ptr->express_dna();
                    }
                    if (round < 2000) {
                        ptr->synapse( food, optimal);
                    } else {
                        ptr->synapse(food, optimal);
                    }
                }
            }
        }

        sort(creatures.begin(), creatures.end(), compare);

        int count = 0;
        if ((round + 1) % 100 == 0) {
            system("clear");
            creatures.at(0)->pr();
            creatures.at(0)->print();
        }
        for (int i = 0; i < 50; i++) {
            for (int k = 0; k < 2; k++) {
                int d = uniformTest(0,24, 1);
                int t = uniformTest(0,45, 1);
                organism* child = new organism(asexual_reproduce(*creatures.at(i)->getGenome()), main);
                //organism* child = new organism(reproduce(*creatures.at(t)->getGenome(), *creatures.at(d)->getGenome()), main);
                creatures.push_back(child);
            }
        }
        if (round == 0) {
            for (int i = 0; i < 1000; i++) {
                delete creatures.at(0);
                creatures.erase(creatures.begin());
            }
        }

        else {
            for (int i = 0; i < 100; i++) {
                delete creatures.at(0);
                creatures.erase(creatures.begin());
            }
        }
        main->clear(0);
    }

    return 0;
}
