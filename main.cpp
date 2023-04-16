
//
// Disclaimer:
// ----------
//
// This code will work only if you selected window, graphics and audio.
//
// Note that the "Run Script" build phase will copy the required frameworks
// or dylibs to your application bundle so you can execute it on any OS X
// computer.
//
// Your resource files (images, sounds, fonts, ...) are also copied to your
// application bundle. To get the path to these resources, use the helper
// function `resourcePath()` from ResourcePath.hpp
//

// #include <SFML/Audio.hpp>
// #include <SFML/Graphics.hpp>

// Here is a small helper for you! Have a look.
#include <iostream>

#include "Network.hpp"
#include "Node.hpp"
#include "organism.hpp"

using namespace std;


int main(int, char const**) {
   

    vector<organism*> creatures;
    Network<Node>* main = new Network<Node>(1);

    vector<double> food;
    vector<double> optimal;
    for (int i = 0; i < 2; i++) {
        int k = uniformTest(0, 1, 1);
        int mult = uniformTest(1, 10, 1);
        // int k=23;
        food.push_back(k);
        optimal.push_back(k);
    }

    int array[34] = {55, 55, 55, 55, 55, 55, 22, 22, 22, 33, 22, 55, 33, 2,  3,  3,  3,
                     2,  33, 98, 55, 5,  55, 55, 55, 55, 98, 33, 67, 77, 76, 76, 55, 55};
    for (int i = 0; i < 1000; i++) {
        organism* ptr = new organism(100, main, array, 19);
        creatures.push_back(ptr);
        for (int m = 0; m < 10; m++) {
            mutate_dna(*creatures.at(i)->getGenome());
        }
        ptr->express_dna();
        ptr->synapse(1, food, optimal);
    }
    for (int round = 0; round < 100000; round++) {
        if (round != 0) {
            for (int k = 0; k < 30; k++) {
                food.clear();
                optimal.clear();
                for (int i = 0; i < 1; i++) {
                    // int k=23;
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
                        optimal.push_back(k * mult);  // Technique officailly works on April 12 Wednesdayxr
                        optimal.push_back(k ^ mult);  // Technique officailly works on April 12 Wednesdayxr
                    }
                }
                for (int i = 0; i < 100; i++) {
                    organism* ptr = creatures.at(i);
                    if (k == 0) {
                        ptr->express_dna();
                    }
                    if (round < 2000) {
                        ptr->synapse(1, food, optimal);
                    } else {
                        ptr->synapse(2, food, optimal);
                    }
                }
            }
        }

        sort(creatures.begin(), creatures.end(), compare);

        int count = 0;
        if ((round + 1) % 100 == 0) {
            system("clear");
            creatures.at(1)->pr();
            creatures.at(1)->print();
        }
        for (int i = 0; i < 10; i++) {
            for (int k = 0; k < 8; k++) {
                int d = uniformTest(0, 10, 1);
                organism* child = new organism(asexual_reproduce(*creatures.at(i)->getGenome()), main);
                // organism* child=new organism(reproduce(*creatures.at(i)->getGenome(),
                // *creatures.at(d)->getGenome()),main);

                creatures.push_back(child);
            }
        }

        for (int i = 20; i < 40; i++) {
            for (int k = 0; k < 1; k++) {
                organism* child = new organism(asexual_reproduce(*creatures.at(i)->getGenome()), main);
                int d = uniformTest(0, 10, 1);
                // organism* child=new organism(reproduce(*creatures.at(i)->getGenome(),
                // *creatures.at(d)->getGenome()),main);
                creatures.push_back(child);
            }
        }

        if (round == 0) {
            for (int i = 0; i < 1000; i++) {
                delete creatures.at(count);
                creatures.erase(creatures.begin());
            }
        }

        else {
            for (int i = 0; i < 100; i++) {
                delete creatures.at(count);
                creatures.erase(creatures.begin());
            }
        }
        main->clear(0);
    }

    return 0;
}
