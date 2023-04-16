//
//  Network.hpp
//  openCVtest
//
//  Created by Liyu Zerihun on 8/17/22.
//

#ifndef Network_hpp
#define Network_hpp

#include <stdio.h>

#include <iostream>
#include <vector>
using namespace std;
template <typename T>

class Network {
   public:
    mutex* mutel[8];

    Network(int threadCount) {
        parameters = 0;
        values = 0;
        collective = 0;
        for (int i = 0; i < threadCount; i++) {
            vector<T*>* f = new vector<T*>;
            f->clear();
            Nodes.push_back(*f);
            mutex* m = new mutex;
            mutel[i] = m;
        }
    }
    vector<T*>* getNodes() { return &Nodes; }

    void push(T* num, int threadID) {
        mutel[threadID]->lock();
        Nodes.at(threadID).push_back(num);
        mutel[threadID]->unlock();
    }

    void clear(int threadID) {
        for (int i = 0; i < Nodes.at(threadID).size(); i++) {
            delete Nodes.at(threadID).at(i);
        }

        Nodes.at(threadID).clear();
    }
    ~Network() {}

   private:
    vector<vector<T*>> Nodes;

    int parameters;
    int values;
    int collective;
};

#endif /* Network_hpp */
