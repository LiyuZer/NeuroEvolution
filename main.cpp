
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
// static void world_cycle(){
//     void* ptr;
//     cout<<"Which organism do you want to look at: "<<endl;
//     cin>>ptr;
//     show_organism(((organism*)ptr)->getGenome());
// }

int main(int, char const**) {
    //    organism m(500);
    //    vector<double> x;
    //
    //    cout<<endl;
    //    organism l(500);
    //    vector<double> t;
    //    cout<<endl;
    //
    //    vector<double> h;
    //
    //    organism child(500);
    //    child.set_genome(*reproduce(m.getGenome(), l.getGenome()));
    //    mutate_dna(child.getGenome());
    //    l.express_dna(t);
    //    m.express_dna(x);
    //    m.print();
    //    cout<<endl;
    //    l.print();
    //    cout<<endl;
    //
    //    child.express_dna(h);
    //    child.print();
    //
    //
    //
    //
    //    for(int i=0; i<100; i++){
    //        m.close();
    //        mutate_dna(m.getGenome());
    //        m.express_dna(x);
    //        show_organism(m);
    //
    //
    //    }

    //
    //    thread j(world_cycle);
    //    Universe world(10,100);
    //    world.cycle();

    //
    //    std::stringstream buffer;
    //    std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
    //
    //        organism m(500);
    //        organism l(500);
    //
    //        for(int i=0; i<100; i++){
    //            mutate_dna(*m.getGenome());
    //            mutate_dna(*l.getGenome());
    ////
    //        }
    //
    //        m.express_dna();
    //
    //    l.express_dna();
    //
    //    organism child(reproduce(*m.getGenome(), *l.getGenome()));
    //    child.express_dna();
    //
    //    show_organism(m);
    //    show_organism(l);
    //    show_organism(m);
    //
    //    show_organism(l);

    //    show_organism(child);

    //
    //    m.print();
    //    cout<<"&"<<endl;
    //    l.print();
    //    cout<<"&"<<endl;
    //    child.print();
    //    string s;
    //
    //
    //
    //
    //    std::string text = buffer.str(); // text will now contain "Bla\n"
    //
    //    string t[3]={""};
    //
    //    int indexLeft=0;
    //    int cn=0;
    //    for(int i=0; i<text.size(); i++){
    //        if(text.at(i)=='&'){
    //            t[cn]=text.substr(indexLeft, i-indexLeft);
    //            indexLeft=i;
    //
    //            cn++;
    //        }
    //
    //    }
    //    t[2]=text.substr(indexLeft, text.size()-indexLeft);
    //    std::cout.rdbuf(old);
    //
    //
    //    show_text(t);
    //
    //
    //
    //    vector<double> input;
    //
    //    for(int i=0; i<10; i++){
    //        input.push_back(i*10);
    //    }
    //
    //    vector<double>* ptr= m.synapse(input);
    //
    //    int d=ptr->size();
    //    for(int i=0; i<d; i++){
    //        cout<<ptr->at(i)<<endl;
    //    }

    //    Network<Node>* main= new Network<Node>(1);
    //
    ////    collection_node* ptr=new collection_node(4,main,true,0);
    ////    ptr->getInputVector()->push_back(1);
    ////    value_node* ptr1= new value_node(main, true,0);
    ////    ptr1->getInputVector()->push_back(0);
    ////
    ////    (*ptr)>>(*ptr1);
    ////    ptr->special_activation();
    ////    cout<<ptr1->getInput(0)<<endl;
    //    int array[5]={50,100,90,100,50};
    //    organism k(500, main, array,5);
    //    for(int m=0; m<10; m++){
    //        mutate_dna(*k.getGenome());
    //            }
    //    k.express_dna();
    //    k.synapse();
    //    k.print();
    //
    /*
      vector<organism*> creatures;
      Network<Node>* main= new Network<Node>(1);

      for(int i=0; i<1000; i++){
          organism* ptr= new organism(100, main);
          creatures.push_back(ptr);
          for(int m=0; m<10; m++){
              mutate_dna(*creatures.at(i)->getGenome());
          }
          ptr->express_dna();
          ptr->synapse(1);
      }
      for(int round=0; round<100000; round++){



          if(round!=0){
              for(int i=0; i<100; i++){
                  organism* ptr= creatures.at(i);
                  ptr->express_dna();
                  for(int k=0; k<1;k++){
                  ptr->synapse(1);
                  }
              }
          }


      sort(creatures.begin(), creatures.end(), compare);

          int count=0;
          if((round+1)%100==0){
          creatures.at(1)->pr();
              creatures.at(1)->print();

          }
          for(int i=0; i<5; i++){
              for(int k=0; k<18;k++){
                  int d=uniformTest(0, 10, 1);
              organism* child=new organism(reproduce(*creatures.at(i)->getGenome(),
      *creatures.at(d)->getGenome()),main); creatures.push_back(child);
              }
          }


          for(int i=5; i<10; i++){
              for(int k=0; k<2;k++){
              organism* child=new organism(reproduce(*creatures.at(i)->getGenome(), *creatures.at(uniformTest(5, 20,
      1))->getGenome()),main); creatures.push_back(child);
              }
          }

          if(round==0)
          {        for(int i=0; i< 1000; i++){
              delete creatures.at(count);
              creatures.erase(creatures.begin());
              }
          }

          else{
          for(int i=0; i< 100; i++){
              delete creatures.at(count);
              creatures.erase(creatures.begin());
          }
          }
          main->clear(0);

      }
   */

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
