////
////  Universe.hpp
////  Neuroevolution_project1
////
////  Created by Liyu Zerihun on 1/25/23.
////  Copyright © 2023 Liyu Zerihun. All rights reserved.
////
//
//#ifndef Universe_hpp
//#define Universe_hpp
//
//#include "organism.hpp"
//#include <SFML/Audio.hpp>
//#include <SFML/Graphics.hpp>
//#include "/Users/liyuzerihun/Documents/MLT/Utilities.cpp"
//
//
//
//const int max_row_count=1000;
//const int max_col_count=1000;
//const int momentum_organism=0.9;
//const double organism_life=20000;
//struct obj{
//    double* life;
//
//    organism* object;
//    food* food_object;
//    sf::Shape* shape;
//    double prev_direction;
//    
//    obj(int x, int y){
//        sf::CircleShape* circle=new sf::CircleShape;
//        circle->setFillColor(sf::Color(0, 255, 0));
//        circle->setPosition(x, y);
//        circle->setRadius(10);
//        shape=circle;
//        organism* ptr=new organism(10000);
//        object=ptr;
//
//
//    }
//    
//    obj(int x, int y, bool val){
//        sf::RectangleShape* s=new sf::RectangleShape;
//        s->setFillColor(sf::Color(200, 0, 100));
//        s->setPosition(x, y);
//        s->setSize(sf::Vector2f(10,10));
//        shape=s;
//        food* ptr=new food(uniformTest(0, 10, 1));
//        food_object=ptr;
//        object=nullptr;
//
//
//    }
//    
//    obj(int x, int y, dna* genome){
//        sf::CircleShape* circle=new sf::CircleShape;
//        circle->setFillColor(sf::Color(0, 255, 0));
//        circle->setPosition(x, y);
//        circle->setRadius(10);
//        shape=circle;
//        organism* ptr=new organism( genome);
//        object=ptr;
//        
//
//
//
//    }
//    
//    
//    void set_pos(int x, int y){
//        shape->setPosition(x, y);
//        if(object!=nullptr){
//            int color= double((organism_life-object->get_life())/organism_life)*255;
//            shape->setFillColor(sf::Color(0,255,color ));
//
//        }
//    }
//
//    
//};
//class Universe{
//    
//public:
//    Universe(int number_organisms, int number_food){
//        
//        
//        
//        for(int i=0; i<max_row_count; i++){
//            for(int m=0; m<max_col_count; m++){
//                universe_objects[i][m]=nullptr;
//            }
//            
//        }
//        for(int i=0; i< number_organisms; i++){
//
//                int r=uniformTest(0, max_row_count-1,1);
//                int c=uniformTest(0, max_col_count-1,1);
//            bool found=false;
//            while(!found){
//                 r=uniformTest(0, max_row_count-1,1);
//                 c=uniformTest(0, max_col_count-1,1);
//                if(universe_objects[r][c]==nullptr){
//                    found=true;
//                    universe_objects[r][c]=new obj(r,c);
//                    cout<<universe_objects[r][c]<<endl;
//                }
//            }
//        }
//    
//    
//    for(int i=0; i< number_food; i++){
//
//            int r=uniformTest(0, max_row_count-1,1);
//            int c=uniformTest(0, max_col_count-1,1);
//        bool found=false;
//        while(!found){
//             r=uniformTest(0, max_row_count-1,1);
//             c=uniformTest(0, max_col_count-1,1);
//            if(universe_objects[r][c]==nullptr){
//                found=true;
//                universe_objects[r][c]=new obj(r,c,false);
//            }
//            
//            }
//            
//        }
//        
//    }
//        
//    
//    void move( obj* ptr, int x, int y){
//        int rand_direction=uniformTest(0, 50, 1);
//        int speed=1;
//        
//        if(rand_direction==45){
//            rand_direction=uniformTest(0, 3, 1);
//            ptr->prev_direction=rand_direction;
//        }
//        else{
//            rand_direction=ptr->prev_direction;
//        }
//        
//        
//        if(rand_direction==0&& x-speed>0 ){
//            if(universe_objects[x-speed][y]==nullptr){
//                ptr->set_pos(x-speed, y);
//                universe_objects[x][y]=nullptr;
//                universe_objects[x-speed][y]=ptr;
//            }
//            else if(universe_objects[x-speed][y]->object!=nullptr){
//                int r=uniformTest(0, max_row_count-1,1);
//                int c=uniformTest(0, max_col_count-1,1);
//                bool found=false;
//                if(can_reproduce(universe_objects[x-speed][y]->object,universe_objects[x][y]->object, organism_life)){
//                while(!found){
//                r=uniformTest(0, max_row_count-1,1);
//                c=uniformTest(0, max_col_count-1,1);
//                if(universe_objects[r][c]==nullptr){
//                    found=true;
//                    dna* ptr=reproduce(*universe_objects[x][y]->object->getGenome(), *universe_objects[x-speed][y]->object->getGenome());
//                    universe_objects[r][c]=new obj(r,c, ptr);
//                }
//                
//            }
//                }
//            
//        }
//        }
//        else if(rand_direction==1 && x+speed<max_row_count ) {
//            
//            if( universe_objects[x+speed][y]==nullptr){
//                ptr->set_pos(x+speed, y);
//                universe_objects[x][y]=nullptr;
//                universe_objects[x+speed][y]=ptr;
//            }
//            else if(universe_objects[x+speed][y]->object!=nullptr){
//                
//                int r=uniformTest(0, max_row_count-1,1);
//                int c=uniformTest(0, max_col_count-1,1);
//            bool found=false;
//                if(can_reproduce(universe_objects[x+speed][y]->object,universe_objects[x][y]->object, organism_life)){
//            while(!found){
//                 r=uniformTest(0, max_row_count-1,1);
//                 c=uniformTest(0, max_col_count-1,1);
//                if(universe_objects[r][c]==nullptr){
//                    found=true;
//                    dna* ptr=reproduce(*universe_objects[x][y]->object->getGenome(), *universe_objects[x+speed][y]->object->getGenome());
//                    universe_objects[r][c]=new obj(r,c, ptr);
//                }
//                
//                
//            }
//                }
//            
//        }
//    }
//        else if(rand_direction==2&& y-speed>0){
//            
//            if( universe_objects[x][y-speed]==nullptr){
//                ptr->set_pos(x, y-speed);
//                universe_objects[x][y]=nullptr;
//                universe_objects[x][y-speed]=ptr;
//            }
//            else if(universe_objects[x][y-speed]->object!=nullptr){
//                
//                int r=uniformTest(0, max_row_count-1,1);
//                int c=uniformTest(0, max_col_count-1,1);
//            bool found=false;
//                if(can_reproduce(universe_objects[x][y-speed]->object,universe_objects[x][y]->object, organism_life)){
//            while(!found){
//                 r=uniformTest(0, max_row_count-1,1);
//                 c=uniformTest(0, max_col_count-1,1);
//                if(universe_objects[r][c]==nullptr){
//                    found=true;
//                    dna* ptr=reproduce(*universe_objects[x][y]->object->getGenome(), *universe_objects[x][y-speed]->object->getGenome());
//                    universe_objects[r][c]=new obj(r,c, ptr);
//                }
//                
//                
//            }
//                }
//            
//        }
//        }
//        else if(rand_direction==3&& y+speed<max_row_count){
//            
//            if( universe_objects[x][y+speed]==nullptr){
//                ptr->set_pos(x, y+speed);
//                universe_objects[x][y]=nullptr;
//                universe_objects[x][y+speed]=ptr;
//            }
//            else if(universe_objects[x][y+speed]->object!=nullptr){
//                
//                int r=uniformTest(0, max_row_count-1,1);
//                int c=uniformTest(0, max_col_count-1,1);
//            bool found=false;
//            if(can_reproduce(universe_objects[x][y+speed]->object,universe_objects[x][y]->object, organism_life)){
//            while(!found){
//                 r=uniformTest(0, max_row_count-1,1);
//                 c=uniformTest(0, max_col_count-1,1);
//                if(universe_objects[r][c]==nullptr){
//                    found=true;
//                    dna* ptr=reproduce(*universe_objects[x][y]->object->getGenome(), *universe_objects[x][y+speed]->object->getGenome());
//                    universe_objects[r][c]=new obj(r,c, ptr);
//                }
//                
//                
//            }
//            }
//            
//        }
//        }
//        
//    }
//    
//    
//
//   void cycle(){
//       sf::RenderWindow window(sf::VideoMode(1000, 1000), "SFML Rectangle");
//       sf::Text text;
//       while (window.isOpen())
//       {
//           count=0;
//
//           
//           sf::Event event;
//           while (window.pollEvent(event))
//           {
//               if (event.type == sf::Event::Closed)
//                   window.close();
//           }
//           window.clear();
//           
//           
//           for (int i=0; i<max_row_count ; i++)
//               {
//                   for (int m=0; m<max_col_count ; m++){
//                       if(universe_objects[i][m]!=nullptr && universe_objects[i][m]->object!=nullptr){
//                           if((universe_objects[i][m]->object->get_life())<=0){
//                           delete universe_objects[i][m];
//                           universe_objects[i][m]=nullptr;
//                           }
//                           else{
//                               universe_objects[i][m]->object->decrease_life();
//                           }
//                       }
//                       if(universe_objects[i][m]!=nullptr && universe_objects[i][m]->object!=nullptr){
//                           move(universe_objects[i][m], i, m);
//                       }
//                   }
//               }
//           for (int i=0; i<max_row_count ; i++)
//               {
//                   for (int m=0; m<max_col_count ; m++){
//
//                       if(universe_objects[i][m]!=nullptr){
//                           window.draw(*universe_objects[i][m]->shape);
//                           if(universe_objects[i][m]->object!=nullptr){
//                               count++;
//
//                           }
//                       }
//                   }
//               }
//           
////
////           if(event.type==sf::Event::MouseButtonPressed){
////               if (event.mouseButton.button == sf::Mouse::Left){
////                   cout<<"I am here"<<endl;
////                   sf::Vector2i mousePos = sf::Mouse::getPosition(window);
////                   if(universe_objects[mousePos.x][mousePos.y]!=nullptr){
////
////                       cout<<universe_objects[mousePos.x][mousePos.y]->object->return_print()<<endl;;
////
////                   }
////
////               }
////
////           }
//           window.display();
//          // cout<<"The number of organisms is: "<<count<<endl;
//
//       }
//    }
//    
//    
//    
//private:
//    obj* universe_objects[max_row_count][max_col_count];
//    int count;
//    
//    
//};
//
//
//
//
//#endif /* Universe_hpp */
