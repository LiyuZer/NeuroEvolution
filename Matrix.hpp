//
//  Matrix.hpp
//  openCVtest
//
//  Created by Liyu Zerihun on 5/25/22.
//

//#include <omp.h>
#include <string>
//#include "/opt/homebrew/Cellar/libomp/15.0.6/include/omp.h"

#include "Node.hpp"
#include "Network.hpp"
#include <vector>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <limits>
#include <thread>
#ifndef Matrix_hpp
#define Matrix_hpp

    

using namespace std;

template <typename T>  class Matrix{
    
 public:

    Matrix(int array[], int sizeOfArray, Network<Node>* n, bool isDynamic, int thrID){
        double size=1;
        int count=0;
        string s_temporary="";
        bool found=false;
        dynamic=isDynamic;

        for(int i=0; i<sizeOfArray; i=i+1){
            size=size*array[i];
            numbers.push_back(array[i]);
        }

        for(int i=0; i<size; i++){
            headers.push_back(T(n, false, thrID));
        }
        threadID=thrID;
    }
    
    
//    Matrix(int array[], T, int sizeOfArray, Network<Node>* n, bool isDynamic){
//        double size=1;
//        int count=0;
//        string s_temporary="";
//        bool found=false;
//        dynamic=isDynamic;
//
//        for(int i=0; i<sizeOfArray; i=i+1){
//            size=size*array[i];
//            numbers.push_back(array[i]);
//        }
//
//        for(int i=0; i<size; i++){
//            headers.push_back(T());
//        }
//
//
//    }
    
   
    
    
    
       Matrix(string s, Network<Node>* n, bool isDynamic, int threadI);


    Matrix( Matrix const& right){
        double size=1;
        int count=0;
        threadID=right.getThreadID();
        for(int i=0; i<right.getNumbers()->size(); i=i+1){
            size=size*right.getNumbers()->at(i);
            numbers.push_back(right.getNumbers()->at(i));
        }

        for(int i=0; i<size; i++){
            headers.push_back(T());
        }
        
        
        
        for(int i=0; i<headers.size(); i++){
            headers.at(i)=right.getHeaders()->at(i);
        }
    }

     int calc(vector<int> v, double index);
     bool checker(vector<int> first, vector<int> second);
     void set(string s, T& num);
     T& get(string s);
     void setA(int s[], int sizeOfArray, T& num);
    void setA_non_reference(int s[], int sizeOfArray, T num);
    void update(bool timeForUpdate, int threadID);
     T& getA(int s[], int sizeOfArray);
     vector<int>* returnDimensions();
    template <typename U> Matrix& correlate(Matrix<U>& kernel, int chan, Matrix<U>& bias, Matrix& returnMatrix);
    template <typename U> Matrix<U>& linearTransform(Matrix<U>& m, bool hidden);

    Matrix<value_node>& attach_collective_nodes( int indexFunc);//Attaches collectuved nodes to vecors
    Matrix<value_node>& addBias( parameter_node* m);//Adds Bias to the nodes

     Matrix* maxPool(Matrix pool);
    template <typename U> T& correlateMultiplication( int startingRow, int startingColoumn, int channel, Matrix<U>& kernel, int chan, Matrix<U>& bias);
    template <typename U>  Matrix<T>* convolution(Matrix<U>& kernel, Matrix<U>& bias);
    template <typename U>   vector<Matrix<T>*>& convolutionTest(Matrix<U>& kernel);

   template <typename U> void  addMatrixSlice(Matrix<U>& matrixSlice, int channel);

     double max( int startingRow, int startingColoumn, Matrix kernel);
    Matrix& operator+(Matrix& a);//Additon of MAtrices
    Matrix& operator^(Matrix& a);// Dot product of matrices

     void printMatrix();
    void printMatrixChannel(int channel);

     void fill2DMatrix();
     void reLu();
   void connectNode(Node& forward); // I added U for this special case where a matrix of value ndoes can be connected to collective nodes
    
    Matrix& convBranch(Matrix& Forward, int dimension, int channel);

    
    void connectNodeWithChannel(Node& forward, int n); // I added U for this special case where a matrix of value ndoes can be connected to collective nodes
    
    template <typename U> void connect(Matrix<U>& forwardMatrix); // I added U for this special case where a matrix of value ndoes can be connected to collective nodes
    template <typename U> void directConnect(Matrix<U>& forwardMatrix); // I added U for this special case where a matrix of value ndoes can be connected to corresponding nodes

    vector<T>* getHeaders() const{
        return &headers;
    }
    void choosefunction(int index){
        for(int i=0; i<headers.size(); i=i+1){
            headers.at(i).setIndex(index);
        }
    }
    void randomize(){
        for(int i=0; i<headers.size(); i++){
            T* n= new T(uniformRand(100), nullptr, this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
            headers.at(i)=*n;
        }
    }
    vector<int>* getNumbers() const{
        return &numbers;
    }
    void push_back(T element);
    void propogate();// Each eleemnet in the matrix will propogate it's value;
    ~Matrix();
    void clearMatrix();
    bool getDynamic(){
        return dynamic;
    }
    // void setSlice(int dimensions[], int dimensionIdex);
     // void getSlice(int dimensions[], int dimensionIdex);

    int getThreadID()const{
        return threadID;
    }
    void setThreadID(int k){
         threadID=k;
    }

 private:
     mutable vector<T> headers;
     mutable vector<int> numbers;
    bool dynamic;
    int threadID;

 };



template <typename T> Matrix<T>:: ~Matrix(){
    }

template <typename T> void Matrix<T>:: clearMatrix(){
    
    for(int i=0; i<headers.size(); i++){
        headers.at(i).zero();
    }
}

template <typename T> Matrix<T>:: Matrix (string s, Network<Node>* n,  bool isDynamic, int thrID){
   double size=1;
   int count=0;
   string s_temporary="";
   bool found=false;
   for(int i=0; i<s.size(); i=i+1){
       if(isdigit(s.at(i))){
           found=true;
           s_temporary.push_back(s.at(i));
       }
       else if(found==true){
           size=size*stoi(s_temporary);
           found=false;
           numbers.push_back(stoi(s_temporary));
           s_temporary="";
       }
   }

   for(int i=0; i<size; i++){
       headers.push_back(T(n, false, thrID));
   }
    threadID=thrID;
    
}


template <typename T>  void Matrix<T>:: propogate(){
    for(int i=0; i< headers.size(); i++){
        headers.at(i).activate();
    }
    
}
template <typename T> template <typename U>  void Matrix<T>:: connect(Matrix<U>& forwardMatrix){
    for(int m=0; m<headers.size(); m++){// Loops through the first node matrix(just a vector)
        for(int i=0; i< forwardMatrix.getHeaders()->size(); i++){// Then connects each node to all of teh forward nodes in the network(this is between two layers)
            headers.at(m)>>forwardMatrix.getHeaders()->at(i);
            }
        }
}
template <typename T>   void Matrix<T>:: connectNode(Node& forward){
    for(int m=0; m<headers.size(); m++){// Loops through the first node matrix(just a vector)
            headers.at(m)>>=forward;
        }
}

template <typename T>   void Matrix<T>:: connectNodeWithChannel(Node& forward, int channel){
    for(int m=0; m<numbers.at(0); m++){// Loops through the first node matrix(just a vector)
        for(int i=0; i<numbers.at(0); i++){// Loops through the first node matrix(just a vector)
            int array[3]={m,i,channel};
            getA(array,3)>>=forward;

           }
        }
}
template <typename T> template <typename U>  void Matrix<T>:: directConnect(Matrix<U>& forwardMatrix){
    for(int m=0; m<headers.size(); m++){// Loops through the first node matrix(just a vector)
            headers.at(m)>>forwardMatrix.getHeaders()->at(m);

        }
}

template <typename T>  Matrix<value_node>& Matrix<T>:: addBias( parameter_node* m){
    int array[1]={this->getNumbers()->at(0)};
    Matrix* n= new Matrix(array,1, m->getNetworkPointer(), true, threadID);
    
    int size1=this->getHeaders()->size();
    Network<Node>* main=m->getNetworkPointer();
    for(int i=0; i<size1 ; i++){
        int arr[1]={i};
        value_node* v= new value_node( main, true,threadID);
        *v=*m+this->getHeaders()->at(i);
        n->setA(arr, 1, *v);
    }
    return *n;
}



//template <typename T, typename U>  void mult(Matrix<T>* largeMatrix, Matrix<U>* m, int x, int location[]){
//
//
//    for(int x=0; x<largeMatrix->getNumbers()->at(1); x++){
//        int array[2]={i,x};
//        int arr[1]={x};
//        (m->getA(arr,1)*this->getA(array, 2))>>=*c;
//    }
//
//}

//template <typename T, typename U> void threadedlinearTransform(Matrix<T>* largeMatrix,Matrix<U>* m, int threadID, Matrix<U>* returnMatrix,  int* finished, vector<collection_node*>* coll){
    
    template <typename T, typename U> void threadedlinearTransform(Matrix<T>* largeMatrix,Matrix<U>* m, int threadID, Matrix<U>* returnMatrix){

    static std::mutex mutex2;


    int bound1=largeMatrix->getNumbers()->at(0);
    int start=threadID*((double)1/8)*bound1;
    int end=(threadID+1)*((double)1/8)*bound1;
    int bound2=largeMatrix->getNumbers()->at(1);
    Network<Node>* main=m->getHeaders()->at(0).getNetworkPointer();
    for(int i=start; i<end; i++){
        collection_node* c= new collection_node( main, true,largeMatrix->getThreadID());
        int location[1]={i};
        for(int x=0; x<bound2; x++){
            int array[2]={i,x};
            int arr[1]={x};
            (m->getA(arr,1)*largeMatrix->getA(array, 2))>>=*c;
        }


        U* temp= new U(main, true, largeMatrix->getThreadID());
        *c>>*temp;
        c->activate();
        returnMatrix->setA(location, 1, *temp);
 
 /*
    if(threadID<4){
    int bound1=largeMatrix->getNumbers()->at(0);
    int bound2=largeMatrix->getNumbers()->at(1);
    int start=threadID*((double)1/4)*bound2;
    int end=(threadID+1)*((double)1/4)*bound2;
    for(int i=0; i<bound1; i++){
        for(int x=start; x<end; x++){
            int array[2]={i,x};
            int arr[1]={x};
            (m->getA(arr,1)*largeMatrix->getA(array, 2))>>=*(coll->at(i));
        }
        mutex2.lock();
        finished[i]++;
        mutex2.unlock();

//        U* temp= new U(main, true);
//        *c>>*temp;
//        c->activate();
//        returnMatrix->setA(location, 1, *temp);
    }
    
    }
    
    else{
        threadID=threadID-4;
        int bound1=largeMatrix->getNumbers()->at(0);
        int bound2=largeMatrix->getNumbers()->at(1);
        int start=threadID*((double)1/4)*bound1;
        int end=(threadID+1)*((double)1/4)*bound1;
        Network<Node>* main=m->getHeaders()->at(0).getNetworkPointer();
        for(int i=start; i<end; i++){
            int location[1]={i};

            while(finished[i]<4){
            }
                   U* temp= new U(main, true);
                    *(coll->at(i))>>=*temp;
                    coll->at(i)->activate();
                    returnMatrix->setA(location, 1, *temp);

        }


//    }
//  */
    /*
    if(threadID<5){
    int bound1=largeMatrix->getNumbers()->at(0);
    int bound2=largeMatrix->getNumbers()->at(1);
    int start=threadID*((double)1/4)*bound2;
    int end=(threadID+1)*((double)1/4)*bound2;
    Network<Node>* main=m->getHeaders()->at(0).getNetworkPointer();
    for(int i=0; i<bound1; i++){
        collection_node* c= new collection_node( main, true);
        int location[1]={i};
        for(int x=start; x<end; x++){
            int array[2]={i,x};
            int arr[1]={x};
            (m->getA(arr,1)*largeMatrix->getA(array, 2))>>=*c;
        }
        bounds[i]++;


    }
        */
    
    
    }
    }

template <typename T> template <typename U> Matrix<U>& Matrix<T>:: linearTransform(Matrix<U>& m, bool hidden){
    int array[1]={this->getNumbers()->at(0)};
      Matrix<U>* returnMatrix=new Matrix<U>(array,1, m.getHeaders()->at(0).getNetworkPointer(), true, this->getThreadID());
      vector<thread> v;
      
      int bound1=this->getNumbers()->at(0);
      int bound2=this->getNumbers()->at(1);
      Network<Node>* main=m.getHeaders()->at(0).getNetworkPointer();
#pragma omp parallel
{
      for(int i=0; i<bound1; i++){
          collection_node* c= new collection_node( main, true, this->getThreadID());
          int location[1]={i};
#pragma omp critical
          {
          for(int x=0; x<bound2; x++){
              int array[2]={i,x};
              int arr[1]={x};
              (m.getA(arr,1)*this->getA(array, 2))>>=*c;
          }
      }
          
          U* temp= new U(main, true, this->getThreadID());
          *c>>*temp;
          c->activate();
          returnMatrix->setA(location, 1, *temp);
      }
}
      return *returnMatrix;
    
}
    
   /*
    if(hidden){
    int array[1]={this->getNumbers()->at(0)};
    Matrix<U>* returnMatrix=new Matrix<U>(array,1, m.getHeaders()->at(0).getNetworkPointer(), true, m.getThreadID());
    vector<thread> v;

    int bound1=this->getNumbers()->at(0);
    int bound2=this->getNumbers()->at(1);
    Network<Node>* main=m.getHeaders()->at(0).getNetworkPointer();
    
    for(int i=0; i<bound1; i++){
        collection_node* c= new collection_node( main, true, m.getThreadID());
        int location[1]={i};
        for(int x=0; x<bound2; x++){
            int array[2]={i,x};
            int arr[1]={x};
            (m.getA(arr,1)*this->getA(array, 2))>>=*c;
        }


        U* temp= new U(main, true, m.getThreadID());
        *c>>*temp;
        c->activate();
        returnMatrix->setA(location, 1, *temp);
    }
    return *returnMatrix;
   
    }
    else{
   
    int array[1]={this->getNumbers()->at(0)};

    Matrix<U>* returnMatrix=new Matrix<U>(array,1, m.getHeaders()->at(0).getNetworkPointer(), true, m.getThreadID());

    vector<thread> vec;
    for(int i=0; i<4; i=i+1){
        vec.push_back(thread(threadedlinearTransform<T,U>,this, &m, i, returnMatrix));
    }
    for(int i=0; i<4; i=i+1){
        vec.at(i).join();
    }
    return *returnMatrix;
    

//    int a[2]={this->getNumbers()->at(0),7};
//        Matrix<value_node> arr(a,2,m.getHeaders()->at(0).getNetworkPointer(), true);
//    vector<collection_node*> coll;
    
    Network<Node>* main=m.getHeaders()->at(0).getNetworkPointer();
    int array[1]={this->getNumbers()->at(0)};
//    for(int i=0; i<array[0]; i++ ){
//        collection_node* c= new collection_node(main, true);
//        coll.push_back(c);
//    }
//    int* finished= new int[this->getNumbers()->at(0)];
//    memset(finished, 0, this->getNumbers()->at(0) * sizeof(int));

    Matrix<U>* returnMatrix=new Matrix<U>(array,1, m.getHeaders()->at(0).getNetworkPointer(), true, threadID);
    vector<thread> vec;
    for(int i=0; i<8; i=i+1){
        vec.push_back(thread(threadedlinearTransform<T,U>,this, &m, i, returnMatrix));
    }
    for(int i=0; i<8; i=i+1){
        vec.at(i).join();
    }
    return *returnMatrix;
     */
       
    





template <typename T, typename U>  Matrix<U>&  linearTransform1D(Matrix<T>& largeMAtrix, Matrix<U>& m){
    int array[1]={largeMAtrix.getNumbers()->at(0)};
    Matrix<U>* returnMatrix=new Matrix<U>(array,1, m.getHeaders()->at(0).getNetworkPointer(), true, largeMAtrix.getThreadID());
    vector<thread> v;
    int rows=0;
    vector<collection_node*> collect;
    collection_node* c= new collection_node( m.getHeaders()->at(0).getNetworkPointer(), true);
    collect.push_back(c);
    for(int i=0; i<largeMAtrix.getHeaders()->size(); i++){
        int cols=(i+1)%largeMAtrix.getNumbers()->at(1);
        int array[2]={rows,cols};
        int arr[1]={cols};
        
        m.getA(arr,1)*largeMAtrix.getA(array, 2)>>=* collect.at(collect.size()-1);
        if(cols==0){
                int location[1]={rows};
                    U* temp= new U(m.getHeaders()->at(0).getNetworkPointer(), true);
                    *collect.at(collect.size()-1)>>*temp;
                    collect.at(collect.size()-1)->activate();
                    returnMatrix->setA(location, 1, *temp);
                    collection_node* f= new collection_node( m.getHeaders()->at(0).getNetworkPointer(), true);
                    collect.push_back(f);
                    rows++;

                }

      
            
        }
    return *returnMatrix;
    
}






template <typename T> Matrix<value_node>& Matrix<T>:: attach_collective_nodes( int indexFunc){
    int array[1]={this->getNumbers()->at(0)};

    Matrix<T>* matrix= new Matrix<T>(array,1, this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
    int size1=headers.size();
    for(int i=0; i<size1; i++){
        int array[1]={i};
        collection_node* c= new collection_node(indexFunc, this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
        value_node* v= new value_node(this->getHeaders()->at(0).getNetworkPointer(),true, threadID);
        headers.at(i)>>=*c;
        *c>>*v;
        c->activate();
        matrix->setA(array, 1, *v);
    }
    
    return *matrix;
    
    
}


template <typename T> Matrix<T>& Matrix<T>:: convBranch(Matrix& forward, int dimension, int channel){
    int array[3]={forward.getNumbers()->at(0),forward.getNumbers()->at(1),1};
    Matrix<T>* ptr=new Matrix(array,this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
    for(int i=0; i<forward.getNumbers()->at(0); i++){
        for(int m=0; m<forward.getNumbers()->at(1); m++){
            int array1[3]={forward.getNumbers()->at(0),forward.getNumbers()->at(1),1};
            int array2[3]={forward.getNumbers()->at(0),forward.getNumbers()->at(1),1};
            
        
    }
    
}
}


template <typename T> Matrix<T>& Matrix<T>:: operator+(Matrix& right){
    int array[numbers.size()];
    for(int i=0; i<numbers.size(); i++){
        array[i]=numbers.at(i);
    }
    Matrix * m= new Matrix(array, numbers.size());
    for(int i=0; i<headers.size(); i++){
        m->getHeaders()->at(i)=right.getHeaders()->at(i)+headers.at(i);
    }
    
    return *m;
}

template <typename T> Matrix<T>& Matrix<T>:: operator^(Matrix& right){
    int array[numbers.size()];
    for(int i=0; i<numbers.size(); i++){
        array[i]=numbers.at(i);
    }
    Matrix<T> * m= new Matrix<T>(array, numbers.size());
    collection_node* c= new collection_node(0, true, threadID);
    value_node* v= new value_node(this->getHeaders()->at(0).getNetworkPointer(), true);
    *c>>*v;
    for(int i=0; i<headers.size(); i++){
        m->getHeaders()->at(i)=right.getHeaders()->at(i)*this->getHeaders()->at(i);
        m->getHeaders()->at(i)>>=*c;

    }

    c->activate();

    
    Matrix<T>* oneDMatrix= new Matrix<T>("1");
    oneDMatrix->getHeaders()->at(0)=*v;

    return *oneDMatrix;
}

template <typename T>  int Matrix<T>:: Matrix :: calc(vector<int> v, double index){
   int sum=1;
   for(int i=index; i<v.size(); i++){
       int number=0;
       if(i==index){
           number=v.at(i)-1;
       }
       else{
           number=numbers.at(i);
       }
       sum=sum*number;
   }
   return sum;
}
template <typename T> bool Matrix<T>:: checker(vector<int> first, vector<int> second){
   for(int i=0; i<first.size(); i++){
       if(first.at(i)<second.at(i)){
           return false;
       }
   }
   return true;
}

template <typename T> void Matrix<T>::set(string s, T& num){
   vector<int> dimensions;
   string s_temporary="";
   bool found=false;
   for(int i=0; i<s.size(); i=i+1){

       if(isdigit(s.at(i))){
           found=true;
           s_temporary=s_temporary+s.at(i);
       }
       else if(found==true){
           dimensions.push_back(stoi(s_temporary)+1);
           found=false;
           s_temporary="";

       }
       
       
       
   }
   
   if(checker(numbers,dimensions)==false){
       throw invalid_argument("Invalid index");
   }
   int location=0;
   for(int i=0; i<dimensions.size(); i=i+1){
       location=location+calc(dimensions, i);
   }


   headers.at(location)=num;
}
template <typename T> T& Matrix<T>:: get(string s){
   vector<int> dimensions;
   string s_temporary="";
   bool found=false;

   for(int i=0; i<s.size(); i=i+1){
       
       if(isdigit(s.at(i))){
           found=true;
           s_temporary.push_back(s.at(i));
       }
       else if(found==true){
           dimensions.push_back(stoi(s_temporary)+1);
           found=false;
           s_temporary="";

       }
   }
   
   int location=0;
   if(checker(numbers,dimensions)==false){
       throw invalid_argument("Invalid index");
   }
   for(int i=0; i<dimensions.size(); i=i+1){

       location=location+calc(dimensions, i);
   }
   return headers.at(location);
}

template <typename T> void Matrix<T>:: setA_non_reference(int s[], int sizeOfArray, T num){
    vector<int> dimensions;
    for(int i=0; i<sizeOfArray; i=i+1){

        dimensions.push_back(s[i]+1);
    }
    
    if(checker(numbers,dimensions)==false){
        throw invalid_argument("Invalid index");
    }
    int location=0;
    for(int i=0; i<dimensions.size(); i=i+1){
        location=location+calc(dimensions, i);
    }
    headers.at(location)=num;
    
}

template <typename T> void Matrix<T>:: setA(int s[], int sizeOfArray, T& num){
   
   vector<int> dimensions;
   for(int i=0; i<sizeOfArray; i=i+1){

       dimensions.push_back(s[i]+1);
   }
   
   if(checker(numbers,dimensions)==false){
       throw invalid_argument("Invalid index");
   }
   int location=0;
   for(int i=0; i<dimensions.size(); i=i+1){
       location=location+calc(dimensions, i);
   }
   headers.at(location)=num;
}

template <typename T> void Matrix<T>:: update(bool timeForUpdate, int threadID){
   
    for(int i=0; i<headers.size(); i++){
        headers.at(i).update(timeForUpdate,threadID);
    }
}

template <typename T>  T& Matrix<T>:: getA(int s[], int sizeOfArray){
   vector<int> dimensions;

   for(int i=0; i<sizeOfArray; i=i+1){
       

           dimensions.push_back(s[i]+1);
   }
   
   int location=0;
   if(checker(numbers,dimensions)==false){
       throw invalid_argument("Invalid index");
   }
   for(int i=0; i<dimensions.size(); i=i+1){

       location=location+calc(dimensions, i);
   }
   return headers.at(location);
}
template <typename T> vector<int>* Matrix<T>:: returnDimensions(){
   return &numbers;
}
template <typename T> template <typename U> T& Matrix<T> :: correlateMultiplication( int startingRow,  int startingColoumn, int channel, Matrix<U>& kernel, int chan, Matrix<U>& bias){

    int array[2]={kernel.returnDimensions()->at(0), kernel.returnDimensions()->at(1)};
    T* returnSum= new value_node(this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
    collection_node* c= new collection_node(0,this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
   int iKernel=0;
   for(int i=startingRow; i<(startingRow+kernel.returnDimensions()->at(0));i=i+1){
       int mKernel=0;
       for(int m=startingColoumn; m<(startingColoumn+kernel.returnDimensions()->at(1));m=m+1){
           int array[3]={i,m,channel};
           int array2[3]={iKernel,mKernel, chan};
           (kernel.getA(array2, 3)*(getA(array, 3)))>>=*c;
           mKernel++;

       }
       iKernel++;

   }
    *c>>*returnSum;
    c->activate();
    collection_node* activationFunc= new collection_node(1 ,this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
    int biasArray[3]={0,0,chan};

    (*returnSum+bias.getA(biasArray, 3))>>=*activationFunc;
    T* finalSum= new value_node(this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
    *activationFunc>>*finalSum;
    activationFunc->activate();
    
    return (*finalSum);
}

template <typename T, typename U> void specialCorrelation( Matrix<T>* largeMatrix,int startingRow,  int startingColoumn, int channel, Matrix<U>& kernel, int chan, Matrix<U>& bias, collection_node* k){

    int array[2]={kernel.returnDimensions()->at(0), kernel.returnDimensions()->at(1)};
    T* returnSum= new value_node(largeMatrix->getHeaders()->at(0).getNetworkPointer(), true);
    collection_node* c= new collection_node(0,largeMatrix->getHeaders()->at(0).getNetworkPointer(), true);
   int iKernel=0;
   for(int i=startingRow; i<(startingRow+kernel.returnDimensions()->at(0));i=i+1){
       int mKernel=0;
       for(int m=startingColoumn; m<(startingColoumn+kernel.returnDimensions()->at(1));m=m+1){
           int array[3]={i,m,channel};
           int array2[3]={iKernel,mKernel, chan};
           (kernel.getA(array2, 3)*(largeMatrix->getA(array, 3)))>>=*c;
           mKernel++;

       }
       iKernel++;

   }
    *c>>*returnSum;
    c->activate();
    collection_node* activationFunc= new collection_node(1 ,largeMatrix->getHeaders()->at(0).getNetworkPointer(), true);
    int biasArray[3]={0,0,chan};

    (*returnSum+bias.getA(biasArray, 3))>>=*activationFunc;
    T* finalSum= new value_node(largeMatrix->getHeaders()->at(0).getNetworkPointer(), true);
    *activationFunc>>*finalSum;
    activationFunc->activate();
    
     (*finalSum)>>=*k;
}

template <typename T> template <typename U> Matrix<T>& Matrix<T>:: correlate(Matrix<U>& kernel, int chan, Matrix<U>& bias, Matrix& returnMatrix){
   int rowLength=(numbers.at(0) - kernel.returnDimensions()->at(0))+1;
   int coloumnLenght=(numbers.at(1) - kernel.returnDimensions()->at(1))+1;
  // Matrix* returnMatrix=new Matrix("("+to_string(rowLength)+"*"+to_string(coloumnLenght)+")", this->getHeaders()->at(0).getNetworkPointer(), true);
    //vector<thread> threads;
   for(int i=0; i< rowLength; i=i+1 ){
       for(int m=0; m< coloumnLenght; m=m+1){
           collection_node* c= new collection_node(0,this->getHeaders()->at(0).getNetworkPointer(), true);
           for(int channel=0; channel<this->getNumbers()->at(2); channel++)
           {
               //threads.push_back(specialCorrelation<T,U>,this, i, m, channel, kernel, chan, bias, c);
               
               T* sum= new value_node(nullptr, this->getHeaders()->at(0).getNetworkPointer(), true);
                *sum=correlateMultiplication(i, m, channel, kernel, chan, bias);
               *sum>>=*c;
//               if(threads.size()>2){
//                   for(int i=0; i<threads.size(); i++){
//                       threads.at(i).join();
//                   }
//                   threads.clear();
//               }
               
           }
          
//           for(int i=0; i<threads.size(); i++){
//               threads.at(i).join();
//           }
           
           T* returnSum= new value_node(this->getHeaders()->at(0).getNetworkPointer(), true);
           *c>>*returnSum;
           c->activate();
            int array[3]={i,m,chan};
           returnMatrix.setA(array, 3, *returnSum);
       }
   }
            return returnMatrix;
}




template <typename T> template <typename U> void Matrix<T>:: addMatrixSlice(Matrix<U>& matrixSlice, int channel){
    
    for(int i=0; i<this->getNumbers()->at(0); i++){
        for(int m=0; m<this->getNumbers()->at(1); m++){
            int array[3]={i,m,channel};
            int array2[3]={i,m};

            this->setA(array, 3, matrixSlice.getA(array2,2));
        }
        
    }
    delete &matrixSlice;
    
}


template <typename T,typename U> void correlateNonStatic(Matrix<T>* largeMatrix, Matrix<U>& kernel, int chan, Matrix<U>& bias, Matrix<T>& returnMatrix){
   int rowLength=(largeMatrix->getNumbers()->at(0) - kernel.returnDimensions()->at(0))+1;
   int coloumnLenght=(largeMatrix->getNumbers()->at(1) - kernel.returnDimensions()->at(1))+1;
  // Matrix* returnMatrix=new Matrix("("+to_string(rowLength)+"*"+to_string(coloumnLenght)+")", this->getHeaders()->at(0).getNetworkPointer(), true);
#pragma clang loop vectorize(assume_safety)
   for(int i=0; i< rowLength; i=i+1 ){
#pragma clang loop unroll_count(4)
       for(int m=0; m< coloumnLenght; m=m+1){
           collection_node* c= new collection_node(0,largeMatrix->getHeaders()->at(0).getNetworkPointer(), true);
#pragma clang loop unroll_count(4)
           for(int channel=0; channel<largeMatrix->getNumbers()->at(2); channel++)
           {
               
               T* sum= new value_node(nullptr, largeMatrix->getHeaders()->at(0).getNetworkPointer(), true);
                *sum=largeMatrix->correlateMultiplication(i, m, channel, kernel, chan, bias);
               *sum>>=*c;
           }
          
          
           T* returnSum= new value_node(largeMatrix->getHeaders()->at(0).getNetworkPointer(), true, largeMatrix->getThreadID());
           *c>>*returnSum;
           c->activate();
            int array[3]={i,m,chan};
           returnMatrix.setA(array, 3, *returnSum);
       }
   }
            
}


template <typename T> template <typename U>  Matrix<T>* Matrix<T>:: convolution(Matrix<U>& kernel, Matrix<U>& bias){
    int rowLength=(numbers.at(0) - kernel.returnDimensions()->at(0))+1;
    int coloumnLenght=(numbers.at(1) - kernel.returnDimensions()->at(1))+1;
    int array[3]={rowLength,coloumnLenght,kernel.getNumbers()->at(2)};
    
    Matrix<T>* v= new Matrix( array, 3, this->getHeaders()->at(0).getNetworkPointer(), true, threadID);
    vector<thread> threads;
    int numberThreads=0;
    for(int channel=0; channel<kernel.getNumbers()->at(2); channel++)
    {

      threads.push_back(thread(correlateNonStatic<T,U>, this, ref(kernel), channel, ref(bias), ref(*v)));
        numberThreads++;
        if(numberThreads>8){
            for(int i=0; i<threads.size(); i++){
                threads.at(i).join();
            }
            threads.clear();
            numberThreads=0;
        }
        
    }
    for(int i=0; i<threads.size(); i++){
        threads.at(i).join();
    }
    
    return v;
    
    
}
template <typename T> template <typename U>  vector<Matrix<T>*>& Matrix<T>:: convolutionTest(Matrix<U>& kernel){
   
    
    vector< Matrix<T>*> v;;
    for(int channel=0; channel<kernel.getNumbers()->at(2); channel++)
    {

        v.push_back(&correlate(kernel, channel));
    }

    
    return v;
    
    
}

template <typename T, typename U> Matrix<T> batchNorm(Matrix<T> largeMatrix, Matrix<U>& alpha, Matrix<U>& beta, vector <double> hiddenVals, Matrix<U>* finalMatrix){
    for(int i=0; i< largeMatrix.getHeaders()->size(); i++){
        hiddenVals.at(i)= largeMatrix.getHeaders()->at(i).getInput(0)+hiddenVals.at(i);
    }
    double mean;
    double variance;
    calculateMeandVariance(hiddenVals, mean, variance);
    
    double varianceNoise=sqrt(variance-0.01);
    for(int i=0; i< largeMatrix.getHeaders()->size(); i++){
        value_node numerator(largeMatrix.getHeaders()->at(i)->getInput(0),nullptr, largeMatrix.getHeaders()->at(0).getNetworkPointer(), true);
        value_node denominator(varianceNoise,nullptr, largeMatrix.getHeaders()->at(0).getNetworkPointer(), true);
        
       value_node* ptr = &(numerator/denominator);
        finalMatrix->getHeaders()->at(i)= (*ptr)*alpha.getHeaders()->at(i)+ beta.getHeaders()->at(i);
    }
    
    
}

template <typename T> void Matrix<T>:: printMatrix(){
   for(int i=0; i<returnDimensions()->at(0); i++){
       int n;
       if( returnDimensions()->size()<2){
           n=1;
       }
       else{
           n=returnDimensions()->at(1);
       }
       for(int m=0; m<n; m++){
           if(n==1){
               int array[1]={i};
               
             // cout<<get("("+to_string(i)+"*"+to_string(m)+")");
               cout<<getA(array, 1).getInput(0)<<" ";
           }
           else{
               int array[2]={i,m};
               
             // cout<<get("("+to_string(i)+"*"+to_string(m)+")");
               
               cout<<getA(array, 2).getInput(0)<<" ";
           }
           

   }
       cout<<endl;
   
}
}
template <typename T> void Matrix<T>:: printMatrixChannel(int channel){
   for(int i=0; i<returnDimensions()->at(0); i++){
       int n;
       if( returnDimensions()->size()<2){
           n=1;
       }
       else{
           n=returnDimensions()->at(1);
       }
       for(int m=0; m<n; m++){
          
           
               int array[3]={i,m, channel};
               
             // cout<<get("("+to_string(i)+"*"+to_string(m)+")");
               
               cout<<getA(array, 3).getInput(0)<<" ";
           }
           
       cout<<endl;

   }
   
}



template <typename T> void Matrix<T>:: fill2DMatrix(){
  
   
   for(int i=0; i<numbers.at(0); i++){
       for(int m=0; m<numbers.at(1); m++){
           std::random_device dev;
              std::mt19937 rng(dev());
              std::uniform_int_distribution<std::mt19937::result_type> dist6(-1,1); // distribution in range [1, 6]
             
           
           int array[2]={i,m};
           setA(array, 2,  dist6(rng));
       }
   }
   
}
template <typename T> double Matrix<T> :: max( int startingRow, int startingColoumn, Matrix<T> kernel){
   int initialArray[2]={startingRow,startingColoumn};
   double greatest=getA(initialArray, 2);
   for(int i=startingRow; i<(startingRow+kernel.returnDimensions()->at(0));i=i+1){
       for(int m=startingColoumn; m<(startingColoumn+kernel.returnDimensions()->at(1));m=m+1){
           int array[2]={i,m};
           int tempNum= getA(array, 2);
           if( tempNum >greatest){
               greatest=tempNum;
           }

       }

   }
   return greatest;
}


template <typename T> Matrix<T>* Matrix<T>:: maxPool(Matrix<T> pool){
   int rowLength=(numbers.at(0) - pool.returnDimensions()->at(0))+1;
   int coloumnLenght=(numbers.at(1) - pool.returnDimensions()->at(1))+1;
   Matrix* returnMatrix=new Matrix("("+to_string(rowLength)+"*"+to_string(coloumnLenght)+")");
   for(int i=0; i< rowLength; i=i+1 ){
       for(int m=0; m< coloumnLenght; m=m+1){
           int greatest=max(i, m, pool);
           int array[2]={i,m};
           returnMatrix->setA(array, 2, greatest);
       }
   }
            return returnMatrix;
   
}

template <typename T> void Matrix<T>:: reLu(){
   for(int i=0; i< numbers.at(0); i=i+1 ){
       for(int m=0; m< numbers.at(1); m=m+1){
           int array[2]={i,m};
           if(getA(array, 2)<0){
               setA(array, 2, 0);
           }
       }
   }
   
   
}
template <typename T> void test (Matrix<T>& kernel){
    for(int i=0; i<kernel.getHeaders()->size(); i++){
        cout<<kernel.getHeaders()->at(i).getInput(0)<<endl;
    }
}

template <typename T> void specialPrint (Matrix<T>& kernel){
    int count=0;
    for(int i=0; i<kernel.getHeaders()->size(); i++){
        
        if(kernel.getHeaders()->at(i)==nullptr){
            cout<<"0"<<" ";
        }
        
        else{
            cout<<kernel.getHeaders()->at(i)->getInput(0)<<" ";

        }
        count++;
        if(count==kernel.getNumbers()->at(1)){
            count=0;
            cout<<"|||"<<endl;
        }
    }
    cout<<endl;
}

template <typename T, typename U> Matrix<U*>& sparseMatrixForm(Matrix<T>& kernel, int array[], int size,vector<vector<int> >& positions, int valuesPerRows){
    int numberOfWaits=0;
    int sapceBetweenVal=sqrt(array[1]) - (kernel.getNumbers()->at(1));
    int countVals=0;
    int indexOfKernel=0;
    int kernelSize=kernel.getHeaders()->size();
    Matrix<value_node*>* sparseMatrix= new Matrix<value_node*>(array,nullptr,size,kernel.getHeaders()->at(0).getNetworkPointer(),true);

    for(int rows=0; rows<array[0]; rows++){
        vector<int> m;
        positions.push_back(m);
        for(int cols=numberOfWaits; cols<array[1]; cols++){
            int array2[2]={rows, cols};
            value_node* val=new value_node(kernel.getHeaders()->at(0).getNetworkPointer(),true);
            kernel.getHeaders()->at(indexOfKernel)|*val;
            sparseMatrix->setA_non_reference(array2, 2, val);
            positions.at(positions.size()-1).push_back(cols);
            countVals++;
            indexOfKernel++;
            if(indexOfKernel==kernelSize){
                if(cols==array[1]-1){
                    rows=array[0];
                }
                indexOfKernel=0;
                cols=array[1];
            }
            if(countVals==kernel.getNumbers()->at(1)){
                cols=cols+sapceBetweenVal;
                countVals=0;
                
            }
        }
        if((rows+1)%valuesPerRows==0 && rows!=0){
            numberOfWaits=numberOfWaits+(kernel.getNumbers()->at(1));

        }
        else{
            numberOfWaits++;
        }
        
    }
    return *sparseMatrix;
    
    
    
}

template <typename T> void Matrix<T>:: push_back(T element){
   
   
   
}

template <typename T, typename U> void linearThread(Matrix<T>* largeMatrix, Matrix<U>* m,  Matrix<T>* finalMatrix, Matrix<T*>* sparseMatrix, int start,int end, int index,vector<vector<int> >& positions  ){

    
    
    for(int rows=start; rows< end; rows++){
        
        collection_node* c= new collection_node(0, m->getHeaders()->at(0).getNetworkPointer(), true);
        for(int cols=0; cols<positions.at(rows).size(); cols++){
            int sparseCol=positions.at(rows).at(cols);
            int arraySparse[2]={rows, sparseCol};
            ((largeMatrix->getHeaders()->at(sparseCol)*(*sparseMatrix->getA(arraySparse, 2))))>>=*c;
        }

        
        *c>>finalMatrix->getHeaders()->at(index);
        c->activate();
        index++;
        
    }

    
    
}


template <typename T, typename U> void reverse(Matrix<T>* largeMatrix, Matrix<U>* m,  Matrix<T>* finalMatrix, Matrix<T*>* sparseMatrix, int start,int end, int index,vector<vector<int> >& positions  ){

    
    
    for(int rows=start; rows< end; rows++){
        
        collection_node* c= new collection_node(0, m->getHeaders()->at(0).getNetworkPointer(), true);
        for(int cols=positions.at(rows).size()-1; cols>=0; cols=cols-1){
            int sparseCol=positions.at(rows).at(cols);
            int arraySparse[2]={rows, sparseCol};
            ((largeMatrix->getHeaders()->at(sparseCol)*(*sparseMatrix->getA(arraySparse, 2))))>>=*c;
        }

        
        *c>>finalMatrix->getHeaders()->at(index);
        c->activate();
        index++;
        
    }

    
    
}



template <typename T, typename U> Matrix<T>& linearSpecial(Matrix<T>& largeMatrix, Matrix<U>& m){
    
    int rowLength=(largeMatrix.getNumbers()->at(0) - m.returnDimensions()->at(0))+1;
    int coloumnLenght=(largeMatrix.getNumbers()->at(1) - m.returnDimensions()->at(1))+1;
    int array[2]={rowLength*coloumnLenght, static_cast<int>(largeMatrix.getHeaders()->size())};
    int dimension[2]={rowLength, coloumnLenght};
    vector<vector<int> > positions;
    
    
   Matrix<T*>* sparseMatrix= &sparseMatrixForm<parameter_node, value_node>(m,array,2,positions,dimension[1]);
    Matrix<T>* finalMatrix=new Matrix<T>(dimension,2,m.getHeaders()->at(0).getNetworkPointer(), true);
   
    
    
    int val=array[0];
    int div=8;
    vector<thread> threads;
    for(int i=0; i<12; i++){
        threads.push_back(thread(linearThread<T,U>, &largeMatrix,&m, finalMatrix, sparseMatrix, val*(i/div), val*((i+1)/div), val*(i/div), ref( positions)));

    }
    
   


    
    for(int i=0; i<threads.size(); i++){
        threads.at(i).join();
    }
    delete sparseMatrix;
    
    return *finalMatrix;
    
    
  
}

 template <typename T> void million (Network<Node>* main,collection_node* c){
    for(double i=0; i<(20000000/8)+1; i++){
        int k=sqrt(2);

        value_node* ptr= new value_node(1, nullptr, main, false);
        *ptr>>=*c;

        
    }

    
}

template <typename T> void no (Network<Node>* main,collection_node* c){
   for(double i=0; i<20000000; i++){
       int k=sqrt(2);

       value_node* ptr= new value_node(1, nullptr, main, false);
       *ptr>>=*c;

       
   }

   
}


template <typename T> void billion (Network<Node>* main,collection_node* c){
   for(double i=0; i<(20000000/8)+1; i++){
       int k=sqrt(2);
       value_node* ptr= new value_node(1, nullptr, main, false);
       *ptr>>=*c;

       
   }
   
}
 template <typename T> double fun(Network<Node>* main){
    collection_node* c= new collection_node(0,main, true);
     collection_node* c1= new collection_node(0,main, true);
     collection_node* c2= new collection_node(0,main, true);
     collection_node* c3= new collection_node(0,main, true);
     collection_node* c4= new collection_node(0,main, true);
     collection_node* c5= new collection_node(0,main, true);
     collection_node* c6= new collection_node(0,main, true);
     collection_node* c7= new collection_node(0,main, true);

    thread t1(million<T>,main,c);
     thread t2(billion<T>,main,c1);
     thread t3(billion<T>,main,c2);
     thread t4(billion<T>,main,c3);
     thread t5(billion<T>,main,c4);
     thread t6(billion<T>,main,c5);
     thread t7(billion<T>,main,c6);
     thread t8(billion<T>,main,c7);
     
    
     t1.join();
     t2.join();
     t3.join();
     t4.join();
     t5.join();
     t6.join();
     t7.join();
     t8.join();

     collection_node specialc(0,main, true);
     
     *c>>specialc;
     *c1>>specialc;
     *c2>>specialc;
     *c3>>specialc;
     *c4>>specialc;
     *c5>>specialc;
     *c6>>specialc;
     *c7>>specialc;
     
     c->activate();
     c1->activate();
     c2->activate();
     c3->activate();
     c4->activate();
     c5->activate();
     c6->activate();
     c7->activate();
     specialc.activate();
    return specialc.getOutput(0);
}


///
///
///

#endif /* Matrix_hpp */



