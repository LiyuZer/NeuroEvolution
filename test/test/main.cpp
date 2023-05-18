//
//  main.cpp
//  test
//
//  Created by Liyu Zerihun on 5/12/23.
//

#include <iostream>
#include <string>
#include <cmath>

using namespace std;
//int main(int argc, const char * argv[]) {
//    // insert code here...
//
//    int array[5]={1,2,3,4,5};
//    int num=500;
//    int dp[500];
//
//    for(int i=0; i<500; i++){
//        for(int m=0; m<5; m++){
//            if(dp[]-array[m])
//        }
//
//    }
//
//    return 0;
//}
//
//int dynamic_distance(string x, string y){
//
//    int** array= new int*[x.size()+1];
//
//    for(int i=0; i<x.size()+1; i++){
//        array[i]=new int[y.size()+1];
//    }
//
//
//    for(int i=0; i<y.size()+1; i++){
//        array[0][i]=i;
//    }
//    for(int i=0; i<x.size()+1; i++){
//        array[i][0]=i;
//    }
//
//    for(int r=1; r<x.size()+1; r++){
//        for(int c=1; c<y.size()+1; c++){
//            array[r][c]=0;
//        }
//    }
//
//
//    for(int r=1; r<x.size()+1; r++){
//        for(int c=1; c<y.size()+1; c++){
//
//            if(x.at(r-1)!=y.at(c-1)){
//                int replace=array[r-1][c-1]+1;
//                int delete_letter=array[r-1][c]+1;
//                int insert=array[r][c-1]+1;
//
//                if(replace<delete_letter){
//
//                    if(replace<insert){
//
//                        array[r][c]=replace;
//                    }
//                    else{
//                    array[r][c]=insert;
//                    }
//                }
//                else{
//                    if(delete_letter<insert){
//                        array[r][c]=delete_letter;
//                    }
//                    else{
//                        array[r][c]=insert;
//                    }
//                }
//            }
//            else{
//                array[r][c]=array[r-1][c-1];
//            }
//
//        }
//
//    }
//
//
////    for(int r=0; r<x.size()+1; r++){
////        for(int c=0; c<y.size()+1; c++){
////            cout<<array[r][c]<<",";
////        }
////        cout<<endl;
////    }
//
//
//    return array[x.size()][y.size()];
//
//
//}
//
//int naive_recursse(string x, string y){
//
//    if(y=="" && x==""){
//        return 0;
//    }
//    else if(x==""){
//        return y.size();
//    }
//    else if(y==""){
//        return x.size();
//    }
//
//    //last character equal part
//    if(x.at(x.size()-1)==y.at(y.size()-1)){
//        return naive_recursse(x.substr(0,x.size()-1), y.substr(0,y.size()-1));
//    }
//
//    //last charcter is not equal
//
//    //replace part
//    int replace=naive_recursse(x.substr(0,x.size()-1),  y.substr(0,y.size()-1));
//    int delete_word=naive_recursse(x.substr(0,x.size()-1), y);
//    int insert_word=naive_recursse(x, y.substr(0,y.size()-1));
//
//    if(replace<delete_word){
//        if(replace<insert_word){
//            return 1+replace;
//        }
//        return 1+insert_word;
//    }
//    else{
//        if(delete_word<insert_word){
//            return 1+delete_word;
//        }
//        else{
//            return 1+insert_word;
//        }
//    }
//
//
//
//
//}
//
//
//int main(){
//
//    string word1="acscnaslcndjacndscklndsjcndcn;dsklcnasdncasklcndjcadjclsac";
//    string word2="ncdspijbcds78chadwcuhewfnc9ewjdfjncdiusclkjsdbadviaubsdvdsjcndlscjknsd89cdjcnldkccdcdlkcd";
//
//    //cout<<word1.substr(0,word1.size()-1)<<endl;
//    cout<<dynamic_distance(word1, word2)<<endl;
//    cout<<naive_recursse(word1, word2)<<endl;
//
//}




//
//
//string naive(string x[], int size, int absolute_pos, int& cost){
//
//    if(size==0){
//        return "";
//    }
//    else{
//        bool line=false;//line has been counted over
//        int count=0;
//        int line_size=0;
//        int min=5000000;
//        int loc=0;
//        string addition="";
//        while(!line){
//            line_size=line_size+ x[count].size()+1;
//            int val=pow((5-line_size),3);
//
//            cout<<x[count]<<" "<<count<<" "<<val<<endl;
//
//            if(line_size>5 || line_size>size){
//                line=true;
//                break;
//            }
//            count++;
//            int sub_size=cost+val;
//            string potential=naive(x+count, size-line_size, absolute_pos+count, sub_size);
//
//
//            if(sub_size<min){
//
//                loc=count;
//                min=sub_size;
//                addition=potential;
//                cost=sub_size;
//            }
//        }
//        string return_string=to_string(absolute_pos+loc)+" "+addition;
//        return return_string;
//    }
//
//}
//
//int main(){
//    string* x= new string[30];
//    string sentence="a a a a a bbbbb";
//
//    string word;
//    int count=0;
//    for(int i=0; i<sentence.size(); i++){
//        if(sentence.at(i)!=' '){
//            word=word+sentence.at(i);
//        }
//        else if(i!=0){
//            if(count==30){
//                break;
//            }
//            x[count]=word;
//            count++;
//            word="";
//        }
//    }
//    x[count]=word;
//
//    int cost=0;
//    cout<<naive(x, sentence.size(),0, cost)<<endl;;
//
//
//
//
//
//
//}

//
//int knapsack(int numbers[], int weights[], int size, int max_size){
//    int** array= new int*[size];
//
//    for(int i=0; i<size; i++){
//        array[i]=new int[max_size];
//    }
//
//    for(int i=0; i<size; i++){
//        array[i][0]=0;
//    }
//
//    for(int r=0; r<size; r++){
//        for(int c=0; c<max_size; c++){
//            array[r][c]=0;
//
//        }
//
//    }
//         //r is the index
//         //c is the remaining size
//
//    for(int c=0; c<max_size; c++){
//        for(int r=size-1; r>=0;r--){
//
//
//            for(int r=0; r<size; r++){
//                for(int c=0; c<max_size; c++){
//                    cout<<array[r][c]<<",";
//
//                }
//                cout<<endl;
//
//            }
//
//            cout<<endl;
//
//            if(r==size-1){
//                if(c-weights[r]>=0){
//                    array[r][c]=weights[r];
//                }
//                else{
//                    array[r][c]=0;
//
//                }
//            }
//            else{
//                if(c-weights[r]>=0){
//                    int addition=array[r+1][c-weights[r]]+numbers[r];
//                    int no_addition=array[r+1][c];
//                    if(addition>no_addition){
//                        array[r][c]=addition;
//                    }
//                    else{
//                        array[r][c]=no_addition;
//                    }
//
//                }
//                else{
//                    array[r][c]=array[r+1][c];
//                }
//            }
//
//
//        }
//
//    }
//
//
//
//
//    return array[1][max_size-1];
//}


//
//int recurse_naive(int profit[], int size[], int indexes_left, int num_sizes, int space_left){
//    int return_val=0;
//
//    if(space_left==0 || indexes_left==0 ){
//        return 0;
//    }
//
//    else{
//
//        int index=0;
//        index=num_sizes-indexes_left;
//            if(space_left-size[index]>=0){
//                int addition = profit[index]+recurse_naive(profit, size, indexes_left-1, num_sizes, space_left-size[index]);
//                int no_addition= recurse_naive(profit, size, indexes_left-1, num_sizes, space_left);
//                if(addition>no_addition){
//                    return_val=addition;
//
//                }
//                else{
//                    return_val=no_addition;
//                }
//            }
//    }
//
//    return return_val;
//
//}
//
//
//int dynamic(int profit[], int size[], int indexes_left, int num_sizes, int space_left, int** dp){
//    int return_val=0;
//
//    if(space_left==0 || indexes_left==0 ){
//        return 0;
//    }
//
//    else{
//
//        if(dp[num_sizes-indexes_left][space_left-1]==-1){
//
//            int index=0;
//            index=num_sizes-indexes_left;
//                if(space_left-size[index]>=0){
//                    int addition = profit[index]+dynamic(profit, size, indexes_left-1, num_sizes, space_left-size[index],dp);
//                    int no_addition= dynamic(profit, size, indexes_left-1, num_sizes, space_left,dp);
//                    if(addition>no_addition){
//                        return_val=addition;
//                        dp[num_sizes-indexes_left][space_left-1]=addition;
//
//                    }
//                    else{
//                        dp[num_sizes-indexes_left][space_left-1]=no_addition;
//                        return_val=no_addition;
//
//                    }
//                }
//        }
//        else{
//            return dp[num_sizes-indexes_left][space_left-1];
//        }
//
//    }
//    return return_val;
//
//}
//
//int main(){
//
//    int profits[10]={1, 6, 10, 16,15,22,21,25,22,12};
//    int sizes[10]={1, 2, 3, 5,1,7,8,1,2,10};
//
//    int** dp=new int*[10];
//    for(int i=0; i<10; i++){
//        dp[i]=new int[30];
//    }
//
//    for(int i=0; i<10; i++){
//        for(int m=0; m<30; m++){
//            dp[i][m]=-1;
//
//        }
//    }
//
//    cout<<dynamic(profits, sizes, 10, 10, 30, dp)<<endl;
//
//
//}


bool recursios(string a, string b, string word, int a_last_index, int b_last_index){

    if(word.size()==0){
        return true;
    }
    bool found_a=false;
    bool found_b=false;
    if(a.at(a_last_index)!=word.at(0) && b.at(b_last_index)!=word.at(0)){
        return false;
    }
    
    else{

        if(a.at(a_last_index)==word.at(0)){
            found_a=recursios(a, b, word.substr(1,word.size()-1), (a_last_index+1)%a.size(), (b_last_index)%b.size());
        }
        if(b.at(b_last_index)==word.at(0)){
            found_b=recursios(a, b, word.substr(1,word.size()-1), (a_last_index)%a.size(), (b_last_index+1)%b.size());
        }
        
        
    }


    return (found_a||found_b);
}


bool dynamic(string a, string b, string word, int a_last_index, int b_last_index, int*** dp){

    if(word.size()==0){
        return true;
    }
    bool found_a=false;
    bool found_b=false;
    if(dp[a_last_index][b_last_index][word.size()-1]==-1){
        if(a.at(a_last_index)!=word.at(0) && b.at(b_last_index)!=word.at(0)){
            dp[a_last_index][b_last_index][word.size()-1]=false;
            return false;
        }
        else{
            if(a.at(a_last_index)==word.at(0)){
                found_a=dynamic(a, b, word.substr(1,word.size()-1), (a_last_index+1)%a.size(), (b_last_index)%b.size(),dp);
            }
            if(b.at(b_last_index)==word.at(0)){
                found_b=dynamic(a, b, word.substr(1,word.size()-1), (a_last_index)%a.size(), (b_last_index+1)%b.size(),dp);
            }
            
            
        }
    }
    else{
        return dp[a_last_index][b_last_index][word.size()-1];
    }
    dp[a_last_index][b_last_index][word.size()-1]=(found_a||found_b);

    return (found_a||found_b);
}

int main(){
    string s="101";
    string m="010";
    
    string inter_leaved="";
    int*** dp=new int**[s.size()];
    
    int count_s=0;
    int count_m=0;
    for(int i=0; i<10; i++){
        if(rand()%2==0){
            inter_leaved=inter_leaved+s.at((count_s%s.size()));
            count_s++;

        }
        else{
            inter_leaved=inter_leaved+m.at(count_m%m.size());
            count_m++;

        }
    }
    
    //inter_leaved="1111111";
    for(int i=0; i<s.size(); i++){
        dp[i]=new int*[m.size()];
        for(int k=0; k<m.size(); k++){
            dp[i][k]=new int[inter_leaved.size()];
        }
    }
    for(int i=0; i<s.size(); i++){
        for(int k=0; k<m.size(); k++){
            for(int l=0; l<inter_leaved.size(); l++){
                dp[i][k][l]=-1;
            }
        }
    }

    
    
//    inter_leaved="inter_leaved";
    inter_leaved="inter_leaved";
   // cout<<recursios(s, m, inter_leaved, 0, 0)<<endl;
    cout<<dynamic(s, m, inter_leaved, 0, 0, dp)<<endl;

    
}
