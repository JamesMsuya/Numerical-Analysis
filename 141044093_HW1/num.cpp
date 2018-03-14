
/* 
 * File:   main.cpp
 * Author: Gme Jemsy
 *
 * Created on March 11, 2017, 9:21 AM
 */

#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cmath>

using namespace std;



double bisect(double a, double b, string str, double e);

double funct(double num,int questNum);



int main(int argc, char** argv){
    
    
    bisect(atof(argv[1]),atof(argv[2]),argv[3],atof(argv[4]));
   
    return 0;
}



double bisect(double a,double b, string str, double e){
	int i=1;
	ofstream out;
    out.open("otput.txt");
    double p_o;
    double p_n;
    double a_1;
    double b_1;
    p_o=a;
    a_1=a;
    b_1=b;
    p_n=(a_1+b_1)/2;
  
    while(!((fabs(p_n-p_o) < e) && (fabs(funct(p_n,4)) < e) && (fabs((p_n-p_o)/p_n) < e ))){
        
          p_o=p_n;
          p_n=(a_1+b_1)/2;
      
       
        if(funct(p_n,4)*funct(a_1,4)>0){     
            a_1 = p_n;

        } 
        if(funct(p_n,4)*funct(a_1,4)<0){     
           b_1 = p_n; 
        } 
        
        
        cout<< i <<"        " << p_n <<"                "<< fabs(p_n-p_o) <<"                " <<fabs((p_n-p_o)/p_n)<<endl;
        out<< i <<"        " << p_n <<"                "<< fabs(p_n-p_o) <<"                " <<fabs((p_n-p_o)/p_n)<<endl;
        if(i>100){     
            cout <<"error "<<endl;
            break;
            exit(1);
        } 
        i++;
    }
    out<<"Root: "<< p_n <<", iterations :" << i<< ", Iterations(Theory) : "<<ceil(log((b-a)/e)/log(2))  <<endl;
    out.close();
    cout<<"Root: "<< p_n <<", iterations :" << i<< ", Iterations(Theory) : "<<ceil(log((b-a)/e)/log(2))  <<endl;

    return(p_n);
}

double funct(double num,int questNum){
    double ans;
    
    if(questNum==1){
        
        ans=(3*num - pow(2.7,num));
        return(ans);
    }
    
    if(questNum==2){
        
       ans = 2 * num + cos(num) * 3 - (pow(2.7,num));
        return(ans);
    }
    
    if(questNum==3){
       
        ans=((num * num) - (4 * num) - log(num));
        return(ans);
    } 
    
   if(questNum==4){
      
       ans = num + 1 - 2*sin(3.14*num);
       return ans;
    }
   return 0; 
}