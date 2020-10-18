#include <bits/stdc++.h>
#include <math.h>
#include <cmath>
#include "rtree.h"
using namespace std;

/*Variables globales*/
RTree *rt;
vector<Nodo*> puntos;
bool ial = false;
bool showCircle=false;
int dim;

Data point2D(float x,float y){
    Data rsp;
    rsp.push_back(x);
    rsp.push_back(y);
    return rsp;
}

int main(){
    cout<<"INIT_PROGRAM"<<endl;
    int m = 2,M = 3;
    dim = 2;
    rt = new RTree(dim,m,M); //(dimension,m,M)
    vData vrd;
	float xmax=-INFINITY,xmin=INFINITY;
	float ymax=-INFINITY,ymin=INFINITY;
    
//    for(int i=0;i<2;i++){
        Data uno = point2D(1,9);
        Data dos = point2D(2,10);
        Data tres = point2D(4,8);
        Data cuatro = point2D(6,7);
        Data cinco = point2D(9,10);
        Data seis = point2D(7,5);
        Data siete = point2D(5,6);
        Data ocho = point2D(4,3);
        Data nueve = point2D(3,2);
        
        vrd.push_back(uno);	
        vrd.push_back(dos);
        vrd.push_back(tres);
        vrd.push_back(cuatro);
        vrd.push_back(cinco);
        vrd.push_back(seis);	
        vrd.push_back(siete);	
//	}


    for(int i=0;i<vrd.size();i++){
        Nodo *tmp = new Nodo(dim,vrd[i]);		
        rt->insert(tmp);
        rt->print();
    }

    rt->getRectangles();

    return 0;
}