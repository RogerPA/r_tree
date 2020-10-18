#include <bits/stdc++.h>
#include "Rrectangle.h"
using namespace std;

//template<typename T,typename L>
class Nodo{
  public:
    bool isLeaf,wasSplit,isData;
    vData tuples; 
    vData I; 
    vector<Nodo*>child; 
    Data rPunto;
    int dim;
    Nodo *parent;
    float areac;
    Nodo(int n_dim,bool leaf=false) { 
        dim = n_dim;
        isLeaf = leaf;
        isData = !leaf;
        parent = NULL;
        areac = 0;
    }
    Nodo(int n_dim,Data dt){ 
        dim = n_dim;
        isLeaf = false;
        isData = true;
        I = makeRectangleFromData(dt);
        rPunto = dt;
        areac = 0;
    }
    
    bool overlap(vData pI){
        for(size_t i=0;i<dim;i++){
            if(I[i][0]>pI[i][0] or pI[i][1]>I[i][1]){
                return false;
            }
        }
        return true;
    }

    int exist(vData pI){
        for(size_t i=0;i<child.size();i++){
            if(child[i]->I==pI){
                return i;
            }
        }
        return -1;
    }

    bool deleteChild(Nodo*H){
        for(size_t i=0;i<child.size();i++)
            if(child[i]==H){
                child.erase(child.begin()+i);
                return true;
            }   
        return false;
    }

    void addEntry(Nodo *E){
        int ext=exist(E->I);
        if(ext==-1){
            child.push_back(E);
            E->parent=this;
            updateRectangleI();
        }
        else{
            //cout<<"Elemento existe"<<endl;
        }     
    }

    bool searchData(Data T){
        if ( std::find(tuples.begin(), tuples.end(), T) != tuples.end() )
            return true;
        return false;
    }
    
    void updateRectangleI(){
        if(child.size()==1){
            I=child[0]->I;
            areac=0;
        }
        else{
            for(size_t i=0;i<child.size();i++){
                for(size_t j=0;j<dim;j++){
                    if(child[i]->I[j][0] < I[j][0] ){
                        I[j][0]=child[i]->I[j][0];
                    }
                    if(child[i]->I[j][1] > I[j][1] ){
                        I[j][1]=child[i]->I[j][1];
                    }
                }
                
            }
            areac = area(I);
        }
    }
    void print(){   
        printVData(I);
    }

};

//template<typename T,typename D>
class RTree{
  public:
    int M,m;
    int dim;
    Nodo *root;
    vector<vData> allRectangles;
    vector<Data> allPoints;
    vector<Data> radioPoints;
    RTree(int n_dim,int n_m,int n_M){
        M=n_M;
        m=n_m;
        root=new Nodo(n_dim,true);
        dim=n_dim;
    }

    bool search(Nodo *&p,vData pI) { //Buscar un punto o un rectangulo
        // S2 /////////
        if(p->isLeaf) {
            return true;
        }
        // S1 ///////
        for(size_t x=0;x<p->child.size();x++) {
            Nodo*currChild=p->child[x];
            if(currChild->overlap(pI)){
                search(currChild,pI);
            }
        }
        return false;
    }

    bool insert(Nodo *E) { // se requiere insertar una nueva entrada E 
        allPoints.push_back(E->rPunto);
        Nodo* L,*LL;
        LL=NULL;       
        // I1 invocamos a chooseLeaf->para seleccionar la hoja
        chooseLeaf(E->I,L);
        
        L->addEntry(E);
        if(L->child.size() > M) {
            splitNode(L,LL);            
        }
        adjustTree(L,LL); // SI no hubo split LL==NULL
        if(L==root and LL!=NULL) {
            Nodo* tempRoot=new Nodo(dim,false);
            tempRoot->addEntry(L);
            tempRoot->addEntry(LL);
            root=tempRoot;
        }
        return true;
    }

    void chooseLeaf(vData E,Nodo *&N) { // selecciona la hoja donde entra E
        //  inicializar como raiz
        N = root;
        //   si estamos en hoja devuelve N
        while(!N->isLeaf){
            float tempArea = INFINITY;
            Nodo *TN;
            int childSize=N->child.size();
            for(size_t i=0;i<childSize;i++) {
                float nTempArea;
                makeRectangle(E,N->child[i]->I,nTempArea);
                //  /////
                nTempArea=nTempArea- N->child[i]->areac;
                if(tempArea>nTempArea){
                    TN=N->child[i];
                    tempArea=nTempArea;
                }
            }
            //  repetimos en 
            N=TN;
        }
        return;
    }

    void adjustTree(Nodo* &N,Nodo*&NN){ // expandimos
        while(N!=root){
            Nodo *P=N->parent;
            P->isLeaf=false;
            Nodo *PP=NULL;
            if(NN!=NULL){
                P->addEntry(NN);
                if(P->child.size()>M){
                    splitNode(P,PP);
                }
            }
            P->updateRectangleI();
            N = P;
            NN = PP;
            
            
        }
        root->updateRectangleI();
        return;
    }

    //Split cuadratico 
    void splitNode(Nodo* &G1,Nodo* &G2){ 
       // nuestras semilaas para el pickseeds 
        Nodo *E1,*E2;
        float Ed1,Ed2;
        vector<Nodo*>LP;
        LP=G1->child;
        bool checkLeaf=G1->isLeaf;
        pickSeeds(LP,E1,E2);
        G1->child.clear();
        G2 = new Nodo(dim,checkLeaf); //Nodos como grupos G1, G2;
        G1->addEntry(E1);
        G2->addEntry(E2); 
        while(LP.size()>0){   
            if( (G1->child.size()+LP.size())==m ){
                G1->child.insert(G1->child.end(), LP.begin(), LP.end());
                G1->updateRectangleI();
                LP.resize(0);
            }
            else if( (G2->child.size()+LP.size())==m ){
                G2->child.insert(G2->child.end(), LP.begin(), LP.end());
                G2->updateRectangleI();
                LP.resize(0);
            }
            if(LP.size()>0)
                pickNext(LP,G1,G2);
        }   
        return;
    }

    void pickSeeds(vector<Nodo*>&LP,Nodo* &E1,Nodo* &E2){
        float d=-INFINITY;
        int indxE1,indxE2;
        // Pentradas para los grupos
        for(int x1=0;x1<LP.size();x1++){
            for(int x2=0;x2<LP.size();x2++){
                if(x1!=x2){
                    float areaJ;
                    vData J=makeRectangle(LP[x1]->I,LP[x2]->I,areaJ);
                    float a1 = LP[x1]->areac;
                    float a2 = LP[x2]->areac;
                    //float td=areaJ;
                    float td = areaJ-a1-a2;//areaJ*a1*a2;
                    if(td >= d){
                        d = td;
                        indxE1 = x1;
						//cout<<indxE1<<endl;
                        indxE2 = x2;
                    }
                }
            }
        }
        E1 = LP[indxE1];
        E2 = LP[indxE2];
                   
        /*eliminamos elementos del grupo*/
        if(indxE1>indxE2){
            swap(indxE1,indxE2);
        }
		//cout<<indxE1<<endl;
        LP.erase(LP.begin()+indxE2);
        LP.erase(LP.begin()+indxE1); 
    }

    void pickNext(vector<Nodo*>&LP,Nodo* &G1,Nodo* &G2){    
        int indxE1,indxE2;
        float dG1=INFINITY,dG2=INFINITY;
        //vData G1,G2;
        float aG1=G1->areac;
        float aG2=G2->areac;
        vData ntG1,ntG2;
        // para cada entrada area
        for(size_t i=0;i<LP.size();i++){
            float areaG1,areaG2;
            vData tG1 = makeRectangle(G1->I,LP[i]->I,areaG1);
            vData tG2 = makeRectangle(G2->I,LP[i]->I,areaG2);
            float tdG1 = areaG1-aG1;
            float tdG2 = areaG2-aG2;
            if(tdG1 < dG1){
                indxE1 = i; 
                dG1 = tdG1; 
				//cout<<dG1<<tdG1<<endl;
                ntG1 = tG1; //nuevo Rectangulo ->  I
            }
            if(tdG2 < dG2){
                indxE2 = i;
                dG2 = tdG2;
				//cout<<dG1<<tdG1<<endl;
                ntG2=tG2;           
            }
            
        }
        // estrada con maxima diferencia
        if(dG1<dG2){
            G1->addEntry(LP[indxE1]);
            LP.erase(LP.begin()+indxE1);
        }
        else{   
            G2->addEntry(LP[indxE2]);
            LP.erase(LP.begin()+indxE2);
        }
    }

	/****************************************************************
	 *                  VISUALIZAR DATOS DEL ARBOL
	 *****************************************************************
	 */

    void print(){
        cout<<"========================================="<<endl;
        cout<<"      Imprimiendo Arbol  "<<endl;
        cout<<"========================================="<<endl;
        printR(root);    
    }

    void printR(Nodo *P){
        cout<<"_________";
        cout<<" Nodo ";      
        cout<<"_________"<<endl;
        P->print();
        if(P->parent==NULL)
            cout<<" - Padre null - ";
        else{
            cout<<" - Padre -";
            P->parent->print();
        }
        cout<<endl;
        if(P->isLeaf){
            cout<<"----------";
            cout<<" Hoja ";
            cout<<"----------"<<endl;
            for(int i=0;i<P->child.size();i++){
                cout<<" --> ";
                P->child[i]->print();
            }
            return;
        }
        
        for(int i=0;i<P->child.size();i++){
            printR(P->child[i]);
        }
    }
    void getRectangles(){
        allRectangles.clear();
        getRectanglesR(root);
    }

    void getRectanglesR(Nodo *P){
        for(size_t i=0;i<P->child.size();i++){
            allRectangles.push_back(P->I);
            getRectanglesR(P->child[i]);
        }
        if(P->isLeaf){
            return;
        }
    }


};
