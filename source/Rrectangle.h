#include <bits/stdc++.h>
using namespace std;
typedef vector<float> Data;
typedef vector<Data> vData;

float area(vData mm){
    float areaCalc=1;
    for(size_t i=0;i<mm.size();i++){
        areaCalc=areaCalc*abs(mm[i][1]-mm[i][0]);
    }
    return areaCalc;
}

float distP2P(Data P1,Data P2){
    float dist=0;
    for(size_t i=0;i<P1.size();i++){
        dist+=abs(P1[i]-P2[i]);
    }
    return dist;
}


vData makeRectangle(vData E1,vData E2,float &area){
    vData R;
    float areaT=1;
    int E1size=E1.size();
    for(size_t i=0;i<E1size;i++){
        float tempMax,tempMin;
        if(E1[i][0]>E2[i][0])
            tempMin=E2[i][0];
        else
            tempMin=E1[i][0];

        if(E1[i][1]<E2[i][1])
            tempMax=E2[i][1];
        else
            tempMax=E1[i][1];
        R.push_back({tempMin,tempMax});
        areaT=areaT*abs(tempMax-tempMin);
    }  
    return R;
}

vData makeRectangleFromData(Data E1){
    vData R;
    for(size_t i=0;i<E1.size();i++){
        vector<float>tempDim(2);
        tempDim[0]=tempDim[1]=E1[i];
        R.push_back(tempDim);
    }
    return R;
}


void printData(Data d){
    for(size_t i=0;i<d.size();i++){
        cout<<d[i]<<", ";
    }
    cout<<endl;
}

//************************************************
//************************************************
void printVData(vData d){
    for(int i=0;i<d.size();i++){
        cout<<"["<<d[i][0]<<","<<d[i][1]<<"] ";
    }
    cout<<endl;
}
