#include <iostream>
#include "Boundary.h"

using namespace std;

int main(){
    Boundary bc1 = Boundary();
    if(bc1.kind() == BC::DBC)
        cout<<"Default OK!"<<endl;
    cout<<bc1.alpha()<<bc1.beta()<<bc1.gamma()<<endl;

    Boundary bc2 = Boundary(1, 0, 0);
    if(bc2.kind() == BC::NBC)
        cout<<"NBC OK!"<<endl;
    cout<<bc2.alpha()<<bc2.beta()<<bc2.gamma()<<endl;

    Boundary bc3 = Boundary(1, 2, 0);
    if(bc3.kind() == BC::RBC)
        cout<<"RBC OK!"<<endl;
    cout<<bc3.alpha()<<bc3.beta()<<bc3.gamma()<<endl;

    Boundary bc4 = Boundary(0, 1, 0);
    if(bc4.kind() == BC::DBC)
        cout<<"DBC OK!"<<endl;
    cout<<bc4.alpha()<<bc4.beta()<<bc4.gamma()<<endl;

    return 0;
}
