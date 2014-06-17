#include <iostream>
#include "Boundary.h"

using namespace std;

int main(){
    Boundary bc1 = Boundary();  // default: PBC
    if(bc1.kind() == BC::PBC)
        cout<<"Default OK!"<<endl;
    cout<<bc1.alpha()<<bc1.beta()<<bc1.gamma()<<endl;

    Boundary bc2 = Boundary(1, 0, 0);  // NBC
    if(bc2.kind() == BC::NBC)
        cout<<"NBC OK!"<<endl;
    cout<<bc2.alpha()<<bc2.beta()<<bc2.gamma()<<endl;

    Boundary bc3 = Boundary(1, 2, 0);  // RBC
    if(bc3.kind() == BC::RBC)
        cout<<"RBC OK!"<<endl;
    cout<<bc3.alpha()<<bc3.beta()<<bc3.gamma()<<endl;

    Boundary bc4 = Boundary(0, 1, 0);  // DBC
    if(bc4.kind() == BC::DBC)
        cout<<"DBC OK!"<<endl;
    cout<<bc4.alpha()<<bc4.beta()<<bc4.gamma()<<endl;

    Boundary bc5 = Boundary(0, 0, 0);  // PBC
    if(bc5.kind() == BC::PBC)
        cout<<"PBC OK!"<<endl;
    cout<<bc5.alpha()<<bc5.beta()<<bc5.gamma()<<endl;

    // Copy constructor
    Boundary bc6(bc3);  // RBC
    if(bc6.kind() == BC::RBC)
        cout<<"Copy constructor OK!"<<endl;
    cout<<bc6.alpha()<<bc6.beta()<<bc6.gamma()<<endl;

    // Asignment
    Boundary bc7;
    bc7 = bc3;
    if(bc7.kind() == BC::RBC)
        cout<<"Asignment OK!"<<endl;
    cout<<bc7.alpha()<<bc7.beta()<<bc7.gamma()<<endl;

    return 0;
}
