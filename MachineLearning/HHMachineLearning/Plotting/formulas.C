#include "TLorentzVector.h"

double LorentzVector2M(float px1,float py1,float pz1,float E1,float px2,float py2,float pz2,float E2){
    auto a = TLorentzVector();
    auto b = TLorentzVector();
    a.SetPxPyPzE(px1,py1,pz1,E1);
    b.SetPxPyPzE(px2,py2,pz2,E2);
    return (a+b).M();
}
double LorentzVector2E(float px1,float py1,float pz1,float E1,float px2,float py2,float pz2,float E2){
    auto a = TLorentzVector();
    auto b = TLorentzVector();
    a.SetPxPyPzE(px1,py1,pz1,E1);
    b.SetPxPyPzE(px2,py2,pz2,E2);
    return (a+b).E();
}
double LorentzVector2Pt(float px1,float py1,float pz1,float E1,float px2,float py2,float pz2,float E2){
    auto a = TLorentzVector();
    auto b = TLorentzVector();
    a.SetPxPyPzE(px1,py1,pz1,E1);
    b.SetPxPyPzE(px2,py2,pz2,E2);
    return (a+b).Pt();
}
double LorentzVector3M(float px1,float py1,float pz1,float E1,float px2,float py2,float pz2,float E2,float px3,float py3,float pz3,float E3){
    auto a = TLorentzVector();
    auto b = TLorentzVector();
    auto c = TLorentzVector();
    a.SetPxPyPzE(px1,py1,pz1,E1);
    b.SetPxPyPzE(px2,py2,pz2,E2);
    c.SetPxPyPzE(px3,py3,pz3,E3);
    return (a+b+c).M();
}
double LorentzVector3E(float px1,float py1,float pz1,float E1,float px2,float py2,float pz2,float E2,float px3,float py3,float pz3,float E3){
    auto a = TLorentzVector();
    auto b = TLorentzVector();
    auto c = TLorentzVector();
    a.SetPxPyPzE(px1,py1,pz1,E1);
    b.SetPxPyPzE(px2,py2,pz2,E2);
    c.SetPxPyPzE(px3,py3,pz3,E3);
    return (a+b+c).E();

}double LorentzVector3Pt(float px1,float py1,float pz1,float E1,float px2,float py2,float pz2,float E2,float px3,float py3,float pz3,float E3){
    auto a = TLorentzVector();
    auto b = TLorentzVector();
    auto c = TLorentzVector();
    a.SetPxPyPzE(px1,py1,pz1,E1);
    b.SetPxPyPzE(px2,py2,pz2,E2);
    c.SetPxPyPzE(px3,py3,pz3,E3);
    return (a+b+c).Pt();
}
