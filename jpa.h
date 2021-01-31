#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/PxPyPzE4D.h"

namespace HHbbWWJPA {
  // Wjj_simple: sum of p4
  // dR_HadW_bjet = lambda bP4,j1P4,j2P4 : op.deltaR(self.Wjj_simple(j1P4,j2P4), bP4)
  //
  // motherPx = lambda ob1p4,ob2p4 : (ob1p4.Px()+ob2p4.Px())
  // motherPy = lambda ob1p4,ob2p4 : (ob1p4.Py()+ob2p4.Py())
  // motherPz = lambda ob1p4,ob2p4 : (ob1p4.Pz()+ob2p4.Pz())
  // motherE  = lambda ob1p4,ob2p4 : (ob1p4.E()+ob2p4.E())
  // betaX    = lambda ob1p4,ob2p4 : motherPx(ob1p4,ob2p4)/motherE(ob1p4,ob2p4)
  // betaY    = lambda ob1p4,ob2p4 : motherPy(ob1p4,ob2p4)/motherE(ob1p4,ob2p4)
  // betaZ    = lambda ob1p4,ob2p4 : motherPz(ob1p4,ob2p4)/motherE(ob1p4,ob2p4)
  // beta2    = lambda ob1p4,ob2p4 : op.pow(betaX(ob1p4,ob2p4),2)+op.pow(betaY(ob1p4,ob2p4),2)+op.pow(betaZ(ob1p4,ob2p4),2)
  // betap    = lambda ob1p4,ob2p4 : betaX(ob1p4,ob2p4)*motherPx(ob1p4,ob2p4) + betaY(ob1p4,ob2p4)*motherPy(ob1p4,ob2p4) + betaZ(ob1p4,ob2p4)*motherPz(ob1p4,ob2p4)
  // gamma2   = lambda ob1p4,ob2p4 : op.switch(beta2(ob1p4,ob2p4) > 0, (gamma(ob1p4,ob2p4)-1)/beta2(ob1p4,ob2p4) , op.c_float(0.0))
  // boostPx  = lambda ob1p4,ob2p4 : ob1p4.Px() + gamma2(ob1p4,ob2p4)*betap(ob1p4,ob2p4)*betaX(ob1p4,ob2p4) + gamma(ob1p4,ob2p4)*betaX(ob1p4,ob2p4)*ob1p4.E()
  // boostPy  = lambda ob1p4,ob2p4 : ob1p4.Px() + gamma2(ob1p4,ob2p4)*betap(ob1p4,ob2p4)*betaY(ob1p4,ob2p4) + gamma(ob1p4,ob2p4)*betaY(ob1p4,ob2p4)*ob1p4.E()
  // boostPz  = lambda ob1p4,ob2p4 : ob1p4.Pz() + gamma2(ob1p4,ob2p4)*betap(ob1p4,ob2p4)*betaZ(ob1p4,ob2p4) + gamma(ob1p4,ob2p4)*betaZ(ob1p4,ob2p4)*ob1p4.E()
  // boostP   = lambda ob1p4,ob2p4 : op.sqrt(op.pow(boostPx(ob1p4,ob2p4),2) + op.pow(boostPy(ob1p4,ob2p4),2) + op.pow(boostPz(ob1p4,ob2p4),2))
  // comp_cosThetaS: = lambda ob1p4,ob2p4 : op.abs(boostPz(ob1p4,ob2p4)/boostP(ob1p4,ob2p4))
  //
  float cosThetaS(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> pA, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> pB) {
    const auto mother = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>(pA + pB);
    const auto mBeta = mother.BoostToCM();
    const auto beta2 = mBeta.Mag2();
    const auto betap = mBeta.Dot(mother.Vect());
    const auto gamma = 1./std::sqrt(1.-beta2);
    const auto gamma2 = ( beta2 > 0. ) ? (gamma-1.)/beta2 : 0.;
    const auto boostP = pA.Vect() + (gamma2*betap + gamma*pA.E())*mBeta;
    return std::abs(boostP.Z()/std::sqrt(boostP.Mag2()));
  }

  // ax = lambda visP4, met : 125.18*125.18 - op.pow(visP4.M(), 2) + 2.*visP4.Px()*met.p4.Px() + 2.*visP4.Py()*met.p4.Py()
  // A  = lambda visP4 : 4.0 * op.pow(visP4.E(), 2) - op.pow(visP4.Pz(), 2)
  // B  = lambda visP4, met : -4.0 * ax(visP4, met) * visP4.Pz()
  // C  = lambda visP4, met : 4.0*op.pow(visP4.E(), 2)*(op.pow(met.p4.Px(), 2) + op.pow(met.p4.Py(), 2)) - op.pow(ax(visP4, met), 2)
  // D  = lambda visP4, met : (op.pow(B(visP4, met), 2) - 4.0*A(visP4)*C(visP4, met))
  // pos   = lambda visP4, met : (-B(visP4, met) + op.sqrt(D(visP4, met)))/(2.*A(visP4))
  // neg   = lambda visP4, met : (-B(visP4, met) - op.sqrt(D(visP4, met)))/(2.*A(visP4))
  // neuPz = lambda visP4, met : (op.switch(D(visP4, met) < 0., -B(visP4, met)/(2*A(visP4)), op.switch(op.abs(pos(visP4, met)) < op.abs(neg(visP4, met)), pos(visP4, met), neg(visP4, met))))
  // self.neuP4 = lambda visP4, met : op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", (met.p4.Px(), met.p4.Py(), neuPz(visP4, met), op.sqrt(op.pow(met.p4.Px(),2)+op.pow(met.p4.Py(),2)+op.pow(neuPz(visP4, met),2)))).result
  // self.Wlep_simple = lambda wj1P4,wj2P4,lepP4,met : lepP4 + self.neuP4(wj1P4+wj2P4+lepP4, met)
  // self.Wjj_simple  = lambda j1P4,j2P4 : j1P4 + j2P4
  // self.HWW_simple  = lambda wj1P4,wj2P4,lepP4,met : self.Wjj_simple(wj1P4,wj2P4) + self.Wlep_simple(wj1P4,wj2P4,lepP4,met)
  
  float HWW_SimpleMassWithRecoNeutrino(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> wj1P4, 
				       const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> wj2P4,
				       const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> lepP4,
				       const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> metP4) {
    const auto visP4 = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>(wj1P4+wj2P4+lepP4);
    double a = std::pow(125.18,2) - visP4.M()*visP4.M() + 2.0*visP4.Px()*metP4.Px() + 2.0*visP4.Py()*metP4.Py();
    double A = 4.0*(std::pow(visP4.E(),2) - std::pow(visP4.Pz(),2));
    double B = -4.0*a*visP4.Pz();
    double C = 4.0*std::pow(visP4.E(),2)*(std::pow(metP4.Px(),2) + std::pow(metP4.Py(),2)) - a*a;
    double delta = B*B - 4.*A*C;

    double metPz = 0.;
    if ( delta < 0. ) {
      metPz = -B/(2.*A);
    } else {
      double pos = (-B + std::sqrt(delta))/(2.*A);
      double neg = (-B - std::sqrt(delta))/(2.*A);
      if ( std::fabs(pos) < std::fabs(neg) ) metPz = pos;
      else metPz = neg;
    }
    double neuPx = metP4.Px();
    double neuPy = metP4.Py();
    double neuPz = metPz;
    double neuEn = std::sqrt(neuPx*neuPx + neuPy*neuPy + neuPz*neuPz); 
    auto neuP4 = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >(neuPx, neuPy, neuPz, neuEn);
    auto HWWP4 = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>(wj1P4+wj2P4+lepP4+neuP4);
    return HWWP4.M();
  }
};
