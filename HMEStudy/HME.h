#include <numeric>
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/PxPyPzE4D.h"
#include "TVector2.h"

#include "TRandom3.h" 
#include "bamboorandom.h" 

class Sampler {
    public:
        std::vector<float> xval;
        std::vector<float> yval;
        std::vector<float> cdf;
        TRandom3 rg;

        void set_xval(std::vector<float> _xval){
            xval = _xval;
        }

        void set_yval(std::vector<float> _yval){
            yval = _yval;
        }

        void computeSampler(){
            float contsum = 0.;
            for (auto y : yval){
                contsum += y;
                cdf.push_back(contsum);
            }
        }
        Sampler(TRandom3 _rg){
            rg = _rg;
        }
        float getRandom(){
            float y = rg.Uniform(0,cdf.back());
            float x = 0;
            for (int i = 0; i < cdf.size() ; i++){
                if (cdf[i] <= y && y <= cdf[i+1]){
                    x = xval[i] + (y-cdf[i]) * (xval[i+1] - xval[i]) / (cdf[i+1]-cdf[i]);
                    break;
                }
            }
            return x;
        }
};

class OnShellW: public Sampler {
    public:
        OnShellW(TRandom3 _rg) : Sampler(_rg) {
            set_xval({40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5, 50.5, 51.5, 52.5, 53.5, 54.5, 55.5, 56.5, 57.5, 58.5, 59.5, 60.5, 61.5, 62.5, 63.5, 64.5, 65.5, 66.5, 67.5, 68.5, 69.5, 70.5, 71.5, 72.5, 73.5, 74.5, 75.5, 76.5, 77.5, 78.5, 79.5, 80.5, 81.5, 82.5, 83.5, 84.5, 85.5, 86.5, 87.5, 88.5, 89.5, 90.5, 91.5, 92.5, 93.5, 94.5, 95.5, 96.5, 97.5, 98.5, 99.5});
            set_yval({784.0, 896.0, 861.0, 1050.0, 1036.0, 1099.0, 1127.0, 1246.0, 1491.0, 1547.0, 1806.0, 1729.0, 2170.0, 2177.0, 2576.0, 2982.0, 3038.0, 3773.0, 3976.0, 4522.0, 4725.0, 5705.0, 6027.0, 6405.0, 6622.0, 7077.0, 6958.0, 8134.0, 8302.0, 9492.0, 9842.0, 11312.0, 12957.0, 16044.0, 19208.0, 25683.0, 35189.0, 54467.0, 100597.0, 217462.0, 308560.0, 152964.0, 58289.0, 26145.0, 14161.0, 8498.0, 5341.0, 3801.0, 2156.0, 1547.0, 1015.0, 889.0, 651.0, 518.0, 343.0, 273.0, 147.0, 84.0, 91.0, 56.0});
            computeSampler();
        }
};

class BjetRescale: public Sampler {
    public:
        BjetRescale(TRandom3 _rg) : Sampler(_rg) {
            set_xval({0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52, 1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78, 1.8, 1.82, 1.84, 1.86, 1.88, 1.9, 1.92, 1.94, 1.96, 1.98, 2.0, 2.02, 2.04, 2.06, 2.08, 2.1, 2.12, 2.14, 2.16, 2.18, 2.2, 2.22, 2.24, 2.26, 2.28, 2.3, 2.32, 2.34, 2.36, 2.38, 2.4, 2.42, 2.44, 2.46, 2.48, 2.5, 2.52, 2.54, 2.56, 2.58, 2.6, 2.62, 2.64, 2.66, 2.68, 2.7, 2.72, 2.74, 2.76, 2.78, 2.8, 2.82, 2.84, 2.86, 2.88, 2.9, 2.92, 2.94, 2.96, 2.98, 3.0, 3.02, 3.04, 3.06, 3.08, 3.1, 3.12, 3.14, 3.16, 3.18, 3.2, 3.22, 3.24, 3.26, 3.28, 3.3, 3.32, 3.34, 3.36, 3.38, 3.4, 3.42, 3.44, 3.46, 3.48, 3.5, 3.52, 3.54, 3.56, 3.58, 3.6, 3.62, 3.64, 3.66, 3.68, 3.7, 3.72, 3.74, 3.76, 3.78, 3.8, 3.82, 3.84, 3.86, 3.88, 3.9, 3.92, 3.94, 3.96, 3.98, 4.0, 4.02, 4.04, 4.06, 4.08, 4.1, 4.12, 4.14, 4.16, 4.18, 4.2, 4.22, 4.24, 4.26, 4.28, 4.3, 4.32, 4.34, 4.36, 4.38, 4.4, 4.42, 4.44, 4.46, 4.48, 4.5, 4.52, 4.54, 4.56, 4.58, 4.6, 4.62, 4.64, 4.66, 4.68, 4.7, 4.72, 4.74, 4.76, 4.78, 4.8, 4.82, 4.84, 4.86, 4.88, 4.9, 4.92, 4.94, 4.96, 4.98, 5.0, 5.02, 5.04, 5.06, 5.08, 5.1, 5.12, 5.14, 5.16, 5.18, 5.2, 5.22, 5.24, 5.26, 5.28, 5.3, 5.32, 5.34, 5.36, 5.38, 5.4, 5.42, 5.44, 5.46, 5.48, 5.5, 5.52, 5.54, 5.56, 5.58, 5.6, 5.62, 5.64, 5.66, 5.68, 5.7, 5.72, 5.74, 5.76, 5.78, 5.8, 5.82, 5.84, 5.86, 5.88, 5.9, 5.92, 5.94, 5.96, 5.98, 6.0});
            set_yval({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 4.0, 6.0, 7.0, 4.0, 4.0, 4.0, 9.0, 6.0, 16.0, 8.0, 7.0, 8.0, 5.0, 6.0, 5.0, 4.0, 8.0, 14.0, 7.0, 21.0, 9.0, 7.0, 14.0, 15.0, 16.0, 9.0, 19.0, 17.0, 28.0, 24.0, 40.0, 51.0, 58.0, 73.0, 88.0, 126.0, 173.0, 269.0, 371.0, 474.0, 594.0, 695.0, 702.0, 777.0, 735.0, 742.0, 636.0, 593.0, 467.0, 458.0, 392.0, 383.0, 341.0, 319.0, 293.0, 270.0, 239.0, 204.0, 184.0, 154.0, 151.0, 153.0, 133.0, 127.0, 101.0, 104.0, 120.0, 77.0, 70.0, 61.0, 57.0, 74.0, 57.0, 73.0, 59.0, 56.0, 47.0, 30.0, 24.0, 38.0, 46.0, 33.0, 32.0, 21.0, 29.0, 30.0, 21.0, 18.0, 25.0, 20.0, 17.0, 19.0, 6.0, 11.0, 14.0, 14.0, 9.0, 12.0, 4.0, 10.0, 11.0, 7.0, 5.0, 7.0, 4.0, 5.0, 4.0, 8.0, 3.0, 2.0, 0.0, 2.0, 8.0, 6.0, 5.0, 0.0, 2.0, 2.0, 6.0, 2.0, 1.0, 1.0, 1.0, 0.0, 2.0, 4.0, 0.0, 1.0, 2.0, 0.0, 2.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
            computeSampler();
        }

};



ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> nuFromOnshellW(
                    float nu_eta, 
                     float nu_phi,
                     const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> l,
                     float mw){
    float deta = nu_eta - l.Eta();
    float dphi = nu_phi - l.Phi();
    float nu_pt = mw*mw/(2*l.Pt()*(cosh(deta)-cos(dphi)));
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> nu;
    if (isinf(nu_pt) || isnan(nu_pt))
        nu.SetCoordinates(0.,0.,0.,0.);
    else
        nu.SetCoordinates(nu_pt,nu_eta,nu_phi,0.);
    return nu;
}

ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> nuFromOffshellW(
                     const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> l1,
                     const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> l2,
                     const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> nu1,
                     const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> met,
                     int control,
                     float mh){
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> nu2;
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> tmplorentz = l1+l2+nu1;
    float nu_tmp_px = met.Px()-nu1.Px();
    float nu_tmp_py = met.Py()-nu1.Px();
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> tmp2lorentz;

    tmp2lorentz.SetPxPyPzE(sqrt(pow(tmplorentz.Pt(),2)+pow(tmplorentz.M(),2)),0,tmplorentz.Pz(),tmplorentz.E());

    TVector2 nu_pxpy(nu_tmp_px, nu_tmp_py);
    float nu_tmp_pt = nu_pxpy.Mod();

    float chdeta = (pow(mh,2)+2*(nu_pxpy.Px()*tmplorentz.Px()+nu_pxpy.Py()*tmplorentz.Py())-pow(tmplorentz.M(),2))/(2*tmp2lorentz.Pt()*nu_tmp_pt);//cosh(nu_offshellW_eta-tmp2lorentz_eta)

    if (chdeta < 1.0) {     
        //delete tmplorentz;     delete tmp2lorentz;     nu2lorentz.SetCoordinates(0, 0, 0, 0);     return false;   }
        nu2.SetCoordinates(0.,0.,0.,0.);
        return nu2;
    }
    float nu_tmp_phi = nu_pxpy.Phi_mpi_pi(nu_pxpy.Phi());
    float deltaeta = acosh(chdeta);
    float nu_tmp_eta = (control == 1) ? (tmp2lorentz.Eta()-deltaeta) : (tmp2lorentz.Eta()+deltaeta);

     if (fabs(nu_tmp_eta) > 7) {
        nu2.SetCoordinates(0.,0.,0.,0.);
        return nu2;
    }
    
    nu2.SetCoordinates(nu_tmp_pt, nu_tmp_eta, nu_tmp_phi, 0);
    return nu2;
    
}

std::pair<float,float> bjetCorrectionsSimple(
    const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> b1,
    const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> b2,
    float mh){

    float c = mh/(b1+b2).M();
    return std::make_pair(c,c);
}

std::pair<float,float> bjetCorrectionsSliding(
    const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> b1,
    const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> b2,
    Sampler gen,
    float mh){
    
    float c1 = 0;
    float c2 = 0;
    
    /* Loop to find solution until found */
    while(true){
        /* Get c1 from pdf */
        c1 = gen.getRandom();
        /* Get c2 from mbb = mh */
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> leadingbjet;
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> trailingbjet;
        if (b1.Pt() > b2.Pt()){
            leadingbjet  = b1;
            trailingbjet = b2;
        }
        else{
            leadingbjet  = b2;
            trailingbjet = b1;
        }
        float x1 = trailingbjet.M2();
        float x2 = 2*c1*(leadingbjet.Dot(trailingbjet));
        float x3 = c1*c1*leadingbjet.M2() - mh*mh;
        /* 
            [c1*P1 + c2*P2]^2 = mH^2
            c2^2 * (m2^2) + c2 * (2*c1*P1.P2) + c1*m1 - mh^2 = 0
            x1 = m2^2
            x2 = 2*c1*P1.P2
            x3 = c1^2 * m1^2 - mh^2
        */
        
        /* Check solutions */
        if (x2 < 0.)
            // P1.P2 > 0 because positive in CM frame and lorentz invariant
            continue; 
        
        if ((x2*x2 - 4*x1*x3) < 0 or x1 == 0)
            // To be solution of the second order polynom
            continue;

        c2 = (-x2 + sqrt(x2*x2 - 4*x1*x3))/(2*x1);

        if (c2 < 0.)
            continue;
    
        if (b1.Pt() > b2.Pt())
            return std::make_pair(c1,c2);
        else
            return std::make_pair(c2,c1);
   }
}




std::pair<float,float> runHME(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> l1,
             const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> l2,
             const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> b1,
             const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> b2,
             const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> met,
             const std::uint32_t eventnr){

    /* Initialize */
    //auto& rg = rdfhelpers::getTRandom3(eventnr);
    TRandom3 rg = rdfhelpers::getTRandom3(eventnr);
    double met_sigma = 25.2;
    int iterations = 10000;
    int converged = 0;
    
    float xmin = 0.;
    float xmax = 5000.;
    float xdelta = 1.;
    std::vector<float> bins (int((xmax-xmin)/xdelta),0.);

    /* PDF generator */
    OnShellW mw_gen(rg);
    BjetRescale rescale_gen(rg);

    /* Iteration loop */
    int it = 0;
    while (it < iterations){
        it += 1;

        /* Random draws */
        float eta_gen = rg.Uniform(-6, 6);
        float phi_gen = rg.Uniform(-3.1415926, 3.1415926);
        float mh = rg.Gaus(125.03, 0.004);
        float mw = mw_gen.getRandom();
        float met_dpx = rg.Gaus(0.0, met_sigma);
        float met_dpy = rg.Gaus(0.0, met_sigma);
//        std::cout<<"-----------------------------------"<<std::endl;
//        std::cout <<"Iteration "<<it<<" - eta = "<<eta_gen<<" - phi = "<<phi_gen<<std::endl;

        /* MW onshell */

        /* bjets corrections */

//        std::cout<<"  Before corrections"<<std::endl;
//        std::cout<<"    MET : "<<met<<std::endl;
//        std::cout<<"    b1  : "<<b1<<std::endl;
//        std::cout<<"    b2  : "<<b2<<std::endl;

        std::pair<float,float> c = bjetCorrectionsSliding(b1,b2,rescale_gen,mh);
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> tmp_b1;
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> tmp_b2;
        tmp_b1 = b1 * c.first;
        tmp_b2 = b2 * c.second;

        /* MET corrections */
        float met_Px_corr = - (c.first-1) * b1.Px() - (c.second-1) * b2.Px();
        float met_Py_corr = - (c.first-1) * b1.Py() - (c.second-1) * b2.Py();
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> tmp_met;
        tmp_met.SetPxPyPzE(met.Px() + met_dpx + met_Px_corr,
                           met.Px() + met_dpx + met_Px_corr,
                           met.Pz(),
                           met.E());
//        std::cout<<"  Corrections : c1 = "<<c.first<<" - c2 = "<<c.second<<std::endl;
//        std::cout<<"  After corrections"<<std::endl;
//        std::cout<<"    MET : "<<tmp_met<<std::endl;
//        std::cout<<"    b1  : "<<tmp_b1<<std::endl;
//        std::cout<<"    b2  : "<<tmp_b2<<std::endl;

        /* H->WW */
        float hme = 0.;
        int solutions = 0; // number of valid solutions

        for (int j = 0; j < 4; j++){
//            std::cout<<"  Solution "<<j<<"/4"<<std::endl;
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> l_offshell;
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> l_onshell;
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> nu_offshell;
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> nu_onshell;

            int is_onshell = j/2;
            /* On-shell */
            if (is_onshell == 0){
                l_onshell = l1;
                l_offshell = l2;
            }
            else{
                l_onshell = l2;
                l_offshell = l1;
            }
            nu_onshell = nuFromOnshellW(eta_gen,phi_gen,l_onshell,mw);
            if (nu_onshell.Pt() == 0. && nu_onshell.Eta() == 0. && nu_onshell.Phi() == 0.){
//                std::cout<<"\t-> No solution for on-shell"<<std::endl;
                continue;   // no solution to on shell
            }

            /* Off-shell */
            int is_offshell = j%2;
            nu_offshell = nuFromOffshellW(l1,l2,nu_onshell,tmp_met,is_offshell,mh);

//            std::cout<<"\tOn-shell  lepton  : "<<l_onshell<<std::endl;
//            std::cout<<"\tOn-shell  neutrino: "<<nu_onshell<<std::endl;
//            std::cout<<"\tOff-shell lepton  : "<<l_offshell<<std::endl;
//            std::cout<<"\tOff-shell neutrino: "<<nu_offshell<<std::endl;

            /* Validity checks */
            if (nu_offshell.Pt() == 0. && nu_offshell.Eta() == 0. && nu_offshell.Phi() == 0.){
//                std::cout<<"\t-> No solution for off-shell"<<std::endl;
                continue;   // no solution to off shell
            }

            /* Produce W p4 */
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> w_onshell;
            ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> w_offshell;

            w_onshell = nu_onshell + l_onshell;
            w_offshell = nu_offshell + l_offshell;

//            std::cout<<"\tOn-shell W  : "<<w_onshell<<std::endl;
//            std::cout<<"\tOff-shell W : "<<w_offshell<<std::endl;


            /* Produce Higgs p4 */
            auto hToww = w_onshell + w_offshell;
            auto hTobb = tmp_b1 + tmp_b2;
            auto h2Tohh = hTobb + hToww;

//            std::cout<<"\th->bb : "<<hTobb<<std::endl;
//            std::cout<<"\th->WW : "<<hToww<<std::endl;
//            std::cout<<"\tH->hh : "<<h2Tohh<<std::endl;

            /* Physics checks */
            if (fabs(hToww.M() - mh) > 1.){
//                std::cout<<"\t-> WW mass "<<hToww.M()<<" too far from mH "<<mh<<std::endl;
                continue; // H->WW on shell condition
            }
            if (w_offshell.M() > hToww.M()/2){
//                std::cout<<"W off-shell mass "<<w_offshell.M()<<" > mHH/2"<<std::endl;
                continue; // W offshell should not be larger than mH/2
            }
//            std::cout<<"\tmW (on-shell)  = "<<w_onshell.M()<<std::endl;
//            std::cout<<"\tmW (off-shell) = "<<w_offshell.M()<<std::endl;
//            std::cout<<"\tmH (h->WW)     = "<<hToww.M()<<std::endl;
//            std::cout<<"\tmH (h->bb)     = "<<hTobb.M()<<std::endl;
//            std::cout<<"\tmH (H->hh)     = "<<h2Tohh.M()<<std::endl;

            /* If valid solution, record value */
            hme += h2Tohh.M();
            solutions++;
        }
        if (solutions > 0){
            hme /= solutions; // average over all possible solutions        
        }
        if (hme != 0.){
            // If we have a value, put it in the binning
            bins[int(std::min(hme,xmax)/xdelta)] ++;   
            converged ++;
        }
    }
    /* Efficiency */
    float eff = float(converged)/iterations *100;
    /* Find max value in the binning */
    std::vector<float>::iterator binmax = std::max_element(bins.begin(),bins.end());       
    float hme = std::distance(bins.begin(), binmax) * xdelta + xmin;

    return std::make_pair(hme,eff);
}
