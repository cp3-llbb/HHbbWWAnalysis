

void testHLT(std::string path){
    TFile* f = new TFile(path.c_str(),"READ");
    TTree* t = (TTree*)f->Get("Events");

    TObjArray* br = (TObjArray*)t->GetListOfBranches()->Clone();
    for(int i = 0; i < br->GetEntries(); ++i) { 
        std::string name(br->At(i)->GetName()); 
        if (name.rfind("HLT_Mu",0) == 0 || name.rfind("HLT_Ele",0) == 0 ||name.rfind("HLT_Iso",0) == 0)
            std::cout<<name<<std::endl;
    }
    //t->Draw("HLT_IsoMu24");
    //t->Draw("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL");
    //t->Draw("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
    //t->Draw("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    //t->Draw("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    //t->Draw("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ");

    t->Draw("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
    //t->Draw("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL");

}
