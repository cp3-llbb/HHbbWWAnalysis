void compareSync(){
    TFile* fT = new TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/sync_bbww_Tallinn_2016_v5.root","READ");
    TTree* tT = (TTree*)fT->Get("syncTree");

    TFile* fL = new TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Synchronization/results/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750.root","READ");
    TTree* tL = (TTree*)fL->Get("Events");
    tL->Draw("PU_weight");
    TObjArray *branchList; 
    branchList  = tL->GetListOfBranches();
    std::string branch;
    auto nBranch     = tL->GetNbranches();

    // Get vector of hist //
    std::vector<std::pair<TH1F*,TH1F*>> hist_list;
    for(int i=0;i<nBranch;i++){
        branch = branchList->At(i)->GetName();

        try {
            tL->Draw((branch+">>hL(100)").c_str(),(branch+"!=-9999").c_str());
            tT->Draw((branch+">>hT(100)").c_str(),(branch+"!=-9999").c_str());
            TH1F* hL = (TH1F*)gROOT->FindObject("hL");
            TH1F* hT = (TH1F*)gROOT->FindObject("hT");
            hL->SetTitle(branch.c_str());
            hL->SetLineWidth(2);
            hT->SetLineWidth(2);
            hL->SetLineColor(kRed-2);
            hT->SetLineColor(kBlue-3);
            hist_list.push_back(std::make_pair(hL,hT));

        } 
        catch (const std::exception& e) { // reference to the base of a polymorphic object
            std::cout << e.what()<<std::endl;; // information from length_error printed
        }
    }

    TCanvas *C = new TCanvas("C","C",800,600); 
    C->Print("syncComparison.pdf[");
    for ( auto hp : hist_list ){
        C->Clear();
        hp.first->Draw();
        hp.second->Draw("same");
        C->Print("syncComparison.pdf");
    }
    C->Print("syncComparison.pdf]");
}
