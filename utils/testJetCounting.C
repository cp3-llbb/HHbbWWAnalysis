
void testJetCounting(std::string filename){
    // Test the jet selection Resolved and Boosted 
    // Commands :
    //      root -l 
    //      .L testJetCounting.C
    //      testJetCounting("path/to/file.root")
    TFile* f = new TFile(filename.c_str(),"OPEN");
    TH1D* BoostedCase = (TH1D*)f->Get("BoostedCase");  
    TH1D* ResolvedCase = (TH1D*)f->Get("ResolvedCase");  
    TH1D* ExclusiveBoostedCase = (TH1D*)f->Get("ExclusiveBoostedCase"); 
    TH1D* ExclusiveResolvedCase = (TH1D*)f->Get("ExclusiveResolvedCase"); 
    TH1D* BoostedAndResolvedCase = (TH1D*)f->Get("BoostedAndResolvedCase");

    auto b = BoostedCase->GetBinContent(2);      
    auto nb = BoostedCase->GetBinContent(1);     
    auto r = ResolvedCase->GetBinContent(2);
    auto nr = ResolvedCase->GetBinContent(1);
    auto xr = ExclusiveResolvedCase->GetBinContent(2); 
    auto nxr = ExclusiveResolvedCase->GetBinContent(1); 
    auto xb = ExclusiveBoostedCase->GetBinContent(2); 
    auto nxb = ExclusiveBoostedCase->GetBinContent(1); 
    auto br = BoostedAndResolvedCase->GetBinContent(2);   
    auto nbr = BoostedAndResolvedCase->GetBinContent(1);

    std::cout<<"Boosted            case : IN = "<<b<<"  OUT = "<<nb<<" / Total = "<<b+nb<<std::endl;
    std::cout<<"Resolved           case : IN = "<<r<<"  OUT = "<<nr<<" / Total = "<<r+nr<<std::endl;
    std::cout<<"Exclusive Resolved case : IN = "<<xr<<"  OUT = "<<nxr<<" / Total = "<<xr+nxr<<std::endl;
    std::cout<<"Exclusive Boosted  case : IN = "<<xb<<"  OUT = "<<nxb<<" / Total = "<<xb+nxb<<std::endl;
    std::cout<<"Boosted + Resolved case : IN = "<<br<<"  OUT = "<<nbr<<" / Total = "<<br+nbr<<std::endl;
    std::cout<<std::endl;
    if (xr == r-br)
        cout<<"Exclusive Resolved behaves correctly"<<std::endl;
    else
        cout<<"[ERROR] Exclusive Resolved behaves incorrectly"<<std::endl;
    if (xb == b-br)
        cout<<"Exclusive Boosted behaves correctly"<<std::endl;
    else
        cout<<"[ERROR] Exclusive Boosted behaves incorrectly"<<std::endl;
        

}
