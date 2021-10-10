
void replaceAll( string &s, const string &search, const string &replace ) {
    for(size_t pos = 0; ; pos += replace.length()) {
        // Locate the substring to replace
        pos = s.find( search, pos );
        if( pos == string::npos ) break;
        // Replace by erasing and inserting
        s.erase( pos, search.length() );
        s.insert( pos, replace );
    }
}


void rename(std::string filepath, std::string input, std::string output, bool debug=false){
    TFile* F = new TFile(filepath.c_str(),"UPDATE");
    unsigned int count = 0;
    for (auto&& keyAsObj : *F->GetListOfKeys()){
        auto key = (TKey*) keyAsObj;
        std::string histname = std::string(key->GetName());
        std::string classname = std::string(key->GetClassName());
        if (classname.rfind("TH", 0) != 0) // is TH so histogram
            continue;
        if (histname.find(input) != std::string::npos){
            count++;
            if (!debug){
                std::string outname = histname;
                replaceAll(outname,input,output);
                auto h = F->Get(histname.c_str());
                h->Write(outname.c_str(),TObject::kOverwrite);
                F->Delete((histname+";1").c_str());
            }
        }
    }
    if (count > 0)
        std::cout<<"Processing "<<filepath<<" : renamed "<<count<<std::endl;
    F->Close();
}
