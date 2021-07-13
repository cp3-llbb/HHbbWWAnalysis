template <typename T>
std::vector<std::vector<float>> getContentFromTH1(const T& h)
{
    std::vector<float> values;
    std::vector<float> errors;
    std::vector<float> edges;
    int Nx = h.GetNbinsX();
    values.reserve(Nx);
    errors.reserve(Nx);
    edges.reserve(Nx+1);
    for (int i = 1 ; i <= Nx ; i++)
    {
        values.emplace_back(h.GetBinContent(i));
        errors.emplace_back(h.GetBinError(i));
        edges.emplace_back(h.GetBinLowEdge(i));
    }
    edges.emplace_back(h.GetBinLowEdge(Nx+1));
    
    return std::vector<std::vector<float>> {edges,values,errors};
}

//std::vector<std::vector<std::vector<float>>> getContentFromTH2D(TH2D& h)
//{
//    std::vector<std::vector<float>> values;
//    std::vector<std::vector<float>> errors;
//    std::vector<std::vector<float>> edges {std::vector<float>(),std::vector<float>()};
//    int Nx = h.GetNbinsX();
//    int Ny = h.GetNbinsY();
//    values.reserve(Nx);
//    errors.reserve(Nx);
//    edges[0].reserve(Nx+1);
//    edges[1].reserve(Ny+1);
//    for (int x = 1 ; x <= Nx; x++)
//    {
//        edges[0].emplace_back(h.GetXaxis()->GetBinLowEdge(x));
//        values.emplace_back(std::vector<float>());
//        errors.emplace_back(std::vector<float>());
//        values.back().reserve(Ny);
//        errors.back().reserve(Ny);
//        for (int y = 1 ; y <= Ny; y++)
//        {
//            values.back().emplace_back(h.GetBinContent(x,y));
//            errors.back().emplace_back(h.GetBinError(x,y));
//            if (x == 1)
//                edges[1].emplace_back(h.GetYaxis()->GetBinLowEdge(y));
//        }
//        if (x == 1)
//            edges[1].emplace_back(h.GetYaxis()->GetBinLowEdge(Ny+1));
//    }
//    edges[0].emplace_back(h.GetXaxis()->GetBinLowEdge(Nx+1));
//    std::cout<<"finished"<<std::endl;
//    return std::vector<std::vector<std::vector<float>>> {edges,values,errors};
//}
template <typename T>
std::vector<std::vector<float>> getContentFromTH2(const T& h)
{
    std::vector<float> values;
    std::vector<float> errors;
    std::vector<float> edges;
    int Nx = h.GetNbinsX();
    int Ny = h.GetNbinsY();
    values.reserve(Nx*Ny);
    errors.reserve(Nx*Ny);
    edges.reserve(Nx+Ny+2);
    for (int x = 1 ; x <= Nx; x++)
    {
        for (int y = 1 ; y <= Ny; y++)
        {
            values.emplace_back(h.GetBinContent(x,y));
            errors.emplace_back(h.GetBinError(x,y));
            if (x == 1)
                edges.emplace_back(h.GetYaxis()->GetBinLowEdge(y));
        }
        if (x == 1)
            edges.emplace_back(h.GetYaxis()->GetBinLowEdge(Ny+1));
        edges.emplace_back(h.GetXaxis()->GetBinLowEdge(x));
    }
    edges.emplace_back(h.GetXaxis()->GetBinLowEdge(Nx+1));
    return std::vector<std::vector<float>> {edges,values,errors};
}

TH1F fillTH1(const float* edges, const float* values, const float * errors, int N, std::string name)
{
    TH1F h = TH1F(name.c_str(),name.c_str(),N,edges);
    for (int i = 0 ; i < N ; i++)
    {
        h.SetBinContent(i+1,values[i]);
        h.SetBinError(i+1,errors[i]);
    }    
    return h;
}
TH1D fillTH1(const double* edges, const double* values, const double * errors, int N, std::string name)
{
    TH1D h = TH1D(name.c_str(),name.c_str(),N,edges);
    for (int i = 0 ; i < N ; i++)
    {
        h.SetBinContent(i+1,values[i]);
        h.SetBinError(i+1,errors[i]);
    }    
    return h;
}


TH2F fillTH2(const float* xedges, const float* yedges, const float* values, const float* errors, int Nx, int Ny, std::string name)
{
    TH2F h = TH2F(name.c_str(),name.c_str(),Nx,xedges,Ny,yedges);
    for (int x = 0 ; x < Nx ; x++)
    {
        for (int y = 0 ; y < Ny ; y++)
        {
            h.SetBinContent(x+1,y+1,values[y+x*Nx]);
            h.SetBinError(x+1,y+1,errors[y+x*Nx]);
        }
    }    
    return h;
}

TH2D fillTH2(const double* xedges, const double* yedges, const double* values, const double* errors, int Nx, int Ny, std::string name)
{
    TH2D h = TH2D(name.c_str(),name.c_str(),Nx,xedges,Ny,yedges);
    for (int x = 0 ; x < Nx ; x++)
    {
        for (int y = 0 ; y < Ny ; y++)
        {
            h.SetBinContent(x+1,y+1,values[y+x*Nx]);
            h.SetBinError(x+1,y+1,errors[y+x*Nx]);
        }
    }    
    return h;
}


