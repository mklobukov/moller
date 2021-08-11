double A(int p, int m) 
{
    double a = ((double)p - (double)m ) / ((double)p + (double)m);
    cout << "Calculating asymmetry from pos/neg: " << p << ", " << m << " -> " << a << endl;
    return a;
}

double U(int p, int m)
{
    double u = 2.0 * sqrt( ((double)p*(double)m) / pow((double)(p+m), 3));
    cout << "uncertainty = " << u << "\n" << endl;
    return u;
}

void FillAsymmetryHistogram(int nbins, TH1D* a, TH1D* plus, TH1D* minus)
{
    for (int i = 0; i < nbins; i++)
    {
        double numPlus = (double)plus->GetBinContent(i);
        double numMinus = (double)minus->GetBinContent(i);
        
        double asymmetry = A(numPlus, numMinus);
        double asymmetryUncertainty = U(numPlus, numMinus);

        if (numPlus + numMinus > 0)
        {
            a->SetBinContent(i, asymmetry);
            a->SetBinError(i, asymmetryUncertainty);
        }
        else 
        {
            a->SetBinContent(i, 0);
            a->SetBinError(i, 0);   
        }
    }
}


void plot2cmSlabAsymmetriesSeparately_LogBin() 
{
    TCanvas *c1 = new TCanvas();
    TCanvas *c2 = new TCanvas();
    // c1->SetLogx();
    // c2->SetLogx();

    TFile *input1 = new TFile("pol8plusCutoffSqrt1.root", "read"); 
    TFile *input2 = new TFile("pol8minusCutoffSqrt1.root", "read"); 
    TFile *input3 = new TFile("pol50plusCutoffSqrt1.root", "read"); 
    TFile *input4 = new TFile("pol50minusCutoffSqrt1.root", "read"); 
    TFile *input5 = new TFile("pol100plusCutoffSqrt1.root", "read"); 
    TFile *input6 = new TFile("pol100minusCutoffSqrt1.root", "read"); 

    TTree *tree1 = (TTree*)input1->Get("Particles emerging from target");
    TTree *tree2 = (TTree*)input2->Get("Particles emerging from target");
    TTree *tree3 = (TTree*)input3->Get("Particles emerging from target");
    TTree *tree4 = (TTree*)input4->Get("Particles emerging from target");
    TTree *tree5 = (TTree*)input5->Get("Particles emerging from target");
    TTree *tree6 = (TTree*)input6->Get("Particles emerging from target");

        // Create bins:
    int nbins = 10;
    int xmax = 11000;
    int xmin = 1;
    Double_t *xbins    = new Double_t[nbins+1];
    Double_t xlogmin = TMath::Log10(xmin);
    Double_t xlogmax = TMath::Log10(xmax);
    Double_t dlogx   = (xlogmax-xlogmin)/((Double_t)nbins);
    for (int i = 0; i <= nbins; i++)
    { 
        Double_t xlog = xlogmin+ i*dlogx;
        xbins[i] = TMath::Exp( TMath::Log(10) * xlog ); 
        cout << "current bin: " << xbins[i] << endl;
    }

    #pragma region hist declaration
    auto h8e_plus = new TH1D("h8e_plus","",nbins,xbins);
    auto h8e_minus = new TH1D("h8e_minus","",nbins,xbins);
    auto h8p_plus = new TH1D("h8p_plus","",nbins,xbins);
    auto h8p_minus = new TH1D("h8p_minus","",nbins,xbins);

    auto h50e_plus = new TH1D("h50e_plus","",nbins,xbins);
    auto h50e_minus = new TH1D("h50e_minus","",nbins,xbins);
    auto h50p_plus = new TH1D("h50p_plus","",nbins,xbins);
    auto h50p_minus = new TH1D("h50p_minus","",nbins,xbins);
    
    auto h100e_plus = new TH1D("h100e_plus","",nbins,xbins);
    auto h100e_minus = new TH1D("h100e_minus","",nbins,xbins);
    auto h100p_plus = new TH1D("h100p_plus","",nbins,xbins);
    auto h100p_minus = new TH1D("h100p_minus","",nbins,xbins);

    auto asymmetries8_e = new TH1D("asymmetries8_e","asymmetries8_e",nbins,xbins);
    auto asymmetries8_p = new TH1D("asymmetries8_p","asymmetries8_p",nbins,xbins);
    auto asymmetries50_e = new TH1D("asymmetries50_e","asymmetries50_e",nbins,xbins);
    auto asymmetries50_p = new TH1D("asymmetries50_p","asymmetries50_p",nbins,xbins);
    auto asymmetries100_e = new TH1D("asymmetries100_e","asymmetries100_e",nbins,xbins);
    auto asymmetries100_p = new TH1D("asymmetries100_p","asymmetries100_p",nbins,xbins);
    #pragma endregion

    h8e_plus->Print("all");


    char *pName = new char[10];
    tree1->SetBranchAddress("ParticleName", pName);

    double e;
    tree1->SetBranchAddress("Energy", &e);

    #pragma region Fill Histos 
    // fill histos for pol8 plus
    for (int i = 0; i < tree1->GetEntries(); i++) 
    {
        tree1->GetEntry(i);
        if (!strcmp(pName, "e-")) h8e_plus->Fill(e);
        else if (!strcmp(pName, "gamma")) h8p_plus->Fill(e);
    }
    h8e_plus->Print("all");
    h8p_plus->Print("all");

    // fill histos for pol8 minus
    tree2->SetBranchAddress("ParticleName", pName);
    tree2->SetBranchAddress("Energy", &e);
    for (int i = 0; i < tree2->GetEntries(); i++) 
    {
        tree2->GetEntry(i);
        if (!strcmp(pName, "e-")) h8e_minus->Fill(e);
        else if (!strcmp(pName, "gamma")) h8p_minus->Fill(e);
    }

    // fill histos for pol50 plus
    tree3->SetBranchAddress("ParticleName", pName);
    tree3->SetBranchAddress("Energy", &e);
    for (int i = 0; i < tree3->GetEntries(); i++) 
    {
        tree3->GetEntry(i);
        if (!strcmp(pName, "e-")) h50e_plus->Fill(e);
        else if (!strcmp(pName, "gamma")) h50p_plus->Fill(e);
    }

    // fill histos for pol50 minus
    tree4->SetBranchAddress("ParticleName", pName);
    tree4->SetBranchAddress("Energy", &e);
    for (int i = 0; i < tree4->GetEntries(); i++) 
    {
        tree4->GetEntry(i);
        if (!strcmp(pName, "e-")) h50e_minus->Fill(e);
        else if (!strcmp(pName, "gamma")) h50p_minus->Fill(e);
    }

    // fill histos for pol100 plus
    tree5->SetBranchAddress("ParticleName", pName);
    tree5->SetBranchAddress("Energy", &e);
    for (int i = 0; i < tree5->GetEntries(); i++) 
    {
        tree5->GetEntry(i);
        if (!strcmp(pName, "e-")) h100e_plus->Fill(e);
        else if (!strcmp(pName, "gamma")) h100p_plus->Fill(e);
    }

    // fill histos for pol50 minus
    tree6->SetBranchAddress("ParticleName", pName);
    tree6->SetBranchAddress("Energy", &e);
    for (int i = 0; i < tree6->GetEntries(); i++) 
    {
        tree6->GetEntry(i);
        if (!strcmp(pName, "e-")) h100e_minus->Fill(e);
        else if (!strcmp(pName, "gamma")) h100p_minus->Fill(e);
    }

    #pragma endregion 

    #pragma region Fill asymmetries bin by bin
    FillAsymmetryHistogram(nbins, asymmetries8_e, h8e_plus, h8e_minus);
    FillAsymmetryHistogram(nbins, asymmetries8_p, h8p_plus, h8p_minus);
    FillAsymmetryHistogram(nbins, asymmetries50_e, h50e_plus, h50e_minus);
    FillAsymmetryHistogram(nbins, asymmetries50_p, h50p_plus, h50p_minus);
    FillAsymmetryHistogram(nbins, asymmetries100_e, h100e_plus, h100e_minus);
    FillAsymmetryHistogram(nbins, asymmetries100_p, h100p_plus, h100p_minus);
    #pragma endregion

    gStyle->SetEndErrorSize(3);	
    gStyle->SetEndErrorSize(2);
    gStyle->SetErrorX(0);



    c1->cd();
    asymmetries8_e->SetMarkerStyle(20);
    asymmetries50_e->SetMarkerStyle(20);
    asymmetries100_e->SetMarkerStyle(20);

    asymmetries8_e->SetMarkerSize(1);
    asymmetries50_e->SetMarkerSize(1);
    asymmetries100_e->SetMarkerSize(1);

    asymmetries8_e->SetMarkerColor(kBlue);
    asymmetries50_e->SetMarkerColor(kRed);
    asymmetries100_e->SetMarkerColor(kGreen);


    asymmetries8_e->Draw("e1 same");
    asymmetries50_e->Draw("e1 same");
    asymmetries100_e->Draw("e1 same");

    TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
    legend2->AddEntry("asymmetries8_e","8% pol","p");
    legend2->AddEntry("asymmetries50_e","50% pol","p");
    legend2->AddEntry("asymmetries100_e","100% pol","p");
    legend2->SetTextSize(0.025);
    legend2->Draw();


    c2->cd();
    asymmetries8_p->SetMarkerStyle(20);
    asymmetries50_p->SetMarkerStyle(20);
    asymmetries100_p->SetMarkerStyle(20);

    asymmetries8_p->SetMarkerSize(1);
    asymmetries50_p->SetMarkerSize(1);
    asymmetries100_p->SetMarkerSize(1);

    asymmetries8_p->SetMarkerColor(kBlue);
    asymmetries50_p->SetMarkerColor(kRed);
    asymmetries100_p->SetMarkerColor(kGreen);

    asymmetries8_p->Draw("e1 same");
    asymmetries50_p->Draw("e1 same");
    asymmetries100_p->Draw("e1 same");

    TLegend *legend3 = new TLegend(0.1,0.7,0.48,0.9);
    legend3->AddEntry("asymmetries8_p","8% pol","p");
    legend3->AddEntry("asymmetries50_p","50% pol","p");
    legend3->AddEntry("asymmetries100_p","100% pol","p");
    legend3->SetTextSize(0.025);
    legend3->Draw();

    // save graphs
    TFile *f = new TFile("asymmBinsasdfasfffdf29.root", "recreate");
    // add multigraphs
    f->GetList()->Add(asymmetries8_e);
    f->GetList()->Add(asymmetries8_p);
    f->GetList()->Add(asymmetries50_e);
    f->GetList()->Add(asymmetries50_p);
    f->GetList()->Add(asymmetries100_e);
    f->GetList()->Add(asymmetries100_p);

    // Write to file
    f->Write();
}