// #include <string>
// #include <fstream>
// #include <iomanip>

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


void plot2cmSlabAsymmetriesSeparately() 
{
    TCanvas *c1 = new TCanvas();
    TCanvas *c2 = new TCanvas();
    c1->SetLogx();
    c2->SetLogx();

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

    int n = 6;
    double _e1[n];
    double _e2[n];
    double _e3[n];

#pragma region initailize arrays holding asymmetries and uncertainties
    _e1[0] = sqrt(1);
    _e1[1] = sqrt(10);
    _e1[2] = sqrt(100);
    _e1[3] = sqrt(1000);
    _e1[4] = sqrt(10000);
    _e1[5] = sqrt(100000); 

    _e2[0] = _e1[0] * 1.05;
    _e2[1] = _e1[1] * 1.05;
    _e2[2] = _e1[2] * 1.05;
    _e2[3] = _e1[3] * 1.05;
    _e2[4] = _e1[4] * 1.05;
    _e2[5] = _e1[5] * 1.05;

    _e3[0] = _e1[0] * 0.95;
    _e3[1] = _e1[1] * 0.95;
    _e3[2] = _e1[2] * 0.95;
    _e3[3] = _e1[3] * 0.95;
    _e3[4] = _e1[4] * 0.95;
    _e3[5] = _e1[5] * 0.95;

    double errorsX[n];
    errorsX[0] = 0;
    errorsX[1] = 0;
    errorsX[2] = 0;
    errorsX[3] = 0;
    errorsX[4] = 0;
    errorsX[5] = 0;

    // 8%
    double _a8e[n];
    double _u8e[n];

    double _a8p[n];
    double _u8p[n];

    // 50%
    double _a50e[n];
    double _u50e[n];

    double _a50p[n];
    double _u50p[n];

    // 100%
    double _a100e[n];
    double _u100e[n];

    double _a100p[n];
    double _u100p[n];

    #pragma endregion 

#pragma region initialize counts
    // plus and minus count electrons
    int pcSqrt1 = 0;
    int pcSqrt10 = 0;
    int pcSqrt100 = 0;
    int pcSqrt1000 = 0;
    int pcSqrt10000 = 0;
    int pcSqrt100000 = 0;
    int ncSqrt1 = 0;
    int ncSqrt10 = 0;
    int ncSqrt100 = 0;
    int ncSqrt1000 = 0;
    int ncSqrt10000 = 0;
    int ncSqrt100000 = 0;

    // plus and minus count photons
    int gammapcSqrt1 = 0;
    int gammapcSqrt10 = 0;
    int gammapcSqrt100 = 0;
    int gammapcSqrt1000 = 0;
    int gammapcSqrt10000 = 0;
    int gammapcSqrt100000 = 0;
    int gammancSqrt1 = 0;
    int gammancSqrt10 = 0;
    int gammancSqrt100 = 0;
    int gammancSqrt1000 = 0;
    int gammancSqrt10000 = 0;
    int gammancSqrt100000 = 0;
    #pragma endregion 

    char *pName = new char[10];
    tree1->SetBranchAddress("ParticleName", pName);

    double e;
    tree1->SetBranchAddress("Energy", &e);

    // count e- and gamma in pol8plus
    for (int i = 0; i < tree1->GetEntries(); i++)
    {
        tree1->GetEntry(i);
        if (!strcmp(pName, "e-")) {
            if (e > sqrt(100000)) { pcSqrt100000++; }
            if (e > sqrt(10000)) { pcSqrt10000++; }
            if (e > sqrt(1000)) { pcSqrt1000++; }
            if (e > sqrt(100)) { pcSqrt100++; }
            if (e > sqrt(10)) { pcSqrt10++; }
            if (e > sqrt(1)) { pcSqrt1++; }
        }
        else if (!strcmp(pName, "gamma")) {
            if (e > sqrt(100000)) { gammapcSqrt100000++; }
            if (e > sqrt(10000)) { gammapcSqrt10000++; }
            if (e > sqrt(1000)) { gammapcSqrt1000++; }
            if (e > sqrt(100)) { gammapcSqrt100++; }
            if (e > sqrt(10)) { gammapcSqrt10++; }
            if (e > sqrt(1)) { gammapcSqrt1++; }
        }
    }
    input1->Close();

    // count e- and gamma in pol8minus
    tree2->SetBranchAddress("ParticleName", pName);
    tree2->SetBranchAddress("Energy", &e);
    for (int i = 0; i < tree2->GetEntries(); i++)
    {
        tree2->GetEntry(i);

        if (!strcmp(pName, "e-")) {
            if (e > sqrt(100000)) { ncSqrt100000++; }
            if (e > sqrt(10000)) { ncSqrt10000++; }
            if (e > sqrt(1000)) { ncSqrt1000++; }
            if (e > sqrt(100)) { ncSqrt100++; }
            if (e > sqrt(10)) { ncSqrt10++; }
            if (e > sqrt(1)) { ncSqrt1++; }
        }
        else if (!strcmp(pName, "gamma")) {
            if (e > sqrt(100000)) { gammancSqrt100000++; }
            if (e > sqrt(10000)) { gammancSqrt10000++; }
            if (e > sqrt(1000)) { gammancSqrt1000++; }
            if (e > sqrt(100)) { gammancSqrt100++; }
            if (e > sqrt(10)) { gammancSqrt10++; }
            if (e > sqrt(1)) { gammancSqrt1++; }
        }
    }
    input2->Close();

#pragma region fill asymmetry and uncertainty for pol = 8%
    // fill asymmetry and uncertainty for electron with pol8
    _a8e[0] = A(pcSqrt1, ncSqrt1);
    _a8e[1] = A(pcSqrt10, ncSqrt10);
    _a8e[2] = A(pcSqrt100, ncSqrt100);
    _a8e[3] = A(pcSqrt1000, ncSqrt1000);
    _a8e[4] = A(pcSqrt10000, ncSqrt10000);
    _a8e[5] = A(pcSqrt100000, ncSqrt100000);

    _u8e[0] = U(pcSqrt1, ncSqrt1);
    _u8e[1] = U(pcSqrt10, ncSqrt10);
    _u8e[2] = U(pcSqrt100, ncSqrt100);
    _u8e[3] = U(pcSqrt1000, ncSqrt1000);
    _u8e[4] = U(pcSqrt10000, ncSqrt10000);
    _u8e[5] = U(pcSqrt100000, ncSqrt100000);

    // fill asymmetry and uncertainty for gamma with pol8
    _a8p[0] = A(gammapcSqrt1, gammancSqrt1);
    _a8p[1] = A(gammapcSqrt10, gammancSqrt10);
    _a8p[2] = A(gammapcSqrt100, gammancSqrt100);
    _a8p[3] = A(gammapcSqrt1000, gammancSqrt1000);
    _a8p[4] = A(gammapcSqrt10000, gammancSqrt10000);
    _a8p[5] = A(gammapcSqrt100000, gammancSqrt100000);

    _u8p[0] = U(gammapcSqrt1, gammancSqrt1);
    _u8p[1] = U(gammapcSqrt10, gammancSqrt10);
    _u8p[2] = U(gammapcSqrt100, gammancSqrt100);
    _u8p[3] = U(gammapcSqrt1000, gammancSqrt1000);
    _u8p[4] = U(gammapcSqrt10000, gammancSqrt10000);
    _u8p[5] = U(gammapcSqrt100000, gammancSqrt100000);
#pragma endregion 

    // Reset the counts:
#pragma region reset counts
    // plus and minus count electrons
     pcSqrt1 = 0;
     pcSqrt10 = 0;
     pcSqrt100 = 0;
     pcSqrt1000 = 0;
     pcSqrt10000 = 0;
     pcSqrt100000 = 0;
     ncSqrt1 = 0;
     ncSqrt10 = 0;
     ncSqrt100 = 0;
     ncSqrt1000 = 0;
     ncSqrt10000 = 0;
     ncSqrt100000 = 0;

    // plus and minus count photons
     gammapcSqrt1 = 0;
     gammapcSqrt10 = 0;
     gammapcSqrt100 = 0;
     gammapcSqrt1000 = 0;
     gammapcSqrt10000 = 0;
     gammapcSqrt100000 = 0;
     gammancSqrt1 = 0;
     gammancSqrt10 = 0;
     gammancSqrt100 = 0;
     gammancSqrt1000 = 0;
     gammancSqrt10000 = 0;
     gammancSqrt100000 = 0;
#pragma endregion

    // positive 50% pol counts: 
    tree3->SetBranchAddress("ParticleName", pName);
    tree3->SetBranchAddress("Energy", &e);

    // count e- and gamma in pol50plus
    for (int i = 0; i < tree3->GetEntries(); i++)
    {
        tree3->GetEntry(i);
        if (!strcmp(pName, "e-")) {
            if (e > sqrt(100000)) { pcSqrt100000++; }
            if (e > sqrt(10000)) { pcSqrt10000++; }
            if (e > sqrt(1000)) { pcSqrt1000++; }
            if (e > sqrt(100)) { pcSqrt100++; }
            if (e > sqrt(10)) { pcSqrt10++; }
            if (e > sqrt(1)) { pcSqrt1++; }
        }
        else if (!strcmp(pName, "gamma")) {
            if (e > sqrt(100000)) { gammapcSqrt100000++; }
            if (e > sqrt(10000)) { gammapcSqrt10000++; }
            if (e > sqrt(1000)) { gammapcSqrt1000++; }
            if (e > sqrt(100)) { gammapcSqrt100++; }
            if (e > sqrt(10)) { gammapcSqrt10++; }
            if (e > sqrt(1)) { gammapcSqrt1++; }
        }
    }
    input3->Close();

    // count e- and gamma in pol50minus
    tree4->SetBranchAddress("ParticleName", pName);
    tree4->SetBranchAddress("Energy", &e);
    for (int i = 0; i < tree4->GetEntries(); i++)
    {
        tree4->GetEntry(i);

        if (!strcmp(pName, "e-")) {
            if (e > sqrt(100000)) { ncSqrt100000++; }
            if (e > sqrt(10000)) { ncSqrt10000++; }
            if (e > sqrt(1000)) { ncSqrt1000++; }
            if (e > sqrt(100)) { ncSqrt100++; }
            if (e > sqrt(10)) { ncSqrt10++; }
            if (e > sqrt(1)) { ncSqrt1++; }
        }
        else if (!strcmp(pName, "gamma")) {
            if (e > sqrt(100000)) { gammancSqrt100000++; }
            if (e > sqrt(10000)) { gammancSqrt10000++; }
            if (e > sqrt(1000)) { gammancSqrt1000++; }
            if (e > sqrt(100)) { gammancSqrt100++; }
            if (e > sqrt(10)) { gammancSqrt10++; }
            if (e > sqrt(1)) { gammancSqrt1++; }
        }
    }
    input4->Close();

#pragma region fill asymmetry and uncertainty for pol = 50%
    // fill asymmetry and uncertainty for electron with pol50
    _a50e[0] = A(pcSqrt1, ncSqrt1);
    _a50e[1] = A(pcSqrt10, ncSqrt10);
    _a50e[2] = A(pcSqrt100, ncSqrt100);
    _a50e[3] = A(pcSqrt1000, ncSqrt1000);
    _a50e[4] = A(pcSqrt10000, ncSqrt10000);
    _a50e[5] = A(pcSqrt100000, ncSqrt100000);

    _u50e[0] = U(pcSqrt1, ncSqrt1);
    _u50e[1] = U(pcSqrt10, ncSqrt10);
    _u50e[2] = U(pcSqrt100, ncSqrt100);
    _u50e[3] = U(pcSqrt1000, ncSqrt1000);
    _u50e[4] = U(pcSqrt10000, ncSqrt10000);
    _u50e[5] = U(pcSqrt100000, ncSqrt100000);

    // fill asymmetry and uncertainty for gamma with pol50
    _a50p[0] = A(gammapcSqrt1, gammancSqrt1);
    _a50p[1] = A(gammapcSqrt10, gammancSqrt10);
    _a50p[2] = A(gammapcSqrt100, gammancSqrt100);
    _a50p[3] = A(gammapcSqrt1000, gammancSqrt1000);
    _a50p[4] = A(gammapcSqrt10000, gammancSqrt10000);
    _a50p[5] = A(gammapcSqrt100000, gammancSqrt100000);

    _u50p[0] = U(gammapcSqrt1, gammancSqrt1);
    _u50p[1] = U(gammapcSqrt10, gammancSqrt10);
    _u50p[2] = U(gammapcSqrt100, gammancSqrt100);
    _u50p[3] = U(gammapcSqrt1000, gammancSqrt1000);
    _u50p[4] = U(gammapcSqrt10000, gammancSqrt10000);
    _u50p[5] = U(gammapcSqrt100000, gammancSqrt100000);
#pragma endregion 

  // Reset the counts:
#pragma region reset counts
    // plus and minus count electrons
     pcSqrt1 = 0;
     pcSqrt10 = 0;
     pcSqrt100 = 0;
     pcSqrt1000 = 0;
     pcSqrt10000 = 0;
     pcSqrt100000 = 0;
     ncSqrt1 = 0;
     ncSqrt10 = 0;
     ncSqrt100 = 0;
     ncSqrt1000 = 0;
     ncSqrt10000 = 0;
     ncSqrt100000 = 0;

    // plus and minus count photons
     gammapcSqrt1 = 0;
     gammapcSqrt10 = 0;
     gammapcSqrt100 = 0;
     gammapcSqrt1000 = 0;
     gammapcSqrt10000 = 0;
     gammapcSqrt100000 = 0;
     gammancSqrt1 = 0;
     gammancSqrt10 = 0;
     gammancSqrt100 = 0;
     gammancSqrt1000 = 0;
     gammancSqrt10000 = 0;
     gammancSqrt100000 = 0;
#pragma endregion

 // positive 100% pol counts: 
    tree5->SetBranchAddress("ParticleName", pName);
    tree5->SetBranchAddress("Energy", &e);

    // count e- and gamma in pol100plus
    for (int i = 0; i < tree5->GetEntries(); i++)
    {
        tree5->GetEntry(i);
        if (!strcmp(pName, "e-")) {
            if (e > sqrt(100000)) { pcSqrt100000++; }
            if (e > sqrt(10000)) { pcSqrt10000++; }
            if (e > sqrt(1000)) { pcSqrt1000++; }
            if (e > sqrt(100)) { pcSqrt100++; }
            if (e > sqrt(10)) { pcSqrt10++; }
            if (e > sqrt(1)) { pcSqrt1++; }
        }
        else if (!strcmp(pName, "gamma")) {
            if (e > sqrt(100000)) { gammapcSqrt100000++; }
            if (e > sqrt(10000)) { gammapcSqrt10000++; }
            if (e > sqrt(1000)) { gammapcSqrt1000++; }
            if (e > sqrt(100)) { gammapcSqrt100++; }
            if (e > sqrt(10)) { gammapcSqrt10++; }
            if (e > sqrt(1)) { gammapcSqrt1++; }
        }
    }
    input5->Close();

    // count e- and gamma in pol100minus
    tree6->SetBranchAddress("ParticleName", pName);
    tree6->SetBranchAddress("Energy", &e);
    for (int i = 0; i < tree6->GetEntries(); i++)
    {
        tree6->GetEntry(i);

        if (!strcmp(pName, "e-")) {
            if (e > sqrt(100000)) { ncSqrt100000++; }
            if (e > sqrt(10000)) { ncSqrt10000++; }
            if (e > sqrt(1000)) { ncSqrt1000++; }
            if (e > sqrt(100)) { ncSqrt100++; }
            if (e > sqrt(10)) { ncSqrt10++; }
            if (e > sqrt(1)) { ncSqrt1++; }
        }
        else if (!strcmp(pName, "gamma")) {
            if (e > sqrt(100000)) { gammancSqrt100000++; }
            if (e > sqrt(10000)) { gammancSqrt10000++; }
            if (e > sqrt(1000)) { gammancSqrt1000++; }
            if (e > sqrt(100)) { gammancSqrt100++; }
            if (e > sqrt(10)) { gammancSqrt10++; }
            if (e > sqrt(1)) { gammancSqrt1++; }
        }
    }
    input6->Close();

#pragma region fill asymmetry and uncertainty for pol = 100%
    // fill asymmetry and uncertainty for electron with pol100
    _a100e[0] = A(pcSqrt1, ncSqrt1);
    _a100e[1] = A(pcSqrt10, ncSqrt10);
    _a100e[2] = A(pcSqrt100, ncSqrt100);
    _a100e[3] = A(pcSqrt1000, ncSqrt1000);
    _a100e[4] = A(pcSqrt10000, ncSqrt10000);
    _a100e[5] = A(pcSqrt100000, ncSqrt100000);

    _u100e[0] = U(pcSqrt1, ncSqrt1);
    _u100e[1] = U(pcSqrt10, ncSqrt10);
    _u100e[2] = U(pcSqrt100, ncSqrt100);
    _u100e[3] = U(pcSqrt1000, ncSqrt1000);
    _u100e[4] = U(pcSqrt10000, ncSqrt10000);
    _u100e[5] = U(pcSqrt100000, ncSqrt100000);

    // fill asymmetry and uncertainty for gamma with pol100
    _a100p[0] = A(gammapcSqrt1, gammancSqrt1);
    _a100p[1] = A(gammapcSqrt10, gammancSqrt10);
    _a100p[2] = A(gammapcSqrt100, gammancSqrt100);
    _a100p[3] = A(gammapcSqrt1000, gammancSqrt1000);
    _a100p[4] = A(gammapcSqrt10000, gammancSqrt10000);
    _a100p[5] = A(gammapcSqrt100000, gammancSqrt100000);

    _u100p[0] = U(gammapcSqrt1, gammancSqrt1);
    _u100p[1] = U(gammapcSqrt10, gammancSqrt10);
    _u100p[2] = U(gammapcSqrt100, gammancSqrt100);
    _u100p[3] = U(gammapcSqrt1000, gammancSqrt1000);
    _u100p[4] = U(gammapcSqrt10000, gammancSqrt10000);
    _u100p[5] = U(gammapcSqrt100000, gammancSqrt100000);
#pragma endregion 

    // Plot the asymmetries again cutoff energy
    auto g8e = new TGraphErrors(n, _e1, _a8e, errorsX, _u8e);
    auto g8p = new TGraphErrors(n, _e1, _a8p, errorsX, _u8p);

    auto g50e = new TGraphErrors(n, _e2, _a50e, errorsX, _u50e);
    auto g50p = new TGraphErrors(n, _e2, _a50p, errorsX, _u50p);

    auto g100e = new TGraphErrors(n, _e3, _a100e, errorsX, _u100e);
    auto g100p = new TGraphErrors(n, _e3, _a100p, errorsX, _u100p);

    TMultiGraph *mge = new TMultiGraph();
    mge->SetTitle("Electron Asymmetry vs Cutoff Energy;Cutoff Energy (MeV);Asymmetry");

    TMultiGraph *mgp = new TMultiGraph();
    mgp->SetTitle("Photon Asymmetry vs Cutoff Energy;Cutoff Energy (MeV);Asymmetry");
    

    g8e->SetTitle("e- asymm, 8% target pol");
    g8e->SetMarkerColor(kRed);
    g8e->SetMarkerStyle(20);
    g8p->SetTitle("photon asymm, 8% target pol");
    g8p->SetMarkerColor(kRed);
    g8p->SetMarkerStyle(20);


    g50e->SetTitle("e- asymm, 50% target pol");
    g50e->SetMarkerColor(kGreen);
    g50e->SetMarkerStyle(20);
    g50p->SetTitle("photon asymm, 50% target pol");
    g50p->SetMarkerColor(kGreen);
    g50p->SetMarkerStyle(20);

    g100e->SetTitle("e- asymm, 100% target pol");
    g100e->SetMarkerColor(kBlue);
    g100e->SetMarkerStyle(20);
    g100p->SetTitle("photon asymm, 100% target pol");
    g100p->SetMarkerColor(kBlue);
    g100p->SetMarkerStyle(20);

    // Plot a line on zero
    auto g3 = new TGraph();
    g3->SetPoint(0,0,0);
    g3->SetPoint(1,sqrt(100000), 0);
    g3->SetLineColor(kRed);
    g3->SetLineStyle(1);
    g3->SetMarkerSize(0);
    g3->SetTitle("Asymm = 0 for reference");

    TLegend *legend = new TLegend(0.1,0.25,0.25,0.9);
    legend->AddEntry("g8e","8% target pol","p");
    legend->AddEntry("g50e","50% target pol","p");
    legend->AddEntry("g100e","100% target pol","p");
    legend->SetTextSize(0.025);
    
    mge->Add(g8e);
    mge->Add(g50e);
    mge->Add(g100e);
    mge->Add(g3);

    mgp->Add(g8p);
    mgp->Add(g50p);
    mgp->Add(g100p);
    mgp->Add(g3);

    
    c1->cd();
    mge->Draw("ALP");
    c1->BuildLegend();
    c1->Update();

    c2->cd();
    TLegend *legend2 = new TLegend(0.1,0.25,0.25,0.9);
    legend2->AddEntry("g8p","8% target pol","p");
    legend2->AddEntry("g50p","50% target pol","p");
    legend2->AddEntry("g100p","100% target pol","p");
    legend2->SetTextSize(0.025);
    mgp->Draw("ALP");
    c2->BuildLegend();
    c2->Update();

    // save graphs
    TFile *f = new TFile("asymmGraphsWithOffset.root", "recreate");
    // add multigraphs
    f->GetList()->Add(mge);
    f->GetList()->Add(mgp);

    // add electron graphs only
    f->GetList()->Add(g8e);
    f->GetList()->Add(g50e);
    f->GetList()->Add(g100e);

    // add photon graphs only
    f->GetList()->Add(g8p);
    f->GetList()->Add(g50p);
    f->GetList()->Add(g100p);

    // Write to file
    f->Write();
}