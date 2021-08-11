void brehm() 
{
    TCanvas *c1 = new TCanvas();
    TCanvas *c2 = new TCanvas();

    c1->cd();

    TFile *input1 = new TFile("brehm.root", "read"); 
    TTree *tree1 = (TTree*)input1->Get("Brehm");

    tree1->Draw("PhotonE>>hph", "PhotonPol"); // weighted hist
    TH1D* hph = (TH1D*)gDirectory->Get("hph");
    tree1->Draw("PhotonE>>photonE");
    TH1D* photonE = (TH1D*)gDirectory->Get("photonE"); // photon energy
    hph->Divide(photonE);
    hph->ClearUnderflowAndOverflow();

    for (int i = 0; i < hph->GetNbinsX(); i++) hph->SetBinError(i, 0);
    hph->SetMarkerStyle(20);
    hph->SetMarkerSize(0.5);
    hph->SetMarkerColor(kGreen);
    hph->ClearUnderflowAndOverflow();
    auto hph_gr = new TGraph(hph);
    hph_gr->SetTitle("Mean Photon Circular Polarization; Energy (MeV); Mean Ciruclar Polarization");



    TTree *t = new TTree("t", "tree from digitized.csv");
    t->ReadFile("dig.csv", "x:y:z");
    t->SetMarkerSize(1);
    t->Draw("z:x");
    TGraph *gr = (TGraph*)gPad->GetPrimitive("Graph");

    //gr->SetMarkerSize(2);

    c2->cd();
    hph_gr->SetMarkerStyle(20);
    hph_gr->SetMarkerSize(0.8);
    TAxis *axis = hph_gr->GetXaxis();
    axis->SetLimits(0.,11000.);               

    hph_gr->Draw("AC* same");
    gr->Draw("same");

    TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
    legend2->AddEntry(gr, "Olsen and Maximon","l");
    legend2->AddEntry(hph_gr,"Geant4 Simulation","p");
    legend2->SetTextSize(0.025);
    legend2->Draw();
}