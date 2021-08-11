#include <math.h> 

void billShapeCheck_Abs()
{

	TCanvas *c1 = new TCanvas("c1", "");
	TCanvas *c2 = new TCanvas("c2", "");
	bool drawAsymmetry = false;
	c1->cd();

	//  TFile *input = new TFile("backup_3B_fullyPolTarget\\dataTESTING.root", "read"); // plus polarization
//	  TFile *input2 = new TFile("backup_3B_fullyPolTarget\\dataTESTINGminus.root", "read"); // minus polarization polarization

	 TFile *input = new TFile("4B_halfPolarized\\data.root", "read"); // plus polarization
	 TFile *input2 = new TFile("4B_halfPolarized\\data-minus.root", "read"); // minus polarization polarization

	int numBins = 50;
	int startBin = 0;
	int endBin = 49;

	if (drawAsymmetry)
	{
		numBins = 20;
		startBin = 7;
		endBin = 20;
	}
	
	try
	{	
		TTree *tree = (TTree*)input->Get("Polarized Ionization >1 GeV");
		TTree *tree2 = (TTree*)input2->Get("Polarized Ionization >1 GeV");
		int entries = tree->GetEntries(); 
		int entries2 = tree2->GetEntries(); 
		
		cout << "Entries: " << endl;
		cout << std::to_string(entries) << endl;
		cout << std::to_string(entries2) << endl;
		
			// primaries 
		double cosa, cosa2;
		tree->SetBranchAddress("CosineCMAngle", &cosa);
		tree2->SetBranchAddress("CosineCMAngle", &cosa2);
		
		// secondaries 
		double sec_cosa, sec_cosa2; // plus and minus
		tree->SetBranchAddress("SecondaryCOMAngleCosine", &sec_cosa);
		tree2->SetBranchAddress("SecondaryCOMAngleCosine", &sec_cosa2);
		
		// TGraph* g = new TGraph(entries);
		// string title = "a";
		// g->SetTitle("Cosine (COM angle) distribution from both helicities;cos(COM Angle);Counts");
		// g->SetMarkerStyle(21);
		// g->SetMarkerSize(0.5);
		
		auto h1 = new TH1D("h1", "Absolute value of cosine (COM angle) for 50% polarized target;| cos(COM Angle) |;Counts", numBins, 0, 1);
		auto h2 = new TH1D("h2", "Just positive;cos(COM Angle);Counts", numBins, -1, 1);
		auto h3 = new TH1D("h3", "Just negative;cos(COM Angle);Counts", numBins, -1, 1);
		auto h4= new TH1D("h4", "Asymmetry vs absolute value of cosine (COM Angle) for 50% polarized target;| cos(COM Angle) |;Asymmetry", numBins, -1, 1);
		
		auto h5= new TH1D("h5", "Asymmetry vs Cosine (COM Angle) for secondaries;cos(COM Angle);Asymmetry", numBins, -1, 1);
		
		
		// positives 
		for (int i = 0; i < entries; i++) 
		{
			tree->GetEntry(i);
			
			// fill the cosine distribution for primaries 
			h1->Fill(abs(cosa));
			
			// fill the cosine distribution for secondaries 
			h5->Fill(abs(sec_cosa));
			
			// fill the asymmetry histogram 
			h4->Fill(abs(cosa));
			
			// fill just the positives histogram 
			h2->Fill(abs(cosa));
		}
		
		// negatives 
		for (int i = 0; i < entries2; i++) 
		{
			tree2->GetEntry(i);
			
			// fill the cosine distribution for primaries 
			h1->Fill(abs(cosa2));
			
			// fill the cosine distribution for secondaries 
			h5->Fill(abs(sec_cosa2));
			
			// fill the asymmetry histogram 
			h4->Fill(abs(cosa));
			
			// fill just the negatives histogram 
			h3->Fill(abs(cosa2));
		}
	
	    h1 ->Print("all");
		h2 ->Print("all");
		h3 ->Print("all");
								
		
		// fill the asymmetry histogram bin by bin
		auto g = new TGraph();
		
		// for (int i = startBin, j = 0; i <=endBin; i++, j++) {
			// double numPlus = (double)h2->GetBinContent(i);
			// double numMinus = (double)h3->GetBinContent(i);
			// double error = 2. * sqrt(numPlus * numMinus / pow((numPlus + numMinus),3));
			// cout << "current bin: " << i << endl;
			// cout << "curr error: " << error << endl;
			// cout << numPlus << " " << numMinus << endl;
			
			// if (error < minErr) { minErr = error; }
			// else if (error >= maxErr) { maxErr = error;}
			
			
			// // double asymmetry = (numPlus - numMinus) / (numPlus + numMinus);
			// // double asymmetryUncertainty = 2. * sqrt(numPlus * numMinus / pow((numPlus + numMinus),3));
			
			// // cout << "Asymmetry and uncertainty: " << asymmetry << " " << asymmetryUncertainty << endl;
			// // cout << "GetBIn: " << h2->GetBin(i) << endl;
			// // g->SetPoint(j, h2->GetBinCenter(i), asymmetry);

			 // //h1->SetBinError(i, error);
		// }
		
			for (int i = 0; i <=h2->GetNbinsX(); i++) {
				double numPlus = (double)h2->GetBinContent(i);
				double numMinus = (double)h3->GetBinContent(i);
				double error = 2. * sqrt(numPlus * numMinus / pow((numPlus + numMinus),3));

				cout << numPlus << " " << numMinus << endl;
				
				double asymmetry = (numPlus - numMinus) / (numPlus + numMinus);
				double asymmetryUncertainty = 2. * sqrt(numPlus * numMinus / pow((numPlus + numMinus),3));
				
				if (numPlus + numMinus > 0)
				{
					h4->SetBinContent(i, asymmetry);
					h4->SetBinError(i, asymmetryUncertainty);
				}
				else {
					h2->SetBinContent(i, 0);
					h4->SetBinError(i, 0);
				}
				cout << "Asymmetry and uncertainty: " << asymmetry << " " << asymmetryUncertainty << endl;
		}
		
		gStyle->SetEndErrorSize(3);	
	    gStyle->SetEndErrorSize(2);
		gStyle->SetErrorX(0);

		
		cout << "integral from bin 0 to bin " << h1->GetNbinsX() << " = " << h1->Integral(0,h1->GetNbinsX(), "width")<< endl;
		
		TF1 *f3 = new TF1("f3", "(3+x^2)^2 / (1 - x^2)^2", -1, 1);
		double integralOfCurve = f3->Integral(-.7, .7);


		double integral = h1->Integral(0,h1->GetNbinsX(), "width");
		cout << "dir int: " << integral << endl;
		// double factorForFunc = integral / 25.3545;
		double factorForFunc = integral / integralOfCurve;
	
		cout << "Factor: " << factorForFunc << endl; // 1.65441
		
		TF1 *f1 = new TF1("f1", "[coeff]* (3+x^2)^2 / (1 - x^2)^2", 0, 1); // (3 + Cos\[Theta]^2)^2/(1 - Cos\[Theta]^2)^2 from 0 to pi
		f1->SetParameter("coeff", factorForFunc);
		f1->SetLineColor(kRed + 2);
		f1->SetLineStyle(2);

		
		// Draw cosine COM distribution for primary 
		h1->SetMarkerColor(kBlue);
		h1->SetMinimum(0);
		h1->SetMarkerStyle(20);
		h1->SetMarkerSize(0.8);
		h1->Draw("P");
		h1->Draw("E1 SAME");
		
		f1->SetLineStyle(2);
	   	TGraph* g2 = new TGraph(f1);
		g2->Draw("same");

		// Legend for the distribution of cosines 
		 TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
	   legend->AddEntry("h1","Simulation data","p");
		legend->AddEntry("f1","William Henry","l");
	   legend->SetTextSize(0.025);
	   legend->Draw();
		
		// // Draw cosine COM distribution for secondary 
		// h5->SetMarkerColor(kGreen);
		// h5->SetMinimum(0);
		// //h1->GetXaxis()->SetRangeUser(0, 1);
		// h5->SetMarkerStyle(20);
		// h5->SetMarkerSize(0.8);
		// h5->Draw("P SAME");
		// h5->Draw("E1 SAME");
		
		
		c2->cd();
		//Draw asymmetry parabola 
		h4->SetFillColor(kBlue);
		h4->SetMarkerStyle(20);
		h4->SetMarkerSize(0.8);
		h4->GetYaxis()->SetRangeUser(-1, 1);
		h4->GetXaxis()->SetRangeUser(0,1);
		h4->Draw("P");
		h4->Print("all");
		h4->Draw("E1 SAME");
		
		// expected asymmetry curve
		TF1 *f2 = new TF1("f2", "-0.5 * ((7 + x^2) * (1-x^2) ) / ( 3 + x^2 )^2",  0, 1);
		f2->SetLineColor(kRed + 2);
		f2->SetLineStyle(2);
		TGraph* g3 = new TGraph(f2);
		g3->Draw("same");
		
		// Legend for asymmetry
	   TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
	   legend2->AddEntry("h4","Calculated asymmetry","p");
		legend2->AddEntry("f2","William Henry","l");
	   legend2->SetTextSize(0.025);
	   legend2->Draw();
		//g->SetMarkerStyle(22);
		//g->Draw("AP");
		 // TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   // legend->AddEntry("g","Simulation","ap");
   // legend->AddEntry("g2","Theory","l");
   // legend->SetTextSize(0.025);
    // legend->Draw();
		//g->SetMaximum(0.1);
		// TF1 *f1 = new TF1("f1", "11000/  2 * (1 + cos(x))", 0, 3.14);
		// f1->SetLineColor(kRed + 2);
		// f1->SetLineStyle(1);
		
		
		// TGraph* g2 = new TGraph(f1);
		// g2->Draw("same");
	
	}
	catch (...)
	{
		input->Close();
	}

	input->Close();
}