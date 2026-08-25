//Make a toy discrete profiling envelope plot

void envelope_plot() {

  //number of points to test and scanning range
  const int npoints = 20;
  const double xmin(-3.), xmax(3.);

  //store the results for each point
  double xs[npoints], y1s[npoints], y2s[npoints], y3s[npoints], yenv[npoints];

  //scan each point
  for(int ipoint = 0; ipoint < npoints; ++ipoint) {
    //x for current point
    const double x = xmin + (ipoint)*(xmax-xmin)/(npoints-1);
    xs[ipoint] = x;
    //make different parabolas for different toy functions
    y1s[ipoint] = 2.0*std::pow((x - 0.0),2) + 0.0;
    y2s[ipoint] = 1.8*std::pow((x - 0.5),2) + 0.5;
    y3s[ipoint] = 2.5*std::pow((x + 0.2),2) + 1.0;
    //store the minimum for the envelope
    yenv[ipoint] = std::min(y1s[ipoint], std::min(y2s[ipoint], y3s[ipoint]));
  }

  //Create a canvas for the plot
  TCanvas* c = new TCanvas();

  //Create a graph for each scan and the envelope
  TGraph* g1 = new TGraph(npoints, xs, y1s);
  TGraph* g2 = new TGraph(npoints, xs, y2s);
  TGraph* g3 = new TGraph(npoints, xs, y3s);
  TGraph* genv = new TGraph(npoints, xs, yenv);

  //configure the graph styles
  g1->SetLineColor(kAzure+6);
  g2->SetLineColor(kRed    );
  g3->SetLineColor(kGreen+1);
  g1->SetMarkerColor(kAzure+6);
  g2->SetMarkerColor(kRed    );
  g3->SetMarkerColor(kGreen+1);
  g1->SetLineWidth(3);
  g2->SetLineWidth(3);
  g3->SetLineWidth(3);
  g1->SetMarkerStyle(20);
  g2->SetMarkerStyle(20);
  g3->SetMarkerStyle(20);

  genv->SetLineColor(kBlack);
  genv->SetLineWidth(3);
  genv->SetLineStyle(kDashed);

  //plot the results
  g1->SetTitle(";#mu_{S}; -2*#Deltaln(L)");
  g1->Draw("APL");
  g2->Draw("PL");
  g3->Draw("PL");
  genv->Draw("L");

  g1->GetYaxis()->SetRangeUser(-0.5, 20.);
  g1->GetYaxis()->SetTitleSize(0.05);
  g1->GetYaxis()->SetTitleOffset(0.85);
  g1->GetXaxis()->SetTitleSize(0.05);
  g1->GetXaxis()->SetTitleOffset(0.85);

  //add a legend with toy function names
  TLegend* leg = new TLegend(0.1, 0.7, 0.4, 0.9);
  leg->AddEntry(g1  , "Polynomial" , "PL");
  leg->AddEntry(g2  , "Exponential", "PL");
  leg->AddEntry(g3  , "Powerlaw"   , "PL");
  leg->AddEntry(genv, "Envelope"   , "L" );
  leg->Draw();

  c->SaveAs("figures/discrete_profiling_cartoon.png");
}
