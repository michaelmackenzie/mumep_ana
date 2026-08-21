#ifndef __MUMEP_ANA_TOOLS_FAMILIES__
#define __MUMEP_ANA_TOOLS_FAMILIES__
// Return a given order of a function family

//------------------------------------------------------------------------------------------------------------------
//Create an exponential PDF sum
RooAbsPdf* create_exponential(RooRealVar& obs, const int order, TString name) {
  if(order <= 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  vector<RooRealVar*> coeffs;
  vector<RooExponential*> exps;
  RooArgList pdfs;
  RooArgList coefficients;
  for(int i = 1; i <= order; ++i) {
    TString base = Form("%s_exp_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "_c", base + " power", -0.1, -3., 3.));
    exps.push_back(new RooExponential(base + "_pdf", base + " pdf", obs, *vars.back()));
    pdfs.add(*exps.back());
    if(i < order) {
      coeffs.push_back(new RooRealVar(base + "_n", base + " norm", 0.1, 0., 1.));
      coefficients.add(*coeffs.back());
    }
  }
  if(order == 1) {
    pdfs.at(0)->SetTitle(Form("Exponential PDF, order %i", order));
    return ((RooAbsPdf*) pdfs.at(0));
  }
  return new RooAddPdf(Form("%s_exp_pdf_order_%i", name.Data(), order), Form("Exponential PDF, order %i", order), pdfs, coefficients, true);
}

//------------------------------------------------------------------------------------------------------------------
//Create an power law PDF sum
RooAbsPdf* create_powerlaw(RooRealVar& obs, int order, TString name) {
  if(order <= 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  vector<RooRealVar*> coeffs;
  vector<RooPowerLaw*> pwrs;
  RooArgList pdfs;
  RooArgList coefficients;
  for(int i = 1; i <= order; ++i) {
    TString base = Form("%s_pwr_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "c", base + " power", 1., -100., 1.));
    pwrs.push_back(new RooPowerLaw(base + "pdf", base + " pdf", obs, *vars.back()));
    pdfs.add(*pwrs.back());
    if(i < order) {
      coeffs.push_back(new RooRealVar(base + "n", base + " norm", 0.1, 0., 1.));
      coefficients.add(*coeffs.back());
    }
  }
  if(order == 1) {
    pdfs.at(0)->SetTitle(Form("Power law PDF, order %i", order));
    return ((RooAbsPdf*) pdfs.at(0));
  }
  return new RooAddPdf(Form("%s_pwr_pdf_order_%i", name.Data(), order), Form("Power law PDF, order %i", order), pdfs, coefficients, true);
}

//------------------------------------------------------------------------------------------------------------------
//Create an power law PDF sum from RooGenericPdf
RooAbsPdf* create_generic_powerlaw(RooRealVar& obs, int order, TString name) {
  vector<RooRealVar*> vars;
  RooArgList var_list;
  var_list.add(obs);
  TString formula = "";
  for(int i = 1; i <= order; ++i) {
    TString base = Form("%s_pwr_order_%i_%i_", name.Data(), order, i);
    if(i == 1) {
      formula = "@1*TMath::Power(@0,@2)";
      vars.push_back(new RooRealVar(base + "n", base + " norm", 1., 0., 1.e8));
      var_list.add(*vars.back());
      vars.push_back(new RooRealVar(base + "c", base + " power", -1., -100., 0.));
      var_list.add(*vars.back());
    } else {
      formula += Form(" + @%i*TMath::Power(@0,@%i)", 2*i+1, 2*i+2);
      vars.push_back(new RooRealVar(base + "n", base + " norm", 1., 0., 1.e8));
      var_list.add(*vars.back());
      vars.push_back(new RooRealVar(base + "c", base + " power", -1., -100., 0.));
      var_list.add(*vars.back());
    }
  }
  RooGenericPdf* pdf = new RooGenericPdf(Form("%s_pwr_pdf_order_%i", name.Data(), order), formula.Data(), var_list);
  pdf->SetTitle(Form("Power law PDF, order %i", order));
  return pdf;
}

//------------------------------------------------------------------------------------------------------------------
//Create a Laurent series PDF sum
RooAbsPdf* create_laurent(RooRealVar& obs, int order, TString name) {
  if(order <= 0 || order > 6) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  RooArgList var_list;
  var_list.add(obs);
  TString formula = "";
  //take the default Laurent series with x0 = 0 from H->Zgamma analysis
  if(order == 1) formula = "TMath::Power(@0, -3)"; //Add a default 1st order
  if(order == 2) formula = "TMath::Power(@0, -3) + @1*TMath::Power(@0, -4)";
  if(order == 3) formula = "TMath::Power(@0, -3) + @1*TMath::Power(@0, -4) + @2*TMath::Power(@0, -5)";
  if(order == 4) formula = "TMath::Power(@0, -3) + @1*TMath::Power(@0, -4) + @2*TMath::Power(@0, -5) + @3*TMath::Power(@0, -6)";
  if(order == 5) formula = "TMath::Power(@0, -2) + @1*TMath::Power(@0, -3) + @2*TMath::Power(@0, -4) + @3*TMath::Power(@0, -5) + @4*TMath::Power(@0, -6)";
  if(order == 6) formula = "TMath::Power(@0, -2) + @1*TMath::Power(@0, -3) + @2*TMath::Power(@0, -4) + @3*TMath::Power(@0, -5) + @4*TMath::Power(@0, -6) + @5*TMath::Power(@0, -7)";
  for(int i = 1; i < order; ++i) {
    TString base = Form("%s_lrt_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "n", base + " norm", 1., 0., 1.e8));
    var_list.add(*vars.back());
  }
  RooGenericPdf* pdf = new RooGenericPdf(Form("%s_lrt_pdf_order_%i", name.Data(), order), formula.Data(), var_list);
  pdf->SetTitle(Form("Laurent series PDF, order %i", order));
  return pdf;
}

//Create a 1/(Polynomial) PDF
RooGenericPdf* create_inv_polynomial(RooRealVar& obs, int order, TString name) {
  if(order <= 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  RooArgList var_list;
  var_list.add(obs);
  TString formula = "";
  //define the formula for fixed orders
  if     (order == 1) formula = "1/(@1*@0 + @2)";
  else if(order == 2) formula = "1/(@1*@0*@0 + @2*@0 + @3)";
  else if(order == 3) formula = "1/(@1*@0*@0*@0 + @2*@0*@0 + @3*@0 + @4)";
  else return nullptr;
  for(int i = 0; i < order+1; ++i) { //N(params) = order + 1
    TString base = Form("%s_inv_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "n", base + " norm", 1., -1.e3, 1.e3));
    var_list.add(*vars.back());
  }
  RooGenericPdf* pdf = new RooGenericPdf(Form("%s_inv_pdf_order_%i", name.Data(), order), formula.Data(), var_list);
  pdf->SetTitle(Form("Inverse polynomial series PDF, order %i", order));
  return pdf;
}

//Create a Gaussian + polynomial(order = order) PDF
RooGenericPdf* create_gaus_poly_pdf(RooRealVar& obs, int order, TString name) {
  if(order < 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  RooArgList var_list;
  var_list.add(obs);
  TString formula = "";
  //define the formula for fixed orders
  if     (order == -1) formula = "TMath::Gaus(@0, @1, @2)";
  else if(order ==  0) formula = "TMath::Gaus(@0, @1, @2) + @3";
  else if(order ==  1) formula = "TMath::Gaus(@0, @1, @2) + @3 + @4*@0/90";
  else if(order ==  2) formula = "TMath::Gaus(@0, @1, @2) + @3 + @4*@0/90 + @5*@0*@0/90/90";
  else if(order ==  3) formula = "TMath::Gaus(@0, @1, @2) + @3 + @4*@0/90 + @5*@0*@0/90/90 + @6*@0*@0*@0/90/90/90";
  else return nullptr;
  //add the Gaussian parameters
  vars.push_back(new RooRealVar(Form("%s_gaus_poly_order_%i_g_0", name.Data(), order),
                                Form("%s_gaus_poly_order_%i_g_0", name.Data(), order),
                                60., 50., 70.)); //mean
  vars.push_back(new RooRealVar(Form("%s_gaus_poly_order_%i_g_1", name.Data(), order),
                                Form("%s_gaus_poly_order_%i_g_1", name.Data(), order),
                                11., 5., 20.)); //sigma

  //add the polynomial parameters
  for(int i = 0; i < order+1; ++i) { //N(params) = order + 1 = a +bx + ...
    TString base = Form("%s_gaus_poly_order_%i_p_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "p",
                                  base + "p",
                                  (i == 0) ? 0.5 :  0. ,
                                  (i == 0) ? -3. : -1,
                                  (i == 0) ?  3. :  1));
  }
  for(auto var : vars) var_list.add(*var);
  RooGenericPdf* pdf = new RooGenericPdf(Form("%s_gaus_poly_pdf_order_%i", name.Data(), order), formula.Data(), var_list);
  pdf->SetTitle(Form("Gaussian polynomial series PDF, order %i", order));
  return pdf;
}

//Create a Gaussian + exponential(order = order) PDF
RooGenericPdf* create_gaus_expo_pdf(RooRealVar& obs, int order, TString name) {
  if(order < 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  RooArgList var_list;
  var_list.add(obs);
  TString formula = "";
  //define the formula for fixed orders
  if     (order ==  0) formula = "TMath::Gaus(@0, @1, @2)";
  else if(order ==  1) formula = "TMath::Gaus(@0, @1, @2) + TMath::Exp(@3 + @4*@0/90)";
  else if(order ==  2) formula = "TMath::Gaus(@0, @1, @2) + TMath::Exp(@3 + @4*@0/90) + TMath::Exp(@5 + @6*@0/90)";
  else return nullptr;
  //add the Gaussian parameters
  vars.push_back(new RooRealVar(Form("%s_gaus_expo_order_%i_g_0", name.Data(), order),
                                Form("%s_gaus_expo_order_%i_g_0", name.Data(), order),
                                60., 50., 70.)); //mean
  vars.push_back(new RooRealVar(Form("%s_gaus_expo_order_%i_g_1", name.Data(), order),
                                Form("%s_gaus_expo_order_%i_g_1", name.Data(), order),
                                11., 5., 20.)); //sigma

  //add the exponential parameters
  for(int i = 0; i < order; ++i) {
    TString base = Form("%s_gaus_expo_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "n",
                                  base + "n",
                                  1., -10., 10.));
    vars.push_back(new RooRealVar(base + "p",
                                  base + "p",
                                  -3., -10., 5.));
  }
  for(auto var : vars) var_list.add(*var);
  RooGenericPdf* pdf = new RooGenericPdf(Form("%s_gaus_expo_pdf_order_%i", name.Data(), order), formula.Data(), var_list);
  pdf->SetTitle(Form("Gaussian exponential series PDF, order %i", order));
  return pdf;
}

//Create a Gaussian + power(order = order) PDF
RooGenericPdf* create_gaus_power_pdf(RooRealVar& obs, int order, TString name) {
  if(order < 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  RooArgList var_list;
  var_list.add(obs);
  TString formula = "";
  //define the formula for fixed orders
  if     (order ==  0) formula = "TMath::Gaus(@0, @1, @2)";
  else if(order ==  1) formula = "TMath::Gaus(@0, @1, @2) + @3*TMath::Power(@0/90, @4)";
  else if(order ==  2) formula = "TMath::Gaus(@0, @1, @2) + @3*TMath::Power(@0/90, @4) + @5*TMath::Power(@0/90, @6)";
  else return nullptr;
  //add the Gaussian parameters
  vars.push_back(new RooRealVar(Form("%s_gaus_power_order_%i_g_0", name.Data(), order),
                                Form("%s_gaus_power_order_%i_g_0", name.Data(), order),
                                60., 50., 70.)); //mean
  vars.push_back(new RooRealVar(Form("%s_gaus_power_order_%i_g_1", name.Data(), order),
                                Form("%s_gaus_power_order_%i_g_1", name.Data(), order),
                                11., 5., 20.)); //sigma

  //add the power law parameters
  for(int i = 0; i < order; ++i) {
    TString base = Form("%s_gaus_power_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "n",
                                  base + "n",
                                  1., 0., 3.));
    vars.push_back(new RooRealVar(base + "p",
                                  base + "p",
                                  -1., -10., 3.));
  }
  for(auto var : vars) var_list.add(*var);
  RooGenericPdf* pdf = new RooGenericPdf(Form("%s_gaus_power_pdf_order_%i", name.Data(), order), formula.Data(), var_list);
  pdf->SetTitle(Form("Gaussian power series PDF, order %i", order));
  return pdf;
}

//Create a Chebychev polynomial PDF
RooChebychev* create_chebychev(RooRealVar& obs, int order, TString name) {
  if(order <= 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  RooArgList list;
  std::map<int, std::vector<double>> initial_params;
  initial_params[0] = {1.};
  initial_params[1] = {1., 0.1};
  initial_params[2] = {-1., 0.1, 0.1};
  initial_params[3] = {-0.7, 0.3, -0.07, -0.001};
  initial_params[4] = {-0.7, 0.3, -0.07, -0.001, 0.015};
  initial_params[5] = {-0.7, 0.3, -0.07, -0.001, 0.015, -0.02};
  initial_params[6] = {-0.7, 0.3, -0.07, -0.001, 0.015, -0.02, -0.01};
  for(int i = 1; i <= order; ++i) {
    const bool has_params = initial_params.count(order) && ((int) initial_params[order].size()) > i;
    TString base = Form("%s_chb_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "c", base + " c",
                                  (has_params) ? initial_params[order][i] : 1./pow(10.,i), (i == 0) ? -25. : -1., (i == 0) ? 25. : 1.));
    list.add(*vars.back());
  }
  return new RooChebychev(Form("%s_chb_pdf_order_%i", name.Data(), order), Form("Chebychev PDF, order %i", order), obs, list);
}

//------------------------------------------------------------------------------------------------------------------
RooAbsPdf* create_generic_bernstein(RooRealVar& obs, int order, TString name) {
  if(order <= 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> variables;
  RooArgSet var_set;
  var_set.add(obs);
  TString formula;
  double xmin = obs.getMin();
  double xmax = obs.getMax();
  TString var_form = Form("(%s - %.3f)/%.3f", obs.GetName(), xmin, (xmax - xmin));
  double vals[] = {1.404, 2.443e-1, 5.549e-1, 3.675e-1, 0., 0., 0., 0., 0., 0., 0.};

  for(int ivar = 1; ivar <= order; ++ivar) {
    TString base = Form("%s_bst_order_%i_%i_", name.Data(), order, ivar);
    RooRealVar* v = new RooRealVar(base + "c", base + " c", 1./pow(10,ivar), -5., 5.);
    formula += Form("%.0f*(%s)^%i*(1-%s)^%i*%s", TMath::Binomial(order, ivar), var_form.Data(), ivar, var_form.Data(), order-ivar, v->GetName());
    if(ivar < order) formula += " + ";
    var_set.add(*v);
    variables.push_back(v);
  }
  cout << "######################\n"
       << "#### Bernstein order " << order << " formula: " << formula.Data() << endl
       << "######################\n";
  RooAbsPdf* pdf = new RooGenericPdf(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), formula.Data(), var_set);
  return pdf;
}

//------------------------------------------------------------------------------------------------------------------
//Create a Combine fast Bernstein polynomial PDF
RooAbsPdf* create_fast_bernstein(RooRealVar& obs, const int order, TString name) {
  if(order <= 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  RooArgList list;
  for(int i = 1; i <= order; ++i) {
    TString base = Form("%s_bst_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "c", base + " c", 1./pow(10.,i), -25., 25.));
    list.add(*vars.back());
  }
  RooAbsPdf* pdf;
  //FIXME have a better way to set the template
  switch(order) {
  case 1 : pdf = new RooBernsteinFast<1 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  case 2 : pdf = new RooBernsteinFast<2 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  case 3 : pdf = new RooBernsteinFast<3 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  case 4 : pdf = new RooBernsteinFast<4 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  case 5 : pdf = new RooBernsteinFast<5 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  case 6 : pdf = new RooBernsteinFast<6 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  case 7 : pdf = new RooBernsteinFast<7 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  // case 8 : pdf = new RooBernsteinFast<8 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  // case 9 : pdf = new RooBernsteinFast<9 >(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  // case 10: pdf = new RooBernsteinFast<10>(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list); break;
  default: return nullptr;
  }
  if(!pdf) return nullptr;
  switch(order) {
  case 1 : ((RooBernsteinFast<1 >*) pdf)->protectSubRange(true); break;
  case 2 : ((RooBernsteinFast<2 >*) pdf)->protectSubRange(true); break;
  case 3 : ((RooBernsteinFast<3 >*) pdf)->protectSubRange(true); break;
  case 4 : ((RooBernsteinFast<4 >*) pdf)->protectSubRange(true); break;
  case 5 : ((RooBernsteinFast<5 >*) pdf)->protectSubRange(true); break;
  case 6 : ((RooBernsteinFast<6 >*) pdf)->protectSubRange(true); break;
  case 7 : ((RooBernsteinFast<7 >*) pdf)->protectSubRange(true); break;
  default: break;
  }

  return pdf;
}

//------------------------------------------------------------------------------------------------------------------
//Create a Bernstein polynomial PDF
RooAbsPdf* create_bernstein(RooRealVar& obs, const int order, TString name) {
  if(order <= 0) {
    cout << __func__ << ": Can't create order " << order << " PDF!\n";
    return nullptr;
  }
  vector<RooRealVar*> vars;
  RooArgList list;
  for(int i = 1; i <= order; ++i) {
    TString base = Form("%s_bst_order_%i_%i_", name.Data(), order, i);
    vars.push_back(new RooRealVar(base + "c", base + " c", 1./pow(10.,i), -5., 5.));
    list.add(*vars.back());
  }
  RooAbsPdf* pdf = new RooBernstein(Form("%s_bst_pdf_order_%i", name.Data(), order), Form("Bernstein PDF, order %i", order), obs, list);
  return pdf;
}

#endif
