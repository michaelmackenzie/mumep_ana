// Useful physics functions/data
#ifndef __mumep_ana_analysis_tools_physics__
#define __mumep_ana_analysis_tools_physics__

//----------------------------------------------------------------------------------------------------------------------------
// mass info, in MeV
double mass_table(int z, int a) {
  const double u(931.49410242);
  // From Matthew:

  // masses relevant for C-12 stopping target
  if (z ==  6 && a ==  12) return  12.         *u; // C (12)  - mother
  if (z ==  5 && a ==  12) return  12.0143526  *u; // B (12)  - no nucleon knockout
  if (z ==  5 && a ==  11) return  11.009305167*u; // B (11)  - 1 neutron knockout
  if (z ==  5 && a ==  10) return  10.012936862*u; // B (10)  - 2 neutron knockout
  if (z ==  4 && a ==  12) return  12.0269221  *u; // Be(12)  - no nucleon knockout, mu- --> e+
  if (z ==  4 && a ==  11) return  11.02166108 *u; // Be(11)  - 1 proton knockout
  if (z ==  4 && a ==  10) return  10.01353470 *u; // Be(10)  - 1 proton + 1 neutron knockout

  // masses relevant for O-16 stopping target
  if (z ==  8 && a ==  16) return  15.99491461960*u; // O(16)  - mother
  if (z ==  7 && a ==  16) return  16.00610190000*u; // N(16)  - no nucleon knockout
  if (z ==  7 && a ==  15) return  15.00010889890*u; // N(15)  - 1 neutron knockout
  if (z ==  7 && a ==  14) return  14.00307400446*u; // N(14)  - 2 neutron knockout
  if (z ==  6 && a ==  15) return  15.01059930000*u; // C(15)  - 1 proton knockout
  if (z ==  6 && a ==  14) return  14.00324198800*u; // C(14)  - 1 proton + 1 neutron knockout

  // masses relevant for Al-27 stopping target
  if (z == 13 && a ==  27) return  26.98153841*u; // Al(27) - mother
  if (z == 12 && a ==  27) return  26.98434063*u; // Mg(27) - no nucleon knockout
  if (z == 12 && a ==  26) return  25.98259297*u; // Mg(26) - 1 neutron knockout
  if (z == 12 && a ==  25) return  24.98583700*u; // Mg(25) - 2 neutron knockout
  if (z == 11 && a ==  26) return  25.99263500*u; // Na(26) - 1 proton knockout
  if (z == 11 && a ==  25) return  24.98995400*u; // Na(25) - 1 proton + 1 neutron knockout

  // masses relevant for Si-28 stopping target
  if (z == 14 && a ==  28) return  27.976926535*u; // Si(28) - mother
  if (z == 13 && a ==  28) return  27.981910090*u; // Al(28) - no nucleon knockout
  if (z == 13 && a ==  27) return  26.981538410*u; // Al(27) - 1 neutron knockout
  if (z == 13 && a ==  26) return  25.986891860*u; // Al(26) - 2 neutron knockout
  if (z == 12 && a ==  27) return  26.984340630*u; // Mg(27) - 1 proton knockout
  if (z == 12 && a ==  26) return  25.982592970*u; // Mg(26) - 1 proton + 1 neutron knockout

  // masses relevant for Ca-40 stopping target
  if (z == 20 && a ==  40) return  39.962590866*u; // Ca(40) - mother
  if (z == 19 && a ==  40) return  39.963998170*u; // K(40)  - no nucleon knockout
  if (z == 19 && a ==  39) return  38.963706487*u; // K(39)  - 1 neutron knockout
  if (z == 19 && a ==  38) return  37.969081120*u; // K(38)  - 2 neutron knockout
  if (z == 18 && a ==  39) return  38.964313000*u; // Ar(39) - 1 proton knockout
  if (z == 18 && a ==  38) return  37.962732100*u; // Ar(38) - 1 proton + 1 neutron knockout

  // masses relevant for Ti-48 stopping target
  if (z == 22 && a ==  48) return  47.94794093*u; // Ti(48) - mother
  if (z == 21 && a ==  48) return  47.95222300*u; // Sc(48) - no nucleon knockout
  if (z == 21 && a ==  47) return  46.95240270*u; // Sc(47) - 1 neutron knockout
  if (z == 21 && a ==  46) return  45.95516750*u; // Sc(46) - 2 neutron knockout
  if (z == 20 && a ==  47) return  46.95454140*u; // Ca(47) - 1 proton knockout
  if (z == 20 && a ==  46) return  45.95368800*u; // Ca(46) - 1 proton + 1 neutron knockout

  // masses relevant for Ni-58 stopping target
  if (z == 28 && a ==  58) return  57.93534180*u; // Ni(58) - mother
  if (z == 27 && a ==  58) return  57.93575140*u; // Co(58) - no nucleon knockout
  if (z == 27 && a ==  57) return  56.93628990*u; // Co(57) - 1 neutron knockout
  if (z == 27 && a ==  56) return  55.93983820*u; // Co(56) - 2 neutron knockout
  if (z == 26 && a ==  57) return  56.93539210*u; // Fe(57) - 1 proton knockout
  if (z == 26 && a ==  56) return  55.93493560*u; // Fe(56) - 1 proton + 1 neutron knockout

  // masses relevant for Ni-60 stopping target
  if (z == 28 && a ==  60) return  59.93078530*u; // Ni(60) - mother
  if (z == 27 && a ==  60) return  59.93381570*u; // Co(60) - no nucleon knockout
  if (z == 27 && a ==  59) return  58.93319370*u; // Co(59) - 1 neutron knockout
  if (z == 27 && a ==  58) return  57.93575140*u; // Co(58) - 2 neutron knockout
  if (z == 26 && a ==  59) return  58.93487360*u; // Fe(59) - 1 proton knockout
  if (z == 26 && a ==  58) return  57.93327370*u; // Fe(58) - 1 proton + 1 neutron knockout

  // masses relevant for Ni-62 stopping target
  if (z == 28 && a ==  62) return  61.92834490*u; // Ni(62) - mother
  if (z == 27 && a ==  62) return  61.93405800*u; // Co(62) - no nucleon knockout
  if (z == 27 && a ==  61) return  60.93247610*u; // Co(61) - 1 neutron knockout
  if (z == 27 && a ==  60) return  59.93381570*u; // Co(60) - 2 neutron knockout
  if (z == 26 && a ==  61) return  60.93674620*u; // Fe(61) - 1 proton knockout
  if (z == 26 && a ==  60) return  59.93407000*u; // Fe(60) - 1 proton + 1 neutron knockout

  // masses relevant for Zr-90 stopping target
  if (z == 40 && a ==  90) return  89.904698760*u; // Zr(90) - mother
  if (z == 39 && a ==  90) return  89.907144800*u; // Y(90)  - no nucleon knockout
  if (z == 39 && a ==  89) return  88.905841200*u; // Y(89)  - 1 neutron knockout
  if (z == 39 && a ==  88) return  87.909501300*u; // Y(88)  - 2 neutron knockout
  if (z == 38 && a ==  89) return  88.907450810*u; // Sr(89) - 1 proton knockout
  if (z == 38 && a ==  88) return  87.905612256*u; // Sr(88) - 1 proton + 1 neutron knockout

  // masses relevant for Mo-98 stopping target
  if (z == 42 && a ==  98) return  97.90540360*u; // Mo(98) - mother
  if (z == 41 && a ==  98) return  97.91033300*u; // Nb(98) - no nucleon knockout
  if (z == 41 && a ==  97) return  96.90809800*u; // Nb(97) - 1 neutron knockout
  if (z == 41 && a ==  96) return  95.90810159*u; // Nb(96) - 2 neutron knockout
  if (z == 40 && a ==  97) return  96.91095740*u; // Zr(97) - 1 proton knockout
  if (z == 40 && a ==  96) return  95.90827762*u; // Zr(96) - 1 proton + 1 neutron knockout

  // masses relevant for Ag-107 stopping target
  if (z == 47 && a == 107) return 106.90509150*u; // Ag(107) - mother
  if (z == 46 && a == 107) return 106.90512810*u; // Pd(107) - no nucleon knockout
  if (z == 46 && a == 106) return 105.90348030*u; // Pd(106) - 1 neutron knockout
  if (z == 46 && a == 105) return 104.90507950*u; // Pd(105) - 2 neutron knockout
  if (z == 45 && a == 106) return 105.90728600*u; // Rh(106) - 1 proton knockout
  if (z == 45 && a == 105) return 104.90568780*u; // Rh(105) - 1 proton + 1 neutron knockout

  // masses relevant for Sn-120 stopping target
  if (z == 50 && a == 120) return 119.90220190*u; // Sn(120) - mother
  if (z == 49 && a == 120) return 119.90797000*u; // In(120) - no nucleon knockout
  if (z == 49 && a == 119) return 118.90585100*u; // In(119) - 1 neutron knockout
  if (z == 49 && a == 118) return 117.90635700*u; // In(118) - 2 neutron knockout
  if (z == 48 && a == 119) return 118.90985000*u; // Cd(119) - 1 proton knockout
  if (z == 48 && a == 118) return 117.90692200*u; // Cd(118) - 1 proton + 1 neutron knockout

  // masses relevant for Pb-208 stopping target
  if (z == 82 && a == 208) return 207.97665190*u; // Pb(208) - mother
  if (z == 81 && a == 208) return 207.98201800*u; // Tl(208) - no nucleon knockout
  if (z == 81 && a == 207) return 206.97741900*u; // Tl(207) - 1 neutron knockout
  if (z == 81 && a == 206) return 205.97611000*u; // Tl(206) - 2 neutron knockout
  if (z == 80 && a == 207) return 206.98230000*u; // Hg(207) - 1 proton knockout
  if (z == 80 && a == 206) return 205.97751400*u; // Hg(206) - 1 proton + 1 neutron knockout

  if(a == 208) {
    if(z == 82) return 207.9766519*u; // Pb(208)
    if(z == 81) return 207.9820180*u; // Tl(208)
  }
  if(a == 207) {
    if(z == 82) return 206.9758967*u; // Pb(207)
    if(z == 81) return 206.977419 *u; // Tl(207)
    if(z == 80) return 206.982300 *u; // Hg(207)
  }
  if(a == 206) {
    if(z == 82) return 205.9744651*u; // Pb(206)
    if(z == 81) return 205.9761100*u; // Tl(206)
    if(z == 80) return 205.977514 *u; // Hg(206)
    if(z == 79) return 205.984740 *u; // Au(206)
  }
  if(a == 197) {
    if(z == 79) return 196.9665701*u; // Au(197) From Wang et al., 2017
    if(z == 78) return 196.9673431*u; // Pt(197)
    if(z == 77) return 196.969657 *u; // Ir(197)
  }
  if(a == 196) {
    if(z == 79) return 195.966571 *u; // Au(196)
    if(z == 78) return 195.9649547*u; // Pt(196)
    if(z == 77) return 195.968400 *u; // Ir(196)
  }
  if(a == 195) {
    if(z == 79) return 194.9650379*u; // Au(195)
    if(z == 78) return 194.9647944*u; // Pt(195)
    if(z == 77) return 194.9659770*u; // Ir(195)
  }
  if(a == 119) {
    if(z == 50) return 118.9033112*u; // Sn(119)
    if(z == 49) return 118.905851 *u; // In(119)
  }
  if(a == 118) {
    if(z == 50) return 117.9016066*u; // Sn(118)
    if(z == 49) return 117.906357 *u; // In(118)
    if(z == 48) return 117.906922 *u; // Cd(118)
  }
  if(a == 117) {
    if(z == 50) return 116.9029540*u; // Sn(117)
    if(z == 49) return 116.904516 *u; // In(117)
    if(z == 48) return 116.9072260*u; // Cd(117)
  }
  if(a == 107) {
    if(z == 47) return 106.9050915*u; // Ag(107)
    if(z == 46) return 106.9051281*u; // Pd(107)
  }
  if(a == 106) {
    if(z == 47) return 105.906664 *u; // Ag(106)
    if(z == 46) return 105.9034803*u; // Pd(106)
    if(z == 45) return 105.907286 *u; // Rh(106)
  }
  if(a == 105) {
    if(z == 47) return 104.906526 *u; // Ag(105)
    if(z == 46) return 104.9050795*u; // Pd(105)
    if(z == 45) return 104.9056878*u; // Rh(105)
  }
  if(a == 96) {
    if(z == 42) return 95.90467477*u; // Mo(96)
    if(z == 41) return 95.90810159*u; // Nb(96)
  }
  if(a == 95) {
    if(z == 42) return 94.90583744*u; // Mo(95)
    if(z == 41) return 94.9068311 *u; // Nb(95)
    if(z == 40) return 94.9080403 *u; // Zr(95)
  }
  if(a == 94) {
    if(z == 42) return 93.90508359*u; // Mo(94)
    if(z == 41) return 93.9072790 *u; // Nb(94)
    if(z == 40) return 93.90631252*u; // Zr(94)
  }
  if(a == 90) {
    if(z == 40) return 89.90469876*u; // Zr(90)
    if(z == 39) return 89.9071448 *u; // Y(90)
  }
  if(a == 89) {
    if(z == 40) return 88.908882  *u; // Zr(89)
    if(z == 39) return 88.9058412 *u; // Y(89)
    if(z == 38) return 88.90745081*u; // Sr(89)
  }
  if(a == 88) {
    if(z == 40) return 87.910221   *u; // Zr(88)
    if(z == 39) return 87.9095013  *u; // Y(88)
    if(z == 38) return 87.905612256*u; // Sr(88)
  }
  if(a == 62) {
    if(z == 28) return 61.9283449*u; // Ni(62)
    if(z == 27) return 61.934058 *u; // Co(62)
  }
  if(a == 61) {
    if(z == 28) return 60.9310549*u; // Ni(61)
    if(z == 27) return 60.9324761*u; // Co(61)
    if(z == 26) return 60.9367462*u; // Fe(61)
  }
  if(a == 60) {
    if(z == 28) return 59.9307853*u; // Ni(60)
    if(z == 27) return 59.9338157*u; // Co(60)
    if(z == 26) return 59.934070 *u; // Fe(60)
  }
  if(a == 59) {
    if(z == 28) return 58.9343456*u; // Ni(59)
    if(z == 27) return 58.9331937*u; // Co(59)
    if(z == 26) return 58.9348736*u; // Fe(59)
  }
  if(a == 58) {
    if(z == 28) return 57.9353418*u; // Ni(58)
    if(z == 27) return 57.9357514*u; // Co(58)
    if(z == 26) return 57.9332737*u; // Fe(58)
  }
  if(a == 57) {
    if(z == 28) return 56.9397915*u; // Ni(57)
    if(z == 27) return 56.9362899*u; // Co(57)
    if(z == 26) return 56.9353921*u; // Fe(57)
  }
  if(a == 56) {
    if(z == 28) return 55.9421279*u; // Ni(56)
    if(z == 27) return 55.9398382*u; // Co(56)
    if(z == 26) return 55.9349356*u; // Fe(56)
  }
  if(a == 48) {
    if(z == 22) return 47.94794093*u; // Ti(48)
    if(z == 21) return 47.952223  *u; // Sc(48)
  }
  if(a == 47) {
    if(z == 22) return 46.95175775*u; // Ti(47)
    if(z == 21) return 46.9524027 *u; // Sc(47)
    if(z == 20) return 46.9545414 *u; // Ca(47)
  }
  if(a == 46) {
    if(z == 22) return 45.95262686*u; // Ti(46)
    if(z == 21) return 45.9551675 *u; // Sc(46)
  }
  if(a == 40) {
    if(z == 20) return 39.962590866*u; // Ca(40)
    if(z == 19) return 39.96399817 *u; // K(40)
  }
  if(a == 39) {
    if(z == 20) return 38.9707108  *u; // Ca(39)
    if(z == 19) return 38.963706487*u; // K(39)
    if(z == 18) return 38.964313   *u; // Ar(39)
  }
  if(a == 38) {
    if(z == 20) return 37.97631923*u; // Ca(38)
    if(z == 19) return 37.96908112*u; // K(38)
  }
  if(a == 28) {
    if(z == 14) return 27.9769265350*u; // Si(28)
    if(z == 13) return 27.98191009  *u; // Al(28)
  }
  if(a == 27) {
    if(z == 14) return 26.98670469*u; // Si(27)
    if(z == 13) return 26.98153841*u; // Al(27)
    if(z == 12) return 26.98434063*u; // Mg(27)
    if(z == 11) return 26.994076  *u; // Na(27)
  }
  if(a == 26) {
    if(z == 14) return 25.99233380*u; // Si(26)
    if(z == 13) return 25.98689186*u; // Al(26)
    if(z == 12) return 25.98259297*u; // Mg(26)
    if(z == 11) return 25.992635  *u; // Na(26)
  }
  if(a == 25) {
    if(z == 14) return 25.004109  *u; // Si(25)
    if(z == 13) return 24.99042831*u; // Al(25)
    if(z == 12) return 24.98583696*u; // Mg(25)
    if(z == 11) return 24.9899540 *u; // Na(25)
    if(z == 10) return 24.997810  *u; // Ne(25)
  }
  if(a == 16) {
    if(z ==  8) return 15.99491461960*u; // O(16)
    if(z ==  7) return 16.0061019    *u; // N(16)
  }
  if(a == 15) {
    if(z ==  8) return 15.0030656   *u; // O(15)
    if(z ==  7) return 15.0001088989*u; // N(15)
    if(z ==  6) return 15.0105993   *u; // C(15)
  }
  if(a == 14) {
    if(z ==  8) return 14.008596706  *u; // O(14)
    if(z ==  7) return 14.00307400446*u; // N(14)
  }
  // cout << __func__ << ": Unknown atom (" << z << ", " << a << ")\n";
  return -1.;
}

//----------------------------------------------------------------------------------------------------------------------------
// Muon 1s binding energy, in MeV
double binding_table(int z, int a) {
  if(a == 197) {
    if(z == 79) return 10.08123; // Au(197)
  }
  if(a == 27) {
    if(z == 13) return 0.464; // Al(27)
  }
  // approximate the binding energy
  // constexpr double m_e(0.511), m_mu(105.66), r_g(13.6e-6);
  // const double e_1s = m_mu / m_e * r_g * z*z;
  constexpr double m_e(0.511), m_mu(105.66), alpha(1./137.036);
  // const double e_1s_base = m_mu * alpha * alpha /2. * z*z;
  const double e_1s_base = m_mu * (1. - sqrt(1. - z*z*alpha*alpha));

  // Evaluate corrections based on the nuclear charge distribution
  // Taken from https://arxiv.org/pdf/2311.16855
  double correction = 0.;
  double val_below(0.), val_above(0.);
  int    z_below(0), z_above(0);
  for(int z_i = 6; z_i <= 92; ++z_i) {
    double corr = 0.;
    //                                             size        dE pert      dE var       dE num
    if     (/*12C  */ z_i ==  6) corr = -m_mu * (3.8967e-6 + -2.3079e-8 + -2.3290e-8 + -2.3727e-8);
    else if(/*16O  */ z_i ==  8) corr = -m_mu * (1.4057e-5 + -9.4133e-8 + -9.4931e-8 + -9.6493e-8);
    else if(/*20Ne */ z_i == 10) corr = -m_mu * (4.0175e-5 + -2.7638e-7 + -2.7854e-7 + -2.8240e-7);
    else if(/*28Si */ z_i == 14) corr = -m_mu * (1.5229e-4 + -1.2805e-6 + -1.2920e-6 + -1.3090e-6);
    else if(/*38Ar */ z_i == 18) corr = -m_mu * (4.4039e-4 + -3.8294e-6 + -3.8660e-6 + -3.9098e-6);
    else if(/*40Ca */ z_i == 20) corr = -m_mu * (6.6509e-4 + -5.9452e-6 + -6.0055e-6 + -6.0708e-6);
    else if(/*66Zn */ z_i == 30) corr = -m_mu * (3.2385e-3 + -2.8141e-5 + -2.8487e-5 + -2.8730e-5);
    else if(/*86Kr */ z_i == 36) corr = -m_mu * (6.3388e-3 + -5.2290e-5 + -5.2994e-5 + -5.3395e-5);
    else if(/*90Zr */ z_i == 40) corr = -m_mu * (9.1096e-3 + -7.3846e-5 + -7.4906e-5 + -7.5446e-5);
    else if(/*120Sn*/ z_i == 50) corr = -m_mu * (1.9954e-2 + -1.3942e-4 + -1.4157e-4 + -1.4241e-4);
    else if(/*136Xe*/ z_i == 54) corr = -m_mu * (2.5930e-2 + -1.6995e-4 + -1.7262e-4 + -1.7357e-4);
    else if(/*142Nd*/ z_i == 60) corr = -m_mu * (3.6374e-2 + -2.2449e-4 + -2.2817e-4 + -2.2935e-4);
    else if(/*176Yb*/ z_i == 70) corr = -m_mu * (6.0941e-2 + -3.0765e-4 + -3.1265e-4 + -3.1397e-4);
    else if(/*185Re*/ z_i == 75) corr = -m_mu * (7.5168e-2 + -3.6513e-4 + -3.7125e-4 + -3.7277e-4);
    else if(/*208Pb*/ z_i == 82) corr = -m_mu * (9.9579e-2 + -4.4041e-4 + -4.4789e-4 + -4.4958e-4);
    else if(/*209Bi*/ z_i == 83) corr = -m_mu * (1.0346e-1 + -4.5137e-4 + -4.5904e-4 + -4.6076e-4);
    else if(/*212Rn*/ z_i == 86) corr = -m_mu * (1.1588e-1 + -4.8286e-4 + -4.9108e-4 + -4.9284e-4);
    else if(/*238U */ z_i == 92) corr = -m_mu * (1.4530e-1 + -5.2566e-4 + -5.3428e-4 + -5.3598e-4);
    else continue;
    if(z < z_i) { // reached the element
      z_above = z_i;
      val_above = corr;
      break;
    } else if(z >= z_i) {
      z_below = z_i;
      val_below = corr;
    }
  }
  if(z_below == z) correction = val_below;
  else if(z_below == 0) correction = val_above; // below the edge of the table
  else if(z_above == 0) correction = val_below; // above the edge of the table
  else {
    correction = val_below + (z - z_below) * (val_above - val_below) / (z_above - z_below); // linear interpolation between points
  }

  const double e_1s = e_1s_base + correction;
  // cout << __func__ << ": Approximating E_1s for Z " << z << " as " << e_1s << endl;
  return e_1s;
}

//----------------------------------------------------------------------------------------------------------------------------
// OMC fractions
double omc_fraction(TString target) {
  if(target.EndsWith("O"   )) return 0.1844; // Taken from 1999 paper
  if(target.EndsWith("Al"  )) return 0.6095;
  if(target.EndsWith("Si"  )) return 0.6587;
  if(target.EndsWith("Ca"  )) return 0.8508;
  if(target.EndsWith("Ti"  )) return 0.8530;
  if(target.EndsWith("Zr"  )) return 0.9529;
  if(target.EndsWith("Ag"  )) return 0.9634;
  if(target.EndsWith("Mo"  )) return 0.9576; // Taken from 1992 paper
  if(target.EndsWith("Sn"  )) return 0.9615;
  if(target.EndsWith("Pb"  )) return 0.9710;
  return -1.;
}

//----------------------------------------------------------------------------------------------------------------------------
// calculate the RMC endpoint
double calculate_endpoint(int z = 13, int a = 27, int dn = 0, int dp = 0) {
  double mass_i = mass_table(z, a);
  double mass_o = mass_table(z-1-dp, a-dn-dp);
  double binding = binding_table(z,a);
  if(mass_i < 0.) {
    // cout << __func__ << ": Unknown initial mass, z = " << z << " a = " << a << endl;
    return -1.;
  }
  if(binding < 0.) {
    // cout << __func__ << ": Unknown binding energy, z = " << z << " a = " << a << endl;
    return -1.;
  }
  if(mass_o < 0.) {
    // cout << __func__ << ": Unknown final mass, z = " << z-1-dp << " a = " << a - dn - dp << endl;
    return -1.;
  }

  const double m_mu(105.6583745), m_e(0.5109989461), m_p(938.272), m_n(939.565);

  // subtract the electron masses to get the nuclear masses (ignore binding energies)
  mass_i -= m_e*z;
  mass_o -= m_e*(z-1-dp);

  // subtract muon binding and add muon mass to the input state
  mass_i += m_mu - binding;

  // add the additional neutrons/proton to the outgoing state
  mass_o += dn*m_n + dp*m_p;

  // calculate the recoil
  const double recoil = (mass_i*mass_i + mass_o*mass_o)/(2.*mass_i) - mass_o;

  // calculate the endpoint
  const double endpoint = mass_i - mass_o - recoil; // mass_i includes muon mass and binding energy already
  return endpoint;
}

//----------------------------------------------------------------------------------------------------------------------------
// calculate the mu- --> e+ endpoint
double calculate_ep_endpoint(int z = 13, int a = 27, int dn = 0, int dp = 0) {
  double mass_i = mass_table(z, a);
  double mass_o = mass_table(z-2-dp, a-dn-dp);
  double binding = binding_table(z,a);
  if(mass_i < 0.) {
    return -1.;
  }
  if(binding < 0.) {
    return -1.;
  }
  if(mass_o < 0.) {
    return -1.;
  }

  const double m_mu(105.6583745), m_e(0.5109989461), m_p(938.272), m_n(939.565);

  // subtract the electron masses to get the nuclear masses (ignore binding energies)
  mass_i -= m_e*z;
  mass_o -= m_e*(z-2-dp);

  // subtract muon binding and add muon mass to the input state
  mass_i += m_mu - binding;

  // add the additional neutrons/proton to the outgoing state
  mass_o += dn*m_n + dp*m_p;

  // calculate the recoil
  const double recoil = (mass_i*mass_i + mass_o*mass_o)/(2.*mass_i) - mass_o;

  // calculate the endpoint
  const double endpoint = mass_i - mass_o - recoil; // mass_i includes muon mass and binding energy already
  return endpoint;
}

//----------------------------------------------------------------------------------------------------------------------------
double closure(const double k, const double kmax) {
  if(k <= 0. || k >= kmax || kmax <= 0.) return 0.;
  const double x = k/kmax;
  const double w = 20./kmax*(1. - 2.*x + 2.*x*x)*x*(1.-x)*(1.-x);
  return w;
}

//----------------------------------------------------------------------------------------------------------------------------
double closure_integral(double xmin, double xmax) {
  double x = xmin;
  double y = xmax;
  double integral = 1./3.*x*x*(-20.*pow(x,4)+72.*pow(x,3) -105.*x*x + 80.*x - 30.);
  integral       -= 1./3.*y*y*(-20.*pow(y,4)+72.*pow(y,3) -105.*y*y + 80.*y - 30.);
  return integral;
}

//----------------------------------------------------------------------------------------------------------------------------
double plestid_closure_integral(double k_1, double k_2, double kmax, int knockout) {
  if(kmax <= 0.) return 0.;
  if(knockout < 0) return 0.;
  k_1 = std::max(0., std::min(kmax, k_1));
  k_2 = std::max(0., std::min(kmax, k_2));
  if(k_1 >= k_2) return 0.;
  const double power = 2. + 1.5*knockout;
  const double x_1 = k_1 / kmax;
  const double x_2 = k_2 / kmax;
  const double val_1 = (x_1 - 1.)*std::pow(1-x_1, power)*(power*x_1+x_1+1.);
  const double val_2 = (x_2 - 1.)*std::pow(1-x_2, power)*(power*x_2+x_2+1.);
  const double integral = val_2 - val_1;
  return integral;
}

//----------------------------------------------------------------------------------------------------------------------------
double plestid_closure(const double k, const double kmax, const int knockout) {
  if(k <= 0. || k >= kmax || kmax <= 0.) return 0.;
  const double power = 2. + 1.5*knockout;
  const double norm = (power + 1.) * (power + 2.) / kmax;
  const double x = k / kmax;
  const double weight = norm * x * std::pow(1. - x, power);
  return weight;
}


//----------------------------------------------------------------------------------------------------------------------------
void z_from_name(TString dataset, int& z, int& a) {
  z = -1; a = -1;
  if(dataset.EndsWith("_C"   )) {z =  6; a =  12; return; }
  if(dataset.EndsWith("_O"   )) {z =  8; a =  16; return; }
  if(dataset.EndsWith("_Al"  )) {z = 13; a =  27; return; }
  if(dataset.EndsWith("_Si"  )) {z = 14; a =  28; return; }
  if(dataset.EndsWith("_Ca"  )) {z = 20; a =  40; return; }
  if(dataset.EndsWith("_Ti"  )) {z = 22; a =  48; return; }
  if(dataset.EndsWith("_Ni58")) {z = 28; a =  58; return; }
  if(dataset.EndsWith("_Ni60")) {z = 28; a =  60; return; }
  if(dataset.EndsWith("_Ni62")) {z = 28; a =  62; return; }
  if(dataset.EndsWith("_Zr"  )) {z = 40; a =  90; return; }
  if(dataset.EndsWith("_Mo"  )) {z = 42; a =  98; return; }
  if(dataset.EndsWith("_Ag"  )) {z = 47; a = 107; return; }
  if(dataset.EndsWith("_Sn"  )) {z = 50; a = 120; return; }
  if(dataset.EndsWith("_Au"  )) {z = 79; a = 197; return; }
  if(dataset.EndsWith("_Pb"  )) {z = 82; a = 208; return; }
}

#endif
