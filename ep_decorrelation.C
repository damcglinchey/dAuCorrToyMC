////////////////////////////////////////////////////////////////////////////////
//
// Study the effect of an event plane decorrelation on two particle
// correlations
//
// Ported heavily from the python test_epdecorrelation.py and corr_funcs.py
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 2 Oct 2016
//
////////////////////////////////////////////////////////////////////////////////

#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

#include <iostream>
#include <vector>
#include <utility>

using namespace std;

typedef pair<double, double> ValErr;

ValErr vn_3sub(const ValErr &C_AB, const ValErr &C_AC, const ValErr &C_BC);
ValErr calc_cn(TH1D* hcorr, int order);

///-----------------------------------------------------------------------------
///
/// Main function
///
/// param[in] Nevent Number of events to simulate
/// param[in] sigd   Sigma value of the gaussian event plane decorrelation
///
void ep_decorrelation(int Nevent = 200000,
                      double sigd = 0.10)
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  //==========================================================================//
  // SET RUNNING CONDITIONS
  //==========================================================================//

  const int NDET = 5; //< Number of detectors (seperate eta ranges)
  double eta[NDET] = { -3.25, -2.0, 0.0, 2.0, 3.25}; //< eta values for each det

  enum DET {BBCS = 0, FVTXS, CNT, FVTXN, BBCN};
  const char* ldet[NDET] = {"BBCS", "FVTXS", "CNT", "FVTXN", "BBCN"};

  int Npart[NDET] = {74, 17, 5, 9, 10}; //< multiplicity in each detector
  // int Npart[NDET] = {5, 5, 5, 5, 5}; //< multiplicity in each detector

  /// True v2 in each detector
  double v2true[NDET] = {0.82 * 0.07, 0.95 * 0.07, 0.07, 0.52 * 0.07, 0.28 * 0.07};

  const int NV2 = 5;
  int dv2[NV2][3] =
  {
    {CNT, FVTXS, BBCS},
    {CNT, FVTXN, BBCN},
    {CNT, FVTXS, FVTXN},
    {CNT, BBCS, FVTXN},
    {CNT, BBCS, BBCN},
  };

  /// measured v2 for different 3 sub-events
  // v2(CNTBBCSFVTXS) = 0.0721258 +/- 0.000145314
  // v2(CNTBBCNFVTXN) = 0.0544472 +/- 0.000486198
  // v2(CNTFVTXNFVTXS) = 0.0873551 +/- 0.000213309
  // v2(CNTBBCNBBCS) = 0.0936908 +/- 0.000899697
  // v2(CNTFVTXNBBCS) = 0.085891 +/- 0.000238102
  double v2data[5][2] =
  {
    {0.0721, 0.0001},
    {0.0544, 0.0005},
    {0.0874, 0.0002},
    {0.0859, 0.0002},
    {0.0937, 0.0009},
  };


  //==========================================================================//
  // DECLARE VARIABLES
  //==========================================================================//

  //-- Event-wise variables (reused)
  double Psi0, m, Psi, dphi;
  vector<double> phi[NDET];

  TRandom3 *rand = new TRandom3();
  TF1 *fphi = new TF1("fphi",
                      "1+2*[0]*TMath::Cos(2*(x-[1]))",
                      -1 * TMath::Pi(), TMath::Pi());

  //-- Event-summed variables
  TH1D* hdphi[NDET][NDET];
  for (int i = 0; i < NDET; i++)
  {
    for (int j = 0; j < NDET; j++)
    {
      hdphi[i][j] = new TH1D(Form("hdphi_%i_%i", i, j),
                             ";#Delta#phi",
                             96, -1.0 * TMath::Pi(), TMath::Pi());
      hdphi[i][j]->Sumw2();
    }
  }
  ValErr c2[NDET][NDET];

  ValErr v2[NV2];

  TH1D* hphi = new TH1D("hphi",
                        ";#Delta#phi",
                        96, -1.0 * TMath::Pi(), TMath::Pi());

  //-- TGraphs
  TGraphErrors* gv2_3sub;
  TGraphErrors* gv2_data;


  //==========================================================================//
  // GENERATE EVENTS
  //==========================================================================//
  cout << endl;
  cout << "--> Generating Events" << endl;

  for (int ievent = 0; ievent < Nevent; ievent++)
  {
    if (ievent % 10000 == 0)
      cout << "----> Generating event " << ievent << endl;

    // get Psi0 and m
    Psi0 = rand->Uniform(-1.0 * TMath::Pi(), TMath::Pi());
    m = sigd > 0 ? rand->Gaus(0, sigd) : 0;

    fphi->SetParameters(v2true[2], 0);
    for (int j = 0; j < 1000; j++)
      hphi->Fill(fphi->GetRandom());

    //-- Loop over detectors and get phi values of particles
    for (int idet = 0; idet < NDET; idet++)
    {
      Psi = m * eta[idet] + Psi0;

      // get random number of particles
      int N = rand->Poisson((double)Npart[idet]);

      // get phi values
      phi[idet].clear();
      fphi->SetParameters(v2true[idet], Psi);
      for (int j = 0; j < N; j++)
        phi[idet].push_back(fphi->GetRandom());

      // cout << " Nphi[" << idet << "]=" << phi[idet].size() << endl;
      // cout << "     phi: ";
      // for (auto i: phi[idet])
      //   cout << i << " ";
      // cout << endl;

    } // idet

    //-- Loop over pairs of detectors and calculate delta phi
    for (int i = 0; i < NDET; i++)
    {
      for (int j = 0; j < NDET; j++)
      {
        if (i == j) continue;

        for (int ip = 0; ip < phi[i].size(); ip++)
        {
          for (int jp = 0; jp < phi[j].size(); jp++)
          {
            dphi = phi[i].at(ip) - phi[j].at(jp);
            if (dphi < -1 * TMath::Pi())
              dphi += TMath::Pi();
            if (dphi > 1 * TMath::Pi())
              dphi -= TMath::Pi();
            hdphi[i][j]->Fill(dphi);
          } // jp
        } // ip
      } // j
    } // i

  } // ievent

  //==========================================================================//
  // CALCULATE C2
  //==========================================================================//
  cout << endl;
  cout << "--> Calculate c2 values" << endl;

  for (int i = 0; i < NDET; i++)
  {
    for (int j = 0; j < NDET; j++)
    {
      c2[i][j] = calc_cn(hdphi[i][j], 2);
      cout << " c2(" << ldet[i] << "-" << ldet[j] << ") = "
           << c2[i][j].first << " +/- " << c2[i][j].second
           << endl;
    } // j
  } // i

  //==========================================================================//
  // CALCULATE v2
  //==========================================================================//
  cout << endl;
  cout << "--> Calculating v2" << endl;

  gv2_3sub = new TGraphErrors();
  gv2_3sub->SetName("gv2_3sub");
  gv2_3sub->SetTitle(";;v_{2}{2}");
  gv2_3sub->SetMarkerStyle(20);
  gv2_3sub->SetMarkerColor(kBlue);
  gv2_3sub->SetLineColor(kBlue);

  for (int iv2 = 0; iv2 < NV2; iv2++)
  {
    v2[iv2] = vn_3sub(c2[dv2[iv2][0]][dv2[iv2][1]],
                      c2[dv2[iv2][0]][dv2[iv2][2]],
                      c2[dv2[iv2][1]][dv2[iv2][2]]);
    cout << "    calc:"
         << " vn_3sub("
         << "c2[" << dv2[iv2][0] << "][" << dv2[iv2][1] << "]" << ","
         << "c2[" << dv2[iv2][0] << "][" << dv2[iv2][2] << "]" << ","
         << "c2[" << dv2[iv2][1] << "][" << dv2[iv2][2] << "]" << ")"
         << endl;

    cout << " v_2{2}("
         << ldet[dv2[iv2][0]] << "-"
         << ldet[dv2[iv2][1]] << "-"
         << ldet[dv2[iv2][2]] << ") = "
         << v2[iv2].first
         << " +/- " << v2[iv2].second
         << endl;

    gv2_3sub->SetPoint(iv2, iv2, v2[iv2].first);
    gv2_3sub->SetPointError(iv2, 0, v2[iv2].second);
  } // iv2

  //==========================================================================//
  // PLOT OBJECTS
  //==========================================================================//
  cout << endl;
  cout << "--> Plotting" << endl;

  // data TGraph
  gv2_data = new TGraphErrors();
  gv2_data->SetName("gv2_data");
  gv2_data->SetTitle(";;v_{2}{2}");
  gv2_data->SetMarkerStyle(21);
  gv2_data->SetMarkerColor(kGreen+2);
  gv2_data->SetLineColor(kGreen+2);
  for (int i = 0; i < 5; i++)
  {
    gv2_data->SetPoint(i, i, v2data[i][0]);
    gv2_data->SetPointError(i, 0, v2data[i][1]);
  }

  //axis
  TH1D* haxis = new TH1D("haxis", ";;v_{2}{2}", 5, -0.5, 4.5);
  haxis->SetMinimum(0.0);
  haxis->SetMaximum(0.14);
  haxis->GetYaxis()->SetTitleFont(63);
  haxis->GetYaxis()->SetTitleSize(20);
  haxis->GetYaxis()->SetTitleOffset(1.5);
  haxis->GetYaxis()->SetLabelFont(63);
  haxis->GetYaxis()->SetLabelSize(14);
  haxis->GetXaxis()->SetTitleFont(63);
  haxis->GetXaxis()->SetTitleSize(20);
  haxis->GetXaxis()->SetTitleOffset(1.7);
  haxis->GetXaxis()->SetLabelFont(63);
  haxis->GetXaxis()->SetLabelSize(14);


  char alabel[500];
  for (int ibin = 1; ibin <= haxis->GetNbinsX(); ibin++)
  {
    sprintf(alabel, "%s-%s-%s", 
            ldet[dv2[ibin-1][0]], ldet[dv2[ibin-1][1]], ldet[dv2[ibin-1][2]]);
    haxis->GetXaxis()->SetBinLabel(ibin, alabel);
  }

  //legend
  TLegend *legv2 = new TLegend(0.5, 0.85, 0.9, 0.9);
  legv2->SetFillStyle(0);
  // legv2->SetBorderSize(0);
  legv2->SetTextFont(63);
  legv2->SetTextSize(14);
  legv2->SetNColumns(2);
  legv2->AddEntry(gv2_data, "Data", "P");
  legv2->AddEntry(gv2_3sub, "Toy MC", "P");


  //other
  TLine lv2;
  lv2.SetLineStyle(2);
  lv2.SetLineColor(kRed);

  TLatex t;
  t.SetNDC();
  t.SetTextFont(63);
  t.SetTextSize(14);

  //==========================================================================//
  // PLOT
  //==========================================================================//

  // TCanvas *cphi = new TCanvas("cphi", "phi", 800, 800);
  // cphi->cd(1);
  // hphi->Scale(1. / hphi->Integral());
  // hphi->Draw();
  // fphi->SetParameters(v2true[2], 0);
  // fphi->Draw("same");
  // ValErr c2test = calc_cn(hphi, 2);
  // cout << " c2: " << c2test.first << " +/- " << c2test.second << endl;



  TCanvas *cv2 = new TCanvas("cv2", "v2", 800, 600);

  cv2->cd(1);
  haxis->Draw();

  gv2_3sub->Draw("P");
  gv2_data->Draw("P");

  lv2.DrawLine(-0.5, v2true[2], 4.5, v2true[2]);

  legv2->Draw("same");

  t.SetTextColor(kBlack);
  t.DrawLatex(0.15, 0.85, Form("N_{event}=%i", Nevent));

  t.SetTextColor(kBlack);
  t.DrawLatex(0.15, 0.80, Form("#sigma_{#Delta}=%.3f", sigd));

  t.SetTextColor(kRed);
  t.DrawLatex(0.15, 0.75, Form("v_{2}^{true}=%.3f", v2true[2]));



  //==========================================================================//
  // PRINT
  //==========================================================================//
  cout << endl;
  cout << "--> Printing plot" << endl;

  cv2->Print("pdfs/ep_decorrelation.pdf");


}
///-----------------------------------------------------------------------------


///-----------------------------------------------------------------------------
///
/// Calculate v_n using the 3 sub-event method
///
ValErr vn_3sub(const ValErr &C_AB, const ValErr &C_AC, const ValErr &C_BC)
{
  // return values;
  double vn = 0;
  double err = 0;

  // Calculate the value
  vn = TMath::Sqrt( (C_AB.first * C_AC.first) / C_BC.first );

  // Calculate the uncertainty using error propogation
  double dC_AB = +0.5 * (C_AC.first / C_BC.first) / vn;
  double dC_AC = +0.5 * (C_AB.first / C_BC.first) / vn;
  double dC_BC = -0.5 * (C_AB.first * C_AC.first) / (C_BC.first * C_BC.first) / vn;

  err += TMath::Power(dC_AB * C_AB.second, 2);
  err += TMath::Power(dC_AC * C_AC.second, 2);
  err += TMath::Power(dC_BC * C_BC.second, 2);
  err = TMath::Sqrt(err);

  // check for nan's (imaginary numbers)
  if (vn != vn)
  {
    vn = 0;
    err = 0;
  }

  return make_pair(vn, err);

}
///-----------------------------------------------------------------------------



///-----------------------------------------------------------------------------
///
/// Calculate C_n coefficients from the histogram numerically
///
ValErr calc_cn(TH1D* hcorr, int order)
{
  if (order < 0 || order > 5)
  {
    cout << "WARNING!! calc_cn() - order parameter " << order
         << " is out of range, returnning 0's" << endl;
    return make_pair(0, 0);
  }
  if (!hcorr)
  {
    cout << "ERROR!! calc_cn() - hcorr not a valid pointer! Returning 0's!" << endl;
    return make_pair(0, 0);
  }

  double val = 0;
  double err = 0;
  double tot = 0;
  for (int ibin = 1; ibin <= hcorr->GetNbinsX(); ibin++)
  {
    double dphi = hcorr->GetBinCenter(ibin);
    double bc = hcorr->GetBinContent(ibin);
    double e = hcorr->GetBinError(ibin);

    val += bc * TMath::Cos((float)order * dphi);
    err += TMath::Power(e * TMath::Cos((float)order * dphi), 2);
    tot += bc;

  }
  val = val / tot;
  err = TMath::Sqrt(err) / tot;

  return make_pair(val, err);
}
///-----------------------------------------------------------------------------

