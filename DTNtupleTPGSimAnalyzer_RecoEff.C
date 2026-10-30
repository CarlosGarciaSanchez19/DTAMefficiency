#include "DTNtupleTPGSimAnalyzer.h"
#include "TVector2.h"
#include "TF1.h"
#include "DTAcceptances.h"
#include <cmath>

constexpr float PI_F = 3.1415927f;
constexpr double PI = 3.141592653589793;

DTNtupleTPGSimAnalyzer::DTNtupleTPGSimAnalyzer(const TString & inFileName,
  const TString & outFileName,
  const TString & quality = "",
  const bool & acc = false,
  const bool & timeCut = false
):
m_outFile(outFileName,"RECREATE"), DTNtupleBaseAnalyzer(inFileName), quality_(quality), acc_(acc), timeCut_(timeCut)
{ 
  m_minMuPt = 20; //DM jp 0 pg 20 
  m_maxDPhi = 0.15; //DM 0.1 or 0.2 (m_maxMuTrigDPhi)? // CGS 0.15
  m_maxMuSegDPhi = 0.1;
  m_maxMuSegDEta = 0.15; // Condition comented out below. Unccomment if needed
  m_minSegHits = 4; // CGS later modified depending on quality
  m_minZSegHits = 0; // Condition comented out below. Unccomment if needed
  m_maxSegT0 = 12.5; // CGS 999 no cut on t0
  Tag = "RECO"; // chamb + algo (Mostly so plotting script works a bit easier)
}


DTNtupleTPGSimAnalyzer::~DTNtupleTPGSimAnalyzer()
{
  
}



void DTNtupleTPGSimAnalyzer::Loop()
{
  book();
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  //for (Long64_t jentry = 0; jentry < 100; jentry++)
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    //if (jentry >= 10000) continue;
    thisEntry = jentry;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEvent(jentry);   nbytes += nb;
    
    //cout << jentry << endl;
    if(jentry % 100 == 0)
    std::cout << "[DTNtupleTPGSimAnalyzer::Loop] processed : "
    << jentry << " entries\r" << std::flush;
    fill();
  }
  
  std::cout << std::endl;
  
  endJob();
}



void DTNtupleTPGSimAnalyzer::book()
{
  m_outFile.cd();
 
  std::vector<std::string> Type = {"matched", "total"};
  std::vector<std::string> chambTags = {"MB1", "MB2", "MB3", "MB4"};
  
  for (const auto & it : Type)
  {
    for (const auto & chamb : chambTags)
    {
       m_plots["EffPhiEta_" + chamb + "_" + Tag + "_" + it] = new TH2D(("hEffPhiEta_" + chamb + "_" + Tag + "_" + it).c_str(),
         "Efficiency vs Phi vs Eta; #phi; #eta",
       150, -3.2, +3.2, 200, -1.2, +1.2); // 64, -3.2, +3.2, 75, -1.5, +1.5
       m_plots["EffEta_" + chamb + "_" + Tag + "_" + it] = new TH1D(("hEffEta_" + chamb + "_" + Tag + "_" + it).c_str(),
          ("Efficiency vs Eta for " + chamb + " " + Tag + "; #eta; Efficiency").c_str(),
        80, -1.2, +1.2);
        m_plots["EffPhi_" + chamb + "_" + Tag + "_" + it] = new TH1D(("hEffPhi_" + chamb + "_" + Tag + "_" + it).c_str(),
          ("Efficiency vs Phi for " + chamb + " " + Tag + "; #phi; Efficiency").c_str(),
        32, -3.2, +3.2);
        std::vector<double> pt_bins = {20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 100.0, 120.0, 160.0};
        m_plots["EffPt_" + chamb + "_" + Tag + "_" + it] = new TH1D(("hEffPt_" + chamb + "_" + Tag + "_" + it).c_str(),
          ("Efficiency vs Pt for " + chamb + " " + Tag + "; pt; Efficiency").c_str(),
	      pt_bins.size()-1, pt_bins.data()); //DM pg 80, 20, 200);  jp 160, 0, 200); dy6 40, 0, 1000);
        m_plots["EffPhiEta_" + chamb + "_" + Tag + "_" + it]->Sumw2();
        m_plots["EffEta_" + chamb + "_" + Tag + "_" + it]->Sumw2();
        m_plots["EffPhi_" + chamb + "_" + Tag + "_" + it]->Sumw2();
        m_plots["EffPt_" + chamb + "_" + Tag + "_" + it]->Sumw2();
    }
    if (!acc_)
      m_plots["hEffPhiEta_COMBINED_" + Tag + "_" + it] = new TH2D(("hEffPhiEta_COMBINED_" + Tag + "_" + it).c_str(),
           "Efficiency vs Phi vs Eta; #phi; #eta", 150, -3.2, +3.2, 200, -1.2, +1.2);
  }
}

Double_t DTNtupleTPGSimAnalyzer::bend_phi(Double_t gen_phi, Double_t gen_pt, int pid, int st)
{
  Double_t a = 0.162; // station 1, 20-200 GeV, PG
  Double_t b = -0.032;
  Double_t c = 0.016;

  if (st == 2) a += 0.01;
  else if (st == 3) a += 0.02;
  else if (st == 4) a += 0.03;

  Double_t Dphi = a * exp(b * gen_pt) + c;

  if (pid == 13) gen_phi += Dphi; // muon
  else gen_phi -= Dphi; // anti-muon

  return gen_phi;
}


bool DTNtupleTPGSimAnalyzer::is_in_accmap(std::string chamb, Double_t genPhi, Double_t genEta)
{
  std::vector<std::string> whTags = {"whNeg2", "whNeg1", "wh0", "wh1", "wh2"};
  std::vector<std::string> secTags = {"sec1", "sec2", "sec3", "sec4", "sec5", "sec6", "sec7", "sec8", "sec9", "sec10", "sec11", "sec12", "sec13", "sec14"};
  if (chamb == "MB1")
  {
    for (const auto & wh : whTags)
    {
      for (const auto & sec : secTags)
      {
        if (sec == "sec13" || sec == "sec14") continue;
        Double_t phi1 = MB1_phi1[wh + "_" + sec];
        Double_t phi2 = MB1_phi2[wh + "_" + sec];
        Double_t eta1 = MB1_eta1[wh + "_" + sec];
        Double_t eta2 = MB1_eta2[wh + "_" + sec];
        if (sec == "sec7" && genPhi <= phi2) genPhi += 2*PI;
        if (sec == "sec7") phi2 += 2*PI;
        if ((genPhi <= phi2 && genPhi >= phi1)
          && (genEta <= eta2 && genEta >= eta1))
          return true;
      }
    }
  }
  if (chamb == "MB2")
  {
    for (const auto & wh : whTags)
    {
      for (const auto & sec : secTags)
      {
        if (sec == "sec13" || sec == "sec14") continue;
        Double_t phi1 = MB2_phi1[wh + "_" + sec];
        Double_t phi2 = MB2_phi2[wh + "_" + sec];
        Double_t eta1 = MB2_eta1[wh + "_" + sec];
        Double_t eta2 = MB2_eta2[wh + "_" + sec];
        if (sec == "sec7" && genPhi <= phi2) genPhi += 2*PI;
        if (sec == "sec7") phi2 += 2*PI;
        if ((genPhi <= phi2 && genPhi >= phi1)
          && (genEta <= eta2 && genEta >= eta1))
          return true;
      }
    }
  }
  if (chamb == "MB3")
  {
    for (const auto & wh : whTags)
    {
      for (const auto & sec : secTags)
      {
        if (sec == "sec13" || sec == "sec14") continue;
        Double_t phi1 = MB3_phi1[wh + "_" + sec];
        Double_t phi2 = MB3_phi2[wh + "_" + sec];
        Double_t eta1 = MB3_eta1[wh + "_" + sec];
        Double_t eta2 = MB3_eta2[wh + "_" + sec];
        if (sec == "sec7" && genPhi <= phi2) genPhi += 2*PI;
        if (sec == "sec7") phi2 += 2*PI;
        if ((genPhi <= phi2 && genPhi >= phi1)
          && (genEta <= eta2 && genEta >= eta1))
          return true;
      }
    }
  }
  if (chamb == "MB4")
  {
    for (const auto & wh : whTags)
    {
      for (const auto & sec : secTags)
      {
        Double_t phi1 = MB4_phi1[wh + "_" + sec];
        Double_t phi2 = MB4_phi2[wh + "_" + sec];
        Double_t eta1 = MB4_eta1[wh + "_" + sec];
        Double_t eta2 = MB4_eta2[wh + "_" + sec];
        if (sec == "sec7" && genPhi <= phi2) genPhi += 2*PI;
        if (sec == "sec7") phi2 += 2*PI;
        if ((genPhi <= phi2 && genPhi >= phi1)
          && (genEta <= eta2 && genEta >= eta1))
          return true;
      }
    }
  }
  return false;
}

void DTNtupleTPGSimAnalyzer::fill()
{
  std::vector<std::string> chambTags = {"MB1", "MB2", "MB3", "MB4"};

  for (std::size_t iGenPart = 0; iGenPart < gen_nGenParts; ++iGenPart)
  {
    Double_t genEta = gen_eta->at(iGenPart);
    Double_t genPhi = gen_phi->at(iGenPart);
    Double_t genPt = gen_pt->at(iGenPart);

    if (std::abs(gen_pdgId->at(iGenPart)) != 13 || genPt < m_minMuPt || abs(genEta) > 1.2) continue;

    std::vector<Int_t> genMuSts = {-1, -1, -1, -1};
    std::vector<Double_t> genBendedPhis = {-1, -1, -1, -1};

    if (!acc_)
      m_plots["hEffPhiEta_COMBINED_" + Tag + "_total"]->Fill(genPhi, genEta);
    for (std::size_t st = 1; st <= chambTags.size(); st++)
    {
      const auto & chamb = "MB" + std::to_string(st);
      genBendedPhis[st - 1] = bend_phi(genPhi, genPt, gen_pdgId->at(iGenPart), st);
      // genBendedPhis[st - 1] = genPhi;
      // if (chamb == "MB4" && abs(genEta) >= 0.75 && acc_) continue;
      if (!is_in_accmap(chamb, genBendedPhis[st - 1], genEta) && acc_) continue;// Aplying acceptance
      m_plots["EffPhiEta_" + chamb + "_" + Tag + "_total"]->Fill(genBendedPhis[st - 1], genEta);
      m_plots["EffEta_" + chamb + "_" + Tag + "_total"]->Fill(genEta);
      m_plots["EffPhi_" + chamb + "_" + Tag + "_total"]->Fill(genBendedPhis[st - 1]);
      m_plots["EffPt_" + chamb + "_" + Tag + "_total"]->Fill(genPt);
      if (acc_) genMuSts[st - 1] = (Int_t)st;
    }

    std::vector<Int_t> bestRMsegs = {-1, -1, -1, -1}; // RM1, RM2, RM3, RM4
    std::vector<Int_t> bestSegNHits = {0, 0, 0, 0};
    // Int_t RM1 = 0;
    // Int_t RM2 = 0;
    // Int_t RM3 = 0;
    // Int_t RM4 = 0;
    // Int_t segSt = -1;

    int minQuality = -99;  //DM took this quality approach from old eff class
    bool qualityCorrelated = false; 

    if      (quality_ == "nothreehits"){ m_minSegHits = 4; }
    else if (quality_ == "correlated"){  m_minSegHits = 6;}
    else if (quality_ == "All" ){        m_minSegHits = 3; }
    else { 
      cout << "Error: Efficiency quality category not found" << endl; 
    }

    for (std::size_t iSeg = 0; iSeg < seg_nSegments; ++iSeg)
    {
      Int_t segNHits = seg_phi_nHits->at(iSeg);
      Int_t segZNHits = seg_z_nHits->at(iSeg);
      Int_t segSt = seg_station->at(iSeg);

      if (abs(seg_phi_t0->at(iSeg)) > m_maxSegT0 && timeCut_) continue; // t0 cut

      Double_t segPhi = seg_posGlb_phi->at(iSeg);
      Double_t segEta = seg_posGlb_eta->at(iSeg);

      // std::cout << "genEta: " << genEta << " genPhi: " << genPhi << std::endl;

      Double_t muSegDPhi = std::abs(acos(cos(genBendedPhis[segSt - 1] - segPhi)));
      // Double_t muSegDPhi = std::abs(genPhi - seg_posGlb_phi->at(iSeg));
      // Double_t muSegDEta = std::abs(genEta - segEta);
      if (muSegDPhi < m_maxMuSegDPhi && segNHits >= m_minSegHits) bestRMsegs[segSt - 1] = segSt;
      //   // (segZNHits >= m_minZSegHits  || segSt==4) &&
      //   segNHits >= bestSegNHits.at(segSt - 1))
      // {
      //   bestSegNHits[segSt - 1] = segNHits;
      //   bestRMsegs[segSt - 1] = segSt;
      // }
      // if (muSegDPhi >= m_maxMuSegDPhi && segNHits < m_minSegHits) continue;

      // bestRMsegs[segSt - 1] = segSt;
    }

    if (std::count_if(bestRMsegs.begin(), bestRMsegs.end(), [](int x){ return x != -1; }) > 1 && !acc_)
      m_plots["hEffPhiEta_COMBINED_" + Tag + "_matched"]->Fill(genPhi, genEta);
    for (std::size_t i = 0; i < chambTags.size(); ++i)
    {
      if (bestRMsegs[i] == -1) continue; // No RM for this station
      if (bestRMsegs[i] != genMuSts[i] && acc_) continue;
      const auto & chamb = "MB" + std::to_string(i + 1);
      m_plots["EffPhiEta_" + chamb + "_" + Tag + "_matched"]->Fill(genBendedPhis[i], genEta);
      m_plots["EffEta_" + chamb + "_" + Tag + "_matched"]->Fill(genEta);
      m_plots["EffPhi_" + chamb + "_" + Tag + "_matched"]->Fill(genBendedPhis[i]);
      m_plots["EffPt_" + chamb + "_" + Tag + "_matched"]->Fill(genPt);
    }
  }
}
 

void DTNtupleTPGSimAnalyzer::endJob()
{
  m_outFile.cd();
  
  m_outFile.Write();
  m_outFile.Close();
}
