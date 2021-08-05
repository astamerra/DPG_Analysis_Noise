//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  7 11:14:58 2021 by ROOT version 6.22/00
// from TTree MuDPGTree/Mu DPG Tree
// found on file: MuDPGNTuple_MWGR3_2021_EXP_run341343_99.root
//////////////////////////////////////////////////////////

#ifndef Analysis_Noisy_Strip_h
#define Analysis_Noisy_Strip_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using namespace std;

class Analysis_Noisy_Strip {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TString 	    fileName; // output filename set in Analysis.cpp
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event_runNumber;
   Int_t           event_lumiBlock;
   Long64_t        event_eventNumber;
   UInt_t          gemDigi_nDigis;
   vector<short>   *gemDigi_station;
   vector<short>   *gemDigi_region;
   vector<short>   *gemDigi_roll;
   vector<short>   *gemDigi_bx;
   vector<short>   *gemDigi_strip;
   vector<short>   *gemDigi_pad;
   vector<float>   *gemDigi_g_r;
   vector<float>   *gemDigi_g_phi;
   vector<float>   *gemDigig_eta;
   vector<float>   *gemDigi_g_x;
   vector<float>   *gemDigi_g_y;
   vector<float>   *gemDigi_g_z;
   UInt_t          gemRecHit_nRecHits;
   vector<int>     *gemRecHit_cluster_size;
   vector<int>     *gemRecHit_firstClusterStrip;
   vector<int>     *gemRecHit_bx;
   vector<int>     *gemRecHit_region;
   vector<int>     *gemRecHit_chamber;
   vector<int>     *gemRecHit_layer;
   vector<int>     *gemRecHit_etaPartition;
   vector<float>   *gemRecHit_loc_r;
   vector<float>   *gemRecHit_loc_phi;
   vector<float>   *gemRecHit_loc_x;
   vector<float>   *gemRecHit_loc_y;
   vector<float>   *gemRecHit_loc_z;
   vector<float>   *gemRecHit_loc_errX;
   vector<float>   *gemRecHit_loc_errY;
   vector<float>   *gemRecHit_g_r;
   vector<float>   *gemRecHit_g_phi;
   vector<float>   *gemRecHit_g_x;
   vector<float>   *gemRecHit_g_y;
   vector<float>   *gemRecHit_g_z;
   UInt_t          gemSegment_nSegments;
   vector<short>   *gemSegment_region;
   vector<short>   *gemSegment_ring;
   vector<short>   *gemSegment_station;
   vector<float>   *gemSegment_posLoc_x;
   vector<float>   *gemSegment_posLoc_y;
   vector<float>   *gemSegment_posLoc_z;
   vector<float>   *gemSegment_dirLoc_x;
   vector<float>   *gemSegment_dirLoc_y;
   vector<float>   *gemSegment_dirLoc_z;
   vector<float>   *gemSegment_posGlb_x;
   vector<float>   *gemSegment_posGlb_y;
   vector<float>   *gemSegment_posGlb_z;
   vector<float>   *gemSegment_time;
   vector<float>   *gemSegment_time_err;
   vector<double>  *gemSegment_chi2;
   vector<float>   *gemSegment_posGlb_phi;
   vector<float>   *gemSegment_posGlb_eta;
   vector<float>   *gemSegment_dirGlb_phi;
   vector<float>   *gemSegment_dirGlb_eta;
   UInt_t          mu_nMuons;
   vector<float>   *mu_pt;
   vector<float>   *mu_phi;
   vector<float>   *mu_eta;
   vector<short>   *mu_charge;
   vector<bool>    *mu_isGlobal;
   vector<bool>    *mu_isStandalone;
   vector<bool>    *mu_isTracker;
   vector<bool>    *mu_isGEM;
   vector<bool>    *mu_isCSC;
   vector<bool>    *mu_isME11;
   vector<bool>    *mu_isLoose;
   vector<bool>    *mu_isMedium;
   vector<bool>    *mu_isTight;
   vector<bool>    *mu_propagated_isME11;
   vector<float>   *mu_propagated_TrackNormChi2;
   vector<float>   *mu_propagated_numberOfValidPixelHits;
   vector<float>   *mu_propagated_innerTracker_ValidFraction;
   vector<float>   *mu_propagated_numberOfValidTrackerHits;
   Float_t         mu_path_length;
   vector<bool>    *mu_isinsideout;
   vector<bool>    *mu_isincoming;
   vector<int>     *mu_propagated_region;
   vector<int>     *mu_propagated_layer;
   vector<int>     *mu_propagated_chamber;
   vector<int>     *mu_propagated_etaP;
   vector<float>   *mu_propagated_pt;
   vector<float>   *mu_propagated_phi;
   vector<float>   *mu_propagated_eta;
   vector<float>   *mu_propagated_charge;
   vector<float>   *mu_propagatedLoc_x;
   vector<float>   *mu_propagatedLoc_y;
   vector<float>   *mu_propagatedLoc_z;
   vector<float>   *mu_propagatedLoc_r;
   vector<float>   *mu_propagatedLoc_phi;
   vector<float>   *mu_propagatedLoc_errX;
   vector<float>   *mu_propagatedLoc_errY;
   vector<float>   *mu_propagatedLoc_dirX;
   vector<float>   *mu_propagatedGlb_x;
   vector<float>   *mu_propagatedGlb_y;
   vector<float>   *mu_propagatedGlb_z;
   vector<float>   *mu_propagatedGlb_r;
   vector<float>   *mu_propagatedGlb_phi;
   vector<float>   *mu_propagatedGlb_errX;
   vector<float>   *mu_propagatedGlb_errY;
   vector<float>   *mu_propagatedGlb_errR;
   vector<float>   *mu_propagatedGlb_errPhi;
   vector<float>   *mu_propagated_EtaPartition_centerX;
   vector<float>   *mu_propagated_EtaPartition_centerY;
   vector<float>   *mu_propagated_EtaPartition_rMax;
   vector<float>   *mu_propagated_EtaPartition_rMin;
   vector<float>   *mu_propagated_EtaPartition_phiMax;
   vector<float>   *mu_propagated_EtaPartition_phiMin;
   vector<float>   *mu_propagated_Innermost_x;
   vector<float>   *mu_propagated_Innermost_y;
   vector<float>   *mu_propagated_Innermost_z;
   vector<float>   *mu_propagated_Outermost_x;
   vector<float>   *mu_propagated_Outermost_y;
   vector<float>   *mu_propagated_Outermost_z;

   // List of branches
   TBranch        *b_event_runNumber;   //!
   TBranch        *b_event_lumiBlock;   //!
   TBranch        *b_event_eventNumber;   //!
   TBranch        *b_gemDigi_nDigis;   //!
   TBranch        *b_gemDigi_station;   //!
   TBranch        *b_gemDigi_region;   //!
   TBranch        *b_gemDigi_roll;   //!
   TBranch        *b_gemDigi_bx;   //!
   TBranch        *b_gemDigi_strip;   //!
   TBranch        *b_gemDigi_pad;   //!
   TBranch        *b_gemDigi_g_r;   //!
   TBranch        *b_gemDigi_g_phi;   //!
   TBranch        *b_gemDigig_eta;   //!
   TBranch        *b_gemDigi_g_x;   //!
   TBranch        *b_gemDigi_g_y;   //!
   TBranch        *b_gemDigi_g_z;   //!
   TBranch        *b_gemRecHit_nRecHits;   //!
   TBranch        *b_gemRecHit_cluster_size;   //!
   TBranch        *b_gemRecHit_firstClusterStrip;   //!
   TBranch        *b_gemRecHit_bx;   //!
   TBranch        *b_gemRecHit_region;   //!
   TBranch        *b_gemRecHit_chamber;   //!
   TBranch        *b_gemRecHit_layer;   //!
   TBranch        *b_gemRecHit_etaPartition;   //!
   TBranch        *b_gemRecHit_loc_r;   //!
   TBranch        *b_gemRecHit_loc_phi;   //!
   TBranch        *b_gemRecHit_loc_x;   //!
   TBranch        *b_gemRecHit_loc_y;   //!
   TBranch        *b_gemRecHit_loc_z;   //!
   TBranch        *b_gemRecHit_loc_errX;   //!
   TBranch        *b_gemRecHit_loc_errY;   //!
   TBranch        *b_gemRecHit_g_r;   //!
   TBranch        *b_gemRecHit_g_phi;   //!
   TBranch        *b_gemRecHit_g_x;   //!
   TBranch        *b_gemRecHit_g_y;   //!
   TBranch        *b_gemRecHit_g_z;   //!
   TBranch        *b_gemSegment_nSegments;   //!
   TBranch        *b_gemSegment_region;   //!
   TBranch        *b_gemSegment_ring;   //!
   TBranch        *b_gemSegment_station;   //!
   TBranch        *b_gemSegment_posLoc_x;   //!
   TBranch        *b_gemSegment_posLoc_y;   //!
   TBranch        *b_gemSegment_posLoc_z;   //!
   TBranch        *b_gemSegment_dirLoc_x;   //!
   TBranch        *b_gemSegment_dirLoc_y;   //!
   TBranch        *b_gemSegment_dirLoc_z;   //!
   TBranch        *b_gemSegment_posGlb_x;   //!
   TBranch        *b_gemSegment_posGlb_y;   //!
   TBranch        *b_gemSegment_posGlb_z;   //!
   TBranch        *b_gemSegment_time;   //!
   TBranch        *b_gemSegment_time_err;   //!
   TBranch        *b_gemSegment_chi2;   //!
   TBranch        *b_gemSegment_posGlb_phi;   //!
   TBranch        *b_gemSegment_posGlb_eta;   //!
   TBranch        *b_gemSegment_dirGlb_phi;   //!
   TBranch        *b_gemSegment_dirGlb_eta;   //!
   TBranch        *b_mu_nMuons;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_isGlobal;   //!
   TBranch        *b_mu_isStandalone;   //!
   TBranch        *b_mu_isTracker;   //!
   TBranch        *b_mu_isGEM;   //!
   TBranch        *b_mu_isCSC;   //!
   TBranch        *b_mu_isME11;   //!
   TBranch        *b_mu_isLoose;   //!
   TBranch        *b_mu_isMedium;   //!
   TBranch        *b_mu_isTight;   //!
   TBranch        *b_mu_propagated_isME11;   //!
   TBranch        *b_mu_propagated_TrackNormChi2;   //!
   TBranch        *b_mu_propagated_numberOfValidPixelHits;   //!
   TBranch        *b_mu_propagated_innerTracker_ValidFraction;   //!
   TBranch        *b_mu_propagated_numberOfValidTrackerHits;   //!
   TBranch        *b_mu_path_length;   //!
   TBranch        *b_mu_isinsideout;   //!
   TBranch        *b_mu_isincoming;   //!
   TBranch        *b_mu_propagated_region;   //!
   TBranch        *b_mu_propagated_layer;   //!
   TBranch        *b_mu_propagated_chamber;   //!
   TBranch        *b_mu_propagated_etaP;   //!
   TBranch        *b_mu_propagated_pt;   //!
   TBranch        *b_mu_propagated_phi;   //!
   TBranch        *b_mu_propagated_eta;   //!
   TBranch        *b_mu_propagated_charge;   //!
   TBranch        *b_mu_propagatedLoc_x;   //!
   TBranch        *b_mu_propagatedLoc_y;   //!
   TBranch        *b_mu_propagatedLoc_z;   //!
   TBranch        *b_mu_propagatedLoc_r;   //!
   TBranch        *b_mu_propagatedLoc_phi;   //!
   TBranch        *b_mu_propagatedLoc_errX;   //!
   TBranch        *b_mu_propagatedLoc_errY;   //!
   TBranch        *b_mu_propagatedLoc_dirX;   //!
   TBranch        *b_mu_propagatedGlb_x;   //!
   TBranch        *b_mu_propagatedGlb_y;   //!
   TBranch        *b_mu_propagatedGlb_z;   //!
   TBranch        *b_mu_propagatedGlb_r;   //!
   TBranch        *b_mu_propagatedGlb_phi;   //!
   TBranch        *b_mu_propagatedGlb_errX;   //!
   TBranch        *b_mu_propagatedGlb_errY;   //!
   TBranch        *b_mu_propagatedGlb_errR;   //!
   TBranch        *b_mu_propagatedGlb_errPhi;   //!
   TBranch        *b_mu_propagated_EtaPartition_centerX;   //!
   TBranch        *b_mu_propagated_EtaPartition_centerY;   //!
   TBranch        *b_mu_propagated_EtaPartition_rMax;   //!
   TBranch        *b_mu_propagated_EtaPartition_rMin;   //!
   TBranch        *b_mu_propagated_EtaPartition_phiMax;   //!
   TBranch        *b_mu_propagated_EtaPartition_phiMin;   //!
   TBranch        *b_mu_propagated_Innermost_x;   //!
   TBranch        *b_mu_propagated_Innermost_y;   //!
   TBranch        *b_mu_propagated_Innermost_z;   //!
   TBranch        *b_mu_propagated_Outermost_x;   //!
   TBranch        *b_mu_propagated_Outermost_y;   //!
   TBranch        *b_mu_propagated_Outermost_z;   //!

   Analysis_Noisy_Strip(TTree *tree, TString fname);
   virtual ~Analysis_Noisy_Strip();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analysis_Noisy_Strip_cxx
Analysis_Noisy_Strip::Analysis_Noisy_Strip(TTree *tree, TString fname) : fChain(0) 
{
   Init(tree);
   fileName = fname;
}

Analysis_Noisy_Strip::~Analysis_Noisy_Strip()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis_Noisy_Strip::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis_Noisy_Strip::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analysis_Noisy_Strip::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   gemDigi_station = 0;
   gemDigi_region = 0;
   gemDigi_roll = 0;
   gemDigi_bx = 0;
   gemDigi_strip = 0;
   gemDigi_pad = 0;
   gemDigi_g_r = 0;
   gemDigi_g_phi = 0;
   gemDigig_eta = 0;
   gemDigi_g_x = 0;
   gemDigi_g_y = 0;
   gemDigi_g_z = 0;
   gemRecHit_cluster_size = 0;
   gemRecHit_firstClusterStrip = 0;
   gemRecHit_bx = 0;
   gemRecHit_region = 0;
   gemRecHit_chamber = 0;
   gemRecHit_layer = 0;
   gemRecHit_etaPartition = 0;
   gemRecHit_loc_r = 0;
   gemRecHit_loc_phi = 0;
   gemRecHit_loc_x = 0;
   gemRecHit_loc_y = 0;
   gemRecHit_loc_z = 0;
   gemRecHit_loc_errX = 0;
   gemRecHit_loc_errY = 0;
   gemRecHit_g_r = 0;
   gemRecHit_g_phi = 0;
   gemRecHit_g_x = 0;
   gemRecHit_g_y = 0;
   gemRecHit_g_z = 0;
   gemSegment_region = 0;
   gemSegment_ring = 0;
   gemSegment_station = 0;
   gemSegment_posLoc_x = 0;
   gemSegment_posLoc_y = 0;
   gemSegment_posLoc_z = 0;
   gemSegment_dirLoc_x = 0;
   gemSegment_dirLoc_y = 0;
   gemSegment_dirLoc_z = 0;
   gemSegment_posGlb_x = 0;
   gemSegment_posGlb_y = 0;
   gemSegment_posGlb_z = 0;
   gemSegment_time = 0;
   gemSegment_time_err = 0;
   gemSegment_chi2 = 0;
   gemSegment_posGlb_phi = 0;
   gemSegment_posGlb_eta = 0;
   gemSegment_dirGlb_phi = 0;
   gemSegment_dirGlb_eta = 0;
   mu_pt = 0;
   mu_phi = 0;
   mu_eta = 0;
   mu_charge = 0;
   mu_isGlobal = 0;
   mu_isStandalone = 0;
   mu_isTracker = 0;
   mu_isGEM = 0;
   mu_isCSC = 0;
   mu_isME11 = 0;
   mu_isLoose = 0;
   mu_isMedium = 0;
   mu_isTight = 0;
   mu_propagated_isME11 = 0;
   mu_propagated_TrackNormChi2 = 0;
   mu_propagated_numberOfValidPixelHits = 0;
   mu_propagated_innerTracker_ValidFraction = 0;
   mu_propagated_numberOfValidTrackerHits = 0;
   mu_isinsideout = 0;
   mu_isincoming = 0;
   mu_propagated_region = 0;
   mu_propagated_layer = 0;
   mu_propagated_chamber = 0;
   mu_propagated_etaP = 0;
   mu_propagated_pt = 0;
   mu_propagated_phi = 0;
   mu_propagated_eta = 0;
   mu_propagated_charge = 0;
   mu_propagatedLoc_x = 0;
   mu_propagatedLoc_y = 0;
   mu_propagatedLoc_z = 0;
   mu_propagatedLoc_r = 0;
   mu_propagatedLoc_phi = 0;
   mu_propagatedLoc_errX = 0;
   mu_propagatedLoc_errY = 0;
   mu_propagatedLoc_dirX = 0;
   mu_propagatedGlb_x = 0;
   mu_propagatedGlb_y = 0;
   mu_propagatedGlb_z = 0;
   mu_propagatedGlb_r = 0;
   mu_propagatedGlb_phi = 0;
   mu_propagatedGlb_errX = 0;
   mu_propagatedGlb_errY = 0;
   mu_propagatedGlb_errR = 0;
   mu_propagatedGlb_errPhi = 0;
   mu_propagated_EtaPartition_centerX = 0;
   mu_propagated_EtaPartition_centerY = 0;
   mu_propagated_EtaPartition_rMax = 0;
   mu_propagated_EtaPartition_rMin = 0;
   mu_propagated_EtaPartition_phiMax = 0;
   mu_propagated_EtaPartition_phiMin = 0;
   mu_propagated_Innermost_x = 0;
   mu_propagated_Innermost_y = 0;
   mu_propagated_Innermost_z = 0;
   mu_propagated_Outermost_x = 0;
   mu_propagated_Outermost_y = 0;
   mu_propagated_Outermost_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_runNumber", &event_runNumber, &b_event_runNumber);
   fChain->SetBranchAddress("event_lumiBlock", &event_lumiBlock, &b_event_lumiBlock);
   fChain->SetBranchAddress("event_eventNumber", &event_eventNumber, &b_event_eventNumber);
   fChain->SetBranchAddress("gemDigi_nDigis", &gemDigi_nDigis, &b_gemDigi_nDigis);
   fChain->SetBranchAddress("gemDigi_station", &gemDigi_station, &b_gemDigi_station);
   fChain->SetBranchAddress("gemDigi_region", &gemDigi_region, &b_gemDigi_region);
   fChain->SetBranchAddress("gemDigi_roll", &gemDigi_roll, &b_gemDigi_roll);
   fChain->SetBranchAddress("gemDigi_bx", &gemDigi_bx, &b_gemDigi_bx);
   fChain->SetBranchAddress("gemDigi_strip", &gemDigi_strip, &b_gemDigi_strip);
   fChain->SetBranchAddress("gemDigi_pad", &gemDigi_pad, &b_gemDigi_pad);
   fChain->SetBranchAddress("gemDigi_g_r", &gemDigi_g_r, &b_gemDigi_g_r);
   fChain->SetBranchAddress("gemDigi_g_phi", &gemDigi_g_phi, &b_gemDigi_g_phi);
   fChain->SetBranchAddress("gemDigig_eta", &gemDigig_eta, &b_gemDigig_eta);
   fChain->SetBranchAddress("gemDigi_g_x", &gemDigi_g_x, &b_gemDigi_g_x);
   fChain->SetBranchAddress("gemDigi_g_y", &gemDigi_g_y, &b_gemDigi_g_y);
   fChain->SetBranchAddress("gemDigi_g_z", &gemDigi_g_z, &b_gemDigi_g_z);
   fChain->SetBranchAddress("gemRecHit_nRecHits", &gemRecHit_nRecHits, &b_gemRecHit_nRecHits);
   fChain->SetBranchAddress("gemRecHit_cluster_size", &gemRecHit_cluster_size, &b_gemRecHit_cluster_size);
   fChain->SetBranchAddress("gemRecHit_firstClusterStrip", &gemRecHit_firstClusterStrip, &b_gemRecHit_firstClusterStrip);
   fChain->SetBranchAddress("gemRecHit_bx", &gemRecHit_bx, &b_gemRecHit_bx);
   fChain->SetBranchAddress("gemRecHit_region", &gemRecHit_region, &b_gemRecHit_region);
   fChain->SetBranchAddress("gemRecHit_chamber", &gemRecHit_chamber, &b_gemRecHit_chamber);
   fChain->SetBranchAddress("gemRecHit_layer", &gemRecHit_layer, &b_gemRecHit_layer);
   fChain->SetBranchAddress("gemRecHit_etaPartition", &gemRecHit_etaPartition, &b_gemRecHit_etaPartition);
   fChain->SetBranchAddress("gemRecHit_loc_r", &gemRecHit_loc_r, &b_gemRecHit_loc_r);
   fChain->SetBranchAddress("gemRecHit_loc_phi", &gemRecHit_loc_phi, &b_gemRecHit_loc_phi);
   fChain->SetBranchAddress("gemRecHit_loc_x", &gemRecHit_loc_x, &b_gemRecHit_loc_x);
   fChain->SetBranchAddress("gemRecHit_loc_y", &gemRecHit_loc_y, &b_gemRecHit_loc_y);
   fChain->SetBranchAddress("gemRecHit_loc_z", &gemRecHit_loc_z, &b_gemRecHit_loc_z);
   fChain->SetBranchAddress("gemRecHit_loc_errX", &gemRecHit_loc_errX, &b_gemRecHit_loc_errX);
   fChain->SetBranchAddress("gemRecHit_loc_errY", &gemRecHit_loc_errY, &b_gemRecHit_loc_errY);
   fChain->SetBranchAddress("gemRecHit_g_r", &gemRecHit_g_r, &b_gemRecHit_g_r);
   fChain->SetBranchAddress("gemRecHit_g_phi", &gemRecHit_g_phi, &b_gemRecHit_g_phi);
   fChain->SetBranchAddress("gemRecHit_g_x", &gemRecHit_g_x, &b_gemRecHit_g_x);
   fChain->SetBranchAddress("gemRecHit_g_y", &gemRecHit_g_y, &b_gemRecHit_g_y);
   fChain->SetBranchAddress("gemRecHit_g_z", &gemRecHit_g_z, &b_gemRecHit_g_z);
   fChain->SetBranchAddress("gemSegment_nSegments", &gemSegment_nSegments, &b_gemSegment_nSegments);
   fChain->SetBranchAddress("gemSegment_region", &gemSegment_region, &b_gemSegment_region);
   fChain->SetBranchAddress("gemSegment_ring", &gemSegment_ring, &b_gemSegment_ring);
   fChain->SetBranchAddress("gemSegment_station", &gemSegment_station, &b_gemSegment_station);
   fChain->SetBranchAddress("gemSegment_posLoc_x", &gemSegment_posLoc_x, &b_gemSegment_posLoc_x);
   fChain->SetBranchAddress("gemSegment_posLoc_y", &gemSegment_posLoc_y, &b_gemSegment_posLoc_y);
   fChain->SetBranchAddress("gemSegment_posLoc_z", &gemSegment_posLoc_z, &b_gemSegment_posLoc_z);
   fChain->SetBranchAddress("gemSegment_dirLoc_x", &gemSegment_dirLoc_x, &b_gemSegment_dirLoc_x);
   fChain->SetBranchAddress("gemSegment_dirLoc_y", &gemSegment_dirLoc_y, &b_gemSegment_dirLoc_y);
   fChain->SetBranchAddress("gemSegment_dirLoc_z", &gemSegment_dirLoc_z, &b_gemSegment_dirLoc_z);
   fChain->SetBranchAddress("gemSegment_posGlb_x", &gemSegment_posGlb_x, &b_gemSegment_posGlb_x);
   fChain->SetBranchAddress("gemSegment_posGlb_y", &gemSegment_posGlb_y, &b_gemSegment_posGlb_y);
   fChain->SetBranchAddress("gemSegment_posGlb_z", &gemSegment_posGlb_z, &b_gemSegment_posGlb_z);
   fChain->SetBranchAddress("gemSegment_time", &gemSegment_time, &b_gemSegment_time);
   fChain->SetBranchAddress("gemSegment_time_err", &gemSegment_time_err, &b_gemSegment_time_err);
   fChain->SetBranchAddress("gemSegment_chi2", &gemSegment_chi2, &b_gemSegment_chi2);
   fChain->SetBranchAddress("gemSegment_posGlb_phi", &gemSegment_posGlb_phi, &b_gemSegment_posGlb_phi);
   fChain->SetBranchAddress("gemSegment_posGlb_eta", &gemSegment_posGlb_eta, &b_gemSegment_posGlb_eta);
   fChain->SetBranchAddress("gemSegment_dirGlb_phi", &gemSegment_dirGlb_phi, &b_gemSegment_dirGlb_phi);
   fChain->SetBranchAddress("gemSegment_dirGlb_eta", &gemSegment_dirGlb_eta, &b_gemSegment_dirGlb_eta);
   fChain->SetBranchAddress("mu_nMuons", &mu_nMuons, &b_mu_nMuons);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_isGlobal", &mu_isGlobal, &b_mu_isGlobal);
   fChain->SetBranchAddress("mu_isStandalone", &mu_isStandalone, &b_mu_isStandalone);
   fChain->SetBranchAddress("mu_isTracker", &mu_isTracker, &b_mu_isTracker);
   fChain->SetBranchAddress("mu_isGEM", &mu_isGEM, &b_mu_isGEM);
   fChain->SetBranchAddress("mu_isCSC", &mu_isCSC, &b_mu_isCSC);
   fChain->SetBranchAddress("mu_isME11", &mu_isME11, &b_mu_isME11);
   fChain->SetBranchAddress("mu_isLoose", &mu_isLoose, &b_mu_isLoose);
   fChain->SetBranchAddress("mu_isMedium", &mu_isMedium, &b_mu_isMedium);
   fChain->SetBranchAddress("mu_isTight", &mu_isTight, &b_mu_isTight);
   fChain->SetBranchAddress("mu_propagated_isME11", &mu_propagated_isME11, &b_mu_propagated_isME11);
   fChain->SetBranchAddress("mu_propagated_TrackNormChi2", &mu_propagated_TrackNormChi2, &b_mu_propagated_TrackNormChi2);
   fChain->SetBranchAddress("mu_propagated_numberOfValidPixelHits", &mu_propagated_numberOfValidPixelHits, &b_mu_propagated_numberOfValidPixelHits);
   fChain->SetBranchAddress("mu_propagated_innerTracker_ValidFraction", &mu_propagated_innerTracker_ValidFraction, &b_mu_propagated_innerTracker_ValidFraction);
   fChain->SetBranchAddress("mu_propagated_numberOfValidTrackerHits", &mu_propagated_numberOfValidTrackerHits, &b_mu_propagated_numberOfValidTrackerHits);
   fChain->SetBranchAddress("mu_path_length", &mu_path_length, &b_mu_path_length);
   fChain->SetBranchAddress("mu_isinsideout", &mu_isinsideout, &b_mu_isinsideout);
   fChain->SetBranchAddress("mu_isincoming", &mu_isincoming, &b_mu_isincoming);
   fChain->SetBranchAddress("mu_propagated_region", &mu_propagated_region, &b_mu_propagated_region);
   fChain->SetBranchAddress("mu_propagated_layer", &mu_propagated_layer, &b_mu_propagated_layer);
   fChain->SetBranchAddress("mu_propagated_chamber", &mu_propagated_chamber, &b_mu_propagated_chamber);
   fChain->SetBranchAddress("mu_propagated_etaP", &mu_propagated_etaP, &b_mu_propagated_etaP);
   fChain->SetBranchAddress("mu_propagated_pt", &mu_propagated_pt, &b_mu_propagated_pt);
   fChain->SetBranchAddress("mu_propagated_phi", &mu_propagated_phi, &b_mu_propagated_phi);
   fChain->SetBranchAddress("mu_propagated_eta", &mu_propagated_eta, &b_mu_propagated_eta);
   fChain->SetBranchAddress("mu_propagated_charge", &mu_propagated_charge, &b_mu_propagated_charge);
   fChain->SetBranchAddress("mu_propagatedLoc_x", &mu_propagatedLoc_x, &b_mu_propagatedLoc_x);
   fChain->SetBranchAddress("mu_propagatedLoc_y", &mu_propagatedLoc_y, &b_mu_propagatedLoc_y);
   fChain->SetBranchAddress("mu_propagatedLoc_z", &mu_propagatedLoc_z, &b_mu_propagatedLoc_z);
   fChain->SetBranchAddress("mu_propagatedLoc_r", &mu_propagatedLoc_r, &b_mu_propagatedLoc_r);
   fChain->SetBranchAddress("mu_propagatedLoc_phi", &mu_propagatedLoc_phi, &b_mu_propagatedLoc_phi);
   fChain->SetBranchAddress("mu_propagatedLoc_errX", &mu_propagatedLoc_errX, &b_mu_propagatedLoc_errX);
   fChain->SetBranchAddress("mu_propagatedLoc_errY", &mu_propagatedLoc_errY, &b_mu_propagatedLoc_errY);
   fChain->SetBranchAddress("mu_propagatedLoc_dirX", &mu_propagatedLoc_dirX, &b_mu_propagatedLoc_dirX);
   fChain->SetBranchAddress("mu_propagatedGlb_x", &mu_propagatedGlb_x, &b_mu_propagatedGlb_x);
   fChain->SetBranchAddress("mu_propagatedGlb_y", &mu_propagatedGlb_y, &b_mu_propagatedGlb_y);
   fChain->SetBranchAddress("mu_propagatedGlb_z", &mu_propagatedGlb_z, &b_mu_propagatedGlb_z);
   fChain->SetBranchAddress("mu_propagatedGlb_r", &mu_propagatedGlb_r, &b_mu_propagatedGlb_r);
   fChain->SetBranchAddress("mu_propagatedGlb_phi", &mu_propagatedGlb_phi, &b_mu_propagatedGlb_phi);
   fChain->SetBranchAddress("mu_propagatedGlb_errX", &mu_propagatedGlb_errX, &b_mu_propagatedGlb_errX);
   fChain->SetBranchAddress("mu_propagatedGlb_errY", &mu_propagatedGlb_errY, &b_mu_propagatedGlb_errY);
   fChain->SetBranchAddress("mu_propagatedGlb_errR", &mu_propagatedGlb_errR, &b_mu_propagatedGlb_errR);
   fChain->SetBranchAddress("mu_propagatedGlb_errPhi", &mu_propagatedGlb_errPhi, &b_mu_propagatedGlb_errPhi);
   fChain->SetBranchAddress("mu_propagated_EtaPartition_centerX", &mu_propagated_EtaPartition_centerX, &b_mu_propagated_EtaPartition_centerX);
   fChain->SetBranchAddress("mu_propagated_EtaPartition_centerY", &mu_propagated_EtaPartition_centerY, &b_mu_propagated_EtaPartition_centerY);
   fChain->SetBranchAddress("mu_propagated_EtaPartition_rMax", &mu_propagated_EtaPartition_rMax, &b_mu_propagated_EtaPartition_rMax);
   fChain->SetBranchAddress("mu_propagated_EtaPartition_rMin", &mu_propagated_EtaPartition_rMin, &b_mu_propagated_EtaPartition_rMin);
   fChain->SetBranchAddress("mu_propagated_EtaPartition_phiMax", &mu_propagated_EtaPartition_phiMax, &b_mu_propagated_EtaPartition_phiMax);
   fChain->SetBranchAddress("mu_propagated_EtaPartition_phiMin", &mu_propagated_EtaPartition_phiMin, &b_mu_propagated_EtaPartition_phiMin);
   fChain->SetBranchAddress("mu_propagated_Innermost_x", &mu_propagated_Innermost_x, &b_mu_propagated_Innermost_x);
   fChain->SetBranchAddress("mu_propagated_Innermost_y", &mu_propagated_Innermost_y, &b_mu_propagated_Innermost_y);
   fChain->SetBranchAddress("mu_propagated_Innermost_z", &mu_propagated_Innermost_z, &b_mu_propagated_Innermost_z);
   fChain->SetBranchAddress("mu_propagated_Outermost_x", &mu_propagated_Outermost_x, &b_mu_propagated_Outermost_x);
   fChain->SetBranchAddress("mu_propagated_Outermost_y", &mu_propagated_Outermost_y, &b_mu_propagated_Outermost_y);
   fChain->SetBranchAddress("mu_propagated_Outermost_z", &mu_propagated_Outermost_z, &b_mu_propagated_Outermost_z);
   Notify();
}

Bool_t Analysis_Noisy_Strip::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis_Noisy_Strip::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analysis_Noisy_Strip::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_Noisy_Strip_cxx
