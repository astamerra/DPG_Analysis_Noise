#define Analysis_Noisy_Strip_cxx
#include "Analysis_Noisy_Strip.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

//initialise output minitree
void outTreeInit(TTree *tree_, Int_t &run, Int_t &lumi, Int_t &evt, Int_t cout[2][2][36]){
    tree_->Branch("run", &run);
    tree_->Branch("evt", &evt);
    tree_->Branch("lumi", &lumi);
    tree_->Branch("count", cout, "count[2][2][36]/I");
}

//initialise output minitree
void outTreeLumiInit(TTree *tree_, Int_t &nEvents,
		      Int_t &lumi, 
		      Int_t nHits_perEndcap[2],
		      Int_t nHits_perLayer[2][2],
		      Int_t nHits_perChamber[2][2][36],
		      Int_t nHits_perEtaP[2][2][36][8],
		      Int_t nHits_perEtaP_without_noise[2][2][36][8],
		      Int_t nHits_perStrip[2][2][36][8][384]){
    tree_->Branch("nEvents", &nEvents, "nEvents/I");
    tree_->Branch("lumi", &lumi);
    tree_->Branch("nHits_perEndcap", nHits_perEndcap, "nHits_perEndcap[2]/I");
    tree_->Branch("nHits_perLayer", nHits_perLayer, "nHits_perLayer[2][2]/I");
    tree_->Branch("nHits_perChamber", nHits_perChamber, "nHits_perChamber[2][2][36]/I");
    tree_->Branch("nHits_perEtaP", nHits_perEtaP, "nHits_perEtaP[2][2][36][8]/I");
    tree_->Branch("nHits_perEtaP_without_noise", nHits_perEtaP_without_noise, "nHits_perEtaP_without_noise[2][2][36][8]/I");
    tree_->Branch("nHits_perStrip", nHits_perStrip, "nHits_perStrip[2][2][36][8][384]/I");
}



void Analysis_Noisy_Strip::Loop()
{
//   In a ROOT session, you can do:
//      root> .L analysisClass.C
//      root> analysisClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   // Creation of output file & final tree
   TString root_fileName = fileName;
   TFile *fout = new TFile(root_fileName, "RECREATE");
   fout->cd();
   TTree *outTree = new TTree("outputTree","outputTree");
   TTree *outLumiTree = new TTree("outputLumiTree","outputLumiTree");
   TTree *outNoisyLumiTree = new TTree("outputNoisyLumiTree","outputNoisyLumiTree");
   //initialise variables
   Int_t run = 0, lumi = 0, evt = 0;
   int counts[2][2][36] = {0};
   int region_index[2] = {-1, +1};
   int layer_index[2] = {1, 2};

   
   vector<Int_t> lumi_list;
   Int_t nevents = 0, nevents1 = 0;
   Int_t nHits_perEndcap[2] = {0};
   Int_t nHits_perLayer[2][2] = {0};
   Int_t nHits_perChamber[2][2][36] = {0};
   Int_t nHits_perEtaP[2][2][36][8] = {0};
   Int_t nHits_perEtaP_without_noise[2][2][36][8] = {0};
   Int_t nHits_perStrip[2][2][36][8][384] = {0};
   

   
      //initialise output minitree
   outTreeInit(outTree, run, lumi, evt, counts);
  
   outTreeLumiInit(outLumiTree, nevents, lumi, nHits_perEndcap, nHits_perLayer, nHits_perChamber, nHits_perEtaP, nHits_perEtaP_without_noise, nHits_perStrip);
   
      if (fChain == 0) return;
   bool isVerbose = false;
   Long64_t nentries = fChain->GetEntriesFast();
  //std::cout << "a qui ci arrivo" << std::endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(isVerbose) cout<<"\n==================\nrun "<<event_runNumber<<" lumi "<<event_lumiBlock<<" evt "<<event_eventNumber<<endl;

      if(isVerbose) cout<<"n. GEM RecHit = "<<gemRecHit_region->size()<<endl;
      //loop on GEM rechits
      for(int j=0; j<gemRecHit_region->size(); j++){
          if(isVerbose) cout<<j<<" RecHit region|chamber|layer= "<<gemRecHit_region->at(j)<<" | "<<gemRecHit_chamber->at(j)<<" | "<<gemRecHit_layer->at(j)<<endl;
          //loop on endcaps
          for(int e=0; e<2; e++){
              //loop on layers
              for(int l=0; l<2 && gemRecHit_region->at(j) == region_index[e]; l++){
                  //loop on chambers
                  for(int c=0; c<36 && gemRecHit_layer->at(j) == layer_index[l]; c++){
                      if(gemRecHit_chamber->at(j) == (c+1)) counts[e][l][c]++;
                  }
              }
          }
      }
      //for each event, fill output Tree
      run = event_runNumber;
      evt = event_eventNumber;
      lumi = event_lumiBlock;
      outTree->Fill();
      //reset counts to zero
      memset(counts, 0, sizeof(counts));

      //per lumi analysis
      lumi_list.push_back(event_lumiBlock);
   } //loop on events
   
        //new loop on events for per lumi analysis
   nentries = fChain->GetEntriesFast();
   nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
            if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(isVerbose) cout<<"\n==================\nrun "<<event_runNumber<<" lumi "<<event_lumiBlock<<" evt "<<event_eventNumber<<endl;
      if(isVerbose) cout<<"jentry "<<jentry<<"event_lumiBlock "<<event_lumiBlock<<" lumi_list.at(jentry) "<<lumi_list.at(jentry)<<endl;
      if((jentry < nentries-1) && event_lumiBlock == lumi_list.at(jentry+1)){
          nevents++;
          if(isVerbose) cout<<"n. GEM RecHit = "<<gemRecHit_region->size()<<endl;
          //loop on GEM rechits
          for(int j=0; j<gemRecHit_region->size(); j++){
              if(isVerbose) cout<<j<<" RecHit region|chamber|layer= "<<gemRecHit_region->at(j)<<" | "<<gemRecHit_chamber->at(j)<<" | "<<gemRecHit_layer->at(j)<<endl;
              //loop on endcaps
              for(int e=0; e<2; e++){
             	 if(gemRecHit_region->at(j) == region_index[e]) nHits_perEndcap[e]++;
                  //loop on layers
                  for(int l=0; l<2 && gemRecHit_region->at(j) == region_index[e]; l++){
                  	if(gemRecHit_layer->at(j) == layer_index[l]) nHits_perLayer[e][l]++;
                      //loop on chambers
                      for(int c=0; c<36 && gemRecHit_layer->at(j) == layer_index[l]; c++){
                      	   if(gemRecHit_chamber->at(j) == c+1) nHits_perChamber[e][l][c]++;
                          //loop on etaP
                          for(int p=0; p<8 && gemRecHit_chamber->at(j) == (c+1); p++){
                          	if(gemRecHit_etaPartition->at(j) == p+1){ nHits_perEtaP[e][l][c][p]++;
                          		 if(not((e==0 && l==0 && c==19 && p==0) ||
                                                 (e==0 && l==0 && c==26 && p==1) ||
                                                 (e==0 && l==0 && c==26 && p==3) ||
                                                 (e==0 && l==0 && c==29 && p==6) ||
                                                 (e==0 && l==1 && c==2 && p==3) ||
                                                 (e==0 && l==1 && c==7 && p==2) ||
                                                 (e==0 && l==1 && c==19 && p==1) ||
                                                 (e==0 && l==1 && c==19 && p==2) ||
                                                 (e==0 && l==1 && c==24) ||
                                                 (e==1 && l==0 && c==27 && p==3) ||
                                                 (e==1 && l==1 && c==6 && p==3) ||
                                                 (e==1 && l==1 && c==25 && p==1) ||
                                                 (e==1 && l==1 && c==20 && p==6) ||
                                                 (e==1 && l==1 && c==27 && p==0))){
                                                 nHits_perEtaP_without_noise[e][l][c][p]++;
                                                 }
                                        }
                          	//loop on strips
                          	for(int s=0; s<384 && gemRecHit_etaPartition->at(j) == (p+1); s++){
                          		//loop on cluster size
                          		for(int cl=0; cl<gemRecHit_cluster_size->at(j) && gemRecHit_firstClusterStrip->at(j) == s; cl++){
                              		 nHits_perStrip[e][l][c][p][s+cl]++;
                                    	}
                              }
                          }
                      }
                  }
              }
          }      
       }else{ //if(jentry>0 && event_lumiBlock != lumi_list.at(jentry-1)){
          if(isVerbose){
              cout<<"event_lumiBlock "<<event_lumiBlock<<" nEvents="<<nevents;
              //cout<<"hits in chamber [1][1][20][21] "<<nhits[1][1][20][21]<<endl; 
              }
          lumi = event_lumiBlock;

	
	outLumiTree->Fill();
               
          memset(nHits_perEndcap, 0, sizeof(nHits_perEndcap));
          memset(nHits_perLayer, 0, sizeof(nHits_perLayer));
          memset(nHits_perChamber, 0, sizeof(nHits_perChamber));
          memset(nHits_perEtaP_without_noise, 0, sizeof(nHits_perEtaP_without_noise));
          memset(nHits_perStrip, 0, sizeof(nHits_perStrip));
          memset(nHits_perEtaP, 0, sizeof(nHits_perEtaP));

          nevents=0;  
   }

   } //loop on events 
   

   
      fout->cd();
   //Write and close the file
   fout->Write();
   fout->Close();

	
  }

     

   




