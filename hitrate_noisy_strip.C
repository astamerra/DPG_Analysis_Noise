#include "TH1F.h"
#include <cmath>
#include <string>


#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLine.h>

void hitrate_noisy_strip()
{
    TString inputfile = "Analysis_miniTree.root"; 
    TChain *t1 = new TChain("outputLumiTree");
    t1->Add(inputfile);

    
    
   TFile saving("plots_etaP.root", "recreate");
    
    TH1F *hrate_etaP[8];
    TH1F *hrate_etaP_without_noise[8];
 


    
        Float_t media_hit_rate[8] = {0.};
       Float_t rms_hit_rate[8] = {0.};
        
        

	
    

					for(int p=0; p<8; p++){
						
         						
        						TString plot2 = "nHits_perEtaP_without_noise[][][]["+std::to_string(p)+"]/nEvents";
						//TString plot2 = "nHits_perEtaP_without_noise["+std::to_string(e)+"]["+std::to_string(l)+"]["+std::to_string(c)+"]["+std::to_string(p)+"]/nEvents";
		    					TString plot2_name = "hrate_etaP_without_noise["+std::to_string(p)+"]";
                    		
                    					t1->Draw(plot2+">>"+plot2_name);
                
                     					hrate_etaP_without_noise[p] = (TH1F *)gDirectory->Get(plot2_name);
                     					hrate_etaP_without_noise[p]->SetTitle(Form("run342955 iEta %d", p+1));
                     					hrate_etaP_without_noise[p]->GetXaxis()->SetTitle("nRecHits / nEvent per lumisection");
                     					media_hit_rate[p] = hrate_etaP_without_noise[p]->GetMean();  // take mean and stddev directly from the plots
                     					rms_hit_rate[p] = hrate_etaP_without_noise[p]->GetRMS();
                     			}
        
 
                
    for(int v=0; v<8; v++){
   	hrate_etaP[v]->Write();
        hrate_etaP_without_noise[v]->Write();            		
     	cout << "Histo " << v+1 << " correctly written" << endl;
    }
    
    
saving.Close();
                                        
    	Int_t nHits_perLayer[2][2] = {0};
       Int_t nHits_perEtaP[2][2][36][8] = {0};
       Int_t nHits_perEtaP_without_noise[2][2][36][8] = {0};
       Int_t nHits_perStrip[2][2][36][8][384] = {0};
       Int_t maximum_nhits = 0;
       
       
       Float_t temp_nHits_perEtaP_without_noise[2][2][36][8] = {0.};
       Float_t temp_nHits_perStrip[2][2][36][8][384] = {0.};
       Int_t nEvents = 0;
       Float_t temp_nEvents = 0.;
       

       Float_t hitrate_perEtaP_without_noise[2][2][36][8] = {0.};
       Float_t hitrate_perStrip[2][2][36][8][384] = {0.};

       
       Float_t sum[8] = {0.};
       Float_t sum2[8] = {0.};
       
       Float_t discrimine[8] = {0.};
       
       Float_t den = 144.0*384.0;
       
       t1->SetBranchAddress("nEvents",&nEvents);
       t1->SetBranchAddress("nHits_perLayer",&nHits_perLayer);
       t1->SetBranchAddress("nHits_perEtaP",&nHits_perEtaP);
       t1->SetBranchAddress("nHits_perStrip",&nHits_perStrip);
       t1->SetBranchAddress("nHits_perEtaP_without_noise",&nHits_perEtaP_without_noise);
       Int_t nentries = (Int_t)t1->GetEntries();
  for (Int_t i=0; i<nentries; i++) {
      t1->GetEntry(i);
       }
       
//Calcolo medie per etaP
       
                        for(int e=0; e<2; e++){
                        for(int l=0; l<2; l++){
                                for(int c=0; c<36; c++){
                                	for(int p=0; p<8; p++){
                                	
                                			 temp_nHits_perEtaP_without_noise[e][l][c][p]=static_cast<Float_t>(nHits_perEtaP_without_noise[e][l][c][p]);
                                                 	 temp_nEvents=static_cast<Float_t>(nEvents);
                                	
                                			 hitrate_perEtaP_without_noise[e][l][c][p] = temp_nHits_perEtaP_without_noise[e][l][c][p] / temp_nEvents;
                                       	         sum[p] += hitrate_perEtaP_without_noise[e][l][c][p];
                                       	     //    sum2[p] += pow(hitrate_perEtaP_without_noise[e][l][c][p], 2);
                                       		// cout <<  "Eta_P " << p+1 << " Somma " << sum[p] << std::endl;
                                	}
                        	}
                	}
        	}
        	
        	
        for(int s=0;s<8;s++){
              //  media_hit_rate[s] = (sum[s]/den);
               // rms_hit_rate[s] = sqrt ( sum2[s]/den );
   
              //  std::cout <<  "Eta_P " << s+1 << " Media " << media_hit_rate[s] << " RMS " << rms_hit_rate[s] << std::endl;
        }

        for(int e=0; e<2; e++){
                        for(int l=0; l<2; l++){
                                for(int c=0; c<36; c++){
                                	for(int p=0; p<8; p++){
                                	//	sum2[p] += pow(hitrate_perEtaP_without_noise[e][l][c][p]-media_hit_rate[p], 2);
                                	}
                                }
                         }
        }
        
        
        for(int s=0;s<8;s++){
               
             //   rms_hit_rate[s] = sqrt ( sum2[s]/den );
   		discrimine[s] = ((media_hit_rate[s] + 5*rms_hit_rate[s]));
                std::cout <<  "Eta_P " << s+1 << " Media " << media_hit_rate[s] << " StdDev " << rms_hit_rate[s] << std::endl;
                std::cout << "discrimine " << discrimine[s] << endl;
        }
    

    
    ofstream SaveFile("Noisy_strip_Run342955.txt");
    ofstream SaveFile2("Normal_strip_Run342955.txt");
    SaveFile << "Endcap " << "Layer " << "Chamber " << "iEta " << "Strip" << endl;
    SaveFile2 << "Endcap " << "Layer " << "Chamber " << "iEta " << "Strip" << " hit rate " << endl;
 for(int e=0;e<2;e++){
                for(int l=0;l<2;l++){
                        for(int c=0;c<36;c++){
                                for(int p=0;p<8;p++){
                                        for(int s=0;s<384;s++){
                                        

                                        	  temp_nHits_perStrip[e][l][c][p][s]=static_cast<Float_t>(nHits_perStrip[e][l][c][p][s]);
                                                 	 temp_nEvents=static_cast<Float_t>(nEvents);
                                	
                                			 hitrate_perStrip[e][l][c][p][s] = temp_nHits_perStrip[e][l][c][p][s] / temp_nEvents;
                                               // for(int i=0;i<9;i++){
                                                if(hitrate_perStrip[e][l][c][p][s] > discrimine[p]) { //noisy strips
                                                      //  n_noisy_strip[i]++;
                                                      SaveFile <<  e << " " <<  l+1 <<  " " << c+1 << " " <<  p+1 << " " << s << "	" << hitrate_perStrip[e][l][c][p][s] <<endl;
                                                      
                                                       }
                                                       if(hitrate_perStrip[e][l][c][p][s] < discrimine[p]) { 
                                                      SaveFile2 <<  e << " " <<  l+1 <<  " " << c+1 << " " <<  p+1 << " " << s << "	" << hitrate_perStrip[e][l][c][p][s] << endl;
					 		}
				          }
				  }
			  }
		 }
	}
	
	     
     

    

    
    
    
    }
