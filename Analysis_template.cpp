#include "Analysis_Noisy_Strip.C"
#include <TROOT.h>
#include <stdio.h>
#include <iostream>
#include <string.h>

// ###### INSTRUCTIONS (input options)
// // (after having compiled the code [.L Analysis.cpp])
//
 using namespace std;

 int Analysis_template(){
    TString fileout = "Analysis_miniTree_342955.root";

    TChain* chain = new TChain("muNtupleProducer/MuDPGTree");
           //AddFile
    chain->Add("/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/2021/GEMCommonNtuples/CRUZET/Run_342955/MuDPGNtuple_2021_CRUZET_342955.root");
/*  for(int i=1;i<10;i++){
    //for(int i=1;i<2;i++){
          chain->Add(Form("/lustre/cms/store/user/gmilella/ExpressCosmics/CRAB3_gem_dpg_ntuple_mwgr3_2021_run341343_express/210507_215122/0000/MuDPGNTuple_MWGR3_2021_EXP_run341343_%d.root", i));        
    }       
 *///OutFile
     chain->Merge(fileout);

     Analysis_Noisy_Strip class_noise(chain, fileout);
     class_noise.Loop();

      return 0;
    }
       
