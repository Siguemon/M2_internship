//
// SimpleAna.cxx in
// /sps/ricochet/goyc/neutronSim/Code/SharedCodes/SimulationAnalysis/SimuAnalysis/
//
// Standard includes
#include <iostream>
#include <string>
#include <functional>
#include <limits>
#include <cmath>
#include "inttypes.h"

/// ROOT includes
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TRint.h"
#include "TStyle.h"
#include "TGraph.h"

/// Ricochet includes
#include "SimuOutput/SimuRunInfo.hh"
#include "SimuOutput/SimuTreeReader.hh"
#include "SimuOutput/InfoHeaderReader.hh"
#include "SimuOutput/PDGConverter.hh"
#include "SimuOutput/VolumeNames.hh"
#include "SimAnaCuts/MaxEnergyDetectorType.hh"
#include "SimAnaCuts/PartitionedEnergy.hh"

/// CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

/// Project includes
#include "ROOT_Tools.hh"
#include "CPP_Tools.hh"
#include "SimuHisto.hh"
#include "Counter.hh"
#include "FullSimuEvent.hh"
#include "Cuts.hh"
#include "BranchReader.hh"
#include "ParseNormOption.hh"
#include "UpdateSimulatedTime.hh"
#include "VolumesReader.hh"
#include "TypedWithThreshold.hh"
#include "SensitiveVetoCut.hh"
//#include "MaxEnergyDetectorType.hh"
// #include "PartitionedEnergy.hh"

void Usage()
{
    std::cout << "Usage: AnalyseSimu -i <SimuOuputFilename_1.root> ... <SimuOuputFilename_N.root> -o <ouput_filename> [-b]\n"
         << "Options:\n"
         << "   -b                            Batch mode, no plot display, exit root-appli \n"
	  << "   -o <filename>                 Choose name of the ROOT output\n" ;
}

int main(int argc, char* argv[])
{
    if(CmdOptionExists(argc, argv, "-h"))
    {
        Usage();
        return EXIT_SUCCESS;
    }

    std::vector<std::string> inputFilenames;
    if(CmdOptionExists(argc, argv, "-i"))
    {
        inputFilenames = GetCmdOptionArray<std::string>(argc, argv, "-i");
    }
    else return EXIT_FAILURE;

    // defautl output rootfile name, overwritten by option-o
    std::string outputFilename = "AnalyseSimu.root";
    if(CmdOptionExists(argc, argv, "-o"))
        outputFilename = GetCmdOption<std::string>(argc, argv, "-o");
    bool batchMode = CmdOptionExists(argc, argv, "-b");



    TH1D::AddDirectory(false);


    std::cout << "List of volumeNames"  << std::endl;
    auto volumeNames = VolumesReader::GetVolumeNames(inputFilenames.front());
    for (auto it = volumeNames.begin() ; it != volumeNames.end() ; ++it) {
      std::cout << it->first << "   "  << it->second << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TFile* outfile = nullptr;
    if(!outputFilename.empty()) outfile = new TFile(outputFilename.c_str(), "RECREATE");


    TTree *tTree = new TTree("ESO","EventSummaryOutput");
    Int_t ev_num;
    tTree->Branch("EvNum",&ev_num,"ev_num/I");
    Float_t time_sim;
    tTree->Branch("IntTime",&time_sim);
    Int_t gen_nb;
    tTree->Branch("Ngen",&gen_nb,"gen_nb/I");
    std::vector<int> gen_pdgid;
    std::vector<float> gen_px;
    std::vector<float> gen_py;
    std::vector<float> gen_pz;
    std::vector<float> gen_e;
    std::vector<float> gen_ekin;
    std::vector<float> gen_x;
    std::vector<float> gen_y;
    std::vector<float> gen_z;
    tTree->Branch("GenPdg",&gen_pdgid);
    tTree->Branch("GenPx",&gen_px);
    tTree->Branch("GenPy",&gen_py);
    tTree->Branch("GenPz",&gen_pz);
    tTree->Branch("GenE",&gen_e);
    tTree->Branch("GenEkin",&gen_ekin);
    tTree->Branch("GenX",&gen_x);
    tTree->Branch("GenY",&gen_y);
    tTree->Branch("GenZ",&gen_z);
    Int_t bolo_nb;
    Int_t topscint_nb;
    Int_t sidescint_nb;
    Int_t cryoscint_nb;
    tTree->Branch("Nbolo",&bolo_nb,"bolo_nb/I");
    tTree->Branch("Ntop",&topscint_nb,"topscint_nb/I");
    tTree->Branch("Nside",&sidescint_nb,"sidescint_nb/I");
    tTree->Branch("Ncryo",&cryoscint_nb,"cryoscint_nb/I");
    // Bolo info
    std::vector<int> bolo_num;
    std::vector<float> bolo_e;
    std::vector<float> bolo_nuclear;
    std::vector<float> bolo_elec;
    std::vector<int> bolo_SID;
    tTree->Branch("BoloNum",&bolo_num);
    tTree->Branch("BoloE",&bolo_e);
    tTree->Branch("BoloNR",&bolo_nuclear);
    tTree->Branch("BoloER",&bolo_elec);
    tTree->Branch("SID",&bolo_SID);
    // Hit info
    std::vector<double> hit_x;
    std::vector<double> hit_y;
    std::vector<double> hit_z;
    std::vector<double> local_hit_r;
    std::vector<double> local_hit_z;
    std::vector<double> deposited;
    std::vector<int> hit_pdg;
    tTree->Branch("HitX",&hit_x);
    tTree->Branch("HitY",&hit_y);
    tTree->Branch("HitZ",&hit_z);
    tTree->Branch("LocalHitR",&local_hit_r);
    tTree->Branch("LocalHitZ",&local_hit_z);
    tTree->Branch("Edep",&deposited);
    tTree->Branch("HitPdg",&hit_pdg);
    // Top Scint info
    std::vector<int> topscint_layer;
    std::vector<int> topscint_num;
    std::vector<float> topscint_e;
    tTree->Branch("TopNum",&topscint_num);
    tTree->Branch("TopLay",&topscint_layer);
    tTree->Branch("TopE",&topscint_e);
    std::vector<int> sidescint_layer;
    std::vector<int> sidescint_num;
    std::vector<float> sidescint_e;
    tTree->Branch("SideNum",&sidescint_num);
    tTree->Branch("SideLay",&sidescint_layer);
    tTree->Branch("SideE",&sidescint_e);
    std::vector<int> cryoscint_num;
    std::vector<int> cryoscint_e;
    tTree->Branch("CryoNum",&cryoscint_num);
    tTree->Branch("CryoE",&cryoscint_e);

    TH1I *h_ngen = new TH1I("h_ngen","  ",50,0.,50.);
    TH1I *h_nbolo = new TH1I("h_nbolo","   ",28,0.,28.);
    TH1I *h_ntop = new TH1I("h_ntop","   ",13,0.,13.);
    TH1I *h_nside = new TH1I("h_nside","   ",23,0.,23.);
    TH1I *h_ncryo = new TH1I("h_ncryo","   ",3,0.,3.);

    double totalSimulatedTime = 0.; /// in second
    unsigned long totalNumberOfEvents = 0;
    //
    // Compute cumulated time  (cosmic muon simulation)
    //
    for(const std::string& filename : inputFilenames) {
      std::cout << " Compute time - File Name:  " << filename.c_str() << std::endl;
      TFile *file = new TFile(filename.c_str(),"READ") ;
      if (file->IsZombie()) continue;
      // Original method; not possible to use IsZombie()
      //      std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "READ"));
      file->cd("GeneratorInfo");
      TTree *genTree = nullptr;
      genTree = LoadRootObject<TTree>("ILLFileM");
      if (!genTree) {
	std::cout << " Not an ILL  cosmic generator : Simulation time not computed " << std::endl;
	continue;
      }
      double simTime;
      genTree->Scan("SimulatedTime");
      genTree->SetBranchAddress("SimulatedTime",&simTime);
      int nb = genTree->GetEntries(); // Case of concatenated files
      for (int i=0;i<nb;i++) {
	genTree->GetEntry(i);
	totalSimulatedTime += simTime;
      }
    }

    time_sim = (float) totalSimulatedTime;

    int nfile = 0;
    for(const std::string& filename : inputFilenames)
    {
      TFile *file = new TFile(filename.c_str(),"READ") ;
      if (file->IsZombie()) continue;
      std::cout << "Loop on Files:  " << filename.c_str() << "  "  <<  nfile << std::endl;
      nfile++;
      TTree* infoTree = LoadRootObject<TTree>("SimuRunInfo");

      if (infoTree) {
	if (infoTree->GetBranch("NbOfEvents")) {
	  unsigned long long nb_of_ev;
	  infoTree->Scan("NbOfEvents");
	  infoTree->SetBranchAddress("NbOfEvents",&nb_of_ev);
	  infoTree->GetEntry(0);
	  std::cout << "main ::  Number of events simulated " << nb_of_ev << std::endl;
	}
      }

      TTree* mainTree = LoadRootObject<TTree>("SimuTree");
      TTree* stepTree = LoadRootObject<TTree, 0>("StepTree");
      std::cout << " Trees ? " << mainTree << " "  << stepTree <<  std::endl;

      // FullSimuEvent eventOri;
      // SimuTreeDev
      FullSimuEventV2 event;
      auto simuEventPtr = &(event.Event);
      mainTree->SetBranchAddress("SimuEvent", &simuEventPtr);
      auto stepsPtr = &(event.StepEvents);
      if(stepTree) stepTree->SetBranchAddress("StepEvent", &stepsPtr);
      uint32_t nbEntries = mainTree->GetEntries();
      totalNumberOfEvents += nbEntries;

      std::cout << "main:: Entries " << nfile << "   "  << nbEntries << "  "  << totalNumberOfEvents << std::endl;

      //nbEntries = 100.;

      for(uint32_t iEntry = 0; iEntry < nbEntries; ++iEntry)
        {
	  ev_num = iEntry;
	  // std::cout << " Event Num " << iEntry << std::endl;
	  // Reset outputTree variables
	  gen_nb = 0;
	  gen_pdgid.clear();
	  gen_px.clear();
	  gen_py.clear();
	  gen_pz.clear();
	  gen_x.clear();
	  gen_y.clear();
	  gen_z.clear();
	  gen_e.clear();
	  gen_ekin.clear();
	  //
	  bolo_nb = 0;
	  bolo_num.clear();
	  bolo_e.clear();
	  bolo_nuclear.clear();
	  bolo_elec.clear();
    bolo_SID.clear();
	  //
    hit_x.clear();
    hit_y.clear();
    hit_z.clear();
    local_hit_r.clear();
    local_hit_z.clear();
    deposited.clear();
    hit_pdg.clear();
    //
	  topscint_nb = 0;
	  topscint_num.clear();
	  topscint_layer.clear();
	  topscint_e.clear();
	  //
	  sidescint_nb = 0;
	  sidescint_num.clear();
	  sidescint_layer.clear();
	  sidescint_e.clear();
	  //
	  cryoscint_nb = 0;
	  cryoscint_num.clear();
	  cryoscint_e.clear();
	  //
	  // Get mainTree event
	  //
	  mainTree->GetEntry(iEntry);
	  //
	  // if applicable get stepTree;
	  //
	  if (stepTree) stepTree->GetEntry(iEntry);
	  int n_steps = (event.StepEvents).size();

	  //
	  // Get Primary Event information
	  //
	  gen_nb = event.nPrim;
	  h_ngen->Fill(gen_nb);
	  for(int i=0;i<event.nPrim;i++)
	    {
	      int pdg_id = (event.PrimaryEvents)[i].PDG;
	      // Gen Position en mm ?
	      float x = (event.PrimaryEvents)[i].Position.X/CLHEP::m;
	      float y = (event.PrimaryEvents)[i].Position.Y/CLHEP::m;
	      float z = (event.PrimaryEvents)[i].Position.Z/CLHEP::m;
	      float px = (event.PrimaryEvents)[i].Px;
	      float py = (event.PrimaryEvents)[i].Py;
	      float pz = (event.PrimaryEvents)[i].Pz;
	      float egev = (event.PrimaryEvents)[i].E/CLHEP::GeV;
	      float ekin = (event.PrimaryEvents)[i].KineticEnergy();
	      //
	      gen_px.emplace_back(px);
	      gen_py.emplace_back(py);
	      gen_pz.emplace_back(pz);
	      gen_e.emplace_back(egev);  //GeV
	      gen_ekin.emplace_back(ekin);  // MeV
	      gen_x.emplace_back(x);
	      gen_y.emplace_back(y);
	      gen_z.emplace_back(z);
	      gen_pdgid.emplace_back(pdg_id);
	    } // end of loop on primaries

	  //
	  // Get  bolometers informations
	  //
	  bolo_nb = event.nDet;
	  h_nbolo->Fill(bolo_nb);
	  for(int i=0;i<event.nDet;i++) {
	    unsigned int detId = (event.DetectorEvents)[i].DetectorID;
	    std::string volName;
	    for (auto it = volumeNames.begin() ; it != volumeNames.end() ; ++it) {
	      if (detId == it->first) volName = it->second;
	    }

	    int detSID = (event.DetectorEvents)[i].DetectorSID;
	    int num = detSID/100;
	    num = num%100;
	    std::cout << " VolName "  <<  volName << "   "  <<  (event.DetectorEvents)[i].DetectorSID<< "  "  << num << std::endl;
	    // Energy per bolometer
	    float e =  (event.DetectorEvents)[i].Energy();  // MeV - NB Energy() est une methode
	    float energy = (event.DetectorEvents)[i].DetectorE;
	    float enucl = (event.DetectorEvents)[i].DetectorNR;
	    float eelec = (event.DetectorEvents)[i].DetectorER;
	    //
	    //  Read Hits if needed
	    //

	    int nhits = (event.DetectorEvents)[i].HitEvents.size();
	    for (const auto& hit : (event.DetectorEvents)[i].HitEvents) {
	      double X = hit.Position.X/CLHEP::m;
	      double Y = hit.Position.Y/CLHEP::m;
              double Z = hit.Position.Z/CLHEP::m;
        double lR = pow((hit.LocalPosition.X/CLHEP::m)*(hit.LocalPosition.X/CLHEP::m) + (hit.LocalPosition.Y/CLHEP::m)*(hit.LocalPosition.Y/CLHEP::m), .5); //positions locales dans le volume du crystal !
        double lZ = hit.LocalPosition.Z/CLHEP::m;
        int PDGH = hit.PDG;
        hit_x.emplace_back(X);
        hit_y.emplace_back(Y);
        hit_z.emplace_back(Z);
        local_hit_r.emplace_back(lR);
        local_hit_z.emplace_back(lZ);
        deposited.emplace_back(hit.Edep);
        hit_pdg.emplace_back(PDGH);
	    }
	    //
	    // Fill EventSumaryOutput
	    //

	    bolo_num.emplace_back(num);
	    bolo_e.emplace_back(energy);
	    bolo_nuclear.emplace_back(enucl);
	    bolo_elec.emplace_back(eelec);
      bolo_SID.emplace_back(detSID);
	  }

	  // Loop on TopScintillators
	  topscint_nb = event.nTop;
	  h_ntop->Fill(topscint_nb);
	  for(int i=0;i<topscint_nb;i++)
	    {
	      int detSID = (event.TopScintEvents)[i].DetectorSID;
	      int scintLayer = detSID%100;
	      int scintModule = ((detSID-scintLayer)/100)%100;
	      float energy =   (event.TopScintEvents)[i].DetectorE; 	      // MeV
	      topscint_num.emplace_back(scintModule);
	      topscint_layer.emplace_back(scintLayer);
	      topscint_e.emplace_back(energy);
	    }// End TopScint

	  // Loop on SideScintillators
	  sidescint_nb = event.nSide;
	  h_nside->Fill(sidescint_nb);
	  for(int i=0;i<sidescint_nb;i++)
	    {
	      int detSID = (event.SideScintEvents)[i].DetectorSID;
	      int scintLayer = detSID%100;
	      int scintModule = ((detSID-scintLayer)/100)%100;
	      float energy =  (event.SideScintEvents)[i].DetectorE;	      // MeV
	      sidescint_num.emplace_back(scintModule);
	      sidescint_layer.emplace_back(scintLayer);
	      sidescint_e.emplace_back(energy);
	    } // End of SideScint
	  //
	  // Loop on Cryo Scintillators
	  //
	  cryoscint_nb = event.nCryo;
	  h_ncryo->Fill(cryoscint_nb);
	  for(int i=0;i<cryoscint_nb;i++)
	    {
	      int detSID = (event.CryoScintEvents)[i].DetectorSID;
	      int cryoModule = detSID%100;
	      float energy =  (event.CryoScintEvents)[i].DetectorE; 	      // MeV
	      cryoscint_num.emplace_back(cryoModule);
	      cryoscint_e.emplace_back(energy);
	    } // End of CryoScint

	  tTree->Fill();
	}

     }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

     auto Process = [&] ()
       {
	 outfile->cd();
	 tTree->Write();
	 h_ngen->Write();
	 h_nbolo->Write();
	 h_ntop->Write();
	 h_nside->Write();
	 h_ncryo->Write();
        if(batchMode) delete outfile;
       };

    if(batchMode) Process();
    else RunRootApp(Process);

    return EXIT_SUCCESS;
    }
