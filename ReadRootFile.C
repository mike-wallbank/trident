// Load dynamic library containing the definitions of the event classes
R__LOAD_LIBRARY(dunend/build/src/libEvent)

void ReadRootFile()
{
  // Load file and tree
  TFile file("ArSMmumuDec11.root");
  TTree* tree = (TTree*) file.Get("Event");

  // Sync tree with empty event
  Event* evt = 0;
  tree->SetBranchAddress("Event", &evt);

  // // Test that everything worked as it is supposed to
  // tree->GetEntry(100);
  // std::cout << "Event no. " << evt->GetMCID()
  //           << " contains " << evt->MCParticleContainer.size()
  //           << " MCParticles." << std::endl;

  //for (size_t event = 0; event < tree->GetEntriesFast(); ++event) {
  for (size_t event = 0; event < 10; ++event) {

    tree->GetEntry(event);

    // std::cout << "Event number " << evt->GetMCID() << " contains " << evt->MCParticleContainer.size()
    // 	      << " particles and " << evt->MCTrackContainer.size() << " tracks" << std::endl;
    // for (std::vector<MCTrack>::iterator track = evt->MCTrackContainer.begin(); track != evt->MCTrackContainer.end(); ++track) {
    //   std::cout << "  Track " << std::distance(evt->MCTrackContainer.begin(),track) << " has " << track->MCHitContainer.size() << " hits" << std::endl;
    //   for (std::vector<MCHit>::iterator hit = track->MCHitContainer.begin(); hit != track->MCHitContainer.end(); ++hit) {
    // 	TLorentzVector xyzt = hit->GetPositionAndTime();
    // 	std::cout << "    Hit " << std::distance(track->MCHitContainer.begin(),hit) << " has "
    // 		  << "(x,y,z) = (" << xyzt.X() << ", " << xyzt.Y() << ", " << xyzt.Z() << ")" << std::endl;
    //   }
    // }

    for (std::vector<MCTrack>::iterator track = evt->MCTrackContainer.begin(); track != evt->MCTrackContainer.end(); ++track) {
      std::cout << "  Track " << std::distance(evt->MCTrackContainer.begin(),track) << " has " << track->MCHitContainer.size() << " hits" << std::endl;
      for (std::vector<MCHit>::iterator hit = track->MCHitContainer.begin(); hit != track->MCHitContainer.end(); ++hit) {
    	TLorentzVector xyzt = hit->GetPositionAndTime();
    	std::cout << "    Hit " << std::distance(track->MCHitContainer.begin(),hit) << " has "
    		  << "(x,y,z) = (" << xyzt.X() << ", " << xyzt.Y() << ", " << xyzt.Z() << ")" << std::endl;
      }
    }
    
  }
  
  file.Close();
}
