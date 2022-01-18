{ 
  TFile f("root/genstudy_pt.root");
  TTree* t = (TTree*)f.Get("tree");
  std::cout << "entries: " << t->GetEntries() << std::endl;
  nbins = 13;
  TH2F denom("denom","denom",nbins, 0. ,nbins*1., nbins, 0. ,nbins*1. );
  int n = tree->Draw("gen_e1_pt:gen_e2_pt>>denom","","goff");
  std::cout << "processed: " << n << std::endl;
  TFile fw("root/denom.root","RECREATE");
  denom.Write();
  fw.Close();
}
