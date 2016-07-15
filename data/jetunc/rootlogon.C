{
   std::cout << "Loading rootlogon.C" << std::endl;
   gROOT->ProcessLine(".L tdrstyle.C");
   setTDRStyle();
   gStyle->SetOptStat(1111111);
   gStyle->SetHistLineWidth(2);
   //gStyle->SetHistFillStyle(1);
   gStyle->SetMarkerStyle(23);
   gStyle->SetMarkerSize(0.3);
   gStyle->SetErrorX(0.5);
   gStyle->SetPalette(1);

   //gSystem->Load("libFWCoreFWLite");
   //AutoLibraryLoader::enable();
   gSystem->Load("libCondFormatsJetMETObjects.so");

   //gSystem->Load("/storage/lhuertas/uerj-storage2/antoniov/Workspace/Analysis/CMSTotem/CMSTotem/TOTEMdataFormat/lib/libTOTEMdataFormat.so");
   //gSystem->Load("/storage/lhuertas/uerj-storage2/antoniov/Workspace/Analysis/CMSTotem/CMSTotem/CMSdataFormat/lib/libCMSdataFormat.so");
   gSystem->Load("/storage/lhuertas/uerj-1/CMSTOTEM/data/CMSSW_5_3_4/src/CMSTotem/TOTEMdataFormat/lib/libTOTEMdataFormat.so");
   gSystem->Load("/storage/lhuertas/uerj-1/CMSTOTEM/data/CMSSW_5_3_4/src/CMSTotem/CMSdataFormat/lib/libCMSdataFormat.so");

   gSystem->AddIncludePath(" -I/storage/lhuertas/uerj-storage2/antoniov/Workspace/Analysis/CMSTotem/CMSTotem/CMSdataFormat/src ");
   gSystem->AddIncludePath(" -I/storage/lhuertas/uerj-storage2/antoniov/Workspace/Analysis/CMSTotem/CMSTotem/TOTEMdataFormat/src ");


}
