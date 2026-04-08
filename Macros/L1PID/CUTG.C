{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.32.06
   
   TCutG *cutg = new TCutG("CUTG",6);
   cutg->SetVarX("L1 PID cut in TLvsQvoxels");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetLineWidth(2);
   cutg->SetMarkerSize(1.2);
   cutg->SetPoint(0,11.1632,746.689);
   cutg->SetPoint(1,65.5506,424.419);
   cutg->SetPoint(2,75.0995,678.843);
   cutg->SetPoint(3,74.2692,967.19);
   cutg->SetPoint(4,6.59629,1408.19);
   cutg->SetPoint(5,11.1632,746.689);
   cutg->Draw("./Cuts/p_TLvsQvoxels_11Li.root");
}
