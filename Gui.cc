#include "Gui.hh"
#include "Simulator.hh"

Gui::Gui(const TGWindow *p, UInt_t w, UInt_t h, TApplication* APP)
  : TGMainFrame(p, w, h), sim(nullptr), zoom(1), app(APP){

  Canvas = new TRootEmbeddedCanvas("EventCanvas", this, hsize, vsize);
  AddFrame(Canvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  
  TGHorizontalFrame *controls = new TGHorizontalFrame(this);
  
  NumEntry = new TGNumberEntry(controls, 1, 6, -1, TGNumberFormat::kNESInteger);
  controls->AddFrame(NumEntry, new TGLayoutHints(kLHintsCenterY, 5, 5, 5, 5));
  particlesLabel = new TGLabel(controls, "Photons");
  controls->AddFrame(particlesLabel, new TGLayoutHints(kLHintsCenterY, 2, 10, 5, 5));
  
  TGTextButton *genBtn = new TGTextButton(controls, "&Generate");
  genBtn->Connect("Clicked()", "Gui", this, "DoRun()");
  controls->AddFrame(genBtn, new TGLayoutHints(kLHintsCenterY, 5, 15, 5, 5));
  
  TGLabel *fileNameLabel = new TGLabel(controls, "File name");
  controls->AddFrame(fileNameLabel, new TGLayoutHints(kLHintsCenterY, 5, 2, 5, 5));
  
  FileNameEntry = new TGTextEntry(controls, new TGTextBuffer(20));
  FileNameEntry->SetText("event_display.png");
  controls->AddFrame(FileNameEntry, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 2, 5, 5, 5));
  
  TGTextButton *saveBtn = new TGTextButton(controls, "&Save");
  saveBtn->Connect("Clicked()", "Gui", this, "DoSave()");
  controls->AddFrame(saveBtn, new TGLayoutHints(kLHintsCenterY, 5, 15, 5, 5));
  
  TGTextButton *zoomBtn = new TGTextButton(controls, "&Zoom");
  zoomBtn->Connect("Clicked()", "Gui", this, "DoZoom()");
  controls->AddFrame(zoomBtn, new TGLayoutHints(kLHintsCenterY, 5, 5, 5, 5));

  TGTextButton *unzoomBtn = new TGTextButton(controls, "&Unzoom");
  unzoomBtn->Connect("Clicked()", "Gui", this, "DoUnzoom()");
  controls->AddFrame(unzoomBtn, new TGLayoutHints(kLHintsCenterY, 5, 15, 5, 5));

  TGTextButton *clearBtn = new TGTextButton(controls, "&Clear");
  clearBtn->Connect("Clicked()", "Gui", this, "DoClear()");
  controls->AddFrame(clearBtn, new TGLayoutHints(kLHintsCenterY, 5, 15, 5, 5));  

  onoffCheck = new TGCheckButton(controls, "Show Info.");
  onoffCheck->SetState(kButtonUp);
  onoffCheck->Connect("Toggled(Bool_t)", "Gui", this, "OnToggle(Bool_t)");
  controls->AddFrame(onoffCheck, new TGLayoutHints(kLHintsCenterY, 5, 5, 5, 5));
  
  TGTextButton *exitBtn = new TGTextButton(controls, "&Exit");
  exitBtn->Connect("Clicked()", "Gui", this, "DoExit()");
  controls->AddFrame(exitBtn, new TGLayoutHints(kLHintsCenterY, 5, 5, 5, 5));
  
  AddFrame(controls, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));

  SetWindowName("Optical simulation GUI");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();

  label[0] = new TLatex(0.01, 0.97, "Material color");
  label[1] = new TLatex(0.01, 0.97-0.03*1, "Normal medium");
  label[2] = new TLatex(0.01, 0.97-0.03*2, "Converter");
  label[3] = new TLatex(0.01, 0.97-0.03*3, "Mirror");
  label[4] = new TLatex(0.01, 0.97-0.03*4, "Diffuser");
  label[5] = new TLatex(0.01, 0.97-0.03*5, "Absorber");
  label[6] = new TLatex(0.01, 0.97-0.03*6, "Detector");
  label[7] = new TLatex(0.01, 0.97-0.03*7, "Mixture");
  label[8] = new TLatex(0.57, 0.97, "Simulation result");
  label[9] = new TLatex(0.57, 0.97-0.03*1, "Go out of world volume");
  label[10] = new TLatex(0.57, 0.97-0.03*2, "Exceed limit of reflections");
  label[11] = new TLatex(0.57, 0.97-0.03*3, "Absorved in normal medium");
  label[12] = new TLatex(0.57, 0.97-0.03*4, "Absorbed by absorber");
  label[13] = new TLatex(0.57, 0.97-0.03*5, "Detected by detector");
  label[14] = new TLatex(0.92, 0.97-0.03*1, "0");
  label[15] = new TLatex(0.92, 0.97-0.03*2, "0");
  label[16] = new TLatex(0.92, 0.97-0.03*3, "0");
  label[17] = new TLatex(0.92, 0.97-0.03*4, "0");
  label[18] = new TLatex(0.92, 0.97-0.03*5, "0");
  label[1]->SetTextColor(2);
  label[2]->SetTextColor(6);
  label[3]->SetTextColor(800);
  label[4]->SetTextColor(8);
  label[5]->SetTextColor(kCyan+2);
  label[6]->SetTextColor(4);
  label[7]->SetTextColor(28);
  for(int i=0;i<19;i++){
    label[i]->SetNDC();
    label[i]->SetTextSize(0.03);
  }
}

Gui::~Gui() {
  Cleanup();
}

void Gui::SetSimulator(Simulator *SIM) {
  sim = SIM;
  if (sim->GetSource()->ChargedMode()) {
    particlesLabel->SetText("Particles");
  }
  c1 = Canvas->GetCanvas();
  c1->cd();
  view = TView::CreateView(1);
  SetGeometry();
  DrawGeometry();
  c1->Update();  
}

void Gui::DoRun() {
  c1->Clear();
  c1->cd();
  view = TView::CreateView(1);
  DrawGeometry();
  view->ZoomView(nullptr,zoom);
  if (sim) {
    sim->SetNevt(NumEntry->GetNumber());
    sim->Run();
    count=sim->GetResult();
    for(int i=0;i<5;i++){
      std::string text = std::to_string(count[i]);
      label[14+i]->SetText(0.92, 0.97 - 0.03 * (i+1), text.c_str());
    }
    DrawLabel();
  }
  c1->Update();  
}

void Gui::DoSave() {
  Canvas->GetCanvas()->Print(FileNameEntry->GetText());
  std::cout << "Saved display to " << FileNameEntry->GetText() << std::endl;
}

void Gui::DoZoom() {
  view->ZoomView(nullptr,1.25);
  zoom*=1.25;
  c1->Update();
}

void Gui::DoUnzoom() {
  view->ZoomView(nullptr,1./1.25);
  zoom/=1.25;
  c1->Update();
}

void Gui::DoClear() {
  c1->Clear();
  c1->cd();
  view = TView::CreateView(1);
  DrawGeometry();
  view->ZoomView(nullptr,zoom);
  for(int i=0;i<5;i++){
    label[14+i]->SetText(0.92, 0.97 - 0.03 * (i+1), "0");
  }
  DrawLabel();
  c1->Update();
}

void Gui::DoExit() {
  CloseWindow();
}
void Gui::OnToggle(Bool_t state) {
  showLabel = state;
  if (showLabel) {
    DrawLabel();
  }
  else { 
    DoClear();
  }
  c1->Update();
}
void Gui::UpdateDisplay() {
  Canvas->GetCanvas()->Modified();
  Canvas->GetCanvas()->Update();
}
void Gui::DrawTrack(const Position& p0, const Position& p1, bool charged){//Draw track of optical photon in the event display
  TPolyLine3D *l = new TPolyLine3D(2);
  l->SetPoint(0,p0[0],p0[1],p0[2]);
  l->SetPoint(1,p1[0],p1[1],p1[2]);
  if(charged){
    l->SetLineColor(633);
    l->SetLineWidth(3);
  }
  else{
    l->SetLineColor(1);
    l->SetLineWidth(2);
  }
  l->Draw();
}
void Gui::DrawLabel(){
  if (showLabel) {
    for(int i=0;i<19;i++){
      label[i]->Draw();
    }
  }
}
bool Gui::ComparePosition(const Position& p0, const Position& p1){
    for (int i = 0; i < 3; ++i) {
        if (p0[i] != p1[i]) return p0[i] < p1[i];
    }
    return false;
}
bool Gui::Parallel(const Direction& n0, const Direction& n1){
  return fabs(n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2]) > 0.999;
}
void Gui::SetGeometry(){
  x_min = y_min = z_min = world;
  x_max = y_max = z_max = -world;
  std::map<std::pair<Position, Position>, std::vector<Direction>> edge_to_norm;  
  for(int m=0;m<sim->GetNMat();m++){//loop for Materials
    edge_to_norm.clear();
    for(int t=0;t<sim->GetMaterial(m)->NTriangle();t++){//loop for Triangles
      for(int j=0;j<3;j++){
        Position p1 = { sim->GetMaterial(m)->GetTriangle(t).X(j),       sim->GetMaterial(m)->GetTriangle(t).Y(j),       sim->GetMaterial(m)->GetTriangle(t).Z(j) };
        Position p2 = { sim->GetMaterial(m)->GetTriangle(t).X((j+1)%3), sim->GetMaterial(m)->GetTriangle(t).Y((j+1)%3), sim->GetMaterial(m)->GetTriangle(t).Z((j+1)%3) };
        if(ComparePosition(p1,p2))std::swap(p1, p2);
        edge_to_norm[{p1, p2}].push_back(sim->GetMaterial(m)->GetTriangle(t).GetNormal());
        ::Compare(x_max, x_min, p1[0]);
        ::Compare(y_max, y_min, p1[1]);
        ::Compare(z_max, z_min, p1[2]);
      }//end of for
    }//end of for
    for (auto it = edge_to_norm.begin(); it != edge_to_norm.end(); ++it) {
      if (it->second.size() != 2) continue;
      if (!Parallel(it->second[0], it->second[1])){
	geo.push_back(std::make_tuple(it->first.first,it->first.second,sim->GetMaterial(m)->TType()));	
      }
    }
  }//end of for
  r_max = (x_max-x_min > y_max-y_min) ? x_max-x_min : y_max-y_min;
  r_max = (z_max-z_min > r_max) ? z_max-z_min : r_max;
  view->SetRange((x_max+x_min-r_max)/2, (y_max+y_min-r_max)/2, (z_max+z_min-r_max)/2,
		 (x_max+x_min+r_max)/2, (y_max+y_min+r_max)/2, (z_max+z_min+r_max)/2);
}
void Gui::DrawGeometry(){
  int type;
  Position p0, p1;
  for(unsigned int i=0;i<geo.size();i++){
    p0=std::get<0>(geo[i]);
    p1=std::get<1>(geo[i]);
    type=std::get<2>(geo[i]);
    TPolyLine3D *l = new TPolyLine3D(2);
    l->SetPoint(0,p0[0],p0[1],p0[2]);
    l->SetPoint(1,p1[0],p1[1],p1[2]);
    if(type==0)     l->SetLineColor(2);
    else if(type==1)l->SetLineColor(6);
    else if(type==2)l->SetLineColor(800);
    else if(type==3)l->SetLineColor(8);
    else if(type==4)l->SetLineColor(kCyan+2);
    else if(type==5)l->SetLineColor(4);
    else if(type>=6)l->SetLineColor(28);
    l->Draw();
  }
  view->SetRange((x_max+x_min-r_max)/2, (y_max+y_min-r_max)/2, (z_max+z_min-r_max)/2,
		 (x_max+x_min+r_max)/2, (y_max+y_min+r_max)/2, (z_max+z_min+r_max)/2);
}
void Gui::CloseWindow(){
  std::cout<<"GUI window was closed."<<std::endl;
  app->Terminate();
}
