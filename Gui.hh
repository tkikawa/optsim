#ifndef GUI_H
#define GUI_H

#include <vector>
#include "Global.hh"
#include <TCanvas.h>
#include <TView.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGNumberEntry.h>
#include <TGTextEntry.h>
#include <TGLabel.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TLatex.h>
#include <RQ_OBJECT.h>
#include <TApplication.h>

class Simulator;

class Gui : public TGMainFrame {
  
public:
  Gui(const TGWindow *p, UInt_t w, UInt_t h, TApplication* APP);
  virtual ~Gui();

  void SetSimulator(Simulator* sim);
  void DoRun();
  void DoSave();
  void DoZoom();
  void DoUnzoom();
  void DoClear();
  void OnToggle(Bool_t state);
  void DoExit();
  void UpdateDisplay();
  void DrawTrack(const Position& p0, const Position& p1, bool charged = false);
  void DrawLabel();
  bool ComparePosition(const Position& p0, const Position& p1);
  bool Parallel(const Direction& p0, const Direction& p1);
  void DrawGeometry();
  void SetGeometry();
  void CloseWindow();
  
private:
  TRootEmbeddedCanvas *Canvas;
  TGNumberEntry *NumEntry;
  TGTextEntry *FileNameEntry;
  TGCheckButton* onoffCheck;
  bool showLabel = false;
  Simulator *sim;
  double x_min, x_max, y_min, y_max, z_min, z_max, r_max;
  TCanvas *c1;
  TView *view;
  double zoom;
  std::vector<std::tuple<Position, Position, int>> geo;
  TGLabel *particlesLabel;
  TLatex* label[19];
  TApplication* app;
  std::string info;
  std::array<int, 5> count;

};
#endif
