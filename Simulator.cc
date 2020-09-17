#include "Simulator.hh"

Simulator::Simulator(std::mt19937 MT, Config config, std::string OUTPUT)
{
  mt = MT;
  output = OUTPUT;
  nevt=10000;
  index0=1;
  std::istringstream(config["Global"]["Number"]) >> nevt;
  std::istringstream(config["Global"]["Index"]) >> index0;
  displaymode=false;
}
Simulator::~Simulator()
{
}
void Simulator::AddMaterial(Material *MAT)
{
  mat.push_back(*MAT);
}
void Simulator::SetSource(Source *SRC)
{
  src = SRC;
}
void Simulator::Run()
{
  file = new TFile(output.c_str(), "recreate");
  tree = new TTree("tree","tree");
  tree->Branch("ipos[3]",&ipos,"ipos[3]/D");
  tree->Branch("fpos[3]",&fpos,"fpos[3]/D");
  tree->Branch("ivec[3]",&ivec,"ivec[3]/D");
  tree->Branch("fvec[3]",&fvec,"fvec[3]/D");
  tree->Branch("time",&time,"time/D");
  tree->Branch("length",&length,"length/D");
  tree->Branch("imat",&imat,"imat/I");
  tree->Branch("fmat",&fmat,"fmat/I");
  tree->Branch("ftype",&ftype,"ftype/I");
  tree->Branch("nref",&nref,"nref/I");
  tree->Branch("npas",&npas,"npas/I");
  for(int j=0;j<5;j++)count[j]=0;
  std::cout<<"Generate "<<nevt<<" optical photons."<<std::endl;
  for(int i=0;i<nevt;i++){
    src->Generate(pos, vec);//Determine the initial position and direction
    mn=PointMaterial(pos);  //Check the material of initial position
    if(mn!=-1){//Initial position is in the defined material.
      matid=mat[mn].ID();
      index=mat[mn].Index();
      attlen=mat[mn].AttLen();
    }
    else{//Initial position is not in the defined materials and is in surrounding material (air).
      matid=0;
      index=index0;
      attlen=0;
    }
    Initialize();
    while(1){
      btype = -2;
      for(int j=0;j<3;j++)cross[j]=pos[j]+vec[j]*world;//far point on the straight-line of the track
      for(unsigned int m=0;m<mat.size();m++){//loop for Materials
	if(matid!=0&&matid<mat[m].ID()&&btype>-1){
	  goto ENDLOOP;
	}
	else if(matid!=0&&matid<mat[m].ID()&&btype==-1){
	  if(mat[m].InSolid(cross)){
	    newmatid=mat[m].ID();
	    newindex=mat[m].Index();
	    newattlen=mat[m].AttLen();
	    btype=mat[m].Type();
	    newmn=m;
	    goto ENDLOOP;
	  }
	  continue;
	}
	for(int t=0;t<mat[m].NTriangle();t++){//loop for Triangles
	  if(mat[m].GetTriangle(t).Collision(pos,cross,cand)){
	    for(int j=0;j<3;j++){
	      newpos[j]=cand[j];
	      if(matid!=0){
		cross[j]=cand[j];
	      }
	      else if(matid==mat[m].ID()){
		cross[j]=cand[j]+vec[j]*micro;
	      }
	      else{
		cross[j]=cand[j]-vec[j]*micro;
	      }
	    }
	    mat[m].GetTriangle(t).GetNormal(normal);
	    if(matid==mat[m].ID()){
	      newmatid=0;
	      newindex=index0;
	      newattlen=0;
	      btype=-1;
	      newmn=-1;
	    }
	    else{
	      newmatid=mat[m].ID();
	      newindex=mat[m].Index();
	      newattlen=mat[m].AttLen();
	      btype=mat[m].Type();
	      newmn=m;
	    }
	  }//end of if
	}//end of for
      }//end of for
      
    ENDLOOP:
      if(btype!=2){
	pl=sqrt(pow(newpos[0]-pos[0],2)+pow(newpos[1]-pos[1],2)+pow(newpos[2]-pos[2],2));
	if(attlen>0){//When attlen==0, absorption in the medium does not occur.
	  apl=-attlen*log(unirand(mt));
	  if(apl<pl){
	    pl=apl;
	    btype=-3;
	    for(int j=0;j<3;j++){
	      newpos[j]=pos[j]+vec[j]*pl;
	    }
	  }
	}
	length+=pl;
	time+=pl/c_0*index;
	if(displaymode){
	  Track(pos,newpos);
	}
      }
      if(btype >= -1 && btype <=2){
	if(btype == -1 || btype == 0){
	  if(Fresnel(vec,newvec,normal,index,newindex)){
	    mn=newmn;
	    matid=newmatid;
	    index=newindex;
	    attlen=newattlen;
	    npas++;
	  }
	  else{
	    nref++;
	  }
	}
	else if(btype == 1){
	  Specular(vec,newvec,normal);//Specular reflection
	  nref++;
	}
	else if(btype == 2){
	  Lambert(vec,newvec,normal);//Diffusion following the Lambert's cosine law
	  nref++;
	}
	for(int j=0;j<3;j++){
	  pos[j]=newpos[j]+newvec[j]*micro;
	  vec[j]=newvec[j];
	}
      }
      else break;//Absorbed, detected or go out of world volume.
      
    }//end of while

    for(int j=0;j<3;j++){
      fpos[j]=newpos[j];
      fvec[j]=newvec[j];
    }
    fmat=newmatid;
    ftype=FType(btype);//Conversion of type ID.
    tree->Fill();
    count[ftype]++;
  }//end of for
  Summary();  
  file->Write();
  file->Close();  
}
int Simulator::PointMaterial(double p[3]){
  for(unsigned int m=0;m<mat.size();m++){
    if(mat[m].InSolid(p)){
      return m;
    }
  }
  return -1;
}
void Simulator::Initialize(){
  for(int j=0;j<3;j++){
    ipos[j]=pos[j];
    ivec[j]=vec[j];
    newpos[j]=pos[j];
    newvec[j]=vec[j];
  }
  imat=matid;
  nref=0;
  npas=0;
  length=0;
  time=0;
}
int Simulator::FType(int bt){
  if(bt==-2)return 1;//Go out of world volume.
  else if(bt==-3)return 2;//Absorbed in normal medium.
  else return bt;//Absorbed by absorber(ftype=3) or detector(ftype=4).
}
void Simulator::Summary(){
  std::cout<<"Go out of world volume   : "<<count[1]<<std::endl;
  std::cout<<"Absorved in normal medium: "<<count[2]<<std::endl;
  std::cout<<"Absorbed by absorber     : "<<count[3]<<std::endl;
  std::cout<<"Detected by detector     : "<<count[4]<<std::endl;
}
bool Simulator::Fresnel(double v[3], double *newv, double norm[3], double idx_in, double idx_out){
  double idx_ratio=idx_out/idx_in;
  double cos_in=v[0]*norm[0]+v[1]*norm[1]+v[2]*norm[2];
  double sin_in=sqrt(1-cos_in*cos_in);
  if(cos_in>0){
    for(int i=0;i<3;i++){
      norm[i]*=-1;
    }
  }
  else{
    cos_in*=-1;
  }
  if(sin_in>idx_ratio){//Total reflection
    Specular(v,newv,norm);
    return false;
  }
  else{
    double sin_out=sin_in/idx_ratio;//Snell's law
    double cos_out=sqrt(1-sin_out*sin_out);
    double ref_prob=0.5*(pow(cos_in-idx_ratio*cos_out,2)/pow(cos_in+idx_ratio*cos_out,2)+pow(idx_ratio*cos_in-cos_out,2)/pow(idx_ratio*cos_in+cos_out,2));//Fresnel equations
    if(unirand(mt)<ref_prob){//Reflection
      Specular(v,newv,norm);
      return false;
    }
    else{//Transmission with refraction
      double para[3];
      for(int i=0;i<3;i++){
	para[i]=v[i]+norm[i]*cos_in;
	if(sin_in!=0){
	  newv[i]=-norm[i]*cos_out+para[i]/sin_in*sin_out;
	}
	else{
	  newv[i]=v[i];
	}
      }
      return true;
    }
  } 
}
void Simulator::Specular(double v[3], double *newv, double norm[3]){
  double cp=v[0]*norm[0]+v[1]*norm[1]+v[2]*norm[2];
  for(int i=0;i<3;i++){
    newv[i]=v[i]-2*norm[i]*cp;
  }
}
void Simulator::Lambert(double v[3], double *newv, double norm[3]){
  double cp=v[0]*norm[0]+v[1]*norm[1]+v[2]*norm[2];
  if(cp>0){
    for(int i=0;i<3;i++){
      norm[i]*=-1;
    }
  }

  double cost,sint,cosp,sinp;
  double ranv[3],vr,ang;
  
  cost=norm[2];
  sint=sqrt(1-cost*cost);
  if(sint!=0){
    cosp=norm[0]/sint;
  }
  else{
    cosp=1;
  }
  if(norm[1]>0){
    sinp=sqrt(1-cosp*cosp);
  }
  else{
    sinp=-sqrt(1-cosp*cosp);
  }

  ranv[2]=sqrt(unirand(mt));//sqrt is give because of the cosine law.
  ang=unirand(mt)*2*pi;
  vr=sqrt(1-ranv[2]*ranv[2]);
  ranv[0]=vr*cos(ang);
  ranv[1]=vr*sin(ang);
  newv[0]=-sinp*ranv[0] +cost*cosp*ranv[1] +sint*cosp*ranv[2];
  newv[1]=cosp*ranv[0]  +cost*sinp*ranv[1] +sint*sinp*ranv[2];
  newv[2]=              -sint*ranv[1]      +cost*ranv[2];
}
void Simulator::Display()
{
  displaymode=true;
  double vtx[3][3];
  TCanvas *c1 = new TCanvas("Event display","Event display",800,800);
  TView *view = TView::CreateView(1);
  x_min = y_min = z_min = world;
  x_max = y_max = z_max = -world;
  for(unsigned int m=0;m<mat.size();m++){//loop for Materials
    for(int t=0;t<mat[m].NTriangle();t++){//loop for Triangles
	for(int j=0;j<3;j++){
	  vtx[j][0]=mat[m].GetTriangle(t).X(j);
	  vtx[j][1]=mat[m].GetTriangle(t).Y(j);
	  vtx[j][2]=mat[m].GetTriangle(t).Z(j);
	  Compare(x_max, x_min, vtx[j][0]);
	  Compare(y_max, y_min, vtx[j][1]);
	  Compare(z_max, z_min, vtx[j][2]);
	}//end of for
	Draw(vtx,mat[m].Type());
    }//end of for
  }//end of for
  r_max = (x_max-x_min > y_max-y_min) ? x_max-x_min : y_max-y_min;
  r_max = (z_max-z_min > r_max) ? z_max-z_min : r_max;
  view->SetRange((x_max-x_min-r_max)/2, (y_max-y_min-r_max)/2, (z_max-z_min-r_max)/2,
		 (x_max-x_min+r_max)/2, (y_max-y_min+r_max)/2, (z_max-z_min+r_max)/2);

  c1->Update();
  std::cout<<"******** Display mode ********"<<std::endl;
  std::cout<<"Red:    normal medium material"<<std::endl;
  std::cout<<"Magenta:mirror material"<<std::endl;
  std::cout<<"Green:  diffuser material"<<std::endl;
  std::cout<<"Cyan:   absorber material"<<std::endl;
  std::cout<<"Blue:   detector material"<<std::endl;
  std::string msg;
  while(1){
    std::cout<<"Waiting for input."<<std::endl;
    std::cout<<"\"quit\" or \"exit\": Terminate the display mode."<<std::endl;
    std::cout<<"Positive integer: Generate photons and display tracks."<<std::endl;
    std::cin >> msg;
    if(msg=="quit" || msg=="exit")break;
    if(std::all_of(msg.cbegin(), msg.cend(), isdigit)){
      if(stoi(msg)>0){
	nevt=stoi(msg);
	Run();
	c1->Update();
      }
    }
  }
  c1->Close();
}
void Simulator::Draw(double vtx[3][3], int type){
  TPolyLine3D *l = new TPolyLine3D(3);
  l->SetPoint(0,vtx[0][0],vtx[0][1],vtx[0][2]);
  l->SetPoint(1,vtx[1][0],vtx[1][1],vtx[1][2]);
  l->SetPoint(2,vtx[2][0],vtx[2][1],vtx[2][2]);
  l->SetPoint(3,vtx[0][0],vtx[0][1],vtx[0][2]);
  if(type==0)     l->SetLineColor(2);
  else if(type==1)l->SetLineColor(6);
  else if(type==2)l->SetLineColor(8);
  else if(type==3)l->SetLineColor(7);
  else if(type==4)l->SetLineColor(4);
  l->Draw();
}
void Simulator::Track(double p0[3],double p1[3]){
  TPolyLine3D *l = new TPolyLine3D(2);
  l->SetPoint(0,p0[0],p0[1],p0[2]);
  l->SetPoint(1,p1[0],p1[1],p1[2]);
  l->SetLineColor(1);
  l->SetLineWidth(2);
  l->Draw();
}
void Simulator::Compare(double &A_max, double &A_min, double A){
  if(A_min > A) A_min = A;
  if(A_max < A) A_max = A;
}