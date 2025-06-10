#include "Charged.hh"

Charged::Charged(std::mt19937 MT, std::vector<Material*> &MAT)
  :mt(MT),
   mat(MAT),
   act_sci(false),
   act_che(false)
{
  std::uniform_int_distribution<int> seed_dist(0, std::numeric_limits<int>::max());
  int seed = seed_dist(mt);
  rng.SetSeed(seed);
}
Charged::~Charged()
{
}
void Charged::Simulate(const Position& pos, const Direction& vec)
{
  particle.clear();
  for(int j=0;j<3;j++)far[j]=pos[j]+vec[j]*world;//far point on the straight-line of the track
  for(unsigned int m=0;m<mat.size();m++){//loop for Materials
    if(mat[m]->TType()<2||mat[m]->TType()==6){//Check if the material type is medium, converter or mixture including medium or converter
      if(act_sci&&mat[m]->Is_Scinti())p_sci=1;
      else p_sci=0;
      if(act_che&&mat[m]->Index()*beta>1)p_che=CalcCheProb(mat[m]->Index());
      else p_che=0;
      if(p_sci==0&&p_che==0)continue;

      //Count the number of intersections between the object and the charged particle's trajectory
      cross.clear();
      for(int t=0;t<mat[m]->NTriangle();t++){//loop for Triangles
	if(mat[m]->GetTriangle(t).Collision(pos,far,crs)){//Check if a line between two points collides with the triangle or not
	  cross.push_back(crs);	  	 
	}
      }

      //Sort the intersection points in order of distance from the particle's origin.
      std::sort(cross.begin(), cross.end(), [this, &pos](const Position& a, const Position& b) {
					      return Distance(pos, a) < Distance(pos, b);
					    });

      //If the number of intersections is odd, add the interaction point (an odd count means the interaction point is inside the object).
      if(cross.size()%2==1){
	cross.insert(cross.begin(),pos);
      }

      for(unsigned int t=0;t<cross.size()/2;t++){
	dist=Distance(cross[2*t],cross[2*t+1]);//Construct the trajectory segments inside the object from pairs of intersection points
	n_sci = NumPhoton(dist,p_sci,true);
	n_che = NumPhoton(dist,p_che,false);
	for(int p=0;p<n_sci;p++){//Generate scintillation photons along the path
	  s = unirand(mt);
	  ptn = Interpolate(cross[2*t], cross[2*t+1], s);
	  dir = Isotropic();
	  time = Distance(pos,ptn)/(c_0*beta)+ScintiDelay(mt);
	  particle.push_back(std::make_tuple(ptn, dir, time));
	}
	for(int p=0;p<n_che;p++){//Generate Cherenkov photons along the path
	  s = unirand(mt);
	  ptn = Interpolate(cross[2*t], cross[2*t+1], s);
	  dir = CherenkovDir(mat[m]->Index(),vec);
	  time = Distance(pos,ptn)/(c_0*beta);
	  particle.push_back(std::make_tuple(ptn, dir, time));
	}
      }      
    }
  }
}
double Charged::Generate(Position& pos, Position& vec, int n)
{
  pos=std::get<0>(particle[n]);
  vec=std::get<1>(particle[n]);
  return std::get<2>(particle[n]);
}
int Charged::GetNPhotons()
{
  return particle.size();
}
int Charged::NumPhoton(double d, double prob, bool is_scinti=true)
{
  if(prob==0)return 0;
  double n_exp;
  if(is_scinti){//Scintillation case
    double mean_dE = d * dedx;
    double width = d * width_pp;
    double dE_sampled = rng.Landau(mean_dE, width);//Landau distribution
    if (dE_sampled < 0) dE_sampled = 0;
    n_exp = dE_sampled * yield;
  }
  else{//Cherenkov case
    n_exp = d*prob;
  }
  std::poisson_distribution<int> poisson_dist(n_exp);
  return poisson_dist(mt);
}
double Charged::Distance(const Position& p1, const Position& p2)
{
    return std::sqrt(
        std::pow(p1[0] - p2[0], 2) +
        std::pow(p1[1] - p2[1], 2) +
        std::pow(p1[2] - p2[2], 2)
    );
}
Position Charged::Interpolate(const Position& p1, const Position& p2, double t) {
    return {
        p1[0] + t * (p2[0] - p1[0]),
        p1[1] + t * (p2[1] - p1[1]),
        p1[2] + t * (p2[2] - p1[2])
    };
}
Direction Charged::Isotropic(){//Randomely and isotropically determine the direction
  Direction v;
  v[2]=1-2*unirand(mt);
  double ang=unirand(mt)*2*pi;
  double vr=sqrt(1-v[2]*v[2]);
  v[0]=vr*cos(ang);
  v[1]=vr*sin(ang);
  return v;
}

Direction Charged::CherenkovDir(double index, const Direction& vec){//Determine the direction of Cherenkov light
  Direction v,vche;
  double sint,cost,sinp,cosp;
  cost=vec[2];
  sint=sqrt(1-cost*cost);
  if(sint!=0){
    cosp=vec[0]/sint;
  }
  else{
    cosp=1;
  }
  if(vec[1] >= 0){
    sinp=sqrt(1-cosp*cosp);
  }
  else{
    sinp=-sqrt(1-cosp*cosp);
  }

  v[2]=1/index/beta;//Cosine of Cherenkov angle
  double ang=unirand(mt)*2*pi;
  double vr=sqrt(1-v[2]*v[2]);
  v[0]=vr*cos(ang);
  v[1]=vr*sin(ang);
  vche[0]=-sinp*v[0] +cost*cosp*v[1] +sint*cosp*v[2];
  vche[1]=cosp*v[0]  +cost*sinp*v[1] +sint*sinp*v[2];
  vche[2]=           -sint*v[1]      +cost*v[2];  
  return vche;
}
void Charged::SetBeta(double BETA){
  beta=BETA;
}
void Charged::SetScinti(std::string sci_type, double YIELD){
  act_sci=true;
  yield=YIELD;
  std::map<std::string, std::tuple<double, double, double, double, double>> scintillators = {
											     {"organic", {5.7, 11.4, 1.032, 64e-6, 10000}},
											     {"NaI",     {32.0, 74.5, 3.67, 528e-6, 38000}},
											     {"CsI",     {54.0, 132.5, 4.51, 560e-6, 54000}},
											     {"BGO",     {75.0, 208.0, 7.13, 710e-6, 9000}},
											     {"LYSO",    {66.0, 174.0, 7.10, 730e-6, 32000}},
											     {"LaBr3",   {47.0, 102.0, 5.08, 470e-6, 63000}}
  };
  Z = std::get<0>(scintillators[sci_type]);
  A = std::get<1>(scintillators[sci_type]);
  density = std::get<2>(scintillators[sci_type]);
  I = std::get<3>(scintillators[sci_type]);
  if(yield==0)yield=std::get<4>(scintillators[sci_type]);
  CalcSciPar();
}
void Charged::SetCherenkov(double WLMIN, double WLMAX){
  act_che=true;
  wlmin=WLMIN;
  wlmax=WLMAX;
}

double Charged::CalcCheProb(double n){
  const double alpha = 1.0 / 137.035999; // Fine-structure constant
  const double nm_to_cm = 1e-7;
  double wlmin_cm = wlmin * nm_to_cm;
  double wlmax_cm = wlmax * nm_to_cm;
  double coeff = 2.0 * pi * alpha;
  double factor = 1.0 - 1.0 / (beta * beta * n * n);
  double integral = (1.0 / wlmin_cm) - (1.0 / wlmax_cm);
  double photons_per_cm = coeff * factor * integral;
  return photons_per_cm / 10.0;
}
void Charged::CalcSciPar(){
  const double K = 0.307075; // MeV mol^-1 cm^2
  const double me = 0.511;   // MeV/c^2
  double beta2 = beta * beta;
  double gamma= 1 / sqrt(1-beta2);
  double gamma2 = gamma * gamma;
  double Wmax = (2.0 * me * beta2 * gamma2) /
                  (1.0 + 2.0 * gamma * me / 105.658 + std::pow(me / 105.658, 2));
  double arg = (2.0 * me * beta2 * gamma2 * Wmax) / (I * I);
  dedx = K * (Z / A) * density * (1.0 / beta2) * (0.5 * std::log(arg) - beta2) / 10.; // dE/dx (MeV/cm)
  width_pp = 0.5 * K * (Z / A) * density * (1.0 / beta2) / 10.; // scaling parameter per path (MeV/mm)
}
