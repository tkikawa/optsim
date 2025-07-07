#include "Charged.hh"

Charged::Charged(std::mt19937& MT, std::vector<Material*> &MAT)
  :mt(MT),
   mat(MAT),
   act_sci(false),
   act_che(false),
   is_electron(false)
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
  double p_sci,p_che,dist,s, time;
  int n_sci,n_che;
  Position crs,ptn,far;
  Direction dir, plr;
  std::vector<Position> cross;
  
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
	  ScintiDir(dir,plr);
	  time = Distance(pos,ptn)/(c_0*beta)+ScintiDelay();
	  particle.push_back(std::make_tuple(ptn, dir, plr, time));
	}
	for(int p=0;p<n_che;p++){//Generate Cherenkov photons along the path
	  s = unirand(mt);
	  ptn = Interpolate(cross[2*t], cross[2*t+1], s);
	  CherenkovDir(mat[m]->Index(),vec,dir,plr);
	  time = Distance(pos,ptn)/(c_0*beta);
	  particle.push_back(std::make_tuple(ptn, dir, plr, time));
	}
      }      
    }
  }
}

void Charged::Generate(Position& pos, Position& vec, Position& pol, double& t, int n)
{
  pos=std::get<0>(particle[n]);
  vec=std::get<1>(particle[n]);
  pol=std::get<2>(particle[n]);
  t=std::get<3>(particle[n]);
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
    if (dE_sampled > 5.5 * mean_dE) dE_sampled = rng.Landau(mean_dE, width);
    if (dE_sampled < 0) dE_sampled = 0;
    n_exp = dE_sampled * yield;
  }
  else{//Cherenkov case
    n_exp = d*prob;
  }
  std::poisson_distribution<int> poisson_dist(n_exp);
  return poisson_dist(mt);
}
Position Charged::Interpolate(const Position& p1, const Position& p2, double t) {
    return {
        p1[0] + t * (p2[0] - p1[0]),
        p1[1] + t * (p2[1] - p1[1]),
        p1[2] + t * (p2[2] - p1[2])
    };
}

void Charged::ScintiDir(Direction& vsci, Direction& psci){
  Isotropic(mt,vsci);// Generate isotropic photon direction vector (vsci)
  RandomPolarization(mt,vsci, psci); // Generate random linear polarization vector in the plane orthogonal to vsci
}

void Charged::CherenkovDir(double index, const Direction& vec, Direction& vche, Direction& pche){
  Direction v;
  double sint, cost, sinp, cosp;
  cost = vec[2];
  sint = sqrt(1.0 - cost * cost);

  if(sint != 0){
    cosp = vec[0] / sint;
  } else {
    cosp = 1.0;
  }

  if(vec[1] >= 0){
    sinp = sqrt(1.0 - cosp * cosp);
  } else {
    sinp = -sqrt(1.0 - cosp * cosp);
  }

  v[2] = 1.0 / (index * beta); // cos(theta_c)
  double ang = unirand(mt) * 2 * pi;
  double vr = sqrt(1.0 - v[2] * v[2]);
  v[0] = vr * cos(ang);
  v[1] = vr * sin(ang);

  // --- Apply rotation to obtain global Cherenkov direction (vche) ---
  vche[0] = -sinp * v[0] + cost * cosp * v[1] + sint * cosp * v[2];
  vche[1] =  cosp * v[0] + cost * sinp * v[1] + sint * sinp * v[2];
  vche[2] =                 -sint * v[1]      + cost * v[2];

  Normalize(vche);

  // --- Polarization vector is orthogonal to both vche and vec (pche = vec × vche) ---
  Cross(vec, vche, pche);
  Normalize(pche);
}

void Charged::SetBeta(double BETA){
  beta=BETA;
}

void Charged::SetElectron(bool is_e){
  is_electron = is_e;
}

void Charged::SetScinti(std::string sci_type, double YIELD, double DELAY){
  act_sci=true;
  yield = YIELD;
  scinti_lifetime = DELAY;
  std::map<std::string, std::tuple<double, double, double, double, double, double>> scintillators
    = {//Name, Z, A, density, I, light_yield, life_time
	{"organic", {5.7, 11.4, 1.032, 64.e-6, 10000., 2.1}},
	{"NaI",     {32.0, 74.5, 3.67, 528.e-6, 38000., 230.}},
	{"CsI",     {54.0, 132.5, 4.51, 560.e-6, 65000., 1000.}},
	{"BGO",     {75.0, 208.0, 7.13, 710.e-6, 8200., 300.}},
	{"LYSO",    {66.0, 174.0, 7.10, 730.e-6, 25000., 50.}},
	{"LaBr3",   {47.0, 102.0, 5.08, 470.e-6, 63000., 26.}}
  };
  Z = std::get<0>(scintillators[sci_type]);
  A = std::get<1>(scintillators[sci_type]);
  density = std::get<2>(scintillators[sci_type]);
  I = std::get<3>(scintillators[sci_type]);
  if(yield==0)yield = std::get<4>(scintillators[sci_type]);
  if(scinti_lifetime==-99999)scinti_lifetime = std::get<5>(scintillators[sci_type]);
  CalcSciPar();
}

void Charged::SetCherenkov(double WLMIN, double WLMAX){
  act_che=true;
  wlmin=WLMIN;
  wlmax=WLMAX;
}
double Charged::CalcCheProb(double n){// Calculate expected Chrenkov photon yield following Frank–Tamm formula
  const double alpha = 1.0 / 137.035999; // Fine-structure constant
  const double nm_to_cm = 1e-7;
  double wlmin_cm = wlmin * nm_to_cm;
  double wlmax_cm = wlmax * nm_to_cm;
  double coeff = 2.0 * pi * alpha;
  double factor = 1.0 - 1.0 / (beta * beta * n * n);
  double integral = (1.0 / wlmin_cm) - (1.0 / wlmax_cm);
  return coeff * factor * integral / 10.0;// Photon yield per unit length (photons/mm)
}

void Charged::CalcSciPar(){
  const double K = 0.307075; // MeV mol^-1 cm^2
  const double me = 0.511;   // MeV/c^2
  double beta2 = beta * beta;
  double gamma= 1 / sqrt(1-beta2);
  double X = log10(beta*gamma);
  double delta = 0.0;

  //Sternheimer paramters
  double X0 = 0.2;
  double X1 = 3.0;
  double C = 4.0;
  
  if(!is_electron){//Calculate dE/dx following Bethe-Broch formula.
  double gamma2 = gamma * gamma;
  double Wmax = (2.0 * me * beta2 * gamma2) /
                  (1.0 + 2.0 * gamma * me / 105.658 + std::pow(me / 105.658, 2));
  double arg = (2.0 * me * beta2 * gamma2 * Wmax) / (I * I);

  // Density effect correction
  if (X < X0) {
    delta = 0.0;
  } else if (X < X1) {
    // Example: parametrization (refer to ICRU or PDG recommended data)
    delta = 2 * log(10) * (X - X0); // Example approximation
  } else {
    delta = 2 * log(10) * X - C; // High-energy asymptotic behavior
  }

  // Shell correction C/Z (approx: ~0.2 for light elements)
  double shell_correction = 0.2; // Should be adjusted using empirical formula by Z
  
  dedx = K * (Z / A) * density * (1.0 / beta2) * (0.5 * std::log(arg) - beta2 - delta / 2.0 - shell_correction / Z) / 10.; // dE/dx (MeV/mm)
  }
  else{//Calculate dE/dx following Berger-Seltzer formula.
    double tau = gamma - 1.0;
    double T = tau * me;
    double F = 1.0 - beta2 + std::log(1.0 + tau * tau / 2.0); // F(tau): correction due to indistinguishability

    // Density effect correction
    if (X < X0) {
      delta = 0.0;
    } else if (X < X1) {
      delta = 2.0 * log(10.0) * pow(X - X0, 2);
    } else {
      delta = 2.0 * log(10.0) * X - C;
    }
    
    dedx = K * (Z / A) * density * (1.0 / beta2) * (std::log(T / I) + F - delta - 2.0) /10.;// dE/dx (MeV/mm)
  }

  width_pp = 0.5 * K * (Z / A) * density * (1.0 / beta2) / 10.; // scaling parameter per path (MeV/mm)
}

double Charged::ScintiDelay(){
  if(scinti_lifetime<=0) return 0;
  else{
    int skip = mt() % 7 + 3;
    mt.discard(skip);
    return exprand(mt)*scinti_lifetime;
  }
}
