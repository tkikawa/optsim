#include "Charged.hh"

Charged::Charged(std::mt19937 MT, std::vector<Material*> &MAT)
  :mt(MT),
   mat(MAT)
{
}
Charged::~Charged()
{
}
void Charged::Simulate(const Position& pos, const Direction& vec)
{
  particle.clear();
  for(int j=0;j<3;j++)far[j]=pos[j]+vec[j]*world;//far point on the straight-line of the track
  for(unsigned int m=0;m<mat.size();m++){//loop for Materials
    if(mat[m]->Type()<2){
      cross.clear();
      p_sci=mat[m]->SciProb();
      p_che=mat[m]->CheProb();
      if(p_sci==0&&p_che==0)continue;
      for(int t=0;t<mat[m]->NTriangle();t++){//loop for Triangles
	if(mat[m]->GetTriangle(t).Collision(pos,far,crs)){//Check if a line between two points collides with the triangle or not
	  cross.push_back(crs);	  	 
	}
      }
      std::sort(cross.begin(), cross.end(), [this, &pos](const Position& a, const Position& b) {
					      return Distance(pos, a) < Distance(pos, b);
					    });
      if(cross.size()%2==1){
	cross.insert(cross.begin(),pos);
      }
      for(unsigned int t=0;t<cross.size()/2;t++){
	dist=Distance(cross[2*t],cross[2*t+1]);
        std::poisson_distribution<int> poisson_sci(dist*p_sci);
	std::poisson_distribution<int> poisson_che(dist*p_che);
        n_sci = poisson_sci(mt);
	n_che = poisson_che(mt);
	for(int p=0;p<n_sci;p++){
	  s = unirand(mt);
	  ptn = Interpolate(cross[2*t], cross[2*t+1], s);
	  dir = Isotropic();
	  particle.push_back(std::make_pair(ptn, dir));
	}
	for(int p=0;p<n_che;p++){
	  s = unirand(mt);
	  ptn = Interpolate(cross[2*t], cross[2*t+1], s);
	  dir = CherenkovDir(mat[m]->Index(),vec);
	  particle.push_back(std::make_pair(ptn, dir));
	}
      }      
    }
  }
}
void Charged::Generate(Position& pos, Position& vec, int n)
{
  pos=particle[n].first;
  vec=particle[n].second;
}
int Charged::GetNPhotons()
{
  return particle.size();
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

  v[2]=1/index;//Cosine of Cherenkov angle assuming beta~1.
  double ang=unirand(mt)*2*pi;
  double vr=sqrt(1-v[2]*v[2]);
  v[0]=vr*cos(ang);
  v[1]=vr*sin(ang);
  vche[0]=-sinp*v[0] +cost*cosp*v[1] +sint*cosp*v[2];
  vche[1]=cosp*v[0]  +cost*sinp*v[1] +sint*sinp*v[2];
  vche[2]=           -sint*v[1]      +cost*v[2];
  
  return vche;
}
