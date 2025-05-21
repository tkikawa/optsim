#include "Geometry.hh"

Geometry::Geometry(std::mt19937 MT)
  : mt(MT)
{
}
Geometry::~Geometry()
{
}
void Geometry::LoadCAD(std::string name){//Read out the input CAD file as a assembly of trianbles

  std::ifstream cadfile(name);
  if(cadfile.good()){
    cadfile.close();
  }
  else{
    std::cerr<<"Error: CAD file "<<name.c_str()<<" is not found."<<std::endl;
    exit(1);
  }

  Assimp::Importer importer;
  const aiScene* scene;
  aiMesh* aim;
  
  scene = importer.ReadFile(name,
			    aiProcess_Triangulate           |
			    aiProcess_JoinIdenticalVertices |
			    aiProcess_CalcTangentSpace);

  aim = scene->mMeshes[0];

  double vertex[3][3];
  for(unsigned int i=0; i < aim->mNumFaces; i++){
    const aiFace& face = aim->mFaces[i];
    for(int j=0;j<3;j++){
      vertex[j][0]=aim->mVertices[face.mIndices[j]].x*cadunit;
      vertex[j][1]=aim->mVertices[face.mIndices[j]].y*cadunit;
      vertex[j][2]=aim->mVertices[face.mIndices[j]].z*cadunit;
    }
    Triangle *tri = new Triangle(vertex);
    triangle.push_back(*tri);
  }
  
}
bool Geometry::InSolid(const Position& pos){//Check if the point is inside or outside of the geometry.
  Position far, tmp;
  far[0]=pos[0]; far[1]=pos[1]; far[2]=world;
  int ncol = 0;
  for(unsigned int i=0; i < triangle.size(); i++){
    if(triangle[i].Collision(pos,far,tmp))ncol++;
  }
  if(ncol%2 == 1) return true;
  else            return false;
}
