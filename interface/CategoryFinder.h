#ifndef CategoryFinder_H
#define CategoryFinder_H



using namespace edm;
using namespace std;


class CategoryFinder {
  public :
    
    CategoryFinder() {};
  ~CategoryFinder() {};
  
  double    pMin;
  double    pMax;
  double    etaMin;
  double    etaMax;
  int       nHitsMin;
  int       nHitsMax;
  int       nPixelHitsMin;
  int       nPixelHitsMax;
  double    chiMin;
  double    chiMax;
  int       withFirstPixel;


};


#endif
