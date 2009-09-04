#ifndef PFTools_h
#define PFTools_h

/**_________________________________________________________________
   class:   PFTools.h

 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)

 version $Id: PFTools.h,v 1.1 2008/08/11 06:06:21 yumiceva Exp $

________________________________________________________________**/


class WorkingPoint{
 public:
	WorkingPoint() {}
		//WorkingPoint(edm::InputTag t,  std::string n, double min, double max=0, std::map< std::string, double> list = 0) : intag_(t), alias_(n), min_(min), max_(max), wpmap_(list) {}
  WorkingPoint(edm::InputTag t,  std::string n, double min, double max=0, std::map< std::string, double> list = 0) : intag_(t), alias_(n), min_(min), max_(max), wpmap_(list) { if (list!=0) wpmap_ = list;}
	
	std::map<std::string, double > list() const {return wpmap_;}
	edm::InputTag inputTag() const {return intag_;}
	std::string alias () const {return alias_;}
	std::string name () const {return alias_;}
	double cut() const { return min_; }
	double Minimum() const { return min_;}
	double Maximum() const { return max_;}
	void print () const  {
		std::cout <<" Working point "<<name()<<" Input Tag "<<inputTag()<<" Cut "<< cut()<<std::endl;
	}
	
 private:
  edm::InputTag intag_;
  std::string alias_;
  double min_;
  double max_;
  std::map< std::string, double > wpmap_;
};


class PFTools {

  public:
	PFTools() {}
	virtual ~PFTools() {}
	
  private:
	

};

#endif
