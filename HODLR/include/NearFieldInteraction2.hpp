class NearFieldInteraction2: public Quadrature {
  /*
  this class doesn't depend on the tree data
  */
public:
  int degreeOfBases; //p(one minus)
  int NumberOfBases;
  std::vector<orderedPair> degreeOfBasesVec;
  int pMax;
  int nChebNodes;
  int rank;
  double epsilon;
  std::vector<double> standardChebNodes1D;
  std::vector<pts2D> standardChebNodes;
  int NLEVELS; //upper limit for FMM tree hight
  integrandInputs II;

  NearFieldInteraction2(int nChebNodes, integrandInputs II, int degreeOfBases): Quadrature() {
    this->nChebNodes = nChebNodes;
    this->rank = nChebNodes*nChebNodes;
    this->II = II;
    nNodes = 10;
    this->epsilon = 1e-4;
    this->NLEVELS = 10;
    this->degreeOfBases = degreeOfBases; //p(one minus)
    this->NumberOfBases = int((degreeOfBases)*(degreeOfBases+1)/2);
    for (size_t d = 0; d < degreeOfBases; d++) {
      for (size_t a = 0; a < d+1; a++) {
        int b = d-a;
        orderedPair op;
        op.x = a;
        op.y = b;
        degreeOfBasesVec.push_back(op);
      }
    }
    this->pMax = 60;
    set_Standard_Cheb_Nodes();
  };

  double factorial(int n) {
    return (n==1 || n==0) ? 1: double(n) * factorial(n - 1);
  }

  double ap(int p, double z) {
    return pow(-1, p)*pow(pow(z, p)/double(factorial(p)), 2.0);
  }

  double gp(int p, double z) {
    double temp = GAMMA + log(z) - Hp(p);
    temp *= ap(p, z);
    return temp;
  }

  std::complex<double> dp(int p, double z) {
    return 2.0*I*ap(p, z)/PI;
  }

  std::complex<double> cp(int p, double z) {
    return ap(p, z) * 2.0*I*gp(p, z)/PI;
  }

  double Hp(int p) {
    return (p==0) ? 0 : Hp(p-1)+1.0/double(p);
  }

  void evaluateV2(Mat &V2) {
    //A is the leaf box, whose neighbor is B
    //level_B: level of box B;
    //B is neighbor of unit box which is a leaf,A
    //leaf can be at any level; can be nLevels or anything less
    //Make set of operators corresponding to leaves at various levels
    //p is index of summation term
    V2 = Mat::Zero(rank, NumberOfBases);
    std::vector<pts2D> Nodes_A;
    shift_scale_Nodes(standardChebNodes, Nodes_A, 0.0, 0.0, II.radius_A);
    for (size_t i = 0; i < rank; i++) {
      for (size_t l = 0; l < NumberOfBases; l++) {
        for (size_t p = 0; p <= pMax; p++) {
          II.p = p;
          II.X.x = (Nodes_A[i].x-II.center_B.x)/II.radius_B;
          II.X.y = (Nodes_A[i].y-II.center_B.y)/II.radius_B;
          II.l = l;
          xcenter = 0.0; //of Quadratire object
          ycenter = 0.0; //of Quadratire object
          Lx = 1.0; //of Quadratire object
          Ly = 1.0; //of Quadratire object
          //obtain_Quadrature(epsilon);//when an integrand is evaluated gauss Quadrature tree is constructed fresh
					obtain_Quadrature_Uniform();//when an integrand is evaluated gauss Quadrature tree is constructed fresh
          V2(i,l) += dp(p, kappa*II.radius_B) * integral;
        }
      }
    }
  }

  dtype integrand(pts2D P) {
    double R2 = (II.X.x-P.x)*(II.X.x-P.x) + (II.X.y-P.y)*(II.X.y-P.y);
    double R = sqrt(R2);
    orderedPair op = degreeOfBasesVec[II.l];
    int a = op.x;
    int b = op.y;
    return pow(R/2.0, 2.0*II.p)*log(R/2.0)*get_ChebPoly(P.x, a)*get_ChebPoly(P.y, b);
	};

  //	get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}

  //	set_Standard_Cheb_Nodes
	void set_Standard_Cheb_Nodes() {
		for (int k=0; k<nChebNodes; ++k) {
			standardChebNodes1D.push_back(-cos((k+0.5)/nChebNodes*PI));
		}
		pts2D temp1;
		for (int j=0; j<nChebNodes; ++j) {
			for (int k=0; k<nChebNodes; ++k) {
				temp1.x	=	standardChebNodes1D[k];
				temp1.y	=	standardChebNodes1D[j];
				standardChebNodes.push_back(temp1);
			}
		}
  }

  void shift_scale_Nodes(std::vector<pts2D>& Nodes, std::vector<pts2D>& shifted_scaled_Nodes, double xShift, double yShift, double radius) {
		for (int k=0; k < Nodes.size(); ++k) {
			pts2D temp;
			temp.x	=	radius*Nodes[k].x+xShift;
			temp.y	=	radius*Nodes[k].y+yShift;
			shifted_scaled_Nodes.push_back(temp);
		}
	}
};
