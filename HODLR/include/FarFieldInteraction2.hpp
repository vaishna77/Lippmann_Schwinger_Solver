class FarFieldInteraction2: public Quadrature {
  /*
  this class doesn't depend on the tree data
  */
public:
  int degreeOfBases; //p(one minus)
  int NumberOfBases;
  std::vector<orderedPair> degreeOfBasesVec;
  int seriesLength;
  int nChebNodes;
  int rank;
  double epsilon;
  std::vector<double> standardChebNodes1D;
  std::vector<pts2D> standardChebNodes;
  int NLEVELS; //upper limit for FMM tree hight
  std::vector<Mat > C; //hold the matrix C for different levels
  double L; //level 0 semi length of box
  std::vector<double> boxRadius;
  Eigen::MatrixXd Q;
  Eigen::MatrixXd Q_pinv;
	CIntegrandInputs CII;

  FarFieldInteraction2(int nChebNodes, double L, int degreeOfBases): Quadrature() {
    this->nChebNodes = nChebNodes;
    this->rank = nChebNodes*nChebNodes;
    nNodes = 10;
    this->epsilon = 1e-3;
    this->NLEVELS = 10;
    this->L = L;
    this->degreeOfBases = degreeOfBases; //p(one minus)
    this->NumberOfBases = int((degreeOfBases)*(degreeOfBases+1)/2);
    boxRadius.push_back(L);
    for (size_t j = 1; j <= NLEVELS; j++) {
      boxRadius.push_back(boxRadius[j-1]/2.0);
    }
    for (size_t d = 0; d < degreeOfBases; d++) {
      for (size_t a = 0; a < d+1; a++) {
        int b = d-a;
        orderedPair op;
        op.x = a;
        op.y = b;
        degreeOfBasesVec.push_back(op);
      }
    }
    seriesLength = 15;
		set_Standard_Cheb_Nodes();
  };
  ~FarFieldInteraction2() {
  };


  /*dtype integrand(pts2D P) {
    int n = CII.n;
    int l = CII.l;
    int level = CII.level;
    double R2 = P.x*P.x + P.y*P.y;
    double normX = sqrt(R2);
    double argumentJ = kappa*boxRadius[level]*normX;
		double theta_y = atan2(P.y, P.x);
    orderedPair op = degreeOfBasesVec[l];
    int a = op.x;
    int b = op.y;
		return besselJ(n, argumentJ)*exp(-I*double(n)*theta_y)*get_ChebPoly(P.x, a)*get_ChebPoly(P.y, b);
	};*/



	dtype integrand(pts2D P) {
    int l = CII.l;
    int level = CII.level;
		pts2D P2;
		P2.x = CII.xDummy.x - P.x;
		P2.y = CII.xDummy.y - P.y;

		pts2D center = CII.center;
    double R2 = P2.x*P2.x + P2.y*P2.y;
    double normX = sqrt(R2);
    double argumentH = kappa*normX;
		pts2D argumentB;
		argumentB.x = (P.x-center.x)/boxRadius[level];
		argumentB.y = (P.y-center.y)/boxRadius[level];
		orderedPair op = degreeOfBasesVec[l];
    int a = op.x;
    int b = op.y;

    dtype H0_true =  besselJ(0, argumentH)+I*besselY(0, argumentH);
    return H0_true*get_ChebPoly(argumentB.x, a)*get_ChebPoly(argumentB.y, b);

    pts2D RHn;
    RHn.x = CII.xDummy.x - CII.center.x;
    RHn.y = CII.xDummy.y - CII.center.y;
    double RHn2 = RHn.x*RHn.x + RHn.y*RHn.y;
    double normRHn = kappa*sqrt(RHn2);
    pts2D RJn;
    RJn.x = P.x - CII.center.x;
    RJn.y = P.y - CII.center.y;
    double RJn2 = RJn.x*RJn.x + RJn.y*RJn.y;
    double normRJn = kappa*sqrt(RJn2);
    double theta_x = atan2(RHn.y, RHn.x);
    double theta_y = atan2(RJn.y, RJn.x);
    dtype H0_series = 0.0 + 0.0*I;
    for (int n = -seriesLength; n <= seriesLength; n++) {
      H0_series += (besselJ(n, normRHn)+I*besselY(n, normRHn))*besselJ(n, normRJn)*exp(I*double(n)*(theta_x-theta_y));
    }
    //cout << "normRJn: " << normRJn << endl;
    return H0_series*get_ChebPoly(argumentB.x, a)*get_ChebPoly(argumentB.y, b);
    //cout << (besselJ(0, argumentH)+I*besselY(0, argumentH))*get_ChebPoly(argumentB.x, a)*get_ChebPoly(argumentB.y, b) << endl;
	};

	std::complex<double> evaluateEntry(pts2D xDummy, int level, pts2D center, int j) {
		xcenter = center.x;
		ycenter = center.y;
		Lx = boxRadius[level];
		Ly = boxRadius[level];
		CII.xDummy = xDummy;
		CII.center = center;
		CII.level = level;
		std::complex<double> temp = 0.0 + 0.0*I;
		for (size_t l = 0; l < NumberOfBases; l++) {
			//cout << "l: " << l << endl;
			CII.l = l;
			obtain_Quadrature(epsilon);
			//obtain_Quadrature_Uniform();
			temp += Q_pinv(l,j)*integral;
		}
		return temp*I/4.0;
	}
	/*
	dtype integrand(pts2D P) {
    return P.x;
	};

	void evaluateC() {
		xcenter = 0.0; //of Quadratire object
		ycenter = 0.0; //of Quadratire object
		Lx = 1.0; //of Quadratire object
		Ly = 1.0; //of Quadratire object
		for (size_t i = 0; i < 5; i++) {
			obtain_Quadrature(epsilon);//when an integrand is evaluated gauss Quadrature tree is constructed fresh
			//obtain_Quadrature_Uniform();//when an integrand is evaluated gauss Quadrature tree is constructed fresh
			cout << integral << endl;
		}
		exit(0);
	}
	*/

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

  void evaluateC() {
		double min;
		double max;
		for (size_t j = 0; j <= NLEVELS; j++) {
			xcenter = 0.0; //of Quadratire object
      ycenter = 0.0; //of Quadratire object
      Lx = 1.0; //of Quadratire object
      Ly = 1.0; //of Quadratire object
      CII.level = j; //of Quadratire object
      C.push_back(Mat::Zero(2*seriesLength+1, NumberOfBases));
			for (int n = -seriesLength; n <= seriesLength; n++) {
        for (size_t l = 0; l < NumberOfBases; l++) {
          CII.n = n;
          CII.l = l;//for each n and l combination the integrand is evaluated;
          obtain_Quadrature(epsilon);//when an integrand is evaluated gauss Quadrature tree is constructed fresh
					//obtain_Quadrature_Uniform();//when an integrand is evaluated gauss Quadrature tree is constructed fresh
					C[j](n+seriesLength,l) = boxRadius[j]*boxRadius[j]*integral;
        }
      }
    }
  }


  void evaluateQ(Eigen::MatrixXd &Q_pinv_buf) {
    Q = Eigen::MatrixXd::Zero(rank, NumberOfBases);
		for (size_t i = 0; i < rank; i++) {
      for (size_t l = 0; l < NumberOfBases; l++) {
        orderedPair op = degreeOfBasesVec[l];
        int a = op.x;
        int b = op.y;
        Q(i,l) = get_ChebPoly(standardChebNodes[i].x, a)*get_ChebPoly(standardChebNodes[i].y, b);
      }
    }
		Q_pinv = Q.completeOrthogonalDecomposition().pseudoInverse();
		Q_pinv_buf = Q_pinv;
  }

  void evaluateM(std::vector<Mat > &M) { //Matrices C and M are invarinat among boxes of same level
		Eigen::MatrixXd Q_pinv_buf;
		evaluateQ(Q_pinv_buf);
		//cout << "Q_pinv_buf: " << endl << Q_pinv_buf << endl;
		evaluateC();
		for (size_t j = 0; j <= NLEVELS; j++) {
			M.push_back(C[j]*Q_pinv);
    }
  }

};
