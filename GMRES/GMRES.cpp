//#include "FMM2DTree_v2.hpp"
#include "G_FMM2DTree.hpp"

int Qchoice;
double gaussian(const pts2D r) {
	double R2 = r.x*r.x + r.y*r.y;
	return 1.5*exp(-160.0*R2);
}

double multipleGaussians(const pts2D r) {
	int noOfGaussians = 20;
	double centersX[20] = {-1.17182014706864, //randomly but well-separated gaussians
-0.861066625441127,
-1.14533939124057,
-0.856639124568525,
-0.620066350133170,
-0.457064031294315,
-0.482782006538245,
-0.635070123186999,
0.108594058866164,
-0.0315722590552472,
-0.177274727610547,
0.0342987980936911,
0.369662097773613,
0.591444264226524,
0.401225715100565,
0.666822699314234,
1.10327781188617,
0.827519268441487,
0.873411278211136,
1.09482906183464};
double centersY[20] = {-0.361066625441127,
0.0719520938231629,
0.354660608759432,
1.14336087543147,
-1.12006635013317,
-0.457064031294315,
0.364929876813001,
0.802261227250738,
-0.394084977497864,
-0.0315722590552472,
0.322725272389453,
1.03429879809369,
-1.13033790222639,
0.0137164199764117,
0.401225715100565,
1.16682269931423,
-0.345187707108795,
-0.172480731558513,
0.373411278211136,
1.09482906183464};
	double q = 0.0;
	double a = 0.0013;
	for (size_t i = 0; i < noOfGaussians; i++) {
		double R2 = (r.x-centersX[i])*(r.x-centersX[i]) + (r.y-centersY[i])*(r.y-centersY[i]);
		q += 1.5*exp(-R2/a);
	}
	return q;
}

double flatBump(const pts2D r) {
	double R2 = r.x*r.x + r.y*r.y;
	return 0.5*erfc(5.0*(R2-1.0));
}

double plasma(const pts2D r) {
	double C = 0.4987;
	double psi = 1.0-(r.x-0.15*(1-r.x*r.x))*(r.x-0.15*(1-r.x*r.x)) - C*(1.0+0.3*r.x)*(1.0+0.3*r.x)*r.y*r.y;
	double a[5] = {0.45, 0.195, 0.51, 0.195, 0.63};
	double x[5] = {0.8, 0.54, -0.14, -0.5, 0.18};
	double y[5] = {0.0, -0.28, 0.7, -0.01, 0.8};
	if (psi <= 0.05) {
		return 0.0;
	}
	else {
		double g = 0.0;
		for (size_t i = 0; i < 5; i++) {
			double temp = (r.x-x[i])*(r.x-x[i])+(r.y-y[i])*(r.y-y[i]);
			g += a[i]*exp(-temp/0.01);
		}
		return -1.5*(psi-0.05)-g*cos(0.9*r.y);
	}
}

double cavity(const pts2D r) {
	double x1 = r.x;
	double x2 = r.y;
	double R = sqrt(x1*x1 + x2*x2);
	double theta = atan2(x2, x1);
	return (1.0-pow(sin(0.5*theta),500))*exp(-2000*(0.1-R*R)*(0.1-R*R));
}

double lens(const pts2D r) {
	double x1 = r.x;
	double x2 = r.y;
	double R = sqrt(x1*x1 + x2*x2);
	return 4*(x2-0.1)*(1.0-erf(25.0*(R-0.3)));
}

class userkernel{
public:
	userkernel() {
	};
	double ContrastFunction(const pts2D r) {
		if (Qchoice == 0)
			return gaussian(r);
		else if (Qchoice == 1)
			return multipleGaussians(r);
		else if (Qchoice == 2)
			return plasma(r);
		else if (Qchoice == 3)
			return flatBump(r);
		else if (Qchoice == 4)
			return cavity(r);
		else
			return lens(r);
	};
	std::complex<double> IncidenceFunction(const pts2D r) {
		std::complex<double> q = exp(I*kappa*r.x);
		return q;
	};
	std::complex<double> RHSFunction(const pts2D r) {
		std::complex<double> q = -1.0*kappa*kappa*ContrastFunction(r)*IncidenceFunction(r);
		return q;
	};
	~userkernel() {};
};

struct inputsToDFMM {
  int nCones_LFR;
	int nChebNodes;
	double L;
	int yes2DFMM;
	int TOL_POW;
	int degreeOfBases;
	int treeAdaptivity;
	int nLevelsUniform;
};

class DFMM {
public:
  inputsToDFMM inputs;
  std::vector<pts2D> particles_X;//locations
	std::vector<pts2D> particles_Y;//dummy values
	std::vector<pts2D> gridPoints;//same as particles_X, particles_Y
	std::vector<int> row_indices;
	std::vector<int> col_indices;
	userkernel* mykernel;
	FMM2DTree<userkernel>* K;
  std::string filename;

  DFMM(inputsToDFMM inputs) {
    this->inputs = inputs;
    filename = "DFMMtree.tex";
    mykernel =	new userkernel();
    K	=	new FMM2DTree<userkernel>(mykernel, inputs.nCones_LFR, inputs.nChebNodes, inputs.L, inputs.yes2DFMM, inputs.TOL_POW, particles_X, particles_Y, row_indices, col_indices, inputs.degreeOfBases, inputs.treeAdaptivity, inputs.nLevelsUniform);
  };

  ~DFMM() {};

  void FMMTreeInitialisation() {
    K->set_Standard_Cheb_Nodes();
  	K->get_Transfer_Matrix();
  	// K->createUniformTree();
		K->createAdaptiveTree();
		K->makeLevelRestriction_Tree();
		string filename = "tree.tex";
		K->outputAdaptiveGrid(filename);
  	K->createCones();
  	K->assign_Tree_Interactions();
		K->assignLeafChargeLocations(gridPoints);//chargeLocations
    K->assignNonLeafChargeLocations();
    K->assignLeafChebNodes();
		// std::string filename3 = "M/M";

		currentDirectory = std::filesystem::current_path();
		char final [512];
		sprintf (final, "%s/precomputations",currentDirectory.c_str());
		mkdir(final, 0775);

		K->evaluatePrecomputations(); //precomputations
		K->writeMToFile(); //precomputations
		K->getNeighborInteractions(); //precomputations


		K->readMFromFile();
    K->getNeighbors();
    K->readNeighborInteractions();

		K->getNodes_LFR();
		K->getNodes_HFR();
		K->getUserCheckPoints();
		K->Assemble_HFR_M2M();
		K->Assemble_HFR_M2L();
		K->Assemble_LFR_M2L();
		K->Assemble_LFR_L2L();
		K->Assemble_HFR_L2L();
		K->collectBoundary();//to plot
  }

  void MatVecProduct(Vec &b, Vec &output) {
		//double start, end;
		//start		=	omp_get_wtime();
		/*if (!K->findPhi) {
			K->getNodes_LFR();
			K->getNodes_HFR();
			K->getUserCheckPoints();
			K->Assemble_HFR_M2M();
			K->Assemble_HFR_M2L();
			K->Assemble_LFR_M2L();
			K->Assemble_LFR_L2L();
			K->Assemble_HFR_L2L();
		}*/
		K->assignLeafCharges(b);
    K->assign_NonLeaf_Charges();
		// cout << "assign_NonLeaf_Charges done" << endl;
		K->LFR_M2M();
		// cout << "LFR_M2M done" << endl;
		K->HFR_M2M();
		// cout << "HFR_M2M done" << endl;
		K->LFR_M2L();
		// cout << "LFR_M2L done" << endl;
		K->HFR_M2L();
		// K->collectPotential(output);
		// cout << "M2L: " << output.norm() << endl;
		// cout << "HFR_M2L done" << endl;
		// cout << "collectPotential done" << endl;
		K->HFR_L2L();
		// cout << "HFR_L2L done" << endl;
		K->LFR_L2L();
		// K->collectPotential(output);
		// cout << "L2L: " << output.norm() << endl;
		//K->evaluate_NearField_Using_Precomputations_H();
		// cout << "LFR_L2L done" << endl;
		// cout << "collectPotential done" << endl;
		K->evaluate_NearField_Using_Precomputations_LS();
		// cout << "evaluate_NearField_Using_Precomputations_LS done" << endl;
		K->collectPotential(output);
		// cout << "NF: " << output.norm() << endl;
		//K->evaluate_NearField();//without precomputations
		// cout << "collectPotential done" << endl;
		//MatVecProduct error
		/*
		int t	=	rand()%K->leafNodes.size();
		double errC = K->perform_Error_Check2(t);
		int j = K->leafNodes[t].x;
		int k = K->leafNodes[t].y;
		std::cout << "j: " << j << "	boxNumber: " << K->tree[j][k].boxNumber << "		Error is: " << errC << endl;// << "	parILA: " << C->tree[nLevels-1][nBox/4].ILActive << std::endl;
		*/
		//end		=	omp_get_wtime();
	  //double timeMatVecProduct = end-start;
		//cout << "timeMatVecProduct: " << timeMatVecProduct << endl;

  }
};

template <typename SolverType>
void GMRES(Vec &x0, SolverType *S, int m, inputsToDFMM inputs, Vec &b, Vec &solution) {
  /*
  N: size of matrix; considering only square matrices
  m: dimension of search space
  b: RHS vector
  The matrix is implicitly defined through the solver
  */
  int N = b.size();
  Vec output;
	S->MatVecProduct(x0, output);
  Vec r0 = b - output;
  double beta = r0.norm();
	Mat V = Mat::Zero(N,m+1);
  Mat W = Mat::Zero(N,m);
  Mat H = Mat::Zero(m+1,m);
  V.col(0) = r0/beta;
  for (size_t j = 0; j < m; j++) {
		Vec lhs = V.col(j);
    Vec rhs;
		S->MatVecProduct(lhs, rhs);
    W.col(j) = rhs;
		for (size_t i = 0; i <= j; i++) {
      H(i,j) = V.col(i).dot(W.col(j));
      W.col(j) = W.col(j) - H(i,j)*V.col(i);
    }
    H(j+1,j) = W.col(j).norm();
		if (fabs(H(j+1,j)) < 1e-8) {
      m = j+1;
      cout << "early exit at m: " << m << endl;
      break;
    }
    V.col(j+1) = W.col(j)/H(j+1,j);
  }
  Vec beta_e1 = Vec::Zero(m+1);
  beta_e1(0) = beta;
  Mat updatedH = H.block(0,0,m+1,m);
  Mat updatedV = V.block(0,0,N,m);
	// cout << "updatedH.n(): " << updatedH.norm() << endl;
	// cout << "updatedH: " << endl << updatedH << endl;
	// cout << "updatedV.n(): " << updatedV.norm() << endl;

  Vec ym = updatedH.colPivHouseholderQr().solve(beta_e1);
	// cout << "ym.n(): " << ym.norm() << endl;
	// cout << "ym: " << endl << ym << endl;

	Vec update = updatedV*ym;
	solution = x0 + update;
}

template <typename SolverType>
void RestartedGMRES(int restart, SolverType *S, int m, inputsToDFMM inputs, Vec &b, Vec &solution) {
	int N = b.size();
  Vec x0 = Vec::Zero(N);
  for (size_t r = 0; r < restart; r++) {
		GMRES<SolverType>(x0, S, m, inputs, b, solution);
    x0 = solution;
		// Vec output;
    // S->MatVecProduct(solution, output);
    // Vec r0 = b - output;
    // cout << "Err in GMRES after restart: " << r0.norm()/b.norm() << endl << endl;
  }
}

int main(int argc, char* argv[]) {
	int nCones_LFR			=	atoi(argv[1]);
	int nChebNodes			=	atoi(argv[2]);
	int treeAdaptivity	=	atoi(argv[3]);
	double L						=	atof(argv[4]);
	kappa 							= atof(argv[5]);
	int yes2DFMM				=	atoi(argv[6]);
	int degreeOfBases 	= atoi(argv[7]);
	int TOL_POW					= atoi(argv[8]);
	int m 							= atoi(argv[9]);
	int restart					= atoi(argv[10]);
	int nLevelsUniform  = atoi(argv[11]);
	Qchoice		    			=	atoi(argv[12]);

  cout << "Wavenumber:		" << kappa << endl;
	cout << "Wavelength:		" << 2*PI/kappa << endl;
	cout << "no. of full cycles:	" << L*kappa/2*PI << endl;
  inputsToDFMM inputs;
  inputs.nCones_LFR = nCones_LFR;
  inputs.nChebNodes = nChebNodes;
  inputs.L = L;
  inputs.yes2DFMM = yes2DFMM;
  inputs.TOL_POW = TOL_POW;
	inputs.degreeOfBases = degreeOfBases;
  inputs.treeAdaptivity = treeAdaptivity;
	inputs.nLevelsUniform = nLevelsUniform;
  Vec Phi;
	double timeIn_getMatrixEntry_Offset = 0.0;
	double timeIn_getMatrixEntry;
	double start, end;

	start		=	omp_get_wtime();
  DFMM *S = new DFMM(inputs);

	S->FMMTreeInitialisation();
	end		=	omp_get_wtime();
	double timeAssemble =	(end-start);
	std::cout << std::endl << "Time taken to assemble: " << timeAssemble << std::endl;

  int N = S->K->N;
  Vec b(N);// = Vec::Random(N);//incidence field
	for (size_t i = 0; i < N; i++) {
		b(i) = S->mykernel->RHSFunction(S->gridPoints[i]); //exp(I*kappa*S->gridPoints[i].x);
	}
	// cout << "b.n(): " << b.lpNorm<Infinity>() << endl;

	//mkdir
	string result_filename = "result";
	string currPath = std::filesystem::current_path();
	char final[256];
	sprintf (final, "%s/%s", currPath.c_str(), result_filename.c_str());
	mkdir(final, 0775);

	result_filename = "result/result_" + std::to_string(Qchoice) + "_" + std::to_string(nChebNodes) + "_" + std::to_string(treeAdaptivity) + "_" + std::to_string(int(kappa)) + "_" + std::to_string(m);
	currPath = std::filesystem::current_path();
	sprintf (final, "%s/%s", currPath.c_str(), result_filename.c_str());
	mkdir(final, 0775);

	Vec x0 = Vec::Zero(N);
	Vec APhi;
	Vec r;
	//GMRES<DFMM>(x0, S, m, inputs, b, Phi);
	start		=	omp_get_wtime();
  RestartedGMRES<DFMM>(restart, S, 10, inputs, b, Phi);
	end		=	omp_get_wtime();
	double timeGMRES =	(end-start);

	string filename0;
	filename0 = result_filename + "/Phi_1_10";
	std::ofstream myfile0;
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 1,10: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;
	// exit(0);

	start		=	omp_get_wtime();
  RestartedGMRES<DFMM>(restart, S, 50, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_1_50";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 1,50: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;

	start		=	omp_get_wtime();
	RestartedGMRES<DFMM>(restart, S, 100, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_1_100";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 1,100: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;

	start		=	omp_get_wtime();
  RestartedGMRES<DFMM>(restart, S, 200, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_1_200";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 1,200: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;

	start		=	omp_get_wtime();
  RestartedGMRES<DFMM>(restart, S, 500, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_1_500";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 1,500: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;

	start		=	omp_get_wtime();
	RestartedGMRES<DFMM>(2, S, 10, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_2_10";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 2,10: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;

	start		=	omp_get_wtime();
	RestartedGMRES<DFMM>(2, S, 50, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_2_50";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 2,50: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;

	start		=	omp_get_wtime();
	RestartedGMRES<DFMM>(2, S, 100, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_2_100";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 2,100: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;

	start		=	omp_get_wtime();
	RestartedGMRES<DFMM>(2, S, 200, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_2_200";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 2,200: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;

	start		=	omp_get_wtime();
	RestartedGMRES<DFMM>(2, S, 500, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_2_50";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << std::endl << "Time taken by GMRES Solver 2,500: " << timeGMRES << std::endl;
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl << endl;


	string filenameB;
	filenameB = result_filename + "/B";


	//true matrix

	/*
	Mat B(N, N);
	#pragma omp parallel for
	for (int j=0; j < N; ++j) {
		#pragma omp parallel for
		for (int k=0; k < N; ++k) {
				B(j,k) = S->K->getMatrixEntry2(j, k);
		}
	}

	std::ofstream myfileB;
	myfileB.open(filenameB.c_str());
	myfileB << B << endl;
	myfileB.close();
	cout << "B formed" << endl;


	Mat BFromFile;
	S->K->ReadFromTextFile(BFromFile, filenameB);
	Vec APhi;
	APhi = BFromFile*Phi;
  Vec r0 = b - APhi;
  cout << "Err in GMRES Phi: " << r0.norm()/b.norm() << endl;


	string filenamePhi = "../../HODLR/build/result/Phi";
	Mat HODLR_Phi_Mat;
	Vec A_H_Phi, HODLR_Phi;
	S->K->ReadFromTextFile(HODLR_Phi_Mat, filenamePhi);
	HODLR_Phi = HODLR_Phi_Mat.col(0);
	A_H_Phi = B*HODLR_Phi;
	//S->MatVecProduct(HODLR_Phi, A_H_Phi);
	double errInHODLRSol = (A_H_Phi - b).norm()/b.norm();
	std::cout << "Err in HODLR Phi : " << errInHODLRSol << std::endl;
	*/

	//////////////////// FIND SCATTERED FIELD FROM PHI ///////////////////////////////////
	Vec U_incidence(N);
	DFMM *S2 = new DFMM(inputs);
	S2->K->findPhi = false;
	S2->FMMTreeInitialisation();
	for (size_t i = 0; i < N; i++) {
		U_incidence(i) = S2->mykernel->IncidenceFunction(S2->gridPoints[i]); //exp(I*kappa*S->gridPoints[i].x);
	}

	Vec U_scattered;
	start		=	omp_get_wtime();
  S2->MatVecProduct(Phi, U_scattered);
	end		=	omp_get_wtime();
	double timeMatVec =	(end-start);
	std::cout << std::endl << "Time taken for MatVec product: " << timeMatVec << std::endl;
	std::cout << std::endl << "Time taken for field computation: " << timeGMRES+timeMatVec << std::endl;

	Vec U_total = U_incidence + U_scattered;
	////////////// write result to file ///////////////////
	string filename;
	filename = result_filename + "/gridPointsX";
	std::ofstream myfile1;
	myfile1.open(filename.c_str());
	for (size_t l = 0; l < S->gridPoints.size(); l++) {
		myfile1 << S->gridPoints[l].x << endl;
	}
	myfile1.close();

	filename = result_filename + "/gridPointsY";
	std::ofstream myfile2;
	myfile2.open(filename.c_str());
	for (size_t l = 0; l < S->gridPoints.size(); l++) {
		myfile2 << S->gridPoints[l].y << endl;
	}
	myfile2	.close();

	filename = result_filename + "/solutionR";
	std::ofstream myfile3;
	myfile3.open(filename.c_str());
	for (size_t l = 0; l < U_total.size(); l++) {
		myfile3 << U_total(l).real() << endl;
	}
	myfile3.close();

	filename = result_filename + "/solutionI";
	std::ofstream myfile4;
	myfile4.open(filename.c_str());
	for (size_t l = 0; l < U_total.size(); l++) {
		myfile4 << U_total(l).imag() << endl;
	}
	myfile4.close();

	filename = result_filename + "/leftRightBoundary"; //to plot
	std::ofstream myfile5;
	myfile5.open(filename.c_str());
	for (size_t l = 0; l < S2->K->leftRightBoundary.size(); l++) {
		myfile5 << S2->K->leftRightBoundary[l].x << "	" << S2->K->leftRightBoundary[l].y << endl;
	}
	myfile5.close();

	filename = result_filename + "/bottomTopBoundary"; //to plot
	std::ofstream myfile6;
	myfile6.open(filename.c_str());
	for (size_t l = 0; l < S2->K->bottomTopBoundary.size(); l++) {
		myfile6 << S2->K->bottomTopBoundary[l].x << "	" << S2->K->bottomTopBoundary[l].y << endl;
	}
	myfile6.close();

	delete S;
	delete S2;
	////////////////////////////////////////////
}
