#include "G_FMM2DTree.hpp"

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
		std::string filename3 = "M/M";
    // K->check();

		K->evaluatePrecomputations(); //precomputations
		K->writeMToFile(filename3); //precomputations
		K->getNeighborInteractions(); //precomputations


		K->readMFromFile(filename3);
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


  void FMMTreeInitialisationNew() {
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
		std::string filename3 = "M/M";
    // K->check();

		K->evaluatePrecomputations(); //precomputations
		K->writeMToFile(filename3); //precomputations
		K->getNeighborInteractions(); //precomputations


		K->readMFromFile(filename3);
    K->getNeighbors();
    K->readNeighborInteractions();
  }

  void assemble() {
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

  void MatVecProduct(Vec b, Vec &output) {
  // void MatVecProduct(Vec b, Vec &output, Vec &trueRHS, Mat &trueA) {
  // void MatVecProduct(Vec b, Vec &output) {
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
		K->LFR_M2M();
		K->HFR_M2M();
    // std::cout << K->tree[2][5].incoming_potential << std::endl << std::endl;
		K->LFR_M2L();
		K->HFR_M2L();
    // std::cout << K->tree[2][5].incoming_potential << std::endl << std::endl;
		// K->collectPotential(output);
    // std::cout << "output: " << output << std::endl;
		K->HFR_L2L();
		K->LFR_L2L();
		//K->evaluate_NearField_Using_Precomputations_H();
    // std::cout << K->tree[2][5].incoming_potential << std::endl << std::endl;

		// K->collectPotential(output);
    // std::cout << "output: " << output << std::endl;
		K->evaluate_NearField_Using_Precomputations_LS();
  	//K->evaluate_NearField();//without precomputations
		K->collectPotential(output);
    // std::cout << "output: " << output << std::endl;
    //MatVecProduct error
    ////////////////////////////////////////////
    // for (size_t i = 0; i < 10; i++) {
    // // for (size_t t = 0; t < K->leafNodes.size(); t++) {
    //   int t	=	rand()%K->leafNodes.size();
    //   double errC = K->perform_Error_Check2(t);
  	// 	int j = K->leafNodes[t].x;
  	// 	int k = K->leafNodes[t].y;
  	// 	std::cout << "j: " << j << "	boxNumber: " << K->tree[j][k].boxNumber << "		Error is: " << errC << endl;// << "	parILA: " << C->tree[nLevels-1][nBox/4].ILActive << std::endl;
    // }
    // Mat trueA;
    // Vec trueRHS;
    // std::cout << "findTotalError: " << K->findTotalError(trueRHS, trueA) << std::endl;
    // std::cout << "-----------------------------" << std::endl;
    ////////////////////////////////////////////
    // exit(0);
		//end		=	omp_get_wtime();
	  //double timeMatVecProduct = end-start;
		//cout << "timeMatVecProduct: " << timeMatVecProduct << endl;

  }
};

class DFMM2 {
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

  DFMM2(inputsToDFMM inputs) {
    this->inputs = inputs;
    filename = "DFMMtree.tex";
    mykernel =	new userkernel();
    K	=	new FMM2DTree<userkernel>(mykernel, inputs.nCones_LFR, inputs.nChebNodes, inputs.L, inputs.yes2DFMM, inputs.TOL_POW, particles_X, particles_Y, row_indices, col_indices, inputs.degreeOfBases, inputs.treeAdaptivity, inputs.nLevelsUniform);
  };

  ~DFMM2() {};

  void FMMTreeInitialisation() {
    K->set_Standard_Cheb_Nodes();
  	K->get_Transfer_Matrix();
  	K->createUniformTree();
		// K->createAdaptiveTree();
		// K->makeLevelRestriction_Tree();
		string filename = "tree.tex";
		K->outputAdaptiveGrid(filename);
  	K->createCones();
  	K->assign_Tree_Interactions();
		K->assignLeafChargeLocations(gridPoints);//chargeLocations
    K->assignNonLeafChargeLocations();
    K->assignLeafChebNodes();
		std::string filename3 = "M/M";
    // K->check();

		// K->evaluatePrecomputations(); //precomputations
		// K->writeMToFile(filename3); //precomputations
		// K->getNeighborInteractions(); //precomputations


		K->readMFromFile(filename3);
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

  void MatVecProduct(Vec b, Vec &output) {
  // void MatVecProduct(Vec b, Vec &output, Vec &trueRHS, Mat &trueA) {
  // void MatVecProduct(Vec b, Vec &output) {
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
		K->LFR_M2M();
		K->HFR_M2M();
    // std::cout << K->tree[2][5].incoming_potential << std::endl << std::endl;
		K->LFR_M2L();
		K->HFR_M2L();
    // std::cout << K->tree[2][5].incoming_potential << std::endl << std::endl;
		// K->collectPotential(output);
    // std::cout << "output: " << output << std::endl;
		K->HFR_L2L();
		K->LFR_L2L();
		//K->evaluate_NearField_Using_Precomputations_H();
    // std::cout << K->tree[2][5].incoming_potential << std::endl << std::endl;

		// K->collectPotential(output);
    // std::cout << "output: " << output << std::endl;
		K->evaluate_NearField_Using_Precomputations_LS();
  	//K->evaluate_NearField();//without precomputations
		K->collectPotential(output);
    // std::cout << "output: " << output << std::endl;
    //MatVecProduct error
    ////////////////////////////////////////////
    // for (size_t i = 0; i < 10; i++) {
    // // for (size_t t = 0; t < K->leafNodes.size(); t++) {
    //   int t	=	rand()%K->leafNodes.size();
    //   double errC = K->perform_Error_Check2(t);
  	// 	int j = K->leafNodes[t].x;
  	// 	int k = K->leafNodes[t].y;
  	// 	std::cout << "j: " << j << "	boxNumber: " << K->tree[j][k].boxNumber << "		Error is: " << errC << endl;// << "	parILA: " << C->tree[nLevels-1][nBox/4].ILActive << std::endl;
    // }
    // Mat trueA;
    // Vec trueRHS;
    // std::cout << "findTotalError: " << K->findTotalError(trueRHS, trueA) << std::endl;
    // std::cout << "-----------------------------" << std::endl;
    ////////////////////////////////////////////
    // exit(0);
		//end		=	omp_get_wtime();
	  //double timeMatVecProduct = end-start;
		//cout << "timeMatVecProduct: " << timeMatVecProduct << endl;

  }
};

template <typename SolverType>
void GMRES(Vec &x0, SolverType *S, HODLR *T, int m, inputsToDFMM inputs, Vec &b_old, Vec &solution) {
  /*
  N: size of matrix; considering only square matrices
  m: dimension of search space
  b: RHS vector
  The matrix is implicitly defined through the solver
  */
  int N = b_old.size();
  Vec b = T->solve(b_old);
  Vec temp, output;
	S->MatVecProduct(x0, temp);//temp = Ax0
	output = T->solve(temp);
  Vec r0 = b - output;
  double beta = r0.norm();
	Mat V = Mat::Zero(N,m+1);
  Mat W = Mat::Zero(N,m);
  Mat H = Mat::Zero(m+1,m);
  V.col(0) = r0/beta;
  for (size_t j = 0; j < m; j++) {
		Vec lhs = V.col(j);
    Vec rhs;
    S->MatVecProduct(lhs, temp);
		rhs = T->solve(temp);
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
  Vec ym = updatedH.colPivHouseholderQr().solve(beta_e1);
	Vec update = updatedV*ym;
	solution = x0 + update;
}

template <typename SolverType>
void RestartedGMRES(int restart, SolverType *S, HODLR *T, int m, inputsToDFMM inputs, Vec &b, Vec &solution) {
	int N = b.size();
  Vec x0 = Vec::Zero(N);
  for (size_t r = 0; r < restart; r++) {
		GMRES<SolverType>(x0, S, T, m, inputs, b, solution);
    x0 = solution;
		// Vec output, temp;
    // S->MatVecProduct(solution, temp);
		// output = T->solve(temp);
    // Vec r0 = b - output;
    // cout << "Err in GMRES after restart: " << r0.norm()/b.norm() << endl << endl;
  }
}

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
