#include "HODLR_Matrix.hpp"
#include "HODLR.hpp"
#include "KDTree.hpp"
#include "common.hpp"
#include "LS_GMRES.cpp"
#include "LS_HODLR.cpp"

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
	int nLevelsUniform  = 4;
	double pre_conditioner_tolerance  = atof(argv[11]); // rank or tolerance
	int preconditioner_target_rank = atof(argv[12]);
	Qchoice					    =	atoi(argv[13]); //global variable

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

	///////////////////////// GMRES initialisation /////////////////////////////////
	cout << "Begining GMRES initialisation" << endl;
	start		=	omp_get_wtime();
  DFMM *S = new DFMM(inputs);
	S->FMMTreeInitialisation();
	end		=	omp_get_wtime();
	double timeAssemble =	(end-start);
	std::cout << std::endl << "Time taken to assemble: " << timeAssemble << std::endl;
  // defining RHS
  int N = S->K->N;
  Vec b(N);// = Vec::Random(N);//incidence field
	for (size_t i = 0; i < N; i++) {
		b(i) = S->mykernel->RHSFunction(S->gridPoints[i]); //exp(I*kappa*S->gridPoints[i].x);
	}
	cout << "GMRES initialisation done" << endl;
	/////////////////////////////////////////////////////////////////////////



	/////////////////////// HODLR pre-conditioner initialisation /////////////////////
  // double pre_conditioner_tolerance  	= pow(10, -4); //preconditioner tolerance
	std::vector<pts2D> particles_X;//locations
  std::vector<pts2D> particles_Y;//dummy values
	userkernel* mykernel		=	new userkernel();
  H_FMM2DTree<userkernel>* F	=	new H_FMM2DTree<userkernel>(mykernel, nCones_LFR, nChebNodes, L, yes2DFMM, TOL_POW, particles_X, particles_Y, kappa, degreeOfBases, treeAdaptivity, nLevelsUniform);
  int H_N, M, dim;
  H_N = F->gridPoints.size();
  M = F->rank;
  dim = 2;
	Kernel* K = new Kernel(H_N, F);
  bool is_sym = false;
  bool is_pd = false;

	string result_filename;
	std::ofstream myfile0;
	Vec b_tilde;
	double timeGMRES;
	Vec APhi;
	Vec r;
	string filename0;
	for (double pre_tol = 1e-2; pre_tol <= 1e-1; pre_tol*=1e+1)
{
	int pre_rank = 5;
	cout << endl << "Begining HODLR pre-conditioner initialisation; pre_tol: " << pre_tol << endl;
	HODLR* T = new HODLR(H_N, M, pre_tol, pre_rank);
	// HODLR* T = new HODLR(H_N, M, pre_conditioner_tolerance, pre_rank);
  T->assemble(K, "rookPivoting", is_sym, is_pd);
	start		=	omp_get_wtime();
  T->factorize();
	end		=	omp_get_wtime();
	double timePreCondFact = end-start;
	cout << "Time taken to preconditioner factorization: " << timePreCondFact << endl;
	cout << "HODLR pre-conditioner initialisation done" << endl;
	/////////////////////////////////////////////////////////////////////////
	//mkdir
	result_filename = "result";
	string currPath = std::filesystem::current_path();
	char final[256];
	sprintf (final, "%s/%s", currPath.c_str(), result_filename.c_str());
	mkdir(final, 0775);

	result_filename = "result/result_" + std::to_string(Qchoice) + "_" + std::to_string(nChebNodes) + "_" + std::to_string(treeAdaptivity) + "_" + std::to_string(int(kappa)) + "_" + std::to_string(m);
	currPath = std::filesystem::current_path();
	sprintf (final, "%s/%s", currPath.c_str(), result_filename.c_str());
	mkdir(final, 0775);

  // running GMRES
	start		=	omp_get_wtime();
  b_tilde = T->solve(b);
  RestartedGMRES<DFMM>(restart, S, T, 10, inputs, b_tilde, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	filename0 = result_filename + "/Phi_1_10_P";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << "Time taken by GMRES Solver 1,10: " << timeGMRES << std::endl << endl;
	if (r.norm()/b.norm() < 1e-10) {
		continue;
	}
	//
	start		=	omp_get_wtime();
  b_tilde = T->solve(b);
  RestartedGMRES<DFMM>(2, S, T, 10, inputs, b_tilde, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	filename0 = result_filename + "/Phi_2_10_P";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << "Time taken by GMRES Solver 2,10: " << timeGMRES << std::endl << endl;
	if (r.norm()/b.norm() < 1e-10) {
		continue;
	}
	// ///////////////////////// GMRES without preconditioner ///////////////////////
	// /////////////////////////////////////////////////////////////
	//
	start		=	omp_get_wtime();
  b_tilde = T->solve(b);
  RestartedGMRES<DFMM>(restart, S, T, 50, inputs, b_tilde, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	filename0 = result_filename + "/Phi_1_50_P";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << "Time taken by GMRES Solver 1,50: " << timeGMRES << std::endl << endl;
	if (r.norm()/b.norm() < 1e-10) {
		continue;
	}
	// ///////////////////////// GMRES without preconditioner ///////////////////////
	// /////////////////////////////////////////////////////////////
	//
	// start		=	omp_get_wtime();
  // b_tilde = T->solve(b);
  // RestartedGMRES<DFMM>(2, S, T, 50, inputs, b_tilde, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	// filename0 = result_filename + "/Phi_2_50_P";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// std::cout << "Time taken by GMRES Solver 2,50: " << timeGMRES << std::endl << endl;
	// if (r.norm()/b.norm() < 1e-10) {
	// 	continue;
	// }
	///////////////////////// GMRES without preconditioner ///////////////////////
	/////////////////////////////////////////////////////////////


	start		=	omp_get_wtime();
  b_tilde = T->solve(b);
  RestartedGMRES<DFMM>(restart, S, T, 100, inputs, b_tilde, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	filename0 = result_filename + "/Phi_1_100_P";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	std::cout << "Time taken by GMRES Solver 1,100: " << timeGMRES << std::endl << endl;
	if (r.norm()/b.norm() < 1e-10) {
		continue;
	}


	// start		=	omp_get_wtime();
  // b_tilde = T->solve(b);
  // RestartedGMRES<DFMM>(2, S, T, 100, inputs, b_tilde, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	// filename0 = result_filename + "/Phi_2_100_P";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// std::cout << "Time taken by GMRES Solver 2,100: " << timeGMRES << std::endl << endl;
	// if (r.norm()/b.norm() < 1e-10) {
	// 	continue;
	// }

	// start		=	omp_get_wtime();
  // b_tilde = T->solve(b);
  // RestartedGMRES<DFMM>(restart, S, T, 200, inputs, b_tilde, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	// filename0 = result_filename + "/Phi_1_200_P";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// std::cout << "Time taken by GMRES Solver 1,200: " << timeGMRES << std::endl << endl;
	delete T;
}
// exit(0);
cout << "------------------------------------------------------" << endl;
	// start		=	omp_get_wtime();
	// b_tilde = T->solve(b);
	// RestartedGMRES<DFMM>(2, S, T, 100, inputs, b_tilde, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	// filename0 = result_filename + "/Phi_2_100_P";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// std::cout << "Time taken by GMRES Solver 2,100: " << timeGMRES << std::endl << endl;

	///////////////////////// GMRES without preconditioner ///////////////////////
	///////////////////////// GMRES without preconditioner ///////////////////////
	start		=	omp_get_wtime();
  RestartedGMRES<DFMM>(restart, S, 10, inputs, b, Phi);
	end		=	omp_get_wtime();
	timeGMRES =	(end-start);
	filename0 = result_filename + "/Phi_1_10";
	myfile0.open(filename0.c_str());
	for (size_t l = 0; l < Phi.size(); l++) {
		myfile0 << Phi(l) << endl;
	}
	myfile0.close();
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	std::cout << "Time taken by GMRES Solver without precond 1,10: " << timeGMRES << std::endl << endl;
	//
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
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	std::cout << "Time taken by GMRES Solver without precond 2,10: " << timeGMRES << std::endl << endl;
	//
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
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	std::cout << "Time taken by GMRES Solver without precond 1,50: " << timeGMRES << std::endl << endl;
	//
	// start		=	omp_get_wtime();
  // RestartedGMRES<DFMM>(2, S, 50, inputs, b, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// filename0 = result_filename + "/Phi_2_50";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	// std::cout << "Time taken by GMRES Solver without precond 2,50: " << timeGMRES << std::endl << endl;

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
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	std::cout << "Time taken by GMRES Solver without precond 1,100: " << timeGMRES << std::endl << endl;

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
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	std::cout << "Time taken by GMRES Solver without precond 2,100: " << timeGMRES << std::endl << endl;
	/////////////////////////////////////////////////////////////


	///////////////////////// GMRES without preconditioner ///////////////////////
	/////////////////////////////////////////////////////////////

	// start		=	omp_get_wtime();
  // b_tilde = T->solve(b);
  // RestartedGMRES<DFMM>(restart, S, T, 200, inputs, b_tilde, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	// filename0 = result_filename + "/Phi_1_200_P";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// std::cout << "Time taken by GMRES Solver 1,200: " << timeGMRES << std::endl << endl;

	///////////////////////// GMRES without preconditioner ///////////////////////
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
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	std::cout << "Time taken by GMRES Solver without precond 1,200: " << timeGMRES << std::endl << endl;

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
	S->MatVecProduct(Phi, APhi);
	r = b - APhi;
	cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	std::cout << "Time taken by GMRES Solver without precond 2,200: " << timeGMRES << std::endl << endl;
	//
	// start		=	omp_get_wtime();
  // RestartedGMRES<DFMM>(restart, S, 500, inputs, b, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// filename0 = result_filename + "/Phi_1_500";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES without precond: " << r.norm()/b.norm() << endl;
	// std::cout << "Time taken by GMRES Solver without precond 1,500: " << timeGMRES << std::endl << endl;
	/////////////////////////////////////////////////////////////

	// start		=	omp_get_wtime();
  // b_tilde = T->solve(b);
  // RestartedGMRES<DFMM>(2, S, T, 200, inputs, b_tilde, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	// filename0 = result_filename + "/Phi_2_200_P";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// std::cout << "Time taken by GMRES Solver 2,200: " << timeGMRES << std::endl << endl;

	///////////////////////// GMRES without preconditioner ///////////////////////
	/////////////////////////////////////////////////////////////

	// start		=	omp_get_wtime();
  // b_tilde = T->solve(b);
  // RestartedGMRES<DFMM>(restart, S, T, 500, inputs, b_tilde, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	// filename0 = result_filename + "/Phi_1_500_P";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// std::cout << "Time taken by GMRES Solver 1,500: " << timeGMRES << std::endl << endl;

	///////////////////////// GMRES without preconditioner ///////////////////////
	/////////////////////////////////////////////////////////////

	// start		=	omp_get_wtime();
  // b_tilde = T->solve(b);
  // RestartedGMRES<DFMM>(2, S, T, 500, inputs, b_tilde, Phi);
	// end		=	omp_get_wtime();
	// timeGMRES =	(end-start);
	// S->MatVecProduct(Phi, APhi);
	// r = b - APhi;
	// cout << "Err in GMRES: " << r.norm()/b.norm() << endl;
	// filename0 = result_filename + "/Phi_2_500_P";
	// myfile0.open(filename0.c_str());
	// for (size_t l = 0; l < Phi.size(); l++) {
	// 	myfile0 << Phi(l) << endl;
	// }
	// myfile0.close();
	// std::cout << "Time taken by GMRES Solver 2,500: " << timeGMRES << std::endl << endl;

	///////////////////////// GMRES without preconditioner ///////////////////////
	/////////////////////////////////////////////////////////////
	// exit(0);
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
