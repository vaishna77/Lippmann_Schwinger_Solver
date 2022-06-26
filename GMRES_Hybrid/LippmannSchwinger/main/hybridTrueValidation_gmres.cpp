#include "HODLR_Matrix.hpp"
#include "HODLR.hpp"
#include "KDTree.hpp"
double timeEntry;
#include "common.hpp"

#include "LS_GMRES.cpp"
#include "LS_HODLR.cpp"
typedef std::complex<double> Dtype;
#include "gmres.hpp"
//	get_ChebPoly
double get_ChebPoly(double x, int n) {
	return cos(n*acos(x));
}

void writeMatrixToFile(Mat& matrix, std::string filename) {
	//create directory
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile << std::setprecision(16);
	myfile << matrix << endl;
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
	double pre_conditioner_tolerance  = atof(argv[11]); // rank or tolerance
	int preconditioner_target_rank = atof(argv[12]);
	Qchoice					    =	atoi(argv[13]); //global variable
	int nLevelsUniform  = atoi(argv[14]);

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
	double timeIn_getMatrixEntry_Offset = 0.0;
	double timeIn_getMatrixEntry;
	double start, end;
	timeEntry = 0.0;
	///////////////////////// GMRES initialisation /////////////////////////////////
	// cout << "Begining GMRES initialisation" << endl;
	start		=	omp_get_wtime();
  DFMM *S = new DFMM(inputs);
	S->FMMTreeInitialisationNew();
	// S->FMMTreeInitialisation();
	end		=	omp_get_wtime();
	double timeInitialise =	(end-start);
	std::cout << std::endl << "Time taken to initialise: " << timeInitialise << std::endl;
	start		=	omp_get_wtime();
	S->assemble();
	end		=	omp_get_wtime();
	double timeAssemble =	(end-start);
	std::cout << std::endl << "Time taken to assemble: " << timeAssemble << std::endl;
	std::cout << std::endl << "Time taken to assemble without matrixEntry: " << timeAssemble-timeEntry << std::endl;
  // defining RHS
  int N = S->K->N;
  Vec b(N);// = Vec::Random(N);//incidence field
	for (size_t i = 0; i < N; i++) {
		b(i) = S->mykernel->RHSFunction(S->gridPoints[i]); //exp(I*kappa*S->gridPoints[i].x);
	}
	// cout << "GMRES initialisation done" << endl;
	////////////////////// SINGLE DAFMM //////////////////////////////////////////////////
	// Vec output;
	// start		=	omp_get_wtime();
	// S->MatVecProduct(b,output);
	// end		=	omp_get_wtime();
	// double timeApplySingleDAFMM =	(end-start);
	// std::cout << std::endl << "Time taken to apply single DAFMM: " << timeApplySingleDAFMM << std::endl;
	// exit(0);
	////////////////////////////////////////////////////////////////////////
	string result_filename;
	std::ofstream myfile0;
	Vec b_tilde;
	Vec APhi;
	Vec r;
	string filename0;
	/////////////////////////////////////////////////////////////////////////
  // running GMRES
	classGMRES* G = new classGMRES();
	int maxIterations = 400;
	double gmres_tolerance = 1.0e-12;
	double errGMRES;
	int noOfIterations;
	Vec Validation_Hybrid_Phi, Validation_GMRES_Phi;
	/////////////// GMRES Solver //////////////////////////////
	cout << "begining GMRES solver..." << endl;
	std::vector<double> e_gmres;
	start		=	omp_get_wtime();
	G->gmres(S, b, maxIterations, gmres_tolerance, Validation_GMRES_Phi, errGMRES, noOfIterations, e_gmres);
	end		=	omp_get_wtime();
	cout << "done" << endl;
	double timeGMRES = end-start;
	cout << "Time taken by GMRES solver: " << timeGMRES << endl;
	cout << "e_GMRES" << endl;
	for (size_t i = 0; i < noOfIterations; i++) {
		cout << e_gmres[i] << " ";
	}
	cout << endl;
	cout << "errGMRES: " << errGMRES << "	noOfIterations: " << noOfIterations << endl;
	/////////////////////////////////////////////
	/////////////////////// HODLR pre-conditioner initialisation /////////////////////
	// std::vector<pts2D> particles_X;//locations
	// std::vector<pts2D> particles_Y;//dummy values
	// userkernel* mykernel		=	new userkernel();
	// H_FMM2DTree<userkernel>* F	=	new H_FMM2DTree<userkernel>(mykernel, nCones_LFR, nChebNodes, L, yes2DFMM, TOL_POW, particles_X, particles_Y, kappa, degreeOfBases, treeAdaptivity, nLevelsUniform);
	// int H_N, M, dim;
	// H_N = F->gridPoints.size();
	// M = F->rank;
	// dim = 2;
	// Kernel* K = new Kernel(H_N, F);
	// bool is_sym = false;
	// bool is_pd = false;
	//
	// // for (size_t preconditioner_target_rank = 5; preconditioner_target_rank <= 25; preconditioner_target_rank+=10) {
	// 	HODLR* T = new HODLR(H_N, M, pre_conditioner_tolerance, int(preconditioner_target_rank));
	//
	// 	T->assemble(K, "rookPivoting", is_sym, is_pd);
	// 	start		=	omp_get_wtime();
	// 	T->factorize();
	// 	end		=	omp_get_wtime();
	// 	double timePreCondFact = end-start;
	// 	cout << "Time taken to preconditioner factorization: " << timePreCondFact << endl;
	// 	cout << "HODLR pre-conditioner initialisation done" << endl;
	//
	// 	cout << "begining Hybrid solver..." << endl;
	// 	std::vector<double> e_hybrid;
	// 	start		=	omp_get_wtime();
	// 	G->gmres(S, T, b, maxIterations, gmres_tolerance, Validation_Hybrid_Phi, errGMRES, noOfIterations, e_hybrid);
	// 	end		=	omp_get_wtime();
	// 	cout << "done" << endl;
	// 	double timeHybrid = end-start;
	// 	cout << "Time taken by Hybrid solver: " << timeHybrid << endl;
	// 	cout << "e_Hybrid" << endl;
	// 	for (size_t i = 0; i < noOfIterations; i++) {
	// 		cout << e_hybrid[i] << " ";
	// 	}
	// 	cout << endl;
	// 	cout << "errHybrid: " << errGMRES << "	noOfIterations: " << noOfIterations << endl;
	//
	// 	delete T;
	///////////////////////////////////////////
// }
	// exit(0);
	//////////////////// FIND SCATTERED FIELD FROM PHI ///////////////////////////////////
	DFMM2 *S2 = new DFMM2(inputs);
	S2->K->findPhi = false;
	cout << "here" << endl;
	S2->FMMTreeInitialisation();
	cout << "FMMTreeInitialisation done" << endl;
	int N2 = S2->K->N;
	Vec U_incidence(N2);

	for (size_t i = 0; i < N2; i++) {
		U_incidence(i) = S2->mykernel->IncidenceFunction(S2->gridPoints[i]); //exp(I*kappa*S->gridPoints[i].x);
	}
	Vec contrastVec(N2);
	for (size_t i = 0; i < N2; i++) {
		contrastVec(i) = S2->mykernel->ContrastFunction(S2->gridPoints[i]); //exp(I*kappa*S->gridPoints[i].x);
	}
	//////////////////////////////////////////
	std::vector<orderedPair> degreeOfBasesVec;
	for (size_t d = 0; d < degreeOfBases; d++) {
		for (size_t a = 0; a < d+1; a++) {
			int b = d-a;
			orderedPair op;
			op.x = a;
			op.y = b;
			degreeOfBasesVec.push_back(op);
		}
	}
	Vec Psi_non_grid = Vec::Zero(N2);
	int NumberOfBases = int((S->K->degreeOfBases)*(S->K->degreeOfBases+1)/2);
	for (size_t i = 0; i < N2; i++) {// iterate over S2->gridPoints[i]
		int j,k,t;//leaf box that contains S2->gridPoints[i]
		pts2D point = S2->gridPoints[i];
		for (t=0; t<S->K->leafNodes.size(); ++t) {
			j = S->K->leafNodes[t].x;
			k = S->K->leafNodes[t].y;
			if((S->K->tree[j][k].center.x-S->K->boxRadius[j] <= point.x && point.x < S->K->tree[j][k].center.x+S->K->boxRadius[j]) && (S->K->tree[j][k].center.y-S->K->boxRadius[j] <= point.y && point.y < S->K->tree[j][k].center.y+S->K->boxRadius[j])) {
				break;
			}
		}
		// std::cout << "i: " << i << "	" << S->K->tree[j][k].center.x-S->K->boxRadius[j] << ", " << point.x << ", " << S->K->tree[j][k].center.x+S->K->boxRadius[j] << std::endl;
		// std::cout << "i: " << i << "	" << S->K->tree[j][k].center.y-S->K->boxRadius[j] << ", " << point.y << ", " << S->K->tree[j][k].center.y+S->K->boxRadius[j] << std::endl;
		for (size_t l = 0; l < NumberOfBases; l++) {
			orderedPair op = degreeOfBasesVec[l];
	    int a = op.x;
	    int b = op.y;
			pts2D P;
			P.x = (S2->gridPoints[i].x - S->K->tree[j][k].center.x)/S->K->boxRadius[j];
			P.y = (S2->gridPoints[i].y - S->K->tree[j][k].center.y)/S->K->boxRadius[j];
			// if(P.x > 1 || P.x < -1 ||P.y > 1 || P.y < -1) {
			// 	std::cout << "oops;	" << "P.x	" << P.x << "	P.y	" << P.y << std::endl;
			// }
			Vec leaf_phi = Validation_GMRES_Phi.segment(t*nChebNodes*nChebNodes, nChebNodes*nChebNodes);
			// Vec leaf_phi = Validation_Hybrid_Phi.segment(t*nChebNodes*nChebNodes, nChebNodes*nChebNodes);
			std::complex<double> cB = 0;
			for (size_t r = 0; r < nChebNodes*nChebNodes; r++) {
				cB += S->K->Q_pinv2(l,r)*leaf_phi(r);
			}
			Psi_non_grid(i) += cB * get_ChebPoly(P.x, a)*get_ChebPoly(P.y, b);
		}
	}
	//////////////////////////////////////////
	Vec U_scattered_GMRES, U_scattered_Hybrid;
	cout << "Psi_non_grid.norm: " << Psi_non_grid.norm() << endl;
	////////////////////// GMRES ////////////////////////////////
	cout << "calculating u_scat_gmres from phi..." << endl;
  S2->MatVecProduct(Psi_non_grid, U_scattered_GMRES);
	cout << "done" << endl;
	Vec u_total_GMRES = U_scattered_GMRES + U_incidence;
	Vec errU_scattered_GMRES = Psi_non_grid/kappa/kappa;
	for (size_t i = 0; i < N2; i++) {
		errU_scattered_GMRES(i) += contrastVec(i)*u_total_GMRES(i);
	}
	cout << "error in u_scat_gmres: " << errU_scattered_GMRES.norm() << endl;
	/////////////////////// HYBRID ///////////////////////////////
	// cout << "calculating u_scat_hybrid from phi..." << endl;
	// S2->MatVecProduct(Psi_non_grid, U_scattered_Hybrid);
	// cout << "done" << endl;
	// Vec u_total_Hybrid = U_scattered_Hybrid + U_incidence;
	// Vec errU_scattered_Hybrid = Psi_non_grid/kappa/kappa;
	// for (size_t i = 0; i < N2; i++) {
	// 	errU_scattered_Hybrid(i) += contrastVec(i)*u_total_Hybrid(i);
	// }
	// cout << "error in u_scat_hybrid: " << errU_scattered_Hybrid.norm() << endl;
	//////////////////////////////////////////////////////
	//mkdir
	result_filename = "result";
	string currPath = std::experimental::filesystem::current_path();
	char final[256];
	sprintf (final, "%s/%s", currPath.c_str(), result_filename.c_str());
	mkdir(final, 0775);

	result_filename = "result/result_" + std::to_string(Qchoice) + "_" + std::to_string(nChebNodes) + "_" + std::to_string(treeAdaptivity) + "_" + std::to_string(int(kappa)) + "_" + std::to_string(m);
	currPath = std::experimental::filesystem::current_path();
	sprintf (final, "%s/%s", currPath.c_str(), result_filename.c_str());
	mkdir(final, 0775);
	string filename;
	std::ofstream myfile1;

	filename = result_filename + "/gridPointsX";
	myfile1.open(filename.c_str());
	for (size_t l = 0; l < N2; l++) {
		myfile1 << S2->gridPoints[l].x << endl;
	}
	myfile1.close();

	filename = result_filename + "/gridPointsY";
	myfile1.open(filename.c_str());
	for (size_t l = 0; l < N2; l++) {
		myfile1 << S2->gridPoints[l].y << endl;
	}
	myfile1.close();

	filename = result_filename + "/solutionR";
	myfile1.open(filename.c_str());
	for (size_t l = 0; l < N2; l++) {
		myfile1 << u_total_GMRES(l).real() << endl;
	}
	myfile1.close();

	filename = result_filename + "/solutionI";
	myfile1.open(filename.c_str());
	for (size_t l = 0; l < N2; l++) {
		myfile1 << u_total_GMRES(l).imag() << endl;
	}
	myfile1.close();

	filename = result_filename + "/errU_scattered_GMRES_real";
	myfile1.open(filename.c_str());
	for (size_t l = 0; l < N2; l++) {
		myfile1 << errU_scattered_GMRES(l).real() << endl;
	}
	myfile1.close();

	filename = result_filename + "/errU_scattered_GMRES_imag";
	myfile1.open(filename.c_str());
	for (size_t l = 0; l < N2; l++) {
		myfile1 << errU_scattered_GMRES(l).imag() << endl;
	}
	myfile1.close();

	// filename = result_filename + "/solutionR";
	// myfile1.open(filename.c_str());
	// for (size_t l = 0; l < N2; l++) {
	// 	myfile1 << u_total_Hybrid(l).real() << endl;
	// }
	// myfile1.close();
	//
	// filename = result_filename + "/solutionI";
	// myfile1.open(filename.c_str());
	// for (size_t l = 0; l < N2; l++) {
	// 	myfile1 << u_total_Hybrid(l).imag() << endl;
	// }
	// myfile1.close();

	// filename = result_filename + "/errU_scattered_Hybrid_real";
	// myfile1.open(filename.c_str());
	// for (size_t l = 0; l < N2; l++) {
	// 	myfile1 << errU_scattered_Hybrid(l).real() << endl;
	// }
	// myfile1.close();
	//
	// filename = result_filename + "/errU_scattered_Hybrid_imag";
	// myfile1.open(filename.c_str());
	// for (size_t l = 0; l < N2; l++) {
	// 	myfile1 << errU_scattered_Hybrid(l).imag() << endl;
	// }
	// myfile1.close();

	delete S2;
	delete S;
	delete G;
}
