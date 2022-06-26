#ifndef _FMM2DTree_HPP__
#define _FMM2DTree_HPP__

class FMM2DCone {
public:
	bool outgoing_ILActive;
	bool incoming_ILActive;
	//Directional multipoles and locals of the box in this cone direction
	Vec multipoles;

	//defined only in HFR
	//std::vector<pts2D> &outgoing_chargePoints;//equivalent points {y_{k}^{B,o,l}}
	Vec outgoing_charges;//equivalent densities {f_{k}^{B,o,l}}
	//std::vector<pts2D> &outgoing_checkPoints;//check points {x_{k}^{B,o,l}}
	Vec outgoing_potential;//check potentials {u_{k}^{B,o,l}}

	//outgoing_chargePoints and outgoing_checkPoints are considered to be same as incoming_checkPoints and incoming_chargePoints respectively; hence are not declared
	//This is because of the assumption that the potentials are evaluated at the charge locations
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i,l}}
	Vec incoming_charges;//equivalent densities {f_{k}^{B,i,l}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i,l}}
	std::vector<int> outgoing_chargePoints;//equivalent points {y_{k}^{B,i,l}}
	std::vector<int> outgoing_checkPoints;//check points {x_{k}^{B,i,l}}
	std::vector<int> user_checkPoints;
	Vec incoming_potential;//check potentials {u_{k}^{B,i,l}}

	ColPivHouseholderQR<Mat> outgoing_Atilde_dec;
	ColPivHouseholderQR<Mat> incoming_Atilde_dec;

	std::vector<orderedPair > InteractionList;
	double angle;
	Mat outgoing_Ar;
	Mat L2L[4], M2M, Ktilde;
	std::vector<Mat> M2L;
	FMM2DCone () {}

};

class FMM2DBox {
public:
	double radius;
	int active;
	bool outgoing_ILActive;
	bool incoming_ILActive;
	int level;
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	bool isLeaf;

	std::vector<orderedPair > neighborNumbers;

	int fineNeighbors[12];//12
	int coarseNeighbors[12];//12
	int separatedFineNeighbors[20];//20
	int colleagueNeighbors[9];//9

	std::vector<orderedPair > InteractionList;
	std::vector<int > MFR_IL;
	//defined only in LFR
	//std::vector<pts2D> &outgoing_chargePoints;//equivalent points {y_{k}^{B,o}}
	Vec outgoing_charges;//equivalent densities {f_{k}^{B,o}}
	//std::vector<pts2D> &outgoing_checkPoints;//check points {x_{k}^{B,o}}
	Vec outgoing_potential;//check potentials {u_{k}^{B,o}}

	//outgoing_chargePoints and outgoing_checkPoints are considered to be same as incoming_checkPoints and incoming_chargePoints respectively; hence are not declared
	//This is because of the assumption that the potentials are evaluated at the charge locations
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	Vec incoming_charges;//equivalent densities {f_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	std::vector<int> user_checkPoints;
	std::vector<int> outgoing_chargePoints;//equivalent points {y_{k}^{B,i,l}}
	std::vector<int> outgoing_checkPoints;//check points {x_{k}^{B,i,l}}
	Vec incoming_potential;//check potentials {u_{k}^{B,i}}

	ColPivHouseholderQR<Mat> outgoing_Atilde_dec;
	ColPivHouseholderQR<Mat> incoming_Atilde_dec;

	Mat outgoing_Ar;
	Mat L2L, Ktilde;
	std::vector<Mat> M2L;

	FMM2DBox () {
		boxNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<4; ++l) {
			childrenNumbers[l]	=	-1;
		}
		isLeaf = false;
		active = true;
		for (size_t i = 0; i < 12; i++) {
			fineNeighbors[i] = -1;
		}
		for (size_t i = 0; i < 12; i++) {
			coarseNeighbors[i] = -1;
		}
		for (size_t i = 0; i < 20; i++) {
			separatedFineNeighbors[i] = -1;
		}
		for (size_t i = 0; i < 9; i++) {
			colleagueNeighbors[i] = -1;
		}
	}

	Vec multipoles;

	pts2D center;
	std::vector<int> chargeLocations;

	std::vector<pts2D> chebNodes;
	std::vector<FMM2DCone> ConeTree;
};

template <typename kerneltype>
class FMM2DTree: public G_LowRank {
public:
	double timeIn_getMatrixEntry;//for time profiling
	long NoIn_getMatrixEntry;
	kerneltype* Q;
	int nLevels;			//	Number of levels in the tree.
	int nChebNodes;			//	Number of Chebyshev nodes along one direction.
	int rank;				//	Rank of interaction, i.e., rank = nChebNodes*nChebNodes.
	int N;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).
	double a;				//	Cut-off for self-interaction. This is less than the length of the smallest box size.
	const double A = 3.0; //higher A, higher level_FarField
	const double B = 0.35;//4.0; //smaller B, higher level_LFR
	int level_LFR;
	int level_FarField;
	int nCones_LFR;
	double epsilonTree;//pow(10,-epsilonTree);
	Vec chargesAll;
	double SFN_angles[20];
	double N_angles[9];
	double FN_angles[12];
	double CN_angles[12];

	std::vector<pts2D> bottomTopBoundary;//used for plotting field
	std::vector<pts2D> leftRightBoundary;//used for plotting field
	bool findPhi;
	string NeighborFilename;
	string MFilename;
	std::vector<pts2D> gridPoints;
	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<double> ConeAperture;
	std::vector<int> nCones;
	std::vector<double> boxHomogRadius;			//	Stores the value of boxRadius^{alpha}
	std::vector<double> boxLogHomogRadius;		//	Stores the value of alpha*log(boxRadius)
	std::vector<std::vector<FMM2DBox> > tree;	//	The tree storing all the information.
	std::vector<std::vector<int> > indexTree;
	int TOL_POW;
	int treeAdaptivity;
	//	Different Operators
	int yes2DFMM;
	std::vector<double> standardChebNodes1D;
	std::vector<pts2D> standardChebNodes;
	std::vector<pts2D> standardChebNodesChild;
	std::vector<pts2D> leafChebNodes;
	int nLevelsUniform;
	std::vector<orderedPair> leafNodes;
	//	Different Operators
	Eigen::MatrixXd M2M[4];					//	Transfer from multipoles of 4 children to multipoles of parent.
	Eigen::MatrixXd L2L[4];					//	Transfer from locals of parent to locals of 4 children.
	FarFieldInteraction* FF; //object needed to computeFarField interaction of Lippmann Schwinger problem
	FarFieldInteraction2* FF2;
	std::vector<Mat > M; //M matrices computed during precomputations
	Eigen::MatrixXd Q_pinv;
	Eigen::MatrixXd Q_pinv2;
	int degreeOfBases;
	Mat colleagueNeighborInteraction[9][9];
	Mat fineNeighborInteraction[9][12];
	Mat separatedFineNeighborInteraction[9][20];
	Mat coarseNeighborInteraction[9][12];

// public:
	FMM2DTree(kerneltype* Q, int nCones_LFR, int nChebNodes, double L, int yes2DFMM, int TOL_POW, std::vector<pts2D> particles_X, std::vector<pts2D> particles_Y, std::vector<int> row_indices, std::vector<int> col_indices, int degreeOfBases, int treeAdaptivity, int nLevelsUniform):
	G_LowRank(TOL_POW, particles_X, particles_Y, row_indices, col_indices) {
		this->nLevelsUniform = nLevelsUniform;
		this->findPhi = true;
		this->NoIn_getMatrixEntry = 0;
		this->timeIn_getMatrixEntry = 0.0;
		this->degreeOfBases = degreeOfBases;//p minus 1
		this->treeAdaptivity = treeAdaptivity;
		this->epsilonTree = pow(10,-treeAdaptivity);
		this->FF = new FarFieldInteraction(nChebNodes, L, degreeOfBases);
		this->FF2 = new FarFieldInteraction2(nChebNodes, L, degreeOfBases);
		FF2->evaluateQ(Q_pinv2);
		this->Q = Q;
		this->TOL_POW = TOL_POW;
		this->nChebNodes		=	nChebNodes;
		this->rank			=	nChebNodes*nChebNodes;
		this->L				=	L;
		this->nCones_LFR		=	nCones_LFR;
		this->yes2DFMM = yes2DFMM;
		nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		this->a			=	smallestBoxSize;
		int k;
		if (yes2DFMM == 1) {
			if (kappa==0)
				this->level_LFR	=	2;	//actually should be 2; but for checking the accuracy of DFMM i am making it 3; so that HFR code runs even for LFR; if that gives good result it means DFMM code is perfect
			else {
				this->level_LFR	=	floor(log(kappa*L/B)/log(2.0));
			}
			if (level_LFR < 2) level_LFR = 2;
		}
		else {
			level_LFR = 2;
		}
		//cout << "level_LFR: " << level_LFR << endl;
	}

	void evaluatePrecomputations() {
		FF->evaluateQ(Q_pinv);
		FF->evaluateM(M);
	}

	void getIfNeighbor(int j, int k, int nj, int nk, int& typeOfNeighbor, int& NumOfNeighbor) {
		if ( fabs(tree[j][k].center.x - tree[nj][nk].center.x) >= 3.0*boxRadius[j] + boxRadius[nj]-epsilonRoundOff
		 || fabs(tree[j][k].center.y - tree[nj][nk].center.y) >= 3.0*boxRadius[j] + boxRadius[nj]-epsilonRoundOff) {
			 typeOfNeighbor = -1;
		}
		else {
			double angle = atan2(tree[nj][nk].center.y-tree[j][k].center.y, tree[nj][nk].center.x-tree[j][k].center.x);
			if (tree[j][k].radius == tree[nj][nk].radius) {
				typeOfNeighbor = 0;
				if (fabs(tree[nj][nk].center.y-tree[j][k].center.y)<epsilonRoundOff && fabs(tree[nj][nk].center.x-tree[j][k].center.x)<epsilonRoundOff) {
					NumOfNeighbor = 8;
				} // self
				else if (fabs(N_angles[0] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 0;
				}
				else if (fabs(N_angles[1] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 1;
				}
				else if (fabs(N_angles[2] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 2;
				}
				else if (fabs(N_angles[3] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 3;
				}
				else if (fabs(N_angles[4] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 4;
				}
				else if (fabs(N_angles[5] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 5;
				}
				else if (fabs(N_angles[6] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 6;
				}
				else if (fabs(N_angles[7] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 7;
				}
			}
			else if (tree[j][k].radius < tree[nj][nk].radius) { //coarse neighbors
				typeOfNeighbor = 1;
				if (fabs(CN_angles[0] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 0;
				}
				else if (fabs(CN_angles[1] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 1;
				}
				else if (fabs(CN_angles[2] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 2;
				}
				else if (fabs(CN_angles[3] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 3;
				}
				else if (fabs(CN_angles[4] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 4;
				}
				else if (fabs(CN_angles[5] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 5;
				}
				else if (fabs(CN_angles[6] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 6;
				}
				else if (fabs(CN_angles[7] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 7;
				}
				else if (fabs(CN_angles[8] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 8;
				}
				else if (fabs(CN_angles[9] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 9;
				}
				else if (fabs(CN_angles[10] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 10;
				}
				else if (fabs(CN_angles[11] - angle) < epsilonRoundOff) {
					NumOfNeighbor = 11;
				}
			}
			else { //fine neighbors and separated fine neighbors
				 if (find_distance(j, k, nj, nk) > 2.5*tree[j][k].radius) { //separated fine neighbors
					 typeOfNeighbor = 3;
					 if (fabs(SFN_angles[0] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 0;
					 }
					 else if (fabs(SFN_angles[1] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 1;
					 }
					 else if (fabs(SFN_angles[2] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 2;
					 }
					 else if (fabs(SFN_angles[3] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 3;
					 }
					 else if (fabs(SFN_angles[4] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 4;
					 }
					 else if (fabs(SFN_angles[5] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 5;
					 }
					 else if (fabs(SFN_angles[6] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 6;
					 }
					 else if (fabs(SFN_angles[7] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 7;
					 }
					 else if (fabs(SFN_angles[8] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 8;
					 }
					 else if (fabs(SFN_angles[9] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 9;
					 }
					 else if (fabs(SFN_angles[10] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 10;
					 }
					 else if (fabs(SFN_angles[11] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 11;
					 }
					 else if (fabs(SFN_angles[12] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 12;
					 }
					 else if (fabs(SFN_angles[13] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 13;
					 }
					 else if (fabs(SFN_angles[14] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 14;
					 }
					 else if (fabs(SFN_angles[15] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 15;
					 }
					 else if (fabs(SFN_angles[16] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 16;
					 }
					 else if (fabs(SFN_angles[17] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 17;
					 }
					 else if (fabs(SFN_angles[18] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 18;
					 }
					 else {
						 NumOfNeighbor = 19;
					 }
				 }
				 else {//fine neighbors
					 typeOfNeighbor = 2;
					 if (fabs(FN_angles[0] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 0;
					 }
					 else if (fabs(FN_angles[1] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 1;
					 }
					 else if (fabs(FN_angles[2] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 2;
					 }
					 else if (fabs(FN_angles[3] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 3;
					 }
					 else if (fabs(FN_angles[4] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 4;
					 }
					 else if (fabs(FN_angles[5] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 5;
					 }
					 else if (fabs(FN_angles[6] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 6;
					 }
					 else if (fabs(FN_angles[7] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 7;
					 }
					 else if (fabs(FN_angles[8] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 8;
					 }
					 else if (fabs(FN_angles[9] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 9;
					 }
					 else if (fabs(FN_angles[10] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 10;
					 }
					 else if (fabs(FN_angles[11] - angle) < epsilonRoundOff) {
						 NumOfNeighbor = 11;
					 }
				 }
			}
		}
	}

	kernel_dtype getMatrixEntry2(const unsigned i, const unsigned j) {

		/*
		pts2D ri = particles_X[i];//particles_X is a member of base class FMM_Matrix
		pts2D rj = particles_X[j];//particles_X is a member of base class FMM_Matrix
		double R2 = (ri.x-rj.x)*(ri.x-rj.x) + (ri.y-rj.y)*(ri.y-rj.y);
		double R = sqrt(R2);
		kernel_dtype out = exp(I*kappa*R)/R;

		double end		=	omp_get_wtime();
		timeIn_getMatrixEntry += (end-start);
		//cout << end-start << endl;
		if (NoIn_getMatrixEntry == 1) {
			cout << "avg time to get a matrix entry: " << timeIn_getMatrixEntry << endl;
		}
		if (R < 1e-8) {
			return R;
		}
		else {
			return out;
		}
		*/

		pts2D ri = particles_X[i];//particles_X is a member of base class FMM_Matrix
		unsigned int leaf_i = i/rank;//identify to which leaf jth charge belongs
		unsigned int index_i = i%rank;
		unsigned int level_i = leafNodes[leaf_i].x;//level to which the leaf box that contains jth charge belongs
		unsigned int box_i = leafNodes[leaf_i].y;

		unsigned int leaf_j = j/rank;//identify to which leaf jth charge belongs
		unsigned int index_j = j%rank;
		unsigned int level_j = leafNodes[leaf_j].x;//level to which the leaf box that contains jth charge belongs
		unsigned int box_j = leafNodes[leaf_j].y;

		bool is_ij_neighbor = false;
		int typeOfNeighbor, NumOfNeighbor;
		getIfNeighbor(level_i, box_i, level_j, box_j, typeOfNeighbor, NumOfNeighbor);
		if (typeOfNeighbor != -1) {
			is_ij_neighbor = true;
		}

		if (is_ij_neighbor) {
			kernel_dtype entry;
			if (typeOfNeighbor == 0) { //colleagueNeighbors
				entry = colleagueNeighborInteraction[level_i-2][NumOfNeighbor](index_i, index_j);
			}
			else if (typeOfNeighbor == 1) { //coarseNeighbors
				entry = coarseNeighborInteraction[level_i-2][NumOfNeighbor](index_i, index_j);
			}
			else if (typeOfNeighbor == 2) { //fine neighbor
				entry = fineNeighborInteraction[level_i-2][NumOfNeighbor](index_i, index_j);
			}
			else { // separatedFineNeighbor
				entry = separatedFineNeighborInteraction[level_i-2][NumOfNeighbor](index_i, index_j);
			}
			if (findPhi) {
				entry = entry*kappa*kappa*(Q->ContrastFunction(ri));
				if (i ==j) {
					entry = 1.0 + entry;
				}
			}
			return entry;
		}
		else {
			double polar_angle_i = std::atan2(ri.y-tree[level_j][box_j].center.y, ri.x-tree[level_j][box_j].center.x);
			double R2	=	(ri.y-tree[level_j][box_j].center.y)*(ri.y-tree[level_j][box_j].center.y) + (ri.x-tree[level_j][box_j].center.x)*(ri.x-tree[level_j][box_j].center.x);
			double R	=	kappa*sqrt(R2);
			kernel_dtype tempSum = 0.0+0.0*I;
			int seriesLength = FF->seriesLength;
			for (int n = -seriesLength; n <= seriesLength; n++) {
				kernel_dtype temp = (besselJ(n, R) + I*besselY(n, R)) * exp(I*double(n)*polar_angle_i);
				tempSum += M[level_j](n+seriesLength, index_j)*temp;
			}
			tempSum *= I/4.0;
			if (findPhi) {
				tempSum *= kappa*kappa*(Q->ContrastFunction(ri));
				if (i==j) {
					tempSum += 1.0;
				}
			}
			return tempSum;
		}

	}

	kernel_dtype getMatrixEntry(const unsigned i, const unsigned j) {
		//NoIn_getMatrixEntry += 1;
		//double start	=	omp_get_wtime();

		/*
		pts2D ri = particles_X[i];//particles_X is a member of base class FMM_Matrix
		pts2D rj = particles_X[j];//particles_X is a member of base class FMM_Matrix
		double R2 = (ri.x-rj.x)*(ri.x-rj.x) + (ri.y-rj.y)*(ri.y-rj.y);
		double R = sqrt(R2);
		kernel_dtype out = exp(I*kappa*R)/R;
		double end		=	omp_get_wtime();
		timeIn_getMatrixEntry += (end-start);
		if (NoIn_getMatrixEntry == 1) {
			cout << "avg time to get a matrix entry: " << timeIn_getMatrixEntry << endl;
		}
		if (R < epsilonRoundOff) {
			return R;
		}
		else {
			return out;
		}
		*/

		pts2D ri = particles_X[i];//particles_X is a member of base class FMM_Matrix
		unsigned int l = j/rank;//identify to which leaf jth charge belongs
		unsigned int chargeIndex = j%rank;
		unsigned int tj = leafNodes[l].x;//level to which the leaf box that contains jth charge belongs
		unsigned int tk = leafNodes[l].y;
		double polar_angle_i = std::atan2(ri.y-tree[tj][tk].center.y, ri.x-tree[tj][tk].center.x);
		double R2	=	(ri.y-tree[tj][tk].center.y)*(ri.y-tree[tj][tk].center.y) + (ri.x-tree[tj][tk].center.x)*(ri.x-tree[tj][tk].center.x);
		double R	=	kappa*sqrt(R2);
		kernel_dtype tempSum = 0.0+0.0*I;
		int seriesLength = FF->seriesLength;
		for (int n = -seriesLength; n <= seriesLength; n++) {
			kernel_dtype temp = (besselJ(n, R) + I*besselY(n, R)) * exp(I*double(n)*polar_angle_i);
			tempSum += M[tj](n+seriesLength, chargeIndex)*temp;
		}
		tempSum *= I/4.0;
		if (findPhi) {
			tempSum *= kappa*kappa*(Q->ContrastFunction(ri));
			if (i==j) {
				tempSum += 1.0;
			}
		}
		/*double end		=	omp_get_wtime();
		timeIn_getMatrixEntry += (end-start);
		if (NoIn_getMatrixEntry == 1) {
			cout << "avg time to get a matrix entry: " << timeIn_getMatrixEntry << endl;
		}*/
		return tempSum;

	}

	void shift_scale_Nodes(std::vector<pts2D>& Nodes, std::vector<pts2D>& shifted_scaled_Nodes, double xShift, double yShift, double radius) {
		for (int k=0; k < Nodes.size(); ++k) {
			pts2D temp;
			temp.x	=	radius*Nodes[k].x+xShift;
			temp.y	=	radius*Nodes[k].y+yShift;
			shifted_scaled_Nodes.push_back(temp);
		}
	}

	//	get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}

	//	get_S
	double get_S(double x, double y, int n) {
		double S	=	0.5;
		for (int k=1; k<n; ++k) {
			S+=get_ChebPoly(x,k)*get_ChebPoly(y,k);
		}
		return 2.0/n*S;
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
		//	Left Bottom child, i.e., Child 0
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x-0.5;
				temp1.y	=	0.5*temp1.y-0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Right Bottom child, i.e., Child 1
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x+0.5;
				temp1.y	=	0.5*temp1.y-0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Right Top child, i.e., Child 2
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x+0.5;
				temp1.y	=	0.5*temp1.y+0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Left Top child, i.e., Child 3
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x-0.5;
				temp1.y	=	0.5*temp1.y+0.5;
				standardChebNodesChild.push_back(temp1);
		}
	}

	void get_Transfer_Matrix() {
		for (int l=0; l<4; ++l) {
			L2L[l]	=	Eigen::MatrixXd(rank,rank);
			for (int j=0; j<rank; ++j) {
				for (int k=0; k<rank; ++k) {
					L2L[l](j,k)	=	get_S(standardChebNodes[k].x, standardChebNodesChild[j+l*rank].x, nChebNodes)*get_S(standardChebNodes[k].y, standardChebNodesChild[j+l*rank].y, nChebNodes);
				}
			}
		}
		for (int l=0; l<4; ++l) {
			M2M[l]	=	L2L[l].transpose();
		}
	}
	//////////////////TREE CREATION////////////////////////////////////////////
	double resolveContrast(int j, int k) {
		Eigen::VectorXd contrastAtChebNodes;
		std::vector<pts2D> Nodes;
		shift_scale_Nodes(standardChebNodes, Nodes, tree[j][k].center.x, tree[j][k].center.y, tree[j][k].radius);
		EvaluateContrastAtNodes(contrastAtChebNodes, Nodes);
		Eigen::VectorXd approximateContrastAtChildChebNodes(4*rank);
		for (size_t c = 0; c < 4; c++) {
			approximateContrastAtChildChebNodes.segment(rank*c, rank) = L2L[c]*contrastAtChebNodes;
		}
		std::vector<pts2D> NodesChild;
		shift_scale_Nodes(standardChebNodesChild, NodesChild, tree[j][k].center.x, tree[j][k].center.y, tree[j][k].radius);
		Eigen::VectorXd contrastAtChildChebNodes;
		EvaluateContrastAtNodes(contrastAtChildChebNodes, NodesChild);
		double err = errInApproximation(contrastAtChildChebNodes, approximateContrastAtChildChebNodes);
		return err;
	}

		double resolveRHS(int j, int k) {
			Vec RHSAtChebNodes;
			std::vector<pts2D> Nodes;
			shift_scale_Nodes(standardChebNodes, Nodes, tree[j][k].center.x, tree[j][k].center.y, tree[j][k].radius);
			EvaluateRHSAtNodes(RHSAtChebNodes, Nodes);
			Vec approximateRHSAtChildChebNodes(4*rank);
			for (size_t c = 0; c < 4; c++) {
				approximateRHSAtChildChebNodes.segment(rank*c, rank) = L2L[c]*RHSAtChebNodes;
			}
			std::vector<pts2D> NodesChild;
			shift_scale_Nodes(standardChebNodesChild, NodesChild, tree[j][k].center.x, tree[j][k].center.y, tree[j][k].radius);
			Vec RHSAtChildChebNodes;
			EvaluateRHSAtNodes(RHSAtChildChebNodes, NodesChild);
			double err = errInApproximation(RHSAtChildChebNodes, approximateRHSAtChildChebNodes);
			return err;
		}

		void createTree(int j, int k) {
			//k: index of box in the vector tree[j]
			//b: boxNumber of box, tree[j][k]
			double errContrast = resolveContrast(j, k);
			double errRHS = resolveRHS(j, k);
			//cout << "errContrast: " << errContrast << endl;
			bool condition1Leaf = false;
			bool condition2Leaf = true; //need to be false
			bool condition3Leaf = false;
			if (errContrast < epsilonTree) condition1Leaf = true;
			if (errRHS < epsilonTree) condition2Leaf = true;
			if (kappa*tree[j][k].radius/2/PI <= 1.0) condition3Leaf = true; //each leaf shoould contain less than or equal to 1 wave cycle
			//nChebnodes are the number of points in 1 wavecycle/Wavelength - in 1D.
			if (condition1Leaf && condition2Leaf && condition3Leaf) {
				tree[j][k].isLeaf = true;
				orderedPair op;
				op.x = j;
				op.y = k;
				leafNodes.push_back(op);
			}
			else {
				if (int(tree.size()) == j+1) {
					std::vector<FMM2DBox> level;
					tree.push_back(level);
					std::vector<int> index;
				  indexTree.push_back(index);
				}
				int n	=	tree[j+1].size();
				int b = tree[j][k].boxNumber;
				for (size_t c = 0; c < 4; c++) {
					FMM2DBox box;
					box.level = j+1;
					box.boxNumber		=	b*4+c;
					box.parentNumber	=	b;
					box.radius = 0.5*tree[j][k].radius;
					if (c==0) {
						box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
					}
					else if (c==1) {
						box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
					}
					else if (c==2) {
						box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
					}
					else {
						box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
					}
					tree[j+1].push_back(box);
					indexTree[j+1].push_back(box.boxNumber);
				}
				for (size_t c = 0; c < 4; c++) {
					createTree(j+1, n+c);
				}
			}
		}

		void createTreeU(int j, int k) {
			//k: index of box in the vector tree[j]
			//b: boxNumber of box, tree[j][k]
			double errContrast = resolveContrast(j, k);
			double errRHS = resolveRHS(j, k);
			bool condition1Leaf = false;
			bool condition2Leaf = false;
			bool condition3Leaf = false;
			if (errContrast < epsilonTree) condition1Leaf = true;
			if (errRHS < epsilonTree) condition2Leaf = true;
			if (kappa*tree[j][k].radius <= 2*PI*3.0) condition3Leaf = true;
			if (condition3Leaf && j == nLevelsUniform) {
				tree[j][k].isLeaf = true;
				orderedPair op;
				op.x = j;
				op.y = k;
				leafNodes.push_back(op);
			}
			else {
				if (int(tree.size()) == j+1) {
					std::vector<FMM2DBox> level;
					tree.push_back(level);
					std::vector<int> index;
				  indexTree.push_back(index);
				}
				int n	=	tree[j+1].size();
				int b = tree[j][k].boxNumber;
				for (size_t c = 0; c < 4; c++) {
					FMM2DBox box;
					box.level = j+1;
					box.boxNumber		=	b*4+c;
					box.parentNumber	=	b;
					box.radius = 0.5*tree[j][k].radius;
					if (c==0) {
						box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
					}
					else if (c==1) {
						box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
					}
					else if (c==2) {
						box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
					}
					else {
						box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
					}
					tree[j+1].push_back(box);
					indexTree[j+1].push_back(box.boxNumber);
				}
				for (size_t c = 0; c < 4; c++) {
					createTreeU(j+1, n+c);
				}
			}
		}

		void createAdaptiveTree() {
			FMM2DBox root;
			root.level = 0;
			root.boxNumber		=	0;
			root.parentNumber	=	-1;
			root.radius = L;
			root.center.x = 0.0;
			root.center.y = 0.0;
			std::vector<FMM2DBox> level;
			level.push_back(root);
			tree.push_back(level);

			std::vector<int> index;
			index.push_back(0);
			indexTree.push_back(index);
			createTree(0, 0);
			nLevels = tree.size() - 1;
			cout << "nLevels: " << nLevels << endl;
			for (size_t j = 1; j <= 10; j++) {
				boxRadius.push_back(boxRadius[j-1]/2.0);
			}
			if (level_LFR >= nLevels) {
				level_LFR = nLevels-1;
			}
		}

		void createUniformTree() {
			FMM2DBox root;
			root.level = 0;
			root.boxNumber		=	0;
			root.parentNumber	=	-1;
			root.radius = L;
			root.center.x = 0.0;
			root.center.y = 0.0;
			std::vector<FMM2DBox> level;
			level.push_back(root);
			tree.push_back(level);

			std::vector<int> index;
			index.push_back(0);
			indexTree.push_back(index);
			createTreeU(0, 0);
			nLevels = tree.size() - 1;
			cout << "nLevels: " << nLevels << endl;
			for (size_t j = 1; j <= 10; j++) {
				boxRadius.push_back(boxRadius[j-1]/2.0);
			}
			if (level_LFR >= nLevels) {
				level_LFR = nLevels-1;
			}
		}

		void getFutureGeneration_LFR(int bj, int bk, int pnj, int pni, std::vector<orderedPair>& futureGeneration) {
			if (tree[pnj][pni].isLeaf) {
				if ( fabs(tree[bj][bk].center.x - tree[pnj][pni].center.x) >= 3.0*boxRadius[bj] + boxRadius[pnj]-epsilonRoundOff
				 || fabs(tree[bj][bk].center.y - tree[pnj][pni].center.y) >= 3.0*boxRadius[bj] + boxRadius[pnj]-epsilonRoundOff) {
				}
				else {
					orderedPair op;
					op.x = pnj;
					op.y = pni;
					futureGeneration.push_back(op);
				}
			}
			else {
				for (int nc=0; nc<4; ++nc) { //children of parents neighbors
					int pnn = tree[pnj][pni].boxNumber;//boxNumber
					int boxB = 4*pnn+nc;//its index=?
					std::vector<int>::iterator indx = std::find(indexTree[pnj+1].begin(), indexTree[pnj+1].end(), boxB);
					int boxB_index = indx-indexTree[pnj+1].begin();
					getFutureGeneration_LFR(bj, bk, pnj+1, boxB_index, futureGeneration);
				}
			}
		}

		void assign_Child_Interaction_LFR(int c, int j, int k, std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR) {
			int parentboxNumber = tree[j][k].boxNumber;
			int boxA = 4*parentboxNumber+c;//child box number; its index=?
			std::vector<int>::iterator indx = std::find(indexTree[j+1].begin(), indexTree[j+1].end(), boxA);
			int boxA_index = indx-indexTree[j+1].begin();
			for (int n=0; n<Tree_Neighbors_LFR[j][k].size(); ++n) {//parents neighbors; so u need its index to access it which is k
				//children of neighbors of parent which are not neighbors to child=IL
				int pnj = Tree_Neighbors_LFR[j][k][n].x;//level
				int pni = Tree_Neighbors_LFR[j][k][n].y;//index
				int pnn = tree[pnj][pni].boxNumber;//boxNumber

				std::vector<orderedPair> futureGeneration;
				getFutureGeneration_LFR(j+1, boxA_index, pnj, pni, futureGeneration);
				for (size_t d = 0; d < futureGeneration.size(); d++) {
					Tree_Neighbors_LFR[j+1][boxA_index].push_back(futureGeneration[d]);
				}
			}
		}

		//	Assigns the interactions for the children of a box
		void assign_Box_Interactions_LFR(int j, int k, std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR) {
			if (!tree[j][k].isLeaf) {
				//#pragma omp parallel for
				for (int c=0; c<4; ++c) {
					assign_Child_Interaction_LFR(c,j,k, Tree_Neighbors_LFR);
				}
			}
		}

		//	Assigns the interactions for the children all boxes at a given level
		void assign_Level_Interactions_LFR(int j, std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR) {
			//#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				assign_Box_Interactions_LFR(j,k, Tree_Neighbors_LFR);//k is index number of box in tree[j] vector
			}
		}

		//	Assigns the interactions for the children all boxes in the tree
		void assign_Tree_Interactions_LFR(std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR) {
			int J = 1;
			for (int c=0; c<4; ++c) {
				for (int n=0; n<4; ++n) {
					orderedPair op;
					op.x = J;
					op.y = n;
					Tree_Neighbors_LFR[J][c].push_back(op);
				}
			}
			for (int j=1; j<=nLevels-1; ++j) {
				assign_Level_Interactions_LFR(j, Tree_Neighbors_LFR);
			}
		}

		void check() {
			int t = 0;
			std::cout << "t: " << t << "	boxN: " << tree[2][t].boxNumber << std::endl;
			std::cout << "IL: " << std::endl;
			for (size_t i = 0; i < tree[2][t].InteractionList.size(); i++) {
				std::cout << tree[2][t].InteractionList[i].x << ", " << tree[2][t].InteractionList[i].y << std::endl;
			}
			std::cout << "NN: " << std::endl;
			for (size_t i = 0; i < tree[2][t].neighborNumbers.size(); i++) {
				std::cout << tree[2][t].neighborNumbers[i].x << ", " << tree[2][t].neighborNumbers[i].y << std::endl;
			}
		}

		void makeLevelRestriction_Tree() {
			for (size_t n = 0; n < nLevels; n++) {
				std::vector<std::vector<std::vector<orderedPair> > > Tree_Neighbors_LFR;
				for (size_t j = 0; j <= nLevels; j++) {
					std::vector<std::vector<orderedPair> > level;
					for (size_t k = 0; k < tree[j].size(); k++) {
						std::vector<orderedPair> box;
						level.push_back(box);
					}
					Tree_Neighbors_LFR.push_back(level);
				}
				assign_Tree_Interactions_LFR(Tree_Neighbors_LFR);
				std::vector<orderedPair> leafNodesOld = leafNodes;
				for (size_t l = 0; l < leafNodesOld.size(); l++) {
					makeLevelRestriction_Box(l, leafNodes[l], Tree_Neighbors_LFR, leafNodesOld);
				}

				int cnt = 0;
 				for (size_t l = 0; l < leafNodesOld.size(); l++) {//removing all the non-leaves from leafNodes vector
					if (leafNodesOld[l].x == -1) {
						leafNodes.erase(leafNodes.begin()+l-cnt);
						cnt++;
					}
				}

				for (size_t j = 0; j <= nLevels; j++) {
					for (size_t k = 0; k < Tree_Neighbors_LFR[j].size(); k++) {
						Tree_Neighbors_LFR[j][k].clear();
					}
					Tree_Neighbors_LFR[j].clear();
				}
				Tree_Neighbors_LFR.clear();
			}
		}

		int getIndexOfInLeafNodes(orderedPair op) {
			for (size_t l = 0; l < leafNodes.size(); l++) {
				if (op.x == leafNodes[l].x && op.y == leafNodes[l].y) {
					return l;
				}
			}
		}

		void makeLevelRestriction_Box(int l, orderedPair leaf, std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR, std::vector<orderedPair>& leafNodesOld) {
			//l is index of leaf in leafNodes
			/*
			cout << "leafNodesOld.size(): " << leafNodesOld.size() << endl;
			for (size_t r = 0; r < leafNodesOld.size(); r++) {
				cout << r << ", "<< leafNodesOld[r].x << ", " << leafNodesOld[r].y << endl;
			}
			*/
			int j = leaf.x;
			int k = leaf.y;
			for (size_t i = 0; i < Tree_Neighbors_LFR[j][k].size(); i++) {
				int nj = Tree_Neighbors_LFR[j][k][i].x;
				int nk = Tree_Neighbors_LFR[j][k][i].y;
				std::vector<orderedPair> newLeaves;

				if (j-nj >= 2 || nj-j >= 2) {
					if (j-nj >= 2) {
						int indexOfBox = getIndexOfInLeafNodes(Tree_Neighbors_LFR[j][k][i]);//index of Tree_Neighbors_LFR[j][k][i]
						if (leafNodesOld[indexOfBox].x != -1) {
							refineBox(indexOfBox, Tree_Neighbors_LFR[j][k][i], leafNodesOld);
						}
					}
					else if(nj-j >= 2) {
						if (leafNodesOld[l].x != -1) {
							refineBox(l, leaf, leafNodesOld);
						}
					}
				}
			}
		}

		void refineBox(int l, orderedPair leaf, std::vector<orderedPair> &leafNodesOld) {
			int j = leaf.x; //level
			int k = leaf.y; //box index
			int b = tree[j][k].boxNumber;
			leafNodesOld[l].x = -1;//making it a non-leaf node
			tree[j][k].isLeaf = false;
			for (size_t c = 0; c < 4; c++) {
				FMM2DBox box;
				box.isLeaf = true;
				box.level = j+1;
				box.boxNumber		=	b*4+c;
				box.parentNumber	=	b;
				box.radius = 0.5*tree[j][k].radius;
				if (c==0) {
					box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
					box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
				}
				else if (c==1) {
					box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
					box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
				}
				else if (c==2) {
					box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
					box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
				}
				else {
					box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
					box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
				}
				orderedPair op;
				op.x = j+1;
				op.y = tree[j+1].size();
				leafNodes.push_back(op);
				tree[j+1].push_back(box);
				indexTree[j+1].push_back(box.boxNumber);
			}
		}

		void assignLeafCharges() {
			for (size_t k = 0; k < leafNodes.size(); k++) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				tree[j][b].multipoles	=	0.5*(Vec::Ones(rank));//+Eigen::VectorXd::Random(rank));
			}
		}

		void assignLeafCharges(Vec &charges) {
			chargesAll = charges;
			int start = 0;
			for (size_t k = 0; k < leafNodes.size(); k++) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				tree[j][b].multipoles	=	charges.segment(start, rank);
				start += rank;
			}
		}

		void assignLeafChargeLocations(std::vector<pts2D> &particles_out) {
			for (size_t k = 0; k < leafNodes.size(); k++) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				int startIndex = gridPoints.size();
				for (size_t i = 0; i < rank; i++) {
					tree[j][b].chargeLocations.push_back(startIndex+i);
				}
				shift_scale_Nodes(standardChebNodes, gridPoints, tree[j][b].center.x, tree[j][b].center.y, boxRadius[j]);
			}
			particles_X = gridPoints;//object of base class FMM_Matrix
			particles_Y = gridPoints;
			particles_out = gridPoints;
		}

		void collectBoundary() {
			for (size_t t = 0; t < leafNodes.size(); t++) {
				int j = leafNodes[t].x;
				int k = leafNodes[t].y;
				pts2D leftRightBoundaryTemp, bottomTopBoundaryTemp;
				leftRightBoundaryTemp.x = tree[j][k].center.x - tree[j][k].radius;//left
				leftRightBoundaryTemp.y = tree[j][k].center.x + tree[j][k].radius;//right
				bottomTopBoundaryTemp.x = tree[j][k].center.y - tree[j][k].radius;//bottom
				bottomTopBoundaryTemp.y = tree[j][k].center.y + tree[j][k].radius;//top
				leftRightBoundary.push_back(leftRightBoundaryTemp);
				bottomTopBoundary.push_back(bottomTopBoundaryTemp);
			}
		}

		void assignNonLeafChargeLocations() {
			for (int j=nLevels-1; j>1; --j) {
				int J	=	j+1;
				#pragma omp parallel for
				for (int k=0; k<tree[j].size(); ++k) {
					if (!tree[j][k].isLeaf) {
						int b = tree[j][k].boxNumber;
						int KboxNumber;
						int K[4];
						std::vector<int>::iterator indx;
						KboxNumber = 4*b+0;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[0] = indx-indexTree[J].begin();
						KboxNumber = 4*b+1;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[1] = indx-indexTree[J].begin();
						KboxNumber = 4*b+2;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[2] = indx-indexTree[J].begin();
						KboxNumber = 4*b+3;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[3] = indx-indexTree[J].begin();

						for (int c=0; c<4; ++c) {
							//Block containing n elements, starting at position i: vector.segment(i,n)
							for (int i = 0; i < tree[J][K[c]].chargeLocations.size(); i++) {
								tree[j][k].chargeLocations.push_back(tree[J][K[c]].chargeLocations[i]);
							}
						}
					}
				}
			}
		}

		void assignLeafChebNodes() {
			for (size_t k = 0; k < leafNodes.size(); k++) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				shift_scale_Nodes(standardChebNodes, tree[j][b].chebNodes, tree[j][b].center.x, tree[j][b].center.y, boxRadius[j]);
			}
		}

		double errInApproximation(Eigen::VectorXd& trueValue, Eigen::VectorXd& approximateValue) {
			Eigen::VectorXd error(trueValue.size());
			for (int k=0; k<trueValue.size(); ++k) {
				error(k)	=	fabs(trueValue(k)-approximateValue(k));
				//error(k)	=	fabs((trueValue-approximateValue)(k)/trueValue(k));
			}
			return error.maxCoeff();///trueValue.maxCoeff();
		}

		double errInApproximation(Vec& trueValue, Vec& approximateValue) {
			Vec error;
			error = trueValue - approximateValue;
			VectorXd errorAbs = error.cwiseAbs();
			VectorXd trueValueAbs = trueValue.cwiseAbs();
			return errorAbs.maxCoeff();///trueValueAbs.maxCoeff();
		}

		void EvaluateContrastAtNodes(Eigen::VectorXd& contrastAtChebNodes, std::vector<pts2D>& Nodes) {
			contrastAtChebNodes = Eigen::VectorXd(Nodes.size());
			for (size_t i = 0; i < Nodes.size(); i++) {
				contrastAtChebNodes(i) = Q->ContrastFunction(Nodes[i]);
			}
		}

		void EvaluateRHSAtNodes(Vec& RHSAtChebNodes, std::vector<pts2D>& Nodes) {
			RHSAtChebNodes = Vec(Nodes.size());
			for (size_t i = 0; i < Nodes.size(); i++) {
				RHSAtChebNodes(i) = Q->RHSFunction(Nodes[i]);
			}
		}

		void writeMToFile(std::string filename) {
			//create directory
			MFilename = "M";
			string currPath = std::experimental::filesystem::current_path();
			char final1[256];
		  sprintf (final1, "%s/%s", currPath.c_str(), MFilename.c_str());
		  mkdir(final1, 0775);

			MFilename = "M/M_" + std::to_string(nChebNodes) + "_" + std::to_string(int(kappa)) + "_" + std::to_string(degreeOfBases);
			char final2[256];
		  sprintf (final2, "%s/%s", currPath.c_str(), MFilename.c_str());
		  mkdir(final2, 0775);

			filename = MFilename + "/M_" + std::to_string(nChebNodes) + "_" + std::to_string(L);
			for (size_t l = 0; l < M.size(); l++) {
				string filenameModified = filename + "_" + std::to_string(l);
				std::ofstream myfile;
				myfile.open(filenameModified.c_str());
				myfile << M[l] << endl;
			}
		}

		void ReadFromTextFile(Mat &matrix ,std::string filename) {
			std::ifstream inFile (filename,std::ios::in);
			if(!inFile.good()) {
				std::cout<<"Error: could not open file:\""<<filename<<"\"for reading \n";
				exit (2);
			}
			//find the no of values in file
			std::istream_iterator<std::string> in{inFile};
			std::istream_iterator<std::string> end;
			long numberofWords=std::distance(in,end);
			//find the no of lines in file
			inFile.clear();
			inFile.seekg(0,std::ios::beg);
			long numberofLines=std::count(std::istreambuf_iterator<char>(inFile),std::istreambuf_iterator<char>(),'\n');
			//std::cout<<"no. of words : "<<numberofWords<<"numberofLines: "<<numberofLines<<std::endl;
			long rows=numberofLines;
			long cols=numberofWords/numberofLines;
			if(rows*cols!=numberofWords) {
				std::cout<<"\n Infile"<<filename<<"cannot form a matrix \n";
				exit(2);
			}
			matrix.array().resize(rows,cols);
			//matrix Base does not allow resizing ...hence change array base
			inFile.clear();
			inFile.seekg(0,std::ios::beg);
			for(unsigned int i=0;i<matrix.rows();i++)
			 for(unsigned int j=0;j<matrix.cols();j++)
				 inFile>>matrix(i,j);
			inFile.close();
		}

		void readMFromFile(std::string filename) {
			std::vector<Mat> M2;
			MFilename = "M/M_" + std::to_string(nChebNodes) + "_" + std::to_string(int(kappa)) + "_" + std::to_string(degreeOfBases);
			filename = MFilename + "/M_" + std::to_string(nChebNodes) + "_" + std::to_string(L);
			for (size_t l = 0; l <= nLevels; l++) {
				string filenameModified = filename + "_" + std::to_string(l);
				Mat m;
				ReadFromTextFile(m ,filenameModified);
				M.push_back(m);
			}
		}

		void outputAdaptiveGrid(std::string filename) {
			double xcenter = 0.0;
			double ycenter = 0.0;
			double Lx = L;
			double Ly = L;

			std::ofstream myfile;
			myfile.open(filename.c_str());
			myfile << "\\documentclass{standalone}" << std::endl;
			myfile << "\\usepackage{tikz}" << std::endl;
			myfile << "\\begin{document}" << std::endl;
			myfile << "\\begin{tikzpicture}" << std::endl;
			for (int k=0; k<(int)leafNodes.size(); ++k) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				myfile << "\\draw (" << tree[j][b].center.x-tree[j][b].radius << ",";
				// myfile << tree[j][b].center.y-tree[j][b].radius << ") rectangle (";
				//myfile << tree[j][b].center.y-tree[j][b].radius << ") rectangle node{\\tiny " << b << "} (";
				myfile << tree[j][b].center.y-tree[j][b].radius << ") rectangle node{\\tiny " << tree[j][b].boxNumber << "} (";
				myfile << tree[j][b].center.x+tree[j][b].radius << ",";
				myfile << tree[j][b].center.y+tree[j][b].radius << ");" << std::endl;
			}
			long double push	=	0.125;
			myfile<< "\\node at (" << xcenter-Lx-push << "," << ycenter-Ly-push << ") {\\tiny$(" << xcenter-Lx << "," << ycenter-Ly << ")$};" << std::endl;
			myfile<< "\\node at (" << xcenter-Lx-push << "," << ycenter+Ly+push << ") {\\tiny$(" << xcenter-Lx << "," << ycenter+Ly << ")$};" << std::endl;
			myfile<< "\\node at (" << xcenter+Lx+push << "," << ycenter-Ly-push << ") {\\tiny$(" << xcenter+Lx << "," << ycenter-Ly << ")$};" << std::endl;
			myfile<< "\\node at (" << xcenter+Lx+push << "," << ycenter+Ly+push << ") {\\tiny$(" << xcenter+Lx << "," << ycenter+Ly << ")$};" << std::endl;
			myfile << "\\end{tikzpicture}" << std::endl;
			myfile << "\\end{document}" << std::endl;
			myfile.close();
		}
		//////////////////////////////////////////////////////////////

	double find_distance(int j1, int k1, int j2, int k2) {
		pts2D r1 = tree[j1][k1].center;
		pts2D r2 = tree[j2][k2].center;
		return	sqrt((r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y));
	}

	void assign_Child_Interaction(int c, int j, int k) {
		// j: level
		// k: parent box
		int parentboxNumber = tree[j][k].boxNumber;
		int boxA = 4*parentboxNumber+c;//child box number; its index=?
		std::vector<int>::iterator indx = std::find(indexTree[j+1].begin(), indexTree[j+1].end(), boxA);
		int boxA_index = indx-indexTree[j+1].begin();
		for (int n=0; n<tree[j][k].neighborNumbers.size(); ++n) {//parents neighbors; so u need its index to access it which is k
			//children of neighbors of parent which are not neighbors to child=IL
			int pnj = tree[j][k].neighborNumbers[n].x;//level
			int pni = tree[j][k].neighborNumbers[n].y;//index
			int pnn = tree[pnj][pni].boxNumber;//boxNumber
			if (j+1 >= level_LFR || (j+1 < level_LFR && tree[j+1][boxA_index].isLeaf)) { //LFR
				if (tree[pnj][pni].isLeaf) {
					if ( fabs(tree[j+1][boxA_index].center.x - tree[pnj][pni].center.x) >= 3.0*boxRadius[j+1] + boxRadius[pnj]-epsilonRoundOff
					 || fabs(tree[j+1][boxA_index].center.y - tree[pnj][pni].center.y) >= 3.0*boxRadius[j+1] + boxRadius[pnj]-epsilonRoundOff) {
						 orderedPair op;
						 op.x = pnj;
						 op.y = pni;
						 tree[j+1][boxA_index].InteractionList.push_back(op);
					}
					else {
						orderedPair op;
						op.x = pnj;
						op.y = pni;
						tree[j+1][boxA_index].neighborNumbers.push_back(op);
					}
				}
				else {
					for (int nc=0; nc<4; ++nc) { //children of parents neighbors
						int boxB = 4*pnn+nc;//its index=?
						std::vector<int>::iterator indx = std::find(indexTree[j+1].begin(), indexTree[j+1].end(), boxB);
						int boxB_index = indx-indexTree[j+1].begin();
						if ( fabs(tree[j+1][boxA_index].center.x - tree[pnj+1][boxB_index].center.x) >= 3.0*boxRadius[j+1] + boxRadius[pnj+1]-epsilonRoundOff
						 || fabs(tree[j+1][boxA_index].center.y - tree[pnj+1][boxB_index].center.y) >= 3.0*boxRadius[j+1] + boxRadius[pnj+1]-epsilonRoundOff) {
							 orderedPair op;
							 op.x = pnj+1;
							 op.y = boxB_index;
							 tree[j+1][boxA_index].InteractionList.push_back(op);
						}
						else {
							if (!tree[j+1][boxA_index].isLeaf) {
								orderedPair op;
								op.x = pnj+1;
								op.y = boxB_index;
								tree[j+1][boxA_index].neighborNumbers.push_back(op);
							}
							else {
								int pj = pnj+1;
								int pk = boxB_index;

								std::vector<orderedPair> futureGeneration;
								getFutureGeneration_LFR(j+1, boxA_index, pj, pk, futureGeneration);
								for (size_t d = 0; d < futureGeneration.size(); d++) {
									tree[j+1][boxA_index].neighborNumbers.push_back(futureGeneration[d]);
								}
							}
						}
					}
				}
			}
			else { //HFR
				if (tree[pnj][pni].isLeaf) {
					if ( fabs(tree[j+1][boxA_index].center.x - tree[pnj][pni].center.x) >= kappa*boxRadius[j+1]*boxRadius[j+1] + boxRadius[j+1] + boxRadius[pnj]-epsilonRoundOff
					 || fabs(tree[j+1][boxA_index].center.y - tree[pnj][pni].center.y) >= kappa*boxRadius[j+1]*boxRadius[j+1] + boxRadius[j+1] + boxRadius[pnj]-epsilonRoundOff) {
							double arg = atan2(tree[j+1][boxA_index].center.y-tree[pnj][pni].center.y, tree[j+1][boxA_index].center.x-tree[pnj][pni].center.x);
	 						arg = fmod(arg+2*PI+PI, 2*PI);
	 						int coneNum = int(arg/ConeAperture[j+1]);
							orderedPair op;
							op.x = pnj;
							op.y = pni;
							tree[j+1][boxA_index].ConeTree[coneNum].InteractionList.push_back(op);
					}
					else {
						orderedPair op;
						op.x = pnj;
						op.y = pni;
						tree[j+1][boxA_index].neighborNumbers.push_back(op);
					}
				}
				else {
					for (int nc=0; nc<4; ++nc) { //children of parents neighbors
						int boxB = 4*pnn+nc;//its index=?
						std::vector<int>::iterator indx = std::find(indexTree[j+1].begin(), indexTree[j+1].end(), boxB);
						int boxB_index = indx-indexTree[j+1].begin();
						if ( fabs(tree[j+1][boxA_index].center.x - tree[pnj+1][boxB_index].center.x) >= kappa*boxRadius[j+1]*boxRadius[j+1] + boxRadius[j+1] + boxRadius[pnj+1]-epsilonRoundOff
						 || fabs(tree[j+1][boxA_index].center.y - tree[pnj+1][boxB_index].center.y) >= kappa*boxRadius[j+1]*boxRadius[j+1] + boxRadius[pnj+1]-epsilonRoundOff) {
							 double arg = atan2(tree[j+1][boxA_index].center.y-tree[pnj+1][boxB_index].center.y, tree[j+1][boxA_index].center.x-tree[pnj+1][boxB_index].center.x);
	 	 						arg = fmod(arg+2*PI+PI, 2*PI);
	 	 						int coneNum = int(arg/ConeAperture[j+1]);
	 							orderedPair op;
	 							op.x = pnj+1;
	 							op.y = boxB_index;
	 							tree[j+1][boxA_index].ConeTree[coneNum].InteractionList.push_back(op);
						}
						else {
							orderedPair op;
							op.x = pnj+1;
							op.y = boxB_index;
							tree[j+1][boxA_index].neighborNumbers.push_back(op);
						}
					}
				}
			}
		}
	}

	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k) {
		if (!tree[j][k].isLeaf) {
			#pragma omp parallel for
			for (int c=0; c<4; ++c) {
				assign_Child_Interaction(c,j,k);
			}
		}
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j) {
		#pragma omp parallel for
		for (int k=0; k<tree[j].size(); ++k) {
			assign_Box_Interactions(j,k);//k is index number of box in tree[j] vector
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions() {
		//j=0, no neighbors, no IL
		//j=1, no IL, neighbors yes
		//neighbor includes self
		int j = 1;
		for (int c=0; c<4; ++c) {
			for (int n=0; n<4; ++n) {
				orderedPair op;
				op.x = 1;
				op.y = n;
				tree[j][c].neighborNumbers.push_back(op);
			}
		}
		for (j=1; j<=nLevels-1; ++j) {
			assign_Level_Interactions(j);
		}
	}

	void getNeighbors() {//with respect to unit box at origin; this is assumed to leaf
		//for a leaf at level_A; the centers need to scaled up by boxRadius[level_A]
		FF->evaluateQ(Q_pinv);

		for (size_t i = 0; i < 20; i++) {
			SFN_angles[i] = atan2(SFN_centers[i].y, SFN_centers[i].x);
		}
		for (size_t i = 0; i < 9; i++) {
			N_angles[i] = atan2(N_centers[i].y, N_centers[i].x);
		}
		for (size_t i = 0; i < 12; i++) {
			FN_angles[i] = atan2(FN_centers[i].y, FN_centers[i].x);
		}
		for (size_t i = 0; i < 12; i++) {
			CN_angles[i] = atan2(CN_centers[i].y, CN_centers[i].x);
		}
	  /*
		double SFN_angles[20]= {5.0*PI/4.0,
	                    3.0*PI/2.0 - arctan(1.5, 2.5),
	                    3.0*PI/2.0 - arctan(0.5, 2.5),
	                    3.0*PI/2.0 + arctan(0.5, 2.5),
	                    3.0*PI/2.0 + arctan(1.5, 2.5),
	                    7.0*PI/4.0,
	                    2*PI - arctan(1.5, 2.5),
	                    2*PI - arctan(0.5, 2.5),
	                    arctan(0.5, 2.5),
	                    arctan(1.5, 2.5),
	                    PI/4.0,
	                    PI/2.0 - arctan(1.5, 2.5),
	                    PI/2.0 - arctan(0.5, 2.5),
	                    PI/2.0 + arctan(0.5, 2.5),
	                    PI/2.0 + arctan(1.5, 2.5),
	                    3.0*PI/4.0,
	                    PI - arctan(1.5, 2.5),
	                    PI - arctan(0.5, 2.5),
	                    PI + arctan(0.5, 2.5),
	                    PI + arctan(1.5, 2.5)
	                    };

	  double N_angles[9] = {5.0*PI/4.0, 3.0*PI/2.0, 7.0*PI/4.0, 0.0, PI/4.0, PI/2.0, 3.0*PI/4.0, PI, 0.0};
	  double FN_angles[12] = {5.0*PI/4.0,
	                          3.0*PI/2.0 - arctan(1.0, 3.0),
	                          3.0*PI/2.0 + arctan(1.0, 3.0),
	                          7.0*PI/4.0,
	                          2.0*PI - arctan(1.0, 3.0),
	                          arctan(1.0, 3.0),
	                          PI/4.0,
	                          PI/2.0 - arctan(1.0, 3.0),
	                          PI/2.0 + arctan(1.0, 3.0),
	                          3.0*PI/4.0,
	                          PI - arctan(1.0, 3.0),
	                          PI + arctan(1.0, 3.0)
	                         };
	  double CN_angles[12] = {5.0*PI/4.0,
	                          3.0*PI/2.0 - arctan(0.5, 1.5),
	                          3.0*PI/2.0 + arctan(0.5, 1.5),
	                          7.0*PI/4.0,
	                          2.0*PI - arctan(0.5, 1.5),
	                          arctan(0.5, 1.5),
	                          PI/4.0,
	                          arctan(1.5, 0.5),
	                          PI/2.0 + arctan(0.5, 1.5),
	                          3.0*PI/4.0,
	                          PI - arctan(0.5, 1.5),
	                          PI + arctan(0.5, 1.5)
													};
		*/
	  for (size_t j = 1; j <= nLevels; j++) {
	    for (size_t k = 0; k < tree[j].size(); k++) {
	      for (size_t n = 0; n < tree[j][k].neighborNumbers.size(); n++) {
	        int nj = tree[j][k].neighborNumbers[n].x;
	        int nk = tree[j][k].neighborNumbers[n].y;
	        double angle = atan2(tree[nj][nk].center.y-tree[j][k].center.y, tree[nj][nk].center.x-tree[j][k].center.x);
	        if (tree[j][k].radius == tree[nj][nk].radius) {
						if (fabs(tree[nj][nk].center.y-tree[j][k].center.y)<epsilonRoundOff && fabs(tree[nj][nk].center.x-tree[j][k].center.x)<epsilonRoundOff) {
							tree[j][k].colleagueNeighbors[8] = tree[j][k].neighborNumbers[n].y;
						} // self
	          else if (fabs(N_angles[0] - angle) < epsilonRoundOff) {
	            tree[j][k].colleagueNeighbors[0] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(N_angles[1] - angle) < epsilonRoundOff) {
	            tree[j][k].colleagueNeighbors[1] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(N_angles[2] - angle) < epsilonRoundOff) {
	            tree[j][k].colleagueNeighbors[2] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(N_angles[3] - angle) < epsilonRoundOff) {
	            tree[j][k].colleagueNeighbors[3] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(N_angles[4] - angle) < epsilonRoundOff) {
	            tree[j][k].colleagueNeighbors[4] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(N_angles[5] - angle) < epsilonRoundOff) {
	            tree[j][k].colleagueNeighbors[5] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(N_angles[6] - angle) < epsilonRoundOff) {
	            tree[j][k].colleagueNeighbors[6] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(N_angles[7] - angle) < epsilonRoundOff) {
	            tree[j][k].colleagueNeighbors[7] = tree[j][k].neighborNumbers[n].y;
	          }
	        }
	        else if (tree[j][k].radius < tree[nj][nk].radius) { //coarse neighbors
	          if (fabs(CN_angles[0] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[0] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[1] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[1] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[2] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[2] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[3] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[3] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[4] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[4] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[5] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[5] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[6] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[6] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[7] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[7] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[8] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[8] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[9] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[9] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[10] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[10] = tree[j][k].neighborNumbers[n].y;
	          }
	          else if (fabs(CN_angles[11] - angle) < epsilonRoundOff) {
	            tree[j][k].coarseNeighbors[11] = tree[j][k].neighborNumbers[n].y;
	          }
	        }
	        else { //fine neighbors and separated fine neighbors
	           if (find_distance(j, k, nj, nk) > 2.5*tree[j][k].radius) { //separated fine neighbors
	             if (fabs(SFN_angles[0] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[0] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[1] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[1] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[2] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[2] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[3] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[3] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[4] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[4] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[5] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[5] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[6] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[6] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[7] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[7] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[8] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[8] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[9] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[9] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[10] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[10] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[11] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[11] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[12] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[12] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[13] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[13] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[14] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[14] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[15] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[15] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[16] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[16] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[17] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[17] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(SFN_angles[18] - angle) < epsilonRoundOff) {
	               tree[j][k].separatedFineNeighbors[18] = tree[j][k].neighborNumbers[n].y;
	             }
	             else {
	               tree[j][k].separatedFineNeighbors[19] = tree[j][k].neighborNumbers[n].y;
	             }
	           }
	           else {//fine neighbors
	             if (fabs(FN_angles[0] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[0] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[1] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[1] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[2] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[2] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[3] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[3] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[4] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[4] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[5] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[5] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[6] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[6] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[7] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[7] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[8] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[8] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[9] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[9] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[10] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[10] = tree[j][k].neighborNumbers[n].y;
	             }
	             else if (fabs(FN_angles[11] - angle) < epsilonRoundOff) {
	               tree[j][k].fineNeighbors[11] = tree[j][k].neighborNumbers[n].y;
	             }
	           }
	        }
	      }
	    }
	  }
	}

	/*
	//the following function is for helmholtz kernel(uses getMatrixEntry)
	void getOperator(int level_A, int level_B, pts2D center_B, pts2D center_A, Mat &matOperator) {
		std::vector<pts2D> Nodes_A;
		shift_scale_Nodes(standardChebNodes, Nodes_A, center_A.x, center_A.y, boxRadius[level_A]);
		std::vector<pts2D> Nodes_B;
		shift_scale_Nodes(standardChebNodes, Nodes_B, center_B.x, center_B.y, boxRadius[level_B]);

		matOperator = Mat::Zero(rank,rank);
		for (size_t i = 0; i < rank; i++) {
			for (size_t j = 0; j < rank; j++) {
				double R2 = (Nodes_B[j].x-Nodes_A[i].x)*(Nodes_B[j].x-Nodes_A[i].x) + (Nodes_B[j].y-Nodes_A[i].y)*(Nodes_B[j].y-Nodes_A[i].y);
				double R = sqrt(R2);
				if (R < epsilonRoundOff) {
					matOperator(i,j) = R;
				}
				else {
					matOperator(i,j) = exp(I*kappa*R)/R;
					//matOperator(i,j) = 1.0/R;
				}
			}
		}
	}
	*/

	//the following function is for Lippmann Schwinger problem(uses equations specific for neighbor interactions)
	void getOperator(int level_A, int level_B, pts2D center_B, pts2D center_A, Mat &matOperator) {
		integrandInputs II;
		II.center_B = center_B;
		II.center_A = center_A;
		II.level_A = level_A;
		II.level_B = level_B;//colleague Neighbor; so same level
		II.radius_A = boxRadius[level_A];
		II.radius_B = boxRadius[level_B];
		NearFieldInteraction1 *NFI1 = new NearFieldInteraction1(nChebNodes, II, degreeOfBases);
		NearFieldInteraction2 *NFI2 = new NearFieldInteraction2(nChebNodes, II, degreeOfBases);
		Mat V1;
		Mat V2;

		//the following two lines use no expansion of H(0,x); it uses H(0,x)=J(0,x)+iY(0,x)
		NFI1->evaluateV(V1);
		Mat temp = V1;

		//the following lines use series expansion of H(0,x); it uses series expansion of J(0,x) and Y(0,x)
		//NFI1->evaluateV1(V1);
		//NFI2->evaluateV2(V2);
		//Mat temp = V1+V2;


		Mat V_temp = temp*I*boxRadius[level_B]*boxRadius[level_B]/4.0;
		//cout << "V_temp.size(): " << V_temp.rows() << ", " << V_temp.cols() << endl;
		matOperator = V_temp*Q_pinv;
		//cout << "V1: " << endl << V1 << endl << endl;
		//cout << "V2: " << endl << V2 << endl << endl;
		//cout << "Q_pinv: " << endl << Q_pinv << endl << endl;
		//exit(0);
		delete NFI1;
		delete NFI2;
	}


	void getNeighborInteractions() {
		//assume domain is always [-1,1]^{2}
		for (size_t level_A = 2; level_A <= 10; level_A++) {
			writeNeighborInteractionsForLevel(level_A);
		}
	}

	void writeNeighborInteractionsForLevel(int level_A) {
		//B is neighbor of A
		Mat colleagueNeighborInteractionTemp[9];
		Mat fineNeighborInteractionTemp[12];
		Mat separatedFineNeighborInteractionTemp[20];
		Mat coarseNeighborInteractionTemp[12];
		pts2D center_A;
		center_A.x = 0.0;
		center_A.y = 0.0;
		//cout << "level_A: " << level_A << endl;
		for (size_t i = 0; i < 9; i++) { //colleagueNeighbors
			int level_B = level_A;
			pts2D center_B = N_centers[i];
			center_B.x *= boxRadius[level_A];
			center_B.y *= boxRadius[level_A];
			//cout << "colleagueNeighbors:	" << i;
			getOperator(level_A, level_B, center_B, center_A, colleagueNeighborInteractionTemp[i]);
			//cout << "	done" << endl;
		}
		for (size_t i = 0; i < 12; i++) { //coarseNeighbors
			int level_B = level_A-1;
			pts2D center_B = CN_centers[i];
			center_B.x *= boxRadius[level_A];
			center_B.y *= boxRadius[level_A];
			//cout << "coarseNeighbors:	" << i;
			getOperator(level_A, level_B, center_B, center_A, coarseNeighborInteractionTemp[i]);
			//cout << "	done" << endl;
		}
		for (size_t i = 0; i < 20; i++) { //separatedFineNeighbors
			int level_B = level_A+1;
			pts2D center_B = SFN_centers[i];
			center_B.x *= boxRadius[level_A];
			center_B.y *= boxRadius[level_A];
			//cout << "separatedFineNeighbors:	" << i;
			getOperator(level_A, level_B, center_B, center_A, separatedFineNeighborInteractionTemp[i]);
			//cout << "	done" << endl;
		}
		for (size_t i = 0; i < 12; i++) { //fineNeighbors
			int level_B = level_A+1;
			pts2D center_B = FN_centers[i];
			center_B.x *= boxRadius[level_A];
			center_B.y *= boxRadius[level_A];
			//cout << "fineNeighbors:	" << i;
			getOperator(level_A, level_B, center_B, center_A, fineNeighborInteractionTemp[i]);
			//cout << "	done" << endl;
		}

		//create directory
		MFilename = "Neighbor";
		string currPath = std::experimental::filesystem::current_path();
		char final1[256];
		sprintf (final1, "%s/%s", currPath.c_str(), MFilename.c_str());
		mkdir(final1, 0775);

		MFilename = "Neighbor/Neighbor_" + std::to_string(nChebNodes) + "_" + std::to_string(int(kappa)) + "_" + std::to_string(degreeOfBases);
		char final2[256];
		sprintf (final2, "%s/%s", currPath.c_str(), MFilename.c_str());
		mkdir(final2, 0775);

		//writeToFile
		std::string filename;
		filename = MFilename + "/colleagueNeighborInteraction";
		filename = filename + "_" + std::to_string(level_A);
		for (size_t n = 0; n < 9; n++) {
			string filenameModified = filename + "_" + std::to_string(n);
			std::ofstream myfile;
			myfile.open(filenameModified.c_str());
			myfile << colleagueNeighborInteractionTemp[n] << endl;
		}

		filename = MFilename + "/fineNeighborInteraction";
		filename = filename + "_" + std::to_string(level_A);
		for (size_t n = 0; n < 12; n++) {
			string filenameModified = filename + "_" + std::to_string(n);
			std::ofstream myfile;
			myfile.open(filenameModified.c_str());
			myfile << fineNeighborInteractionTemp[n] << endl;
		}

		filename = MFilename + "/separatedFineNeighborInteraction";
		filename = filename + "_" + std::to_string(level_A);
		for (size_t n = 0; n < 20; n++) {
			string filenameModified = filename + "_" + std::to_string(n);
			std::ofstream myfile;
			myfile.open(filenameModified.c_str());
			myfile << separatedFineNeighborInteractionTemp[n] << endl;
		}

		filename = MFilename + "/coarseNeighborInteraction";
		filename = filename + "_" + std::to_string(level_A);
		for (size_t n = 0; n < 12; n++) {
			string filenameModified = filename + "_" + std::to_string(n);
			std::ofstream myfile;
			myfile.open(filenameModified.c_str());
			myfile << coarseNeighborInteractionTemp[n] << endl;
		}
	}

	void readNeighborInteractions() {
		int numberOfLevelOperators = nLevels-2+1;
		std::string filename;
		NeighborFilename = "Neighbor/Neighbor_" + std::to_string(nChebNodes) + "_" + std::to_string(int(kappa)) + "_" + std::to_string(degreeOfBases);
		for (size_t level_A = 2; level_A <= nLevels; level_A++) {
			filename = NeighborFilename + "/colleagueNeighborInteraction";
			filename = filename + "_" + std::to_string(level_A);
			for (size_t n = 0; n < 9; n++) {
				string filenameModified = filename + "_" + std::to_string(n);
				ReadFromTextFile(colleagueNeighborInteraction[level_A-2][n] ,filenameModified);
			}

			filename = NeighborFilename + "/fineNeighborInteraction";
			filename = filename + "_" + std::to_string(level_A);
			for (size_t n = 0; n < 12; n++) {
				string filenameModified = filename + "_" + std::to_string(n);
				ReadFromTextFile(fineNeighborInteraction[level_A-2][n] ,filenameModified);
			}

			filename = NeighborFilename + "/separatedFineNeighborInteraction";
			filename = filename + "_" + std::to_string(level_A);
			for (size_t n = 0; n < 20; n++) {
				string filenameModified = filename + "_" + std::to_string(n);
				ReadFromTextFile(separatedFineNeighborInteraction[level_A-2][n] ,filenameModified);
			}

			filename = NeighborFilename + "/coarseNeighborInteraction";
			filename = filename + "_" + std::to_string(level_A);
			for (size_t n = 0; n < 12; n++) {
				string filenameModified = filename + "_" + std::to_string(n);
				ReadFromTextFile(coarseNeighborInteraction[level_A-2][n] ,filenameModified);
			}
		}
	}

	void createCones() {
		N = leafNodes.size()*rank;
		cout << "level_LFR: " << level_LFR << endl;
		cout << "Number of particles: " << N << endl;
		nCones.push_back(nCones_LFR*pow(2.0, level_LFR-1));
		ConeAperture.push_back(2*PI/nCones[0]);
		for (size_t j = 1; j <= level_LFR-1; j++) {
			nCones.push_back(nCones[j-1]/2.0);
			ConeAperture.push_back(ConeAperture[j-1]*2.0);
		}
		for (size_t j = 0; j <= level_LFR-1; j++) {
			for (int k=0; k<tree[j].size(); ++k) {
				if (!tree[j][k].isLeaf) {
					tree[j][k].ConeTree.clear();////
					for (int c=0; c<nCones[j]; ++c) {
						FMM2DCone cone;
						cone.angle	=	ConeAperture[j]/2.0 + c*ConeAperture[j];
						tree[j][k].ConeTree.push_back(cone);
					}
				}
			}
		}
	}

	void assign_NonLeaf_Charges() {
		for (int j=nLevels-1; j>1; --j) {
			int J	=	j+1;
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (!tree[j][k].isLeaf) {
					int b = tree[j][k].boxNumber;
					int KboxNumber;
					int K[4];
					std::vector<int>::iterator indx;
					KboxNumber = 4*b+0;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[0] = indx-indexTree[J].begin();
					KboxNumber = 4*b+1;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[1] = indx-indexTree[J].begin();
					KboxNumber = 4*b+2;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[2] = indx-indexTree[J].begin();
					KboxNumber = 4*b+3;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[3] = indx-indexTree[J].begin();

					int NumCharges = tree[J][K[0]].chargeLocations.size() + tree[J][K[1]].chargeLocations.size() +tree[J][K[2]].chargeLocations.size() + tree[J][K[3]].chargeLocations.size();
					tree[j][k].multipoles = Vec::Zero(NumCharges);
					int start = 0;
					for (int c=0; c<4; ++c) {
						//Block containing n elements, starting at position i: vector.segment(i,n)
						int NumElem = tree[J][K[c]].chargeLocations.size();
						tree[j][k].multipoles.segment(start, NumElem) = tree[J][K[c]].multipoles;
						start += NumElem;
					}
				}
			}
		}
	}

	void getParticlesFromChildrenLFR_outgoing_col(int j, int k, std::vector<int>& searchNodes) {
		if (tree[j][k].isLeaf) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			for (int c = 0; c < 4; c++) {
				int KboxNumber = 4*b+c;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (J >= level_LFR || (J<level_LFR && tree[J][K].isLeaf)) { //LFR
					// if (tree[J][K].incoming_checkPoints.size() == 0) {
					// 	cout << "problem LFR_outgoing_row: " << j << ", " << k << ", " << J << ", " << K << endl;
					// }
					searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());
				}
				else { //HFR
					for (size_t cone = 0; cone < nCones[J]; cone++) {
						searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone].outgoing_chargePoints.begin(), tree[J][K].ConeTree[cone].outgoing_chargePoints.end());
					}
				}
			}
		}
	}

	void getParticlesFromChildrenLFR_outgoing_row(int j, int k, std::vector<int>& searchNodes, int j1, int k1) {
		if (tree[j][k].isLeaf) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			for (int c = 0; c < 4; c++) {
				int KboxNumber = 4*b+c;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (J >= level_LFR || (J<level_LFR && tree[J][K].isLeaf)) { //LFR
					searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());
				}
				else { //HFR
					for (size_t cone = 0; cone < nCones[J]; cone++) {
						searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone].incoming_checkPoints.begin(), tree[J][K].ConeTree[cone].incoming_checkPoints.end());
					}
				}
			}
		}
	}

	void getNodes_LFR() {
		for (int j=nLevels; j>=2; j--) {
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].incoming_chargePoints.size() > 0)
					tree[j][k].incoming_chargePoints.clear();
				if (tree[j][k].incoming_checkPoints.size() > 0)
					tree[j][k].incoming_checkPoints.clear();
				if (tree[j][k].outgoing_chargePoints.size() > 0)
					tree[j][k].outgoing_chargePoints.clear();
				if (tree[j][k].outgoing_checkPoints.size() > 0)
					tree[j][k].outgoing_checkPoints.clear();
				if (tree[j][k].user_checkPoints.size() > 0)
					tree[j][k].user_checkPoints.clear();
			}
		}
		for (int j=nLevels; j>=level_LFR; j--) {
			getNodes_LFR_incoming_level(j);
			getNodes_LFR_outgoing_level(j);
		}
	}

	void getNodes_LFR_outgoing_box(int j, int k, int &n_rows, int &n_cols, int &ComputedRank) {
		int ILcheck = 0;
		if (tree[j][k].active == true) {
			std::vector<int> boxA_Nodes;
			getParticlesFromChildrenLFR_outgoing_col(j, k, boxA_Nodes);

			//sort( boxA_Nodes.begin(), boxA_Nodes.end() );
			//boxA_Nodes.erase( unique( boxA_Nodes.begin(), boxA_Nodes.end() ), boxA_Nodes.end() );

			std::vector<int> IL_Nodes;//indices
			for (int l=0; l<tree[j][k].InteractionList.size(); ++l) {
				int jIL = tree[j][k].InteractionList[l].x;
				int kIL = tree[j][k].InteractionList[l].y;
				std::vector<int> chargeLocations;
				getParticlesFromChildrenLFR_outgoing_row(jIL, kIL, chargeLocations, j, k);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}

			//sort( IL_Nodes.begin(), IL_Nodes.end() );
			//IL_Nodes.erase( unique( IL_Nodes.begin(), IL_Nodes.end() ), IL_Nodes.end() );

			n_rows = IL_Nodes.size();
			n_cols = boxA_Nodes.size();
			int tol_pow = TOL_POW;
			double tol_ACA = pow(10,-1.0*tol_pow);
			row_indices = IL_Nodes;
			col_indices = boxA_Nodes;//object of base class G_LowRank
			std::vector<int> row_bases, col_bases;
			if (n_rows != 0 && n_cols != 0) {
				tree[j][k].outgoing_ILActive = true;
				Mat dummy;
				ACA_only_nodes(row_bases, col_bases, ComputedRank, tol_ACA, dummy, tree[j][k].outgoing_Ar);
				int minN = n_rows;
				if (n_rows > n_cols) {
					minN = n_cols;
				}
				// if(k==2) {
				// 	std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
				// 	std::vector<int> IL_charges2;//indices
				// 	for (int l=0; l<tree[j][k].InteractionList.size(); ++l) {
				// 		int jIL = tree[j][k].InteractionList[l].x;
				// 		int kIL = tree[j][k].InteractionList[l].y;
				// 		IL_charges2.insert(IL_charges2.end(), tree[jIL][kIL].chargeLocations.begin(), tree[jIL][kIL].chargeLocations.end());
				// 	}
				// 	Mat Afull2 = getMatrix(IL_charges2, tree[j][k].chargeLocations);
				// 	std::cout << "row_indices.size(): " << row_indices.size() << "	col_indices.size(): " << col_indices.size() << std::endl;
				// 	std::cout << "Afull2: " << std::endl << Afull2 << std::endl;
				// 	std::cout << "Afull2.norm(): " << Afull2.norm() << std::endl;
				// }
				// std::cout << "j: " << j << "	k: " << k << "	outgoing ComputedRank: " << ComputedRank << std::endl;
				if (ComputedRank > 0) {
					for (int r = 0; r < row_bases.size(); r++) {
						tree[j][k].outgoing_checkPoints.push_back(IL_Nodes[row_bases[r]]);
					}
					for (int c = 0; c < col_bases.size(); c++) {
						tree[j][k].outgoing_chargePoints.push_back(boxA_Nodes[col_bases[c]]);
					}
					std::vector<int> row_indices_local;
					for (size_t r = 0; r < row_bases.size(); r++) {
						row_indices_local.push_back(IL_Nodes[row_bases[r]]);
					}
					std::vector<int> col_indices_local;
					for (size_t c = 0; c < col_bases.size(); c++) {
						col_indices_local.push_back(boxA_Nodes[col_bases[c]]);
					}
					Mat Atilde = getMatrix(row_indices_local, col_indices_local);
					tree[j][k].outgoing_Atilde_dec = Atilde.colPivHouseholderQr();
					///////////////////////////////////////
					std::vector<int> IL_charges;//indices
					for (int l=0; l<tree[j][k].InteractionList.size(); ++l) {
						int jIL = tree[j][k].InteractionList[l].x;
						int kIL = tree[j][k].InteractionList[l].y;
						IL_charges.insert(IL_charges.end(), tree[jIL][kIL].chargeLocations.begin(), tree[jIL][kIL].chargeLocations.end());
					}
					// Mat Afull = getMatrix(IL_charges, tree[j][k].chargeLocations);
					// Mat C = getMatrix(IL_charges, tree[j][k].outgoing_chargePoints);
					// Mat R = getMatrix(tree[j][k].outgoing_checkPoints, tree[j][k].chargeLocations);
					// Mat A_cur = C*tree[j][k].outgoing_Atilde_dec.solve(R);
					// double nested_error = (Afull-A_cur).norm()/Afull.norm();
					// cout << "Afull.size(): " << Afull.rows() << ", " << Afull.cols() << endl;
					// cout << "nested size: " << row_indices.size() << ", " << col_indices.size() << endl;
					// std::cout << "j: " << j << "	k: " << k << "	outgoing error: " << nested_error << std::endl;
					///////////////////////////////////////
				}
			}
			// if (n_rows == 0) {
			else {
				tree[j][k].outgoing_ILActive = false;
				getParticlesFromChildrenLFR_outgoing_col(j, k, tree[j][k].outgoing_chargePoints);
			}
		}
		else {
			tree[j][k].outgoing_ILActive = false;
		}
	}

	void getNodes_LFR_outgoing_level(int j) { //LFR; box interactions
	  /*
	  computes
	  1. x^{B,o}
	  2. y^{B,o}
	  3. matrix decomposition: K(x^{B,o}, y^{B,o})
	  */
	  //for (int j=nLevels; j>=2; j--) {
	    int rankPerLevel = 0;
	    int n_rows_checkpoint;
	    int n_cols_checkpoint;
	    std::vector<int> boxA_Particles_checkpoint;
	    std::vector<int> IL_Particles_checkpoint;
			int ComputedRank, n_rows, n_cols;
			//#pragma omp parallel for
	    for (int k=0; k<tree[j].size(); ++k) {
				getNodes_LFR_outgoing_box(j, k, n_rows, n_cols, ComputedRank);
				if (rankPerLevel < ComputedRank) {
					rankPerLevel = ComputedRank;
					n_rows_checkpoint = n_rows;
					n_cols_checkpoint = n_cols;
				}
	    }
	  //}
		cout << "O;	j: " << j << "	Nboxes: " << tree[j].size() << "	rows,cols: " << n_rows_checkpoint << "," << n_cols_checkpoint << "	Crank: " << rankPerLevel << endl;
	}


	void getParticlesFromChildrenLFR_incoming_row(int j, int k, std::vector<int>& searchNodes) {
		if (tree[j][k].isLeaf) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			for (int c = 0; c < 4; c++) {
				int KboxNumber = 4*b+c;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (J >= level_LFR || (J<level_LFR && tree[J][K].isLeaf)) { //LFR
					// if (tree[J][K].incoming_checkPoints.size() == 0) {
					// 	cout << "problem LFR_outgoing_row: " << j << ", " << k << ", " << J << ", " << K << endl;
					// }
					searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());
				}
				else { //HFR
					for (size_t cone = 0; cone < nCones[J]; cone++) {
						searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone].incoming_checkPoints.begin(), tree[J][K].ConeTree[cone].incoming_checkPoints.end());
					}
				}
			}
		}
	}

	void getParticlesFromChildrenLFR_incoming_col(int j, int k, std::vector<int>& searchNodes) {
		if (tree[j][k].isLeaf) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			for (int c = 0; c < 4; c++) {
				int KboxNumber = 4*b+c;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (J >= level_LFR || (J<level_LFR && tree[J][K].isLeaf)) { //LFR
					// if (tree[J][K].incoming_checkPoints.size() == 0) {
					// 	cout << "problem LFR_outgoing_row: " << j << ", " << k << ", " << J << ", " << K << endl;
					// }
					searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());
				}
				else { //HFR
					for (size_t cone = 0; cone < nCones[J]; cone++) {
						searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone].outgoing_chargePoints.begin(), tree[J][K].ConeTree[cone].outgoing_chargePoints.end());
					}
				}
			}
		}
	}

	void getNodes_LFR_incoming_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
		int ILcheck = 0;
		if (tree[j][k].active == true) {
			std::vector<int> boxA_Nodes;
			getParticlesFromChildrenLFR_incoming_row(j, k, boxA_Nodes);

			//sort( boxA_Nodes.begin(), boxA_Nodes.end() );
			//boxA_Nodes.erase( unique( boxA_Nodes.begin(), boxA_Nodes.end() ), boxA_Nodes.end() );

			std::vector<int> IL_Nodes;//indices
			for (int l=0; l<tree[j][k].InteractionList.size(); ++l) {
				int jIL = tree[j][k].InteractionList[l].x;
				int kIL = tree[j][k].InteractionList[l].y;
				std::vector<int> chargeLocations;
				getParticlesFromChildrenLFR_incoming_col(jIL, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}

			//sort( IL_Nodes.begin(), IL_Nodes.end() );
			//IL_Nodes.erase( unique( IL_Nodes.begin(), IL_Nodes.end() ), IL_Nodes.end() );

			n_rows = boxA_Nodes.size();
			n_cols = IL_Nodes.size();
			int tol_pow = TOL_POW;
			double tol_ACA = pow(10,-1.0*tol_pow);
			row_indices = boxA_Nodes;//object of base class G_LowRank
			col_indices = IL_Nodes;
			std::vector<int> row_bases, col_bases;
			if (n_rows != 0 && n_cols != 0) {
				tree[j][k].incoming_ILActive = true;
				Mat dummy1, dummy2;
				ACA_only_nodes(row_bases, col_bases, ComputedRank, tol_ACA, dummy1, dummy2);
				// std::cout << "j: " << j << "	k: " << k << "	incoming ComputedRank: " << ComputedRank << std::endl;
				int minN = n_rows;
				if (n_rows > n_cols) {
					minN = n_cols;
				}
				if(ComputedRank > 0) {
					for (int r = 0; r < row_bases.size(); r++) {
						tree[j][k].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
					}
					for (int c = 0; c < col_bases.size(); c++) {
						tree[j][k].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
					}
					std::vector<int> row_indices_local;
					for (size_t r = 0; r < row_bases.size(); r++) {
						row_indices_local.push_back(boxA_Nodes[row_bases[r]]);
					}
					std::vector<int> col_indices_local;
					for (size_t c = 0; c < col_bases.size(); c++) {
						col_indices_local.push_back(IL_Nodes[col_bases[c]]);
					}
					Mat Atilde = getMatrix(row_indices_local, col_indices_local);
					tree[j][k].incoming_Atilde_dec = Atilde.colPivHouseholderQr();
					///////////////////////////////////////
					std::vector<int> IL_charges;//indices
					for (int l=0; l<tree[j][k].InteractionList.size(); ++l) {
						int jIL = tree[j][k].InteractionList[l].x;
						int kIL = tree[j][k].InteractionList[l].y;
						IL_charges.insert(IL_charges.end(), tree[jIL][kIL].chargeLocations.begin(), tree[jIL][kIL].chargeLocations.end());
					}
					// Mat Afull = getMatrix(tree[j][k].chargeLocations, IL_charges);
					// Mat C = getMatrix(tree[j][k].chargeLocations, tree[j][k].incoming_chargePoints);
					// Mat R = getMatrix(tree[j][k].incoming_checkPoints, IL_charges);
					// if(j==3 && k==5) {
					// 	std::cout << "Atilde.determinant(): " << Atilde.determinant() << std::endl;
					// 	std::cout << "Afull.norm(): " << Afull.norm() << std::endl;
					// 	std::cout << "Atilde: " << Atilde << std::endl;
					// 	std::cout << "C: " << C << std::endl;
					// 	std::cout << "R: " << R << std::endl;
					// }
					// Mat A_cur = C*tree[j][k].incoming_Atilde_dec.solve(R);
					// double nested_error = (Afull-A_cur).norm()/Afull.norm();
					// cout << "Afull.size(): " << Afull.rows() << ", " << Afull.cols() << endl;
					// cout << "nested size: " << row_indices.size() << ", " << col_indices.size() << endl;
					// std::cout << "j: " << j << "	k: " << k << "	incoming error: " << nested_error << std::endl;
					///////////////////////////////////////
				}
			}
			// if (n_cols == 0) {
			else {
				tree[j][k].incoming_ILActive = false;
				getParticlesFromChildrenLFR_incoming_row(j, k, tree[j][k].incoming_checkPoints);
			}
		}
		else {
			tree[j][k].incoming_ILActive = false;
		}
	}

	void getNodes_LFR_incoming_level(int j) { //LFR; box interactions
		/*
		computes
		1. x^{B,o}
		2. y^{B,o}
		3. matrix decomposition: K(x^{B,o}, y^{B,o})
		*/
		//for (int j=nLevels; j>=2; j--) {
			int rankPerLevel = 0;
			int n_rows_checkpoint;
			int n_cols_checkpoint;
			std::vector<int> boxA_Particles_checkpoint;
			std::vector<int> IL_Particles_checkpoint;
			int ComputedRank, n_rows, n_cols;
			//#pragma omp parallel for
	    for (int k=0; k<tree[j].size(); ++k) {
				getNodes_LFR_incoming_box(j, k, n_rows, n_cols, ComputedRank);
				if (rankPerLevel < ComputedRank) {
					rankPerLevel = ComputedRank;
					n_rows_checkpoint = n_rows;
					n_cols_checkpoint = n_cols;
				}
			}
		//}
		cout << "I;	j: " << j << "	Nboxes: " << tree[j].size() << "	rows,cols: " << n_rows_checkpoint << "," << n_cols_checkpoint << "	Crank: " << rankPerLevel << endl;
	}

	void LFR_M2M_ILActive_True(int j, int k) {
		std::vector<int> source_points;// = tree[j][k].chargeLocations//source points
		Vec source_densities;
		if (tree[j][k].isLeaf) {
			int Veclength = tree[j][k].multipoles.size();
			source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int l = 0; l < tree[j][k].multipoles.size(); l++) {
				source_densities(g) = tree[j][k].multipoles(l);
				++g;
			}
			source_points.insert(source_points.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			int KboxNumber;
			int K[4];
			std::vector<int>::iterator indx;
			KboxNumber = 4*b+0;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[0] = indx-indexTree[J].begin();
			KboxNumber = 4*b+1;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[1] = indx-indexTree[J].begin();
			KboxNumber = 4*b+2;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[2] = indx-indexTree[J].begin();
			KboxNumber = 4*b+3;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[3] = indx-indexTree[J].begin();

			int Veclength = tree[J][K[0]].outgoing_charges.size()+tree[J][K[1]].outgoing_charges.size()+tree[J][K[2]].outgoing_charges.size()+tree[J][K[3]].outgoing_charges.size();
			source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {
				for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
					source_densities(g) = tree[J][K[child]].outgoing_charges(l);
					++g;
				}
				source_points.insert(source_points.end(), tree[J][K[child]].outgoing_chargePoints.begin(), tree[J][K[child]].outgoing_chargePoints.end());//outgoing_chargePoints
			}
		}
		tree[j][k].outgoing_potential = tree[j][k].outgoing_Ar*source_densities;//u^{B,o}
		tree[j][k].outgoing_charges = tree[j][k].outgoing_Atilde_dec.solve(tree[j][k].outgoing_potential);//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
	}

	void LFR_M2M_ILActive_False(int j, int k) {
		if (tree[j][k].isLeaf) {
			tree[j][k].outgoing_charges = tree[j][k].multipoles;
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			int KboxNumber;
			int K[4];
			std::vector<int>::iterator indx;
			KboxNumber = 4*b+0;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[0] = indx-indexTree[J].begin();
			KboxNumber = 4*b+1;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[1] = indx-indexTree[J].begin();
			KboxNumber = 4*b+2;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[2] = indx-indexTree[J].begin();
			KboxNumber = 4*b+3;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[3] = indx-indexTree[J].begin();

			int Veclength = tree[J][K[0]].outgoing_charges.size()+tree[J][K[1]].outgoing_charges.size()+tree[J][K[2]].outgoing_charges.size()+tree[J][K[3]].outgoing_charges.size();
			tree[j][k].outgoing_charges = Vec::Zero(Veclength);
			int g = 0;
			for (int child = 0; child < 4; child++) {
				for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
					tree[j][k].outgoing_charges(g) = tree[J][K[child]].outgoing_charges(l);
					++g;
				}
			}
		}
	}

	void LFR_M2M() {//outgoing operations
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
	  x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		A = K(x^{B,o}, y^{B,o})
		x^{B,o}=tree[j][k].outgoing_checkPoints
		y^{B,o}=tree[j][k].outgoing_chargePoints
		*/
		for (int j=nLevels; j>=2; --j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (j<level_LFR && !tree[j][k].isLeaf) {
					continue;
				}
				if (tree[j][k].active == true) {
					if (tree[j][k].outgoing_ILActive) {
						if (tree[j][k].outgoing_chargePoints.size() > 0) {
							LFR_M2M_ILActive_True(j, k);
						}
					}
					else {
						LFR_M2M_ILActive_False(j, k);
					}
				}
			}
		}
	}

	void Assemble_LFR_M2L() {
		/*
		u^{B,i}: tree[j][boxB].ConeTree[coneB].incoming_potential
		x^{B,i}: tree[j][boxB].ConeTree[coneB].incoming_checkPoints
		f^{A,o}: tree[j][boxA].ConeTree[coneA].outgoing_charges
		y^{A,o}: tree[j][boxA].ConeTree[coneB].outgoing_chargePoints
		*/
		#pragma omp parallel for
		for (int j=2; j<=nLevels; ++j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {//BoxA
				if (j<level_LFR && !tree[j][k].isLeaf) {
					continue;
				}
				if (tree[j][k].active == true) {
					if (!tree[j][k].isLeaf) {
						//tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].user_checkPoints.size());
						if (tree[j][k].incoming_ILActive) {
							for (int l = 0; l < tree[j][k].InteractionList.size(); l++) {
								int jIL = tree[j][k].InteractionList[l].x;
								int kIL = tree[j][k].InteractionList[l].y;
								if (jIL<level_LFR && !tree[jIL][kIL].isLeaf) {//HFR
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									// if (n_rows != 0 && n_cols != 0) {
										tree[j][k].M2L.push_back(getMatrix(tree[j][k].user_checkPoints, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints));
										//tree[j][k].incoming_potential += R*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									// }
								}
								else {
									int n_rows = tree[j][k].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									tree[j][k].M2L.push_back(getMatrix(tree[j][k].user_checkPoints, tree[jIL][kIL].outgoing_chargePoints));
									//tree[j][k].incoming_potential += R*tree[jIL][kIL].outgoing_charges;
								}
							}
						}
					}
					else {
						//tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].chargeLocations.size());
						if (tree[j][k].incoming_ILActive) {
							for (int l = 0; l < tree[j][k].InteractionList.size(); l++) {
								int jIL = tree[j][k].InteractionList[l].x;
								int kIL = tree[j][k].InteractionList[l].y;

								if (jIL<level_LFR && !tree[jIL][kIL].isLeaf) {//HFR
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].chargeLocations.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();

									// if (n_rows != 0 && n_cols != 0) {
										tree[j][k].M2L.push_back(getMatrix(tree[j][k].chargeLocations, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints));
										//tree[j][k].incoming_potential += R*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									// }
								}
								else {
									int n_rows = tree[j][k].chargeLocations.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									tree[j][k].M2L.push_back(getMatrix(tree[j][k].chargeLocations, tree[jIL][kIL].outgoing_chargePoints));
									//tree[j][k].incoming_potential += R*tree[jIL][kIL].outgoing_charges;
								}
							}
						}
					}
				}
			}
		}
	}

	void LFR_M2L() {
		/*
		u^{B,i}: tree[j][boxB].ConeTree[coneB].incoming_potential
		x^{B,i}: tree[j][boxB].ConeTree[coneB].incoming_checkPoints
		f^{A,o}: tree[j][boxA].ConeTree[coneA].outgoing_charges
		y^{A,o}: tree[j][boxA].ConeTree[coneB].outgoing_chargePoints
		*/
		#pragma omp parallel for
		for (int j=2; j<=nLevels; ++j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {//BoxA
				if (j<level_LFR && !tree[j][k].isLeaf) {
					continue;
				}
				if (tree[j][k].active == true) {
					if (!tree[j][k].isLeaf) {
						tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].user_checkPoints.size());
						if (tree[j][k].incoming_ILActive) {
							#pragma omp parallel for
							for (int l = 0; l < tree[j][k].InteractionList.size(); l++) {
								int jIL = tree[j][k].InteractionList[l].x;
								int kIL = tree[j][k].InteractionList[l].y;
								if (jIL<level_LFR && !tree[jIL][kIL].isLeaf) {//HFR
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									if (n_rows != 0 && n_cols != 0 && RHS_size != 0) {
										//Mat R = getMatrix(tree[j][k].user_checkPoints, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints);
										tree[j][k].incoming_potential += tree[j][k].M2L[l]*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									}
								}
								else {
									int n_rows = tree[j][k].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									//Mat R = getMatrix(tree[j][k].user_checkPoints, tree[jIL][kIL].outgoing_chargePoints);
									tree[j][k].incoming_potential += tree[j][k].M2L[l]*tree[jIL][kIL].outgoing_charges;
								}
							}
						}
					}
					else {
						tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].chargeLocations.size());
						if (tree[j][k].incoming_ILActive) {
							if (tree[j][k].incoming_checkPoints.size() == 0) {
								continue;
							}
							#pragma omp parallel for
							for (int l = 0; l < tree[j][k].InteractionList.size(); l++) {
								int jIL = tree[j][k].InteractionList[l].x;
								int kIL = tree[j][k].InteractionList[l].y;
								if (jIL<level_LFR && !tree[jIL][kIL].isLeaf) {//HFR
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].chargeLocations.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									if (n_rows != 0 && n_cols != 0 && RHS_size != 0) {
										//Mat R = getMatrix(tree[j][k].chargeLocations, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints);
										tree[j][k].incoming_potential += tree[j][k].M2L[l]*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									}
								}
								else {
									int n_rows = tree[j][k].chargeLocations.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									//Mat R = getMatrix(tree[j][k].chargeLocations, tree[jIL][kIL].outgoing_chargePoints);
									tree[j][k].incoming_potential += tree[j][k].M2L[l]*tree[jIL][kIL].outgoing_charges;
								}
							}
						}
					}
				}
			}
		}
	}

	void Assemble_LFR_L2L_IL_True(int j, int k) {
		//tree[j][k].multipoles//source densities
		//tree[j][k].chargeLocations//source points
		//x^{B,o}=tree[j][k].outgoing_checkPoints
		//kernel evaluation between x^{B,o} and source points
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();

		#pragma omp parallel for
		for (int child = 0; child < 4; child++) {
			//x^{C,i}: tree[J][4*k+c].incoming_checkPoints
			//f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
			//y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
			Mat R;
			if (!tree[J][K[child]].isLeaf) {
				tree[J][K[child]].L2L = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].incoming_chargePoints);
			}
			else {
				if (j >= level_LFR) {//LFR
					tree[J][K[child]].L2L = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].incoming_chargePoints);
				}
				else {
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						int n_rows = tree[J][K[child]].chargeLocations.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							tree[J][K[child]].L2L = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						}
					}
				}
			}
		}
	}

	void LFR_L2L_IL_True(int j, int k) {
		//tree[j][k].multipoles//source densities
		//tree[j][k].chargeLocations//source points
		//x^{B,o}=tree[j][k].outgoing_checkPoints
		//kernel evaluation between x^{B,o} and source points
		tree[j][k].incoming_charges = tree[j][k].incoming_Atilde_dec.solve(tree[j][k].incoming_potential);//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();

		#pragma omp parallel for
		for (int child = 0; child < 4; child++) {
			//x^{C,i}: tree[J][4*k+c].incoming_checkPoints
			//f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
			//y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
			Mat R;
			if (!tree[J][K[child]].isLeaf) {
				tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].incoming_charges;//u^{B,o}
			}
			else {
				if (j >= level_LFR) {//LFR
					tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].incoming_charges;//u^{B,o}
				}
				else {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						int n_rows = tree[J][K[child]].chargeLocations.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
				}
			}
		}
	}

	void LFR_L2L_IL_False(int j, int k) {
		int J = j+1;
		int b = tree[j][k].boxNumber;

		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		int offset = 0;

		for (int child = 0; child < 4; child++) {
			if (!tree[J][K[child]].isLeaf) {
				for (size_t i = 0; i < tree[J][K[child]].user_checkPoints.size(); i++) {
					tree[J][K[child]].incoming_potential(i) += tree[j][k].incoming_potential(offset+i);
				}
				offset += tree[J][K[child]].user_checkPoints.size();
			}
			else {
				if (j >= level_LFR) {//LFR
					for (size_t i = 0; i < tree[J][K[child]].chargeLocations.size(); i++) {
						tree[J][K[child]].incoming_potential(i) += tree[j][k].incoming_potential(offset+i);
					}
					offset += tree[J][K[child]].chargeLocations.size();
				}
				else {
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						for (size_t i = 0; i < tree[J][K[child]].chargeLocations.size(); i++) {
							tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
						}
					}
					offset += tree[J][K[child]].chargeLocations.size();
				}
			}
		}
	}

	void LFR_L2L() {//outgoing operations
		for (int j=2; j<nLevels; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (j<level_LFR && !tree[j][k].isLeaf) {
					continue;
				}
				if (!tree[j][k].isLeaf) {
					if (tree[j][k].active == true) {
						if (tree[j][k].incoming_ILActive) {
							if (tree[j][k].incoming_chargePoints.size() > 0) {
								LFR_L2L_IL_True(j, k);
							}
						}
						else {
							LFR_L2L_IL_False(j, k);
						}
					}
				}
			}
		}
	}


		void Assemble_LFR_L2L() {//outgoing operations
			#pragma omp parallel for
			for (int j=2; j<nLevels; ++j) {//parent
				#pragma omp parallel for
				for (int k=0; k<tree[j].size(); ++k) {
					if (j<level_LFR && !tree[j][k].isLeaf) {
						continue;
					}
					if (!tree[j][k].isLeaf) {
						if (tree[j][k].active == true) {
							if (tree[j][k].incoming_ILActive) {
								Assemble_LFR_L2L_IL_True(j, k);
							}
						}
					}
				}
			}
		}

	void evaluate_NearField_Using_Precomputations_H() {
		for (int t=0; t<leafNodes.size(); ++t) {
			int j = leafNodes[t].x;
			int k = leafNodes[t].y;
			if (tree[j][k].active == true) {
				for (size_t nColleagues = 0; nColleagues < 9; nColleagues++) {
					int nj = j;
					int nk = tree[j][k].colleagueNeighbors[nColleagues];
					if (nk != -1) {
						Mat boxOperator = colleagueNeighborInteraction[j-2][nColleagues];
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				for (size_t nFine = 0; nFine < 12; nFine++) {
					int nj = j+1;
					int nk = tree[j][k].fineNeighbors[nFine];
					if (nk != -1) {
						Mat boxOperator = fineNeighborInteraction[j-2][nFine];
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				for (size_t nSFine = 0; nSFine < 20; nSFine++) {
					int nj = j+1;
					int nk = tree[j][k].separatedFineNeighbors[nSFine];
					if (nk != -1) {
						Mat boxOperator = separatedFineNeighborInteraction[j-2][nSFine];
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				for (size_t nCoarse = 0; nCoarse < 12; nCoarse++) {
					int nj = j-1;
					int nk = tree[j][k].coarseNeighbors[nCoarse];
					if (nk != -1) {
						Mat boxOperator = coarseNeighborInteraction[j-2][nCoarse];
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
			}
		}
	}

	void evaluate_NearField_Using_Precomputations_LS() {
		#pragma omp parallel for
		for (int t=0; t<leafNodes.size(); ++t) {
			int j = leafNodes[t].x;
			int k = leafNodes[t].y;
			if (tree[j][k].active == true) {
				VectorXd contrastVec(rank);
				#pragma omp parallel for
				for (size_t i = 0; i < rank; i++) {
					contrastVec(i) = kappa*kappa*(Q->ContrastFunction(tree[j][k].chebNodes[i]));
				}
				MatrixXd contrastMat = contrastVec.asDiagonal();
				#pragma omp parallel for
				for (size_t nColleagues = 0; nColleagues < 9; nColleagues++) {
					int nj = j;
					int nk = tree[j][k].colleagueNeighbors[nColleagues];
					if (nk != -1) {
						Mat boxOperator = MatrixXd::Zero(rank,rank);;
						if (findPhi) {
							boxOperator = contrastMat*colleagueNeighborInteraction[j-2][nColleagues];
							// boxOperator = colleagueNeighborInteraction[j-2][nColleagues];
							if (nColleagues == 8) //self interaction
								boxOperator = boxOperator + MatrixXd::Identity(rank,rank);
						}
						else {
							boxOperator = colleagueNeighborInteraction[j-2][nColleagues];
						}
						// if (k==5) {
						// 	std::cout << "boxOperator" << std::endl << boxOperator << std::endl;
						// 	std::cout << "tree[nj][nk].multipoles" << std::endl << tree[nj][nk].multipoles << std::endl;
						// 	// std::cout << "boxOperator" << std::endl << boxOperator << std::endl;
						// }
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				#pragma omp parallel for
				for (size_t nFine = 0; nFine < 12; nFine++) {
					int nj = j+1;
					int nk = tree[j][k].fineNeighbors[nFine];
					Mat boxOperator;
					if (nk != -1) {
						if (findPhi) {
							boxOperator = contrastMat*fineNeighborInteraction[j-2][nFine];
						}
						else {
							boxOperator = fineNeighborInteraction[j-2][nFine];
						}
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				#pragma omp parallel for
				for (size_t nSFine = 0; nSFine < 20; nSFine++) {
					int nj = j+1;
					int nk = tree[j][k].separatedFineNeighbors[nSFine];
					Mat boxOperator;
					if (nk != -1) {
						if (findPhi) {
							boxOperator = contrastMat*separatedFineNeighborInteraction[j-2][nSFine];
						}
						else {
							boxOperator = separatedFineNeighborInteraction[j-2][nSFine];
						}
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				#pragma omp parallel for
				for (size_t nCoarse = 0; nCoarse < 12; nCoarse++) {
					int nj = j-1;
					int nk = tree[j][k].coarseNeighbors[nCoarse];
					Mat boxOperator;
					if (nk != -1) {
						if (findPhi) {
							boxOperator = contrastMat*coarseNeighborInteraction[j-2][nCoarse];
						}
						else {
							boxOperator = coarseNeighborInteraction[j-2][nCoarse];
						}
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
			}
		}
	}

	void evaluate_NearField() {
	//we are always making sure that we have a tree high enough to be in LFR;
	//so, near filed needs to be done only in LFR
		//#pragma omp parallel for
		for (int t=0; t<leafNodes.size(); ++t) {
			int j = leafNodes[t].x;
			int k = leafNodes[t].y;
			if (tree[j][k].active == true) {
				//Neighbor Interaction
				for (int l = 0; l < tree[j][k].neighborNumbers.size(); l++) {//excluding self
					int nj = tree[j][k].neighborNumbers[l].x;
					int nk = tree[j][k].neighborNumbers[l].y;
					int n_rows = tree[j][k].chargeLocations.size();
					int n_cols = tree[nj][nk].chargeLocations.size();
					Mat R = getMatrix(tree[j][k].chargeLocations, tree[nj][nk].chargeLocations);
					tree[j][k].incoming_potential += R*tree[nj][nk].multipoles;//u^{B,o}
				}
			}
		}
	}

	double perform_Error_Check(int t) {
		int j = leafNodes[t].x;//boxA
		int k = leafNodes[t].y;
		if (tree[j][k].active == true) {
			Vec potential	=	Vec::Zero(tree[j][k].chargeLocations.size());
			for (int p=0; p<leafNodes.size(); ++p) {
				int nj = leafNodes[p].x;//boxB
				int nk = leafNodes[p].y;
				Mat NaiveOper;
				getOperator(j, nj, tree[nj][nk].center, tree[j][k].center, NaiveOper);

				VectorXd contrastVec(rank);
				for (size_t i = 0; i < rank; i++) {
					contrastVec(i) = kappa*kappa*(Q->ContrastFunction(tree[j][k].chebNodes[i]));
				}
				MatrixXd contrastMat = contrastVec.asDiagonal();
				NaiveOper = contrastMat*NaiveOper;
				if (j==nj && k==nk) {
					NaiveOper = NaiveOper + MatrixXd::Identity(rank,rank);
				}

				potential = potential + NaiveOper*tree[nj][nk].multipoles;
			}
			Eigen::VectorXd error(tree[j][k].chargeLocations.size());
			for (int p=0; p<tree[j][k].chargeLocations.size(); ++p) {
				error(p)	=	abs(potential(p)-tree[j][k].incoming_potential(p));
			}
			VectorXd absTruePotential = potential.cwiseAbs();
			VectorXd absCalcPotential = tree[j][k].incoming_potential.cwiseAbs();
			return error.maxCoeff()/absTruePotential.maxCoeff();
		}
		else {
			return 0.0;
		}
	}

	// double perform_Error_Check2(int t) {
	// 	int lj = leafNodes[t].x;//boxA
	// 	int lk = leafNodes[t].y;
	// 	if (tree[lj][lk].active == true) {
	// 		Vec truePotential	=	Vec::Zero(tree[lj][lk].chargeLocations.size());
	// 		for (size_t t = 0; t < truePotential.size(); t++) {
	// 			int i = tree[lj][lk].chargeLocations[t];
	// 			for (int j = 0; j < N; j++) {
	// 				truePotential(t) = truePotential(t) + getMatrixEntry2(i,j)*chargesAll[j];
	// 			}
	// 		}
	// 		Eigen::VectorXd error(tree[lj][lk].chargeLocations.size());
	// 		for (int p=0; p<tree[lj][lk].chargeLocations.size(); ++p) {
	// 			error(p)	=	abs(truePotential(p)-tree[lj][lk].incoming_potential(p));
	// 		}
	// 		VectorXd absTruePotential = truePotential.cwiseAbs();
	// 		VectorXd absCalcPotential = tree[lj][lk].incoming_potential.cwiseAbs();
	// 		return error.maxCoeff()/absTruePotential.maxCoeff();
	// 	}
	// 	else {
	// 		return 0.0;
	// 	}
	// }

	double perform_Error_Check2(int t) {
		int lj = leafNodes[t].x;//boxA
		int lk = leafNodes[t].y;
		if (tree[lj][lk].active == true) {
			Vec truePotential	=	Vec::Zero(tree[lj][lk].chargeLocations.size());
				// for (size_t y = 0; y < tree[2][t].InteractionList.size(); y++) {//box
				// 	int il = tree[2][t].InteractionList[y].y;
				// 	for (size_t u = 0; u < tree[2][il].chargeLocations.size(); u++) {
				// 		int ind = tree[2][il].chargeLocations[u];
				// 		truePotential(t) = truePotential(t) + getMatrixEntry2(i,ind)*chargesAll[ind];
				// 	}
				// }
				for (size_t y = 0; y < tree[lj][lk].neighborNumbers.size(); y++) { //box
					int il = tree[lj][lk].neighborNumbers[y].y;
					int jl = tree[lj][lk].neighborNumbers[y].x;
					Mat check(tree[lj][lk].chargeLocations.size(),tree[jl][il].chargeLocations.size());
					// if (il==lk) {
						for (size_t t = 0; t < truePotential.size(); t++) {
							int i = tree[lj][lk].chargeLocations[t];
							for (size_t u = 0; u < tree[jl][il].chargeLocations.size(); u++) {
								int ind = tree[jl][il].chargeLocations[u];
								// truePotential(t) = truePotential(t) + getMatrixEntry2(i,ind)*chargesAll[ind];
								check(t,u) = getMatrixEntry2(i,ind);
								// if(lk==5) {
								// 	std::cout << chargesAll[ind] << std::endl;
								// }
							}
						}
						truePotential += check*tree[jl][il].multipoles;
						// if(lk==5) {
						// 	std::cout << "Mat" << std::endl << check << std::endl;
						// }
					// }
				}
				for (size_t y = 0; y < tree[lj][lk].InteractionList.size(); y++) { //box
					int il = tree[lj][lk].InteractionList[y].y;
					int jl = tree[lj][lk].InteractionList[y].x;
					Mat check(tree[lj][lk].chargeLocations.size(),tree[jl][il].chargeLocations.size());
					// if (il==lk) {
						for (size_t t = 0; t < truePotential.size(); t++) {
							int i = tree[lj][lk].chargeLocations[t];
							for (size_t u = 0; u < tree[jl][il].chargeLocations.size(); u++) {
								int ind = tree[jl][il].chargeLocations[u];
								// truePotential(t) = truePotential(t) + getMatrixEntry2(i,ind)*chargesAll[ind];
								check(t,u) = getMatrixEntry2(i,ind);
								// if(lk==5) {
								// 	std::cout << chargesAll[ind] << std::endl;
								// }
							}
						}
						truePotential += check*tree[jl][il].multipoles;
						// if(lk==5) {
						// 	std::cout << "Mat" << std::endl << check << std::endl;
						// }
					// }
				}
				// for (int j = 0; j < N; j++) {
				// 	truePotential(t) = truePotential(t) + getMatrixEntry2(i,j)*chargesAll[j];
				// }
			// if(lk==5) {
			// 	// std::cout << "multipoles" << std::endl << tree[lj][lk].multipoles << std::endl;
			// 	std::cout << "truePotential" << std::endl << truePotential << std::endl;
			// 	// std::cout << "truePotential2" << std::endl << check*tree[lj][lk].multipoles << std::endl;
			// 	std::cout << "tree[lj][lk].incoming_potential" << std::endl << tree[lj][lk].incoming_potential << std::endl;
			// }
			Eigen::VectorXd error(tree[lj][lk].chargeLocations.size());
			for (int p=0; p<tree[lj][lk].chargeLocations.size(); ++p) {
				error(p)	=	abs(truePotential(p)-tree[lj][lk].incoming_potential(p));
			}
			VectorXd absTruePotential = truePotential.cwiseAbs();
			VectorXd absCalcPotential = tree[lj][lk].incoming_potential.cwiseAbs();
			// return error.maxCoeff()/absTruePotential.maxCoeff();
			// if (t==0) {
			// 	std::cout << "truePotential" << std::endl << truePotential << std::endl;
			// 	std::cout << "tree[lj][lk].incoming_potential" << std::endl << tree[lj][lk].incoming_potential << std::endl;
			// }
			return error.norm()/truePotential.norm();
		}
		else {
			return 0.0;
		}
	}

	void collectPotential(Vec &potential) {
		potential = VectorXcd::Zero(N);
		int start = 0;
		for (size_t t = 0; t < leafNodes.size(); t++) {
			int j = leafNodes[t].x;
			int k = leafNodes[t].y;
			potential.segment(start, rank) = tree[j][k].incoming_potential;
			start += rank;
		}
	}

	double findTotalError(Vec &trueRHS, Mat &S) {
		Vec DAFMMpotential;
		collectPotential(DAFMMpotential);
		// Mat S(N,N);
		S = Mat(N,N);
		for (size_t i = 0; i < N; i++) {
			for (size_t j = 0; j < N; j++) {
				S(i,j) = getMatrixEntry2(i,j);
			}
		}
		trueRHS = S*chargesAll;
		Vec trueErr = trueRHS - DAFMMpotential;
		double DAFMMerr = trueErr.norm()/trueRHS.norm();
		return DAFMMerr;
		// std::cout << "DAFMMerr: " << DAFMMerr << std::endl;
	}
//////////////////////////////////////////////////////////////////
/*
void getUserCheckPoints() {
	for (size_t j = level_LFR; j <= nLevels; j++) {//LFR
		for (size_t k = 0; k < tree[j].size(); k++) {
			if (tree[j][k].isLeaf) {
				tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
			}
			else {
				if (tree[j][k].ILActive) {
					tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].incoming_checkPoints.begin(), tree[j][k].incoming_checkPoints.end());
				}
				else {
					int J = j+1;
					int b = tree[j][k].boxNumber;
					int KboxNumber;
					int K[4];
					std::vector<int>::iterator indx;
					KboxNumber = 4*b+0;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[0] = indx-indexTree[J].begin();
					KboxNumber = 4*b+1;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[1] = indx-indexTree[J].begin();
					KboxNumber = 4*b+2;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[2] = indx-indexTree[J].begin();
					KboxNumber = 4*b+3;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[3] = indx-indexTree[J].begin();
					for (size_t ch = 0; ch < 4; ch++) {
						if (tree[J][K[ch]].isLeaf) {
							tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[J][K[ch]].chargeLocations.begin(), tree[J][K[ch]].chargeLocations.end());
						}
						else {
							tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[J][K[ch]].incoming_checkPoints.begin(), tree[J][K[ch]].incoming_checkPoints.end());
						}
					}
				}
			}
		}
	}
	for (size_t j = 2; j < level_LFR; j++) {//HFR
		for (size_t k = 0; k < tree[j].size(); k++) {
			if (tree[j][k].isLeaf) {
				tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
			}
			else {//HFR
				for (size_t cone = 0; cone < nCones[j]; cone++) {
					if (tree[j][k].ConeTree[cone].ILActive) {
						tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[j][k].ConeTree[cone].incoming_checkPoints.begin(), tree[j][k].ConeTree[cone].incoming_checkPoints.end());
					}
					else {
						int J = j+1;
						int b = tree[j][k].boxNumber;
						int KboxNumber;
						int K[4];
						std::vector<int>::iterator indx;
						KboxNumber = 4*b+0;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[0] = indx-indexTree[J].begin();
						KboxNumber = 4*b+1;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[1] = indx-indexTree[J].begin();
						KboxNumber = 4*b+2;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[2] = indx-indexTree[J].begin();
						KboxNumber = 4*b+3;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[3] = indx-indexTree[J].begin();
						int cone_child =  cone/2;
						for (size_t ch = 0; ch < 4; ch++) {
							if (tree[J][K[ch]].isLeaf) {
								tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].chargeLocations.begin(), tree[J][K[ch]].chargeLocations.end());
							}
							else {
								if (J == level_LFR) {
									tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].incoming_checkPoints.begin(), tree[J][K[ch]].incoming_checkPoints.end());
								}
								else {
									tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].ConeTree[cone_child].incoming_checkPoints.begin(), tree[J][K[ch]].ConeTree[cone_child].incoming_checkPoints.end());
								}
							}
						}
					}
				}
			}
		}
	}
}
*/
void getUserCheckPoints() {
	for (size_t j = nLevels; j >= level_LFR; j--) {//LFR
		#pragma omp parallel for
		for (size_t k = 0; k < tree[j].size(); k++) {
			if (tree[j][k].isLeaf) {
				tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
			}
			else {
				if (tree[j][k].incoming_ILActive) {
					tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].incoming_checkPoints.begin(), tree[j][k].incoming_checkPoints.end());
				}
				else {
					int J = j+1;
					int b = tree[j][k].boxNumber;
					int KboxNumber;
					int K[4];
					std::vector<int>::iterator indx;
					KboxNumber = 4*b+0;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[0] = indx-indexTree[J].begin();
					KboxNumber = 4*b+1;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[1] = indx-indexTree[J].begin();
					KboxNumber = 4*b+2;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[2] = indx-indexTree[J].begin();
					KboxNumber = 4*b+3;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[3] = indx-indexTree[J].begin();
					for (size_t ch = 0; ch < 4; ch++) {
						tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[J][K[ch]].user_checkPoints.begin(), tree[J][K[ch]].user_checkPoints.end());
					}
				}
			}
		}
	}
	for (size_t j = level_LFR-1; j >= 2; j--) {//HFR
		#pragma omp parallel for
		for (size_t k = 0; k < tree[j].size(); k++) {
			if (tree[j][k].isLeaf) {
				tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
			}
			else {//HFR
				for (size_t cone = 0; cone < nCones[j]; cone++) {
					if (tree[j][k].ConeTree[cone].incoming_ILActive) {
						tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[j][k].ConeTree[cone].incoming_checkPoints.begin(), tree[j][k].ConeTree[cone].incoming_checkPoints.end());
					}
					else {
						int J = j+1;
						int b = tree[j][k].boxNumber;
						int KboxNumber;
						int K[4];
						std::vector<int>::iterator indx;
						KboxNumber = 4*b+0;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[0] = indx-indexTree[J].begin();
						KboxNumber = 4*b+1;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[1] = indx-indexTree[J].begin();
						KboxNumber = 4*b+2;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[2] = indx-indexTree[J].begin();
						KboxNumber = 4*b+3;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[3] = indx-indexTree[J].begin();
						int cone_child =  cone/2;
						for (size_t ch = 0; ch < 4; ch++) {
							if (tree[J][K[ch]].isLeaf) {
								tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].user_checkPoints.begin(), tree[J][K[ch]].user_checkPoints.end());
							}
							else {
								if (J == level_LFR) {
									tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].user_checkPoints.begin(), tree[J][K[ch]].user_checkPoints.end());
								}
								else {
									tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].ConeTree[cone_child].user_checkPoints.begin(), tree[J][K[ch]].ConeTree[cone_child].user_checkPoints.end());
								}
							}
						}
					}
				}
			}
		}
	}
}

void getNodes_HFR() {
	for (int j=level_LFR-1; j>=2; --j) {
		for (int k=0; k<tree[j].size(); ++k) {
			if (!tree[j][k].isLeaf) {
				for (int cone=0; cone<nCones[j]; ++cone) {
					if (tree[j][k].ConeTree[cone].incoming_chargePoints.size() > 0)
						tree[j][k].ConeTree[cone].incoming_chargePoints.clear();
					if (tree[j][k].ConeTree[cone].incoming_checkPoints.size() > 0)
						tree[j][k].ConeTree[cone].incoming_checkPoints.clear();
					if (tree[j][k].ConeTree[cone].outgoing_chargePoints.size() > 0)
						tree[j][k].ConeTree[cone].outgoing_chargePoints.clear();
					if (tree[j][k].ConeTree[cone].outgoing_checkPoints.size() > 0)
						tree[j][k].ConeTree[cone].outgoing_checkPoints.clear();
					if (tree[j][k].ConeTree[cone].user_checkPoints.size() > 0)
						tree[j][k].ConeTree[cone].user_checkPoints.clear();
				}
			}
		}
	}
	for (int j=level_LFR-1; j>=2; --j) {
		getNodes_HFR_outgoing_level(j);
		getNodes_HFR_incoming_level(j);
	}
}


void getParticlesFromChildrenHFR_outgoing_row(int j, int k, int cone_parent, std::vector<int>& searchNodes) {
	if (tree[j][k].isLeaf) {
		searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
	}
	else {
		int J = j+1;//child
		int b = tree[j][k].boxNumber;
		if (j == level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (tree[J][K].isLeaf) {
					searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());//outgoing_chargePoints
				}
				else {
					searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone_child].incoming_checkPoints.begin(), tree[J][K].ConeTree[cone_child].incoming_checkPoints.end());//outgoing_chargePoints
				}
			}
		}
	}
}

void getParticlesFromChildrenHFR_outgoing_col(int j, int k, int cone_parent, std::vector<int>& searchNodes) {
	if (tree[j][k].isLeaf) {
		searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
	}
	else {
		int J = j+1;//child
		int b = tree[j][k].boxNumber;
		if (j == level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (tree[J][K].isLeaf) {
					searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());//outgoing_chargePoints
				}
				else {
					searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone_child].outgoing_chargePoints.begin(), tree[J][K].ConeTree[cone_child].outgoing_chargePoints.end());//outgoing_chargePoints
				}
			}
		}
	}
}

void getNodes_HFR_outgoing_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
	if (tree[j][k].active == true) {
		for (int cone=0; cone<nCones[j]; ++cone) {
			std::vector<int> boxA_Nodes;
			getParticlesFromChildrenHFR_outgoing_col(j, k, cone, boxA_Nodes);

			//sort( boxA_Nodes.begin(), boxA_Nodes.end() );
			//boxA_Nodes.erase( unique( boxA_Nodes.begin(), boxA_Nodes.end() ), boxA_Nodes.end() );

			std::vector<int> IL_Nodes;
			for (int b = 0; b < tree[j][k].ConeTree[cone].InteractionList.size(); b++) {
					int jB = tree[j][k].ConeTree[cone].InteractionList[b].x;
					int boxB = tree[j][k].ConeTree[cone].InteractionList[b].y;
					double arg = atan2(tree[j][k].center.y-tree[jB][boxB].center.y, tree[j][k].center.x-tree[jB][boxB].center.x);
					double argB = fmod(arg+2*PI, 2*PI);
					int coneB = int(argB/ConeAperture[jB]);//=coneB; direction of nBox wrt k
					std::vector<int> chargeLocations;
					getParticlesFromChildrenHFR_outgoing_row(jB, boxB, coneB, chargeLocations);
					IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}

			//sort( IL_Nodes.begin(), IL_Nodes.end() );
			//IL_Nodes.erase( unique( IL_Nodes.begin(), IL_Nodes.end() ), IL_Nodes.end() );

			n_rows = IL_Nodes.size();
			n_cols = boxA_Nodes.size();
			int tol_pow = TOL_POW;
			double tol_ACA = pow(10,-1.0*tol_pow);
			row_indices = IL_Nodes;
			col_indices = boxA_Nodes;
			std::vector<int> row_bases, col_bases;
			if (n_rows != 0 && n_cols != 0) {
				tree[j][k].ConeTree[cone].outgoing_ILActive = true;
				Mat dummy;
				ACA_only_nodes(row_bases, col_bases, ComputedRank, tol_ACA, dummy, tree[j][k].ConeTree[cone].outgoing_Ar);
				if (ComputedRank > 0) {
					for (int r = 0; r < row_bases.size(); r++) {
						tree[j][k].ConeTree[cone].outgoing_checkPoints.push_back(IL_Nodes[row_bases[r]]);
					}
					for (int c = 0; c < col_bases.size(); c++) {
						tree[j][k].ConeTree[cone].outgoing_chargePoints.push_back(boxA_Nodes[col_bases[c]]);
					}
					std::vector<int> row_indices_local;
					for (size_t r = 0; r < row_bases.size(); r++) {
						row_indices_local.push_back(IL_Nodes[row_bases[r]]);
					}
					std::vector<int> col_indices_local;
					for (size_t c = 0; c < col_bases.size(); c++) {
						col_indices_local.push_back(boxA_Nodes[col_bases[c]]);
					}
					Mat Atilde = getMatrix(row_indices_local, col_indices_local);
					tree[j][k].ConeTree[cone].outgoing_Atilde_dec = Atilde.colPivHouseholderQr();
				}
			}
			// if (n_rows == 0) {
			else {
				tree[j][k].ConeTree[cone].outgoing_ILActive = false;
				getParticlesFromChildrenHFR_outgoing_col(j, k, cone, tree[j][k].ConeTree[cone].outgoing_chargePoints);
			}
		}
	}
}

void getNodes_HFR_outgoing_level(int j) { //HFR; cone interactions
	//for (int j=level_LFR-1; j>=2; --j) {
		int rankPerLevel = 0;
		int n_rows_checkpoint;
		int n_cols_checkpoint;
		int kMax;
		int n_rows, n_cols, ComputedRank;
		std::vector<int> boxA_Particles_checkpoint;
		std::vector<int> IL_Particles_checkpoint;
		for (int k=0; k<tree[j].size(); ++k) {
			if (tree[j][k].isLeaf) {
				getNodes_LFR_outgoing_box(j, k, n_rows, n_cols, ComputedRank);
			}
			else {
				getNodes_HFR_outgoing_box(j, k, n_rows, n_cols, ComputedRank);
			}
			if (rankPerLevel < ComputedRank) {
				rankPerLevel = ComputedRank;
				n_rows_checkpoint = n_rows;
				n_cols_checkpoint = n_cols;
				kMax = k;
			}
		}
	//}
		cout << "O;	j: " << j << "	Nboxes: " << tree[j].size() << "	k: " << kMax << "	rows,cols: " << n_rows_checkpoint << "," << n_cols_checkpoint << "	Crank: " << rankPerLevel << endl;
	}

void getParticlesFromChildrenHFR_incoming_row(int j, int k, int cone_parent, std::vector<int>& searchNodes) {
	if (tree[j][k].isLeaf) {
		searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
	}
	else {
		int J = j+1;//child
		int b = tree[j][k].boxNumber;
		if (j == level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (tree[J][K].isLeaf) {
					searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());//outgoing_chargePoints
				}
				else {
					searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone_child].incoming_checkPoints.begin(), tree[J][K].ConeTree[cone_child].incoming_checkPoints.end());//outgoing_chargePoints
				}
			}
		}
	}
}

void getParticlesFromChildrenHFR_incoming_col(int j, int k, int cone_parent, std::vector<int>& searchNodes) {
	if (tree[j][k].isLeaf) {
		searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
	}
	else {
		int J = j+1;//child
		int b = tree[j][k].boxNumber;
		if (j == level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (tree[J][K].isLeaf) {
					searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());//outgoing_chargePoints
				}
				else {
					searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone_child].outgoing_chargePoints.begin(), tree[J][K].ConeTree[cone_child].outgoing_chargePoints.end());//outgoing_chargePoints
				}
			}
		}
	}
}

void getNodes_HFR_incoming_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
	if (tree[j][k].active == true) {
		for (int cone=0; cone<nCones[j]; ++cone) {
			std::vector<int> boxA_Nodes;
			getParticlesFromChildrenHFR_incoming_row(j, k, cone, boxA_Nodes);

			//sort( boxA_Nodes.begin(), boxA_Nodes.end() );
			//boxA_Nodes.erase( unique( boxA_Nodes.begin(), boxA_Nodes.end() ), boxA_Nodes.end() );

			std::vector<int> IL_Nodes;
			for (int b = 0; b < tree[j][k].ConeTree[cone].InteractionList.size(); b++) {
					int jB = tree[j][k].ConeTree[cone].InteractionList[b].x;
					int boxB = tree[j][k].ConeTree[cone].InteractionList[b].y;
					double arg = atan2(tree[j][k].center.y-tree[jB][boxB].center.y, tree[j][k].center.x-tree[jB][boxB].center.x);
					double argB = fmod(arg+2*PI, 2*PI);
					int coneB = int(argB/ConeAperture[jB]);//=coneB; direction of nBox wrt k
					std::vector<int> chargeLocations;
					getParticlesFromChildrenHFR_incoming_col(jB, boxB, coneB, chargeLocations);
					IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}

			//sort( IL_Nodes.begin(), IL_Nodes.end() );
			//IL_Nodes.erase( unique( IL_Nodes.begin(), IL_Nodes.end() ), IL_Nodes.end() );

			n_rows = boxA_Nodes.size();
			n_cols = IL_Nodes.size();
			int tol_pow = TOL_POW;
			double tol_ACA = pow(10,-1.0*tol_pow);
			row_indices = boxA_Nodes;
			col_indices = IL_Nodes;
			std::vector<int> row_bases, col_bases;
			if (n_rows != 0 && n_cols != 0) {
				tree[j][k].ConeTree[cone].incoming_ILActive = true;
				Mat dummy1, dummy2;
				ACA_only_nodes(row_bases, col_bases, ComputedRank, tol_ACA, dummy1, dummy2);
				if (ComputedRank > 0) {
					for (int r = 0; r < row_bases.size(); r++) {
						tree[j][k].ConeTree[cone].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
					}
					for (int c = 0; c < col_bases.size(); c++) {
						tree[j][k].ConeTree[cone].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
					}
					std::vector<int> row_indices_local;
					for (size_t r = 0; r < row_bases.size(); r++) {
						row_indices_local.push_back(boxA_Nodes[row_bases[r]]);
					}
					std::vector<int> col_indices_local;
					for (size_t c = 0; c < col_bases.size(); c++) {
						col_indices_local.push_back(IL_Nodes[col_bases[c]]);
					}
					Mat Atilde = getMatrix(row_indices_local, col_indices_local);
					tree[j][k].ConeTree[cone].incoming_Atilde_dec = Atilde.colPivHouseholderQr();
				}
			}
			// if (n_cols == 0) {
			else {
				tree[j][k].ConeTree[cone].incoming_ILActive = false;
				getParticlesFromChildrenHFR_incoming_row(j, k, cone, tree[j][k].ConeTree[cone].incoming_checkPoints);
			}
		}
	}
}

void getNodes_HFR_incoming_level(int j) { //HFR; cone interactions
		int rankPerLevel = 0;
		int n_rows_checkpoint;
		int n_cols_checkpoint;
		int kMax;
		std::vector<int> boxA_Particles_checkpoint;
		std::vector<int> IL_Particles_checkpoint;
		int n_rows, n_cols, ComputedRank;
		for (int k=0; k<tree[j].size(); ++k) {
			if (tree[j][k].isLeaf){
				getNodes_LFR_incoming_box(j, k, n_rows, n_cols, ComputedRank);
			}
			else {
				getNodes_HFR_incoming_box(j, k, n_rows, n_cols, ComputedRank);
			}
			if (rankPerLevel < ComputedRank) {
				rankPerLevel 			= ComputedRank;
				n_rows_checkpoint = n_rows;
				n_cols_checkpoint = n_cols;
				kMax 							= k;
			}
		}
	cout << "I;	j: " << j << "	Nboxes: " << tree[j].size() << "	k: " << kMax << "	rows,cols: " << n_rows_checkpoint << "," << n_cols_checkpoint << "	Crank: " << rankPerLevel << endl;
}

	void Assemble_HFR_M2M_ILActiveTrue(int j, int k, int cone_parent) {
		/*
		A = K(x^{B,o}, y^{B,o})
		x^{B,o}=tree[j][k].outgoing_checkPoints
		y^{B,o}=tree[j][k].outgoing_chargePoints
		*/
		std::vector<int> source_points;// = tree[j][k].chargeLocations//source points
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		if (j==level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				source_points.insert(source_points.end(), tree[J][K[child]].outgoing_chargePoints.begin(), tree[J][K[child]].outgoing_chargePoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			int Veclength = 0;
			for (int child = 0; child < 4; child++) {
				if (tree[J][K[child]].isLeaf) {
					source_points.insert(source_points.end(), tree[J][K[child]].outgoing_chargePoints.begin(), tree[J][K[child]].outgoing_chargePoints.end());//outgoing_chargePoints
				}
				else {
					source_points.insert(source_points.end(), tree[J][K[child]].ConeTree[cone_child].outgoing_chargePoints.begin(), tree[J][K[child]].ConeTree[cone_child].outgoing_chargePoints.end());//outgoing_chargePoints
				}
			}
		}
		int n_rows = tree[j][k].ConeTree[cone_parent].outgoing_checkPoints.size();//outgoing_checkPoints
		int n_cols = source_points.size();
		if (n_rows != 0 && n_cols != 0) {
			tree[j][k].ConeTree[cone_parent].M2M = getMatrix(tree[j][k].ConeTree[cone_parent].outgoing_checkPoints, source_points);
		}
	}

	void HFR_M2M_ILActiveTrue(int j, int k, int cone_parent) {
		/*
		A = K(x^{B,o}, y^{B,o})
		x^{B,o}=tree[j][k].outgoing_checkPoints
		y^{B,o}=tree[j][k].outgoing_chargePoints
		*/
		std::vector<int> source_points;// = tree[j][k].chargeLocations//source points
		Vec source_densities;
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		if (j==level_LFR-1) {
			int Veclength = tree[J][K[0]].outgoing_charges.size()+tree[J][K[1]].outgoing_charges.size()+tree[J][K[2]].outgoing_charges.size()+tree[J][K[3]].outgoing_charges.size();
			source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {//child
				for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
					source_densities(g) = tree[J][K[child]].outgoing_charges(l);
					++g;
				}
			}
			/*for (int child = 0; child < 4; child++) {
				source_points.insert(source_points.end(), tree[J][K[child]].outgoing_chargePoints.begin(), tree[J][K[child]].outgoing_chargePoints.end());//outgoing_chargePoints
			}*/
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			int Veclength = 0;
			#pragma omp parallel for
			for (int child = 0; child < 4; child++) {//child
				if (tree[J][K[child]].isLeaf) {
					Veclength += tree[J][K[child]].outgoing_charges.size();
				}
				else {
					Veclength += tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size();
				}
			}
			source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {//child
				if (tree[J][K[child]].isLeaf) {
					for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
						source_densities(g) = tree[J][K[child]].outgoing_charges(l);
						++g;
					}
				}
				else {
					for (int l = 0; l < tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size(); l++) {
						source_densities(g) = tree[J][K[child]].ConeTree[cone_child].outgoing_charges(l);
						++g;
					}
				}
			}
			/*for (int child = 0; child < 4; child++) {
				if (tree[J][K[child]].isLeaf) {
					source_points.insert(source_points.end(), tree[J][K[child]].outgoing_chargePoints.begin(), tree[J][K[child]].outgoing_chargePoints.end());//outgoing_chargePoints
				}
				else {
					source_points.insert(source_points.end(), tree[J][K[child]].ConeTree[cone_child].outgoing_chargePoints.begin(), tree[J][K[child]].ConeTree[cone_child].outgoing_chargePoints.end());//outgoing_chargePoints
				}
			}*/
		}
		int n_rows = tree[j][k].ConeTree[cone_parent].outgoing_checkPoints.size();//outgoing_checkPoints
		int n_cols = source_densities.size();
		if (n_rows != 0 && n_cols != 0) {
			/*Mat R = getMatrix(tree[j][k].ConeTree[cone_parent].outgoing_checkPoints, source_points);
			Mat Err = tree[j][k].ConeTree[cone_parent].outgoing_Ar-R;
			if (Err.norm() != 0.0) {
				cout << "j: " << j << "	k: " << k << "	c: " << cone_parent << "	Er: " << Err.norm() << "	Ar.n: " << tree[j][k].ConeTree[cone_parent].outgoing_Ar.norm() << ", R.n: " << R.norm() << endl;
			}
			tree[j][k].ConeTree[cone_parent].outgoing_potential = R*source_densities;//u^{B,o}*/
			tree[j][k].ConeTree[cone_parent].outgoing_potential = tree[j][k].ConeTree[cone_parent].M2M*source_densities;//u^{B,o}
			tree[j][k].ConeTree[cone_parent].outgoing_charges = tree[j][k].ConeTree[cone_parent].outgoing_Atilde_dec.solve(tree[j][k].ConeTree[cone_parent].outgoing_potential);//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
		}
	}

	void HFR_M2M_ILActiveFalse(int j, int k, int cone_parent) {
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		if (j==level_LFR-1) {
			int Veclength = tree[J][K[0]].outgoing_charges.size()+tree[J][K[1]].outgoing_charges.size()+tree[J][K[2]].outgoing_charges.size()+tree[J][K[3]].outgoing_charges.size();
			tree[j][k].ConeTree[cone_parent].outgoing_charges = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {//child
				for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
					tree[j][k].ConeTree[cone_parent].outgoing_charges(g) = tree[J][K[child]].outgoing_charges(l);
					++g;
				}
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			int Veclength = 0;
			#pragma omp parallel for
			for (int child = 0; child < 4; child++) {//child
				if (tree[J][K[child]].isLeaf) {
					Veclength += tree[J][K[child]].outgoing_charges.size();
				}
				else {
					Veclength += tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size();
				}
			}
			tree[j][k].ConeTree[cone_parent].outgoing_charges = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {//child
				if (tree[J][K[child]].isLeaf) {
					for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
						tree[j][k].ConeTree[cone_parent].outgoing_charges(g) = tree[J][K[child]].outgoing_charges(l);
						++g;
					}
				}
				else {
					for (int l = 0; l < tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size(); l++) {
						tree[j][k].ConeTree[cone_parent].outgoing_charges(g) = tree[J][K[child]].ConeTree[cone_child].outgoing_charges(l);
						++g;
					}
				}
			}
		}
	}

	void HFR_M2M() {//outgoing operations
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
		x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		*/
		for (int j=level_LFR-1; j>=2; --j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active == true) {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						if (tree[j][k].ConeTree[cone_parent].outgoing_ILActive) {
							if (tree[j][k].ConeTree[cone_parent].outgoing_chargePoints.size() > 0) {
								HFR_M2M_ILActiveTrue(j, k, cone_parent);
							}
						}
						else {
							HFR_M2M_ILActiveFalse(j, k, cone_parent);
						}
					}
				}
			}
		}
	}

	void Assemble_HFR_M2M() {//outgoing operations
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
		x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		*/
		for (int j=level_LFR-1; j>=2; --j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active == true) {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						if (tree[j][k].ConeTree[cone_parent].outgoing_ILActive) {
							Assemble_HFR_M2M_ILActiveTrue(j, k, cone_parent);
						}
					}
				}
			}
		}
	}

	void Assemble_HFR_M2L() {
		/*
		l:tree[j][boxA].ConeTree[coneA]; l':tree[j][boxB].ConeTree[coneB]
		u^{B,i,l'}: tree[j][boxB].ConeTree[coneB].incoming_potential
		x^{B,i,l'}: tree[j][boxB].ConeTree[coneB].incoming_checkPoints
		f^{A,o,l}: tree[j][boxA].ConeTree[coneA].outgoing_charges
		y^{A,o,l}: tree[j][boxA].ConeTree[coneA].outgoing_chargePoints
		*/
		#pragma omp parallel for
		for (int j=2; j<level_LFR; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k < tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active) {
					#pragma omp parallel for
					for (int coneB = 0; coneB < nCones[j]; coneB++) {
						//tree[j][k].ConeTree[coneB].incoming_potential = Vec::Zero(tree[j][k].ConeTree[coneB].user_checkPoints.size());
						if (tree[j][k].ConeTree[coneB].incoming_ILActive) {
							int numIL = tree[j][k].ConeTree[coneB].InteractionList.size();
							for (int l = 0; l < numIL; l++) {
								int jIL = tree[j][k].ConeTree[coneB].InteractionList[l].x;
								int kIL = tree[j][k].ConeTree[coneB].InteractionList[l].y;
								if (tree[jIL][kIL].isLeaf) {//LFR
									int n_rows = tree[j][k].ConeTree[coneB].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].outgoing_charges.size();
									// if (n_rows != 0 && n_cols != 0) {
										tree[j][k].ConeTree[coneB].M2L.push_back(getMatrix(tree[j][k].ConeTree[coneB].user_checkPoints, tree[jIL][kIL].outgoing_chargePoints));
										//tree[j][k].ConeTree[coneB].incoming_potential += R*tree[jIL][kIL].outgoing_charges;//u^{B,o}
									// }
								}
								else {
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].ConeTree[coneB].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									// if (n_rows != 0 && n_cols != 0) {
										tree[j][k].ConeTree[coneB].M2L.push_back(getMatrix(tree[j][k].ConeTree[coneB].user_checkPoints, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints));
										//tree[j][k].ConeTree[coneB].incoming_potential += R*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									// }
								}
							}
						}
					}
				}
			}
		}
	}

	void HFR_M2L() {
		/*
		l:tree[j][boxA].ConeTree[coneA]; l':tree[j][boxB].ConeTree[coneB]
		u^{B,i,l'}: tree[j][boxB].ConeTree[coneB].incoming_potential
		x^{B,i,l'}: tree[j][boxB].ConeTree[coneB].incoming_checkPoints
		f^{A,o,l}: tree[j][boxA].ConeTree[coneA].outgoing_charges
		y^{A,o,l}: tree[j][boxA].ConeTree[coneA].outgoing_chargePoints
		*/
		#pragma omp parallel for
		for (int j=2; j<level_LFR; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k < tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active) {
					#pragma omp parallel for
					for (int coneB = 0; coneB < nCones[j]; coneB++) {
						tree[j][k].ConeTree[coneB].incoming_potential = Vec::Zero(tree[j][k].ConeTree[coneB].user_checkPoints.size());
						if (tree[j][k].ConeTree[coneB].incoming_ILActive) {
							int numIL = tree[j][k].ConeTree[coneB].InteractionList.size();
							#pragma omp parallel for
							for (int l = 0; l < numIL; l++) {
								int jIL = tree[j][k].ConeTree[coneB].InteractionList[l].x;
								int kIL = tree[j][k].ConeTree[coneB].InteractionList[l].y;
								if (tree[jIL][kIL].isLeaf) {//LFR
									int n_rows = tree[j][k].ConeTree[coneB].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].outgoing_charges.size();
									if (n_rows != 0 && n_cols != 0 && RHS_size != 0) {
										/*Mat R = getMatrix(tree[j][k].ConeTree[coneB].user_checkPoints, tree[jIL][kIL].outgoing_chargePoints);
										Mat Err = tree[j][k].ConeTree[coneB].M2L[l]-R;
										if (Err.norm() != 0.0) {
											cout << "j: " << j << "	k: " << k << "	coneB: " << coneB << "	Err: " << Err.norm() << endl;
										}*/
										tree[j][k].ConeTree[coneB].incoming_potential += tree[j][k].ConeTree[coneB].M2L[l]*tree[jIL][kIL].outgoing_charges;//u^{B,o}
									}
								}
								else {
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].ConeTree[coneB].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									if (n_rows != 0 && n_cols != 0 && RHS_size != 0) {
										/*Mat R = getMatrix(tree[j][k].ConeTree[coneB].user_checkPoints, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints);
										Mat Err = tree[j][k].ConeTree[coneB].M2L[l]-R;
										if (Err.norm() != 0.0) {
											cout << "j: " << j << "	k: " << k << "	coneB: " << coneB << "	Err: " << Err.norm() << endl;
										}*/
										tree[j][k].ConeTree[coneB].incoming_potential += tree[j][k].ConeTree[coneB].M2L[l]*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	void Assemble_HFR_L2L_ILActiveTrue (int j, int k, int cone_parent) {
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		#pragma omp parallel for
		for (int child = 0; child < 4; child++) {
			//cout << "J: " << J << "	K: " << K[child] << "	iL: " << tree[J][K[child]].isLeaf << endl;
			if (j==level_LFR-1) {
				if (level_LFR != nLevels) {
					if (!tree[J][K[child]].isLeaf) {
						int n_rows = tree[J][K[child]].user_checkPoints.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
							//tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
					else {
						int n_rows = tree[J][K[child]].chargeLocations.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
							//tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
				}
				else {
					int n_rows = tree[J][K[child]].chargeLocations.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						//tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
				}
			}
			else {
				int cone_child = cone_parent/2;
				if (!tree[J][K[child]].isLeaf) {
					//cout << "J: " << J << "	K: " << K[child] << "	rc: " << tree[J][K[child]].ConeTree[cone_child].L2L.rows() << ", " << tree[J][K[child]].ConeTree[cone_child].L2L.cols() << "	rhs: " << tree[j][k].ConeTree[cone_parent].incoming_charges.size() << endl;

					int n_rows = tree[J][K[child]].ConeTree[cone_child].user_checkPoints.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].ConeTree[cone_child].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						//tree[J][K[child]].ConeTree[cone_child].incoming_potential += tree[J][K[child]].ConeTree[cone_child].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
					//cout << "done" << endl;
				}
				else {
					int n_rows = tree[J][K[child]].user_checkPoints.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					//cout << "J: " << J << "	K: " << K[child] << "	rc: " << tree[J][K[child]].L2L.rows() << ", " << tree[J][K[child]].L2L.cols() << "	rhs: " << tree[j][k].ConeTree[cone_parent].incoming_charges.size() << ", " << tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size() << endl;
					if (n_rows != 0 && n_cols != 0) {
						tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						//cout << "R.rc: " << tree[J][K[child]].L2L.rows() << ", " << tree[J][K[child]].L2L.cols() << endl;
						//tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
					//cout << "done" << endl;
				}
			}
		}
	}

	void HFR_L2L_ILActiveTrue (int j, int k, int cone_parent) {
		/*
		x^{C,i}: tree[J][4*k+c].incoming_checkPoints
		f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
		y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
		x^{C,i}: tree[J][4*k+c].incoming_checkPoints
		f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
		y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
		x^{C,i,l'}: tree[J][4*k+c].ConeTree[cone_child].incoming_checkPoints
		f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
		y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
		*/
		tree[j][k].ConeTree[cone_parent].incoming_charges = tree[j][k].ConeTree[cone_parent].incoming_Atilde_dec.solve(tree[j][k].ConeTree[cone_parent].incoming_potential);//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		#pragma omp parallel for
		for (int child = 0; child < 4; child++) {
			if (j==level_LFR-1) {
				if (level_LFR != nLevels) {
					if (!tree[J][K[child]].isLeaf) {
						int n_rows = tree[J][K[child]].user_checkPoints.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							/*Mat R = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
							Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
							if (Err.norm() != 0.0) {
								cout << "J: " << J << "	K[child]: " << K[child] << "	cone_parent: " << cone_parent << "	Err: " << Err.norm() << endl;
							}*/
							tree[J][K[child]].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
					else {
						int n_rows = tree[J][K[child]].chargeLocations.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							/*Mat R = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
							Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
							if (Err.norm() != 0.0) {
								cout << "J: " << J << "	K[child]: " << K[child] << "	cone_parent: " << cone_parent << "	Err: " << Err.norm() << endl;
							}*/
							tree[J][K[child]].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
				}
				else {
					int n_rows = tree[J][K[child]].chargeLocations.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						/*Mat R = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
						if (Err.norm() != 0.0) {
							cout << "J: " << J << "	K[child]: " << K[child] << "	cone_parent: " << cone_parent << "	Err: " << Err.norm() << endl;
						}*/
						tree[J][K[child]].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
				}
			}
			else {
				int cone_child = cone_parent/2;
				if (!tree[J][K[child]].isLeaf) {
					int n_rows = tree[J][K[child]].ConeTree[cone_child].user_checkPoints.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						/*Mat R = getMatrix(tree[J][K[child]].ConeTree[cone_child].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
						if (Err.norm() != 0.0) {
							cout << "J: " << J << "	K[child]: " << K[child] << "	cone_child: " << cone_child << "	Err: " << Err.norm() << endl;
						}*/
						tree[J][K[child]].ConeTree[cone_child].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
				}
				else {
					int n_rows = tree[J][K[child]].user_checkPoints.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						/*Mat R = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
						if (Err.norm() != 0.0) {
							cout << "J: " << J << "	K[child]: " << K[child] << "	cone_child: " << cone_child << "	Err: " << Err.norm() << endl;
						}*/
						tree[J][K[child]].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
				}
			}
		}
	}

	void HFR_L2L_ILActiveFalse(int j, int k, int cone_parent) {
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		int offset = 0;
		for (int child = 0; child < 4; child++) {
			if (j==level_LFR-1) {
				if (level_LFR != nLevels) {
					if (!tree[J][K[child]].isLeaf) {
						for (size_t i = 0; i < tree[J][K[child]].user_checkPoints.size(); i++) {
							tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
						}
						offset += tree[J][K[child]].user_checkPoints.size();
					}
					else {
						for (size_t i = 0; i < tree[J][K[child]].chargeLocations.size(); i++) {
							tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
						}
						offset += tree[J][K[child]].chargeLocations.size();
					}
				}
				else {
					for (size_t i = 0; i < tree[J][K[child]].chargeLocations.size(); i++) {
						tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
					}
					offset += tree[J][K[child]].chargeLocations.size();
				}
			}
			else {
				int cone_child = cone_parent/2;
				if (tree[J][K[child]].isLeaf) {
					for (size_t i = 0; i < tree[J][K[child]].user_checkPoints.size(); i++) {
						tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
					}
					offset += tree[J][K[child]].user_checkPoints.size();
				}
				else {
					for (size_t i = 0; i < tree[J][K[child]].ConeTree[cone_child].user_checkPoints.size(); i++) {
						tree[J][K[child]].ConeTree[cone_child].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
					}
					offset += tree[J][K[child]].ConeTree[cone_child].user_checkPoints.size();
				}
			}
		}
	}

	void Assemble_HFR_L2L() {//outgoing operations; finds locals untill level: level_LFR
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
		x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		*/
		for (int j=2; j<level_LFR; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active == true && !tree[j][k].isLeaf) {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						if (tree[j][k].ConeTree[cone_parent].incoming_ILActive) {
							Assemble_HFR_L2L_ILActiveTrue(j,k,cone_parent);
						}
					}
				}
			}
		}
	}

	void HFR_L2L() {//outgoing operations; finds locals untill level: level_LFR
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
		x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		*/
		for (int j=2; j<level_LFR; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active == true && !tree[j][k].isLeaf) {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						if (tree[j][k].ConeTree[cone_parent].incoming_ILActive) {
							if (tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size() > 0) {
								HFR_L2L_ILActiveTrue(j,k,cone_parent);
							}
						}
						else {
							HFR_L2L_ILActiveFalse(j,k,cone_parent);
						}
					}
				}
			}
		}
	}
//////////////////////////////////////////////////////////////////
};

#endif
