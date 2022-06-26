//#ifndef _singularNodes_HPP__
//#define _singularNodes_HPP__
class FMM_Matrix {
public:
	std::vector<pts2D> particles_X;
	std::vector<pts2D> particles_Y;

	FMM_Matrix(std::vector<pts2D> particles_X, std::vector<pts2D> particles_Y){
			this->particles_X = particles_X;
			this->particles_Y = particles_Y;
	}

	virtual kernel_dtype getMatrixEntry(const unsigned i, const unsigned j) {
		cout << "virtual getInteraction" << endl;
		return 0.0;
	}
	/*
	kernel_dtype getMatrixEntry(const unsigned i, const unsigned j) {
		pts2D r1 = particles_X[i];
		pts2D r2 = particles_Y[j];
		return getInteraction(r1, r2);
	}

	virtual kernel_dtype getInteraction(const pts2D r1, const pts2D r2) {
		cout << "virtual getInteraction" << endl;
		return 0.0;
	}
	*/
	Vec getRow(const int j, std::vector<int> col_indices) {
		int n_cols = col_indices.size();
		Vec row(n_cols);
    //#pragma omp parallel for
    for(int k = 0; k < n_cols; k++) {
        row(k) = this->getMatrixEntry(j, col_indices[k]);
    }
    return row;
  }

  Vec getCol(const int k, std::vector<int> row_indices) {
		int n_rows = row_indices.size();
    Vec col(n_rows);
    //#pragma omp parallel for
    for (int j=0; j<n_rows; ++j) {
			//cout << "j: " << j << "	row_indices[j]: " << row_indices[j] << "	k: " << k << endl;
      col(j) = this->getMatrixEntry(row_indices[j], k);
    }
    return col;
  }

  Mat getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
		int n_rows = row_indices.size();
		int n_cols = col_indices.size();
    Mat mat(n_rows, n_cols);
    //#pragma omp parallel for
    for (int j=0; j < n_rows; ++j) {
        //#pragma omp parallel for
        for (int k=0; k < n_cols; ++k) {
            mat(j,k) = this->getMatrixEntry(row_indices[j], col_indices[k]);
        }
    }
    return mat;
  }
};

class G_LowRank: public FMM_Matrix {
public:
	double tol_ACA_unused;
	double tol_robust;
	/*std::vector<pts2D> particles_X;
	std::vector<pts2D> particles_Y;*/
	std::vector<int> row_indices;
	std::vector<int> col_indices;
	/*
	kernel_dtype getMatrixEntry(const unsigned i, const unsigned j) {
		pts2D r1 = particles_X[row_indices[i]];
		pts2D r2 = particles_Y[col_indices[j]];
		return getInteraction(r1, r2);
	}

	virtual kernel_dtype getInteraction(const pts2D r1, const pts2D r2) {
		cout << "virtual getInteraction" << endl;
		return 0.0;
	}
	*/
	G_LowRank(int tol_pow, std::vector<pts2D> particles_X, std::vector<pts2D> particles_Y, std::vector<int> row_indices, std::vector<int> col_indices):
	FMM_Matrix(particles_X, particles_Y) {
			this->tol_ACA_unused = pow(10,-1.0*tol_pow);
			this->tol_robust = pow(10,-10);

			//this->particles_X = particles_X;
			//this->particles_Y = particles_Y;
			this->row_indices = row_indices;
			this->col_indices = col_indices;
	}
	/*
  Vec getRow(const int j) {
		int n_cols = col_indices.size();
		Vec row(n_cols);
    //#pragma omp parallel for
    for(int k = 0; k < n_cols; k++) {
        row(k) = this->getMatrixEntry(j, k);
    }
    return row;
  }

  Vec getCol(const int k) {
		int n_rows = row_indices.size();
    Vec col(n_rows);
    #pragma omp parallel for
    for (int j=0; j<n_rows; ++j) {
        col(j) = this->getMatrixEntry(j, k);
    }
    return col;
  }

  Mat getMatrix() {
			int n_rows = row_indices.size();
			int n_cols = col_indices.size();
      Mat mat(n_rows, n_cols);
      #pragma omp parallel for
      for (int j=0; j < n_rows; ++j) {
          #pragma omp parallel for
          for (int k=0; k < n_cols; ++k) {
              mat(j,k) = this->getMatrixEntry(j, k);
          }
      }
      return mat;
  }
	*/
	// void maxAbsVector(const Vec& v, const std::set<int>& allowed_indices,
	// 															dtype& max, int& index
	// 														 ) {
	// 	std::set<int>::iterator it;
	// 	index = *allowed_indices.begin();
	// 	max   = v(index);
	//
	// 	for(it = allowed_indices.begin(); it != allowed_indices.end(); it++) {
	// 			if(abs(v(*it)) > abs(max)) {
	// 					index   =   *it;
	// 					max     =   v(index);
	// 			}
	// 	}
	// }

	void maxAbsVector(const Vec& v, const std::set<int>& allowed_indices,
																	double max, int& index
																 ) {
			std::set<int>::iterator it;
			index = *allowed_indices.begin();
			max   = abs(v(index));

			for(it = allowed_indices.begin(); it != allowed_indices.end(); it++) {
					if(abs(v(*it)) > max) {
							index   =   *it;
							max     =   abs(v(index));
					}
			}
	}

	// void maxAbsVector(const Vec& v, const std::set<int>& allowed_indices,
	// 															dtype& max, int& index
	// 														 ) {
	// 	std::set<int>::iterator it;
	// 	index = *allowed_indices.begin();
	// 	max   = v(index);
	//
	// 	for(it = allowed_indices.begin(); it != allowed_indices.end(); it++) {
	// 			if(abs(v(*it)) > abs(max)) {
	// 					index   =   *it;
	// 					max     =   v(index);
	// 			}
	// 	}
	// }

	// void ACA_only_nodes(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, double tol_ACA, Mat &Ac, Mat &Ar) {
	// 	int col_index;
	// 	int row_index;
	// 	int N1 = row_indices.size();
	// 	int N2 = col_indices.size();
	// 	Vec row(N2), col(N1), v(N2), u(N1);
	// 	Vec row_temp, col_temp;
	// 	std::vector<Vec> Uvec;
	// 	std::vector<Vec> Vvec;
	// 	std::vector<Vec> AcVec;
	// 	std::vector<Vec> ArVec;
	// 	computed_rank = 0;
	// 	dtype max;
	// 	int min = N1;
	// 	if (N1 > N2) {
	// 		min = N2;
	// 	}
	// 	// Mat Ac = Mat(N1, min);
	// 	// Mat Ar = Mat(min, N2);
	//
	// 	std::set<int> remaining_row_ind;
	// 	std::set<int> remaining_col_ind;
	//
	// 	for(int l = 0; l < N1; l++) {
	// 			remaining_row_ind.insert(l);
	// 	}
	// 	for(int l= 0; l < N2; l++) {
	// 			remaining_col_ind.insert(l);
	// 	}
	// 	if (N1 < N2) {
	// 		int l_local;
	// 		for(l_local = 0; l_local < N1; l_local++) {
	// 			row_index=l_local;
	// 			row = getRow(row_indices[row_index], col_indices);
	// 			this->maxAbsVector(row, remaining_col_ind, max, col_index);
	// 			if(abs(row(int(col_index))) > tol_robust) {
	// 				break;
	// 			}
	// 		}
	// 		if (l_local == N1) {
	// 			Ac = Mat(N1,computed_rank);
	// 			Ar = Mat(computed_rank,N2);
	// 			// Ac_modified = Ac.block(0,0,N1,computed_rank);
	// 			// Ar_modified = Ar.block(0,0,computed_rank,N2);
	// 			return;
	// 		}
	// 		v=row/row(int(col_index));
	// 		row_bases.push_back(row_index);
	// 		col_bases.push_back(int(col_index));
	// 		col = getCol(col_indices[col_index], row_indices);
	// 		u	=	col;
	// 		Uvec.push_back(u);
	// 		Vvec.push_back(v);
	// 		// Ac.col(computed_rank) = col;
	// 		// Ar.row(computed_rank) = row;
	// 		AcVec.push_back(col);
	// 		ArVec.push_back(row);
	// 		remaining_col_ind.erase(col_index);
	// 		remaining_row_ind.erase(row_index);
	// 		computed_rank = 1;
	//
	// 		double normS	=	0.0;
	// 		double prev_normS = DBL_MAX;
	// 		this->maxAbsVector(col, remaining_row_ind, max, row_index);
	//
	// 		// while (abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
	// 		// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && (abs(prev_normS-normS) >= tol_ACA*prev_normS || abs(prev_normS-normS) >= tol_ACA) && computed_rank < min) {
	// 		while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
	// 			row_bases.push_back(int(row_index));
	// 			row = getRow(row_indices[row_index], col_indices);
	// 			row_temp = row;
	// 			for (int l=0; l<=computed_rank-1; ++l){
	// 				for (size_t s = 0; s < Vvec[l].size(); s++) {
	// 					row(s)	-=	Uvec[l](row_index)*Vvec[l](s);
	// 				}
	// 			}
	// 			this->maxAbsVector(row, remaining_col_ind, max, col_index);
	// 			col_bases.push_back(int(col_index));
	// 			for (size_t l = 0; l < row.size(); l++) {
	// 				v(l) = row(l)/row(int(col_index));
	// 			}
	// 			col = getCol(col_indices[col_index], row_indices);
	// 			col_temp = col;
	// 			for (int l=0; l<=computed_rank-1; ++l){
	// 				for (size_t s = 0; s < Uvec[l].size(); s++) {
	// 					col(s)	-=	Vvec[l](col_index)*Uvec[l](s);
	// 				}
	// 			}
	// 			u	=	col;
	// 			if(u.norm()< tol_robust || v.norm()< tol_robust) {
	// 				row_bases.pop_back();
	// 				col_bases.pop_back();
	// 				break;
	// 			}
	// 			Uvec.push_back(u);
	// 			Vvec.push_back(v);
	// 			// Ac.col(computed_rank) = col_temp;
	// 			// Ar.row(computed_rank) = row_temp;
	// 			AcVec.push_back(col_temp);
	// 			ArVec.push_back(row_temp);
	//
	// 			++computed_rank;
	// 			remaining_col_ind.erase(col_index);
	// 			remaining_row_ind.erase(row_index);
	// 			if (computed_rank != 2)
	// 				prev_normS = normS;
	// 			normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
	// 			for (int l=0; l<=computed_rank-1; ++l){
	// 				normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
	// 			}
	// 			normS	=	sqrt(normS);
	// 			this->maxAbsVector(col, remaining_row_ind, max, row_index);
	// 		}
	// 	}
	// 	else {
	// 		int l_local;
	// 		for(l_local = 0; l_local < N2; l_local++) {
	// 			col_index=l_local;
	// 			col = getCol(col_indices[col_index], row_indices);
	// 			this->maxAbsVector(col, remaining_row_ind, max, row_index);
	// 			if(abs(col(int(row_index))) > tol_robust) {
	// 				break;
	// 			}
	// 		}
	// 		if (l_local == N2) {
	// 			// Ac_modified = Ac.block(0,0,N1,computed_rank);
	// 			// Ar_modified = Ar.block(0,0,computed_rank,N2);
	// 			Ac = Mat(N1,computed_rank);
	// 			Ar = Mat(computed_rank,N2);
	// 			return;
	// 		}
	//
	// 		col_bases.push_back(col_index);
	// 		v=col/col(int(row_index));
	// 		row_bases.push_back(int(row_index));
	// 		row = getRow(row_indices[row_index], col_indices);
	// 		u	=	row;
	// 		Uvec.push_back(u);
	// 		Vvec.push_back(v);
	// 		// Ac.col(computed_rank) = col;
	// 		// Ar.row(computed_rank) = row;
	// 		AcVec.push_back(col);
	// 		ArVec.push_back(row);
	// 		computed_rank = 1;
	//
	// 		remaining_row_ind.erase(row_index);
	// 		remaining_col_ind.erase(col_index);
	//
	// 		double normS	=	0.0;
	// 		double prev_normS = DBL_MAX;
	// 		this->maxAbsVector(row, remaining_col_ind, max, col_index);
	//
	// 		// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && (abs(prev_normS-normS) >= tol_ACA*prev_normS || abs(prev_normS-normS) >= tol_ACA) && computed_rank < min) {
	// 		while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
	// 		// while (abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
	// 			col_bases.push_back(int(col_index));
	// 		  col = getCol(col_indices[col_index], row_indices);
	// 			col_temp = col;
	// 		  for (int l=0; l<=computed_rank-1; ++l){
	// 		    for (size_t s = 0; s < Vvec[l].size(); s++) {
	// 		      col(s)	-=	Uvec[l](col_index)*Vvec[l](s);
	// 		    }
	// 		  }
	// 		  this->maxAbsVector(col, remaining_row_ind, max, row_index);
	// 		  row_bases.push_back(int(row_index));
	// 		  for (size_t l = 0; l < col.size(); l++) {
	// 		    v(l) = col(l)/col(int(row_index));
	// 		  }
	// 		  row = getRow(row_indices[row_index], col_indices);
	// 			row_temp = row;
	// 		  for (int l=0; l<=computed_rank-1; ++l){
	// 		    for (size_t s = 0; s < Uvec[l].size(); s++) {
	// 		      row(s)	-=	Vvec[l](row_index)*Uvec[l](s);
	// 		    }
	// 		  }
	// 		  u	=	row;
	// 		  if(u.norm()< tol_robust || v.norm()< tol_robust) {
	// 		    col_bases.pop_back();
	// 		    row_bases.pop_back();
	// 		    break;
	// 		  }
	// 		  Uvec.push_back(u);
	// 		  Vvec.push_back(v);
	// 			// Ac.col(computed_rank) = col_temp;
	// 			// Ar.row(computed_rank) = row_temp;
	// 			AcVec.push_back(col_temp);
	// 			ArVec.push_back(row_temp);
	//
	// 		  ++computed_rank;
	// 			remaining_row_ind.erase(row_index);
	// 		  remaining_col_ind.erase(col_index);
	// 		  if (computed_rank != 2)
	// 		    prev_normS = normS;
	// 		  normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
	// 		  for (int l=0; l<=computed_rank-1; ++l){
	// 		    normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
	// 		  }
	// 		  normS	=	sqrt(normS);
	// 		  this->maxAbsVector(row, remaining_col_ind, max, col_index);
	// 		}
	// 	}
	// 	// Ac_modified = Ac.block(0,0,N1,computed_rank);
	// 	// Ar_modified = Ar.block(0,0,computed_rank,N2);
	// 	Ac = Mat(N1,computed_rank);
	// 	Ar = Mat(computed_rank,N2);
	// 	for (size_t i = 0; i < computed_rank; i++) {
	// 		Ac.col(i) = AcVec[i];
	// 		Ar.row(i) = ArVec[i];
	// 	}
	// }

	// void ACA_only_nodes(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, double tol_ACA, Mat &Ac, Mat &Ar) {
	// 	int col_index;
	// 	int row_index;
	// 	double tol_robust = 1e-10;
	// 	int N1 = row_indices.size();
	// 	int N2 = col_indices.size();
	// 	Vec row(N2), col(N1), v(N2), u(N1);
	// 	Vec row_temp, col_temp;
	// 	std::vector<Vec> Uvec;
	// 	std::vector<Vec> Vvec;
	// 	std::vector<Vec> AcVec;
	// 	std::vector<Vec> ArVec;
	// 	computed_rank = 0;
	// 	dtype max;
	// 	int min = N1;
	// 	if (N1 > N2) {
	// 		min = N2;
	// 	}
	// 	// Mat Ac = Mat(N1, min);
	// 	// Mat Ar = Mat(min, N2);
	//
	// 	std::set<int> remaining_row_ind;
	// 	std::set<int> remaining_col_ind;
	//
	// 	for(int l = 0; l < N1; l++) {
	// 			remaining_row_ind.insert(l);
	// 	}
	// 	for(int l= 0; l < N2; l++) {
	// 			remaining_col_ind.insert(l);
	// 	}
	// 	if (N1 < N2) {
	// 		int l_local;
	// 		for(l_local = 0; l_local < N1; l_local++) {
	// 			row_index=l_local;
	// 			row = getRow(row_indices[row_index], col_indices);
	// 			this->maxAbsVector(row, remaining_col_ind, max, col_index);
	// 			if(abs(row(int(col_index))) > tol_robust) {
	// 				break;
	// 			}
	// 		}
	// 		if (l_local == N1) {
	// 			Ac = Mat(N1,computed_rank);
	// 			Ar = Mat(computed_rank,N2);
	// 			// Ac_modified = Ac.block(0,0,N1,computed_rank);
	// 			// Ar_modified = Ar.block(0,0,computed_rank,N2);
	// 			return;
	// 		}
	// 		v=row/row(int(col_index));
	// 		row_bases.push_back(row_index);
	// 		col_bases.push_back(int(col_index));
	// 		col = getCol(col_indices[col_index], row_indices);
	// 		u	=	col;
	// 		Uvec.push_back(u);
	// 		Vvec.push_back(v);
	// 		// Ac.col(computed_rank) = col;
	// 		// Ar.row(computed_rank) = row;
	// 		AcVec.push_back(col);
	// 		ArVec.push_back(row);
	// 		remaining_col_ind.erase(col_index);
	// 		remaining_row_ind.erase(row_index);
	// 		computed_rank = 1;
	//
	// 		double normS	=	0.0;
	// 		double prev_normS = DBL_MAX;
	// 		this->maxAbsVector(col, remaining_row_ind, max, row_index);
	//
	// 		// while (abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
	// 		// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && (abs(prev_normS-normS) >= tol_ACA*prev_normS || abs(prev_normS-normS) >= tol_ACA) && computed_rank < min) {
	// 		while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
	// 			row_bases.push_back(int(row_index));
	// 			row = getRow(row_indices[row_index], col_indices);
	// 			row_temp = row;
	// 			for (int l=0; l<=computed_rank-1; ++l){
	// 				for (size_t s = 0; s < Vvec[l].size(); s++) {
	// 					row(s)	-=	Uvec[l](row_index)*Vvec[l](s);
	// 				}
	// 			}
	// 			this->maxAbsVector(row, remaining_col_ind, max, col_index);
	// 			col_bases.push_back(int(col_index));
	// 			for (size_t l = 0; l < row.size(); l++) {
	// 				v(l) = row(l)/row(int(col_index));
	// 			}
	// 			col = getCol(col_indices[col_index], row_indices);
	// 			col_temp = col;
	// 			for (int l=0; l<=computed_rank-1; ++l){
	// 				for (size_t s = 0; s < Uvec[l].size(); s++) {
	// 					col(s)	-=	Vvec[l](col_index)*Uvec[l](s);
	// 				}
	// 			}
	// 			u	=	col;
	// 			if(u.norm()< tol_robust || v.norm()< tol_robust) {
	// 				row_bases.pop_back();
	// 				col_bases.pop_back();
	// 				break;
	// 			}
	// 			Uvec.push_back(u);
	// 			Vvec.push_back(v);
	// 			// Ac.col(computed_rank) = col_temp;
	// 			// Ar.row(computed_rank) = row_temp;
	// 			AcVec.push_back(col_temp);
	// 			ArVec.push_back(row_temp);
	//
	// 			++computed_rank;
	// 			remaining_col_ind.erase(col_index);
	// 			remaining_row_ind.erase(row_index);
	// 			if (computed_rank != 2)
	// 				prev_normS = normS;
	// 			normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
	// 			for (int l=0; l<=computed_rank-1; ++l){
	// 				normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
	// 			}
	// 			normS	=	sqrt(normS);
	// 			this->maxAbsVector(col, remaining_row_ind, max, row_index);
	// 		}
	// 	}
	// 	else {
	// 		int l_local;
	// 		for(l_local = 0; l_local < N2; l_local++) {
	// 			col_index=l_local;
	// 			col = getCol(col_indices[col_index], row_indices);
	// 			this->maxAbsVector(col, remaining_row_ind, max, row_index);
	// 			if(abs(col(int(row_index))) > tol_robust) {
	// 				break;
	// 			}
	// 		}
	// 		if (l_local == N2) {
	// 			// Ac_modified = Ac.block(0,0,N1,computed_rank);
	// 			// Ar_modified = Ar.block(0,0,computed_rank,N2);
	// 			Ac = Mat(N1,computed_rank);
	// 			Ar = Mat(computed_rank,N2);
	// 			return;
	// 		}
	//
	// 		col_bases.push_back(col_index);
	// 		v=col/col(int(row_index));
	// 		row_bases.push_back(int(row_index));
	// 		row = getRow(row_indices[row_index], col_indices);
	// 		u	=	row;
	// 		Uvec.push_back(u);
	// 		Vvec.push_back(v);
	// 		// Ac.col(computed_rank) = col;
	// 		// Ar.row(computed_rank) = row;
	// 		AcVec.push_back(col);
	// 		ArVec.push_back(row);
	// 		computed_rank = 1;
	//
	// 		remaining_row_ind.erase(row_index);
	// 		remaining_col_ind.erase(col_index);
	//
	// 		double normS	=	0.0;
	// 		double prev_normS = DBL_MAX;
	// 		this->maxAbsVector(row, remaining_col_ind, max, col_index);
	//
	// 		// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && (abs(prev_normS-normS) >= tol_ACA*prev_normS || abs(prev_normS-normS) >= tol_ACA) && computed_rank < min) {
	// 		while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
	// 		// while (abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
	// 			col_bases.push_back(int(col_index));
	// 		  col = getCol(col_indices[col_index], row_indices);
	// 			col_temp = col;
	// 		  for (int l=0; l<=computed_rank-1; ++l){
	// 		    for (size_t s = 0; s < Vvec[l].size(); s++) {
	// 		      col(s)	-=	Uvec[l](col_index)*Vvec[l](s);
	// 		    }
	// 		  }
	// 		  this->maxAbsVector(col, remaining_row_ind, max, row_index);
	// 		  row_bases.push_back(int(row_index));
	// 		  for (size_t l = 0; l < col.size(); l++) {
	// 		    v(l) = col(l)/col(int(row_index));
	// 		  }
	// 		  row = getRow(row_indices[row_index], col_indices);
	// 			row_temp = row;
	// 		  for (int l=0; l<=computed_rank-1; ++l){
	// 		    for (size_t s = 0; s < Uvec[l].size(); s++) {
	// 		      row(s)	-=	Vvec[l](row_index)*Uvec[l](s);
	// 		    }
	// 		  }
	// 		  u	=	row;
	// 		  if(u.norm()< tol_robust || v.norm()< tol_robust) {
	// 		    col_bases.pop_back();
	// 		    row_bases.pop_back();
	// 		    break;
	// 		  }
	// 		  Uvec.push_back(u);
	// 		  Vvec.push_back(v);
	// 			// Ac.col(computed_rank) = col_temp;
	// 			// Ar.row(computed_rank) = row_temp;
	// 			AcVec.push_back(col_temp);
	// 			ArVec.push_back(row_temp);
	//
	// 		  ++computed_rank;
	// 			remaining_row_ind.erase(row_index);
	// 		  remaining_col_ind.erase(col_index);
	// 		  if (computed_rank != 2)
	// 		    prev_normS = normS;
	// 		  normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
	// 		  for (int l=0; l<=computed_rank-1; ++l){
	// 		    normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
	// 		  }
	// 		  normS	=	sqrt(normS);
	// 		  this->maxAbsVector(row, remaining_col_ind, max, col_index);
	// 		}
	// 	}
	// 	// Ac_modified = Ac.block(0,0,N1,computed_rank);
	// 	// Ar_modified = Ar.block(0,0,computed_rank,N2);
	// 	Ac = Mat(N1,computed_rank);
	// 	Ar = Mat(computed_rank,N2);
	// 	for (size_t i = 0; i < computed_rank; i++) {
	// 		Ac.col(i) = AcVec[i];
	// 		Ar.row(i) = ArVec[i];
	// 	}
	// }
	void ACA_only_nodes(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, double epsilon, Mat &Ac, Mat &Ar) {
		int col_index;
		int row_index;
		int N1 = row_indices.size();
		int N2 = col_indices.size();
		Vec row(N2), col(N1), v(N2), u(N1);
		Vec row_temp, col_temp;
		std::vector<Vec> Uvec;
		std::vector<Vec> Vvec;
		computed_rank = 0;
		double max;
		int min = N1;
		if (N1 > N2) {
			min = N2;
		}
		Ac = Mat(N1, min);
		Ar = Mat(min, N2);

		std::set<int> remaining_row_ind;
		std::set<int> remaining_col_ind;

		for(int l = 0; l < N1; l++) {
				remaining_row_ind.insert(l);
		}
		for(int l= 0; l < N2; l++) {
				remaining_col_ind.insert(l);
		}
		if (N1 < N2) {
			// row_index=0;
			// row_bases.push_back(row_index);
			// row = getRow(row_indices[row_index], col_indices);
			// this->maxAbsVector(row, remaining_col_ind, max, col_index);
			// v=row/row(int(col_index));
			// col_bases.push_back(int(col_index));
			// col = getCol(col_indices[col_index], row_indices);
			// u	=	col;
			// Uvec.push_back(u);
			// Vvec.push_back(v);
			// Ac.col(computed_rank) = col;
			// Ar.row(computed_rank) = row;
			// remaining_col_ind.erase(col_index);
			// remaining_row_ind.erase(row_index);
			// computed_rank = 1;

			int l_local;
			for(l_local = 0; l_local < N1; l_local++) {
				row_index=l_local;
				row = getRow(row_indices[row_index], col_indices);
				this->maxAbsVector(row, remaining_col_ind, max, col_index);
				if(abs(row(int(col_index))) > pow(10,-1.0*16)) {
					break;
				}
			}
			if (l_local == N1) {
				Ac = Ac.block(0,0,N1,computed_rank);
				Ar = Ar.block(0,0,computed_rank,N2);
				return;
			}
			v=row/row(int(col_index));
			row_bases.push_back(row_index);
			col_bases.push_back(int(col_index));
			col = getCol(col_indices[col_index], row_indices);
			u	=	col;
			Uvec.push_back(u);
			Vvec.push_back(v);
			Ac.col(computed_rank) = col;
			Ar.row(computed_rank) = row;
			remaining_col_ind.erase(col_index);
			remaining_row_ind.erase(row_index);
			computed_rank = 1;

			double normS	=	0.0;
			double prev_normS = DBL_MAX;//1000.0;
			this->maxAbsVector(col, remaining_row_ind, max, row_index);

			// std::cout << abs(row(int(col_index))) << ";	" << abs(col(int(row_index))) << std::endl;
			// while(u.norm()*v.norm() > epsilon*normS && computed_rank < min) {
			// while (abs(row(int(col_index))) > pow(10,-1.0*25) && abs(col(int(row_index))) > pow(10,-1.0*25) && abs(prev_normS-normS) >= epsilon*prev_normS && computed_rank < min) {
			// while (abs(prev_normS-normS) >= epsilon && computed_rank < min) {
			while (abs(prev_normS-normS) >= epsilon*prev_normS && computed_rank < min) {
				row_bases.push_back(int(row_index));
				row = getRow(row_indices[row_index], col_indices);
				row_temp = row;
				for (int l=0; l<=computed_rank-1; ++l){
					for (size_t s = 0; s < Vvec[l].size(); s++) {
						row(s)	-=	Uvec[l](row_index)*Vvec[l](s);
					}
				}
				this->maxAbsVector(row, remaining_col_ind, max, col_index);
				col_bases.push_back(int(col_index));
				for (size_t l = 0; l < row.size(); l++) {
					v(l) = row(l)/row(int(col_index));
				}
				col = getCol(col_indices[col_index], row_indices);
				col_temp = col;
				for (int l=0; l<=computed_rank-1; ++l){
					for (size_t s = 0; s < Uvec[l].size(); s++) {
						col(s)	-=	Vvec[l](col_index)*Uvec[l](s);
					}
				}
				u	=	col;
				// if(u.norm()< epsilon*pow(10,-1.0*16) || v.norm()< epsilon*pow(10,-1.0*16)) {
				// 	row_bases.pop_back();
				// 	col_bases.pop_back();
				// 	// std::cout << "here" << std::endl;
				// 	break;
				// }
				Uvec.push_back(u);
				Vvec.push_back(v);
				Ac.col(computed_rank) = col_temp;
				Ar.row(computed_rank) = row_temp;

				++computed_rank;
				remaining_col_ind.erase(col_index);
				remaining_row_ind.erase(row_index);
				if (computed_rank != 2)
					prev_normS = normS;
				normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
				for (int l=0; l<=computed_rank-1; ++l){
					normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
				}
				normS	=	sqrt(normS);
				this->maxAbsVector(col, remaining_row_ind, max, row_index);
			}
			// exit(0);
		}
		else {
			// col_index=0;
			// col_bases.push_back(col_index);
			// //cout << "row_indices.size(): " << row_indices.size() << endl;
			// //for (size_t d = 0; d < row_indices.size(); d++) {
			// //	cout << row_indices[d] << endl;
			// //}
			// //cout << "col_indices[col_index]: " << col_indices[col_index] << endl;
			// col = getCol(col_indices[col_index], row_indices);
			// //cout << "in ACA" << endl;
			//
			// this->maxAbsVector(col, remaining_row_ind, max, row_index);
			// v=col/col(int(row_index));
			// row_bases.push_back(int(row_index));
			// row = getRow(row_indices[row_index], col_indices);
			// u	=	row;
			// Uvec.push_back(u);
			// Vvec.push_back(v);
			// Ac.col(computed_rank) = col;
			// Ar.row(computed_rank) = row;
			// computed_rank = 1;

			int l_local;
			for(l_local = 0; l_local < N2; l_local++) {
				col_index=l_local;
				col = getCol(col_indices[col_index], row_indices);
				this->maxAbsVector(col, remaining_row_ind, max, row_index);
				if(abs(col(int(row_index))) > pow(10,-1.0*16)) {
					break;
				}
			}
			if (l_local == N2) {
				Ac = Ac.block(0,0,N1,computed_rank);
				Ar = Ar.block(0,0,computed_rank,N2);
				return;
			}

			// col_index=0;
			col_bases.push_back(col_index);
			// col = getCol(col_indices[col_index], row_indices);
			// this->maxAbsVector(col, remaining_row_ind, max, row_index);
			v=col/col(int(row_index));
			row_bases.push_back(int(row_index));
			row = getRow(row_indices[row_index], col_indices);
			u	=	row;
			Uvec.push_back(u);
			Vvec.push_back(v);
			Ac.col(computed_rank) = col;
			Ar.row(computed_rank) = row;
			computed_rank = 1;

			remaining_row_ind.erase(row_index);
			remaining_col_ind.erase(col_index);

			double normS	=	0.0;
			double prev_normS = DBL_MAX;//1000.0;
			this->maxAbsVector(row, remaining_col_ind, max, col_index);

			while (abs(prev_normS-normS) >= epsilon*prev_normS && computed_rank < min) {
			// while (abs(prev_normS-normS) >= epsilon && computed_rank < min) {
			// while(u.norm()*v.norm() > epsilon*normS &&  computed_rank < min) {
			// while (abs(row(int(col_index))) > pow(10,-1.0*25) && abs(col(int(row_index))) > pow(10,-1.0*25) && abs(prev_normS-normS) >= epsilon*prev_normS && computed_rank < min) {
				col_bases.push_back(int(col_index));
			  col = getCol(col_indices[col_index], row_indices);
				col_temp = col;
			  for (int l=0; l<=computed_rank-1; ++l){
			    for (size_t s = 0; s < Vvec[l].size(); s++) {
			      col(s)	-=	Uvec[l](col_index)*Vvec[l](s);
			    }
			  }
			  this->maxAbsVector(col, remaining_row_ind, max, row_index);
			  row_bases.push_back(int(row_index));
			  for (size_t l = 0; l < col.size(); l++) {
			    v(l) = col(l)/col(int(row_index));
			  }
			  row = getRow(row_indices[row_index], col_indices);
				row_temp = row;
			  for (int l=0; l<=computed_rank-1; ++l){
			    for (size_t s = 0; s < Uvec[l].size(); s++) {
			      row(s)	-=	Vvec[l](row_index)*Uvec[l](s);
			    }
			  }
			  u	=	row;
			  // if(u.norm()< epsilon*pow(10,-1.0*16) || v.norm()< epsilon*pow(10,-1.0*16)) {
			  //   col_bases.pop_back();
			  //   row_bases.pop_back();
			  //   break;
			  // }
			  Uvec.push_back(u);
			  Vvec.push_back(v);
				Ac.col(computed_rank) = col_temp;
				Ar.row(computed_rank) = row_temp;

			  ++computed_rank;
				remaining_row_ind.erase(row_index);
			  remaining_col_ind.erase(col_index);
			  if (computed_rank != 2)
			    prev_normS = normS;
			  normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
			  for (int l=0; l<=computed_rank-1; ++l){
			    normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
			  }
			  normS	=	sqrt(normS);
			  this->maxAbsVector(row, remaining_col_ind, max, col_index);
			}
		}
		Ac = Ac.block(0,0,N1,computed_rank);
		Ar = Ar.block(0,0,computed_rank,N2);
		// if(Ac.norm() < epsilon || Ar.norm() < epsilon) {
		// 	cout << "Ac.n(): " << Ac.norm() << "	Ar.n(): " << Ar.norm() << endl;
		// }
	}

};
