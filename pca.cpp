#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <assert.h>

#include <stdio.h>
#include <stdlib.h>

#include "tnt_array1d.h"
#include "tnt_array2d.h"

#include "jama_eig.h"


using namespace std;
using namespace TNT;
using namespace JAMA;

namespace PCA {
	bool debug = false;
	
	template < class T>
	void convert_from_string(T& value, const string& s)
	{
		stringstream ss(s);
		ss >> value;
	}
	
	void load_data_from_file(Array2D<double>& d, char*& file_path) {
		ifstream in(file_path);
        string line;
        int r = 0;
		
        if (in.is_open()) {
            while (in.good()) {
				int col = 0;
                getline(in, line);
                if (line.empty()) continue;
                
				size_t start_pos = 0;
				char space = ' ';
				while (true) {
					size_t pos = line.find(space, start_pos);
					string data(line.substr(start_pos, pos - start_pos));
					if (!data.empty()) {
						double v = 0;
						convert_from_string(v, data);
						
						if (debug)
							cout << "value: " << v << endl;
						d[r][col] = v;
					}
              
					if ((int)pos != -1) {
						start_pos = pos + 1;
					}
					else {
						break;
					}
					col += 1;
				}
				r += 1;
            }
            in.close();
        }
	}
	
	void adjust_data(Array2D<double>& d, Array1D<double>& means) {
	   for (int i=0; i<d.dim2(); ++i) { 
		   double mean = 0;
		   for (int j=0; j<d.dim1(); ++j) {
			   mean += d[j][i];
		   }

		   mean /= d.dim1();

		   // store the mean
		   means[i] = mean;

		   // subtract the mean
		   for (int j=0; j<d.dim1(); ++j) {
			   d[j][i] -= mean;
		   }
	   }
	}

	double compute_covariance(const Array2D<double>& d, int i, int j) {
	   double cov = 0;
	   for (int k=0; k<d.dim1(); ++k) {
		   cov += d[k][i] * d[k][j];
	   }

	   return cov / (d.dim1() - 1);
	}

	void compute_covariance_matrix(const Array2D<double> & d, Array2D<double> & covar_matrix) {
		int dim = d.dim2();
		assert(dim == covar_matrix.dim1());
		assert(dim == covar_matrix.dim2());
		for (int i=0; i<dim; ++i) {
			for (int j=i; j<dim; ++j) {
				covar_matrix[i][j] = compute_covariance(d, i, j);
			}
		}


		// fill the Left triangular matrix
		for (int i=1; i<dim; i++) {
			for (int j=0; j<i; ++j) {
				covar_matrix[i][j] = covar_matrix[j][i];
			}
		}

	}

	// Calculate the eigenvectors and eigenvalues of the covariance
	// matrix
	void eigen(const Array2D<double> & covar_matrix, Array2D<double>& eigenvector, Array2D<double>& eigenvalue) {
		Eigenvalue<double> eig(covar_matrix);
		eig.getV(eigenvector);
		eig.getD(eigenvalue);
	}


	void transpose(const Array2D<double>& src, Array2D<double>& target) {
		for (int i=0; i<src.dim1(); ++i) {
			for (int j=0; j<src.dim2(); ++j) {
				target[j][i] = src[i][j];
			}
		}
	}

	// z = x * y
	void multiply(const Array2D<double>& x, const Array2D<double>& y, Array2D<double>& z) {
		assert(x.dim2() == y.dim1());
		for (int i=0; i<x.dim1(); ++i) {
			for (int j=0; j<y.dim2(); ++j) {
				double sum = 0;
				int d = y.dim1();
				for (int k=0; k<d; k++) {
					sum += x[i][k] * y[k][j];
				}
				z[i][j] = sum;
			}
		}
	}
}

int main(int argc, char* argv[]) {
    if (argc != 2) {        
        cout << "Usage: " << argv[0] << " data_file" << endl;        
        return -1;   
    }
	using namespace PCA;
	
    const int row = 10;
    const int col = 2;

    Array2D<double> d(row, col);
    load_data_from_file(d, argv[1]);
    Array1D<double> means(col);
    adjust_data(d, means);

    Array2D<double> covar_matrix(col, col);
    compute_covariance_matrix(d, covar_matrix);

    int dim = covar_matrix.dim1();
    // get the eigenvectors
    Array2D<double> eigenvector(dim, dim);

    // get the eigenvalues
    Array2D<double> eigenvalue(dim, dim);
    eigen(covar_matrix, eigenvector, eigenvalue);
    cout << "The eigenvectors:" << endl;
    cout << eigenvector << endl;

    cout << "The eigenvalues:" << endl;
    cout << eigenvalue << endl;

    // restore the old data
    // final_data = RowFeatureVector * RowDataAdjust
    Array2D<double> final_data(row, col);
    Array2D<double> transpose_data(col, row);
    transpose(d, transpose_data);
    multiply(eigenvector, transpose_data, final_data);

    cout << "the final data" << endl;
    cout << final_data << endl;

    return 0;
}
