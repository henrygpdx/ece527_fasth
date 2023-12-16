// SVD.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"

// JordanGauss.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <stdio.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <cmath>

#include "xparameters.h"	/* SDK generated parameters */
#include "xil_printf.h"
#include "xil_cache.h"
#include "xplatform_info.h"

// Timer
#include "xtmrctr.h"
#define TMRCTR_DEVICE_ID    XPAR_AXI_TIMER_0_DEVICE_ID

//#include <boost/multiprecision/cpp_bin_float.hpp>
//#include <boost/math/special_functions/gamma.hpp>

//using namespace boost::multiprecision;

//typedef number<backends::cpp_bin_float<4096, backends::digit_base_2, void, boost::int16_t, -16382, 16383>, et_off>  cpp_bin_float_4096;

// (2, 50), (3, 20)
/*#define MATRIX_SIZE 4
#define RESIZE_FACTOR 20
#define ACCURACY 1e-5

using namespace std;



std::vector<std::vector<float>> matrix_i;



template<typename Arg = float>
void matrix_transpose(std::vector<std::vector<Arg>> matrix1,
	std::vector<std::vector<Arg>>& matrix2)
{
	matrix2.resize(matrix1.size());
	for (std::size_t row = 0; row < matrix1.size(); row++)
	{
		matrix2[row].resize(matrix1[row].size());
		for (std::size_t col = 0; col < matrix1[row].size(); col++)
			matrix2[row][col] = matrix1[col][row];
	}
}

template<typename Arg = float>
void matrix_by_matrix(std::vector<std::vector<Arg>> matrix1,
	std::vector<std::vector<Arg>>& matrix2, std::vector<std::vector<Arg>>& matrix3)
{
	matrix3.resize(matrix1.size());
	for (std::size_t row = 0; row < matrix1.size(); row++)
	{
		matrix3[row].resize(matrix1[row].size());
		for (std::size_t col = 0; col < matrix1[row].size(); col++)
		{
			matrix3[row][col] = 0.00;
			for (std::size_t k = 0; k < matrix1[row].size(); k++)
				matrix3[row][col] += matrix1[row][k] * matrix2[k][col];
		}
	}
}

template<typename Arg = float>
void get_hermitian_matrix(std::vector<Arg> eigenvector,
	std::vector<std::vector<Arg>>& h_matrix)
{
	h_matrix.resize(eigenvector.size());
	for (std::size_t row = 0; row < eigenvector.size(); row++)
		h_matrix[row].resize(eigenvector.size());

	h_matrix[0][0] = 1.0 / eigenvector[0];
	for (std::size_t row = 1; row < eigenvector.size(); row++)
		h_matrix[row][0] = -eigenvector[row] / eigenvector[0];

	for (std::size_t row = 1; row < eigenvector.size(); row++)
		h_matrix[row][row] = 1;
}

template<typename Arg = float>
void get_hermitian_matrix_inverse(std::vector<Arg> eigenvector,
	std::vector<std::vector<Arg>>& ih_matrix)
{
	ih_matrix.resize(eigenvector.size());
	for (std::size_t row = 0; row < eigenvector.size(); row++)
		ih_matrix[row].resize(eigenvector.size());

	ih_matrix[0][0] = eigenvector[0];
	for (std::size_t row = 1; row < eigenvector.size(); row++)
		ih_matrix[row][0] = -eigenvector[row];

	for (std::size_t row = 1; row < eigenvector.size(); row++)
		ih_matrix[row][row] = 1;
}

template<typename Arg = float>
void jordan_gaussian_transform(
	std::vector<std::vector<Arg>> matrix, std::vector<Arg>& eigenvector)
{
	const Arg eps = 10e-6; bool eigenv_found = false;
	for (std::size_t s = 0; s < matrix.size() - 1 && !eigenv_found; s++)
	{
		std::size_t col = s; Arg alpha = matrix[s][s];
		while (col < matrix[s].size() && alpha != 0 && alpha != 1)
			matrix[s][col++] /= alpha;

		for (std::size_t col = s; col < matrix[s].size() && !alpha; col++)
			std::swap(matrix[s][col], matrix[s + 1][col]);

		for (std::size_t row = 0; row < matrix.size(); row++)
		{
			Arg gamma = matrix[row][s];
			for (std::size_t col = s; col < matrix[row].size() && row != s; col++)
				matrix[row][col] = matrix[row][col] - matrix[s][col] * gamma;
		}
	}

	std::size_t row = 0;
	while (row < matrix.size())
		eigenvector.push_back(-matrix[row++][matrix.size() - 1]);

	eigenvector[eigenvector.size() - 1] = 1;
}

template<typename Arg = float>
void get_inverse_diagonal_matrix(std::vector<std::vector<Arg>> matrix,
	std::vector<std::vector<Arg>>& inv_matrix)
{
	inv_matrix.resize(matrix.size());
	for (std::size_t index = 0; index < matrix.size(); index++)
	{
		inv_matrix[index].resize(matrix[index].size());
		inv_matrix[index][index] = 1.0 / matrix[index][index];
	}
}

template<typename Arg = float>
void get_reduced_matrix(std::vector<std::vector<Arg>> matrix,
	std::vector<std::vector<Arg>>& r_matrix, std::size_t new_size)
{
	if (new_size > 1)
	{
		r_matrix.resize(new_size);
		std::size_t index_d = matrix.size() - new_size;
		std::size_t row = index_d, row_n = 0;
		while (row < matrix.size())
		{
			r_matrix[row_n].resize(new_size);
			std::size_t col = index_d, col_n = 0;
			while (col < matrix.size())
				r_matrix[row_n][col_n++] = matrix[row][col++];

			row++; row_n++;
		}
	}

	else if (new_size == 1)
	{
		r_matrix.resize(new_size);
		r_matrix[0].resize(new_size);
		r_matrix[0][0] = matrix[1][1];
	}
}

template<typename Arg = float>
void generate_matrix(std::vector<std::vector<float>>& \
	matrix, std::size_t rows, std::size_t cols)
{
	//std::srand((unsigned int)std::time(nullptr));
	matrix.resize(rows);
	for (std::size_t row = 0; row < matrix.size(); row++)
	{
		matrix[row].resize(cols);
		for (std::size_t col = 0; col < matrix[row].size(); col++)
			matrix[row][col] = std::rand() % 10;
	}
}

template<typename Arg = float>
void print_matrix(std::vector<std::vector<float>>	matrix)
{
	for (std::size_t row = 0; row < matrix.size(); row++)
	{
		for (std::size_t col = 0; col < matrix[row].size(); col++) {
			std::cout << std::setprecision(5) << std::fixed << matrix[row][col] << " ";
			//int whole = matrix[row][col];
			//int thousandths = (matrix[row][col]-whole)*1000;
			//xil_printf("%d.%5d ", whole, thousandths);
		}

		std::cout << "\n";
	}

	std::cout << "\n";
}
template<typename Arg = float>
void compute_evd(std::vector<std::vector<Arg>> matrix,
	std::vector<Arg>& eigenvalues, std::vector<std::vector<Arg>>& eigenvectors, std::size_t eig_count)
{
	std::cout << "entering compute_evd...\n";
	std::size_t m_size = matrix.size();
	std::vector<Arg> vec; vec.resize(m_size);
	std::generate(vec.begin(), vec.end(), []() {
		return rand() % 100;
	});

	if (eigenvalues.size() == 0 && eigenvectors.size() == 0)
	{
		eigenvalues.resize(m_size);
		eigenvectors.resize(eigenvalues.size());
		matrix_i = matrix;
	}

	std::vector<std::vector<Arg>> m;
	m.resize(m_size); // correct number of columns
	for (std::size_t row = 0; row < m_size; row++)
		m[row].resize(RESIZE_FACTOR);

	Arg lambda_old = 0;

	std::size_t index = 0; bool is_eval = false;
	std::cout << "entering while loop...\n";
	while (is_eval == false)
	{
		std::cout << "Beginning while loop...";
		for (std::size_t row = 0; row < m_size && (index % RESIZE_FACTOR) == 0; row++)
			m[row].resize(m[row].size() + RESIZE_FACTOR);

		for (std::size_t row = 0; row < m_size; row++)
		{
			m[row][index] = 0;
			for (std::size_t col = 0; col < m_size; col++)
				m[row][index] += matrix[row][col] * vec[col];
		}

		for (std::size_t col = 0; col < m_size; col++)
			vec[col] = m[col][index];

		if (index > 0)
		{
			Arg lambda = ((index - 1) > 0) ? \
				(m[0][index] / m[0][index - 1]) : m[0][index];
			//is_eval = (boost::multiprecision::fabs(lambda - lambda_old) < 10e-10) ? true : false; //10e-15
			is_eval = (std::fabs(lambda - lambda_old) < ACCURACY) ? true : false;

			eigenvalues[eig_count] = lambda; lambda_old = lambda;
		}
		std::cout << "Incrementing index\n";
		index++;
	}

	std::vector<std::vector<Arg>> matrix_new;
	std::cout << "entering if statements...\n";
	if (m_size > 1)
	{
		std::vector<std::vector<Arg>> matrix_target;
		matrix_target.resize(m_size);

		for (std::size_t row = 0; row < m_size; row++)
		{
			matrix_target[row].resize(m_size);
			for (std::size_t col = 0; col < m_size; col++)
				matrix_target[row][col] = (row == col) ? \
				(matrix[row][col] - eigenvalues[eig_count]) : matrix[row][col];
		}

		std::vector<Arg> eigenvector;
		jordan_gaussian_transform(matrix_target, eigenvector);

		std::vector<std::vector<Arg>> hermitian_matrix;
		get_hermitian_matrix(eigenvector, hermitian_matrix);

		std::vector<std::vector<Arg>> ha_matrix_product;
		matrix_by_matrix(hermitian_matrix, matrix, ha_matrix_product);

		std::vector<std::vector<Arg>> inverse_hermitian_matrix;
		get_hermitian_matrix_inverse(eigenvector, inverse_hermitian_matrix);

		std::vector<std::vector<Arg>> iha_matrix_product;
		matrix_by_matrix(ha_matrix_product, inverse_hermitian_matrix, iha_matrix_product);

		get_reduced_matrix(iha_matrix_product, matrix_new, m_size - 1);
	}

	if (m_size <= 1)
	{
		for (std::size_t index = 0; index < eigenvalues.size(); index++)
		{
			Arg lambda = eigenvalues[index];
			std::vector<std::vector<Arg>> matrix_target;
			matrix_target.resize(matrix_i.size());

			for (std::size_t row = 0; row < matrix_i.size(); row++)
			{
				matrix_target[row].resize(matrix_i.size());
				for (std::size_t col = 0; col < matrix_i.size(); col++)
					matrix_target[row][col] = (row == col) ? \
					(matrix_i[row][col] - lambda) : matrix_i[row][col];
			}

			eigenvectors.resize(matrix_i.size());
			jordan_gaussian_transform(matrix_target, eigenvectors[index]);

			Arg eigsum_sq = 0;
			for (std::size_t v = 0; v < eigenvectors[index].size(); v++)
				eigsum_sq += std::pow(eigenvectors[index][v], 2.0);

			for (std::size_t v = 0; v < eigenvectors[index].size(); v++)
				eigenvectors[index][v] /= std::sqrt(eigsum_sq);

			eigenvalues[index] = std::sqrt(eigenvalues[index]);
		}

		return;
	}

	compute_evd(matrix_new, eigenvalues, eigenvectors, eig_count + 1);

	return;
}
template<typename Arg = float>
void svd(std::vector<std::vector<Arg>> matrix, std::vector<std::vector<Arg>>& s,
	std::vector<std::vector<Arg>>& u, std::vector<std::vector<Arg>>& v)
{
	std::vector<std::vector<Arg>> matrix_t;
	matrix_transpose(matrix, matrix_t);

	std::vector<std::vector<Arg>> matrix_product1;
	matrix_by_matrix(matrix, matrix_t, matrix_product1);

	std::vector<std::vector<Arg>> matrix_product2;
	matrix_by_matrix(matrix_t, matrix, matrix_product2);

	std::vector<std::vector<Arg>> u_1;
	std::vector<std::vector<Arg>> v_1;

	std::vector<Arg> eigenvalues;
	compute_evd(matrix_product2, eigenvalues, v_1, 0);

	matrix_transpose(v_1, v);

	s.resize(matrix.size());
	for (std::size_t index = 0; index < eigenvalues.size(); index++)
	{
		s[index].resize(eigenvalues.size());
		s[index][index] = eigenvalues[index];
	}

	std::vector<std::vector<Arg>> s_inverse;
	get_inverse_diagonal_matrix(s, s_inverse);

	std::vector<std::vector<Arg>> av_matrix;
	matrix_by_matrix(matrix, v, av_matrix);
	matrix_by_matrix(av_matrix, s_inverse, u);
}

int main(int argc, char* argv[])
{
	std::size_t matrix_size = MATRIX_SIZE;
	std::vector<std::vector<float>> u, v;
	std::vector<std::vector<float>> matrix, s;
	std::cout << "Singular Value Decomposition (SVD):\nResize Factor: " << RESIZE_FACTOR << "\nAccuracy: " << ACCURACY << "\n";

	//std::cout << "Enter size of matrix N = (50x50 max): "; std::cin >> matrix_size;

	if (matrix_size <= 50)
	{
		generate_matrix(matrix, matrix_size, matrix_size);
		//matrix[0][0] = 5;
		//matrix[0][1] = 5;
		//matrix[1][0] = -1;
		//matrix[1][1] = 7;
		std::cout << "\nA = \n"; print_matrix(matrix);
		std::cout << "Finished printing\n";
		svd(matrix, s, u, v);

		std::cout << "\nS = \n"; print_matrix(s);
		std::cout << "\nU = \n"; print_matrix(u);
		std::cout << "\nV = \n"; print_matrix(v);
	}

	else std::cout << "Wrong matrix size... (matrix decomposition recommended)";

	//std::cin.get();
	//std::cin.get();

	return 0;
}*/

#include <stdio.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <cmath>

#define d 16
#define log2d 1
#define bs 8


void matrix_by_matrix(float matrix1[d][d], float matrix2[d][d], float matrix3[d][d])
{
	//matrix3.resize(matrix1.size());
	for (std::size_t row = 0; row < d; row++)
	{
		//matrix3[row].resize(matrix1[row].size());
		for (std::size_t col = 0; col < d; col++)
		{
			matrix3[row][col] = 0.00;
			for (std::size_t k = 0; k < d; k++)
				matrix3[row][col] += matrix1[row][k] * matrix2[k][col];
		}
	}
}
void print_matrix(float	matrix[d][d], int rows, int cols)
{
	for (std::size_t row = 0; row < rows; row++)
	{
		for (std::size_t col = 0; col < cols; col++)
			std::cout << std::setprecision(5) << std::fixed << matrix[row][col] << " ";

		std::cout << "\n";
	}

	std::cout << "\n";
}
void Q(float v[d][d], float M[d][d]) {
    //float M[d][d] = {0};

    // Set M to be identity matrix
    for (int i=0; i<d; i++) {
        M[i][i] = 1;
    }

    for (int i=0; i<d; i++) { // for each column
        float I[d][d] = {0};
        for (int j=0; j<d; j++) {
            I[j][j] = 1;
        }

        for (int j=0; j<d; j++) {
            for (int k=0; k<d; k++) {
                //std::cout << "j = " << j << "\tk = " << k << "\tComputing " << I[j][k] << " - 2 * " << v[j][i] << " * " << v[k][i] << "\n";
                I[j][k] = I[j][k] - 2 * v[j][i] * v[k][i]; // H
            }
        }
        //std::cout << "H() is returning\n";
        //print_matrix(I, d, d);

        // M = M * H
        float newM[d][d];
        //std::cout << "M = \n";
        //print_matrix(M, d, d);
        matrix_by_matrix(M, I, newM);
        //print_matrix(newM, d, d);

        // Copy the values
        for (int j=0; j<d; j++) {
            for (int k=0; k<d; k++) {
                M[j][k] = newM[j][k];
            }
        }

    }
}

void fasthpp(float v[d][d], float x[d][bs]) {
    // Y = v.T
    float y[d][d];
    float w[d][d];
    for (int i=0; i<d; i++) {
        for (int j=0; j<d; j++) {
            y[i][j] = v[j][i];
            w[i][j] = -2*v[j][i];
        }
    }

    // "Step 1"
    int k=1;
    for (int i=0; i<log2d; i++) {
        int k_2 = k;
        k = k*2;

        float m1[d/2][1];
        for (int row=0; row<d/2; row++) {
            m1[row][0] = 0;
            for (int col=0; col<d; col++) {
                m1[row][0] += y[row*2][col] * w[row*2+1][col];
            }
        }
        //std::cout << "m1 = [" << m1[0][0] << ", " << m1[1][0] << "]\n";

        float m2[d/2][d];
        for (int j=0; j<d/2; j++) {
            for (int row=0; row<d; row++) {
                m2[j][row] = w[j*2][row] * m1[j][0];
            }
        }
        //std::cout << "m2 =\n";
        //print_matrix(m2, d/2, d);

        for (int row=0; row<d/2; row++) {
            for (int col=0; col<d; col++) {
                w[2*row+1][col] += m2[row][col];
            }
        }

        //std::cout << "w =\n";
        //print_matrix(w, d, d);
    }

    // "Step 2"
    //for (int j=d/k; j>-1; j--) { // Variable loop bound
    for (int j=d/k-1; j>-1; j--) {
    	float newX[d][bs] = {0};

    	        //std::cout << "Y = \n";
    	        //print_matrix(y, d, d);

    	       /*std::cout << "X =\n";
    	        for (std::size_t row = 0; row < d; row++)
    	        {
    	            for (std::size_t col = 0; col < bs; col++)
    	                std::cout << std::setprecision(5) << std::fixed << x[row][col] << " ";

    	            std::cout << "\n";
    	        }*/

    	        for (int row=0; row<k; row++) {
    	            for (int col=0; col<bs; col++) {
    	                for (int l=0; l<d; l++) {
    	                    newX[row][col] += y[j*k+row][l] * x[l][col];
    	                }
    	            }
    	        }
    	        /*std::cout << "\nW =\n";
    	        for (std::size_t row = 0; row < d; row++)
    	        {
    	            for (std::size_t col = 0; col < d; col++)
    	                std::cout << std::setprecision(5) << std::fixed << w[row][col] << " ";

    	            std::cout << "\n";
    	        }

    	        std::cout << "\n";*/

    	        float newX1[d][bs] = {0};
    	        for (int row=0; row<d; row++) {
    	            for (int col=0; col<bs; col++) {
    	                for (int l=0; l<k; l++) {
    	                    //std::cout << "w[" << l
    	                    newX1[row][col] += w[j*k+l][row] * newX[l][col];
    	                }
    	            }
    	        }
    	        /*std::cout << "\nnewX1 =\n";
    	        for (std::size_t row = 0; row < d; row++)
    	        {
    	            for (std::size_t col = 0; col < bs; col++)
    	                std::cout << std::setprecision(5) << std::fixed << newX1[row][col] << " ";

    	            std::cout << "\n";
    	        }
    	        std::cout << "\n";*/

    	        for (int row=0; row<d; row++) {
    	            for (int col=0; col<bs; col++) {
    	                x[row][col] += newX1[row][col];
    	            }
    	        }

    	        /*std::cout << "\nX =\n";
    	        for (std::size_t row = 0; row < d; row++)
    	        {
    	            for (std::size_t col = 0; col < bs; col++)
    	                std::cout << std::setprecision(5) << std::fixed << x[row][col] << " ";

    	            std::cout << "\n";
    	        }

    	        std::cout << "\n";*/




    }
}


int main() {

	//Timer variables to measure time
	u32 timerCount_Stop;
	u32 timerCount_Start;

	//AXI device object for AXI Timer
	XTmrCtr timer;
	//Initialize timer object
	int status;
	status = XTmrCtr_Initialize(&timer, TMRCTR_DEVICE_ID);
	if (status != XST_SUCCESS)
	{
		xil_printf("Timer init fail\n\r");
		return XST_FAILURE;
	}

	status = XTmrCtr_SelfTest(&timer, 0);
	if (status != XST_SUCCESS)
	{
		xil_printf("Timer self test fail\n\r");
		return XST_FAILURE;
	}

    float V[d][d] = {0}; // Householder vectors
    float X[d][bs];
    // V[0][0] = 0.8206;
    // V[0][1] = 0.4881;
    // V[1][0] = 0.5715;
    // V[1][1] = 0.8728;
    /*V[0][0] = 0.8395;
    V[0][1] = 0.5717;
    V[0][2] = 0.7215;
    V[0][3] = -0.681;
    V[1][0] = 0.2956;
    V[1][1] = -0.4746;
    V[1][2] = -0.0345;
    V[1][3] = -0.519;
    V[2][0] = -0.3277;
    V[2][1] = 0.6338;
    V[2][2] = -0.3144;
    V[2][3] = -0.454;
    V[3][0] = -0.3171;
    V[3][1] = -0.215;
    V[3][2] = -0.6159;
    V[3][3] = 0.2466;

    X[0][0] = 1.6423;
    X[0][1] = -0.1596;
    X[0][2] = -0.4974;
    X[0][3] = 0.4396;
    X[0][4] = -0.7581;
    X[0][5] = 1.0783;
    X[0][6] = 0.8008;
    X[0][7] = 1.6806;
    X[1][0] = 1.2791;
    X[1][1] = 1.2964;
    X[1][2] = 0.6105;
    X[1][3] = 1.3347;
    X[1][4] = -0.2316;
    X[1][5] = 0.0418;
    X[1][6] = -0.2516;
    X[1][7] = 0.8599;
    X[2][0] = -1.3847;
    X[2][1] = -0.8712;
    X[2][2] = -0.2234;
    X[2][3] = 1.7174;
    X[2][4] = 0.3189;
    X[2][5] = -0.4245;
    X[2][6] = 0.3057;
    X[2][7] = -0.7746;
    X[3][0] = -1.5576;
    X[3][1] = 0.9956;
    X[3][2] = -0.8798;
    X[3][3] = -0.6011;
    X[3][4] = -1.2742;
    X[3][5] = 2.1228;
    X[3][6] = -1.2347;
    X[3][7] = -0.4879;*/
    /*V[0][0] = 0.4854;
    V[0][1] = 0.5008;
    V[0][2] = 0.5791;
    V[0][3] = -0.5974;
    V[0][4] = 0.358;
    V[0][5] = -0.349;
    V[0][6] = -0.0193;
    V[0][7] = -0.5071;
    V[1][0] = -0.1895;
    V[1][1] = 0.5551;
    V[1][2] = -0.2523;
    V[1][3] = -0.3982;
    V[1][4] = -0.3841;
    V[1][5] = -0.1582;
    V[1][6] = -0.3445;
    V[1][7] = 0.2409;
    V[2][0] = 0.4137;
    V[2][1] = -0.0537;
    V[2][2] = -0.3198;
    V[2][3] = 0.1247;
    V[2][4] = -0.4001;
    V[2][5] = 0.3049;
    V[2][6] = 0.3588;
    V[2][7] = 0.5311;
    V[3][0] = 0.3222;
    V[3][1] = 0.4365;
    V[3][2] = 0.3925;
    V[3][3] = 0.3787;
    V[3][4] = -0.1222;
    V[3][5] = 0.0118;
    V[3][6] = -0.1127;
    V[3][7] = 0.2717;
    V[4][0] = -0.3488;
    V[4][1] = -0.2933;
    V[4][2] = -0.1436;
    V[4][3] = 0.4873;
    V[4][4] = 0.1683;
    V[4][5] = -0.12;
    V[4][6] = 0.137;
    V[4][7] = -0.2448;
    V[5][0] = -0.3924;
    V[5][1] = 0.3352;
    V[5][2] = -0.5656;
    V[5][3] = -0.1706;
    V[5][4] = -0.6723;
    V[5][5] = 0.6002;
    V[5][6] = -0.5532;
    V[5][7] = -0.1542;
    V[6][0] = -0.2302;
    V[6][1] = -0.2216;
    V[6][2] = 0.0502;
    V[6][3] = 0.1492;
    V[6][4] = -0.2575;
    V[6][5] = 0.3368;
    V[6][6] = -0.3647;
    V[6][7] = -0.2326;
    V[7][0] = -0.3535;
    V[7][1] = 0.0121;
    V[7][2] = -0.0408;
    V[7][3] = 0.1917;
    V[7][4] = -0.0516;
    V[7][5] = 0.5215;
    V[7][6] = -0.5307;
    V[7][7] = 0.4372;
    X[0][0] = 1.4451;
    X[0][1] = 0.8564;
    X[0][2] = 2.2181;
    X[0][3] = 0.5232;
    X[0][4] = 0.3466;
    X[0][5] = -0.1973;
    X[0][6] = -1.0546;
    X[0][7] = 1.278;
    X[1][0] = -0.1722;
    X[1][1] = 0.5238;
    X[1][2] = 0.0566;
    X[1][3] = 0.4263;
    X[1][4] = 0.575;
    X[1][5] = -0.6417;
    X[1][6] = -2.2064;
    X[1][7] = -0.7508;
    X[2][0] = 0.0109;
    X[2][1] = -0.3387;
    X[2][2] = -1.3407;
    X[2][3] = -0.5854;
    X[2][4] = 0.5362;
    X[2][5] = 0.5246;
    X[2][6] = 1.1412;
    X[2][7] = 0.0516;
    X[3][0] = 0.744;
    X[3][1] = -0.4816;
    X[3][2] = -1.0495;
    X[3][3] = 0.6039;
    X[3][4] = -1.7223;
    X[3][5] = -0.8278;
    X[3][6] = 1.3347;
    X[3][7] = 0.4835;
    X[4][0] = -2.5095;
    X[4][1] = 0.488;
    X[4][2] = 0.7846;
    X[4][3] = 0.0286;
    X[4][4] = 0.6408;
    X[4][5] = 0.5832;
    X[4][6] = 1.0669;
    X[4][7] = -0.4502;
    X[5][0] = -0.1853;
    X[5][1] = 0.7528;
    X[5][2] = 0.4048;
    X[5][3] = 0.1785;
    X[5][4] = 0.2649;
    X[5][5] = 1.2732;
    X[5][6] = -0.0013;
    X[5][7] = -0.3036;
    X[6][0] = -1.457;
    X[6][1] = -0.1023;
    X[6][2] = -0.5992;
    X[6][3] = 0.4771;
    X[6][4] = 0.7262;
    X[6][5] = 0.0912;
    X[6][6] = -0.3891;
    X[6][7] = 0.5279;
    X[7][0] = -0.0127;
    X[7][1] = 0.2408;
    X[7][2] = 0.1325;
    X[7][3] = 0.7642;
    X[7][4] = 1.095;
    X[7][5] = 0.3399;
    X[7][6] = 0.72;
    X[7][7] = 0.4114;*/

    V[0][0] = 0.3373;
    V[0][1] = 0.4428;
    V[0][2] = 0.2279;
    V[0][3] = -0.5619;
    V[0][4] = 0.2712;
    V[0][5] = -0.3712;
    V[0][6] = -0.0153;
    V[0][7] = -0.3291;
    V[0][8] = -0.2048;
    V[0][9] = 0.3833;
    V[0][10] = -0.0882;
    V[0][11] = -0.3784;
    V[0][12] = -0.1769;
    V[0][13] = -0.1449;
    V[0][14] = -0.1875;
    V[0][15] = 0.2076;
    V[1][0] = 0.2875;
    V[1][1] = -0.0475;
    V[1][2] = -0.1258;
    V[1][3] = 0.1173;
    V[1][4] = -0.3031;
    V[1][5] = 0.3242;
    V[1][6] = 0.2844;
    V[1][7] = 0.3447;
    V[1][8] = 0.3482;
    V[1][9] = 0.3014;
    V[1][10] = 0.1372;
    V[1][11] = 0.3599;
    V[1][12] = -0.0563;
    V[1][13] = 0.0108;
    V[1][14] = -0.0613;
    V[1][15] = 0.2341;
    V[2][0] = -0.2424;
    V[2][1] = -0.2594;
    V[2][2] = -0.0565;
    V[2][3] = 0.4583;
    V[2][4] = 0.1275;
    V[2][5] = -0.1276;
    V[2][6] = 0.1086;
    V[2][7] = -0.1588;
    V[2][8] = -0.424;
    V[2][9] = 0.2315;
    V[2][10] = -0.1977;
    V[2][11] = -0.1621;
    V[2][12] = -0.3096;
    V[2][13] = 0.55;
    V[2][14] = -0.3011;
    V[2][15] = -0.1328;
    V[3][0] = -0.16;
    V[3][1] = -0.1959;
    V[3][2] = 0.0197;
    V[3][3] = 0.1403;
    V[3][4] = -0.1951;
    V[3][5] = 0.3582;
    V[3][6] = -0.2891;
    V[3][7] = -0.1509;
    V[3][8] = -0.382;
    V[3][9] = 0.0084;
    V[3][10] = -0.0143;
    V[3][11] = 0.1822;
    V[3][12] = -0.0238;
    V[3][13] = 0.4779;
    V[3][14] = -0.2888;
    V[3][15] = 0.3767;
    V[4][0] = 0.253;
    V[4][1] = 0.255;
    V[4][2] = 0.5611;
    V[4][3] = 0.1396;
    V[4][4] = 0.1386;
    V[4][5] = -0.0593;
    V[4][6] = -0.3746;
    V[4][7] = 0.2621;
    V[4][8] = -0.0469;
    V[4][9] = 0.1218;
    V[4][10] = 0.0127;
    V[4][11] = 0.1149;
    V[4][12] = 0.1397;
    V[4][13] = -0.1663;
    V[4][14] = -0.538;
    V[4][15] = -0.2044;
    V[5][0] = 0.0019;
    V[5][1] = -0.1008;
    V[5][2] = -0.3392;
    V[5][3] = -0.1562;
    V[5][4] = 0.2144;
    V[5][5] = 0.1577;
    V[5][6] = 0.4054;
    V[5][7] = 0.0106;
    V[5][8] = 0.2025;
    V[5][9] = -0.112;
    V[5][10] = -0.2358;
    V[5][11] = 0.1628;
    V[5][12] = -0.4185;
    V[5][13] = -0.2145;
    V[5][14] = 0.3254;
    V[5][15] = 0.1316;
    V[6][0] = -0.4393;
    V[6][1] = 0.1453;
    V[6][2] = 0.1985;
    V[6][3] = 0.0076;
    V[6][4] = 0.2562;
    V[6][5] = 0.1754;
    V[6][6] = 0.379;
    V[6][7] = -0.0923;
    V[6][8] = -0.0504;
    V[6][9] = 0.175;
    V[6][10] = 0.0909;
    V[6][11] = 0.0481;
    V[6][12] = 0.0644;
    V[6][13] = 0.3298;
    V[6][14] = -0.0003;
    V[6][15] = -0.0827;
    V[7][0] = -0.255;
    V[7][1] = -0.0305;
    V[7][2] = -0.1516;
    V[7][3] = 0.1273;
    V[7][4] = 0.2903;
    V[7][5] = 0.0274;
    V[7][6] = -0.1382;
    V[7][7] = 0.1083;
    V[7][8] = -0.0035;
    V[7][9] = 0.056;
    V[7][10] = 0.0298;
    V[7][11] = 0.206;
    V[7][12] = 0.2661;
    V[7][13] = 0.0881;
    V[7][14] = 0.1756;
    V[7][15] = 0.112;
    V[8][0] = 0.338;
    V[8][1] = 0.3012;
    V[8][2] = -0.3634;
    V[8][3] = -0.3015;
    V[8][4] = -0.0544;
    V[8][5] = 0.4917;
    V[8][6] = 0.2326;
    V[8][7] = 0.1181;
    V[8][8] = 0.3108;
    V[8][9] = 0.0043;
    V[8][10] = -0.4057;
    V[8][11] = 0.2495;
    V[8][12] = -0.0912;
    V[8][13] = 0.2677;
    V[8][14] = -0.1674;
    V[8][15] = 0.1734;
    V[9][0] = -0.1703;
    V[9][1] = 0.2853;
    V[9][2] = 0.4096;
    V[9][3] = 0.3871;
    V[9][4] = 0.1077;
    V[9][5] = -0.0632;
    V[9][6] = -0.2603;
    V[9][7] = 0.0214;
    V[9][8] = 0.0949;
    V[9][9] = 0.2249;
    V[9][10] = -0.1046;
    V[9][11] = 0.4327;
    V[9][12] = -0.6027;
    V[9][13] = -0.1082;
    V[9][14] = -0.2915;
    V[9][15] = 0.2212;
    V[10][0] = -0.3327;
    V[10][1] = 0.0681;
    V[10][2] = 0.0063;
    V[10][3] = -0.0923;
    V[10][4] = 0.1147;
    V[10][5] = -0.2197;
    V[10][6] = 0.0621;
    V[10][7] = -0.2243;
    V[10][8] = -0.4362;
    V[10][9] = 0.3145;
    V[10][10] = 0.2896;
    V[10][11] = 0.0141;
    V[10][12] = -0.3759;
    V[10][13] = 0.196;
    V[10][14] = 0.1891;
    V[10][15] = 0.5517;
    V[11][0] = 0.0063;
    V[11][1] = 0.0359;
    V[11][2] = -0.2038;
    V[11][3] = -0.0554;
    V[11][4] = -0.3726;
    V[11][5] = -0.4783;
    V[11][6] = -0.4035;
    V[11][7] = -0.1072;
    V[11][8] = -0.1412;
    V[11][9] = -0.349;
    V[11][10] = -0.4329;
    V[11][11] = 0.0345;
    V[11][12] = 0.2486;
    V[11][13] = -0.144;
    V[11][14] = 0.1717;
    V[11][15] = 0.1933;
    V[12][0] = 0.3106;
    V[12][1] = -0.2744;
    V[12][2] = 0.2435;
    V[12][3] = -0.0899;
    V[12][4] = -0.4699;
    V[12][5] = 0.1076;
    V[12][6] = 0.1701;
    V[12][7] = 0.2776;
    V[12][8] = 0.1432;
    V[12][9] = 0.491;
    V[12][10] = -0.117;
    V[12][11] = -0.2513;
    V[12][12] = 0.045;
    V[12][13] = 0.2769;
    V[12][14] = 0.3186;
    V[12][15] = 0.1252;
    V[13][0] = -0.1426;
    V[13][1] = -0.304;
    V[13][2] = -0.1252;
    V[13][3] = -0.158;
    V[13][4] = 0.0617;
    V[13][5] = 0.1325;
    V[13][6] = -0.0527;
    V[13][7] = -0.4755;
    V[13][8] = -0.1084;
    V[13][9] = 0.2512;
    V[13][10] = -0.4001;
    V[13][11] = 0.4066;
    V[13][12] = 0.0752;
    V[13][13] = -0.1296;
    V[13][14] = 0.2524;
    V[13][15] = 0.46;
    V[14][0] = -0.0008;
    V[14][1] = 0.4962;
    V[14][2] = 0.0389;
    V[14][3] = -0.2829;
    V[14][4] = -0.229;
    V[14][5] = 0.0251;
    V[14][6] = 0.142;
    V[14][7] = 0.4079;
    V[14][8] = -0.0196;
    V[14][9] = -0.2106;
    V[14][10] = -0.4603;
    V[14][11] = -0.2915;
    V[14][12] = 0.0043;
    V[14][13] = 0.0203;
    V[14][14] = 0.0471;
    V[14][15] = 0.1115;
    V[15][0] = -0.1626;
    V[15][1] = 0.0822;
    V[15][2] = -0.1363;
    V[15][3] = 0.1234;
    V[15][4] = -0.3486;
    V[15][5] = -0.0082;
    V[15][6] = -0.1255;
    V[15][7] = 0.3002;
    V[15][8] = 0.3418;
    V[15][9] = -0.1662;
    V[15][10] = 0.1919;
    V[15][11] = 0.1383;
    V[15][12] = 0.1312;
    V[15][13] = 0.1465;
    V[15][14] = 0.1233;
    V[15][15] = 0.0606;
    X[0][0] = -0.6855;
    X[0][1] = 0.5636;
    X[0][2] = -1.5072;
    X[0][3] = -1.6107;
    X[0][4] = -1.479;
    X[0][5] = 0.4323;
    X[0][6] = -0.125;
    X[0][7] = 0.7821;
    X[1][0] = -1.5988;
    X[1][1] = -0.1091;
    X[1][2] = 0.7152;
    X[1][3] = 0.0391;
    X[1][4] = 1.3059;
    X[1][5] = 0.2466;
    X[1][6] = -1.9776;
    X[1][7] = 0.0179;
    X[2][0] = -1.3793;
    X[2][1] = 0.6258;
    X[2][2] = -2.585;
    X[2][3] = -0.024;
    X[2][4] = -0.1222;
    X[2][5] = -0.747;
    X[2][6] = 1.7093;
    X[2][7] = 0.0579;
    X[3][0] = 1.193;
    X[3][1] = 1.9373;
    X[3][2] = 0.7287;
    X[3][3] = 0.9809;
    X[3][4] = 0.4146;
    X[3][5] = 1.1566;
    X[3][6] = 0.2691;
    X[3][7] = -0.0366;
    X[4][0] = 0.9733;
    X[4][1] = -1.0151;
    X[4][2] = -0.5419;
    X[4][3] = -0.441;
    X[4][4] = -0.3136;
    X[4][5] = -0.1293;
    X[4][6] = -0.715;
    X[4][7] = -0.0476;
    X[5][0] = 2.0207;
    X[5][1] = 0.2539;
    X[5][2] = 0.9364;
    X[5][3] = 0.7122;
    X[5][4] = -0.0318;
    X[5][5] = 0.1016;
    X[5][6] = 1.3433;
    X[5][7] = 0.7133;
    X[6][0] = 0.4038;
    X[6][1] = -0.714;
    X[6][2] = 0.8337;
    X[6][3] = -0.9585;
    X[6][4] = 0.4536;
    X[6][5] = 1.2461;
    X[6][6] = -2.3065;
    X[6][7] = -1.2869;
    X[7][0] = 0.1799;
    X[7][1] = -2.1268;
    X[7][2] = -0.1341;
    X[7][3] = -1.0408;
    X[7][4] = -0.7647;
    X[7][5] = -0.0553;
    X[7][6] = 1.2049;
    X[7][7] = -0.9825;
    X[8][0] = 0.4334;
    X[8][1] = -0.7172;
    X[8][2] = 1.0554;
    X[8][3] = -1.4534;
    X[8][4] = 0.4652;
    X[8][5] = 0.3714;
    X[8][6] = -0.0047;
    X[8][7] = 0.0795;
    X[9][0] = 0.3782;
    X[9][1] = 0.7051;
    X[9][2] = -1.7237;
    X[9][3] = -0.8435;
    X[9][4] = 0.4351;
    X[9][5] = 0.2659;
    X[9][6] = -0.5871;
    X[9][7] = 0.0827;
    X[10][0] = 0.8854;
    X[10][1] = 0.1824;
    X[10][2] = 0.7864;
    X[10][3] = -0.0579;
    X[10][4] = 0.5667;
    X[10][5] = -0.7098;
    X[10][6] = -0.4875;
    X[10][7] = 0.0501;
    X[11][0] = 0.6084;
    X[11][1] = 1.6309;
    X[11][2] = -0.0847;
    X[11][3] = 1.0844;
    X[11][4] = 0.9478;
    X[11][5] = -0.6766;
    X[11][6] = -0.573;
    X[11][7] = -0.3303;
    X[12][0] = -0.7939;
    X[12][1] = 0.3752;
    X[12][2] = 0.0879;
    X[12][3] = -1.2415;
    X[12][4] = -0.3203;
    X[12][5] = -0.8444;
    X[12][6] = -0.5513;
    X[12][7] = 1.989;
    X[13][0] = 1.9003;
    X[13][1] = 1.6951;
    X[13][2] = 0.0281;
    X[13][3] = -0.1754;
    X[13][4] = -1.7735;
    X[13][5] = -0.7046;
    X[13][6] = -0.3947;
    X[13][7] = 1.8868;
    X[14][0] = -0.2184;
    X[14][1] = 0.1663;
    X[14][2] = 2.1442;
    X[14][3] = 1.7046;
    X[14][4] = 0.3459;
    X[14][5] = 0.6425;
    X[14][6] = -0.204;
    X[14][7] = 0.6854;
    X[15][0] = -0.1397;
    X[15][1] = -1.1808;
    X[15][2] = -1.2829;
    X[15][3] = 0.4485;
    X[15][4] = -0.5907;
    X[15][5] = 0.8541;
    X[15][6] = -0.4901;
    X[15][7] = -0.3595;

    print_matrix(V, d, d);

    float M[d][d] = {0};
    print_matrix(M, d, d);
    float out[d][bs] = {0};

    //Reset the timer
	XTmrCtr_Reset(&timer, 0);
	XTmrCtr_ClearStats(&timer);
	XTmrCtr_Start(&timer, 0);
	// Start timer
	timerCount_Start = XTmrCtr_GetValue(&timer, 0);

    Q(V, M); // Naive Householder Multiplication
    // Perform the multiplication with the input
    for (int row=0; row<d; row++) {
        for (int col=0; col<bs; col++) {
            for (int k=0; k<d; k++) {
                out[row][col] += M[row][k] * X[k][col];
            }
        }
    }

    timerCount_Stop = XTmrCtr_GetValue(&timer, 0);
	//Display time
	xil_printf("Time Accelerator: %d\n\r", (timerCount_Stop - timerCount_Start));

    std::cout << "\nX =\n";
    for (std::size_t row = 0; row < d; row++)
    {
        for (std::size_t col = 0; col < bs; col++)
            std::cout << std::setprecision(5) << std::fixed << out[row][col] << " ";

        std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "Beginning fasth\n";

    //Reset the timer
	XTmrCtr_Reset(&timer, 0);
	XTmrCtr_ClearStats(&timer);
	XTmrCtr_Start(&timer, 0);
	// Start timer
	timerCount_Start = XTmrCtr_GetValue(&timer, 0);

    fasthpp(V, X);
    timerCount_Stop = XTmrCtr_GetValue(&timer, 0);
	xil_printf("Time Accelerator: %d\n\r", (timerCount_Stop - timerCount_Start));


    std::cout << "\nX =\n";
	for (std::size_t row = 0; row < d; row++)
	{
		for (std::size_t col = 0; col < bs; col++)
			std::cout << std::setprecision(5) << std::fixed << X[row][col] << " ";

		std::cout << "\n";
	}
	std::cout << "\n";
    return 0;
}


