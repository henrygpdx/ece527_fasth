#include <stdio.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <string.h>

#define d 129
#define log2d 1
#define bs 128
using namespace std;

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

void Q(float v[d][d], float M[d][d]) {
    //float M[d][d] = {0};

    // Set M to be identity matrix
    for (int i=0; i<d; i++) {
        M[i][i] = 1;
    }

    for (int i=0; i<d; i++) { // for each column
        float I[d][d] = {0};
        for (int j=0; j<d; j++) {
            // #pragma HLS pipeline
            I[j][j] = 1;
        }

        // M = M * H
        float newM[d][d];
        matrix_by_matrix(M, I, newM);

        // Copy the values
        for (int j=0; j<d; j++) {
            // #pragma HLS pipeline
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
            #pragma HLS pipeline
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
            #pragma HLS pipeline
            m1[row][0] = 0;
            for (int col=0; col<d; col++) {
                m1[row][0] += y[row*2][col] * w[row*2+1][col];
            }
        }

        float m2[d/2][d];
        for (int j=0; j<d/2; j++) {
            #pragma HLS pipeline
            for (int row=0; row<d; row++) {
                m2[j][row] = w[j*2][row] * m1[j][0];
            }
        }

        for (int row=0; row<d/2; row++) {
            #pragma HLS pipeline
            for (int col=0; col<d; col++) {
                w[2*row+1][col] += m2[row][col];
            }
        }
    }

    // "Step 2"
    //for (int j=d/k; j>-1; j--) { // Variable loop bound
    for (int j=d/k-1; j>-1; j--) {
        float newX[d][bs] = {0};

        for (int row=0; row<k; row++) {
            for (int col=0; col<bs; col++) {
                #pragma HLS pipeline
                for (int l=0; l<d; l++) {
                    newX[row][col] += y[j*k+row][l] * x[l][col];
                }
            }
        }

        float newX1[d][bs] = {0};
        for (int row=0; row<d; row++) {
            #pragma HLS pipeline
            for (int col=0; col<bs; col++) {
                for (int l=0; l<k; l++) {
                    newX1[row][col] += w[j*k+l][row] * newX[l][col];
                }
            }
        }

        for (int row=0; row<d; row++) {
            #pragma HLS pipeline
            for (int col=0; col<bs; col++) {
                x[row][col] += newX1[row][col];
            }
        }
    }
}

void read_input(float in[d][bs], float in_buf[d][bs]){
    memcpy(in_buf, in, d*bs*sizeof(float));
}
void write_output(float in[d][bs], float in_buf[d][bs]){
    memcpy(in_buf, in, d*bs*sizeof(float));
}
int main_func(float X_t[d][bs]) {
    #pragma HLS INTERFACE m_axi depth=24*24 port=X_t offset=slave bundle=input
    #pragma HLS INTERFACE s_axilite register port=return bundle=CTL
    float X[d][bs];
    read_input(X_t, X);
    float V[d][d] = {0}; // Householder vectors


    V[0][0] = 0.4854;
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
    X[7][7] = 0.4114;

    float M[d][d] = {0};
    float out[d][bs] = {0};

    Q(V, M);
    for (int row=0; row<d; row++) {
        // #pragma HLS unroll factor=5

        for (int col=0; col<bs; col++) {
            // #pragma HLS pipeline
            for (int k=0; k<d; k++) {
                out[row][col] += M[row][k] * X[k][col];
            }
        }
    }

    // fasthpp(V, X);
    write_output(out, X_t);
    return 0;
}
