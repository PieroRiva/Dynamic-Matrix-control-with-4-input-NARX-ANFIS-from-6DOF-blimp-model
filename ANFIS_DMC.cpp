/* ANFISblimpCtrl.cpp */

//-------------------------------| 
//          LIBRARIES            |
//-------------------------------| 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//-------------------------------| 
//       GLOBAL VARIABLES        |
//-------------------------------| 
#define N 9
const double pi = 3.14159265359;
const double tol = 1e-6;
const int mInputs = 3;
const int mStates = 6;
const int nData = 6000;

const int nInputs = 4;
const int nFuzzy = 5;
const int nRules = 625;

const int Np = 4;
const int Nc = 3;

double tt[nData];
double UT[mInputs][nData];
double U[mInputs][nData];
double X[mStates][nData];
double Y[mStates][nData];
double Vearth[3][nData];
double Xearth[3][nData];
double lambda[3][3];
double dX[mStates][nData];
double Vtot[nData];
double Wtot[nData];
double G[mStates][nData];
double A[mStates][nData];
double Fd[mStates][nData];
double P[mStates][nData];

double OUTPUT[nData];
double INPUT1[nData];
double INPUT2[nData];
double INPUT3[nData];
double INPUT4[nData];

double O5[nData];
double En[nData];
double muIn[nFuzzy][nInputs];
double w[nRules];
double wn[nRules];
double fi[nRules];
double fi_wn[nRules];
double sumW;

double aIne[nFuzzy][nInputs];
double cIne[nFuzzy][nInputs];
double aOute[nRules];
double cOute[nRules];

double XX[nInputs];
double ukk[mInputs][nData];
double dukk[mInputs * Nc];
double xkk[mStates][nData];
double ykk[mStates][nData];
double yee[nData];

double ref[nData];
double wref[nData];

double Q[Np][Np];
double R[Nc * mInputs][Nc * mInputs];

double gi[Np];
double Gij[Np][mInputs * Nc];
double Wref[Np];
double Yfree[Np];
double Ystep[Np][mInputs];

double Gijt[mInputs * Nc][Np];
double Gijt_Q[mInputs * Nc][Np];
double Gijt_Q_Gij[mInputs * Nc][mInputs * Nc];
double Gijt_Q_Gij_R[mInputs * Nc][mInputs * Nc];
double Gijt_Q_Gij_R_1[mInputs * Nc][mInputs * Nc];
double Gijt_Q_Gij_R_1_Gijt[mInputs * Nc][Np];
double Gijt_Q_Gij_R_1_Gijt_Q[mInputs * Nc][Np];
double Wref_Yfree[Np];



//-------------------------------| 
//    FUNCTION DECLARATION       |
//-------------------------------| 

double h(double t);
double r(double t);
void fill(double* arr, int n, int m, double val);
double sign(double x);
double map(double input, double minIn, double maxIn, double minOut, double maxOut);

double reference(double t);

double gauss(double x, double a, double c);
double sigmoid(double x, double a, double c);
double invSigmoid(double x, double a, double c);
double dGauss_da(double x, double a, double c);
double dGauss_dc(double x, double a, double c);
double dinvSigmoid_da(double x, double a, double c);
double dinvSigmoid_dc(double x, double a, double c);
void matrixMultiplication(int a, double* A, int c, double* B, int b, double* C);
void matrixTranspose(int a, double* A, int b, double* B);
void matrixSum(int n, int m, double* A, double* B, double* C);
void matrixSubstract(int n, int m, double* A, double* B, double* C);
void set_ANFIS();

void getCfactor(double M[N][N], double t[N][N], int p, int q, int n);
double DET(double M[N][N], int n);
void ADJ(double M[N][N], double adj[N][N]);
bool INV(double M[N][N], double inv[N][N]);


void mINVERSE(double M[81], double inverse[81]);







//-------------------------------| 
//            MAIN()             |
//-------------------------------| 
int main()
{
    printf("CPP_ANFISblimpCtrl.exe");
    printf("\n\nPiero A. Riva Riquelme, pieroriva@udec.cl");
    printf("\nUndergraduate student at Universidad de Concepcion, Chile");
    printf("\nJanuary of 2021");
    // Step 1: Loading ANFIS
    printf("\n\nStep 1: Loading ANFIS parameters...");
    set_ANFIS();

    // Step 2: Loading blimp parameters
    printf("\n\nStep 2: Loading blimp parameters...");
    // a) Physical constants
    printf("\n\ta) Physical constants");
    const double rhoAir = 1.205;            // Density of air at NTP (20°C, 1atm)
    const double rhoHe = 0.1664;           // Density of Helium at NTP (20°C, 1atm)
    const double g_acc = 9.80665;           // Acceleration of gravity
    const double deg2rad = pi / 180;          // Degrees to radians conversion 
    const double rad2deg = pow(deg2rad, -1.0); // Radians to degrees conversion

    // b) Vehicle geometry and parameters
    printf("\n\tb) Vehicle geometry and parameters");
    const double blimp_a = 0.9;                                     // Blimp's x length 
    const double blimp_b = 0.45;                                    // Blimp's y length 
    const double blimp_c = 0.45;                                    // Blimp's z length 
    const double blimp_volume = 4.0 * pi * blimp_a * blimp_b * blimp_c / 3.0;    // Blimp's volume 
    const double blimp_area = pow(blimp_volume, 0.6666666666);      // Blimp's area 
    const double blimp_mHe = blimp_volume * rhoHe;                    // Blimp's mass of helium 
    const double blimp_mAir = blimp_volume * rhoAir;                  // Blimp's mass of air
    const double blimp_mass = blimp_mAir - blimp_mHe;               // Blimp's mass (chosen for 0 buoyancy)
    const double blimp_mTotal = blimp_mass + blimp_mHe;               // Blimp's total mass
    const double blimp_dx = blimp_a / 8.0;                              // Blimp's x axis distace from CV to propellers
    const double blimp_dy = blimp_b / 2.0;                              // Blimp's y axis distace from CV to propellers
    const double blimp_dz = blimp_c;                                // Blimp's z axis distace from CV to propellers
    const double blimp_ax = 0.0;                                      // Blimp's x axis distance from center of gravity CG and center of volume CV
    //const double blimp_ay = 0.0;                                      // Blimp's y axis distance from center of gravity CG and center of volume CV
    const double blimp_az = -(0.2 * blimp_mass) * blimp_b / blimp_mTotal;   // Blimp's z axis distance from center of gravity CG and center of volume CV

    // c) Masses and inertias
    printf("\n\tc) Masses and inertias");
    const double blimp_Ix = blimp_mTotal * (blimp_b * blimp_b + blimp_c * blimp_c) / 5.0;
    const double blimp_Iy = blimp_mTotal * (blimp_c * blimp_c + blimp_a * blimp_a) / 5.0;
    const double blimp_Iz = blimp_mTotal * (blimp_a * blimp_a + blimp_b * blimp_b) / 5.0;

    // c.1) Tuckerman fo a prolate ellipsoid
    const double tuckerman_e = sqrt(1.0 - blimp_c * blimp_c / (blimp_a * blimp_a));
    const double tuckerman_alpha = (1.0 - tuckerman_e * tuckerman_e) * (log((1.0 + tuckerman_e) / (1.0 - tuckerman_e)) - 2.0 * tuckerman_e) / (tuckerman_e * tuckerman_e * tuckerman_e);
    const double tuckerman_beta = (1.0 - tuckerman_e * tuckerman_e) * ((tuckerman_e / (1.0 - tuckerman_e * tuckerman_e)) - 0.5 * log((1.0 + tuckerman_e) / (1.0 - tuckerman_e))) / (tuckerman_e * tuckerman_e * tuckerman_e);
    const double tuckerman_gamma = tuckerman_beta;

    const double tuckerman_K1 = blimp_volume * (tuckerman_alpha / (2.0 - tuckerman_alpha));
    const double tuckerman_K2 = blimp_volume * (tuckerman_beta / (2.0 - tuckerman_beta));
    const double tuckerman_K3 = blimp_volume * (tuckerman_gamma / (2.0 - tuckerman_gamma));
    const double tuckerman_K1_ = blimp_volume * blimp_Ix * (pow((blimp_b * blimp_b - blimp_c * blimp_c) / (blimp_b * blimp_b + blimp_c * blimp_c), 2.0) * ((tuckerman_gamma - tuckerman_beta) / (2.0 * ((blimp_b * blimp_b - blimp_c * blimp_c) / (blimp_b * blimp_b + blimp_c * blimp_c)) - (tuckerman_gamma - tuckerman_beta))));
    const double tuckerman_K2_ = blimp_volume * blimp_Iy * (pow((blimp_c * blimp_c - blimp_a * blimp_a) / (blimp_c * blimp_c + blimp_a * blimp_a), 2.0) * ((tuckerman_alpha - tuckerman_gamma) / (2.0 * ((blimp_c * blimp_c - blimp_a * blimp_a) / (blimp_c * blimp_c + blimp_a * blimp_a)) - (tuckerman_alpha - tuckerman_gamma))));
    const double tuckerman_K3_ = blimp_volume * blimp_Iz * (pow((blimp_a * blimp_a - blimp_b * blimp_b) / (blimp_a * blimp_a + blimp_b * blimp_b), 2.0) * ((tuckerman_beta - tuckerman_alpha) / (2.0 * ((blimp_a * blimp_a - blimp_b * blimp_b) / (blimp_a * blimp_a + blimp_b * blimp_b)) - (tuckerman_beta - tuckerman_alpha))));

    // c.2) Virtual masses and inertias
        // Tuckerman
    const double blimp_Xu = -tuckerman_K1 * rhoAir;
    const double blimp_Yv = -tuckerman_K2 * rhoAir;
    const double blimp_Zw = -tuckerman_K3 * rhoAir;
    const double blimp_Lp = 0.0;
    const double blimp_Mq = -tuckerman_K2_ * rhoAir;
    const double blimp_Nr = -tuckerman_K3_ * rhoAir;

    // Gomes
    const double blimp_Mu = 0.0;
    const double blimp_Lv = 0.0;
    const double blimp_Nv = 0.0;
    const double blimp_Mw = 0.0;
    const double blimp_Yp = 0.0;
    const double blimp_Xq = 0.0;
    const double blimp_Zq = 0.0;
    const double blimp_Yr = 0.0;

    // Groups
    const double blimp_mx = blimp_mTotal - blimp_Xu;
    const double blimp_my = blimp_mTotal - blimp_Yv;
    const double blimp_mz = blimp_mTotal - blimp_Zw;
    const double blimp_Jx = blimp_Ix - blimp_Lp;
    const double blimp_Jy = blimp_Iy - blimp_Mq;
    const double blimp_Jz = blimp_Iz - blimp_Nr;
    const double blimp_Jxz = 0.0;


    // d) M matrix
    printf("\n\td) M matrix");
    const double M[6][6] = {
        {blimp_mx                           , 0.0                                   , 0.0                                   , 0.0                                   , blimp_mTotal * blimp_az - blimp_Xq    , 0.0                               },
        {0.0                                , blimp_my                              , 0.0                                   , -blimp_mTotal * blimp_az - blimp_Yp   , 0.0                                   , blimp_mTotal * blimp_ax - blimp_Yr},
        {0.0                                , 0.0                                   , blimp_mz                              , 0.0                                   , -blimp_mTotal * blimp_ax - blimp_Zq   , 0.0                               },
        {0.0                                , -blimp_mTotal * blimp_az - blimp_Lv   , 0.0                                   , blimp_Ix - blimp_Lp                   , 0.0                                   , -blimp_Jxz                        },
        {blimp_mTotal * blimp_az - blimp_Mu , 0.0                                   , -blimp_mTotal * blimp_ax - blimp_Mw   , 0.0                                   , blimp_Iy - blimp_Mq                   , 0.0                               },
        {0.0                                , blimp_mTotal * blimp_ax - blimp_Nv    , 0.0                                   , -blimp_Jxz                            , 0.0                                   , blimp_Iz - blimp_Nr               }
    };

    const double invM[6][6] = {
        {0.916844157881075  , 0.0                   , 0.0               , 0.0                   , 0.287823601986896 , 0.0               },
        {0.0                , 0.666945041828523     , 0.0               , -0.638717492340345    , 0.0               , 0.0               },
        {0.0                , 0.0                   , 0.637872073637673 , 0.0                   , 0.0               , 0.0               },
        {0.0                , -0.638717492340345    , 0.0               , 14.032280169794683    , 0.0               , 0.0               },
        {0.287823601986896  , 0.0                   , 0.0               , 0.0                   , 4.489659394858919 , 0.0               },
        {0.0                , 0.0                   , 0.0               , 0.0                   , 0.0               , 4.399303334727424 }
    };

    // Step 3: Control configuration
    printf("\n\nStep 3: Control configuration...");

    // a) Time definition
    printf("\n\ta) Time definition...");
    double ti = 0.1;
    double step = 0.1;
    double tf = 600.0;
    //double tt[nData];
    for (int i = 0; i < nData; i++) {
        tt[i] = ti + step * i;
    }

    // b) Configuration
    printf("\n\tb) Configuration...");
    // b.1) DMC
    printf("\n\t\tb.1) Configuring DMC options...");
    const double alpha = 0.9000;
    const double Q_par = 0001.0000;
    const double R1 = 0005.0000; // 1st Control action weight
    const double R2 = 1000.0000; // 2nd Control action weight
    const double R3 = 0005.0000; // 3rd Control action weight
    const double Pu1 = 0000.1000; // 1st Control action proportional gain
    const double Pu2 = 0000.0000; // 2nd Control action proportional gain
    const double Pu3 = 0000.1000; // 3rd Control action proportional gain

    // Error and control action weight
    for (int i = 0; i < Np; i++) {
        for (int j = 0; j < Np; j++) {
            if (i == j) {
                Q[i][j] = Q_par;
            }
            else {
                Q[i][j] = 0.0;
            }
        }
    }
    for (int i = 0; i < Nc * mInputs; i++) {
        for (int j = 0; j < Nc * mInputs; j++) {
            if (i == j) {
                if (i < Nc) {
                    R[i][j] = R1;
                }
                else
                    if (Nc <= i && i < 2 * Nc) {
                        R[i][j] = R2;
                    }
                    else
                        if (2 * Nc <= i && i < 3 * Nc) {
                            R[i][j] = R3;
                        }
            }
            else {
                R[i][j] = 0.0;
            }
        }
    }


    // b.2) ANFIS
    printf("\n\t\tb.2) Configuring ANFIS options...");

    const double minOut = 10.0;
    const double maxOut = 90.0;

    const double minUT1 = -0.1;
    const double maxUT1 = 0.1;
    const double minUT2 = 0.45;
    const double maxUT2 = 0.55;
    const double minUT3 = -pi / 2;
    const double maxUT3 = pi / 2;

    const double minX1 = -0.4;
    const double maxX1 = 0.4;
    const double minX2 = -1.0;
    const double maxX2 = 1.0;
    const double minX3 = -1.0;
    const double maxX3 = 1.0;
    const double minX4 = -pi;
    const double maxX4 = pi;
    const double minX5 = -pi;
    const double maxX5 = pi;
    const double minX6 = -pi;
    const double maxX6 = pi;

    const double minY1 = -100;
    const double maxY1 = 100;
    const double minY2 = -100;
    const double maxY2 = 100;
    const double minY3 = -100;
    const double maxY3 = 100;
    const double minY4 = -pi;
    const double maxY4 = pi;
    const double minY5 = -pi;
    const double maxY5 = pi;
    const double minY6 = -pi;
    const double maxY6 = pi;

    // c) Workspace
    printf("\n\tc) Creating workspace...");
    // c.1) DMC
    printf("\n\t\tc.1) DMC workspace...");

    fill((double*)gi, Np, 1, 0.0);
    fill((double*)Gij, Np, mInputs * Nc, 0.0);
    fill((double*)Wref, Np, 1, 0.0);
    fill((double*)Yfree, Np, 1, 0.0);
    fill((double*)Ystep, Np, mInputs, 0.0);

    // c.2) Process
    printf("\n\t\tc.2) Process workspace...");
    fill((double*)UT, mInputs, nData, 0.0);
    fill((double*)UT + nData, 1, nData, 0.5);
    fill((double*)U, mInputs, nData, 0.0);
    fill((double*)X, mStates, nData, 0.0);
    fill((double*)Y, mStates, nData, 0.0);
    fill((double*)Vearth, 3, nData, 0.0);
    fill((double*)Xearth, 3, nData, 0.0);
    fill((double*)lambda, 3, 3, 0.0);
    fill((double*)dX, mStates, nData, 0.0);
    fill((double*)Vtot, nData, 1, 0.0);
    fill((double*)Wtot, nData, 1, 0.0);
    fill((double*)G, mStates, nData, 0.0);
    fill((double*)A, mStates, nData, 0.0);
    fill((double*)Fd, mStates, nData, 0.0);
    fill((double*)P, mStates, nData, 0.0);
    double ddu[3];

    double f1;
    double f2;
    double f3;
    double f4;
    double f5;
    double f6;

    double P1;
    double P2;
    double P3;
    double P4;
    double P5;
    double P6;

    double CD = 0.9;
    double CY = 0.9;
    double CL = 0.9;
    double Cl = 0.9;
    double Cm = 0.9;
    double Cn = 0.9;

    double coefB1;
    double coefB2;
    double coefB3;
    double coefB4;
    double coefB5;
    double coefB6;

    double A1;
    double A2;
    double A3;
    double A4;
    double A5;
    double A6;

    double G1;
    double G2;
    double G3;
    double G4;
    double G5;
    double G6;

    double aux_differential_equation[mStates];


    // c.3) ANFIS
    printf("\n\t\tc.3) ANFIS workspace...");

    fill((double*)XX, nInputs, 1, 50.0);
    fill((double*)ukk, mInputs, nData, 50.0);
    fill((double*)dukk, mInputs * Nc, 1, 0.0);
    fill((double*)xkk, mStates, nData, 50.0);
    fill((double*)ykk, mStates, nData, 50.0);
    fill((double*)yee, nData, 1, 50.0);

    // d) Reference definition
    printf("\n\td) Reference definition...");
    fill((double*)ref, nData, 1, 0.0);
    for (int i = 0; i < nData; i++) {
        ref[i] = reference(tt[i]) * 2.0 - 50;
    }

    // e) Reference filter
    printf("\n\te) Filtering reference...");
    fill((double*)wref, nData, 1, 0.0);
    wref[0] = ref[0];
    for (int i = 1; i < nData; i++) {
        wref[i] = alpha * wref[i - 1] + (1.0 - alpha) * ref[i];
    }

    // f) Reference plot
    printf("\n\tf) Plotting part VI...");

    // g) Index table
    printf("\n\tg) Index table...");
    printf("\n\t\td.2) Index table");
    static int indexTable[nRules][nInputs];
    for (int k = 1; k <= nInputs; k++) {
        int l = 1;
        for (int j = 1; j <= nRules; j = j + (int)pow(nFuzzy, long long(k) - 1)) {
            for (int i = 1; i <= (int)pow(nFuzzy, (long long(k) - 1)); i++) {
                indexTable[j + i - 2][nInputs - k] = l;
            }
            l = l + 1;
            if (l > nFuzzy) {
                l = 1;
            }
        }
    }
    /*
          for(int i=1; i<=nRules;i++){
                printf("\n");
            for(int j=1;j<=nInputs;j++){
                printf(" %d", indexTable[i-1][j-1]);
            }
        }
          */

          // Step 4: Control sequence
    printf("\n\nStep 4: Control sequence starts...");


    double clk_nData = 0.0;
    double clk_blimp = 0.0;
    double clk_ANFIS = 0.0;
    double clk_prediction = 0.0;
    double clk_Yfree = 0.0;
    double clk_Ystep = 0.0;
    double clk_duk = 0.0;
    clock_t clock_nData1;
    clock_t clock_nData2;
    clock_t clock_blimp1;
    clock_t clock_blimp2;
    clock_t clock_ANFIS1;
    clock_t clock_ANFIS2;
    clock_t clock_prediction1;
    clock_t clock_prediction2;
    clock_t clock_Yfree1;
    clock_t clock_Yfree2;
    clock_t clock_Ystep1;
    clock_t clock_Ystep2;
    clock_t clock_duk1;
    clock_t clock_duk2;

    for (int n = 3; n < nData - Np; n++) {
        clock_nData1 = clock();
        // a) OUTPUT reading
        // Blimp (yk)

        clock_blimp1 = clock();
        // a) Dynamics vector, Fd
        f1 = -blimp_mz * X[2][n - 1] * X[4][n - 1] + blimp_my * X[5][n - 1] * X[1][n - 1] + blimp_mTotal * (blimp_ax * (X[4][n - 1] * X[4][n - 1] + X[5][n - 1] * X[5][n - 1]) - blimp_az * X[5][n - 1] * X[3][n - 1]);
        f2 = -blimp_mx * X[0][n - 1] * X[5][n - 1] + blimp_mz * X[3][n - 1] * X[2][n - 1] + blimp_mTotal * (-blimp_ax * X[3][n - 1] * X[4][n - 1] - blimp_az * X[5][n - 1] * X[4][n - 1]);
        f3 = -blimp_my * X[1][n - 1] * X[3][n - 1] + blimp_mx * X[4][n - 1] * X[0][n - 1] + blimp_mTotal * (-blimp_ax * X[5][n - 1] * X[3][n - 1] + blimp_az * (X[4][n - 1] * X[4][n - 1] + X[3][n - 1] * X[3][n - 1]));
        f4 = -(blimp_Jz - blimp_Jy) * X[5][n - 1] * X[4][n - 1] + blimp_Jxz * X[3][n - 1] * X[4][n - 1] + blimp_mTotal * blimp_az * (X[0][n - 1] * X[5][n - 1] - X[3][n - 1] * X[2][n - 1]);
        f5 = -(blimp_Jx - blimp_Jz) * X[3][n - 1] * X[5][n - 1] + blimp_Jxz * (X[5][n - 1] * X[5][n - 1] - X[3][n - 1] * X[3][n - 1]) + blimp_mTotal * (blimp_ax * (X[1][n - 1] * X[3][n - 1] - X[4][n - 1] * X[0][n - 1]) - blimp_az * (X[2][n - 1] * X[4][n - 1] - X[5][n - 1] * X[1][n - 1]));
        f6 = -(blimp_Jy - blimp_Jx) * X[4][n - 1] * X[3][n - 1] - blimp_Jxz * X[4][n - 1] * X[5][n - 1] + blimp_mTotal * (-blimp_ax * (X[0][n - 1] * X[5][n - 1] - X[3][n - 1] * X[2][n - 1]));
        Fd[0][n - 1] = f1;
        Fd[1][n - 1] = f2;
        Fd[2][n - 1] = f3;
        Fd[3][n - 1] = f4;
        Fd[4][n - 1] = f5;
        Fd[5][n - 1] = f6;

        // b) Propulsion vector, P
        U[0][n - 1] = UT[0][n - 1] * UT[1][n - 1];            // Alpha* Tmax
        U[1][n - 1] = UT[0][n - 1] * (1.0 - UT[1][n - 1]);      // (1 - Alpha)* Tmax
        U[2][n - 1] = UT[2][n - 1];

        P1 = (U[0][n - 1] + U[1][n - 1]) * cos(U[2][n - 1]);
        P2 = 0;
        P3 = -(U[0][n - 1] + U[1][n - 1]) * sin(U[2][n - 1]);
        P4 = (U[1][n - 1] - U[0][n - 1]) * sin(U[2][n - 1]) * blimp_dy;
        P5 = (U[0][n - 1] + U[1][n - 1]) * (blimp_dz * cos(U[2][n - 1]) - blimp_dx * sin(U[2][n - 1]));
        P6 = (U[1][n - 1] - U[0][n - 1]) * cos(U[2][n - 1]) * blimp_dy;
        P[0][n - 1] = P1;
        P[1][n - 1] = P2;
        P[2][n - 1] = P3;
        P[3][n - 1] = P4;
        P[4][n - 1] = P5;
        P[5][n - 1] = P6;

        // c) Aerodynamic force vector, A
        Vtot[n - 1] = pow(X[0][n - 1] * X[0][n - 1] + X[1][n - 1] * X[1][n - 1] + X[2][n - 1] * X[2][n - 1], 0.5);
        Wtot[n - 1] = pow(X[3][n - 1] * X[3][n - 1] + X[4][n - 1] * X[4][n - 1] + X[5][n - 1] * X[5][n - 1], 0.5);

        coefB1 = 0.5 * rhoAir * X[0][n - 1] * X[0][n - 1] * sign(X[0][n - 1]) * blimp_area;
        coefB2 = 0.5 * rhoAir * X[1][n - 1] * X[1][n - 1] * sign(X[1][n - 1]) * blimp_area;
        coefB3 = 0.5 * rhoAir * X[2][n - 1] * X[2][n - 1] * sign(X[2][n - 1]) * blimp_area;
        coefB4 = 0.5 * rhoAir * X[3][n - 1] * X[3][n - 1] * sign(X[3][n - 1]) * blimp_volume;
        coefB5 = 0.5 * rhoAir * X[4][n - 1] * X[4][n - 1] * sign(X[4][n - 1]) * blimp_volume;
        coefB6 = 0.5 * rhoAir * X[5][n - 1] * X[5][n - 1] * sign(X[5][n - 1]) * blimp_volume;

        A1 = -CD * coefB1;
        A2 = -CY * coefB2;
        A3 = -CL * coefB3;
        A4 = -Cl * coefB4;
        A5 = -Cm * coefB5;
        A6 = -Cn * coefB6;

        A[0][n - 1] = A1;
        A[1][n - 1] = A2;
        A[2][n - 1] = A3;
        A[3][n - 1] = A4;
        A[4][n - 1] = A5;
        A[5][n - 1] = A6;

        // d) Gravitational force vector, G
        lambda[0][0] = cos(Y[4][n - 1]) * cos(Y[5][n - 1]);
        lambda[0][1] = cos(Y[4][n - 1]) * sin(Y[5][n - 1]);
        lambda[0][2] = sin(Y[4][n - 1]);
        lambda[1][0] = (-cos(Y[3][n - 1]) * sin(Y[5][n - 1]) + sin(Y[3][n - 1]) * sin(Y[4][n - 1]) * cos(Y[5][n - 1]));
        lambda[1][1] = (cos(Y[3][n - 1]) * cos(Y[5][n - 1]) + sin(Y[3][n - 1]) * sin(Y[4][n - 1]) * sin(Y[5][n - 1]));
        lambda[1][2] = sin(Y[3][n - 1]) * cos(Y[4][n - 1]);
        lambda[2][0] = (sin(Y[3][n - 1]) * sin(Y[5][n - 1]) + cos(Y[3][n - 1]) * sin(Y[4][n - 1]) * cos(Y[5][n - 1]));
        lambda[2][1] = (-sin(Y[3][n - 1]) * cos(Y[5][n - 1]) + cos(Y[3][n - 1]) * sin(Y[4][n - 1]) * sin(Y[5][n - 1]));
        lambda[2][2] = cos(Y[3][n - 1]) * cos(Y[4][n - 1]);

        double B = rhoAir * g_acc * blimp_volume;
        double W = blimp_mTotal * g_acc;

        G1 = lambda[2][0] * (W - B);
        G2 = lambda[2][1] * (W - B);
        G3 = lambda[2][2] * (W - B);
        G4 = -lambda[2][1] * blimp_az * W;
        G5 = (lambda[2][0] * blimp_az - lambda[2][2] * blimp_ax) * W;
        G6 = lambda[2][1] * blimp_ax * W;

        G[0][n - 1] = G1;
        G[1][n - 1] = G2;
        G[2][n - 1] = G3;
        G[3][n - 1] = G4;
        G[4][n - 1] = G5;
        G[5][n - 1] = G6;

        // e) Differential equation
        for (int i = 0; i < mStates; i++) {
            aux_differential_equation[i] = P[i][n - 1] + Fd[i][n - 1] + A[i][n - 1] + G[i][n - 1];
        }

        for (int i = 0; i < mStates; i++) {
            for (int j = 0; j < mStates; j++) {
                dX[i][n - 1] = dX[i][n - 1] + invM[i][j] * aux_differential_equation[j];
            }
        }

        // f) Integrate differential equation
        for (int i = 0; i < mStates; i++) {
            for (int j = 0; j < n; j++) {
                X[i][n] = X[i][n] + dX[i][j];
            }
            X[i][n] = X[i][n] + (dX[i][n - 1] - dX[i][0]) * 0.5;
            X[i][n] = X[i][n] * step;
        }

        // g) Calculate vehicle position in terms of displacements in the north, east and vertical directions 
        for (int i = 0; i < mStates; i++) {
            for (int j = 0; j < n; j++) {
                Y[i][n] = Y[i][n] + X[i][j];
            }
            Y[i][n] = Y[i][n] + (X[i][n - 1] - X[i][0]) * 0.5;
            Y[i][n] = Y[i][n] * step;                           // XYZ aligned with NEU
        }


        clock_blimp2 = clock();
        // Scaling to fuzzy
        xkk[0][n] = map(X[0][n], minX1, maxX1, minOut, maxOut);
        xkk[1][n] = map(X[1][n], minX2, maxX2, minOut, maxOut);
        xkk[2][n] = map(X[2][n], minX3, maxX3, minOut, maxOut);
        xkk[3][n] = map(X[3][n], minX4, maxX4, minOut, maxOut);
        xkk[4][n] = map(X[4][n], minX5, maxX5, minOut, maxOut);
        xkk[5][n] = map(X[5][n], minX6, maxX6, minOut, maxOut);

        ykk[0][n] = map(Y[0][n], minY1, maxY1, minOut, maxOut);
        ykk[1][n] = map(Y[1][n], minY2, maxY2, minOut, maxOut);
        ykk[2][n] = map(Y[2][n], minY3, maxY3, minOut, maxOut);
        ykk[3][n] = map(Y[3][n], minY4, maxY4, minOut, maxOut);
        ykk[4][n] = map(Y[4][n], minY5, maxY5, minOut, maxOut);
        ykk[5][n] = map(Y[5][n], minY6, maxY6, minOut, maxOut);

        // b) ANFIS OUTPUT reading
        clock_ANFIS1 = clock();
        // Prelayer
        XX[0] = ukk[0][n - 1];
        XX[1] = ukk[1][n - 1];
        XX[2] = ukk[2][n - 1];
        XX[3] = xkk[0][n - 1];
        // Layer 1: Input fuzzyfication
        for (int k = 0; k < nInputs; k++) {
            for (int j = 0; j < nFuzzy; j++) {
                muIn[j][k] = gauss(XX[k], aIne[j][k], cIne[j][k]);
                if (muIn[j][k] < tol) {
                    muIn[j][k] = tol;
                }
                if (muIn[j][k] > (1.0 - tol)) {
                    muIn[j][k] = 1.0 - tol;
                }
            }
        }
        // Layer 2: Calculation weights
        for (int j = 0; j < nRules; j++) {
            w[j] = 1.0;
            for (int k = 0; k < nInputs; k++) {
                w[j] = w[j] * muIn[indexTable[j][k] - 1][k];
            }
        }
        // Layer 3: Normalizing
        sumW = 0.0;
        for (int j = 0; j < nRules; j++) {
            sumW = sumW + w[j];
        }
        if (sqrt(sumW * sumW) < tol) {
            sumW = tol;
        }
        for (int j = 0; j < nRules; j++) {
            wn[j] = w[j] / sumW;
        }
        // Layer 4: Calculation of wn*fi
        for (int j = 0; j < nRules; j++) {
            fi[j] = invSigmoid(w[j], aOute[j], cOute[j]);
            fi_wn[j] = fi[j] * wn[j];
        }
        // Layer 5: Calculating output
        double f = 0;
        for (int j = 0; j < nRules; j++) {
            f = f + fi_wn[j];
        }
        yee[n] = f;

        clock_ANFIS2 = clock();


        // Wref
        for (int j = 0; j < Np; j++) {
            Wref[j] = wref[n];
        }

        // Error
        double ekk = wref[n] - xkk[0][n];

        // c) Prediction
        clock_prediction1 = clock();
        // ANFIS(Yfree)
        clock_Yfree1 = clock();
        Yfree[0] = xkk[0][n - 1];
        for (int z = 0; z < Np - 1; z++) {
            // Prelayer
            XX[0] = ukk[0][n - 1];
            XX[1] = ukk[1][n - 1];
            XX[2] = ukk[2][n - 1];
            XX[3] = Yfree[z];
            // Layer 1: Input fuzzyfication
            for (int k = 0; k < nInputs; k++) {
                for (int j = 0; j < nFuzzy; j++) {
                    muIn[j][k] = gauss(XX[k], aIne[j][k], cIne[j][k]);
                    if (muIn[j][k] < tol) {
                        muIn[j][k] = tol;
                    }
                    if (muIn[j][k] > (1.0 - tol)) {
                        muIn[j][k] = 1.0 - tol;
                    }
                }
            }
            // Layer 2: Calculation weights
            for (int j = 0; j < nRules; j++) {
                w[j] = 1.0;
                for (int k = 0; k < nInputs; k++) {
                    w[j] = w[j] * muIn[indexTable[j][k] - 1][k];
                }
            }
            // Layer 3: Normalizing
            sumW = 0.0;
            for (int j = 0; j < nRules; j++) {
                sumW = sumW + w[j];
            }
            if (sqrt(sumW * sumW) < tol) {
                sumW = tol;
            }
            for (int j = 0; j < nRules; j++) {
                wn[j] = w[j] / sumW;
            }
            // Layer 4: Calculation of wn*fi
            for (int j = 0; j < nRules; j++) {
                fi[j] = invSigmoid(w[j], aOute[j], cOute[j]);
                fi_wn[j] = fi[j] * wn[j];
            }
            // Layer 5: Calculating output
            double f = 0.0;
            for (int j = 0; j < nRules; j++) {
                f = f + fi_wn[j];
            }
            Yfree[z + 1] = f;
        }
        clock_Yfree2 = clock();
        clock_Ystep1 = clock();

        // ANFIS(Ystep)
            // ukk1
        if (ukk[0][n - 1] > 50) {
            ddu[0] = -35;
        }
        else {
            ddu[0] = +35;
        }

        // ANFIS
        for (int z = 0; z < Np; z++) {
            // Prelayer
            if (z == 0) {
                XX[3] = xkk[0][n - 1];
            }
            else {
                XX[3] = Ystep[z - 1][0];
            }
            XX[0] = ukk[0][n - 1] + ddu[0];
            XX[1] = ukk[1][n - 1];
            XX[2] = ukk[2][n - 1];
            // Layer 1: Input fuzzyfication
            for (int k = 0; k < nInputs; k++) {
                for (int j = 0; j < nFuzzy; j++) {
                    muIn[j][k] = gauss(XX[k], aIne[j][k], cIne[j][k]);
                    if (muIn[j][k] < tol) {
                        muIn[j][k] = tol;
                    }
                    if (muIn[j][k] > (1.0 - tol)) {
                        muIn[j][k] = 1.0 - tol;
                    }
                }
            }
            // Layer 2: Calculation weights
            for (int j = 0; j < nRules; j++) {
                w[j] = 1.0;
                for (int k = 0; k < nInputs; k++) {
                    w[j] = w[j] * muIn[indexTable[j][k] - 1][k];
                }
            }
            // Layer 3: Normalizing
            sumW = 0.0;
            for (int j = 0; j < nRules; j++) {
                sumW = sumW + w[j];
            }
            if (sqrt(sumW * sumW) < tol) {
                sumW = tol;
            }
            for (int j = 0; j < nRules; j++) {
                wn[j] = w[j] / sumW;
            }
            // Layer 4: Calculation of wn*fi
            for (int j = 0; j < nRules; j++) {
                fi[j] = invSigmoid(w[j], aOute[j], cOute[j]);
                fi_wn[j] = fi[j] * wn[j];
            }
            // Layer 5: Calculating output
            double f = 0.0;
            for (int j = 0; j < nRules; j++) {
                f = f + fi_wn[j];
            }
            Ystep[z][0] = f;

        }



        // ukk2
        if (ukk[1][n - 1] > 50) {
            ddu[1] = -35;
        }
        else {
            ddu[1] = +35;
        }

        // ANFIS
        for (int z = 0; z < Np; z++) {
            // Prelayer
            if (z == 0) {
                XX[3] = xkk[0][n - 1];
            }
            else {
                XX[3] = Ystep[z - 1][1];
            }
            XX[0] = ukk[0][n - 1];
            XX[1] = ukk[1][n - 1] + ddu[1];
            XX[2] = ukk[2][n - 1];
            // Layer 1: Input fuzzyfication
            for (int k = 0; k < nInputs; k++) {
                for (int j = 0; j < nFuzzy; j++) {
                    muIn[j][k] = gauss(XX[k], aIne[j][k], cIne[j][k]);
                    if (muIn[j][k] < tol) {
                        muIn[j][k] = tol;
                    }
                    if (muIn[j][k] > (1.0 - tol)) {
                        muIn[j][k] = 1.0 - tol;
                    }
                }
            }
            // Layer 2: Calculation weights
            for (int j = 0; j < nRules; j++) {
                w[j] = 1.0;
                for (int k = 0; k < nInputs; k++) {
                    w[j] = w[j] * muIn[indexTable[j][k] - 1][k];
                }
            }
            // Layer 3: Normalizing
            sumW = 0.0;
            for (int j = 0; j < nRules; j++) {
                sumW = sumW + w[j];
            }
            if (sqrt(sumW * sumW) < tol) {
                sumW = tol;
            }
            for (int j = 0; j < nRules; j++) {
                wn[j] = w[j] / sumW;
            }
            // Layer 4: Calculation of wn*fi
            for (int j = 0; j < nRules; j++) {
                fi[j] = invSigmoid(w[j], aOute[j], cOute[j]);
                fi_wn[j] = fi[j] * wn[j];
            }
            // Layer 5: Calculating output
            double f = 0.0;
            for (int j = 0; j < nRules; j++) {
                f = f + fi_wn[j];
            }
            Ystep[z][1] = f;
        }


        // ukk3
        if (ukk[2][n - 1] > 50) {
            ddu[2] = -35;
        }
        else {
            ddu[2] = +35;
        }

        // ANFIS
        for (int z = 0; z < Np; z++) {
            // Prelayer
            if (z == 0) {
                XX[3] = xkk[0][n - 1];
            }
            else {
                XX[3] = Ystep[z - 1][2];
            }
            XX[0] = ukk[0][n - 1];
            XX[1] = ukk[1][n - 1];
            XX[2] = ukk[2][n - 1] + ddu[2];
            // Layer 1: Input fuzzyfication
            for (int k = 0; k < nInputs; k++) {
                for (int j = 0; j < nFuzzy; j++) {
                    muIn[j][k] = gauss(XX[k], aIne[j][k], cIne[j][k]);
                    if (muIn[j][k] < tol) {
                        muIn[j][k] = tol;
                    }
                    if (muIn[j][k] > (1.0 - tol)) {
                        muIn[j][k] = 1.0 - tol;
                    }
                }
            }
            // Layer 2: Calculation weights
            for (int j = 0; j < nRules; j++) {
                w[j] = 1.0;
                for (int k = 0; k < nInputs; k++) {
                    w[j] = w[j] * muIn[indexTable[j][k] - 1][k];
                }
            }
            // Layer 3: Normalizing
            sumW = 0.0;
            for (int j = 0; j < nRules; j++) {
                sumW = sumW + w[j];
            }
            if (sqrt(sumW * sumW) < tol) {
                sumW = tol;
            }
            for (int j = 0; j < nRules; j++) {
                wn[j] = w[j] / sumW;
            }
            // Layer 4: Calculation of wn*fi
            for (int j = 0; j < nRules; j++) {
                fi[j] = invSigmoid(w[j], aOute[j], cOute[j]);
                fi_wn[j] = fi[j] * wn[j];
            }
            // Layer 5: Calculating output
            double f = 0.0;
            for (int j = 0; j < nRules; j++) {
                f = f + fi_wn[j];
            }
            Ystep[z][2] = f;
        }

        clock_Ystep2 = clock();

        // Gij
        for (int nn = 0; nn < mInputs; nn++) {
            for (int j = 0; j < Np; j++) {
                gi[j] = (Ystep[j][nn] - Yfree[j]) / ddu[nn];
            }
            for (int j = 0; j < Np; j++) {
                for (int k = 0; k < Nc; k++) {
                    if (k > j) {
                        Gij[j][k + (nn)*Nc] = 0.0;
                    }
                    else {
                        Gij[j][k + (nn)*Nc] = gi[j - k];
                    }
                }
            }
        }




        // d) Optimal control sequence calculation
        clock_duk1 = clock();

        fill((double*)Gijt, mInputs * Nc, Np, 0.0);
        fill((double*)Gijt_Q, mInputs * Nc, Np, 0.0);
        fill((double*)Gijt_Q_Gij, mInputs * Nc, mInputs * Nc, 0.0);
        fill((double*)Gijt_Q_Gij_R, mInputs * Nc, mInputs * Nc, 0.0);
        fill((double*)Gijt_Q_Gij_R_1, mInputs * Nc, mInputs * Nc, 0.0);
        fill((double*)Gijt_Q_Gij_R_1_Gijt, mInputs * Nc, Np, 0.0);
        fill((double*)Gijt_Q_Gij_R_1_Gijt_Q, mInputs * Nc, Np, 0.0);
        fill((double*)Wref_Yfree, Np, 1, 0.0);
        fill((double*)dukk, mInputs * Nc, 1, 0.0);

        matrixTranspose(Np, (double*)Gij, mInputs * Nc, (double*)Gijt);
        matrixMultiplication(mInputs * Nc, (double*)Gijt, Np, (double*)Q, Np, (double*)Gijt_Q);
        matrixMultiplication(mInputs * Nc, (double*)Gijt_Q, Np, (double*)Gij, mInputs * Nc, (double*)Gijt_Q_Gij);
        matrixSum(mInputs * Nc, mInputs * Nc, (double*)Gijt_Q_Gij, (double*)R, (double*)Gijt_Q_Gij_R);

        mINVERSE((double*)Gijt_Q_Gij_R, (double*)Gijt_Q_Gij_R_1);

        matrixMultiplication(mInputs * Nc, (double*)Gijt_Q_Gij_R_1, mInputs * Nc, (double*)Gijt, Np, (double*)Gijt_Q_Gij_R_1_Gijt);
        matrixMultiplication(mInputs * Nc, (double*)Gijt_Q_Gij_R_1_Gijt, Np, (double*)Q, Np, (double*)Gijt_Q_Gij_R_1_Gijt_Q);
        matrixSubstract(Np, 1, (double*)Wref, (double*)Yfree, (double*)Wref_Yfree);
        matrixMultiplication(mInputs * Nc, (double*)Gijt_Q_Gij_R_1_Gijt_Q, Np, (double*)Wref_Yfree, 1, (double*)dukk);

        clock_duk2 = clock();

        /*
        printf("\nGijt_Q_Gij:");
        for (int i = 1; i <= mInputs * Nc; i++) {
            printf("\n");
            for (int j = 1; j <= Nc*mInputs; j++) {
                printf("Gijt_Q_Gij(%d,%d) = %06.3f; ", i, j, Gijt_Q_Gij[i - 1][j - 1]);
            }
        }

        printf("\nR:");
        for (int i = 1; i <= Nc*mInputs; i++) {
            printf("\n");
            for (int j = 1; j <= Nc*mInputs; j++) {
                printf("R(%d,%d)=%06.3f; ", i, j, R[i - 1][j - 1]);
            }
        }

        printf("\nGijt_Q_Gij x R:");
        for (int i = 0; i < mInputs*Nc; i++) {
            printf("\n");
            for (int j = 0; j < mInputs*Nc; j++) {
                printf("[%06.3f] ", Gijt_Q_Gij_R[i][j]);
            }
        }*/

        // e) Determination of ut
        // ukk(fuzzy)
        ukk[0][n] = ukk[0][n - 1] + Pu1 * dukk[0];
        ukk[1][n] = ukk[1][n - 1] + Pu2 * dukk[Nc];
        ukk[2][n] = ukk[2][n - 1] + Pu3 * dukk[2 * Nc];

        // Limiting control action
        if (ukk[0][n] < minOut) {
            ukk[0][n] = minOut;
        }
        if (ukk[1][n] < minOut) {
            ukk[1][n] = minOut;
        }
        if (ukk[2][n] < minOut) {
            ukk[2][n] = minOut;
        }
        if (ukk[0][n] > maxOut) {
            ukk[0][n] = maxOut;
        }
        if (ukk[1][n] > maxOut) {
            ukk[1][n] = maxOut;
        }
        if (ukk[2][n] > maxOut) {
            ukk[2][n] = maxOut;
        }

        clock_prediction2 = clock();


        // Scaling to Force(N)
        UT[0][n] = map(ukk[0][n], minOut, maxOut, minUT1, maxUT1);
        UT[1][n] = map(ukk[1][n], minOut, maxOut, minUT2, maxUT2);
        UT[2][n] = map(ukk[2][n], minOut, maxOut, minUT3, maxUT3);

        clock_nData2 = clock();
        // printf

        clk_blimp = clk_blimp + (double(clock_blimp2) - double(clock_blimp1));
        clk_ANFIS = clk_ANFIS + (double(clock_ANFIS2) - double(clock_ANFIS1));
        clk_prediction = clk_prediction + (double(clock_prediction2) - double(clock_prediction1));
        clk_nData = clk_nData + (double(clock_nData2) - double(clock_nData1));
        clk_Yfree = clk_Yfree + (double(clock_Yfree2) - double(clock_Yfree1));
        clk_Ystep = clk_Ystep + (double(clock_Ystep2) - double(clock_Ystep1));
        clk_duk = clk_duk + (double(clock_duk2) - double(clock_duk1));
        //printf("\nnData = %04d. Blimp dynamics = %10.10f s. ANFIS network = %10.10f s. Prediction = %10.10f s. Overall = %10.10f s.", n, (double(clock_blimp2) - double(clock_blimp1)) / CLOCKS_PER_SEC, (double(clock_ANFIS2) - double(clock_ANFIS1)) / CLOCKS_PER_SEC, (double(clock_prediction2) - double(clock_prediction1)) / CLOCKS_PER_SEC, (double(clock_nData2) - double(clock_nData1)) / CLOCKS_PER_SEC);
    }

    printf("\n\tControl simulation completed in %.3f seconds.", clk_nData / (double(CLOCKS_PER_SEC)));
    printf("\n\tShowing average speeds:");
    printf("\n\t\ta) Overall step calculation: \t\t%.9f s", clk_nData / (double(CLOCKS_PER_SEC) * (nData - Np - 3)));
    printf("\n\t\tb) Blimp dynamics calculation: \t\t%.9f s", clk_blimp / (double(CLOCKS_PER_SEC) * (nData - Np - 3)));
    printf("\n\t\tc) C++ ANFIS calculation: \t\t%.9f s", clk_ANFIS / (double(CLOCKS_PER_SEC) * (nData - Np - 3)));
    printf("\n\t\td) Overall prediction calculation: \t%.9f s", clk_prediction / (double(CLOCKS_PER_SEC) * (nData - Np - 3)));
    printf("\n\t\t\td.1) Yfree calculation: \t  %.9f s", clk_Yfree / (double(CLOCKS_PER_SEC) * (nData - Np - 3)));
    printf("\n\t\t\td.2) Ystep calculation: \t  %.9f s", clk_Ystep / (double(CLOCKS_PER_SEC) * (nData - Np - 3)));
    printf("\n\t\t\td.3) Optimizer (duk): \t\t  %.9f s", clk_duk / (double(CLOCKS_PER_SEC) * (nData - Np - 3)));


    // Step 5: Control results
    printf("\n\nStep 5: Control results");
    printf("\n\ta) Writing to ANFISblimpCtrl.m file...");

    FILE* blimpFile;
    blimpFile = fopen("ANFISblimpCtrl.m", "w");
    if (blimpFile != NULL) {
        for (int j = 1; j <= nData; j++) {
            fprintf(blimpFile, "\nttt(%d)\t= %.15f;", j, tt[j - 1]);
        }
        for (int i = 1; i <= mInputs; i++) {
            for (int j = 1; j <= nData; j++) {
                fprintf(blimpFile, "\nukk(%d,%d)\t= %.15f;", i, j, ukk[i - 1][j - 1]);
            }
        }
        for (int i = 1; i <= mStates; i++) {
            for (int j = 1; j <= nData; j++) {
                fprintf(blimpFile, "\nxkk(%d,%d)\t= %.15f;", i, j, xkk[i - 1][j - 1]);
            }
        }
        for (int i = 1; i <= mStates; i++) {
            for (int j = 1; j <= nData; j++) {
                fprintf(blimpFile, "\nykk(%d,%d)\t= %.15f;", i, j, ykk[i - 1][j - 1]);
            }
        }
        for (int j = 1; j <= nData; j++) {
            fprintf(blimpFile, "\nwref(%d)\t= %.15f;", j, wref[j - 1]);
        }
        for (int j = 1; j <= nData; j++) {
            fprintf(blimpFile, "\nyee(%d)\t= %.15f;", j, yee[j - 1]);
        }

        printf("\n\t¡Success!");
    }
    fclose(blimpFile);

    return 0;
}




double h(double t) {
    if (t >= 0.0) {
        return 1.0;
    }
    else {
        return 0.0;
    }
}

double r(double t) {
    return h(t) * t;
}

void fill(double* arr, int n, int m, double val) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            *(arr + (long long(i) * long long(m) + long long(j))) = val;
        }
    }
}

double sign(double x) {
    if (x >= 0.0) {
        return 1.0;
    }
    else {
        return -1.0;
    }
}

double map(double input, double minIn, double maxIn, double minOut, double maxOut) {
    return ((input - minIn) / (maxIn - minIn)) * (maxOut - minOut) + minOut;
}

double reference(double t) {
    return h(t) * 50 - h(t - 30) * 10 + h(t - 80) * 10 + h(t - 130) * 10 - h(t - 180) * 2 - r(t - 230) * 0.1 + r(t - 260) * 0.1 - h(t - 290) * 10 - h(t - 320) * 2 + r(t - 350) * 0.1 - r(t - 450) * 0.1;
}

double gauss(double x, double a, double c) {
    return exp(-(x - c) * (x - c) / (a * a));
}

double sigmoid(double x, double a, double c) {
    return 1.0 / (1.0 + exp(-(x - c) / a));
}

double invSigmoid(double x, double a, double c) {
    return  c - a * log(1.0 / x - 1.0);
}

double dGauss_da(double x, double a, double c) {
    return (2.0 * exp(-(-c + x) * (-c + x) / (a * a)) * (-c + x) * (-c + x)) / (a * a * a);
}

double dGauss_dc(double x, double a, double c) {
    return (2.0 * exp(-(-c + x) * (-c + x) / (a * a)) * (-c + x)) / (a * a);
}

double dinvSigmoid_da(double x, double a, double c) {
    return -log(1.0 / x - 1.0);
}

double dinvSigmoid_dc(double x, double a, double c) {
    return 1.0;
}

void matrixMultiplication(int a, double* A, int c, double* B, int b, double* C) {
    // C = A * B
    fill((double*)C, a, b, 0.0);
    for (int i = 0; i < a; i++) {
        for (int j = 0; j < b; j++) {
            for (int k = 0; k < c; k++) {
                *(C + long long(i) * long long(b) + long long(j)) = *(C + long long(i) * long long(b) + long long(j)) + *(A + long long(i) * long long(c) + long long(k)) * *(B + long long(k) * long long(c) + long long(j));
            }
        }
    }
}

void matrixTranspose(int a, double* A, int b, double* B) {
    // B = A'
    fill((double*)B, b, a, 0.0);
    for (int i = 0; i < a; i++) {
        for (int j = 0; j < b; j++) {
            *(B + long long(j) * long long(a) + long long(i)) = *(A + long long(i) * long long(b) + long long(j));
        }
    }
}

void matrixSum(int n, int m, double* A, double* B, double* C) {
    // C = A + B
    fill((double*)C, n, m, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            *(C + long long(i) * long long(m) + long long(j)) = *(A + long long(i) * long long(m) + long long(j)) + *(B + long long(i) * long long(m) + long long(j));
        }
    }
}

void matrixSubstract(int n, int m, double* A, double* B, double* C) {
    // C = A - B
    fill((double*)C, n, m, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            *(C + long long(i) * long long(m) + long long(j)) = *(A + long long(i) * long long(m) + long long(j)) - *(B + long long(i) * long long(m) + long long(j));
        }
    }
}

void getCfactor(double M[N][N], double t[N][N], int p, int q, int n) {
    int i = 0, j = 0;
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            if (r != p && c != q) {
                t[i][j++] = M[r][c];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double DET(double M[N][N], int n) {
    double D = 0;
    if (n == 1) {
        return M[0][0];
    }
    double t[N][N];
    int s = 1;
    for (int f = 0; f < n; f++) {
        getCfactor(M, t, 0, f, n);
        D += s * M[0][f] * DET(t, n - 1);
        s = -s;
    }
    return D;
}

void ADJ(double M[N][N], double adj[N][N]) {
    if (N == 1) {
        adj[0][0] = 1;
        return;
    }
    int s = 1;
    double t[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            getCfactor(M, t, i, j, N);
            s = ((i + j) % 2 == 0) ? 1 : -1;
            adj[j][i] = (s) * (DET(t, N - 1));
        }
    }
}

bool INV(double M[N][N], double inv[N][N]) {
    double det = DET(M, N);
    if (det == 0) {
        printf("can't find its inverse");
        return false;
    }
    double adj[N][N];
    ADJ(M, adj);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            inv[i][j] = adj[i][j] / double(det);
        }
    }
    return true;
}

void mINVERSE(double M[81], double inverse[81])
{
    double x[81];
    double smax;
    int b_i;
    int i;
    int i2;
    int ix;
    int iy;
    int j;
    int jA;
    int jp1j;
    int k;
    signed char ipiv[9];
    signed char p[9];

    //  Inverse function
    for (i = 0; i < 81; i++) {
        inverse[i] = 0.0;
        x[i] = M[i];
    }

    for (i = 0; i < 9; i++) {
        ipiv[i] = static_cast<signed char>(i + 1);
    }

    for (j = 0; j < 8; j++) {
        int b_tmp;
        int mmj_tmp;
        mmj_tmp = 7 - j;
        b_tmp = j * 10;
        jp1j = b_tmp + 2;
        iy = 9 - j;
        jA = 0;
        ix = b_tmp;
        smax = fabs(x[b_tmp]);
        for (k = 2; k <= iy; k++) {
            double s;
            ix++;
            s = fabs(x[ix]);
            if (s > smax) {
                jA = k - 1;
                smax = s;
            }
        }

        if (x[b_tmp + jA] != 0.0) {
            if (jA != 0) {
                iy = j + jA;
                ipiv[j] = static_cast<signed char>(iy + 1);
                ix = j;
                for (k = 0; k < 9; k++) {
                    smax = x[ix];
                    x[ix] = x[iy];
                    x[iy] = smax;
                    ix += 9;
                    iy += 9;
                }
            }

            i = (b_tmp - j) + 9;
            for (b_i = jp1j; b_i <= i; b_i++) {
                x[b_i - 1] /= x[b_tmp];
            }
        }

        iy = b_tmp + 9;
        jA = b_tmp;
        for (jp1j = 0; jp1j <= mmj_tmp; jp1j++) {
            smax = x[iy];
            if (x[iy] != 0.0) {
                ix = b_tmp + 1;
                i = jA + 11;
                i2 = (jA - j) + 18;
                for (b_i = i; b_i <= i2; b_i++) {
                    x[b_i - 1] += x[ix] * -smax;
                    ix++;
                }
            }

            iy += 9;
            jA += 9;
        }
    }

    for (i = 0; i < 9; i++) {
        p[i] = static_cast<signed char>(i + 1);
    }

    for (k = 0; k < 8; k++) {
        signed char i1;
        i1 = ipiv[k];
        if (i1 > k + 1) {
            iy = p[i1 - 1];
            p[i1 - 1] = p[k];
            p[k] = static_cast<signed char>(iy);
        }
    }

    for (k = 0; k < 9; k++) {
        jp1j = 9 * (p[k] - 1);
        inverse[k + jp1j] = 1.0;
        for (j = k + 1; j < 10; j++) {
            i = (j + jp1j) - 1;
            if (inverse[i] != 0.0) {
                i2 = j + 1;
                for (b_i = i2; b_i < 10; b_i++) {
                    iy = (b_i + jp1j) - 1;
                    inverse[iy] -= inverse[i] * x[(b_i + 9 * (j - 1)) - 1];
                }
            }
        }
    }

    for (j = 0; j < 9; j++) {
        iy = 9 * j;
        for (k = 8; k >= 0; k--) {
            jA = 9 * k;
            i = k + iy;
            smax = inverse[i];
            if (smax != 0.0) {
                inverse[i] = smax / x[k + jA];
                for (b_i = 0; b_i < k; b_i++) {
                    jp1j = b_i + iy;
                    inverse[jp1j] -= inverse[i] * x[b_i + jA];
                }
            }
        }
    }
}








//-------------------------------| 
//        CUDA WRAPPERS          |
//-------------------------------| 





// END




void set_ANFIS() {

    aIne[0][0] = 15.000000000000000000000000000000;
    aIne[1][0] = 15.000000000000000000000000000000;
    aIne[2][0] = 15.000000000000000000000000000000;
    aIne[3][0] = 15.000000000000000000000000000000;
    aIne[4][0] = 15.000000000000000000000000000000;
    aIne[0][1] = 15.000000000000000000000000000000;
    aIne[1][1] = 15.000000000000000000000000000000;
    aIne[2][1] = 15.000000000000000000000000000000;
    aIne[3][1] = 15.000000000000000000000000000000;
    aIne[4][1] = 15.000000000000000000000000000000;
    aIne[0][2] = 15.000000000000000000000000000000;
    aIne[1][2] = 15.000000000000000000000000000000;
    aIne[2][2] = 15.000000000000000000000000000000;
    aIne[3][2] = 15.000000000000000000000000000000;
    aIne[4][2] = 15.000000000000000000000000000000;
    aIne[0][3] = 15.000000000000000000000000000000;
    aIne[1][3] = 15.000000000000000000000000000000;
    aIne[2][3] = 15.000000000000000000000000000000;
    aIne[3][3] = 15.000000000000000000000000000000;
    aIne[4][3] = 15.000000000000000000000000000000;
    cIne[0][0] = 0.000000000000000000000000000000;
    cIne[1][0] = 25.000000000000000000000000000000;
    cIne[2][0] = 50.000000000000000000000000000000;
    cIne[3][0] = 75.000000000000000000000000000000;
    cIne[4][0] = 100.000000000000000000000000000000;
    cIne[0][1] = 0.000000000000000000000000000000;
    cIne[1][1] = 25.000000000000000000000000000000;
    cIne[2][1] = 50.000000000000000000000000000000;
    cIne[3][1] = 75.000000000000000000000000000000;
    cIne[4][1] = 100.000000000000000000000000000000;
    cIne[0][2] = 0.000000000000000000000000000000;
    cIne[1][2] = 25.000000000000000000000000000000;
    cIne[2][2] = 50.000000000000000000000000000000;
    cIne[3][2] = 75.000000000000000000000000000000;
    cIne[4][2] = 100.000000000000000000000000000000;
    cIne[0][3] = 0.000000000000000000000000000000;
    cIne[1][3] = 25.000000000000000000000000000000;
    cIne[2][3] = 50.000000000000000000000000000000;
    cIne[3][3] = 75.000000000000000000000000000000;
    cIne[4][3] = 100.000000000000000000000000000000;
    aOute[0] = -2.084817645120403017955368341063;
    aOute[1] = -5.495744141622523670775990467519;
    aOute[2] = -5.311598485855140161504550633254;
    aOute[3] = -7.828101292786593390360394550953;
    aOute[4] = -10.080141135474692148932263080496;
    aOute[5] = 13.672808484567406495102659391705;
    aOute[6] = -13.176313192280673547429614700377;
    aOute[7] = -5.144128926688597225336252449779;
    aOute[8] = -11.840466607734260406914472696371;
    aOute[9] = -9.495792595406102876154363912065;
    aOute[10] = -2.700203491526818044832225496066;
    aOute[11] = -11.820597543278056562598976597656;
    aOute[12] = -34.735624016845051187374338041991;
    aOute[13] = -18.241571808841460011763047077693;
    aOute[14] = -18.000381680616275303918882855214;
    aOute[15] = -8.439672909591276805940651684068;
    aOute[16] = -18.040665212306397080510578234680;
    aOute[17] = -34.760083116886363541198079474270;
    aOute[18] = -34.760083116886363541198079474270;
    aOute[19] = -9.837327425481268150520008930471;
    aOute[20] = -11.514747143698901510333598707803;
    aOute[21] = -24.427426112197622387611772865057;
    aOute[22] = -34.760083116886363541198079474270;
    aOute[23] = -34.760083116886363541198079474270;
    aOute[24] = -12.004552723516388468283366819378;
    aOute[25] = -4.235680952643297025872470840113;
    aOute[26] = -10.697026580212176227746567747090;
    aOute[27] = -13.572702317246386982674266619142;
    aOute[28] = -10.225655929418566714161897834856;
    aOute[29] = -7.658881454173118719097601569956;
    aOute[30] = -9.872781077978263652994428412057;
    aOute[31] = -13.264654668891010302900212991517;
    aOute[32] = -17.516833574551569085997471120209;
    aOute[33] = -10.252512505138909659763157833368;
    aOute[34] = -14.476676270219252629090078698937;
    aOute[35] = 0.063453063298635592470731126014;
    aOute[36] = -6.574214063560877896463807701366;
    aOute[37] = -18.684634756470735794664506101981;
    aOute[38] = -20.487476057506640358951699454337;
    aOute[39] = -15.528277603167886411483777919784;
    aOute[40] = 1.594299461176739818313308205688;
    aOute[41] = -6.365054156309841282279649021802;
    aOute[42] = -7.584225829042147815073349192971;
    aOute[43] = -30.031819984220732067115022800863;
    aOute[44] = -10.038851707666241352967517741490;
    aOute[45] = -2.720405154799239610952099610586;
    aOute[46] = -7.123976765140384692642783193151;
    aOute[47] = -8.182905506504257076016983774025;
    aOute[48] = -34.189471888858719239578931592405;
    aOute[49] = -8.581865361336651787382834299933;
    aOute[50] = -16.314564480290574266518888180144;
    aOute[51] = -7.651148896944667932018546707695;
    aOute[52] = -2.020412648125384347252975203446;
    aOute[53] = -30.401970233808846444389928365126;
    aOute[54] = -22.222182903459263059176009846851;
    aOute[55] = -6.979409912431629336992955359165;
    aOute[56] = -3.449094763560766274679281195858;
    aOute[57] = -6.418385024024479967863499041414;
    aOute[58] = -34.412064111943564626017177943140;
    aOute[59] = -24.668675884871433368061843793839;
    aOute[60] = 15.835831923959625910924842173699;
    aOute[61] = -16.813192287148467585211619734764;
    aOute[62] = -22.113526369286358175259010749869;
    aOute[63] = -34.760083116886363541198079474270;
    aOute[64] = -20.070639046698001806134925573133;
    aOute[65] = 6.580061235148209952683373558102;
    aOute[66] = -4.532515603161148654010048630880;
    aOute[67] = -8.659041444838168288811175443698;
    aOute[68] = -11.662816366043800897500659630168;
    aOute[69] = -9.247108439690405745636780920904;
    aOute[70] = 3.112140441032947180133305664640;
    aOute[71] = -3.939738701016491706496935876203;
    aOute[72] = -7.677050192419641483354553201934;
    aOute[73] = -14.292570489393868982119784050155;
    aOute[74] = -7.551558679599160051054695941275;
    aOute[75] = 3.929209598645452050647008945816;
    aOute[76] = -6.698280528488928631247745215660;
    aOute[77] = -9.883218980514056895003704994451;
    aOute[78] = -34.537207832335255375255655962974;
    aOute[79] = -30.760760221321497454027849016711;
    aOute[80] = 11.837144424448691637508090934716;
    aOute[81] = 4.808758776898513609410201752326;
    aOute[82] = -9.171613997006781815457543416414;
    aOute[83] = -34.121628901414112533529987558722;
    aOute[84] = -32.061150329724590335445100208744;
    aOute[85] = 22.119515751912324219574657035992;
    aOute[86] = 5.102766733839051660481800354319;
    aOute[87] = -9.616619352807223464196795248426;
    aOute[88] = -21.463541795425435054767149267718;
    aOute[89] = -0.599740133241022599364100642561;
    aOute[90] = 26.511086896219396180640615057200;
    aOute[91] = -0.891165881964516515267860086169;
    aOute[92] = -9.555669467523230764527397695929;
    aOute[93] = -27.848985271225764392966084415093;
    aOute[94] = 0.995154163717485440443510924524;
    aOute[95] = 12.000706520289247336563676071819;
    aOute[96] = -1.267339079273520852098044997547;
    aOute[97] = -31.161101243121080983655701857060;
    aOute[98] = -34.628509329969901386903075035661;
    aOute[99] = -5.150833542821337651673729851609;
    aOute[100] = 10.593806256153664335784014838282;
    aOute[101] = -8.196300357999364294414590403903;
    aOute[102] = -22.637151099664265530009288340807;
    aOute[103] = -10.901623765334365501189495262224;
    aOute[104] = -12.658309242106470904332127247471;
    aOute[105] = 17.778111973846115745345741743222;
    aOute[106] = 5.153351217652725857476525561651;
    aOute[107] = -9.358680062607284355635783867911;
    aOute[108] = -14.399475342070774175340375222731;
    aOute[109] = 6.611984334634814608477881847648;
    aOute[110] = 18.558941075256502273305159178562;
    aOute[111] = 11.441473946441032438769980217330;
    aOute[112] = -9.667385539984092446275099064223;
    aOute[113] = -16.476745576231461853922155569308;
    aOute[114] = 4.285033884930492753539965633536;
    aOute[115] = 16.271168846193518220388796180487;
    aOute[116] = 2.918939127424355017836887782323;
    aOute[117] = -7.767867844336468152732777525671;
    aOute[118] = -9.482328162687183592538531229366;
    aOute[119] = -0.615621151699472290985681866005;
    aOute[120] = 21.681560608423470881689354428090;
    aOute[121] = 6.532353241739095928153346903855;
    aOute[122] = -5.524953323024761608905919274548;
    aOute[123] = -7.892259655151920405558030324755;
    aOute[124] = 6.362210000323178071823804202722;
    aOute[125] = 16.398944062458507886503866757266;
    aOute[126] = -12.092496000334874040049726318102;
    aOute[127] = -3.711999756760710145186976660625;
    aOute[128] = -8.692154712995293408539509982802;
    aOute[129] = -19.734654268804057153374742483720;
    aOute[130] = 14.205383384396290225026859843638;
    aOute[131] = 4.061257614704441465391937526874;
    aOute[132] = -3.711999756760710145186976660625;
    aOute[133] = -14.977006593683620749857254850212;
    aOute[134] = -19.953727544665817106306349160150;
    aOute[135] = 14.687602396733844045684236334637;
    aOute[136] = 7.192944050667620992101092269877;
    aOute[137] = -15.634889606683261220609892916400;
    aOute[138] = -8.846712117729042645919435017277;
    aOute[139] = -6.447907663093634944573295797454;
    aOute[140] = 12.607147859606175899216395919211;
    aOute[141] = -8.499968646715450049100581964012;
    aOute[142] = -12.296743110317438052447869267780;
    aOute[143] = -13.488858203653675715827375825029;
    aOute[144] = -18.027733041686584414264871156774;
    aOute[145] = 10.680955975133072044513937726151;
    aOute[146] = -10.812229545475188885461648169439;
    aOute[147] = -15.951518189833178595904428220820;
    aOute[148] = -26.016825455262679867018960067071;
    aOute[149] = -9.791515972543889390067306521814;
    aOute[150] = -5.354550059388626159773139079334;
    aOute[151] = 3.340696394628563492545936242095;
    aOute[152] = -6.280908608969569151270206930349;
    aOute[153] = -7.405953772716732608216716471361;
    aOute[154] = -7.135499673684098098647154984064;
    aOute[155] = -5.727993799210259417975521500921;
    aOute[156] = 3.594008186149934669373351425747;
    aOute[157] = -10.812048447632799152984262036625;
    aOute[158] = -6.659138688315616505519756174181;
    aOute[159] = -8.600317924069601716041688632686;
    aOute[160] = 9.732669513492721691250153526198;
    aOute[161] = -2.647298166201000668706910801120;
    aOute[162] = -9.518631661704658242229015741032;
    aOute[163] = -31.940890393167869376611633924767;
    aOute[164] = -19.200929315138488107095326995477;
    aOute[165] = 12.500065356879911604437438654713;
    aOute[166] = -2.636290785390785096353738481412;
    aOute[167] = -5.054254645014958136073346395278;
    aOute[168] = -16.321992192999857707036426290870;
    aOute[169] = -8.742036100432812872895738109946;
    aOute[170] = 19.386734436211835230778888217174;
    aOute[171] = -2.773229509339630638464768708218;
    aOute[172] = -6.821679994710295069637595588574;
    aOute[173] = -7.546101381791159745660024782410;
    aOute[174] = -8.662669114146970272827275039162;
    aOute[175] = 0.163478843261819883903029904104;
    aOute[176] = 1.140455832716614104072050395189;
    aOute[177] = -5.238757713039999153181724977912;
    aOute[178] = -24.354507237356877169531799154356;
    aOute[179] = -33.998117511717708794094505719841;
    aOute[180] = 13.062347736236946005305981088895;
    aOute[181] = -0.690674422385158970350005347427;
    aOute[182] = -6.191623537465112292466073995456;
    aOute[183] = -33.116741564762826044443499995396;
    aOute[184] = -34.685359055258913940633647143841;
    aOute[185] = 17.863289305262792794337656232528;
    aOute[186] = 2.230941629134592041339146817336;
    aOute[187] = -6.187222171614940258166370767867;
    aOute[188] = -28.811973511460436725428735371679;
    aOute[189] = -34.760083116886363541198079474270;
    aOute[190] = 11.692909660262126081420319678728;
    aOute[191] = -0.161371175476341277565950349526;
    aOute[192] = -3.250773276653635690536248148419;
    aOute[193] = -22.999034410598561350980162387714;
    aOute[194] = -34.760083116886363541198079474270;
    aOute[195] = 13.903916990598844449777971021831;
    aOute[196] = -0.223926560486721326626025074802;
    aOute[197] = -3.630186397835856837446044664830;
    aOute[198] = -3.689556613508623517105888822698;
    aOute[199] = -23.620267059508954332613939186558;
    aOute[200] = 11.829222576957672430353341042064;
    aOute[201] = 9.671211493756414867561943538021;
    aOute[202] = -0.614006959879394687007447828364;
    aOute[203] = -14.857569257624204439593995630275;
    aOute[204] = -33.460351617957947212289582239464;
    aOute[205] = 12.048916974901837662059733702336;
    aOute[206] = -0.280453683130871989082066875199;
    aOute[207] = -9.432179086886932850575249176472;
    aOute[208] = -13.415127217793331837469850142952;
    aOute[209] = -33.721507825971947625021130079404;
    aOute[210] = 11.252993650728063457222560828086;
    aOute[211] = 5.301652360524331975000222882954;
    aOute[212] = -3.509328009665552805529387114802;
    aOute[213] = -14.136290952591386016479191312101;
    aOute[214] = -19.772046698998934743940480984747;
    aOute[215] = 16.564656283556463023387550492771;
    aOute[216] = -0.585428901412651692126587477105;
    aOute[217] = -6.340443474480793817349422170082;
    aOute[218] = -11.841309968184146939051970548462;
    aOute[219] = -8.145957291403112421335208637174;
    aOute[220] = 16.575451527689910591334410128184;
    aOute[221] = 4.773346151337615417276083462639;
    aOute[222] = -15.243442106429437643555502290837;
    aOute[223] = -16.838150018546226505122831440531;
    aOute[224] = -19.818154495366794520805342472158;
    aOute[225] = 9.776727629035866229401108284947;
    aOute[226] = 2.449078743718190231959397351602;
    aOute[227] = -16.794156970876549905824504094198;
    aOute[228] = -8.329641295075409956893963681068;
    aOute[229] = -3.149405450182843235751306565362;
    aOute[230] = 11.414923882509038222110575588886;
    aOute[231] = -0.803499110010593264874501073791;
    aOute[232] = -5.194010843425130197203998250188;
    aOute[233] = -5.512066239414147084119122155244;
    aOute[234] = 3.251298259389269862396076860023;
    aOute[235] = 6.888610688118261826673460745951;
    aOute[236] = 3.777706637628776231707661281689;
    aOute[237] = -0.752706456005254542773741377459;
    aOute[238] = -14.372064503065812957061098131817;
    aOute[239] = 3.435514774472468335630992442020;
    aOute[240] = 17.615241912514676414502901025116;
    aOute[241] = -3.512753068241333043175700368010;
    aOute[242] = -1.854224711266238667661809813580;
    aOute[243] = -10.984353027212023334868717938662;
    aOute[244] = 0.605003073110891631536389922985;
    aOute[245] = 27.168413130356572793289160472341;
    aOute[246] = 26.795161318753208945508959004655;
    aOute[247] = 13.220386618005866452563168422785;
    aOute[248] = -1.057988489557932787477056990610;
    aOute[249] = 13.147036603541140564743727736641;
    aOute[250] = -6.971926812953943120021449431079;
    aOute[251] = -17.666960989521196978557782131247;
    aOute[252] = -19.803143202257391664034003042616;
    aOute[253] = -19.477101563780816917414995259605;
    aOute[254] = -29.931309908884102810588956344873;
    aOute[255] = 11.441631579215309599817373964470;
    aOute[256] = -11.190493802307187820588296744972;
    aOute[257] = -14.636430792258929400873057602439;
    aOute[258] = -14.217280826152753903102166077588;
    aOute[259] = -34.371195292965701639786857413128;
    aOute[260] = 25.718292722143040407445369055495;
    aOute[261] = 18.210666950511249950750425341539;
    aOute[262] = 5.124241089997775944198110664729;
    aOute[263] = -11.020892612687932299309068184812;
    aOute[264] = -20.356667074368353098634543130174;
    aOute[265] = 18.579992880900373819486048887484;
    aOute[266] = 11.276865685231157954149239230901;
    aOute[267] = 0.114956453979764838813082405977;
    aOute[268] = -9.180853108340896184813573199790;
    aOute[269] = -14.242814463566880078815302113071;
    aOute[270] = 13.255539011671197968667001987342;
    aOute[271] = -0.810268602613907984455465793872;
    aOute[272] = -15.607801497160597392621639301069;
    aOute[273] = -19.744000910167049056553878472187;
    aOute[274] = -8.533556939448464717656861466821;
    aOute[275] = -6.648558464879472218456157861510;
    aOute[276] = 2.478972662028056106464646290988;
    aOute[277] = -3.872630733782064016423873908934;
    aOute[278] = -7.514785129869033575289449800039;
    aOute[279] = -20.381108910756900343130837427452;
    aOute[280] = 11.714110434877293087652105896268;
    aOute[281] = 4.402884239709353941805147769628;
    aOute[282] = 7.744915310340108938191860943334;
    aOute[283] = -12.293482675754534128031991713215;
    aOute[284] = -13.390768019026074497901390714105;
    aOute[285] = 12.599365327307127770950501144398;
    aOute[286] = 1.164817517508382138657907489687;
    aOute[287] = -6.870743207692028953204044228187;
    aOute[288] = -9.941405177129459147522538842168;
    aOute[289] = 1.090432655703475361974597035442;
    aOute[290] = 14.447812760924684027941111708060;
    aOute[291] = 10.368230678630878216495148080867;
    aOute[292] = 11.504813157898803765988304803614;
    aOute[293] = -10.008176132009490189034295326564;
    aOute[294] = -15.939152691147155493922582536470;
    aOute[295] = 17.692781628601824905899775330909;
    aOute[296] = 13.705622033209399646125348226633;
    aOute[297] = 6.437440551014769418713967752410;
    aOute[298] = -4.927782908846036313832428277237;
    aOute[299] = -17.310832892573625230170364375226;
    aOute[300] = 11.627130511040213889373262645677;
    aOute[301] = -4.267643793813418184868169191759;
    aOute[302] = 2.162472263507926761860744591104;
    aOute[303] = -16.427307884624866574085899628699;
    aOute[304] = -23.476610645583061653951517655514;
    aOute[305] = 16.638033495690127239186040242203;
    aOute[306] = 10.970748657810164061743307684083;
    aOute[307] = -6.513793953116242541057090420509;
    aOute[308] = -3.107305693837132842816117772600;
    aOute[309] = -14.140328388027127104464852891397;
    aOute[310] = 7.928501312682240786955389921786;
    aOute[311] = 22.984736372060130094041596748866;
    aOute[312] = -0.499108019612376885287829964000;
    aOute[313] = -4.008365034136486038107705098810;
    aOute[314] = -26.361220190771557980724537628703;
    aOute[315] = 9.149611382360744116226669575553;
    aOute[316] = 4.949467905593613714643197454279;
    aOute[317] = -0.838656065481923063131830531347;
    aOute[318] = -0.480784520757674360691424908509;
    aOute[319] = -19.974439850128245410587624064647;
    aOute[320] = 14.982061155410915986863074067514;
    aOute[321] = 5.579683008679146460906395077473;
    aOute[322] = 4.151351616800103450088954559760;
    aOute[323] = -14.640581338083157092455621750560;
    aOute[324] = -16.316648528336010315342718968168;
    aOute[325] = 11.295120909155903277110155613627;
    aOute[326] = 7.949738532931195145181391126243;
    aOute[327] = -8.907136100148015600552753312513;
    aOute[328] = -12.714337855571923086017704918049;
    aOute[329] = -17.875748106283424476714571937919;
    aOute[330] = 11.167180123525954016372452315409;
    aOute[331] = 11.100629994122675370249453408178;
    aOute[332] = 3.226653560342472193411822445341;
    aOute[333] = -7.236304183924755939472106547328;
    aOute[334] = -11.274371278909555371683381963521;
    aOute[335] = 18.281430106034296301231734105386;
    aOute[336] = 5.407313859237697251103327289457;
    aOute[337] = -0.558899517830174330512704727880;
    aOute[338] = 0.620836223238954376668630175118;
    aOute[339] = -18.025212444811792522614268818870;
    aOute[340] = 17.192621157863424485867653856985;
    aOute[341] = 5.851101381139955925902995659271;
    aOute[342] = 1.222732110587017473335436079651;
    aOute[343] = -9.118300352260364505241341248620;
    aOute[344] = -22.611307419716617062022123718634;
    aOute[345] = 17.448661212614393178910177084617;
    aOute[346] = 23.605447980130310980939611908980;
    aOute[347] = 8.993356286474700311828200938180;
    aOute[348] = -6.246983162694403013404098601313;
    aOute[349] = -32.810083116886268328471487620845;
    aOute[350] = 19.931635821815330444906066986732;
    aOute[351] = 14.728679607851235644488951948006;
    aOute[352] = -6.191653394483475736365107877646;
    aOute[353] = -8.180937633772051853497941920068;
    aOute[354] = -16.333246166981982838706244365312;
    aOute[355] = 15.443482386921361992904166982044;
    aOute[356] = 10.011956263483897444643844210077;
    aOute[357] = -2.628568926049349219198347782367;
    aOute[358] = 2.055288988690687101978937789681;
    aOute[359] = -13.405547342315381698085730022285;
    aOute[360] = 25.036375979950616255109707708471;
    aOute[361] = 0.393053360005446161817843631070;
    aOute[362] = 1.113952441460661191285907989368;
    aOute[363] = -11.074548201555769466608580842149;
    aOute[364] = 4.627842264899031832214859605301;
    aOute[365] = 23.398395717337344024144840659574;
    aOute[366] = 8.297519436689649197091966925655;
    aOute[367] = 2.288869452773643953236160086817;
    aOute[368] = -12.832046794299692749063979135826;
    aOute[369] = 4.473974757925892831167402619030;
    aOute[370] = 34.780083116886359562158759217709;
    aOute[371] = 34.780083116886359562158759217709;
    aOute[372] = 22.756826047907203758313698926941;
    aOute[373] = -4.214988120104019131417771859560;
    aOute[374] = 1.909178333677674199719831449329;
    aOute[375] = -3.608035557005876192704363347730;
    aOute[376] = 4.814161468933737531017413857626;
    aOute[377] = 0.980938422988455016593434265815;
    aOute[378] = -1.972023612026655658979734653258;
    aOute[379] = -10.600023899188487419564808078576;
    aOute[380] = -4.226054029003589285196085256757;
    aOute[381] = 6.175104219858911669405188149540;
    aOute[382] = 2.314161468933745524623191158753;
    aOute[383] = -2.116447977708054040135721152183;
    aOute[384] = -11.795283986085307148528045217972;
    aOute[385] = -6.325736183533089196373566664988;
    aOute[386] = 4.456053134425366657467293407535;
    aOute[387] = 8.022464810495897324926772853360;
    aOute[388] = -9.459328990648200985447147104423;
    aOute[389] = -22.791591062995369298960213200189;
    aOute[390] = -3.889895775499108587780483503593;
    aOute[391] = 6.472672358612224563501058582915;
    aOute[392] = 5.590456556798844012234894762514;
    aOute[393] = -1.472432750021260616790641506668;
    aOute[394] = -11.119128448706133482914992782753;
    aOute[395] = -3.767927168717190689051221852424;
    aOute[396] = 8.473502756544565528429302503355;
    aOute[397] = 2.898389517342121024512380245142;
    aOute[398] = 1.058017493764253469379355010460;
    aOute[399] = -10.979621141501064940371179545764;
    aOute[400] = -8.074971876830062456065206788480;
    aOute[401] = 4.418745480391357283167508285260;
    aOute[402] = 7.814161468933737531017413857626;
    aOute[403] = -7.577511030435574390651254361728;
    aOute[404] = -10.208899214554813639210806286428;
    aOute[405] = -5.238426873449856024933524167864;
    aOute[406] = 3.973363019082453817532041284721;
    aOute[407] = 5.084496278291243065439175552456;
    aOute[408] = 0.455095866559146966512372500802;
    aOute[409] = -12.012646945433258949265109549742;
    aOute[410] = -1.812590908830089997394452439039;
    aOute[411] = 13.979531970275294483485595264938;
    aOute[412] = 8.975654751373280504367357934825;
    aOute[413] = -2.093971777623848673499651340535;
    aOute[414] = -16.665253504475796120232189423405;
    aOute[415] = -5.791777053159317567576636065496;
    aOute[416] = 8.314161468933738419195833557751;
    aOute[417] = 5.224630438687330347136139607755;
    aOute[418] = 0.117232366049706923649686984845;
    aOute[419] = -9.956260961529933339875242381822;
    aOute[420] = -0.113153831460361439797956961684;
    aOute[421] = 9.433399329045373704616395116318;
    aOute[422] = 6.814161468933737531017413857626;
    aOute[423] = -6.174464856977839311014122358756;
    aOute[424] = -8.317490572236904355918341025244;
    aOute[425] = 4.101034021699699216867429640843;
    aOute[426] = 16.246668685812682753066837904043;
    aOute[427] = 12.873499626727284095295544830151;
    aOute[428] = 1.607512942379803799752835402614;
    aOute[429] = -19.159351780450588620396956684999;
    aOute[430] = 6.813289565216547494230781012448;
    aOute[431] = 11.289421570354393864477060560603;
    aOute[432] = 3.484658303101352494479669985594;
    aOute[433] = 1.047403347114264526851457048906;
    aOute[434] = -8.438706429572865985733187699225;
    aOute[435] = 25.573794841826217094649109640159;
    aOute[436] = 12.006039912710912531679241510574;
    aOute[437] = 2.713188911473991993261734023690;
    aOute[438] = 0.344623516853865274889301417716;
    aOute[439] = -5.244071056488379234394869854441;
    aOute[440] = -1.476389508129834826632986732875;
    aOute[441] = 3.308061190690648523116124124499;
    aOute[442] = 4.722452731319182639424525405047;
    aOute[443] = 0.692577772131530489119199955894;
    aOute[444] = -7.805319632012657216080242505996;
    aOute[445] = -3.211732735638463065441783328424;
    aOute[446] = 7.828751857639585587378405762138;
    aOute[447] = 13.985904185912749753128991869744;
    aOute[448] = -1.643103652883143084295625158120;
    aOute[449] = -9.888941924416547735177118738648;
    aOute[450] = 27.650056160453015507982854614966;
    aOute[451] = 15.989833805457939419625290611293;
    aOute[452] = 2.558579945993551874039440008346;
    aOute[453] = -5.409160506265659584812510729535;
    aOute[454] = -14.063230080941801602989471575711;
    aOute[455] = 34.707990016511217845618375577033;
    aOute[456] = 2.973412380201165383653005847009;
    aOute[457] = 0.075925297304148722798977644288;
    aOute[458] = -2.821192157624883734001741686370;
    aOute[459] = -5.271903912290406779561635630671;
    aOute[460] = 31.864204038489077674967120401561;
    aOute[461] = 20.456053799188374853201821679249;
    aOute[462] = 5.988718231403898251130613061832;
    aOute[463] = -2.866263024925616686999774174183;
    aOute[464] = -7.244343026119493877956756477943;
    aOute[465] = 31.423323006109900035198734258302;
    aOute[466] = 21.259801781647592378021727199666;
    aOute[467] = 7.486567940204695403849655122031;
    aOute[468] = 1.149710634296620792582643844071;
    aOute[469] = -14.009146664881834354332568182144;
    aOute[470] = 33.518570846886468928005342604592;
    aOute[471] = 34.780083116886359562158759217709;
    aOute[472] = 1.263988456900567936003199065453;
    aOute[473] = -3.388234026458271763715401903028;
    aOute[474] = -22.063658350576208277971090865321;
    aOute[475] = 33.684216588443426587673457106575;
    aOute[476] = 26.578285653354136286452558124438;
    aOute[477] = -6.902950791282751019650731905131;
    aOute[478] = 5.267339003364505778392867796356;
    aOute[479] = 4.679465358175151301622918253997;
    aOute[480] = 34.780083116886359562158759217709;
    aOute[481] = 18.775258925110755114928906550631;
    aOute[482] = -4.355650123965109976609255681979;
    aOute[483] = 1.313388717152813356392471177969;
    aOute[484] = 2.386996918733698791470487776678;
    aOute[485] = 34.620644676230696745733439456671;
    aOute[486] = 21.839293197962323489491609507240;
    aOute[487] = 5.032800861316990115312819398241;
    aOute[488] = -2.666566767462771547769762037206;
    aOute[489] = -1.669849120407292275558575056493;
    aOute[490] = 34.713798955774393562023760750890;
    aOute[491] = 20.460043828076777572277933359146;
    aOute[492] = 9.165239137188647688958553771954;
    aOute[493] = 2.514771846964982238148422766244;
    aOute[494] = 14.637194785518655493206097162329;
    aOute[495] = 34.780083116886359562158759217709;
    aOute[496] = 34.780083116886359562158759217709;
    aOute[497] = 11.894021177220951202002652280498;
    aOute[498] = -4.401372780115357485897220612969;
    aOute[499] = -4.737203054874510321781144739361;
    aOute[500] = 7.388216348648568576606976421317;
    aOute[501] = 23.244371491572053400886943563819;
    aOute[502] = 7.314161468933737531017413857626;
    aOute[503] = 2.186542330918926246852151962230;
    aOute[504] = -5.276664325194610150049356889213;
    aOute[505] = 7.412005015235014759866771782981;
    aOute[506] = 19.954459422203960627939522964880;
    aOute[507] = 7.814161468933737531017413857626;
    aOute[508] = 2.314161468933745524623191158753;
    aOute[509] = -2.941734454447907598506617432577;
    aOute[510] = 3.630788218274721756984035891946;
    aOute[511] = 16.917659984613305113043679739349;
    aOute[512] = 10.180445547083619928230291407090;
    aOute[513] = -0.905211128161025113669779784686;
    aOute[514] = -9.398980328548423912593534623738;
    aOute[515] = 6.954536531347243588641049427679;
    aOute[516] = 23.443636703196133908022602554411;
    aOute[517] = 6.814161468933737531017413857626;
    aOute[518] = 3.026650478529138155181499314494;
    aOute[519] = -8.668668677851615100848903239239;
    aOute[520] = 6.971483746735942510497352486709;
    aOute[521] = 26.968111974240013495318635250442;
    aOute[522] = 6.783616512671087228625310672214;
    aOute[523] = 5.814161468933737531017413857626;
    aOute[524] = -4.096398665603918409772177255945;
    aOute[525] = 8.092950440416014856737092486583;
    aOute[526] = 19.816714931213439854218449909240;
    aOute[527] = 7.807523468180282755213283962803;
    aOute[528] = 7.814161468933737531017413857626;
    aOute[529] = 4.710615458713315106820118671749;
    aOute[530] = 11.502577010982461658272768545430;
    aOute[531] = 14.439840945782995262902659305837;
    aOute[532] = 8.828694164850249492815237317700;
    aOute[533] = -1.751479180374845068257627644925;
    aOute[534] = -7.663080461650141117502244014759;
    aOute[535] = -6.527771597965650762773748283507;
    aOute[536] = 16.745253293421519913408701540902;
    aOute[537] = 17.960496070887995045950447092764;
    aOute[538] = -4.322996293722198402065259870142;
    aOute[539] = -20.283400585181354358610406052321;
    aOute[540] = 11.978164172276118293325453123543;
    aOute[541] = 25.427603912302782163123993086629;
    aOute[542] = 8.266329171131946651485122856684;
    aOute[543] = -0.673726232192641827367651785607;
    aOute[544] = -7.646690234084824844273953203810;
    aOute[545] = 10.882974398618966205276592518203;
    aOute[546] = 26.968111974240013495318635250442;
    aOute[547] = 6.975149952258738572652418952202;
    aOute[548] = 6.168918755702086365033665060764;
    aOute[549] = 6.566657242031994634601232974092;
    aOute[550] = -1.126804259896590387768355867593;
    aOute[551] = 20.881026823462480734860946540721;
    aOute[552] = 16.279029495808035932213897467591;
    aOute[553] = 8.172126219385333456557418685406;
    aOute[554] = -5.774282761019464693674763111630;
    aOute[555] = -0.637609761530236518467518180842;
    aOute[556] = 15.587282915542719052837128401734;
    aOute[557] = 8.196702508927657504500530194491;
    aOute[558] = 3.943211433577485713897203822853;
    aOute[559] = -13.829119552130356041175218706485;
    aOute[560] = -3.394740394659540694277666261769;
    aOute[561] = 12.437803003312735583563153340947;
    aOute[562] = 8.929863270691203780415889923461;
    aOute[563] = -5.306200373189214225533305580029;
    aOute[564] = -8.585218490766322929630405269563;
    aOute[565] = 13.884412592069439895681171037722;
    aOute[566] = 17.338258368833322720092837698758;
    aOute[567] = 12.849639547406763284698172356002;
    aOute[568] = 0.555095014296749456050861226686;
    aOute[569] = -12.935281589071230712306714849547;
    aOute[570] = 13.903571299669930994014066527598;
    aOute[571] = 23.355934306443160153321514371783;
    aOute[572] = 26.585204099718442449784561176784;
    aOute[573] = 15.109178119604555590171912626829;
    aOute[574] = 0.994150837758983052516725820169;
    aOute[575] = 11.945145270735416431762132560834;
    aOute[576] = -9.187910219642652975835517281666;
    aOute[577] = 8.194270359959766381052759243175;
    aOute[578] = -4.181836071509643204535677796230;
    aOute[579] = -17.049722280751165470746855135076;
    aOute[580] = 13.210079495944851757371907297056;
    aOute[581] = -0.476139262145280461879792710533;
    aOute[582] = 20.131290282248539114107188652270;
    aOute[583] = 4.183239231128037083351500768913;
    aOute[584] = -8.367961740005661397390213096514;
    aOute[585] = 15.844642151916056960203604830895;
    aOute[586] = 2.208307894811383231825629991363;
    aOute[587] = 12.944337300683955049862561281770;
    aOute[588] = 10.879175944744439163969218498096;
    aOute[589] = -4.719255243634544960684706893517;
    aOute[590] = 14.041571795631444530272347037680;
    aOute[591] = 7.501511085116463029009992169449;
    aOute[592] = 13.898717791654414455138066841755;
    aOute[593] = 6.660144683427950518250781897223;
    aOute[594] = -10.038025599490644879097089869902;
    aOute[595] = 16.422173037595005240518730715849;
    aOute[596] = 4.447748279538211058081742521608;
    aOute[597] = 18.045495354237768736993530183099;
    aOute[598] = 1.672326746139092046306018346513;
    aOute[599] = -17.300788909329021691974048735574;
    aOute[600] = 8.785632133898539919414361065719;
    aOute[601] = 4.008576225299805173563072457910;
    aOute[602] = 8.870858998661523386886074149515;
    aOute[603] = 6.404193837924704091335570410592;
    aOute[604] = 5.201088308932622972236003988655;
    aOute[605] = 13.949757135511875461020281363744;
    aOute[606] = 8.674559365092852303291692805942;
    aOute[607] = 12.396258228621482189169000776019;
    aOute[608] = 9.161728486709348473482350527775;
    aOute[609] = 8.384205238997635589726087346207;
    aOute[610] = 21.862687762648608469362443429418;
    aOute[611] = 4.925017981937862110441983531928;
    aOute[612] = -1.376527060368659238775990161230;
    aOute[613] = 14.668422397070193241574997955468;
    aOute[614] = 3.339167437472250643537563519203;
    aOute[615] = 20.721361409938509723360766656697;
    aOute[616] = 12.047526627436932145087666867767;
    aOute[617] = 0.799429382126391918639285449899;
    aOute[618] = 16.494099652724152349492214852944;
    aOute[619] = 5.451525861036491171773832320469;
    aOute[620] = 31.374735720583732501154372585006;
    aOute[621] = 34.780083116886359562158759217709;
    aOute[622] = 10.889026523553322789439334883355;
    aOute[623] = 6.486860580427967981620440696133;
    aOute[624] = 2.031195935780393924119380244520;
    cOute[0] = 1.138993039872924661892739095492;
    cOute[1] = 5.666000551878929591964606515830;
    cOute[2] = 5.474641747201421715374181076186;
    cOute[3] = 7.254793953637979520010503620142;
    cOute[4] = 10.697956509877563746613304829225;
    cOute[5] = -13.861117206305278770628319762181;
    cOute[6] = 11.992025621461950279922348272521;
    cOute[7] = 5.887208343929106746372781344689;
    cOute[8] = 12.207059599034081998070178087801;
    cOute[9] = 10.777138536928198675468593137339;
    cOute[10] = 4.176941298645776790010586410062;
    cOute[11] = 15.583733715842335243451088899747;
    cOute[12] = 36.693160039963295560028200270608;
    cOute[13] = 23.180117154146117286472872365266;
    cOute[14] = 20.068071817002966383824968943372;
    cOute[15] = 11.362248710939892859528299595695;
    cOute[16] = 32.360670975530574366985092638060;
    cOute[17] = 37.494442091245325343606964452192;
    cOute[18] = 37.654698501501734142493660328910;
    cOute[19] = 12.600597978279401090162537002470;
    cOute[20] = 15.918538552134576136154464620631;
    cOute[21] = 36.672218286147177934708452085033;
    cOute[22] = 38.295724142527340916331013431773;
    cOute[23] = 38.455980552783763926072424510494;
    cOute[24] = 15.915198135428463288576494960580;
    cOute[25] = 7.913177126138285721879128686851;
    cOute[26] = 12.720289130072464089948880427983;
    cOute[27] = 13.648385283616098817560668976512;
    cOute[28] = 13.190113868061965618494468799327;
    cOute[29] = 11.935337334600271574913676886354;
    cOute[30] = 14.925406102892271320570216630585;
    cOute[31] = 17.106430078276169126638706075028;
    cOute[32] = 14.650061671682260922011664661113;
    cOute[33] = 13.699849431782052278094852226786;
    cOute[34] = 18.717315162041529674752382561564;
    cOute[35] = 5.372093388794739077241047198186;
    cOute[36] = 11.741066987707247548655686841812;
    cOute[37] = 25.239186812223216804795811185613;
    cOute[38] = 25.350295747923965450354444328696;
    cOute[39] = 21.029580883172037886197358602658;
    cOute[40] = 3.808303556646996490542278479552;
    cOute[41] = 12.406416749180191771984027582221;
    cOute[42] = 18.319535319117335347982589155436;
    cOute[43] = 40.966429326561645041238080011681;
    cOute[44] = 17.255655246081786913237010594457;
    cOute[45] = 9.035921138567802657348693173844;
    cOute[46] = 15.455113176971547517268845695071;
    cOute[47] = 29.607315302695663916665580472909;
    cOute[48] = 42.375006250079607639236201066524;
    cOute[49] = 17.151450824819313822899857768789;
    cOute[50] = 27.220777828141830667618705774657;
    cOute[51] = 11.955430202072841794347368704621;
    cOute[52] = 37.522630078350644566853588912636;
    cOute[53] = 35.841100639468052690972399432212;
    cOute[54] = 30.359399044106762488581807701848;
    cOute[55] = 15.288683572414564437735862156842;
    cOute[56] = 14.044087258813513940935990831349;
    cOute[57] = 33.664703256896018501720391213894;
    cOute[58] = 43.261554887216817633088794536889;
    cOute[59] = 33.673482224844036636568489484489;
    cOute[60] = -7.848457490407599657089576794533;
    cOute[61] = 22.688695213028818642442274722271;
    cOute[62] = 34.601345713855700125805014977232;
    cOute[63] = 44.866236963040130092394974781200;
    cOute[64] = 31.880780625984485254775790963322;
    cOute[65] = 3.474541210123223144989879074274;
    cOute[66] = 20.704487457549028306402760790661;
    cOute[67] = 21.706025497571062032875488512218;
    cOute[68] = 25.372231890587645608547973097302;
    cOute[69] = 21.038840982844064342316414695233;
    cOute[70] = 7.219198918545458809603587724268;
    cOute[71] = 18.665667840138461031074257334694;
    cOute[72] = 20.933369357313228675820937496610;
    cOute[73] = 29.160256776271882728224227321334;
    cOute[74] = 20.107283455799915827810764312744;
    cOute[75] = 6.706893644495744766231837274972;
    cOute[76] = 30.672188725222660821145836962387;
    cOute[77] = 20.478470152147384197860446874984;
    cOute[78] = 47.173232988736764070836215978488;
    cOute[79] = 43.347742445441944880712981102988;
    cOute[80] = -0.738724760905446697911713727080;
    cOute[81] = 25.184984782725809537851091590710;
    cOute[82] = 18.230739880496802385323462658562;
    cOute[83] = 47.972442471031200739162159152329;
    cOute[84] = 46.219257096199541479109029751271;
    cOute[85] = -13.379177521362390024250998976640;
    cOute[86] = 6.191732177581128482302119664382;
    cOute[87] = 24.210294457162039094555439078249;
    cOute[88] = 36.771396606745540225347212981433;
    cOute[89] = 13.527141593792199714130219945218;
    cOute[90] = -14.743254242808600196212864830159;
    cOute[91] = 7.319896153877347266814012982650;
    cOute[92] = 24.022849749562769972044407040812;
    cOute[93] = 42.741487312887599614441569428891;
    cOute[94] = 12.435213128785175484836145187728;
    cOute[95] = 1.820880854621270028559365528054;
    cOute[96] = 19.769533397696637422313870047219;
    cOute[97] = 49.698188106758109938709822017699;
    cOute[98] = 50.426293121931898610910138813779;
    cOute[99] = 21.969098818681906948313553584740;
    cOute[100] = 4.268126638284175378146301227389;
    cOute[101] = 26.137891650262879750243882881477;
    cOute[102] = 42.007354527277378508642868837342;
    cOute[103] = 27.735783993438381855867191916332;
    cOute[104] = 29.863216382100539192379073938355;
    cOute[105] = -2.368161818055459111320715237525;
    cOute[106] = 10.402340945002501726435184536967;
    cOute[107] = 28.846220846659107905907148960978;
    cOute[108] = 31.136549513548246181926515419036;
    cOute[109] = 10.746007663532616049906209809706;
    cOute[110] = -2.413394992632387658915149586392;
    cOute[111] = 2.948832549402431890683828896726;
    cOute[112] = 27.807310372761168792976604891010;
    cOute[113] = 34.767883032010253430144075537100;
    cOute[114] = 15.776817665231064324871113058180;
    cOute[115] = 5.459873789632221807721634831978;
    cOute[116] = 13.516573816799963836388087656815;
    cOute[117] = 25.221599302951080545653894660063;
    cOute[118] = 28.351775578293427315657027065754;
    cOute[119] = 19.489955242976918725616997107863;
    cOute[120] = -3.412661212427466672636455768952;
    cOute[121] = 13.058759905007510937480219581630;
    cOute[122] = 25.513104768741598604719911236316;
    cOute[123] = 27.489759987285538045398425310850;
    cOute[124] = 13.578799805700404235153655463364;
    cOute[125] = 3.502273777755944283995859223069;
    cOute[126] = 30.182022248539947639756064745598;
    cOute[127] = 23.907094300158302502268270472996;
    cOute[128] = 28.475130677659734601547825150192;
    cOute[129] = 42.132971889899195616635552141815;
    cOute[130] = 7.128854125236977701263185736025;
    cOute[131] = 16.905130918771828163471582229249;
    cOute[132] = 24.453126299365170126520752091892;
    cOute[133] = 36.301109157786250136723538162187;
    cOute[134] = 41.661373257153194060720124980435;
    cOute[135] = 7.549573328070484201646195288049;
    cOute[136] = 18.999896147856038197687666979618;
    cOute[137] = 36.761544065163391792339098174125;
    cOute[138] = 30.509590924162246494688588427380;
    cOute[139] = 28.195840728776566663782432442531;
    cOute[140] = 9.692984728037755459695290483069;
    cOute[141] = 32.308838654414451241336792008951;
    cOute[142] = 39.536038617377556647625169716775;
    cOute[143] = 37.737494690096561100745020667091;
    cOute[144] = 40.075420681011621582001680508256;
    cOute[145] = 12.463717434365067404655746940989;
    cOute[146] = 35.303088349995761063837562687695;
    cOute[147] = 44.301641374713923937633808236569;
    cOute[148] = 53.972482045432371933202375657856;
    cOute[149] = 33.334953028617107406716968398541;
    cOute[150] = 29.228844686249896511753831873648;
    cOute[151] = 20.083021425880140498065884457901;
    cOute[152] = 29.359160388893052129333227640018;
    cOute[153] = 31.440357124864600990576946060173;
    cOute[154] = 31.094836389980592628035083180293;
    cOute[155] = 32.111284716374356662527134176344;
    cOute[156] = 20.898938689002406476902251597494;
    cOute[157] = 30.147803655556867141740440274589;
    cOute[158] = 31.354697797482398158308569691144;
    cOute[159] = 32.250741759444281342439353466034;
    cOute[160] = 15.835779730369520024169105454348;
    cOute[161] = 29.235344236780953508514357963577;
    cOute[162] = 34.548953751958926261522719869390;
    cOute[163] = 59.868903475767602628820895915851;
    cOute[164] = 44.828072596691455942163884174079;
    cOute[165] = 12.832026914031525066661743039731;
    cOute[166] = 20.632124822772958339101023739204;
    cOute[167] = 28.861318926968607456728932447731;
    cOute[168] = 42.934658224123452896492381114513;
    cOute[169] = 33.905743810336772980917885433882;
    cOute[170] = 5.165666851941877801834834826877;
    cOute[171] = 16.777010199802973033911257516593;
    cOute[172] = 34.486100809569158798240096075460;
    cOute[173] = 35.778411346248475410902756266296;
    cOute[174] = 36.835076730147434886930568609387;
    cOute[175] = 28.220558929852206375699097407050;
    cOute[176] = 27.313024009568948713422287255526;
    cOute[177] = 41.692488268880765645008068531752;
    cOute[178] = 51.691810808802884480428474489599;
    cOute[179] = 62.082938833067345285599003545940;
    cOute[180] = 13.578248252031817955298720335122;
    cOute[181] = 28.946211840899920275660406332463;
    cOute[182] = 30.422817450834443775420368183404;
    cOute[183] = 63.489454356100864629297575447708;
    cOute[184] = 64.257262604065800815078546293080;
    cOute[185] = 12.183087597085473063884819566738;
    cOute[186] = 30.724202973139885131104165338911;
    cOute[187] = 35.585140616183629447277780855075;
    cOute[188] = 61.248432799312688246118341339752;
    cOute[189] = 65.058544655347844809512025676668;
    cOute[190] = 10.811193826852459665133210364729;
    cOute[191] = 25.946059334761400805291486904025;
    cOute[192] = 35.438281483129721038949355715886;
    cOute[193] = 56.200870720713787420663720695302;
    cOute[194] = 65.859826706629931436509650666267;
    cOute[195] = 10.352001738004556585792670375668;
    cOute[196] = 28.922936396677229708984668832272;
    cOute[197] = 35.088439924478237230687227565795;
    cOute[198] = 33.278481055106418295963521813974;
    cOute[199] = 56.977285112758991658665763679892;
    cOute[200] = 17.813475944715005283569553284906;
    cOute[201] = 32.233254971587264492427493678406;
    cOute[202] = 26.963729292107451840365683892742;
    cOute[203] = 48.862328624356017314767086645588;
    cOute[204] = 66.022110137483323910601029638201;
    cOute[205] = 17.641200430583154457053751684725;
    cOute[206] = 29.927514113090118996751698432490;
    cOute[207] = 41.046244212947257778978382702917;
    cOute[208] = 47.615608149601996501587564125657;
    cOute[209] = 68.263672860475978154681797605008;
    cOute[210] = 30.832198829933460615393414627761;
    cOute[211] = 31.711482609024766787797489087097;
    cOute[212] = 38.433628947480542592529673129320;
    cOute[213] = 39.701662950875480362356029218063;
    cOute[214] = 55.549417913617709530171850929037;
    cOute[215] = 24.062924600779208361700511886738;
    cOute[216] = 30.607404518473224186436709715053;
    cOute[217] = 43.469949486758793000262812711298;
    cOute[218] = 45.643716266317227336912765167654;
    cOute[219] = 42.285932300330458133430511225015;
    cOute[220] = 19.013286571663300605905533302575;
    cOute[221] = 35.924640323540046438210993073881;
    cOute[222] = 53.863716256851397190530406078324;
    cOute[223] = 52.467243492934109383440954843536;
    cOute[224] = 57.065406860079768591731408378109;
    cOute[225] = 28.359296540822249710345204221085;
    cOute[226] = 37.539761429528965663848794065416;
    cOute[227] = 57.044569790644160889314662199467;
    cOute[228] = 44.141683435434146076659089885652;
    cOute[229] = 40.143291186419077121172449551523;
    cOute[230] = 28.086265790615669857288594357669;
    cOute[231] = 29.162886189474978237967661698349;
    cOute[232] = 48.858832978817780201552523067221;
    cOute[233] = 40.582181727159145623318181606010;
    cOute[234] = 34.258701740610696617750363657251;
    cOute[235] = 32.950720732918505007091880543157;
    cOute[236] = 34.458088366348064823796448763460;
    cOute[237] = 40.054322229424755619220377411693;
    cOute[238] = 53.757808933957633712452661711723;
    cOute[239] = 35.375767276809554573446803260595;
    cOute[240] = 21.200877166411274288293498102576;
    cOute[241] = 33.457407905843915330024174181744;
    cOute[242] = 38.874459824365196425333124352619;
    cOute[243] = 50.693625994064881012945988913998;
    cOute[244] = 36.931837007287185770110227167606;
    cOute[245] = 10.941550364286507956990135426167;
    cOute[246] = 9.502862820528108045436965767294;
    cOute[247] = 25.091737114163052524418162647635;
    cOute[248] = 40.417548127122628898177936207503;
    cOute[249] = 24.779962485972429675484818289988;
    cOute[250] = 48.199729203369493291120306821540;
    cOute[251] = 60.260298573438852542949462076649;
    cOute[252] = 63.712627706151572226644930196926;
    cOute[253] = 64.460630760926733273663558065891;
    cOute[254] = 72.135452629890181697192019782960;
    cOute[255] = 30.198402358893126518069038866088;
    cOute[256] = 54.464618206754437323979800567031;
    cOute[257] = 59.011089689473379849005141295493;
    cOute[258] = 61.441309846421702900443051476032;
    cOute[259] = 76.227575173213992343335121404380;
    cOute[260] = 14.478118325645615627195184060838;
    cOute[261] = 22.313863200923265139863360673189;
    cOute[262] = 33.928673007022702279300574446097;
    cOute[263] = 53.148365054283154051972815068439;
    cOute[264] = 62.322129594686984432883036788553;
    cOute[265] = 24.231805133967547760676097823307;
    cOute[266] = 31.943116569711257568542350782081;
    cOute[267] = 43.032118540596272282527934294194;
    cOute[268] = 48.829452895641651366531732492149;
    cOute[269] = 52.976395008803926600648992462084;
    cOute[270] = 29.827655832071478414491139119491;
    cOute[271] = 42.330596448624497440960112726316;
    cOute[272] = 61.526624839630720487093640258536;
    cOute[273] = 64.774627505875926658518437761813;
    cOute[274] = 49.785910665907294969656504690647;
    cOute[275] = 51.091813867407850580093509051949;
    cOute[276] = 41.450456465611253520364698488265;
    cOute[277] = 44.183893634670951655607495922595;
    cOute[278] = 48.453930320151663124761398648843;
    cOute[279] = 68.600522016200386588025139644742;
    cOute[280] = 34.127133287128337713056680513546;
    cOute[281] = 41.482414373646015803842601599172;
    cOute[282] = 28.041642087009620354365324601531;
    cOute[283] = 52.899180723480313304207811597735;
    cOute[284] = 55.718919751950821250829903874546;
    cOute[285] = 33.518524072822302173335629049689;
    cOute[286] = 46.768982038819643776150769554079;
    cOute[287] = 53.532635496701658439633320085704;
    cOute[288] = 47.788559406774105298154609045014;
    cOute[289] = 42.264780018779880776946811238304;
    cOute[290] = 32.419484804309220749019004870206;
    cOute[291] = 36.513974837381020677185006206855;
    cOute[292] = 33.814904331134961523730453336611;
    cOute[293] = 58.349592419268830667533620726317;
    cOute[294] = 65.335320150358697333103918936104;
    cOute[295] = 29.409504703422992832884119707160;
    cOute[296] = 33.000540041677680846987641416490;
    cOute[297] = 35.014644357024124587951519060880;
    cOute[298] = 50.995213783413682051559590036049;
    cOute[299] = 66.411514536498501115602266509086;
    cOute[300] = 36.186820120679378476324927760288;
    cOute[301] = 54.192337209955859123056143289432;
    cOute[302] = 46.030749636169687732945021707565;
    cOute[303] = 67.787901481971104544754780363292;
    cOute[304] = 74.741539181352123932811082340777;
    cOute[305] = 31.367579955427181914728862466291;
    cOute[306] = 38.855694572617622384314017836004;
    cOute[307] = 55.149908373935012662059307331219;
    cOute[308] = 62.786288091326063920405431417748;
    cOute[309] = 69.912579886325900702104263473302;
    cOute[310] = 49.750898685983038660651800455526;
    cOute[311] = 23.127240142696390279297702363692;
    cOute[312] = 55.623049109184776739311928395182;
    cOute[313] = 65.393452610236963096213003154844;
    cOute[314] = 82.304887961360108761255105491728;
    cOute[315] = 49.348154659312015724026423413306;
    cOute[316] = 54.909549561711315845968783833086;
    cOute[317] = 50.082281569930934495005203643814;
    cOute[318] = 67.082832845795437037850206252187;
    cOute[319] = 73.823308748693960978926043026149;
    cOute[320] = 36.465584568543611965196760138497;
    cOute[321] = 49.565096170177646683896455215290;
    cOute[322] = 49.189883459427747425252164248377;
    cOute[323] = 67.511584024381050994634279049933;
    cOute[324] = 67.999761381799814330406661611050;
    cOute[325] = 44.740394648471820460144954267889;
    cOute[326] = 48.510825583764443535983446054161;
    cOute[327] = 63.450431216856408411786105716601;
    cOute[328] = 66.232708361801726937301282305270;
    cOute[329] = 71.339761418921014524130441714078;
    cOute[330] = 44.542219700237190238567563937977;
    cOute[331] = 41.571901767883318257190694566816;
    cOute[332] = 58.193153751658201144891791045666;
    cOute[333] = 64.794299621744102068987558595836;
    cOute[334] = 65.317479048942018948764598462731;
    cOute[335] = 32.676281599268264699276187457144;
    cOute[336] = 34.369606289690374012479878729209;
    cOute[337] = 50.565984454213129595245845848694;
    cOute[338] = 64.532651449548481537021871190518;
    cOute[339] = 76.448230496639837383554549887776;
    cOute[340] = 41.978286340736055137767834821716;
    cOute[341] = 34.647938522711882569637964479625;
    cOute[342] = 56.140625802721146442308963742107;
    cOute[343] = 58.790921817871982568703970173374;
    cOute[344] = 79.051707387880369992672058288008;
    cOute[345] = 38.735426227670515686440921854228;
    cOute[346] = 31.064250353062238474421974387951;
    cOute[347] = 44.626771002568453639014478540048;
    cOute[348] = 60.380323517208019268309726612642;
    cOute[349] = 89.699570296373337896511657163501;
    cOute[350] = 38.724224049315758122702391119674;
    cOute[351] = 36.720289107785376359061046969146;
    cOute[352] = 65.199413172280173967010341584682;
    cOute[353] = 64.963284263546853480875142849982;
    cOute[354] = 73.712657513126060848662746138871;
    cOute[355] = 40.942082227654594817067845724523;
    cOute[356] = 38.858753169330050525331898825243;
    cOute[357] = 70.315417161303898296864645089954;
    cOute[358] = 54.130496612405664791367598809302;
    cOute[359] = 71.363469443373759304449777118862;
    cOute[360] = 28.193584004635965811758069321513;
    cOute[361] = 37.900037376817508061321859713644;
    cOute[362] = 49.093356449827751930570229887962;
    cOute[363] = 63.254835292123111400997004238889;
    cOute[364] = 53.072879440527827910045743919909;
    cOute[365] = 30.235340630032855102626854204573;
    cOute[366] = 45.922358244962460105398349696770;
    cOute[367] = 46.145756721178635473279427969828;
    cOute[368] = 71.727022274122916201122279744595;
    cOute[369] = 52.330906845569046481614350341260;
    cOute[370] = 24.524788677985551288429633132182;
    cOute[371] = 24.685045088241974298171044210903;
    cOute[372] = 35.437176043356188870347978081554;
    cOute[373] = 65.527139787657958436284388881177;
    cOute[374] = 55.667082844100107763551932293922;
    cOute[375] = 64.850871284481002021493623033166;
    cOute[376] = 54.439308543349532953925518086180;
    cOute[377] = 62.039670232326493248820042936131;
    cOute[378] = 67.585933185905361142431502230465;
    cOute[379] = 71.678520063678661244921386241913;
    cOute[380] = 66.871240711537524248342378996313;
    cOute[381] = 54.775598007951899148793017957360;
    cOute[382] = 65.276830995904916221661551389843;
    cOute[383] = 68.536752618956541027728235349059;
    cOute[384] = 73.234497957691061742480087559670;
    cOute[385] = 68.215052994951818732261017430574;
    cOute[386] = 58.792592064162001008753577480093;
    cOute[387] = 53.540408959490683571402769302949;
    cOute[388] = 73.865450776805133159541583154351;
    cOute[389] = 85.144883547374135446261789184064;
    cOute[390] = 68.232547511409279650251846760511;
    cOute[391] = 56.049200682623165903351036831737;
    cOute[392] = 67.005635923898353212280198931694;
    cOute[393] = 70.560740980351411621995794121176;
    cOute[394] = 72.252403656408873189320729579777;
    cOute[395] = 69.758091503891719753482902888209;
    cOute[396] = 53.111935744363243827592668822035;
    cOute[397] = 65.764846623435516903555253520608;
    cOute[398] = 76.035642469221045303129358217120;
    cOute[399] = 76.967577078132762835593894124031;
    cOute[400] = 71.870439117577561205507663544267;
    cOute[401] = 59.868378263260105143217515433207;
    cOute[402] = 57.361967756528208894906128989533;
    cOute[403] = 72.481109157006457621719164308161;
    cOute[404] = 76.217001771805556131766934413463;
    cOute[405] = 70.235886792253822363818471785635;
    cOute[406] = 61.745385615736445572565571637824;
    cOute[407] = 55.826587976718165862166642909870;
    cOute[408] = 70.191628851956906487430387642235;
    cOute[409] = 87.739551028674426902398408856243;
    cOute[410] = 67.734712172386380757416191045195;
    cOute[411] = 52.632946435288268105523457052186;
    cOute[412] = 58.966014641825807984787388704717;
    cOute[413] = 69.803121352392253129437449388206;
    cOute[414] = 83.480415247375873377677635289729;
    cOute[415] = 72.322404079279692723503103479743;
    cOute[416] = 57.832870801888006440094613935798;
    cOute[417] = 51.894298857706971261904982384294;
    cOute[418] = 72.327377316208782076500938273966;
    cOute[419] = 88.357986031243939351043081842363;
    cOute[420] = 67.396331852535269035797682590783;
    cOute[421] = 56.356813485729205126517626922578;
    cOute[422] = 60.872766071927344455616548657417;
    cOute[423] = 71.749456177271653700699971523136;
    cOute[424] = 66.114325843371616997501405421644;
    cOute[425] = 63.343799066232669758846896002069;
    cOute[426] = 52.648168798288793368556071072817;
    cOute[427] = 51.443875275908482080922112800181;
    cOute[428] = 64.735223686734840953249658923596;
    cOute[429] = 91.100449400158240109703911002725;
    cOute[430] = 62.677102802594859554119466338307;
    cOute[431] = 58.633449842033058985180105082691;
    cOute[432] = 63.615813046687051723893091548234;
    cOute[433] = 74.709087568168712323313229717314;
    cOute[434] = 77.974458308863518141151871532202;
    cOute[435] = 44.015596017493038516477099619806;
    cOute[436] = 59.714610677715157294187520164996;
    cOute[437] = 57.302307043012149279093137010932;
    cOute[438] = 72.879140840712040017024264670908;
    cOute[439] = 78.728431445733207283410592935979;
    cOute[440] = 73.254909998541748450406885240227;
    cOute[441] = 68.797224958398913940982311032712;
    cOute[442] = 60.641614885752936459084594389424;
    cOute[443] = 73.545530249216326978967117611319;
    cOute[444] = 77.575089095287879104034800548106;
    cOute[445] = 75.700752880550837176087952684611;
    cOute[446] = 68.096865143616298610140802338719;
    cOute[447] = 56.457615183727959617954184068367;
    cOute[448] = 73.081624342450979270324751269072;
    cOute[449] = 82.089705888042502124335442204028;
    cOute[450] = 43.820158174482656932013924233615;
    cOute[451] = 55.245750169975011090173211414367;
    cOute[452] = 66.128540555033808345797297079116;
    cOute[453] = 78.461070396956529293674975633621;
    cOute[454] = 87.481780792978781846613856032491;
    cOute[455] = 38.146583549780359589931322261691;
    cOute[456] = 72.250575518071485703330836258829;
    cOute[457] = 71.388142420933007770145195536315;
    cOute[458] = 63.245855326992021616661077132449;
    cOute[459] = 78.570189200146174357541895005852;
    cOute[460] = 40.973425648725580572317994665354;
    cOute[461] = 57.724268859630058159382315352559;
    cOute[462] = 67.599098724516863967437529936433;
    cOute[463] = 74.747753019613313085756090003997;
    cOute[464] = 80.037228217565811405620479490608;
    cOute[465] = 43.135634466704075862253375817090;
    cOute[466] = 52.078644028807630661503935698420;
    cOute[467] = 75.929445840577344029043160844594;
    cOute[468] = 63.926507655362073023752600420266;
    cOute[469] = 90.121870632599097916681785136461;
    cOute[470] = 41.550335348083024200604995712638;
    cOute[471] = 40.710686113882928793827886693180;
    cOute[472] = 75.347519647929203756575589068234;
    cOute[473] = 78.786186635543771217271569184959;
    cOute[474] = 99.172941173288108984706923365593;
    cOute[475] = 42.192328231276249539405398536474;
    cOute[476] = 49.099406974782162649262318154797;
    cOute[477] = 83.113610952776269868991221301258;
    cOute[478] = 71.360482083929952068501734174788;
    cOute[479] = 72.283876870478181331236555706710;
    cOute[480] = 42.152993806190622194662864785641;
    cOute[481] = 54.260183893213351780104858335108;
    cOute[482] = 82.743323463905781522953475359827;
    cOute[483] = 73.717063821530715017615875694901;
    cOute[484] = 75.481361956171937777071434538811;
    cOute[485] = 43.020560018584646400086057838053;
    cOute[486] = 58.332083186160957666288595646620;
    cOute[487] = 80.661894422295645767917449120432;
    cOute[488] = 74.259013104655878123594447970390;
    cOute[489] = 77.899566821473513300588820129633;
    cOute[490] = 43.755557908754738605239253956825;
    cOute[491] = 58.341019843528833632717578439042;
    cOute[492] = 73.575473933789510283531853929162;
    cOute[493] = 75.359368526781054242746904492378;
    cOute[494] = 63.422483437959265017980214906856;
    cOute[495] = 44.556839960036818126809521345422;
    cOute[496] = 44.717096370293226925696217222139;
    cOute[497] = 68.664890452024735623126616701484;
    cOute[498] = 85.381102695340032937565410975367;
    cOute[499] = 82.916651793023262939641426783055;
    cOute[500] = 72.692095008882361639734881464392;
    cOute[501] = 56.134960422275788971546717220917;
    cOute[502] = 69.379047314506223642638360615820;
    cOute[503] = 77.297727959998255187201721128076;
    cOute[504] = 86.206771701997396917249716352671;
    cOute[505] = 72.944923150640491371632379014045;
    cOute[506] = 58.378553321407508747142856009305;
    cOute[507] = 71.967566164041826937136647757143;
    cOute[508] = 75.591562245406208830900141038001;
    cOute[509] = 78.822359754279801791199133731425;
    cOute[510] = 77.726342068369632443136652000248;
    cOute[511] = 62.082293420644376169548195321113;
    cOute[512] = 70.244791104619508814721484668553;
    cOute[513] = 83.328412048612079843223909847438;
    cOute[514] = 91.299051026264180563885020092130;
    cOute[515] = 75.448343233073188685011700727046;
    cOute[516] = 56.269094749874767558139865286648;
    cOute[517] = 75.937645991740382100942952092737;
    cOute[518] = 79.736167258786551315097312908620;
    cOute[519] = 90.810978116526328562940761912614;
    cOute[520] = 76.240068247248601096544007305056;
    cOute[521] = 54.673506626703670008282642811537;
    cOute[522] = 70.754952415659417397364450152963;
    cOute[523] = 77.498230466999629584279318805784;
    cOute[524] = 87.859472260882157002015446778387;
    cOute[525] = 75.257322032021491509112820494920;
    cOute[526] = 60.862516744076870622848218772560;
    cOute[527] = 75.038478649707201384444488212466;
    cOute[528] = 77.461280194313900437919073738158;
    cOute[529] = 75.439185702515203502116492018104;
    cOute[530] = 70.555551016096501371066551655531;
    cOute[531] = 65.767915596590427185219596140087;
    cOute[532] = 73.150425280480448009257088415325;
    cOute[533] = 77.849603712172836367244599387050;
    cOute[534] = 80.030151714913131399953272193670;
    cOute[535] = 93.565542492623563930465024895966;
    cOute[536] = 68.771484503404423094252706505358;
    cOute[537] = 70.017031721689349410553404595703;
    cOute[538] = 94.980145416706577066179306712002;
    cOute[539] = 106.708154287972163842823647428304;
    cOute[540] = 72.330982575775337295453937258571;
    cOute[541] = 58.430605974477735742311779176816;
    cOute[542] = 73.522831512454203561901522334665;
    cOute[543] = 77.750650422131869277109217364341;
    cOute[544] = 90.560722017563804797646298538893;
    cOute[545] = 75.566481712725803276953229214996;
    cOute[546] = 58.679916883113712344766099704430;
    cOute[547] = 72.376199982647534625357366167009;
    cOute[548] = 77.432896285548167725210078060627;
    cOute[549] = 78.109261117079597624979214742780;
    cOute[550] = 88.777196288847818550493684597313;
    cOute[551] = 67.139471538076236356573645025492;
    cOute[552] = 70.499150001211091876029968261719;
    cOute[553] = 79.317511262396806159813422709703;
    cOute[554] = 95.437969409844612300730659626424;
    cOute[555] = 89.581239521959474814138957299292;
    cOute[556] = 74.604040451150396506818651687354;
    cOute[557] = 82.287703596196976718601945322007;
    cOute[558] = 90.746520708060145921081129927188;
    cOute[559] = 104.472155611867705715667398180813;
    cOute[560] = 95.418704747189039494514872785658;
    cOute[561] = 77.556007354830185818173049483448;
    cOute[562] = 81.553399846605387324416369665414;
    cOute[563] = 99.282262402813600488116207998246;
    cOute[564] = 100.232458421756533084590046200901;
    cOute[565] = 76.326166148501343400312180165201;
    cOute[566] = 73.182663806689845387154491618276;
    cOute[567] = 78.665848298534172045037848874927;
    cOute[568] = 94.798670259896425704937428236008;
    cOute[569] = 104.925929703516374047467252239585;
    cOute[570] = 76.559104285052683280810015276074;
    cOute[571] = 66.015191603215924942560377530754;
    cOute[572] = 62.500316077818979465519078075886;
    cOute[573] = 71.422317569645684898205217905343;
    cOute[574] = 88.029358394304836110677570104599;
    cOute[575] = 80.679799259197707783641817513853;
    cOute[576] = 101.563392126254925074135826434940;
    cOute[577] = 84.129094099643637605367985088378;
    cOute[578] = 96.821797750537186288966040592641;
    cOute[579] = 109.465233806279158557117625605315;
    cOute[580] = 79.287490219349066933318681549281;
    cOute[581] = 94.968386774647328252285660710186;
    cOute[582] = 72.089540611198273722948215436190;
    cOute[583] = 87.620942773377223034003691282123;
    cOute[584] = 102.978757155602252737480739597231;
    cOute[585] = 78.057320759489272177233942784369;
    cOute[586] = 91.793874606967108320532133802772;
    cOute[587] = 78.582982738882577677941299043596;
    cOute[588] = 83.026918121939957018184941262007;
    cOute[589] = 98.551571656425664968992350623012;
    cOute[590] = 83.041568291603425677749328315258;
    cOute[591] = 87.706319866641564431120059452951;
    cOute[592] = 79.800793603055012681579682976007;
    cOute[593] = 88.818645502041590589215047657490;
    cOute[594] = 105.511121785602682621174608357251;
    cOute[595] = 79.017271954077855866671598050743;
    cOute[596] = 91.399076427102002639912825543433;
    cOute[597] = 77.803664188620061281653761398047;
    cOute[598] = 95.426657264353977438986476045102;
    cOute[599] = 114.754378652918816783312649931759;
    cOute[600] = 88.311638644853800883538497146219;
    cOute[601] = 93.381651657871586280634801369160;
    cOute[602] = 86.612837763985382366627163719386;
    cOute[603] = 90.317482845097615040685923304409;
    cOute[604] = 89.973590122006754654648830182850;
    cOute[605] = 83.471417562136906553860171698034;
    cOute[606] = 88.421350421341628589289030060172;
    cOute[607] = 84.646418752859943879229831509292;
    cOute[608] = 88.051426395289439597036107443273;
    cOute[609] = 88.079293876322935830103233456612;
    cOute[610] = 75.227206753133742722638999111950;
    cOute[611] = 92.186353720333386263519059866667;
    cOute[612] = 99.346685048972830145430634729564;
    cOute[613] = 82.552257783567768001375952735543;
    cOute[614] = 95.384495059268218142278783489019;
    cOute[615] = 77.176160126025294516693975310773;
    cOute[616] = 85.982874722011828794165921863168;
    cOute[617] = 97.688965233044129377049102913588;
    cOute[618] = 80.767414089336043048206192906946;
    cOute[619] = 93.522272918860934964868647512048;
    cOute[620] = 67.994238638390655182774935383350;
    cOute[621] = 64.749147652344689163328439462930;
    cOute[622] = 86.089780392842314427070959936827;
    cOute[623] = 93.351058184660317351699632126838;
    cOute[624] = 97.086036434254140203847782686353;
}
