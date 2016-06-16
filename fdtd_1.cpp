#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <direct.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"


//サブルーティン
void file_open(char*);
void file_close();
void parameter(char*);
void initialize_matrix();
void modeling();
void set_epsilon();
void source_func();
void observation_func();
void calc_efield();
void calc_hfield();
void absorpt_bound_condition();
void saving_electric_field();
//void magnitude();
void output_field(char*);
void output_field_write(char *);
void output_model();
void calc_energy();
void calc_power();
void calc_powerHz();
void calc_poynting_power();
void calc_poynting_powerHz();

void mcircle(int, int, double);  //make circle function
void rightquartercircle1(int, int, double);
void leftquartercircle1(int, int, double);
void rightquartercircle2(int, int, double);
void leftquartercircle2(int, int, double);

FILE *model_xy;
FILE *model_yz;
FILE *model_xz;
char *dir_name[] = {"1550"};


#define INT_DIV(X,Y)   ((int)((X)/(Y)))
#define SQ(X)          ((X)*(X))

static void _strtime(char *buf)
{
    time_t     current;
    struct tm  *local;
    time(&current);                     /* 現在の時刻を取得 */
    local = localtime(&current);        /* 地方時の構造体に変換 */
    asctime_r(local, buf);
}

/* MPI global parameter */
#define NODE 2
#define IRANK_MIN 0				// 最小値(この数字の割り当てられた計算機に結果が入る)
static int IRANK_MAX = NODE-1;
static int irank = 0;
//static int isize = 0;

/* Scale */

#define XMAX_ALL    100
#define XMAX (INT_DIV(XMAX_ALL, NODE) + 1)
#define YMAX        50
//#define ZMAX        50
#define RX(x)       ((irank * XMAX) + x)

#define dt       (28e-18) 		// 時間ステップ[s]
#define Nmax     (300) 		// 最終時間ステップ

#define Ncut 500 			// 時間ステップ数を表示させる間隔
#define Tcut 200 			// エネルギーの平均の算出を開始する時間ステップ
#define Fcut 200 			// フィールドを出力する時間ステップ数 (最終計算ステップからの差)

static int xmax, ymax, zmax;
//#define ymax YMAX
int n;
//#define zmax ZMAX

#define CELL_SIZE 100  // セルサイズ[nm]

#define dblCellSize (CELL_SIZE * 1e-9)
#define dx dblCellSize
#define dy dblCellSize
#define dz dblCellSize
#define inv_dx (1/dx)
#define inv_dy (1/dy)
#define inv_dz (1/dz)


#define intCladHeight_top 10
#define intCladHeight_bottom 10
#define intSlabHeight 20
#define int_air 10

#define ZMAX (intCladHeight_top+intCladHeight_bottom+intSlabHeight+int_air)
//#define zmax ZMAX

#define PitchX 400		//[nm]
#define intPitchX  (INT_DIV(PitchX,CELL_SIZE)
#define intPitchY  (INT_DIV((PitchX * sqrt(3.0)/2 + 0.5), CELL_SIZE)); 		// 格子定数のセルサイズ(Y方向)	+0.5は四捨五入のため

#define intPcwLine (INT_DIV(CELL_SIZE*XMAX_ALL,PitchX))
#define intPcwRow  6

#define RADIUS 100	//円孔半径[nm]

#define dblRadius (RADIUS * 1e-9) 					// 円孔半径


/*モデル*/
#define CLAD 0									// クラッド
#define CORE 1									// コア
#define AIR 3									// 空気
#define CIRCLE_REF_INDEX CORE					// 円孔

/*****************************************************************************/
// 物理量[MKSA系]
/*****************************************************************************/

#define PI 3.141592
#define C0 2.997924e8			// 真空中の光速 [m/s]
#define EPSILON0 8.854e-12		// 真空中の誘電率 [F/m]
#define MU0 (PI*(4.0e-7))		// 真空中の透磁率 [N/A^2]

static const double cnstHxyz = dt / MU0; //磁界の計算に使う定数

static const double n_core = 3.43; 		// コアの屈折率(CORE)
static const double n_clad = 1.444; 		// クラッド屈折率(SiO2)

#define epsilon0  EPSILON0 					// 真空の誘電率
#define epsilon1  (EPSILON0 * SQ(n_core)) 		// コアの誘電率(Si)
#define epsilon2  (EPSILON0 * SQ(n_clad)) 		// クラッドの誘電率(SiO2)


// 励振関数定数の設定
//lambda = atof(dir_name) * 1e-9;
//omega0 = 2.0*PI*C0/lambda;
//sigma = omega0 * delta_omega;

//double lambda; 					// 励振関数の波長 [m] (プログラム中で dir_name をdouble化して代入)
//double omega0; 					// 励振関数の角周波数 [s^-1]
// ガウシアンパルス
//double sigma; 					// 広がり幅を決定する定数

static double Ex[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
static double Ey[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
static double Ez[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
static double Hx[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
static double Hy[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
static double Hz[XMAX+1][YMAX+1][ZMAX+1] = {0.0};

// 誘電率 (YMAX+1, ZMAX+1は対称境界条件で必要なため)

static double epsilonx[XMAX][YMAX+1][ZMAX+1] = {0.0};
static double epsilony[XMAX][YMAX+1][ZMAX+1] = {0.0};
static double epsilonz[XMAX][YMAX+1][ZMAX+1] = {0.0};
static int cell[XMAX][YMAX+1][ZMAX+1] = {0};

// 吸収境界条件適用のときに使う電界の配列
double
    Exn2y00[XMAX+1][ZMAX+1] = {0.0},
    Exn1y00[XMAX+1][ZMAX+1+1] = {0.0},
    Exn2y01[XMAX+1][ZMAX+1] = {0.0},
    Exn1y01[XMAX+1][ZMAX+1+1] = {0.0};
double
    Exn2z00[XMAX+1][YMAX+1] = {0.0},
    Exn1z00[XMAX+1][YMAX+1+1] = {0.0},
    Exn2z01[XMAX+1][YMAX+1] = {0.0},
    Exn1z01[XMAX+1][YMAX+1+1] = {0.0};

double
    Eyn2z00[XMAX+1][YMAX] = {0.0},
    Eyn1z00[XMAX+1][YMAX+1] = {0.0},
    Eyn2z01[XMAX+1][YMAX] = {0.0},
    Eyn1z01[XMAX+1][YMAX+1] = {0.0};

double
    Eyn2x00[YMAX][ZMAX+1] = {0.0},
    Eyn1x00[YMAX+1][ZMAX+1+1] = {0.0},
    Eyn2x01[YMAX][ZMAX+1] = {0.0},
    Eyn1x01[YMAX+1][ZMAX+1+1] = {0.0};

double Eyn2xm1[YMAX][ZMAX+1] = {0.0},
    Eyn1xm1[YMAX+1][ZMAX+1+1] = {0.0},
    Eyn2xm0[YMAX][ZMAX+1] = {0.0},
    Eyn1xm0[YMAX+1][ZMAX+1+1] = {0.0};

double
    Ezn2y00[XMAX+1][ZMAX] = {0.0},
    Ezn1y00[XMAX+1][ZMAX+1] = {0.0},
    Ezn2y01[XMAX+1][ZMAX] = {0.0},
    Ezn1y01[XMAX+1][ZMAX+1] = {0.0};

double
    Ezn2x00[YMAX+1][ZMAX] = {0.0},
    Ezn1x00[YMAX+1+1][ZMAX+1] = {0.0},
    Ezn2x01[YMAX+1][ZMAX] = {0.0},
    Ezn1x01[YMAX+1+1][ZMAX+1] = {0.0};

double
    Ezn2xm1[YMAX+1][ZMAX] = {0.0},
    Ezn1xm1[YMAX+1+1][ZMAX+1] = {0.0},
    Ezn2xm0[YMAX+1][ZMAX] = {0.0},
    Ezn1xm0[YMAX+1+1][ZMAX+1] = {0.0};

/*
// 磁界分布出力用
double field_xy[XMAX][YMAX]; 	// Hz-field のファイル出力 (面垂直方向の磁界成分)
double field_yz[YMAX][ZMAX];	// Hz-field のファイル出力 (面方向の磁界成分)
double field_zx_Hz[ZMAX][XMAX];	// Hy-field のファイル出力 (面方向の磁界成分)
double field_zx_Hy[ZMAX][XMAX];	// Hy-field のファイル出力 (面方向の磁界成分)

int x, y, z; 				// 離散座標
int xmax, ymax, zmax; 			// 空間分割の最大値
int xmax_all; //分割前の最大値
//int n, Nmax; 				//n: 時間ステップ，Nmax: 時間ステップの最大値
int n; 					//時間ステップ数
//int icut, jcut, kcut; 		//icut, jcut, kcut: フィールド出力のときに使う
//int PrintStat;
//int PrintEnd;
//int source_i1, source_i2;
//int source_j1, source_j2;
//int source_k;
int x_cen, y_cen, z_cen; 			//x_cen, y_cen, z_cen: 解析空間の中心離散座標
int x_model_cen, y_model_cen; 	//x_model_cen, y_model_cen:モデルの中心離散座標
*/

void file_open(char* dir_name_def){
	char dir_name[40];
	//_mkdir(strcpy(dir_name, dir_name_def)); 		// 振り分けできるかテスト

	model_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xy.txt"), "w"); 		// 振り分けできるかテスト
	model_yz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_yz.txt"), "w");
	model_xz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xz.txt"), "w");

}

/*出力用ファイルを閉じる*/
void file_close(){

	fclose(model_xy);
	fclose(model_yz);
	fclose(model_xz);

}



/*配列の初期化*/
void initialize_matrix()
{
    /*
    //最後のノードだけのりしろ不要なのでx方向に1セル小さい
    xmax = (irank != IRANK_MAX) ? XMAX : (XMAX-1);

    // 解析空間の中心座標
    x_cen = xmax/2;
    y_cen = ymax/2;
    z_cen = zmax/2;

    //モデルの中心と解析空間の中心は１セル分ずれているので要注意
    x_model_cen = x_cen + 1;
    y_model_cen = y_cen + 1;

    for (x = 0; x < xmax_all; x++){
        for(y = 0; y < ymax+1; y++){
            for(z = 0; z < zmax+1; z++){
                // 誘電率
                ALL_epsilonx[x][y][z] = EPSILON0;
                ALL_epsilony[x][y][z] = EPSILON0;
                ALL_epsilonz[x][y][z] = EPSILON0;
            }
        }
    }

    for(x = 0; x < xmax+1; x++){
        for(y = 0; y < ymax+1; y++){
            for(z = 0; z < zmax+1; z++){
                // 誘電率(分割する必要がないような．．．)
                epsilonx[x][y][z] = EPSILON0;
                epsilony[x][y][z] = EPSILON0;
                epsilonz[x][y][z] = EPSILON0;
            }
        }
    }

    for(x = 0; x < xmax_all; x++){
        for(y = 0; y < ymax+1; y++){
            for(z = 0; z < zmax+1; z++){
                // セルの目印
                cell[RX(x)][y][z] = CLAD;
            }
        }
    }

    for(x = 0; x < xmax; x++){
        for(y = 0; y < ymax; y++){
            for(z = 0; z < zmax; z++){
                // セルの目印
                cell[x][y][z] = 0;
            }
        }
    }*/

}

void modeling()
{
  if(irank != IRANK_MAX){
    xmax = XMAX;
    ymax = YMAX;
    zmax = ZMAX;
  }

  //最後のノードだけのりしろ不要なのでx方向に1セル小さい
  if(irank == IRANK_MAX){
    xmax = XMAX - 1;
    ymax = YMAX;
    zmax = ZMAX;
  }

	int n_temp; 		//屈折率の値保存用
	double epsilon_temp; 		//誘電率の値保存用
  int x,y,z;
	/****************************** スラブの形成 ******************************/
	for(x = 0; x < xmax; x++){
		for(y = 0; y < YMAX; y++){
			for(z = 0; z < ZMAX; z++){

				n_temp = CLAD;
				epsilon_temp = epsilon2;

				if(z < intCladHeight_bottom){			//下部クラッド層に設定
				}
				if(z >= intCladHeight_bottom && z < (intSlabHeight + intCladHeight_bottom) ){			//スラブに設定
					n_temp = CORE;
					epsilon_temp = epsilon1;
				}
				if(z >= (intSlabHeight + intCladHeight_bottom) && z < (intSlabHeight + intCladHeight_bottom + intCladHeight_top) ){			//上部クラッド層に設定
				}
				if(z >= (intSlabHeight + intCladHeight_bottom + intCladHeight_top) && z < ZMAX){	// スラブに設定
					n_temp = AIR;
					epsilon_temp = epsilon0;
				}

				cell[x][y][z] = n_temp;
				epsilonx[x][y][z] = epsilon_temp;
				epsilony[x][y][z] = epsilon_temp;
				epsilonz[x][y][z] = epsilon_temp;
			}
		}
	 }

	/****************************** スラブの形成 ******************************/


	/****************************** フォトニック結晶 ******************************/
	/*
	for(x = 0; x < intPcwLine; x++){
		int y2, y_poo;
		if (x < input_NormPcw_Xend){
			if (x == 0){
				y_poo = 0;
				for (y2 = intPcwRow-1; y2 >= 0; y2--){
					Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす
					if(y2 % 2 == 1){		// if y:even
						Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchX * x - 1;	// 配列の引数に使用するので-1
					}
					else{				// if y:odd 0.5Aずらす
						Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
					}
					y_poo++;
				}
			}
			else{
				for (y2 = intPcwWid-1; y2 >= 0; y2--){
					Pnum[x][y2].Y = Pnum[x-1][y2].Y;
					Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;
				}
			}
		}
	}

	double R;

	//半径セル数の計算
	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています

	mcircle(Pnum[x][y].X, Pnum[x][y].Y, R);
*/

	/****************************** フォトニック結晶 ******************************/



	/****************************** 対称境界部分の誘電率の設定 ******************************/
/*
	for(x = 0; x < xmax_all+1; x++){
		for(z = 0; z < zmax_all+1; z++){
			ALL_cell[x][ymax][z] = ALL_cell[x][ymax-1][z];
			ALL_epsilonx[x][ymax][z] = ALL_epsilonx[x][ymax-1][z];
			ALL_epsilony[x][ymax][z] = ALL_epsilony[x][ymax-1][z];
			ALL_epsilonz[x][ymax][z] = ALL_epsilonz[x][ymax-1][z];
		}
	}

	for(x = 0; x < xmax_all+1; x++){
		for(y = 0; y < ymax_all+1; y++){
			ALL_cell[x][y][zmax] = ALL_cell[x][y][zmax-1];
			ALL_epsilonx[x][y][zmax] = ALL_epsilonx[x][y][zmax-1];
			ALL_epsilony[x][y][zmax] = ALL_epsilony[x][y][zmax-1];
			ALL_epsilonz[x][y][zmax] = ALL_epsilonz[x][y][zmax-1];
		}
	}
*/
	/****************************** 対称境界部分の誘電率の設定 ******************************/





}


/****************************** 円孔の生成 ******************************/
/*
void mcircle(int x_circ, int y_circ, double R){
	rightquartercircle1(x_circ, y_circ, R);
	leftquartercircle1(x_circ-1, y_circ, R);
	rightquartercircle2(x_circ, y_circ-1, R);
	leftquartercircle2(x_circ-1, y_circ-1, R);
}

void rightquartercircle1(int x_circ, int y_circ, double R){

	int x, y, Ie, Je;
	double r;

	Ie = (int) (x_circ+R-1);
	Je = (int) (y_circ+R-1);
	for(x = x_circ; x <= Ie; x++){
		for(y = y_circ; y <= Je; y++){
			for(z >= intCladHeight_bottom && z < (intSlabHeight + intCladHeight_bottom){
				r = sqrt(double((x-x_circ+1) * (x-x_circ+1) + (y-y_circ+1) * (y-y_circ+1)) ) - 0.5;
				if(r <= R){
					ALL_cell[x][y][z] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}
		}
	}
}


void leftquartercircle1(int x_circ, int y_circ, double R){

	int x, y, Ie, Je;
	double r;

	Ie = (int) (x_circ-R+1);
	Je = (int) (y_circ+R-1);
	for(x = x_circ; x >= Ie; x--){
		for(y = y_circ; y <= Je; y++){
			for(z >= intCladHeight_bottom && z < (intSlabHeight + intCladHeight_bottom){
			r = sqrt(double((x-x_circ-1) * (x-x_circ-1) + (y-y_circ+1) * (y-y_circ+1)) ) - 0.5;
				if(r <= R){
					ALL_cell[x][y][z] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}
		}
	}
}


void rightquartercircle2(int x_circ, int y_circ, double R){

	int x, y, Ie, Je;
	double r;

	Ie = (int) (x_circ+R-1);
	Je = (int) (y_circ-R+1);
	for(x = x_circ; x <= Ie; x++){
		for(y = y_circ; y >= Je; y--){
			for(z >= intCladHeight_bottom && z < (intSlabHeight + intCladHeight_bottom){
				r = sqrt(double((x-x_circ+1) * (x-x_circ+1) + (y-y_circ-1) * (y-y_circ-1)) ) - 0.5;
				if(r <= R){
					ALL_cell[x][y][z] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}
		}
	}
}

void leftquartercircle2(int x_circ, int y_circ, double R){

	int x, y, Ie, Je;
	double r;

	Ie = (int) (x_circ-R+1);
	Je = (int) (y_circ-R+1);
	for(x = x_circ; x >= Ie; x--){
		for(y = y_circ; y >= Je; y--){
			for(z >= intCladHeight_bottom && z < (intSlabHeight + intCladHeight_bottom){
				r = sqrt(double((x-x_circ-1) * (x-x_circ-1) + (y-y_circ-1) * (y-y_circ-1)) ) - 0.5;
				if(r <= R){
					ALL_cell[x][y][z] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}
		}
	}
}
*/

/****************************** 円孔の生成 ******************************/





// 励振関数
void source_func(int n)
{
    /*
    static const double delta_omega = 0.05; 				// 中心周波数で規格した半値全幅
    static const int Npeak = 500; 							// ピークステップ数

    // 励振点の座標☆☆
    static const int ex_y_st = YMAX -24 ; 	//- intWireWid_2	// 導波路断面始セル数(横) ←解析空間の中間セル座標から導波路幅(既に1/2値になっている)を引いている☆☆
    static const int ex_y_ed = YMAX; 					// 導波路断面終セル数(横)
    static const int ex_z_st = ZMAX - intSlabHeigPer; 	// 導波路断面始セル数(縦) ←解析空間の中間セル座標から導波路幅(の1/2値)を引いている
    static const int ex_z_ed = ZMAX; 			// 導波路断面終セル数(縦)

    if(irank == intExctPortNum){
        // 励振点の設定
        int x = intExctLenPart;
        for(int y = ex_y_st; y < ex_y_ed; y++){
            for(int z = ex_z_st; z < ex_z_ed; z++){
#if _EXITATION_FUNC	// CW励振
                Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt); // 01
#else	// Gaussian励振
                Hz[x][y][z] += 1000 * cos(omega0*(n-Npeak)*dt) * exp(-(SQ(sigma*dt*(n-Npeak))/2)) * cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st));
#endif
            }
        }
    }
     */

    for(int x = 0; x < xmax+1; x++){
        for(int z = 0; z < zmax; z++){
            Hx[x][ymax][z] = Hx[x][ymax-1][z];		// 偶関数
        }
    }
    for(int x = 0; x < xmax; x++){
        for(int z = 0; z < zmax+1; z++){
            Hz[x][ymax][z] = Hz[x][ymax-1][z];		// 偶関数
        }
    }
    for(int x = 0; x < xmax; x++){
        for(int y = 0; y < ymax+1; y++){
            Hy[x][y][zmax] = -Hy[x][y][zmax-1];		// 奇関数
        }
    }
    for(int x = 0; x < xmax+1; x++){
        for(int y = 0; y < ymax; y++){
            Hx[x][y][zmax] = -Hx[x][y][zmax-1];		// 奇関数
        }
    }
}

void calc_efield()
{
    // Ex
    for(int x = 0; x < xmax; x++){
        for(int y = 1; y < ymax+1; y++){		// Exはy軸に対して奇関数
            for(int z = 1; z < zmax+1; z++){
                double cnstEx = dt / epsilonx[x][y][z];
                double dex = ( (Hz[x][y][z] - Hz[x][y-1][z]) / dy) - ( (Hy[x][y][z] - Hy[x][y][z-1]) / dz);
                Ex[x][y][z] = Ex[x][y][z] + cnstEx * dex;
            }
        }
    }

    // Ey
    for(int x = 1; x < xmax; x++){
        for(int y = 0; y < ymax; y++){
            for(int z = 1; z < zmax+1; z++){
                double cnstEy = dt / epsilony[x][y][z];
                double dey = ( (Hx[x][y][z] - Hx[x][y][z-1]) / dz)-( (Hz[x][y][z] - Hz[x-1][y][z]) / dx);
                Ey[x][y][z] = Ey[x][y][z] + cnstEy * dey;
            }
        }
    }

    // Ez
    for(int x = 1; x < xmax; x++){
        for(int y = 1; y < ymax+1; y++){		// Ezはy軸に対して奇関数
            for(int z = 0; z < zmax; z++){
                double cnstEz = dt / epsilonz[x][y][z];
                double dez = ( (Hy[x][y][z] - Hy[x-1][y][z]) / dx) - ( (Hx[x][y][z] - Hx[x][y-1][z]) / dy);
                Ez[x][y][z] = Ez[x][y][z] + cnstEz * dez;
            }
        }
    }

    /****************************** 電界の対称境界条件 ******************************/

    // 境界面で反対称となる電界成分の境界面上の値を0としている
    for(int x = 0; x < xmax; x++){
        for(int z = 0; z < zmax+1; z++){
            Ex[x][ymax][z] = 0.0;		// 奇関数
        }
    }
    for(int x = 0; x < xmax+1; x++){
        for(int z = 0; z < zmax; z++){
            Ez[x][ymax][z] = 0.0;		// 奇関数
        }
    }
    /****************************** 電界の対称境界条件 ******************************/
}



void calc_hfield()
{
    // Hx
    for(int x = 0; x < xmax+1; x++){
        for(int y = 0; y < ymax; y++){
            for(int z = 0; z < zmax; z++){
                double dhx = ( (Ey[x][y][z+1] - Ey[x][y][z]) / dz) - ( (Ez[x][y+1][z] - Ez[x][y][z]) / dy);
                Hx[x][y][z] = Hx[x][y][z] + cnstHxyz * dhx;
            }
        }
    }

    // Hy
    for(int x = 0; x < xmax; x++){
        for(int y = 0; y < ymax+1; y++){
            for(int z = 0; z < zmax; z++){
                double dhy = ( (Ez[x+1][y][z] - Ez[x][y][z]) / dx) - ( (Ex[x][y][z+1] - Ex[x][y][z]) / dz);
                Hy[x][y][z] = Hy[x][y][z] + cnstHxyz * dhy;
            }
        }
    }

    // Hz
    for(int x = 0; x < xmax; x++){
        for(int y = 0; y < ymax; y++){
            for(int z = 0; z < zmax+1; z++){
                double dhz = ((Ex[x][y+1][z] - Ex[x][y][z]) / dy) - ((Ey[x+1][y][z] - Ey[x][y][z]) / dx);
                Hz[x][y][z] = Hz[x][y][z] + cnstHxyz * dhz;
            }
        }
    }
}


// Mur2次，1次の吸収境界条件から端面の計算
void absorpt_bound_condition()
{

    /****************************** 対称境界条件 ******************************/

    // 4回対称
    for(int z = 0; z < zmax+1; z++){
        Eyn1x00[ymax][z] = Eyn1x00[ymax-1][z];
        Eyn1x01[ymax][z] = Eyn1x01[ymax-1][z];
        Eyn1xm0[ymax][z] = Eyn1xm0[ymax-1][z];
        Eyn1xm1[ymax][z] = Eyn1xm1[ymax-1][z];
    }
    for(int x = 0; x < xmax+1; x++){
        Eyn1z00[x][ymax] = Eyn1z00[x][ymax-1];
        Eyn1z01[x][ymax] = Eyn1z01[x][ymax-1];
    }
    for(int x = 0; x < xmax; x++){
        Exn1z00[x][ymax+1] = -Exn1z00[x][ymax-1];
        Exn1z01[x][ymax+1] = -Exn1z01[x][ymax-1];
    }
    for(int z = 0; z < zmax; z++){
        Ezn1x00[ymax+1][z] = -Ezn1x00[ymax-1][z];
        Ezn1x01[ymax+1][z] = -Ezn1x01[ymax-1][z];
        Ezn1xm0[ymax+1][z] = -Ezn1xm0[ymax-1][z];
        Ezn1xm1[ymax+1][z] = -Ezn1xm1[ymax-1][z];
    }

    // 8回対称
    for(int y = 0; y < ymax+1; y++){
        Ezn1x00[y][zmax] = -Ezn1x00[y][zmax-1];
        Ezn1x01[y][zmax] = -Ezn1x01[y][zmax-1];
        Ezn1xm0[y][zmax] = -Ezn1xm0[y][zmax-1];
        Ezn1xm1[y][zmax] = -Ezn1xm1[y][zmax-1];
    }
    for(int x = 0; x < xmax+1; x++){
        Ezn1y00[x][zmax] = -Ezn1y00[x][zmax-1];
        Ezn1y01[x][zmax] = -Ezn1y01[x][zmax-1];
    }
    for(int x = 0; x < xmax; x++){
        Exn1y00[x][zmax+1] = Exn1y00[x][zmax-1];
        Exn1y01[x][zmax+1] = Exn1y01[x][zmax-1];
    }
    for(int y = 0; y <= ymax-1; y++){
        Eyn1x00[y][zmax+1] = Eyn1x00[y][zmax-1];
        Eyn1x01[y][zmax+1] = Eyn1x01[y][zmax-1];
        Eyn1xm0[y][zmax+1] = Eyn1xm0[y][zmax-1];
        Eyn1xm1[y][zmax+1] = Eyn1xm1[y][zmax-1];
    }

    /****************************** 対称境界条件 ******************************/

    /****************************** Murの2次の吸収境界条件(Ex) ******************************/

    if(irank != IRANK_MAX){
        for(int x = 1; x < xmax; x++){
            for(int z = 1; z < zmax+1; z++){
                double velo_dt = (C0 / sqrt(epsilonx[x][0][z]/epsilon0) ) * dt;

                double u1ax1 = (velo_dt - dy) / (velo_dt + dy);
                double u2ax1 = (2.0 * dy) / (velo_dt + dy);
                double u3ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy) );
                double u4ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy) );

                Ex[x][0][z] = -Exn2y01[x][z]
                +u1ax1 * (Ex[x][1][z] + Exn2y00[x][z])
                +u2ax1 * (Exn1y00[x][z] + Exn1y01[x][z])
                +u3ax1 * (Exn1y00[x+1][z] - 2.0 * Exn1y00[x][z] + Exn1y00[x-1][z] + Exn1y01[x+1][z] - 2.0 * Exn1y01[x][z] + Exn1y01[x-1][z])
                +u4ax1 * (Exn1y00[x][z+1] - 2.0 * Exn1y00[x][z] + Exn1y00[x][z-1] + Exn1y01[x][z+1] - 2.0 * Exn1y01[x][z] + Exn1y01[x][z-1]);
            }
        }
        for(int x = 1; x < xmax; x++){
            for(int y = 1; y < ymax+1; y++){
                double velo_dt = (C0 / sqrt(epsilonx[x][y][0]/epsilon0) ) * dt;

                double u1bx1 = (velo_dt - dz) / (velo_dt + dz);
                double u2bx1 = (2.0 * dz) / (velo_dt + dz);
                double u3bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) );
                double u4bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) );

                Ex[x][y][0] = -Exn2z01[x][y]
                +u1bx1 * (Ex[x][y][1] + Exn2z00[x][y])
                +u2bx1 * (Exn1z00[x][y] + Exn1z01[x][y])
                +u3bx1 * (Exn1z00[x+1][y] - 2.0 * Exn1z00[x][y] + Exn1z00[x-1][y] + Exn1z01[x+1][y] - 2.0*Exn1z01[x][y] + Exn1z01[x-1][y])
                +u4bx1 * (Exn1z00[x][y+1] - 2.0 * Exn1z00[x][y] + Exn1z00[x][y-1] + Exn1z01[x][y+1] - 2.0*Exn1z01[x][y] + Exn1z01[x][y-1]);
            }
        }
    }
    else{
        for(int x = 1; x < xmax-1; x++){
            for(int z = 1; z < zmax+1; z++){
                double velo_dt = (C0 / sqrt(epsilonx[x][0][z]/epsilon0) ) * dt;

                double u1ax1 = (velo_dt - dy) / (velo_dt + dy);
                double u2ax1 = (2.0 * dy) /(velo_dt + dy);
                double u3ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy));
                double u4ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy));

                Ex[x][0][z] = -Exn2y01[x][z]
                +u1ax1 * (Ex[x][1][z] + Exn2y00[x][z])
                +u2ax1 * (Exn1y00[x][z] + Exn1y01[x][z])
                +u3ax1 * (Exn1y00[x+1][z] - 2.0 * Exn1y00[x][z] + Exn1y00[x-1][z] + Exn1y01[x+1][z] - 2.0 * Exn1y01[x][z] + Exn1y01[x-1][z])
                +u4ax1 * (Exn1y00[x][z+1] - 2.0 * Exn1y00[x][z] + Exn1y00[x][z-1] + Exn1y01[x][z+1] - 2.0 * Exn1y01[x][z] + Exn1y01[x][z-1]);
            }
        }
        for(int x = 1; x < xmax-1; x++){
            for(int y = 1; y < ymax+1; y++){
                double velo_dt = (C0 / sqrt(epsilonx[x][y][0]/epsilon0) ) * dt;

                double u1bx1 = (velo_dt - dz) / (velo_dt + dz);
                double u2bx1 = (2.0 * dz) / (velo_dt + dz);
                double u3bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) );
                double u4bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) );

                Ex[x][y][0] = -Exn2z01[x][y]
                +u1bx1 * (Ex[x][y][1] + Exn2z00[x][y])
                +u2bx1 * (Exn1z00[x][y] + Exn1z01[x][y])
                +u3bx1 * (Exn1z00[x+1][y] - 2.0 * Exn1z00[x][y] + Exn1z00[x-1][y] + Exn1z01[x+1][y] - 2.0 * Exn1z01[x][y] + Exn1z01[x-1][y])
                +u4bx1 * (Exn1z00[x][y+1] - 2.0 * Exn1z00[x][y] + Exn1z00[x][y-1] + Exn1z01[x][y+1] - 2.0 * Exn1z01[x][y] + Exn1z01[x][y-1]);
            }
        }
    }

    /****************************** Murの2次の吸収境界条件(Ex) ******************************/

    /****************************** Murの1次の吸収境界条件(Ex) ******************************/

    for(int y = 1; y < ymax+1; y++){

        if(irank == IRANK_MIN){
            double velo_dt = (C0 / sqrt(epsilonx[0][y][0]/epsilon0) ) * dt;
            double u2xa = (velo_dt - dz) / (velo_dt + dz);

            Ex[0][y][0] = Exn1z01[0][y] + u2xa * (Ex[0][y][1] - Exn1z00[0][y]);
        }

        if(irank == IRANK_MAX){
            double velo_dt = (C0 / sqrt(epsilonx[xmax-1][y][0]/epsilon0) ) * dt;
            double u2xc = (velo_dt - dz) / (velo_dt + dz);

            Ex[xmax-1][y][0] = Exn1z01[xmax-1][y] + u2xc * (Ex[xmax-1][y][1] - Exn1z00[xmax-1][y]);
        }
    }

    for(int z = 1; z < zmax; z++){

        if(irank == IRANK_MIN){
            double velo_dt = (C0 / sqrt(epsilonx[0][0][z]/epsilon0) ) * dt;
            double u1xa = (velo_dt - dy) / (velo_dt + dy);

            Ex[0][0][z] = Exn1y01[0][z] + u1xa * (Ex[0][1][z] - Exn1y00[0][z]);
        }

        if(irank == IRANK_MAX){
            double velo_dt = (C0 / sqrt(epsilonx[xmax-1][0][z]/epsilon0) ) * dt;
            double u1xc = (velo_dt - dy) / (velo_dt + dy);

            Ex[xmax-1][0][z] = Exn1y01[xmax-1][z] + u1xc * (Ex[xmax-1][1][z] - Exn1y00[xmax-1][z]);
        }
    }


    // 辺(Murの1次の吸収境界条件) -- y平面とz平面からそれぞれ算出される値の平均値を取る
    if (irank != IRANK_MIN){
        for(int x = 0; x < xmax; x++){
            double velo_dt = (C0 / sqrt(epsilonx[x][0][0]/epsilon0) ) * dt;

            double u2xa1 = (velo_dt - dx) / (velo_dt + dx);

            Ex[x][0][0] = 0.5 * (Exn1z01[x][0] + u2xa1 * (Ex[x][0][1] - Exn1z00[x][0])
                                 + Exn1y01[x][0] + u2xa1 * (Ex[x][1][0] - Exn1y00[x][0]) );
        }
    }
    else{
        for(int x = 1; x < xmax; x++){
            double velo_dt = (C0 / sqrt(epsilonx[x][0][0]/epsilon0) ) * dt;

            double u2xa1 = (velo_dt - dx) / (velo_dt + dx);

            Ex[x][0][0] = 0.5 * (Exn1z01[x][0] + u2xa1 * (Ex[x][0][1] - Exn1z00[x][0])
                                 + Exn1y01[x][0] + u2xa1 * (Ex[x][1][0] - Exn1y00[x][0]) );
        }
    }

    /****************************** Murの1次の吸収境界条件(Ex) ******************************/

    /****************************** Murの2次の吸収境界条件(Ey) ******************************/

    for(int x = 1; x < xmax; x++){
        for(int y = 1; y < ymax; y++){
            double velo_dt = (C0 / sqrt(epsilony[x][y][0]/epsilon0) ) * dt;

            double u1by1 = (velo_dt - dz) / (velo_dt + dz);
            double u2by1 = (2.0 * dz) / (velo_dt + dz);
            double u3by1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) );
            double u4by1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) );

            Ey[x][y][0] = -Eyn2z01[x][y]
            + u1by1 * (Ey[x][y][1] + Eyn2z00[x][y])
            + u2by1 * (Eyn1z00[x][y] + Eyn1z01[x][y])
            + u3by1 * (Eyn1z00[x+1][y] - 2.0 * Eyn1z00[x][y] + Eyn1z00[x-1][y] + Eyn1z01[x+1][y] - 2.0 * Eyn1z01[x][y] + Eyn1z01[x-1][y])
            + u4by1 * (Eyn1z00[x][y+1] - 2.0 * Eyn1z00[x][y] + Eyn1z00[x][y-1] + Eyn1z01[x][y+1] - 2.0 * Eyn1z01[x][y] + Eyn1z01[x][y-1]);
        }
    }
    if(irank == IRANK_MIN){
        for(int y = 1; y < ymax; y++){
            for(int z = 1; z < zmax+1; z++){
                double velo_dt = (C0 / sqrt(epsilony[0][y][z]/epsilon0) ) * dt;

                double u1cy1 = (velo_dt - dx) / (velo_dt + dx);
                double u2cy1 = (2.0 * dx) / (velo_dt + dx);
                double u3cy1 = (dx * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) );
                double u4cy1 = (dx * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) );

                Ey[0][y][z] = -Eyn2x01[y][z]
                +u1cy1 * (Ey[1][y][z] + Eyn2x00[y][z])
                +u2cy1 * (Eyn1x00[y][z] + Eyn1x01[y][z])
                +u3cy1 * (Eyn1x00[y+1][z] - 2.0 * Eyn1x00[y][z] + Eyn1x00[y-1][z] + Eyn1x01[y+1][z] - 2.0 * Eyn1x01[y][z] + Eyn1x01[y-1][z])
                +u4cy1 * (Eyn1x00[y][z+1] - 2.0 * Eyn1x00[y][z] + Eyn1x00[y][z-1] + Eyn1x01[y][z+1] - 2.0 * Eyn1x01[y][z] + Eyn1x01[y][z-1]);
            }
        }
    }
    if(irank == IRANK_MAX){
        for(int y = 1; y < ymax; y++){
            for(int z = 1; z < zmax+1; z++){
                double velo_dt = (C0 / sqrt(epsilony[xmax][y][z]/epsilon0) ) * dt;

                double u1cy2 = (velo_dt - dx) / (velo_dt + dx);
                double u2cy2 = (2.0 * dx) / (velo_dt + dx);
                double u3cy2 = (dx * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) );
                double u4cy2 = (dx * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) );

                Ey[xmax][y][z] =  -Eyn2xm1[y][z]
                + u1cy2 * (Ey[xmax-1][y][z] + Eyn2xm0[y][z])
                + u2cy2 * (Eyn1xm0[y][z] + Eyn1xm1[y][z])
                + u3cy2 * (Eyn1xm0[y+1][z] - 2.0 * Eyn1xm0[y][z] + Eyn1xm0[y-1][z] + Eyn1xm1[y+1][z] - 2.0 * Eyn1xm1[y][z] + Eyn1xm1[y-1][z])
                + u4cy2 * (Eyn1xm0[y][z+1] - 2.0 * Eyn1xm0[y][z] + Eyn1xm0[y][z-1] + Eyn1xm1[y][z+1] - 2.0 * Eyn1xm1[y][z] + Eyn1xm1[y][z-1]);
            }
        }
    }

    /****************************** Murの2次の吸収境界条件(Ey) ******************************/

    /****************************** Murの1次の吸収境界条件(Ey) ******************************/

    for(int x = 1; x < xmax; x++){
        double velo_dt = (C0 / sqrt(epsilony[x][0][0]/epsilon0) ) * dt;
        double u2ya = (velo_dt - dz) / (velo_dt + dz);
        Ey[x][0][0] = Eyn1z01[x][0] + u2ya * (Ey[x][0][1] - Eyn1z00[x][0]);
    }

    for(int z = 1; z < zmax+1; z++){
        if(irank == IRANK_MIN){
            double velo_dt = (C0 / sqrt(epsilony[0][0][z]/epsilon0) ) * dt;
            double u3ya = (velo_dt - dx) / (velo_dt + dx);
            Ey[0][0][z] = Eyn1x01[0][z] + u3ya * (Ey[1][0][z] - Eyn1x00[0][z]);
        }
        if(irank == IRANK_MAX){
            double velo_dt = (C0 / sqrt(epsilony[xmax][0][0]/epsilon0) ) * dt;
            double u3yb = (velo_dt - dx) / (velo_dt + dx);
            Ey[xmax][0][z] = Eyn1xm1[0][z] + u3yb * (Ey[xmax-1][0][z] - Eyn1xm0[0][z]);
        }
    }

    // 辺(Murの1次の吸収境界条件) --x平面とz平面からそれぞれ算出される値の平均値を取る
    for(int y = 0; y < ymax; y++){

        if(irank == IRANK_MIN){
            double velo_dt = (C0 / sqrt(epsilonx[0][y][0]/epsilon0) ) * dt;

            double u2ya1 = (velo_dt - dz) / (velo_dt + dz);
            Ey[0][y][0] = 0.5 * (Eyn1z01[0][y] + u2ya1 * (Ey[0][y][1] - Eyn1z00[0][y])
                                 + Eyn1x01[y][0] + u2ya1 * (Ey[1][y][0] - Eyn1x00[y][0]));
        }
        if(irank == IRANK_MAX){
            double velo_dt = (C0 / sqrt(epsilony[xmax][y][0]/epsilon0) ) * dt;

            double u2yc1 = (velo_dt - dz) / (velo_dt + dz);
            Ey[xmax][y][0] = 0.5*(Eyn1z01[xmax][y] + u2yc1 * (Ey[xmax][y][1] - Eyn1z00[xmax][y])
                                  + Eyn1xm1[y][0] + u2yc1 * (Ey[xmax-1][y][0] - Eyn1xm0[y][0]));
        }
    }
    /****************************** Murの1次の吸収境界条件(Ey) ******************************/

    /****************************** Murの2次の吸収境界条件(Ez) ******************************/

    for(int x = 1; x < xmax; x++){
        for(int z = 1; z < zmax; z++){
            double velo_dt = (C0 / sqrt(epsilonz[x][0][z] / epsilon0) ) * dt;

            double u1az1 = (velo_dt - dy) / (velo_dt + dy);
            double u2az1 = (2.0 * dy) / (velo_dt + dy);
            double u3az1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy) );
            double u4az1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy) );

            Ez[x][0][z] = -Ezn2y01[x][z]
            + u1az1 * (Ez[x][1][z] + Ezn2y00[x][z])
            + u2az1 * (Ezn1y00[x][z] + Ezn1y01[x][z])
            + u3az1 * (Ezn1y00[x+1][z] - 2.0 * Ezn1y00[x][z] + Ezn1y00[x-1][z] + Ezn1y01[x+1][z] - 2.0 * Ezn1y01[x][z] + Ezn1y01[x-1][z])
            + u4az1 * (Ezn1y00[x][z+1] - 2.0 * Ezn1y00[x][z] + Ezn1y00[x][z-1] + Ezn1y01[x][z+1] - 2.0 * Ezn1y01[x][z] + Ezn1y01[x][z-1]);
        }
    }

    for(int y = 1; y < ymax+1; y++){
        for(int z = 1; z < zmax; z++){
            if(irank == IRANK_MIN){
                double velo_dt = (C0 / sqrt(epsilonz[0][y][z] / epsilon0) ) * dt;

                double u1cz1 = (velo_dt - dx) / (velo_dt + dx);
                double u2cz1 = (2.0 * dx) / (velo_dt + dx);
                double u3cz1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) );
                double u4cz1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) );

                Ez[0][y][z] = -Ezn2x01[y][z]
                + u1cz1 * (Ezn2x00[y][z] + Ez[1][y][z])
                + u2cz1 * (Ezn1x00[y][z] + Ezn1x01[y][z])
                + u3cz1 * (Ezn1x00[y+1][z] - 2.0 * Ezn1x00[y][z] + Ezn1x00[y-1][z] + Ezn1x01[y+1][z] - 2.0 * Ezn1x01[y][z] + Ezn1x01[y-1][z])
                + u4cz1 * (Ezn1x00[y][z+1] - 2.0 * Ezn1x00[y][z] + Ezn1x00[y][z-1] + Ezn1x01[y][z+1] - 2.0 * Ezn1x01[y][z] + Ezn1x01[y][z-1]);
            }
            if(irank == IRANK_MAX){
                double velo_dt = (C0 / sqrt(epsilonz[xmax][y][z] / epsilon0) ) * dt;

                double u1cz2 = (velo_dt - dx) / (velo_dt + dx);
                double u2cz2 = (2.0 * dx) / (velo_dt + dx);
                double u3cz2 = (dy * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) );
                double u4cz2 = (dz * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) );

                Ez[xmax][y][z] = -Ezn2xm1[y][z]
                + u1cz2 * (Ezn2xm0[y][z] + Ez[xmax-1][y][z])
                + u2cz2 * (Ezn1xm1[y][z] + Ezn1xm0[y][z])
                + u3cz2 * (Ezn1xm1[y+1][z] - 2.0 * Ezn1xm1[y][z] + Ezn1xm1[y-1][z] + Ezn1xm0[y+1][z] - 2.0 * Ezn1xm0[y][z] + Ezn1xm0[y-1][z])
                + u4cz2 * (Ezn1xm1[y][z+1] - 2.0 * Ezn1xm1[y][z] + Ezn1xm1[y][z-1] + Ezn1xm0[y][z+1] - 2.0 * Ezn1xm0[y][z] + Ezn1xm0[y][z-1]);
            }
        }
    }

    /****************************** Murの2次の吸収境界条件(Ez) ******************************/

    /****************************** Murの1次の吸収境界条件(Ez) ******************************/
    for(int x = 1; x < xmax; x++){
        double velo_dt = (C0 / sqrt(epsilonz[x][0][0] / epsilon0) ) * dt;
        double u1za = (velo_dt - dy) / (velo_dt + dy);

        Ez[x][0][0] = Ezn1y01[x][0] + u1za * (Ez[x][1][0] - Ezn1y00[x][0]);
    }

    for(int y = 1; y < ymax+1; y++){

        if(irank == IRANK_MIN){
            double velo_dt = (C0 / sqrt(epsilonz[0][y][0] / epsilon0) ) * dt;
            double u3za = (velo_dt - dx) / (velo_dt + dx);

            Ez[0][y][0] = Ezn1x01[y][0] + u3za * (Ez[1][y][0] - Ezn1x00[y][0]);
        }
        if(irank == IRANK_MAX){
            double velo_dt = (C0 / sqrt(epsilonz[xmax][y][0] / epsilon0) ) * dt;
            double u3zb = (velo_dt - dx) / (velo_dt + dx);

            Ez[xmax][y][0] = Ezn1xm1[y][0] + u3zb * (Ez[xmax-1][y][0] - Ezn1xm0[y][0]);
        }
    }

    // 辺(Murの1次の吸収境界条件) --x平面とy平面からそれぞれ算出される値の平均値を取る
    for(int z = 0; z < zmax+1; z++){

        if(irank == IRANK_MIN){
            double velo_dt = (C0 / sqrt(epsilonz[0][0][z] / epsilon0) ) * dt;
            double u1za1 = (velo_dt - dy) / (velo_dt + dy);

            Ez[0][0][z] = 0.5 * (Ezn1y01[0][z] + u1za1 * (Ez[0][1][z] - Ezn1y00[0][z])
                                 + Ezn1x01[0][z] + u1za1 * (Ez[1][0][z] - Ezn1x00[0][z]) );
        }
        if(irank == IRANK_MAX){
            double velo_dt = (C0 / sqrt(epsilonz[xmax][0][z] / epsilon0) ) * dt;
            double u1zb1 = (velo_dt - dy) / (velo_dt + dy);

            Ez[xmax][0][z] = 0.5 * (Ezn1y01[xmax][z] + u1zb1 * (Ez[xmax][1][z] - Ezn1y00[xmax][z])
                                    + Ezn1xm1[0][z] + u1zb1 * (Ez[xmax-1][0][z] - Ezn1xm0[0][z]));
        }
    }
    /****************************** Murの1次の吸収境界条件(Ez) ******************************/

}

/*電界の保存*/
void saving_electric_field()
{
    // Ex
    for(int x = 0; x < xmax+1; x++){
        for(int z = 0; z < zmax+1; z++){
            Exn2y00[x][z] = Exn1y00[x][z];
            Exn1y00[x][z] = Ex[x][0][z];
            Exn2y01[x][z] = Exn1y01[x][z];
            Exn1y01[x][z] = Ex[x][1][z];
        }
    }
    for(int x = 0; x < xmax+1; x++){
        for(int y = 0; y < ymax+1; y++){
            Exn2z00[x][y] = Exn1z00[x][y];
            Exn1z00[x][y] = Ex[x][y][0];
            Exn2z01[x][y] = Exn1z01[x][y];
            Exn1z01[x][y] = Ex[x][y][1];
        }
    }

    // Ey
    for(int x = 0; x < xmax+1; x++){
        for(int y = 0; y < ymax; y++){
            Eyn2z00[x][y] = Eyn1z00[x][y];
            Eyn1z00[x][y] = Ey[x][y][0];
            Eyn2z01[x][y] = Eyn1z01[x][y];
            Eyn1z01[x][y] = Ey[x][y][1];
        }
    }
    for(int y = 0; y < ymax; y++){
        for(int z = 0; z < zmax+1; z++){
            if(irank == IRANK_MIN){
                Eyn2x00[y][z] = Eyn1x00[y][z];
                Eyn1x00[y][z] = Ey[0][y][z];
                Eyn2x01[y][z] = Eyn1x01[y][z];
                Eyn1x01[y][z] = Ey[1][y][z];
            }
            if(irank == IRANK_MAX){
                Eyn2xm1[y][z] = Eyn1xm1[y][z];
                Eyn1xm1[y][z] = Ey[xmax-1][y][z];
                Eyn2xm0[y][z] = Eyn1xm0[y][z];
                Eyn1xm0[y][z] = Ey[xmax][y][z];
            }
        }
    }

    //Ez
    for(int x = 0; x < xmax+1; x++){
        for(int z = 0; z < zmax; z++){
            Ezn2y00[x][z] = Ezn1y00[x][z];
            Ezn1y00[x][z] = Ez[x][0][z];
            Ezn2y01[x][z] = Ezn1y01[x][z];
            Ezn1y01[x][z] = Ez[x][1][z];
        }
    }
    for(int y = 0; y < ymax+1; y++){
        for(int z = 0; z < zmax; z++){
            if(irank == IRANK_MIN){
                Ezn2x00[y][z] = Ezn1x00[y][z];
                Ezn1x00[y][z] = Ez[0][y][z];
                Ezn2x01[y][z] = Ezn1x01[y][z];
                Ezn1x01[y][z] = Ez[1][y][z];
            }
            if(irank == IRANK_MAX){
                Ezn2xm1[y][z] = Ezn1xm1[y][z];
                Ezn1xm1[y][z] = Ez[xmax-1][y][z];
                Ezn2xm0[y][z] = Ezn1xm0[y][z];
                Ezn1xm0[y][z] = Ez[xmax][y][z];
            }
        }
    }
}

void calc_poynting_powerHz() // 透過スペクトル計算用の観測面中央でのHzの出力
{
    /*
    double pmax_01 = 0;		// 入力パワーの最大値を記録する変数
    double pmax_03 = 0;		// 出力パワーの最大値を記録
    double pmin_01 = 0;		// 出力パワーの最小値を記録
    double pmin_03 = 0;		// 出力パワーの最小値を記録
    int x = 0;
    if(irank == intObseInPortNum){ //入射
        x = intObseLenPart1;
        y = ymax - 1;
        z = zmax - 1;
        fprintf(fpHz1, "%e¥n", Hz[x][y][z]);
    }

    if(irank == intObseOutPortNum){			// 出射
        x = intObseLenPart4;
        y = ymax - 1;
        z = zmax - 1;
        fprintf(fpHz5, "%e¥n", Hz[x][y][z]);
    }
    */
}


void output_field_write(char *dir_name_def){
    FILE *EXY2;
    char  fname3[40];
    int x, y;
    int z = 20;
    sprintf(fname3, "/Field_E_XY_%d_%d.txt", n, irank);
    EXY2 = fopen(strcat(strcpy(dir_name[0], dir_name_def), fname3), "w");
    double E_xy;
    for(int x = 0; x < xmax; x++){
      for(int y = 0; y < ymax; y++){
        E_xy = Hz[x][y][z];
        fprintf(EXY2, "%e¥t", E_xy);
      }
      fprintf(EXY2, "¥n");
    }
    fclose(EXY2);
}


void output_field(char *dir_name_def){
  if(n == Nmax){
			output_field_write (dir_name_def);
	}
}

void output_model(){
  int x, y, z;

  if(irank == IRANK_MIN){
  for(x = 0; x < xmax; x++){
    for(y = 0; y < ymax; y++){
      fprintf(model_xy, "%d¥t", cell[x][y][20]);
    }
    fprintf(model_xy, "¥n");
  }
  fclose(model_xy);
  for(x = 0; x < xmax; x++){
    for(y = 0; z < zmax; y++){
      fprintf(model_xz, "%d¥t", cell[x][ymax/2][z]);
    }
    fprintf(model_xz, "¥n");
  }
  fclose(model_xz);
  for(x = 0; z < zmax; x++){
    for(y = 0; y < ymax; y++){
      fprintf(model_yz, "%d¥t", cell[xmax/2][y][z]);
    }
    fprintf(model_yz, "¥n");
  }
  fclose(model_yz);
}
}



void main_calc(int wave, int left, int right)
{
    double s_time, e_time;
    char time[40] = {0};

    MPI_Status status;
    int tag_send = 0, tag_recv = 0;
    int dir_count = 0;
    initialize_matrix(); 						// 配列の初期化
    modeling(); 								// モデルの設定
    file_open(dir_name[dir_count]); 			// ファイルを開く
    //parameter(dir_name[dir_count]); 			// パラメータの設定と出力

    // 計算開始時刻の出力
    if (irank == IRANK_MIN){
        _strtime(time);
        //fprintf(fpparameter, "Start Time:¥t %s¥n", time);
        s_time = MPI_Wtime();
    }

    // 電磁界計算
    //    observation_func(); 	// 観測点の設定
        output_model(); 		// モデルの出力
    //    set_epsilon(); 			// 誘電率の割り当て

    for(int n = 1 ; n <= Nmax; n++){

        // 時間ステップ数の表示
        if(n % Ncut == 0){
            _strtime(time);
            printf("n = %d, time = %s\n", n, time);
        }

        // 励振関数の設定
        source_func(n);

        // 一度同期をとる(同期はノード間で速度にばらつきが生じる作業)
        MPI_Barrier (MPI_COMM_WORLD);

        // 電界の計算
        calc_efield();

        // 吸収境界条件による端面の計算
        absorpt_bound_condition();

        // 一度同期をとる
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Sendrecv( &Ex[1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_send,
                     &Ex[xmax][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
        MPI_Sendrecv( &Ey[1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_send,
                     &Ey[xmax][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
        MPI_Sendrecv( &Ez[1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_send,
                     &Ez[xmax][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
        MPI_Sendrecv( &Ex[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
                     &Ex[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
        MPI_Sendrecv( &Ey[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
                     &Ey[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
        MPI_Sendrecv( &Ez[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
                     &Ez[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);

        // 電界の保存
        saving_electric_field();

        // 磁界の計算
        calc_hfield();

        MPI_Sendrecv( &Hy[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
                     &Hy[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
        MPI_Sendrecv( &Hz[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
                     &Hz[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);

        // フィールドの出力
        output_field (dir_name[dir_count]);

        // ポインティングパワー計算と出力
        //calc_poynting_power();
        calc_poynting_powerHz();

        // 一度同期をとる
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (irank == IRANK_MIN){
        _strtime(time);
        //fprintf(fpparameter, "End Time:¥t %s¥n", time);
        //時刻の出力
        e_time = MPI_Wtime();
        printf ("¥ntime = %f¥n", e_time - s_time);
    }
    file_close(); 			// ファイルを閉じる
}

int main(int argc, char **argv)
{
    char processor_name[MPI_MAX_PROCESSOR_NAME] = {0};
    int namelen, isize;
    // MPIによる通信の開始
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &irank);
    MPI_Comm_size (MPI_COMM_WORLD, &isize);
    IRANK_MAX = isize - 1;
    MPI_Get_processor_name (processor_name, &namelen);

    if(irank == 0) {
        printf ("Parallel processing on %d nodes\n", isize);
    }
    printf ("Process %d on %s\n", irank, processor_name);

    // 隣の計算機の番号の指定
    int left = (irank == IRANK_MIN) ? MPI_PROC_NULL : irank - 1;
    int right = (irank == IRANK_MAX) ? MPI_PROC_NULL : irank + 1;
    main_calc(1550, left, right);
    MPI_Finalize();
}
