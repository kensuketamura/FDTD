/*****************************************************************************/
// 一つの計算機での解析領域
/*****************************************************************************/
#define XMAX (INT_DIV(XMAX_ALL, NODE) + 1)
#define YMAX YMAX_ALL
#define ZMAX ZMAX_ALL


/*****************************************************************************/
// 励振関数
/*****************************************************************************/
//char *dir_name[] = {"1550", "1580"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1555", "1560", "1570", "1575"}; 		// 励振関数の波長 [nm]
double lambda; 											// 励振関数の波長 [m] (プログラム中で dir_name をdouble化して代入)
double omega0; 											// 励振関数の角周波数 [s^-1]

// ガウシアンパルス
double sigma; 											// 広がり幅を決定する定数
//static const double delta_omega = 0.1; 						// 中心周波数で規格した半値全幅
//static const int Npeak = 3000; 								// ピークステップ数


/*****************************************************************************/
// Model掃き出し用の定数
/*****************************************************************************/
#define CLAD 0									// クラッド
#define CORE 1									// コア
#define GaInAsP	2								// GaInAsP
#define AIR_GaInAsP	3							// CLAD/GaInAsP
#define EXITATION 30							// 励振点//モデル上だと31
#define OBSERVATION 20							// 観測点//モデル上だと21
#define CIRCLE_REF_INDEX	CLAD				//関数mcircleで書き込む数字
#define CIRCLE_REF_INDEX2	CLAD		//関数mcircle2で書き込む数字
#define CIRCLE_REF_INDEX3	2				//場所確認用ドットの色指定


/*****************************************************************************/
// フォトニック結晶モデル用宣言
/*****************************************************************************/
//厚さ方向パラメータ
struct PNUM {int X; int Y; }; 					//円柱の中心座標を与える変数.sankakuで使用


/*****************************************************************************/
// グローバル変数 (MPIの通信ではグローバル変数を使用しないとエラーが生じる模様)
/*****************************************************************************/

#if _CALCULATION_TYPE == _PROPAGATION_CALCULATION

// 電磁界 (XMAX+1は並列計算での"のりしろ"が必要なため，YMAX+1, ZMAX+1は対称境界条件で必要なため)
double Ex[XMAX+1][YMAX+1][ZMAX+1]; 
double Ey[XMAX+1][YMAX+1][ZMAX+1]; 
double Ez[XMAX+1][YMAX+1][ZMAX+1]; 
double Hx[XMAX+1][YMAX+1][ZMAX+1]; 
double Hy[XMAX+1][YMAX+1][ZMAX+1]; 
double Hz[XMAX+1][YMAX+1][ZMAX+1]; 

// 誘電率 (YMAX+1, ZMAX+1は対称境界条件で必要なため)
double ALL_epsilonx[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1]; 
double ALL_epsilony[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1]; 
double ALL_epsilonz[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1]; 
double epsilonx[XMAX][YMAX+1][ZMAX+1]; 
double epsilony[XMAX][YMAX+1][ZMAX+1]; 
double epsilonz[XMAX][YMAX+1][ZMAX+1]; 
int ALL_cell[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];		// 実際はALL_cell[XMAX_ALL+1][YMAX_ALL][ZMAX_ALL]で良いが，プログラムの都合上誘電体の配列と同じ次元にした
int cell[XMAX][YMAX+1][ZMAX+1];							// 実際はcell[XMAX][YMAX][ZMAX]で良いが，プログラムの都合上誘電体の配列と同じ次元にした
double epsilon_xy[XMAX+1][YMAX+1], epsilon_yz[YMAX+1][ZMAX+1], epsilon_zx[XMAX+1][ZMAX+1]; 
double epsilon_zx2[XMAX+1][ZMAX+1]; 
int cell_xy[XMAX][YMAX], cell_yz[YMAX][ZMAX]; 
int cell_zx[ZMAX][XMAX];//☆16/1/6


/*吸収境界条件適用のときに使う電界の配列*/
double Exn2y00[XMAX+1][ZMAX+1], 	Exn1y00[XMAX+1][ZMAX+1+1], 	Exn2y01[XMAX+1][ZMAX+1], 	Exn1y01[XMAX+1][ZMAX+1+1]; 
//double exn2ym1[XMAX+1][ZMAX+1], 	exn1ym1[XMAX+1][ZMAX+1+1], 	exn2ym0[XMAX+1][ZMAX+1], 	exn1ym0[XMAX+1][ZMAX+1+1]; 
double Exn2z00[XMAX+1][YMAX+1], 	Exn1z00[XMAX+1][YMAX+1+1], 	Exn2z01[XMAX+1][YMAX+1], 	Exn1z01[XMAX+1][YMAX+1+1]; 
//double exn2zm1[XMAX+1][YMAX+1], 	exn1zm1[XMAX+1][YMAX+1+1], 	exn2zm0[XMAX+1][YMAX+1], 	exn1zm0[XMAX+1][YMAX+1+1]; 

double Eyn2z00[XMAX+1][YMAX], 	Eyn1z00[XMAX+1][YMAX+1], 	Eyn2z01[XMAX+1][YMAX], 	Eyn1z01[XMAX+1][YMAX+1]; 
//double eyn2zm1[XMAX+1][YMAX], 	eyn1zm1[XMAX+1][YMAX+1], 	eyn2zm0[XMAX+1][YMAX], 	eyn1zm0[XMAX+1][YMAX+1]; 
double Eyn2x00[YMAX][ZMAX+1], 	Eyn1x00[YMAX+1][ZMAX+1+1], 	Eyn2x01[YMAX][ZMAX+1], 	Eyn1x01[YMAX+1][ZMAX+1+1]; 
double Eyn2xm1[YMAX][ZMAX+1], 	Eyn1xm1[YMAX+1][ZMAX+1+1], 	Eyn2xm0[YMAX][ZMAX+1], 	Eyn1xm0[YMAX+1][ZMAX+1+1]; 

double Ezn2y00[XMAX+1][ZMAX], 	Ezn1y00[XMAX+1][ZMAX+1], 	Ezn2y01[XMAX+1][ZMAX], 	Ezn1y01[XMAX+1][ZMAX+1]; 
//double ezn2ym1[XMAX+1][ZMAX], 	ezn1ym1[XMAX+1][ZMAX+1], 	ezn2ym0[XMAX+1][ZMAX], 	ezn1ym0[XMAX+1][ZMAX+1]; 
double Ezn2x00[YMAX+1][ZMAX], 	Ezn1x00[YMAX+1+1][ZMAX+1], 	Ezn2x01[YMAX+1][ZMAX], 	Ezn1x01[YMAX+1+1][ZMAX+1]; 
double Ezn2xm1[YMAX+1][ZMAX], 	Ezn1xm1[YMAX+1+1][ZMAX+1], 	Ezn2xm0[YMAX+1][ZMAX], 	Ezn1xm0[YMAX+1+1][ZMAX+1]; 

// 磁界分布出力用
double field_xy[XMAX][YMAX]; 	// Hz-field のファイル出力 (面垂直方向の磁界成分)
double field_yz[YMAX][ZMAX];	// Hz-field のファイル出力 (面方向の磁界成分)
double field_zx_Hz[ZMAX][XMAX];	// Hy-field のファイル出力 (面方向の磁界成分)
double field_zx_Hy[ZMAX][XMAX];	// Hy-field のファイル出力 (面方向の磁界成分)

int x, y, z; 				// 離散座標
int xmax, ymax, zmax; 			// 空間分割の最大値
int xmax_all, ymax_all, zmax_all; //分割前の最大値
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
//char fname[30]; 
int irank, isize; 

// 電磁界計算に使う定数の宣言
//double dex, dey, dez; 				//電界の計算に使う差分値	dex: Ex, dey: Ey, dez: Ez
//double cnstEx, cnstEy, cnstEz; 	//電界の計算に使う定数		
//double dhx, dhy, dhz; 				//磁界の計算に使う差分値	dhx: Hx, dhy: Hy, dhz: Hz
static const double cnstHxyz = dt / MU0; 					//磁界の計算に使う定数		cnstHxyz: Magnetic field calculation


///*励振用定数の宣言*/
//double omega0; 	//omega0: 励振関数の角周波数
//double sigma; 	//ガウシアンパルスのための広がり幅を決定する定数
//int Npeak; 		//ガウシアンパルスのピークステップ数

/*モデル設定に使う定数の宣言 -- 関数modeling()以外でも使う定数があるのでここで宣言する．*/
double b; //ディスク厚さ
int bk; //ディスク厚さの離散値

double w; 		//支柱幅の半分
int wij; 		//支柱幅の半分の離散値
double h; 		//支柱の高さ
int hk; 			//支柱の高さの離散値
//double np; 		//支柱の屈折率
/*
#define CLAD 0
#define AlGaAs3 1
#define AlGaAs8 2
#define GaAs 3
#define AlAs 4
#define AlOxAs 5

double htslab; 		//スラブの厚さの半分
int ihtslab; 		//スラブの厚さの半分の離散値
double rhole; 		//円孔の半径
int irhole; 			//円孔半径の離散値
double hphole; 		//円孔配列ピッチの半分
int ihphole; 		//円孔配列ピッチの半分の離散値
double hpholev; 		//円孔配列ピッチの半分の√3倍
int ihpholev; 		//円孔配列ピッチの半分の√3倍の離散値
*/
//入出力パワーの最大値と最小値を記録する変数
double powermax_in; 
double powermin_in; 
double powermax_out; 
double powermin_out; 

//double n_core; 		//活性層屈折率
//double n_clad; 		//クラッド層屈折率
double topy1, topy5, tohz1, tohz5; 
double topy1h, topy5h, tohz1h, tohz5h; 
//double epsilon1; 
//double epsilon2; 

void mcircle(int, int , int, int);  //make circle function
void halfcircle(int, int , int, int); 							//make circle function
void rightquartercircle1(int, int, int, int, double); 
void leftquartercircle1(int, int, int, int, double); 
void rightquartercircle2(int, int, int, int, double); 
void leftquartercircle2(int, int, int, int, double); 

/////////////////////////////////////////////////////
/*観測点*/
int inputi, inputj, inputk; 
int outputi, outputj, outputk; 
double pyn_in, pyn_out; 

/*WGチャープ構造*/
#define WG_chirp_gradient 1
#define WG_chirp_gradient2 0.5
#define WG_chirp_gradient3 0.08 //BOUNDARYLINE == 9で使用
#define WG_chirp_gradient4 0.04 //BOUNDARYLINE == 10で使用
#define WG_chirp_gradient5 0.015 //BOUNDARYLINE == 14で使用
#define WG_chirp_gradient6 0.020 //BOUNDARYLINE == 15で使用
#define WG_chirp_gradient7 0.030 //BOUNDARYLINE == 16で使用
#define WG_chirp_gradient8 0.020 //BOUNDARYLINE == 17で使用
#define WG_chirp_gradient9 0.012 //BOUNDARYLINE == 18で使用
#define WG_chirp_gradient10 0.010 //BOUNDARYLINE == 19で使用
#define WG_chirp_gradient11 0.008 //BOUNDARYLINE == 20で使用
#define WG_chirp_gradient12 0.012 //BOUNDARYLINE == 22で使用
#define WG_chirp_gradient13 0.016 //BOUNDARYLINE == 23で使用
#define WG_chirp_gradient14 0.020 //BOUNDARYLINE == 24で使用
#define WG_chirp_gradient15 0.024 //BOUNDARYLINE == 25で使用 1/7
#define WG_chirp_gradient16 0.004 //BOUNDARYLINE == 26で使用 1/7
#define WG_chirp_gradient17 0.002 //BOUNDARYLINE == 27で使用 1/7


#define WG_chirp_off_in_x  187
#define WG_chirp_off_in_x2  187-16
#define WG_chirp_off_in_x21  291 // 1/4
#define WG_chirp_off_in_x33  316

#define WG_chirp_off_in_x_2  0 //BOUNDARYLINE == 9で使用


#define WG_chirp_off_out_x 800+1 //これでテーパが行き過ぎないように制限
#define WG_chirp_off_out_x2 800+1+16
#define WG_chirp_off_out_x9 788+200 //BOUNDARYLINE == 14で使用
#define WG_chirp_off_out_x18 788+200+42 //BOUNDARYLINE == 18で使用
#define WG_chirp_off_out_x19 788+200+63 //BOUNDARYLINE == 19で使用
#define WG_chirp_off_out_x20 788+200+84 //BOUNDARYLINE == 20で使用
#define WG_chirp_off_out_x21 788+200+84 //BOUNDARYLINE == 21で使用

#define WG_chirp_off_out_y 75.5+3+40+16
#define WG_chirp_off_out_y9 75.5+3+40+16+1  //BOUNDARYLINE == 9で使用
#define WG_chirp_off_out_y15 75.5+3+40+16+1  //BOUNDARYLINE == 15で使用
#define WG_chirp_off_out_y16 75.5+3+40+16  //BOUNDARYLINE == 16で使用
#define WG_chirp_off_out_y17 75.5+3+40+16+1  //BOUNDARYLINE == 17で使用
#define WG_chirp_off_out_y18 75.5+3+40+16+1 //BOUNDARYLINE == 18で使用
#define WG_chirp_off_out_y19 75.5+3+40+16+1  //BOUNDARYLINE == 19で使用
#define WG_chirp_off_out_y20 75.5+3+40+16+1-10  //BOUNDARYLINE == 20で使用
#define WG_chirp_off_out_y21 75.5+3+40+16+1  //BOUNDARYLINE == 21で使用



/*ファイルポインタ宣言*/
//モデル			//Ez				//Hz
FILE *model_xy; 		FILE *fpfez_xy; 		FILE *fpfhz_xy; 
FILE *model_yz; 		FILE *fpfez_yz; 		FILE *fpfhz_yz; 
FILE *model_xz; 		FILE *fpfez_zx; 		FILE *fpfhz_zx;//☆ 
FILE *model_xy2; 	FILE *fpfez2_xy; 	FILE *fpfhz2_xy; 
FILE *allmodel_xy; 	 
FILE *allmodel_yz1, *allmodel_yz4, *allmodel_yz7;
FILE *allmodel_zx1, *allmodel_zx4, *allmodel_zx7; //☆16/1/6 
FILE *allmodel_zx; //☆16/1/6 

//誘電率
FILE *fpepsilonx; 
FILE *fpallepsilonx; 
FILE *fpepsilony; 
FILE *fpepsilonz; 
FILE *fpepsilony2; 
FILE *fpepsilonz2; 
FILE *fpAllEpsilon; 
FILE *fpEpsilon; 


FILE *fpex; 		FILE *fphx; 
FILE *fpey; 		FILE *fphy; 
FILE *fpez; 		FILE *fphz; 
/*	FILE *fpez1; 	FILE *fphz1; 
FILE *fpez2; 	FILE *fphz2; 
FILE *fpez3; 	FILE *fphz3; 
FILE *fpez4; 	FILE *fphz4; 
FILE *fpez5; 	FILE *fphz5; 
FILE *fpez6; 	FILE *fphz6; 
*/

//	FILE *fpfhr_xy; 	FILE *fpfhth_xy; 	FILE *fpfer_xy; 		FILE *fpfeth_xy; 	FILE *fpfethhz_xy; 
//	FILE *fpfhr_yz; 	FILE *fpfhth_yz; 	FILE *fpfer_yz; 		FILE *fpfeth_yz; 	FILE *fpfethhz_yz; 
//	FILE *fpfhr_zx; 	FILE *fpfhth_zx; 	FILE *fpfer_zx; 		FILE *fpfeth_zx; 	FILE *fpfethhz_zx; 
//	FILE *fpfhr2_xy; FILE *fpfhth2_xy; 	FILE *fpfer2_xy; 	FILE *fpfeth2_xy; 	FILE *fpfethhz2_xy; 

//	FILE *fphr1; 	FILE *fphth1; 		FILE *fper1; 		FILE *fpeth1; 
//	FILE *fphr2; 	FILE *fphth2; 		FILE *fper2; 		FILE *fpeth2; 
//	FILE *fphr3; 	FILE *fphth3; 		FILE *fper3; 		FILE *fpeth3; 
//	FILE *fphr4; 	FILE *fphth4; 		FILE *fper4; 		FILE *fpeth4; 
//	FILE *fphr5; 	FILE *fphth5; 		FILE *fper5; 		FILE *fpeth5; 
//	FILE *fphr6; 	FILE *fphth6; 		FILE *fper6; 		FILE *fpeth6; 

FILE *fpparameter; 	//計算パラメータ保存ファイル
FILE *fppoynt_para; 		//ポインティングパワーEx*Hy, Ey*Hx各成分保存ファイル
FILE *fppoynt1; 		//ポインティングパワー保存ファイル
FILE *fppoynt2; 		//ポインティングパワー保存ファイル
FILE *fppoynt3; 		//ポインティングパワー保存ファイル
FILE *fppoynt4; 		//ポインティングパワー保存ファイル
FILE *fppoynt5; 		//ポインティングパワー保存ファイル
FILE *fppoynt6; 		//ポインティングパワー保存ファイル
FILE *fppoynt1h; 		//ポインティングパワー保存ファイル
FILE *fppoynt2h; 		//ポインティングパワー保存ファイル
FILE *fppoynt3h; 		//ポインティングパワー保存ファイル
FILE *fppoynt4h; 		//ポインティングパワー保存ファイル
FILE *fppoynt5h; 		//ポインティングパワー保存ファイル
FILE *fppoynt6h; 		//ポインティングパワー保存ファイル
FILE *fppowerHz1; 		//断面全体Hz保存用
FILE *fppowerHz2; 		//断面全体Hz保存用
FILE *fppowerHz3; 		//断面全体Hz保存用
FILE *fppowerHz4; 	//断面全体Hz保存用
FILE *fppowerHz5; 	//断面全体Hz保存用
FILE *fppowerHz6; 	//断面全体Hz保存用
FILE *fppowerHz1h; 	//断面全体Hz保存用
FILE *fppowerHz2h; 	//断面全体Hz保存用
FILE *fppowerHz3h; 	//断面全体Hz保存用
FILE *fppowerHz4h; 	//断面全体Hz保存用
FILE *fppowerHz5h; 	//断面全体Hz保存用
FILE *fppowerHz6h; 	//断面全体Hz保存用
FILE *fpHz1; 		// 入射観測面Hz保存ファイル
FILE *fpHz5; 		// 出射観測面Hz保存ファイル
FILE *avpoynt1; 	
FILE *avpoynt5; 
FILE *avhz1; 	
FILE *avhz5; 	
FILE *avpoynt1h; 	
FILE *avpoynt5h; 
FILE *avhz1h; 	
FILE *avhz5h; 	
//	FILE *fpenergy; 		FILE *fppower; 
//	FILE *fpenergye; 	FILE *fppowerv; 
//	FILE *fpenergyh; 	FILE *fppowerp; 

#elif _CALCULATION_TYPE == _BAND_CALCULATION

#endif
