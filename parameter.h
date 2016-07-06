/*****************************************************************************/
// 並列計算時の計算機の台数
/*****************************************************************************/
#if _FDTD
#define NODE 1
#else
#define NODE 1
#endif

// 並列計算時の計算機のランク
#define IRANK_MAX (NODE-1)		// 最大値
#define IRANK_MIN 0				// 最小値(この数字の割り当てられた計算機に結果が入る)
#define ISIZE NODE				// サイズ

#define TRUE 1
#define FALSE 0
/*****************************************************************************/
// 計算パラメータ		単位はnm
// (配列の定義やプリプロセッサで使用するため，マクロ化)
// 各値はセルサイズの整数倍になっていること!!
/*****************************************************************************/


/****************************** 3列目格子シフト構造 ******************************/
#define CELL_SIZE 21			// セルサイズ
#define PITCH 399				// PC 格子定数
#define PITCH_SHIFT_MAX 399 //480		// 格子定数変化PCWのPC格子定数の最大値

#define SLAB_HEIGHT 210		// スラブ厚
#define CLAD_HEIGHT1 525	// 上部クラッド高さ
#define CLAD_HEIGHT2 0		// 下部クラッド高さ
#define AIR_HEIGHT 0		// 空気層高さ

#define RADIUS 120			// PCの標準円孔半径
//#define SX3 80				// 伝搬(X)方向の3列目格子シフト量(SX2,SX4=0でないと使えない要改善!!mainの1380行目)
//#define SX1 0				// 伝搬(X)方向の1列目格子シフト量
//#define SX2 0				// 伝搬(X)方向の2列目格子シフト量(SX3,SX4=0でないと使えない要改善!!mainの1380行目)
//#define SX4 0				//伝搬(X)方向の4列目格子シフト量(SX2,SX3=0でないと使えない要改善!!mainの1380行目)
//#define SY 0				// 横(Y)方向の導波路全体格子シフト量(対称境界を使用しているので実際にはこの倍)

#define EXCT_LEN 840					// 励振点 (モデルの左端からの距離)
#define EXCT_OBSE_LEN 975				// 励振点から観測面の中心までの距離
#define OBSE_WIRE_LEN 2540				// 観測面の中心から細線導波路端までの距離
#define OBSE_INTER (PITCH)				// 観測面の長さ(透過率，反射率を求めるために使用)
#define WIRE_OUTPUT_LEN (EXCT_LEN + 45)	// 出射細線導波路の長さ(直打ちの数字はセルサイズをNODE数の整数倍にするために使用)
#define WIRE_OUTPUT_OFFSET 0			// 出射細線導波路のスラブ終端の長さ
#define WIRE_WID_OFFSET 180				// 細線幅のオフセット量の半分(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_SiSLAB_TERMINATION_LEN 255	// PCW横のCOREスラブ終端の長さ
#define PCW_SiSLAB_OFFSET 0				// PCW縦のCOREスラブのオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_WIDTH_CHIRP 180				// PCW幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)

// これ以下は周期数	
//#define NORM_PCW_PER 0				// 通常PCW周期数
//#define CHIRP_3RD_LS_PER 0			// 3列目格子シフト量チャープLSPCW周期数． 2~4列目まで対応．※シフト量/周期数<cellsizeだとチャープしない
//#define CHIRP_2ND_LS_PER 0			//導波路幅チャープとシフト量チャープを同時に行う際のチャープLSPCW周期数．PITCH_SHIFT_PERより小さくないとダメ?　2,3列目のみ対応．※シフト量/周期数<cellsizeだとチャープしない
//#define PITCH_SHIFT_PER 5				// 格子定数変化PCWの周期数//修正=>導波路幅チャープの周期数
//#define LSPCW_SHIFT_DESCRETE FALSE			// 格子定数変化PCWのとき，はじめからLSPCWにする場合はFALSE　しない場合はTRUE
//#define PITCH_SHIFT_CHIRP_PER 0		// 格子定数変化チャープPCWの周期数
//#define LSPCW_PER 20				// LSPCW周期数
//#define PCW_WID 6				// PCWの列数
// これ以上は周期数

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// 通常PCW長
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// チャープLSPCW長
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW長
#define WIRE_OUTPUT_OFFSET_PER (WIRE_OUTPUT_OFFSET, CELL_SIZE)	// 出射細線導波路のスラブ終端の長さ

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// 入射観測面の中心座標 (モデルの左端からの距離)
#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)	// 出射観測面の中心座標 (モデルの右端からの距離)

#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)			// 入射細線長 
#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)		// 出射細線長

/****************************** 3列目格子シフト構造 ******************************/

/*****************************************************************************/
// 全解析領域
/*****************************************************************************/
#define XMAX_ALL 200	// SiO2 1340
#define YMAX_ALL 163	//163 PCW_WID:6/163  PCW_WID:8/209  PCW_WID:10/255
#define ZMAX_ALL 42		// CLAD_HEIGHT1:525/42  CLAD_HEIGHT1:750/57  CLAD_HEIGHT1:990/73


/*****************************************************************************/
// セルサイズ [m]
/*****************************************************************************/
static const double dblCellSize = CELL_SIZE * 1e-9; 
static const double dx = dblCellSize; 
static const double dy = dblCellSize; 
static const double dz = dblCellSize; 
static const double inv_dx = 1/dx; 
static const double inv_dy = 1/dy; 
static const double inv_dz = 1/dz; 

/*****************************************************************************/
// 時間ステップ (クーラントの安定条件などに注意)
/*****************************************************************************/
#if _FDTD

/*-------------------- CELL_SIZE:15nm --------------------*/
static const double dt = 2.8e-17; 			// 時間ステップ[s]
static const int Nmax = 100; 				// 最終時間ステップ
/*-------------------- CELL_SIZE:15nm --------------------*/

static const int Ncut = 50000; 				// 時間ステップ数を表示させる間隔
static const int Tcut = 1000; 				// パワーの平均の算出を開始する時間ステップ  (最終計算ステップからの差)
static const int Fcut = 200; 				// フィールドを出力する時間ステップ数 (最終計算ステップからの差)

#else
static const double dt = 2.8e-17; 			// 時間ステップ[s]
static const int Ncut = 5; 				// 時間ステップ数を表示させる間隔
static const int Tcut = 30; 			// エネルギーの平均の算出を開始する時間ステップ
static const int Fcut = 30; 			// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
static const int Nmax = 1; 			// 最終時間ステップ
#endif

static const int Ncheck = 10; 					// 動作確認用のフィールドを出力する時間ステップ
static const int Ncutfield = Ncut; 			// フィールドを出力する時間ステップ
static const int Ncutfield2 = 5; 			// 安定状態でのフィールドを出力する時間ステップ間隔

/*****************************************************************************/
// 物理量[MKSA系]
/*****************************************************************************/
#define PI 3.141592
#define C0 2.997924e8			// 真空中の光速 [m/s]
#define EPSILON0 8.854e-12		// 真空中の誘電率 [F/m]
#define MU0 (PI*(4.0e-7))		// 真空中の透磁率 [N/A^2]


/*****************************************************************************/
// スラブと解析領域
/*****************************************************************************/
static const double dblSlabHeig = SLAB_HEIGHT * 1.0e-9; 				// スラブ厚
static const double dblCladHeight1 = CLAD_HEIGHT1 * 1.0e-9; 				// 上部クラッド層厚
static const int	intSlabHeigPer = INT_DIV (dblSlabHeig/2.0, dblCellSize); 				//スラブ列数(対称境界条件を使用しているので，スラブ厚を1/2)
static const int	intCladHeight1 = INT_DIV (CLAD_HEIGHT1, CELL_SIZE); 				//上部クラッド列数
static const int	air_hc = 0; 					//クラッド上部空気層厚
static const int	intSlabCen = air_hc + intCladHeight1 + intSlabHeigPer; 			//活性層中心セル(+1は[intSlabHeigPer]が奇数セル数のときに中央値にするため)

/*****************************************************************************/
// フォトニック結晶導波路
/*****************************************************************************/
static const double dblPitchCellComp = dblCellSize/2.0; 			// 格子定数の丸め
static const double dblPitch = PITCH * 1e-9; 					// 円孔格子定数

static const double dblRadius = RADIUS * 1e-9; 					// 円孔半径
static const double dblDiamter = 2.0 * dblRadius; 				// 円孔直径

static const int intPitchX = INT_DIV (PITCH, CELL_SIZE); 		// 格子定数のセルサイズ(X方向)
static const int intPitchY = (INT_DIV((PITCH * sqrt(3.0)/2 + 0.5), CELL_SIZE)); 		// 格子定数のセルサイズ(Y方向)	+0.5は四捨五入のため
static const int intRadius = INT_DIV (RADIUS, CELL_SIZE); 		// 円孔半径

static const int Row_x = 10;
static const int Row_y = 9;
static const int Wid = 10;


/*****************************************************************************/
// 材料の屈折率と誘電率
/*****************************************************************************/
static const double nw = 4.5; //量子井戸層の屈折率
static const double nsch1 = 3.265; //SCH1層の屈折率 3.265:バンドギャップ波長1100nm
static const double nsch2 = 3.292; //SCH2層の屈折率 3.292:バンドギャップ波長1150nm
static const double nsch3 = 3.32; //SCH3層の屈折率 3.32:バンドギャップ波長1200nm
static const double np = 3.17; //支柱の屈折率

static const double n_core = 3.43; 		// コアの屈折率(CORE)
static const double n_clad = 1.444; 		// クラッド屈折率(SiO2)

static const double epsilon0 = EPSILON0; 					// 真空の誘電率
static const double epsilon1 = EPSILON0 * SQ(n_core); 		// コアの誘電率
static const double epsilon2 = EPSILON0 * SQ(n_clad); 		// クラッドの誘電率


/*****************************************************************************/
// 励振関数
/*****************************************************************************/
#if _FDTD
#if _EXITATION_FUNC
char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]
#else
char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]
#endif
#else
char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]
#endif

static const double delta_omega = 0.05; 						// 中心周波数で規格した半値全幅
static const int Npeak = 500; 								// ピークステップ数

// 励振点の座標
static const int ex_y_st = YMAX_ALL - 10; 		// 導波路断面始セル数(横) ←解析空間の中間セル座標から導波路幅(既に1/2値になっている)を引いている
static const int ex_y_ed = YMAX_ALL; 					// 導波路断面終セル数(横)
static const int ex_z_st = ZMAX_ALL - intSlabHeigPer; 	// 導波路断面始セル数(縦) ←解析空間の中間セル座標から導波路幅(の1/2値)を引いている
static const int ex_z_ed = ZMAX_ALL; 			// 導波路断面終セル数(縦)


/*****************************************************************************/
// 観測点と励振点の設定 (計算誤差を防ぐために桁上げして計算)
/*****************************************************************************/

// 分割前の励振点，観測面 (分割後はプログラム中で処理)
static const int intExctLen = INT_DIV (EXCT_LEN, CELL_SIZE); //励振点 分割前の励振点

// ポインティングパワーの最大値，最小値の計算区間
static const int intObseInter = INT_DIV (OBSE_INTER, CELL_SIZE); 	// 観測区間

// 分割後の励振点，観測面
int intExctPortNum;	// 励振点がある計算機の番号

int intExctLenPart;	// 励振点

// ポインティングパワーの最大値の算出用

static const int WGlength = (int)((2.00e-6*1e10)/(dx*1e10)); //矩形導波路長(入射側．セル数に変換)
