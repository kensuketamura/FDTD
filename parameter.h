/*****************************************************************************/
// 並列計算時の計算機の台数
/*****************************************************************************/
#if _FDTD
#define NODE 2
#else
#define NODE 2
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

#if _PROGRAM_TEST
/****************************** 本番用(SiO2) ******************************/
/****************************** 3列目格子シフト構造 ******************************/
#define CELL_SIZE 21			// セルサイズ   21
#define PITCH 399				// PC 格子定数   399
#define PITCH_SHIFT_MAX 399 	// 格子定数変化PCWのPC格子定数の最大値
#define PITCH_SHIFT_MAX2 399

#define SLAB_HEIGHT 210		// スラブ厚
#define CLAD_HEIGHT1 1995	// 上部クラッド高さ
#define CLAD_HEIGHT2 0		// 下部クラッド高さ
#define AIR_HEIGHT 0		// 空気層高さ

#define RADIUS 110			// PCの標準円孔半径      120
#define SX3 84				// 伝搬(X)方向の3列目格子シフト量(SX2,SX4=0でないと使えない要改善!!mainの1380行目)
#define SX1 0				// 伝搬(X)方向の1列目格子シフト量
#define SX2 0				// 伝搬(X)方向の2列目格子シフト量(SX3,SX4=0でないと使えない要改善!!mainの1380行目)
#define SX4 0				//伝搬(X)方向の4列目格子シフト量(SX2,SX3=0でないと使えない要改善!!mainの1380行目)
#define SY -588				// 横(Y)方向の導波路全体格子シフト量(対称境界を使用しているので実際にはこの倍)

#define EXCT_LEN 840			//840		// ☆励振点 (モデルの左端からの距離) これは基本変えない．
#define EXCT_OBSE_LEN 3000-420		//975 		// ☆☆励振点から観測面の中心までの距離  これも基本変えない  3000//15/12/25
//OBSE_WIRE_LENを+aしたときにはこれを-aする
#define OBSE_WIRE_LEN 5815		//2540(121 セル) 2750(131セル) 5080(242セル)//  5500(262セル)2/23 //5815(297セル) 35セル増やした 2/24  さらに長尺化(by馬場先生) ☆(出射)観測面の中心から細線導波路端までの距離
//これを+aすると，モニタとPCW端の距離が+a，モニタと出射端の距離が-2aされる．出射の細線距離は-aになる.入射モニタとPCW端との距離も+a
//#define OBSE_INTER (2 * PITCH)		// 観測面の長さ(透過率，反射率を求めるために使用)
#define OBSE_INTER (PITCH)				// 観測面の長さ(透過率，反射率を求めるために使用)
#define WIRE_OUTPUT_LEN (EXCT_LEN + 45)	// 出射細線導波路の長さ(直打ちの数字はセルサイズをNODE数の整数倍にするために使用)
#define WIRE_OUTPUT_OFFSET 0			// 出射細線導波路のスラブ終端の長さ
#define WIRE_WID_OFFSET 378//168				// 細線幅のオフセット量の半分(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define WIRE_WID_OFFSET_OUT 378//168
#define PCW_SiSLAB_TERMINATION_LEN 255	// PCW横のCOREスラブ終端の長さ
#define PCW_SiSLAB_OFFSET 0				// PCW縦のCOREスラブのオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_WIDTH_CHIRP 378//168				// PCW幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_WIDTH_CHIRP_OUT 378//168
#define PCW_WIDTH_Para 1

//☆導波路幅チャープを全て0にするとバグが生じるので便宜的に全て1にする

// これ以下は周期数
#define NORM_PCW_PER 0				// 通常PCW周期数
#define CHIRP_3RD_LS_PER 0			// 3列目格子シフト量チャープLSPCW周期数． 2~4列目まで対応．※シフト量/周期数<cellsizeだとチャープしない
#define CHIRP_2ND_LS_PER 5//5			//導波路幅チャープとシフト量チャープを同時に行う際のチャープLSPCW周期数．PITCH_SHIFT_PERより小さくないとダメ?　2,3列目のみ対応．※シフト量/周期数<cellsizeだとチャープしない
#define CHIRP_2ND_LS_PER_OUT 5//5
#define PITCH_SHIFT_PER 5//5				//導波路幅チャープの周期数
#define PITCH_SHIFT_PER_OUT 5//5
#define PITCH_SHIFT_CHIRP_PER2 5//5		// 格子定数変化チャープPCWの周期数(導波路幅チャープと同時)
#define PITCH_SHIFT_CHIRP_PER2_OUT 5//5
#define LSPCW_SHIFT_DESCRETE FALSE			// 格子定数変化PCWのとき，はじめからLSPCWにする場合はFALSE　しない場合はTRUE
#define PITCH_SHIFT_CHIRP_PER 0		// 格子定数変化チャープPCWの周期数
#define LSPCW_PER 20				// LSPCW周期数
#define PCW_WID 8				// PCWの列数
// これ以上は周期数

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// 通常PCW長
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// チャープLSPCW長
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW長
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// 出射細線導波路のスラブ終端の長さ

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// 入射観測面の中心座標 (モデルの左端からの距離)
//この位置はおそらく図には表示されない　?intObseLenPart4とは向きが逆
//PS.15/12/22 これは機能してない



#define WIRE_LEN1 (OBSE_LEN1 + 3696)			// ☆入射細線長　//2541 長尺 //2961　さらに長尺2/23 //さらに長尺2 3696 2/24
//#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)　//元のプログラム
//#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)

#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN) // ☆☆ 出射細線長 数字はnm　出射端から出射モニタの距離
//モニターを後ろにする場合はオフセット137*21
//#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW長

/****************************** 3列目格子シフト構造 ******************************/

#endif
/****************************** 本番用(SiO2) ******************************/
#endif

#else

/****************************** 解析領域決定用 ******************************/

#define CELL_SIZE 35		// セルサイズ
#define PITCH 455			// PC 格子定数
#define PITCH_SHIFT_MAX 455 // 格子定数変化PCWのPC格子定数の最大値

#define SLAB_HEIGHT 210		// スラブ厚
#define CLAD_HEIGHT1 525	// 上部クラッド高さ
#define CLAD_HEIGHT2 0		// 下部クラッド高さ
#define AIR_HEIGHT 0		// 空気層高さ

#define RADIUS 160			// PCの標準円孔半径
#define SX3 105				// 伝搬(X)方向の3列目格子シフト量
#define SX1 0				// 伝搬(X)方向の3列目格子シフト量
#define SY 0				// 横(Y)方向の導波路全体格子シフト量

#define EXCT_LEN 840					// 励振点 (モデルの左端からの距離)
#define EXCT_OBSE_LEN 980				// 励振点から観測面の中心までの距離
#define OBSE_WIRE_LEN 2520				// 観測面の中心から細線導波路端までの距離
#define OBSE_INTER (2 * PITCH)			// 観測面の長さ(透過率，反射率を求めるために使用)
#define WIRE_OUTPUT_LEN EXCT_LEN		// 出射細線導波路の長さ
#define WIRE_OUTPUT_OFFSET 0			// 出射細線導波路のスラブ終端の長さ
#define WIRE_WID_OFFSET 0				// 細線幅のオフセット量(プラス側が幅広)
#define PCW_SiSLAB_TERMINATION_LEN 105	// PCW横のCOREスラブ終端の長さ
#define PCW_SiSLAB_OFFSET 0				// PCW縦のCOREスラブのオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_WIDTH_CHIRP 0				// PCW幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)

// これ以下は周期数
#define NORM_PCW_PER 0				// 通常PCW周期数
#define CHIRP_3RD_LS_PER 0			// 3列目格子シフト量チャープLSPCW周期数
#define PITCH_SHIFT_PER 5			// 格子定数変化PCWの周期数
#define PITCH_SHIFT_CHIRP_PER 0		// 格子定数変化チャープPCWの周期数
#define LSPCW_PER 20				// LSPCW周期数
#define PCW_WID 8					// PCWの列数
// これ以上は周期数

#define LSPCW_ROW 3					// 導波路から数えて何列目の格子点を伝搬方向にシフトさせるか．

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)		// 通常PCW長
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// チャープLSPCW長
#define LSPCW_LEN (LSPCW_PER * PITCH)			// LSPCW長
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)			// 出射細線導波路のスラブ終端の長さ
// これ以上は周期数

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)	// 入射観測面の中心座標 (モデルの左端からの距離)
#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)		// 出射観測面の中心座標 (モデルの右端からの距離)

#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)						// 入射細線長
#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)					// 出射細線長
#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW長

/****************************** 解析領域決定用 ******************************/

#endif

/*****************************************************************************/
// 全解析領域
/*****************************************************************************/
#if _PROGRAM_TEST

#define XMAX_ALL 1240 //948(2015/11/11) //1250(2015/11/19) 1250-165=1085(2015/12/16) 1250-165-16=1069(2015/12/24) 1069+151-55=1165(2015/12/25)入射長尺化  1220 とすると1206になる(2065/2/23)入射長尺化1
//1222 で 1220 //1240 で 1235 //1260 で 1235 //1270 で 1235 1240より大きくしても意味ないので基本1240
#define YMAX_ALL 172	//10列では145+32(この場合XMAXを50程度減らす必要あり)．8列では145.長尺では172
#define ZMAX_ALL 100	//15/12/16 70(容量不足のため) 従来100

#else

#define XMAX_ALL 616
#define YMAX_ALL 95
#define ZMAX_ALL 18

#endif

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

#if _PROGRAM_TEST

static const double dt = 38e-18; 			// 時間ステップ[s]
static const int Nmax = 150000; //150000 				// ☆最終時間ステップ
static const int Ncut = 1000; //50000				// ☆時間ステップ数を表示させる間隔
static const int Tcut = 1000; 	//1000			//　☆☆パワーの平均の算出を開始する時間ステップ  (最終計算ステップからの差)
static const int Fcut = 200; 				// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
//時間ステップ数が10000だと，
#else

static const double dt = 67e-18; 			// 時間ステップ[s]

/******************** それなり短く (<5min) ********************/
static const int Ncut = 500; 			// 時間ステップ数を表示させる間隔
static const int Tcut = 200; 			// エネルギーの平均の算出を開始する時間ステップ
static const int Fcut = 200; 			// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
static const int Nmax = 1000; 			// 最終時間ステップ

#endif

#else
static const double dt = 2.8e-17; 			// 時間ステップ[s]
static const int Ncut = 5; 				// 時間ステップ数を表示させる間隔
static const int Tcut = 30; //30			// ☆エネルギーの平均の算出を開始する時間ステップ

static const int Fcut = 30; 			// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
static const int Nmax = 1; 			// 最終時間ステップ
#endif

static const int Ncheck = 10; //10					// ☆動作確認用のフィールドを出力する時間ステップ
// これは初めにあるやつ
static const int Ncutfield = Ncut; 			// フィールドを出力する時間ステップ
//static const int Ncutfield2 = 10; 			// 安定状態でのフィールドを出力する時間ステップ間隔
static const int Ncutfield2 = 5; 			// ☆安定状態でのフィールドを出力する時間ステップ間隔
//これは149800～150000まで5刻みで出力する

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
static const int	air_hc = 	(int)((0.0e-6*1e10)/(dz*1e10)); 					//クラッド上部空気層厚
static const int	intSlabCen = air_hc + intCladHeight1 + intSlabHeigPer; 			//活性層中心セル(+1は[intSlabHeigPer]が奇数セル数のときに中央値にするため)


/*****************************************************************************/
// フォトニック結晶導波路
/*****************************************************************************/
static const double dblPitchCellComp = dblCellSize/2.0; 			// 格子定数の丸め
static const double dblPitch = PITCH * 1e-9; 					// 円孔格子定数

static const double dblRadius = RADIUS * 1e-9; 					// 円孔半径
static const double dblRadius2 = 0.18e-6; 					// 大円孔半径
static const double dblRadius3 = 0.12e-6; 					// 入射部円孔の半径
static const double dblRadius4 = 0.04e-6; 					// 場所確認用ドットの半径
static const double dblRadius5 = 0.06e-6; 						//中央円孔列先端円孔半径
static const double dblRadius6 = 0.10e-6; 						//中央円孔列先端円孔半径
static const double dblRadius7 = 0.14e-6; 						//中央円孔列先端円孔半径
static const double dblRadius8 = 0.16e-6; 						//中央円孔列先端円孔半径
static const double dblDiamter = 2.0 * dblRadius; 				// 円孔直径

static const int intPitchX = INT_DIV (PITCH, CELL_SIZE); 		// 格子定数のセルサイズ(X方向)
static const int intPitchY = (INT_DIV((PITCH * sqrt(3.0)/2 + 0.5), CELL_SIZE)); 		// 格子定数のセルサイズ(Y方向)	+0.5は四捨五入のため
static const int intRadius = INT_DIV (RADIUS, CELL_SIZE); 		// 円孔半径

static const int intNormPcwPer = NORM_PCW_PER; 				// PCW列数
static const int intPitchShiftPcwPer = PITCH_SHIFT_PER;				// 格子定数変化PCWの周期数
static const int intPitchShiftPcwPerOut = PITCH_SHIFT_PER_OUT;
static const int intPitchShiftChirpPcwPer = PITCH_SHIFT_CHIRP_PER;	// 格子定数変化チャープPCWの周期数

static const int intChirp3rdLsPer = CHIRP_3RD_LS_PER; 		// チャープLSPCW列数
static const int intChirp2ndLsPer = CHIRP_2ND_LS_PER; 		// チャープLSPCW列数
static const int intChirp2ndLsPerOut = CHIRP_2ND_LS_PER_OUT;
static const int intLspcwPer = LSPCW_PER; 					// LSPCW列数
static const int intPcwPer = 2 * (intNormPcwPer + intChirp3rdLsPer + intPitchShiftChirpPcwPer) + intPitchShiftPcwPer + intPitchShiftPcwPerOut + intLspcwPer; 	// 全PCW列数
static const int intSx4Per = INT_DIV (SX4, CELL_SIZE); 		// 伝搬(X)方向の3列目格子シフト列数
static const int intSx3Per = INT_DIV (SX3, CELL_SIZE); 		// 伝搬(X)方向の3列目格子シフト列数
static const int intSx2Per = INT_DIV (SX2, CELL_SIZE); 		// 伝搬(X)方向の3列目格子シフト列数
static const int intSx1Per = INT_DIV (SX1, CELL_SIZE); 		// 伝搬(X)方向の1列目格子シフト列数
static const int intSyPer = INT_DIV (SY, CELL_SIZE); 		// 横(Y)方向の導波路全体格子シフト列数


//static const int intNormPcwLen = INT_DIV (NORM_PCW_LEN, CELL_SIZE); 				// PCW列数
//static const int intChirpLsLen = INT_DIV (CHIRP_3RD_LS_LEN, CELL_SIZE); 				// チャープLSPCW列数
//static const int intLspcwLen = INT_DIV (LSPCW_LEN, CELL_SIZE); 					// LSPCW列数
//static const int intPcwLen = 2 * (intNormPcwLen + intChirpLsLen) + intLspcwLen; 	// 全PCW列数
static const int intPcwWid = PCW_WID; 											// 横方向のPCW列数
static const int intPcwStartX = intRadius; 										// 円柱配置の開始X列数
static const int intPcwStartY = intRadius + INT_DIV (PCW_SiSLAB_TERMINATION_LEN, CELL_SIZE);	// 円柱配置の開始Y列数


/*****************************************************************************/
// フォトニック結晶導波路接続用 細線導波路
/*****************************************************************************/
static const double dblWireLen1 = WIRE_LEN1 * 1e-9; 		// 入射
static const double dblWireLen2 = WIRE_LEN2 * 1e-9; 		// 出射

static const int intWireLen1 = INT_DIV(WIRE_LEN1, CELL_SIZE); 		// 入射
static const int intWireLen2 = INT_DIV(WIRE_LEN2, CELL_SIZE); 		// 出射
int intWirePer2, intWirePer3; 		// 出射
static const int intWireWid_2 = intPitchY + INT_DIV(WIRE_WID_OFFSET, CELL_SIZE) + INT_DIV(SY, CELL_SIZE); 		// 細線幅の半分のセル数☆
static const int intWireWid_2_Out = intPitchY + INT_DIV(WIRE_WID_OFFSET_OUT, CELL_SIZE) + INT_DIV(SY, CELL_SIZE);
static const int intWireWid = intWireWid_2 * 2; 		// 細線幅☆




/*****************************************************************************/
// 材料の屈折率と誘電率
/*****************************************************************************/
static const double nw = 4.5; //量子井戸層の屈折率
static const double nsch1 = 3.265; //SCH1層の屈折率 3.265:バンドギャップ波長1100nm
static const double nsch2 = 3.292; //SCH2層の屈折率 3.292:バンドギャップ波長1150nm
static const double nsch3 = 3.32; //SCH3層の屈折率 3.32:バンドギャップ波長1200nm
static const double np = 3.17; //支柱の屈折率
//static const double n_core = 3.5; 		// コアの屈折率(CORE)
static const double n_core = 3.43; 		// コアの屈折率(CORE)
#if PCW_Air_Or_SiO2
static const double n_clad = 1.0; 		// クラッド屈折率(真空)
#else
static const double n_clad = 1.444; 		// クラッド屈折率(SiO2)
#endif
static const double epsilon0 = EPSILON0; 					// 真空の誘電率
static const double epsilon1 = EPSILON0 * SQ(n_core); 		// コアの誘電率
static const double epsilon2 = EPSILON0 * SQ(n_clad); 		// クラッドの誘電率


/*****************************************************************************/
// 励振関数//1000
/*****************************************************************************/
#if _FDTD
#if _EXITATION_FUNC
#if _PROGRAM_TEST
#if PCW_Air_Or_SiO2
char *dir_name[] = {"1565","1575"}; 		// 励振関数の波長 [nm]
#else
#if PCW_S1S3_Shift
  char *dir_name[] = {"1570"}; 		// 励振関数の波長 [nm]
#else
  char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]☆計算時
#endif
#endif
#else
  char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]
#endif
#else
  char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]
#endif
#else
  char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]☆吐き出し時
#endif

static const double delta_omega = 0.05; 						// 中心周波数で規格した半値全幅
static const int Npeak = 500; 								// ピークステップ数

// 励振点の座標☆☆
static const int ex_y_st = YMAX_ALL -24 ; 	//- intWireWid_2	// 導波路断面始セル数(横) ←解析空間の中間セル座標から導波路幅(既に1/2値になっている)を引いている☆☆
static const int ex_y_ed = YMAX_ALL; 					// 導波路断面終セル数(横)
static const int ex_z_st = ZMAX_ALL - intSlabHeigPer; 	// 導波路断面始セル数(縦) ←解析空間の中間セル座標から導波路幅(の1/2値)を引いている
static const int ex_z_ed = ZMAX_ALL; 			// 導波路断面終セル数(縦)


static const int y_zx = YMAX_ALL - 1;
/*****************************************************************************/
// 観測点と励振点の設定 (計算誤差を防ぐために桁上げして計算)
/*****************************************************************************/

// 分割前の励振点，観測面 (分割後はプログラム中で処理)
static const int intExctLen = INT_DIV (EXCT_LEN, CELL_SIZE); //励振点 分割前の励振点☆
static const int intObseLen1 = INT_DIV (OBSE_LEN1, CELL_SIZE); //観測点1
//static const int intObseLen5 = XMAX_ALL - INT_DIV (OBSE_LEN5, CELL_SIZE); //☆観測点3 //観測点5の間違えでは？？
//これは多分機能してない

// 分割前の励振点，観測面 (分割後はプログラム中で処理)

/****************************** 観測面の修正(2013/8/20) ******************************/
//static const int intObseWid = 2 * intWireWid_2; // 観測面の幅(Y方向)
static const int intObseWid = 2 * intPitchY; // 観測面の幅(Y方向)
static const int intObseHeig = 2 * intSlabHeigPer; // 観測面の高さ(Z方向)
/****************************** 観測面の修正(2013/8/20) ******************************/

// ポインティングパワーの最大値，最小値の計算区間
static const int intObseInter = INT_DIV (OBSE_INTER, CELL_SIZE); 	// 観測区間

// 分割後の励振点，観測面
int intExctPortNum;	// 励振点がある計算機の番号
int intObseInPortNum;	// 入射 観測点がある計算機の番号
int intObseOutPortNum;	// 出射 観測点がある計算機の番号
int intObseCenPortNum;	// 中央 観測点がある計算機の番号

int intExctLenPart;	// 励振点
int intObseLenPart1; // 入射 観測点(始点)
int intObseLenPart2; // 入射 観測点(終点)
int intObseLenPart3; // 入射 観測点(中点)
int intObseLenPart4; // 出射 観測点(始点)
int intObseLenPart5; // 出射 観測点(終点)
int intObseLenPart6; // 出射 観測点(中点)
int intObseLenPart7; // 中心 観測点(中点)


// 磁界の観測面
int intObseLenPartHz1; 		// 入射 観測点
int intObseLenPartHz5; 		// 出射 観測点


static const int WGlength = (int)((2.00e-6*1e10)/(dx*1e10)); //矩形導波路長(入射側．セル数に変換)
