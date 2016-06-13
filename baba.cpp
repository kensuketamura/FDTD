/*
3次元_FDTD法による電磁界解析 ver. 2.01
From September, 2000; 
Designed by Atsushi SAKAI; 
supported by 
Hiroya DAN (3D-_FDTD), 
Yoshitaka WATANABE (Periodic Boundary Condition in 3D-_FDTD), 
Hiroshi YAMADA(2D-_FDTD, PLRC), 
Tomoki YONEHANA (2D-_FDTD PML, Non-Linear), 
Toshihiko BABA (Photonic crystal bend model: April, 2001), 
Kosuke MORITO (3D_Symmetry_Condition : October, 2003).
Takashi KAWASAKI (PCCW: 2007)
Koichiro YOSHIDA (Observation Area: 2008)
Norihiro ISHIKURA (October, 2012)
*/

#define _FDTD 1		// FDTD計算			0 : モデル掃き出し(プリプロセッサでコンパイルを変更させる)
//										1 : 計算実行

#define _BAND_CALCULATION 0			// 計算の種類 バンド計算
#define _PROPAGATION_CALCULATION 1	// 計算の種類 伝搬計算

#define _CALCULATION_TYPE _PROPAGATION_CALCULATION	// 計算の種類

#define _EXITATION_FUNC 1	// 励振関数の種類		0 : Gaussian 
//													1 : CW

#define _PROGRAM_TEST 1		// プログラムの動作テスト	0: TEST(最終計算ステップ，出力ファイルを短く)
//															1: 本番

#define _MODEL_ALL_EPSILON 0 	// XY断面をモデル全体の出力	0: なし
//																1: あり

#define _CRT_SECURE_NO_WARNINGS //	警告を発生させないようにする

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <direct.h>
#include <stdlib.h>
#include <string>
#include "mpi.h"

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

//#include "grobal_function.cpp"

//#include "parameter.h"
/*****************************************************************************/
// 並列計算時の計算機の台数
/*****************************************************************************/
#if _FDTD
#define NODE 2
#else
#define NODE 1
#endif

// PCWの構造
#define PCW_Air_Or_SiO2	0		// 0:SiO2, 1:Air-brige
#define PCW_S1S3_Shift	0		// 0:3列目格子シフト構造, 1:1,3列目格子シフト構造

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

#if PCW_Air_Or_SiO2
/****************************** 本番用(Air-brige) ******************************/

/*-------------------- CELL_SIZE:15nm --------------------*/
#define CELL_SIZE 15			// セルサイズ
#define PITCH 450				// PC 格子定数
#define PITCH_SHIFT_MAX 450 //480		// 格子定数変化PCWのPC格子定数の最大値

#define SLAB_HEIGHT 210		// スラブ厚
#define CLAD_HEIGHT1 525	// 上部クラッド高さ
#define CLAD_HEIGHT2 0		// 下部クラッド高さ
#define AIR_HEIGHT 0		// 空気層高さ

#define RADIUS 155			// PCの標準円孔半径
#define SX1 0				// 伝搬(X)方向の1列目格子シフト量
#define SX3 90				// 伝搬(X)方向の3列目格子シフト量
#define SY 0				// 横(Y)方向の導波路全体格子シフト量

#define EXCT_LEN 840					// 励振点 (モデルの左端からの距離)
#define EXCT_OBSE_LEN 975				// 励振点から観測面の中心までの距離
#define OBSE_WIRE_LEN 2540				// 観測面の中心から細線導波路端までの距離
//#define OBSE_INTER (2 * PITCH)		// 観測面の長さ(透過率，反射率を求めるために使用)
#define OBSE_INTER (PITCH)				// 観測面の長さ(透過率，反射率を求めるために使用)
#define WIRE_OUTPUT_LEN (EXCT_LEN + 105)	// 出射細線導波路の長さ(直打ちの数字はセルサイズをNODE数の整数倍にするために使用)
#define WIRE_OUTPUT_OFFSET 0			// 出射細線導波路のスラブ終端の長さ
#define WIRE_WID_OFFSET 0				// 細線幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_SiSLAB_TERMINATION_LEN 255	// PCW横のCOREスラブ終端の長さ
#define PCW_SiSLAB_OFFSET 0				// PCW縦のCOREスラブのオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_WIDTH_CHIRP 0				// PCW幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)

// これ以下は周期数
#define NORM_PCW_PER 0				// 通常PCW周期数
#define CHIRP_3RD_LS_PER 0			// 3列目格子シフト量チャープLSPCW周期数
#define PITCH_SHIFT_PER 5			// 格子定数変化PCWの周期数
#define LSPCW_SHIFT_DESCRETE FALSE			// 格子定数変化PCWのとき，はじめからLSPCWにする場合はFALSE　しない場合はTRUE
#define PITCH_SHIFT_CHIRP_PER 0		// 格子定数変化チャープPCWの周期数
#define LSPCW_PER 20				// LSPCW周期数
#define PCW_WID 6					// PCWの列数
// これ以上は周期数
#define LSPCW_ROW 3					// 導波路から数えて何列目の格子点を伝搬方向にシフトさせるか．

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// 通常PCW長
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// チャープLSPCW長
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW長
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// 出射細線導波路のスラブ終端の長さ

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// 入射観測面の中心座標 (モデルの左端からの距離)
#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)	// 出射観測面の中心座標 (モデルの右端からの距離)

#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)			// 入射細線長
#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)		// 出射細線長
#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW長
/*-------------------- CELL_SIZE:15nm --------------------*/

/*-------------------- CELL_SIZE:21nm --------------------*/
//#define CELL_SIZE 21			// セルサイズ
//#define PITCH 441				// PC 格子定数
//#define PITCH_SHIFT_MAX 441 //480		// 格子定数変化PCWのPC格子定数の最大値
//
//#define SLAB_HEIGHT 210		// スラブ厚
//#define CLAD_HEIGHT1 525	// 上部クラッド高さ
//#define CLAD_HEIGHT2 0		// 下部クラッド高さ
//#define AIR_HEIGHT 0		// 空気層高さ
//
//#define RADIUS 147			// PCの標準円孔半径
//#define SX3 105				// 伝搬(X)方向の3列目格子シフト量
//#define SY 0				// 横(Y)方向の導波路全体格子シフト量
//
//#define EXCT_LEN 840					// 励振点 (モデルの左端からの距離)
//#define EXCT_OBSE_LEN 966				// 励振点から観測面の中心までの距離
//#define OBSE_WIRE_LEN 2541				// 観測面の中心から細線導波路端までの距離
////#define OBSE_INTER (2 * PITCH)		// 観測面の長さ(透過率，反射率を求めるために使用)
//#define OBSE_INTER (PITCH)				// 観測面の長さ(透過率，反射率を求めるために使用)
//#define WIRE_OUTPUT_LEN (EXCT_LEN + 0)	// 出射細線導波路の長さ(直打ちの数字はセルサイズをNODE数の整数倍にするために使用)
//#define WIRE_OUTPUT_OFFSET 0			// 出射細線導波路のスラブ終端の長さ
//#define WIRE_WID_OFFSET 0				// 細線幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
//#define PCW_SiSLAB_TERMINATION_LEN 252	// PCW横のCOREスラブ終端の長さ
//#define PCW_SiSLAB_OFFSET 0				// PCW縦のCOREスラブのオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
//#define PCW_WIDTH_CHIRP 0				// PCW幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
//
//// これ以下は周期数
//#define NORM_PCW_PER 0				// 通常PCW周期数
//#define CHIRP_3RD_LS_PER 0			// 3列目格子シフト量チャープLSPCW周期数
//#define PITCH_SHIFT_PER 5			// 格子定数変化PCWの周期数
//#define PITCH_SHIFT_CHIRP_PER 0		// 格子定数変化チャープPCWの周期数
//#define LSPCW_PER 20				// LSPCW周期数
//#define PCW_WID 6					// PCWの列数
//// これ以上は周期数
//
//#define LSPCW_ROW 3					// 導波路から数えて何列目の格子点を伝搬方向にシフトさせるか．
//
//#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// 通常PCW長
//#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// チャープLSPCW長
//#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW長
//#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// 出射細線導波路のスラブ終端の長さ
//
//#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// 入射観測面の中心座標 (モデルの左端からの距離)
//#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)	// 出射観測面の中心座標 (モデルの右端からの距離)
//
//#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)			// 入射細線長
//#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)		// 出射細線長
//#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW長
/*-------------------- CELL_SIZE:21nm --------------------*/

/****************************** 本番用(Air-brige) ******************************/
#else

/****************************** 本番用(SiO2) ******************************/

#if PCW_S1S3_Shift

/****************************** 1,3列目格子シフト構造 ******************************/
#define CELL_SIZE 15			// セルサイズ
#define PITCH 405				// PC 格子定数
#define PITCH_SHIFT_MAX 495 //480		// 格子定数変化PCWのPC格子定数の最大値

#define SLAB_HEIGHT 210		// スラブ厚
#define CLAD_HEIGHT1 750	// 上部クラッド高さ
//#define CLAD_HEIGHT1 990	// 上部クラッド高さ
#define CLAD_HEIGHT2 0		// 下部クラッド高さ
#define AIR_HEIGHT 0		// 空気層高さ

#define RADIUS 110			// PCの標準円孔半径
//#define SX3 90				// 伝搬(X)方向の3列目格子シフト量
//#define SX1 0				// 伝搬(X)方向の3列目格子シフト量
//#define SX3 90				// 伝搬(X)方向の3列目格子シフト量
//#define SX1 120				// 伝搬(X)方向の1列目格子シフト量
#define SX3 0				// 伝搬(X)方向の3列目格子シフト量
#define SX1 0				// 伝搬(X)方向の1列目格子シフト量
#define SY 15				// 横(Y)方向の導波路全体格子シフト量(対称境界を使用しているので実際にはこの倍)

#define EXCT_LEN 840					// 励振点 (モデルの左端からの距離)
#define EXCT_OBSE_LEN 975				// 励振点から観測面の中心までの距離
#define OBSE_WIRE_LEN 2540				// 観測面の中心から細線導波路端までの距離
//#define OBSE_INTER (2 * PITCH)		// 観測面の長さ(透過率，反射率を求めるために使用)
#define OBSE_INTER (PITCH)				// 観測面の長さ(透過率，反射率を求めるために使用)
#define WIRE_OUTPUT_LEN (EXCT_LEN + 60)	// 出射細線導波路の長さ(直打ちの数字はセルサイズをNODE数の整数倍にするために使用)
#define WIRE_OUTPUT_OFFSET 0			// 出射細線導波路のスラブ終端の長さ
#define WIRE_WID_OFFSET 0				// 細線幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_SiSLAB_TERMINATION_LEN 240	// PCW横のCOREスラブ終端の長さ
#define PCW_SiSLAB_OFFSET 0				// PCW縦のCOREスラブのオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)
#define PCW_WIDTH_CHIRP 0				// PCW幅のオフセット量(プラス側が幅広) (負数の場合には，1セル分大きく設定しないと丸め誤差発生)

// これ以下は周期数
#define NORM_PCW_PER 0				// 通常PCW周期数
#define CHIRP_3RD_LS_PER 0			// 3列目格子シフト量チャープLSPCW周期数
#define PITCH_SHIFT_PER 0			// 格子定数変化PCWの周期数
#define PITCH_SHIFT_CHIRP_PER 5	// 格子定数変化チャープPCWの周期数
#define LSPCW_PER 5				// LSPCW周期数
#define PCW_WID 6					// PCWの列数
// これ以上は周期数

//#define LSPCW_ROW 3					// 導波路から数えて何列目の格子点を伝搬方向にシフトさせるか．

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// 通常PCW長
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// チャープLSPCW長
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW長
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// 出射細線導波路のスラブ終端の長さ

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// 入射観測面の中心座標 (モデルの左端からの距離)
#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)	// 出射観測面の中心座標 (モデルの右端からの距離)

#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)			// 入射細線長
#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)		// 出射細線長
//#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW長

/****************************** 1,3列目格子シフト構造 ******************************/

#else

/****************************** 3列目格子シフト構造 ******************************/
#define CELL_SIZE 21			// セルサイズ   15
#define PITCH 399				// PC 格子定数   405
#define PITCH_SHIFT_MAX 399 //480		// 格子定数変化PCWのPC格子定数の最大値
#define PITCH_SHIFT_MAX2 399

#define SLAB_HEIGHT 210		// スラブ厚
#define CLAD_HEIGHT1 1995	// 上部クラッド高さ   750
//#define CLAD_HEIGHT1 990	// 上部クラッド高さ
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

//#define LSPCW_ROW 3					// 導波路から数えて何列目の格子点を伝搬方向にシフトさせるか．

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// 通常PCW長
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// チャープLSPCW長
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW長
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// 出射細線導波路のスラブ終端の長さ

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// 入射観測面の中心座標 (モデルの左端からの距離)
//#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN - 210)	// ☆　出射観測面の中心座標 (モデルの右端からの距離)　これは実際の長さ(nm)
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


/****************************** 動作確認用 ******************************/

//#define CELL_SIZE 100		// セルサイズ
//#define PITCH 400			// PC 格子定数
//
//#define SLAB_HEIGHT 200		// スラブ厚
//#define CLAD_HEIGHT1 700	// 上部クラッド高さ
//#define CLAD_HEIGHT2 0		// 下部クラッド高さ
//#define AIR_HEIGHT 0		// 空気層高さ
//
//#define RADIUS 100			// PCの標準円孔半径
//#define SX3 100				// 伝搬(X)方向の3列目格子シフト量
//#define SY 0				// 横(Y)方向の導波路全体格子シフト量
//
//#define NORM_PCW_PER 5			// 通常PCW周期数
//#define CHIRP_3RD_LS_PER 1			// チャープLSPCW周期数
//#define LSPCW_PER 5				// LSPCW周期数
//#define PCW_WID 6					// PCWの列数
//#define LSPCW_ROW 3					// 導波路から数えて何列目の格子点を伝搬方向にシフトさせるか．
//#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)		// 通常PCW長
//#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// チャープLSPCW長
//#define LSPCW_LEN (LSPCW_PER * PITCH)			// LSPCW長
//
//#define EXCT_LEN 800					// 励振点 (モデルの左端からの距離)
//#define EXCT_OBSE_LEN 1900				// 励振点から観測面の中心までの距離
//#define OBSE_WIRE_LEN EXCT_OBSE_LEN		// 観測面の中心から細線導波路までの距離
//#define OBSE_INTER (2 * PITCH)			// 観測面の長さ(透過率，反射率を求めるために使用)
//#define OBSE_LEN1 (EXCT_LEN + OBSE_INTER / 2 + EXCT_OBSE_LEN)	// 入射観測面の中心座標 (モデルの左端からの距離)
//#define OBSE_LEN5 (EXCT_LEN + OBSE_INTER / 2)		// 出射観測面の中心座標 (モデルの右端からの距離)
//
//#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)	// 入射細線長
//#define WIRE_LEN2 (OBSE_LEN5 + OBSE_WIRE_LEN)	// 出射細線長
//#define WIRE_WID_OFFSET 0		// 細線幅のオフセット量

/****************************** 動作確認用 ******************************/

#endif

/*****************************************************************************/
// 全解析領域
/*****************************************************************************/
#if _PROGRAM_TEST

#if PCW_Air_Or_SiO2
/*-------------------- CELL_SIZE:15nm --------------------*/
#define XMAX_ALL 1440	// SiO2 1340, Air 1440
#define YMAX_ALL 183	// SiO2 163 Air 183
#define ZMAX_ALL 42
/*-------------------- CELL_SIZE:15nm --------------------*/
#else
#if PCW_S1S3_Shift
/*-------------------- CELL_SIZE:15nm --------------------*/
#define XMAX_ALL 1725	// SiO2 1340
#define YMAX_ALL 180	// PCW_WID:6/163  PCW_WID:8/209  PCW_WID:10/255
#define ZMAX_ALL 57		// CLAD_HEIGHT1:525/42  CLAD_HEIGHT1:750/57  CLAD_HEIGHT1:990/73
//#define ZMAX_ALL 73	// CLAD_HEIGHT1:750/57  CLAD_HEIGHT1:990/73
/*-------------------- CELL_SIZE:15nm --------------------*/
#else
/*-------------------- CELL_SIZE:15nm --------------------*/
//☆15nmは気にしなくて良い．
#define XMAX_ALL 1240 //948(2015/11/11) //1250(2015/11/19) 1250-165=1085(2015/12/16) 1250-165-16=1069(2015/12/24) 1069+151-55=1165(2015/12/25)入射長尺化  1220 とすると1206になる(2065/2/23)入射長尺化1
//1222 で 1220 //1240 で 1235 //1260 で 1235 //1270 で 1235 1240より大きくしても意味ないので基本1240
#define YMAX_ALL 172	//10列では145+32(この場合XMAXを50程度減らす必要あり)．8列では145.長尺では172
#define ZMAX_ALL 100	//15/12/16 70(容量不足のため) 従来100
//#define ZMAX_ALL 73	// CLAD_HEIGHT1:750/57  CLAD_HEIGHT1:990/73
/*-------------------- CELL_SIZE:15nm --------------------*/
#endif
#endif

/*-------------------- CELL_SIZE:21nm --------------------*/
//#define XMAX_ALL 1010 // SiO2 1340, Air 695
//#define YMAX_ALL 127 //SiO2 163 Air 183
//#define ZMAX_ALL 30
/*-------------------- CELL_SIZE:21nm --------------------*/

#else

#define XMAX_ALL 616
#define YMAX_ALL 95
#define ZMAX_ALL 18

//#define XMAX_ALL 950
//#define YMAX_ALL 114
//#define ZMAX_ALL 38

//#define XMAX_ALL 162
//#define YMAX_ALL 19
//#define ZMAX_ALL 8

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
#if PCW_Air_Or_SiO2
/*-------------------- CELL_SIZE:15nm --------------------*/
static const double dt = 28e-18; 			// 時間ステップ[s]
static const int Nmax = 300000; 				// 最終時間ステップ
/*-------------------- CELL_SIZE:15nm --------------------*/
#else
/*-------------------- CELL_SIZE:15nm --------------------*/ //☆
static const double dt = 38e-18; 			// 時間ステップ[s]
static const int Nmax = 150000; //150000 				// ☆最終時間ステップ
//static const int Nmaxp = 5000;  				// ☆準最終時間ステップ
/*-------------------- CELL_SIZE:15nm --------------------*/
#endif

/*-------------------- CELL_SIZE:21nm --------------------*/
//static const double dt = 39e-18; 			// 時間ステップ[s]
//static const int Nmax = 130000; 				// 最終時間ステップ
//static const int Nmax = 250000; 				// 最終時間ステップ
/*-------------------- CELL_SIZE:21nm --------------------*/


static const int Ncut = 1000; //50000				// ☆時間ステップ数を表示させる間隔

static const int Tcut = 1000; 	//1000			//　☆☆パワーの平均の算出を開始する時間ステップ  (最終計算ステップからの差)
//static const int Fcut = 500; 				// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
static const int Fcut = 200; 				// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
//時間ステップ数が10000だと，
#else

static const double dt = 67e-18; 			// 時間ステップ[s]

/******************** すごく短く (<<1min) ********************/
//static const int Ncut = 5; 				// 時間ステップ数を表示させる間隔
//static const int Tcut = 30; 			// エネルギーの平均の算出を開始する時間ステップ
//static const int Fcut = 99; 			// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
//static const int Nmax = 100; 			// 最終時間ステップ

/******************** それなり短く (<5min) ********************/
static const int Ncut = 500; 			// 時間ステップ数を表示させる間隔
static const int Tcut = 200; 			// エネルギーの平均の算出を開始する時間ステップ
static const int Fcut = 200; 			// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
static const int Nmax = 1000; 			// 最終時間ステップ

/******************** やや短く (<10min) ********************/
//static const int Ncut = 500; 				// 時間ステップ数を表示させる間隔
//static const int Tcut = 3000; 			// エネルギーの平均の算出を開始する時間ステップ
//static const int Fcut = 500; 				// フィールドを出力する時間ステップ数 (最終計算ステップからの差)
//static const int Nmax = 10000; 			// 最終時間ステップ

/******************** 長く********************/
//static const int Ncut = 20000; 				// 時間ステップ数を表示させる間隔
//static const int Nmax = 100000; 			// 最終時間ステップ
//static const int Tcut = Nmax - 100; 			// エネルギーの平均の算出を開始する時間ステップ
//static const int Fcut = 200; 				// フィールドを出力する時間ステップ数 (最終計算ステップからの差)

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

//static const int intPcwWid = 9; 		// 円孔行数(層厚方向)(Ex.13)
//static const int intPcwLen = 51; 		// 円孔列数(導波路方向)(Ex.10)
//static const int intPcwStartX = 8; 		// 円柱配置の開始X座標(導波路方向，セル数入力)(Ex.105)
//static const int intPcwStartY = 8; 		// 円柱配置の開始Y座標(幅方向，セル数入力)(Ex.4)
//static const int PCmargin = 6; 		//フォトニック結晶領域マージン．この長さだけ最外円孔の中心から伝搬方向にスペースを設ける(セル数入力，"0" = 界面間隔0.08um)

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

//#define taper_x		INT_DIV (0.0e-6, dx)		//テーパ距離(導波路方向) Ex: (int)((0.44e-6*1e10)/(dx*1e10))
//#define taper_y		INT_DIV (0.0e-6, dx)		//テーパ距離(幅方向) Ex: (int)((0.24e-6*1e10)/(dy*1e10))
//#define taper		INT_DIV (0.0e-6, dx)		//テーパ距離(導波路・幅方向に等距離の場合) Ex: (int)((0.44e-6*1e10)/(dx*1e10))


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
//char *dir_name[] = {"1525","1545"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1525","1545"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1525","1545","1565","1585"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1585","1588","1582","1555"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1582","1565","1575"}; 		// 励振関数の波長 [nm]
char *dir_name[] = {"1565","1575"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1565"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1582"}; 		// 励振関数の波長 [nm]
#else
#if PCW_S1S3_Shift
//char *dir_name[] = {"1485","1525"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1575","1565","1555","1545","1535","1525"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1565","1568","1571","1574"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1559","1556","1553","1550"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1547","1544","1541","1538"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1580","1570","1560","1550","1540"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1573","1575","1567"}; 		// 励振関数の波長 [nm]
char *dir_name[] = {"1570"}; 		// 励振関数の波長 [nm]
#else
//char *dir_name[] = {"1530", "1535"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1545", "1550"}; 		// 励振関数の波長 [nm]
/*const*/ char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]☆計算時
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
//static const int ex_y_st = YMAX_ALL - 26; 		// 導波路断面始セル数(横) ←解析空間の中間セル座標から導波路幅(既に1/2値になっている)を引いている
//static const int ex_y_ed = YMAX_ALL; 			// 導波路断面終セル数(横)
//static const int ex_z_st = ZMAX_ALL - 7; 		// 導波路断面始セル数(縦) ←解析空間の中間セル座標から導波路幅(の1/2値)を引いている
//static const int ex_z_ed = ZMAX_ALL; 			// 導波路断面終セル数(縦)



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
//int intObseLenPart8; // 中心 観測点(中点)
//int intObseLenPart9; // 中心 観測点(中点)

// 磁界の観測面
int intObseLenPartHz1; 		// 入射 観測点
int intObseLenPartHz5; 		// 出射 観測点

// ポインティングパワーの最大値の算出用
//static const int intObseLenPart2 = (int)((3.84e-6*1e10)/(dx*1e10)); //intObseLenPart1から格子定数1周期ずれ
//static const int intObseLenPart2 = (int) INT_DIV ((2.7e-6 + dblPitch*2), dx); //intObseLenPart1から格子定数2周期ずれ
//static const int intObseLenPart5 = (int)((1.32e-6*1e10)/(dx*1e10)); //intObseLenPart4から格子定数1周期ずれ
//static const int intObseLenPart5 = (int)((1.3e-6*1e10 + dblPitch*2*1e10)/(dx*1e10)); //intObseLenPart4から格子定数2周期ずれ NODE 6
//static const int intObseLenPart5 = (int) INT_DIV ((9.3e-6 + dblPitch*2), dx); //intObseLenPart4から格子定数2周期ずれ NODE 2

static const int WGlength = (int)((2.00e-6*1e10)/(dx*1e10)); //矩形導波路長(入射側．セル数に変換)


//#include "module0.h"

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
double Ex[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
double Ey[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
double Ez[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
double Hx[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
double Hy[XMAX+1][YMAX+1][ZMAX+1] = {0.0};
double Hz[XMAX+1][YMAX+1][ZMAX+1] = {0.0};

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

int main(int argc, char **argv){

#if _CALCULATION_TYPE == _PROPAGATION_CALCULATION
	double s_time, e_time; 
	char processor_name[MPI_MAX_PROCESSOR_NAME]; 
	char time[9]; 
	int tag_send = 0, tag_recv = 0; 
	int right, left; 
	int namelen; 

	MPI_Status status; 

	// MPIによる通信の開始
	MPI_Init (&argc, &argv); 
	MPI_Comm_rank (MPI_COMM_WORLD, &irank);
    MPI_Comm_size (MPI_COMM_WORLD, &isize);
	MPI_Get_processor_name (processor_name, &namelen);

	if (isize != ISIZE){
		printf ("MPIで設定した計算機の台数(%d)がプログラム中の値と一致しません．\n終了します¥n", ISIZE); 
		return 0; 
	}

	printf ("%d分割並列処理スタート¥n", isize); 
	printf ("Process %d on %s¥n", irank, processor_name); 

	// 隣の計算機の番号の指定
	left = irank - 1; 
	if(irank == IRANK_MIN){
		left = MPI_PROC_NULL; 
	}
	right = irank + 1; 
	if(irank == IRANK_MAX){
		right = MPI_PROC_NULL; 
	}

	// dir_name (励振波長) の配列長だけ繰り返し
	for(int dir_count = 0; dir_count < (sizeof(dir_name) / sizeof(dir_name[0]) ); dir_count++){

		initialize_matrix(); 						// 配列の初期化
		modeling(); 								// モデルの設定
		file_open(dir_name[dir_count]); 			// ファイルを開く
		parameter(dir_name[dir_count]); 			// パラメータの設定と出力


		// 計算開始時刻の出力
		if (irank == IRANK_MIN){
			_strtime(time); 
			fprintf(fpparameter, "Start Time:¥t %s¥n", time); 
			s_time = MPI_Wtime(); 
		}

		// 電磁界計算
		for(n = 1 ; n <= Nmax; n++){

			// 時間ステップ数の表示
			if(n % Ncut == 0){
				_strtime(time); 
				printf("n = %d, ¥t¥t", n); 
				printf("time = %s¥n", time); 
			}

			// 励振関数の設定
			source_func(); 

#if _FDTD

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
#if _EXITATION_FUNC
			calc_poynting_power(); 					
#else
			calc_poynting_powerHz(); 					
#endif

			// 一度同期をとる
			MPI_Barrier(MPI_COMM_WORLD); 
#endif
			if(n == 1) {
				observation_func(); 	// 観測点の設定
				output_model(); 		// モデルの出力
				set_epsilon(); 			// 誘電率の割り当て
			}
		}

		if (irank == IRANK_MIN){
			_strtime(time);
			fprintf(fpparameter, "End Time:¥t %s¥n", time); 	/*計算終了時刻の出力*/
			//時刻の出力
			e_time = MPI_Wtime(); 
			printf ("¥ntime = %f¥n", e_time - s_time); 
		}

		file_close(); 			// ファイルを閉じる
	}

	//MPI_Finalize(); 			// MPIを終了する
#elif _CALCULATION_TYPE == _BAND_CALCULATION

#endif
}



// 出力用ファイルを開く
void file_open(char* dir_name_def){
	char dir_name[40]; 
	//_mkdir(strcpy(dir_name, dir_name_def)); 		// 振り分けできるかテスト

	if (irank == IRANK_MIN){
		fpparameter = fopen (strcat(strcpy(dir_name, dir_name_def), "/Parameter.txt"), "w"); 
		allmodel_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_xy.txt"), "w"); 
		fpallepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/All_Epsilon_xy.txt"), "w"); 
	}

	if(irank == intObseInPortNum){ // 入射
		allmodel_yz1 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_yz1.txt"), "w"); 
	}
	if(irank == intObseOutPortNum){ // 出射
		allmodel_yz4 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_yz4.txt"), "w"); 
	}
	if(irank == intObseCenPortNum){			// 中央
		allmodel_yz7 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_yz7.txt"), "w"); 
	}

	model_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xy.txt"), "w"); 		// 振り分けできるかテスト
	model_xy2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xy2.txt"), "w"); 
	model_yz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_yz.txt"), "w"); 
	model_xz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xz.txt"), "w"); 

	fpepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_xy.txt"), "w"); 
	fpepsilony = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz.txt"), "w"); 
	fpepsilonz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx.txt"), "w"); 
	fpepsilony2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz2.txt"), "w"); 
	fpepsilonz2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx2.txt"), "w"); 

	if(irank == intObseInPortNum){
		fppoynt1 = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower1.txt"), "w");
		avpoynt1 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePoyntingPower1.txt"), "w");
		fpHz1 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Hz01.txt"), "w");
	}

    if(irank == intObseOutPortNum){
		fppoynt5 = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower5.txt"), "w"); 
		avpoynt5 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePoyntingPower5.txt"), "w");
		fpHz5 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Hz05.txt"), "w");
	}
}


/*出力用ファイルを閉じる*/
void file_close(){

	if (irank == IRANK_MIN){
		fclose(fpparameter); 
		//fclose(allmodel_xy); 
		fclose(fpallepsilonx); 
	}

	fclose(model_xy);
	fclose(model_xy2);
	fclose(model_yz); 
	fclose(model_xz); 

	fclose(fpepsilonx); 
	fclose(fpepsilony); 
	fclose(fpepsilonz); 
	fclose(fpepsilony2); 
	fclose(fpepsilonz2); 

	if(irank == intObseInPortNum){
		fclose(fppoynt1);
		fclose(avpoynt1);
		fclose(fpHz1);
	}
	if(irank == intObseOutPortNum){
		fclose(fppoynt5); 
		fclose(avpoynt5);
		fclose(fpHz5);
	}
}


// 計算用パラメータの設定と出力
void parameter(char* dir_name){

	if (irank == IRANK_MIN){

		//printf("Nodes = %d¥n", NODE); 
		//printf("Cell Size[nm] = %d¥n", CELL_SIZE); 
		//printf("Time Step[s] = %e¥n", dt); 
		//printf("Final Time Step: %d¥n", Nmax); 
		fprintf(fpparameter, "XMAX_ALL = %d¥n", XMAX_ALL); 
		fprintf(fpparameter, "YMAX_ALL = %d¥n", YMAX_ALL); 
		fprintf(fpparameter, "ZMAX_ALL = %d¥n", ZMAX_ALL); 
		fprintf(fpparameter, "Nodes = %d¥n", NODE); 
		fprintf(fpparameter, "Cell Size [nm] = %d¥n", CELL_SIZE); 
		fprintf(fpparameter, "Time Step [s] = %e¥n", dt); 
		fprintf(fpparameter, "Final Time Step = %d¥n", Nmax); 
		fprintf(fpparameter, "Final Time [s] = %e¥n", (double) Nmax * dt); 
		fprintf(fpparameter, "¥n"); 

		//printf("Slab Height [nm] = %d¥n", SLAB_HEIGHT); 
		//printf("Hole Pitch [nm] = %d¥n", PITCH); 
		//printf("Hole Diameter [nm] = %d¥n", RADIUS*2); 
		//printf("Normal PCW Period = %d¥n", NORM_PCW_PER); 
		//printf("Chirped LSPCW Period = %d¥n", CHIRP_3RD_LS_PER); 
		//printf("LSPCW Period = %d¥n", LSPCW_PER); 
		//printf("Hole Column[y] = %d¥n", intPcwWid); 
		//printf("Hole Row[x] = %d ¥n", intPcwPer); 
		//printf("Hole Start Coordinate = (%d, %d)¥n", intPcwStartX, intPcwStartY); 

		fprintf(fpparameter, "Upper Clad Index = %lf¥n", n_clad); 
		fprintf(fpparameter, "Slab Index = %lf¥n", n_core); 
		fprintf(fpparameter, "Upper Height [nm] = %d¥n", CLAD_HEIGHT1); 
		fprintf(fpparameter, "Slab Height [nm] = %d¥n", SLAB_HEIGHT); 
		fprintf(fpparameter, "Hole Pitch [nm] = %d¥n", PITCH); 
		fprintf(fpparameter, "Pitch Shift Max [nm] = %d¥n", PITCH_SHIFT_MAX); 
		fprintf(fpparameter, "Hole Diameter [nm] = %d¥n", RADIUS*2);
		fprintf(fpparameter, "%d-row Shift [nm] = %d¥n", 1, SX1);
		fprintf(fpparameter, "%d-row Shift [nm] = %d¥n", 3, SX3);
		fprintf(fpparameter, "%Y-direction All Shift [nm] = %d¥n", SY);
		//fprintf(fpparameter, "%d-row Shift [nm] = %d¥n", LSPCW_ROW, SX3);
		fprintf(fpparameter, "Normal PCW Period = %d¥n", NORM_PCW_PER); 
		fprintf(fpparameter, "Chirped LSPCW Period = %d¥n", CHIRP_3RD_LS_PER); 
		fprintf(fpparameter, "Pitch Shift PCW Period = %d¥n", PITCH_SHIFT_PER); 
		fprintf(fpparameter, "Pitch Shift Chirp PCW Period = %d¥n", PITCH_SHIFT_CHIRP_PER); 
		fprintf(fpparameter, "LSPCW Period = %d¥n", LSPCW_PER); 
		fprintf(fpparameter, "Hole Column[y] = %d ¥n", intPcwWid); 
		fprintf(fpparameter, "Hole Row[x] = %d ¥n", intPcwPer); 
		fprintf(fpparameter, "Hole Start Coordinate = (%d, %d)¥n", intPcwStartX, intPcwStartY); 
		fprintf(fpparameter, "¥n"); 

		fprintf(fpparameter, "Exctation [nm] = %d¥n", EXCT_LEN);
		fprintf(fpparameter, "Exctation to Observation [nm] = %d¥n", EXCT_OBSE_LEN); 
		fprintf(fpparameter, "Observation to Wire [nm] = %d¥n", OBSE_WIRE_LEN);
		fprintf(fpparameter, "Output Wire to Termination [nm] = %d¥n", WIRE_OUTPUT_LEN);
		fprintf(fpparameter, "Termination Length [nm] = %d¥n", WIRE_OUTPUT_OFFSET);
		if(intWireWid < 0){
			fprintf(fpparameter, "Wire Width [nm] = %d¥n", (intWireWid + 1) * CELL_SIZE);
		}
		else{
			fprintf(fpparameter, "Wire Width [nm] = %d¥n", intWireWid * CELL_SIZE);
		}
		fprintf(fpparameter, "PCW Slab Termination Length [nm] = %d¥n", PCW_SiSLAB_TERMINATION_LEN); 
		if(PCW_SiSLAB_OFFSET < 0){
			fprintf(fpparameter, "PCW Slab Offset [nm] = %d¥n", PCW_SiSLAB_OFFSET + CELL_SIZE); 
		}
		else{
			fprintf(fpparameter, "PCW Slab Offset [nm] = %d¥n", PCW_SiSLAB_OFFSET); 
		}
		if(PCW_WIDTH_CHIRP < 0){
			fprintf(fpparameter, "PCW Width Chirp [nm] = %d¥n", PCW_WIDTH_CHIRP + CELL_SIZE); 
		}
		else{
			fprintf(fpparameter, "PCW Width Chirp [nm] = %d¥n", PCW_WIDTH_CHIRP); 
		}

		if(LSPCW_SHIFT_DESCRETE == FALSE){
			fprintf(fpparameter, "LSPCW_SHIFT_DESCRETE = FALSE¥n");
		}
		else{
			fprintf(fpparameter, "LSPCW_SHIFT_DESCRETE = TRUE¥n"); 
		}


		fprintf(fpparameter, "¥n"); 
	}

	// 励振関数定数の設定
	lambda = atof(dir_name) * 1e-9; 		
	omega0 = 2.0*PI*C0/lambda; 
	sigma = omega0 * delta_omega;
}

/*配列の初期化*/
void initialize_matrix(){

	//各ノードの座標
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

	// 解析空間の最大値
	xmax_all = XMAX_ALL; 
	ymax_all = YMAX_ALL; 
	zmax_all = ZMAX_ALL; 

	// 解析空間の中心座標
	x_cen = xmax/2; 
	y_cen = ymax/2; 
	z_cen = zmax/2; 	

	//モデルの中心と解析空間の中心は１セル分ずれているので要注意
	x_model_cen = x_cen + 1; 
	y_model_cen = y_cen + 1; 

	int x, y, z;

	for (x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax+1; z++){
				// 電界
				Ex[x][y][z] = 0.0; 
				Ey[x][y][z] = 0.0;
				Ez[x][y][z] = 0.0; 

				// 磁界
				Hx[x][y][z] = 0.0;
				Hy[x][y][z] = 0.0; 
				Hz[x][y][z] = 0.0; 
			}	
		}	
	}

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
				ALL_cell[x][y][z] = CLAD;
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
	}

}


// モデルの設定
void modeling(){

	int n_temp; 		//屈折率の値保存用
	double epsilon_temp; 		//誘電率の値保存用

	/****************************** スラブの形成 ******************************/

	for(x = 0; x < xmax_all+1; x++){
		for(y = 0; y < ymax_all; y++){
			for(z = 0; z < zmax_all; z++){		

				//n_temp = CLAD; 
				//epsilon_temp = EPSILON0; 
				n_temp = CLAD; 
				epsilon_temp = epsilon2; 

				if(z < air_hc){			//空気層に設定
				}
				if(z >= air_hc && z < (air_hc + intCladHeight1) ){			//上部クラッドに設定
				}
				if(z >= (air_hc + intCladHeight1) && z < (air_hc + intCladHeight1 + intSlabHeigPer)){	// スラブに設定
					n_temp = CORE; 
					epsilon_temp = epsilon1; 
				}

				ALL_cell[x][y][z] = n_temp; 
				ALL_epsilonx[x][y][z] = epsilon_temp; 
				ALL_epsilony[x][y][z] = epsilon_temp; 
				ALL_epsilonz[x][y][z] = epsilon_temp; 
			}
		}
	}
	/****************************** スラブの形成 ******************************/


	/****************************** フォトニック結晶 ******************************/

	int s_x3; 		// チャープLSPCWのシフト量
	int s_x2;
	int s_x4;
	int z_end; 				// 円孔の終了座標

	if (intPcwPer == 0){
		intWirePer2 = intWireLen1 - 1;									// 出射COREスラブの開始点
		intWirePer3 = intWirePer2 + intWireLen2;						// 出射COREスラブの終了点
	}

	else{

		//struct PNUM Pnum[intPcwWid*10][intPcwPer*10]; 	// 円柱の中心座標
		struct PNUM Pnum[100][10]; 	// 円柱の中心座標
		struct PNUM Pnum_Init[1][10]; 		// 円柱の標準格子定数による中心座標

		z_end = zmax_all; 		// 円孔が貫通している場合を考える
		Pnum_Init[0][0].Y = intPcwStartY; 	// 円柱を配置する最初のY座標

		if(intPcwWid % 2 == 1){		// if y:even 
			Pnum_Init[0][intPcwWid-1].X = intWireLen1 + intPcwStartX + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1	
		}
		else{						// if y:odd 0.5Aずらす
			Pnum_Init[0][intPcwWid-1].X = intWireLen1 + intPcwStartX - 1;	// 配列の引数に使用するので-1
		}

		if(y != 0){
			Pnum_Init[0][intPcwWid-1].Y = Pnum_Init[0][0].Y + intPitchY * (intPcwWid - 1); 		//ｙ座標を(root3)/2*intPitchXだけずらす
		}

		for(z = 0; z < z_end; z++){
			for(y = 0; y < intPcwWid; y++){

				//if(y % 2 == 1){		// if y:even 
				//	Pnum[0][y].X = intWireLen1 + intPcwStartX - 1;	// 配列の引数に使用するので-1
				//}
				//else{				// if y:odd 0.5Aずらす
				//	Pnum[0][y].X = intWireLen1 + intPcwStartX + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1	
				//}

				//if(y != 0){
				//	Pnum[0][y].Y = Pnum[0][y-1].Y + intPitchY; 		//ｙ座標を(root3)/2*intPitchXだけずらす
				//}
				int input_NormPcw_Xend;
				int input_PitchShiftPcw_Xend;
				int input_PitchShiftChirpPcw_Xend;
				int input_Chirp_Ls_Xend;
				int Lspcw_Xend;
				int output_Chirp_Ls_Xend;
				int output_PitchShiftChirpPcw_Xend;
				int output_PitchShiftPcw_Xend;
				int output_PCW_Xend;

				if (intNormPcwPer > 0){
					input_NormPcw_Xend = intNormPcwPer;
				}
				else{
					input_NormPcw_Xend = 0;
				}
				/******************** 格子シフト量チャープ(2013/7/19) ********************/
				input_PitchShiftPcw_Xend = input_NormPcw_Xend + intPitchShiftPcwPer + intPitchShiftChirpPcwPer;
				input_PitchShiftChirpPcw_Xend = input_PitchShiftPcw_Xend;
				//input_PitchShiftPcw_Xend = input_NormPcw_Xend + intPitchShiftPcwPer;
				//if (intPitchShiftChirpPcwPer > 0){
				//	input_PitchShiftChirpPcw_Xend = input_PitchShiftPcw_Xend + (intPitchShiftChirpPcwPer - 1);
				//}
				//else{
				//	input_PitchShiftChirpPcw_Xend = input_PitchShiftPcw_Xend;
				//}
				/******************** 格子シフト量チャープ(2013/7/19) ********************/

				if (intChirp3rdLsPer > 0){
					input_Chirp_Ls_Xend = input_PitchShiftChirpPcw_Xend + (intChirp3rdLsPer);
				}
				else{
					input_Chirp_Ls_Xend = input_PitchShiftChirpPcw_Xend;
				}
				Lspcw_Xend = input_Chirp_Ls_Xend + intLspcwPer + 1;
				if (intChirp3rdLsPer > 0){
					output_Chirp_Ls_Xend = Lspcw_Xend + (intChirp3rdLsPer - 1);
				}
				else{
					output_Chirp_Ls_Xend = Lspcw_Xend;
				}
				/******************** 格子シフト量チャープ(2013/7/19) ********************/
				output_PitchShiftChirpPcw_Xend = output_Chirp_Ls_Xend;
				output_PitchShiftPcw_Xend = output_PitchShiftChirpPcw_Xend + intPitchShiftPcwPer + intPitchShiftChirpPcwPer;
				//if (intPitchShiftChirpPcwPer > 0){
				//	output_PitchShiftChirpPcw_Xend = output_Chirp_Ls_Xend + (intPitchShiftChirpPcwPer - 1);
				//}
				//else{
				//	output_PitchShiftChirpPcw_Xend = output_Chirp_Ls_Xend;
				//}
				//output_PitchShiftPcw_Xend = output_PitchShiftChirpPcw_Xend + intPitchShiftPcwPer;
				/******************** 格子シフト量チャープ(2013/7/19) ********************/

				output_PCW_Xend = output_Chirp_Ls_Xend + intNormPcwPer;

				int intPitchShiftX = INT_DIV (PITCH_SHIFT_MAX, CELL_SIZE);
				int intPitchShiftY = (INT_DIV((PITCH_SHIFT_MAX * sqrt(3.0)/2 + 0.5), CELL_SIZE));
				int intPCWwidthChirp = INT_DIV (PCW_WIDTH_CHIRP, CELL_SIZE);
				int intPreviousPCWwidthOffset;
				int intNowPCWwidthOffset;

				/******************** 格子シフト量チャープ(2013/7/19) ********************/
				double dblPitchShiftChirpX, dblPitchShiftChirpY, dblPitchShiftChirpY2;
				int intPitchShiftChirpX, intPitchShiftChirpX2, intPitchShiftChirpY;
				/******************** 格子シフト量チャープ(2013/7/19) ********************/


				/****************************** LSPCW ******************************/
				for(x = 0; x < intPcwPer; x++){
					int y2, y_poo, y_poo2;

					s_x3 = 0;
					s_x2 = 0;
					s_x4 = 0;
					// 入射 通常PCW
					if (x < input_NormPcw_Xend){
						if (x == 0){
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								if(y2 % 2 == 1){		// if y:even 
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

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

					// 入射 格子定数変化PCW
					else if (x < input_PitchShiftPcw_Xend && x >= input_NormPcw_Xend){
						if (x == 0){
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchShiftY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								if(y2 % 2 == 1){		// if y:even 
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchShiftX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchShiftX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路1列目シフト構造(2013/7/12) ********************/
								if (y2 != intPcwWid - 1){
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
									//Pnum[x][y2].X -= INT_DIV( (SX1 + (SX1 * (double)(intPitchShiftX / intPitchX) * (x - intPitchShiftChirpPcwPer) / (double) intPitchShiftChirpPcwPer)), CELL_SIZE);
									Pnum[x][y2].X -= INT_DIV(SX1, CELL_SIZE);
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
								}
								/******************** 導波路1列目シフト構造(2013/7/12) ********************/

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;
							}
						}
						else{
							/******************** 格子シフト量チャープ(2013/7/19) ********************/
							intPitchShiftChirpX = intPitchShiftX;
							intPitchShiftChirpY = 0;

							if (x >= intPitchShiftPcwPer){
								intPitchShiftChirpX -= (int) ((intPitchShiftX - intPitchX) * ((x - intPitchShiftPcwPer - 1) / (double) (intPitchShiftChirpPcwPer - 1)));

								if (x > intPitchShiftPcwPer){
									dblPitchShiftChirpY = (intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 1);
								}
							}

							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								if(x < intPitchShiftPcwPer + intPitchShiftChirpPcwPer - 1){
									if (y2 < intPcwWid-1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchShiftY - dblPitchShiftChirpY * (x - intPitchShiftPcwPer)) * (intPcwWid-1 - y2));	
									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
									}
								}
								else{
									if (y2 < intPcwWid-1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);	
									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
									}
								}
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
							}
							/******************** 格子シフト量チャープ(2013/7/19) ********************/

						}
						if (LSPCW_SHIFT_DESCRETE == FALSE){
							// 2列目格子シフト
							if (y == intPcwWid - 2){
								s_x2 = INT_DIV(SX2, CELL_SIZE);
							}
							// 3列目格子シフト
							if (y == intPcwWid - 3){
								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								s_x3 = INT_DIV(SX3, CELL_SIZE);
								//s_x3 = INT_DIV( (SX3 + (SX3 * (double)(intPitchShiftX / intPitchX) * (x - intPitchShiftChirpPcwPer) / (double) intPitchShiftChirpPcwPer)), CELL_SIZE);
								//s_x3 = INT_DIV( (SX3 * (double)(intPitchShiftX / intPitchX)), CELL_SIZE);
								/******************** 格子シフト量チャープ(2013/7/19) ********************/

							}
							if (y == intPcwWid - 4){
								s_x4 = INT_DIV(SX4, CELL_SIZE);
							}
						}
					}

					// 入射 チャープLSPCW
					else if (x >= input_NormPcw_Xend && x < input_Chirp_Ls_Xend){
						if (x == 0){
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								if(y2 % 2 == 1){		// if y:even 
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;
							}
						}
						else{
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;

							}

							if (y == intPcwWid - 3){
								if (intChirp3rdLsPer == 0){
									s_x3 = 0;
								}
								else{
									s_x3 = INT_DIV (intSx3Per, intChirp3rdLsPer) * (x - input_NormPcw_Xend); 
								}
							}
							if (y == intPcwWid - 2){
								if (intChirp3rdLsPer == 0){
									s_x2 = 0;
								}
								else{
									s_x2 = INT_DIV (intSx2Per, intChirp3rdLsPer) * (x - input_NormPcw_Xend); 
								}
							}
							if (y == intPcwWid - 4){
								if (intChirp3rdLsPer == 0){
									s_x4 = 0;
								}
								else{
									s_x4 = INT_DIV (intSx4Per, intChirp3rdLsPer) * (x - input_NormPcw_Xend); 
								}
							}
						}
					}


					// 出射 格子定数変化PCW
					else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){

						// 格子定数変化PCWとの入射接続部
						if (x == output_PitchShiftChirpPcw_Xend){
							y_poo = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								//Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchShiftY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								//if(y2 % 2 == 1){		// if y:even 
								//	Pnum[x][y2].X = Pnum[x-1][intPcwWid-1].X + intPitchShiftX;	// 配列の引数に使用するので-1
								//}
								//else{				// if y:odd 0.5Aずらす
								//	Pnum[x][y2].X = Pnum[x-1][intPcwWid-1].X + intPitchShiftX + INT_DIV (intPitchX, 2.0);	// 配列の引数に使用するので-1
								//}
								/******************** 格子シフト量チャープ(2013/7/19) ********************/


								/******************** 導波路1列目シフト構造(2013/7/12) ********************/
								if (y2 != intPcwWid - 1){
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
									Pnum[x][y2].X -= INT_DIV(SX1, CELL_SIZE);
									//Pnum[x][y2].X -= INT_DIV( (SX1 + (SX1 * (double)(intPitchShiftX / intPitchX) * (x - output_PitchShiftChirpPcw_Xend) / (double) intPitchShiftChirpPcwPer)), CELL_SIZE);
									//Pnum[x][y2].X -= INT_DIV(SX1, CELL_SIZE);
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
								}
								/******************** 導波路1列目シフト構造(2013/7/12) ********************/

								y_poo++;
							}
						}							
						/******************** 格子シフト量チャープ(2013/7/19) ********************/
						intPitchShiftChirpX = intPitchX;
						intPitchShiftChirpY = 0;
						dblPitchShiftChirpY = 0;

						if (x >= output_PitchShiftChirpPcw_Xend){
							intPitchShiftChirpX += (int) ((intPitchShiftX - intPitchX) * ((x - output_PitchShiftChirpPcw_Xend) / (double) (intPitchShiftChirpPcwPer - 1)));
							dblPitchShiftChirpY = (intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 1);

							//if (x > output_PitchShiftChirpPcw_Xend){
							//	dblPitchShiftChirpY = (intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 1);
							//}
						}

						for (y2 = intPcwWid-1; y2 >= 0; y2--){
							if (x != output_PitchShiftChirpPcw_Xend){
								if (y2 < intPcwWid-1){
									if (y2 % 2 == 1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchY + dblPitchShiftChirpY * (x - output_PitchShiftChirpPcw_Xend)) * (intPcwWid-1 - y2));	
									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchY + dblPitchShiftChirpY * (x - output_PitchShiftChirpPcw_Xend + 1)) * (intPcwWid-1 - y2));	
									}
								}
								else{
									Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
								}
							}
							else{
								if (y2 < intPcwWid-1){
									if (y2 % 2 == 1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);	
									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchY + dblPitchShiftChirpY * (x - output_PitchShiftChirpPcw_Xend + 1)) * (intPcwWid-1 - y2));	
									}
								}
								else{
									Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
								}
							}

							//if(x > output_PitchShiftChirpPcw_Xend){
							//if (x >= output_PitchShiftChirpPcw_Xend){
							//	if (y2 < intPcwWid-1){
							//		Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchY + dblPitchShiftChirpY * (x - output_PitchShiftChirpPcw_Xend + 1)) * (intPcwWid-1 - y2));	
							//	}
							//	else{
							//		Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
							//	}
							//}
							//else{
							//	if (y2 < intPcwWid-1){
							//		Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);	
							//	}
							//	else{
							//		Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
							//	}
							//}
							//}
							//else{
							//	if (y2 < intPcwWid-1){
							//		Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);	
							//	}
							//	else{
							//		Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
							//	}
							//}
							Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
						}

						////else{
						//intPitchShiftChirpX = intPitchX;
						//intPitchShiftChirpX2 = intPitchX;
						//intPitchShiftChirpY = 0;
						//if (x >= output_PitchShiftChirpPcw_Xend){
						//	intPitchShiftChirpX += (int) ((intPitchShiftX - intPitchX) * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (intPitchShiftChirpPcwPer)));
						//	intPitchShiftChirpX2 += (int) ((intPitchShiftX - intPitchX) * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (intPitchShiftChirpPcwPer - 1)));
						//	//if (x > output_PitchShiftChirpPcw_Xend){
						//	dblPitchShiftChirpY = -(intPitchShiftY - intPitchY) / (double) intPitchShiftChirpPcwPer;
						//	dblPitchShiftChirpY2 = -(intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 1);
						//	intPitchShiftChirpY = (int) ((intPitchShiftY - intPitchY) * (output_PitchShiftChirpPcw_Xend-x-1) / (double) intPitchShiftChirpPcwPer) - (int) ((intPitchShiftY - intPitchY) * (output_PitchShiftChirpPcw_Xend-x) / (double) intPitchShiftChirpPcwPer);
						//	//}
						//}

						//for (y2 = intPcwWid-1; y2 >= 0; y2--){
						//	if (y2 < intPcwWid-1){
						//		if (y2 % 2 == 1){
						//			Pnum[x][y2].Y = Pnum[x-1][y2].Y + (int) (dblPitchShiftChirpY * (double) (intPcwWid-1 - y2));	
						//		}	
						//		else{
						//			Pnum[x][y2].Y = Pnum[x-1][y2].Y + intPitchShiftChirpY + (int) (dblPitchShiftChirpY2 * (double) (intPcwWid-1 - y2 - 1));	
						//		}
						//	}
						//	else{
						//		Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
						//	}
						//	if (y2 % 2 == 1){
						//		Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
						//	}
						//	else{
						//		Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX2;
						//	}
						//}

						////for (y2 = intPcwWid-1; y2 >= 0; y2--){
						////	Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
						////	Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftX;
						////}
						////}
						/******************** 格子シフト量チャープ(2013/7/19) ********************/

						if (LSPCW_SHIFT_DESCRETE == FALSE){
							// 3列目格子シフト
							if (y == intPcwWid - 2){
								s_x2 = INT_DIV(SX2, CELL_SIZE);
							}
							if (y == intPcwWid - 4){
								s_x4 = INT_DIV(SX4, CELL_SIZE);
							}
							if (y == intPcwWid - 3){

								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								s_x3 = INT_DIV(SX3, CELL_SIZE);
								//s_x3 = INT_DIV( (SX3 + (SX3 * (double)(intPitchShiftX / intPitchX) * (x - output_PitchShiftChirpPcw_Xend) / (double) intPitchShiftChirpPcwPer)) , CELL_SIZE);
								//s_x3 = INT_DIV( (SX3 * (double)(intPitchShiftX / intPitchX)), CELL_SIZE);
								/******************** 格子シフト量チャープ(2013/7/19) ********************/
							}
						}
					}

					// 出射 チャープLSPCW
					else if (x >= Lspcw_Xend && x < output_Chirp_Ls_Xend){
						for (y2 = intPcwWid-1; y2 >= 0; y2--){
							Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
							Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;

						}
						if (y == intPcwWid - 3){
							if (intChirp3rdLsPer == 0){
								s_x3 = 0;
							}
							else{
								s_x3 = INT_DIV (intSx3Per, intChirp3rdLsPer) * (output_Chirp_Ls_Xend - x); 
							}
						}
						if (y == intPcwWid - 4){
							if (intChirp3rdLsPer == 0){
								s_x4 = 0;
							}
							else{
								s_x4 = INT_DIV (intSx4Per, intChirp3rdLsPer) * (output_Chirp_Ls_Xend - x); 
							}
						}
						if (y == intPcwWid - 2){
							if (intChirp3rdLsPer == 0){
								s_x2 = 0;
							}
							else{
								s_x2 = INT_DIV (intSx2Per, intChirp3rdLsPer) * (output_Chirp_Ls_Xend - x); 
							}
						}
					}


					// 出射 通常PCW
					else if (x >= output_Chirp_Ls_Xend && x < output_PCW_Xend){
						if (x == 0){
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								if(y2 % 2 == 1){		// if y:even 
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

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

					// 通常格子定数PCW or LSPCW
					else{
						if (x != 0){
							y_poo = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (intPitchY * y_poo);	
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;
								y_poo++;
							}
						}
						///******************** 格子シフト量チャープ(2013/7/19) ********************/
						////// 格子定数変化PCWとの出射接続部
						//if (x == output_Chirp_Ls_Xend - 1){
						//	intPitchShiftChirpX = intPitchX;
						//	intPitchShiftChirpY = 0;

						//	intPitchShiftChirpX += (int) ((intPitchShiftX - intPitchX) * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (intPitchShiftChirpPcwPer - 1)));
						//	//if (x > output_PitchShiftChirpPcw_Xend){
						//	dblPitchShiftChirpY = -(intPitchShiftY - intPitchY) / (double) intPitchShiftChirpPcwPer;
						//	dblPitchShiftChirpY2 = -(intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 2);
						//	intPitchShiftChirpY = (int) ((intPitchShiftY - intPitchY) * (intPitchShiftChirpPcwPer-x-1) / (double) intPitchShiftChirpPcwPer) - (int) ((intPitchShiftY - intPitchY) * (intPitchShiftChirpPcwPer-x) / (double) intPitchShiftChirpPcwPer);
						//	//}

						//	for (y2 = intPcwWid-1; y2 >= 0; y2--){
						//		if (y2 < intPcwWid-1){
						//			if (y2 % 2 == 1){
						//				//Pnum[x][y2].Y = Pnum[x-1][y2].Y + (int) (dblPitchShiftChirpY * (double) (intPcwWid-1 - y2));	
						//			}	
						//			else{
						//				Pnum[x][y2].Y = Pnum[x-1][y2].Y + intPitchShiftChirpY + (int) (dblPitchShiftChirpY2 * (double) (intPcwWid-1 - y2 - 1.5));	
						//			}
						//		}
						//		else{
						//			Pnum[x][y2].Y = Pnum[x-1][y2].Y;	
						//		}
						//		Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
						//	}
						//	//	y_poo = 0;
						//	//	for (y2 = intPcwWid-1; y2 >= 0; y2--){
						//	//		if(y2 % 2 == 0){
						//	//			Pnum[x][y2].X += (intPitchShiftX - intPitchX);
						//	//			Pnum[x][y2].Y -= (intPitchShiftY - intPitchY) * y_poo;
						//	//		}
						//	//		y_poo++;
						//	//	}
						//}
						///******************** 格子シフト量チャープ(2013/7/19) ********************/


						// 3列目格子シフト
						if (y == intPcwWid - 3){
							s_x3 = intSx3Per;
						}
						if (y == intPcwWid - 2){
							s_x2 = intSx2Per;
						}
						if (y == intPcwWid - 4){
							s_x4 = intSx4Per;
						}
					}
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// 導波路幅チャープ
					if (PCW_WIDTH_CHIRP != 0){
						if (x < input_PitchShiftPcw_Xend){
							////////入射シフト量チャープここから(田村)
							if(intChirp2ndLsPer != 0){
							//3列目
							if (y == intPcwWid - 3){
								if (intChirp2ndLsPer == 0){
									s_x3 = 0;
								}
								else if (INT_DIV (intSx3Per, intChirp2ndLsPer) * (x +10 - input_PitchShiftPcw_Xend)<INT_DIV(SX3,CELL_SIZE)){
									s_x3 = INT_DIV (intSx3Per, intChirp2ndLsPer) * (x +10 - input_PitchShiftPcw_Xend); 
								}
								else{
									s_x3 = INT_DIV(SX3,CELL_SIZE);
								}
							}
							//2列目
							if (y == intPcwWid - 2){
								if (intChirp2ndLsPer == 0){
									s_x2 = 0;
								}
								else if (INT_DIV (intSx2Per, intChirp2ndLsPer) * (x +10 - input_PitchShiftPcw_Xend)<INT_DIV(SX2,CELL_SIZE)){
									s_x2 = INT_DIV (intSx2Per, intChirp2ndLsPer) * (x +10 - input_PitchShiftPcw_Xend); 
								}
								else{
									s_x2 = INT_DIV(SX2,CELL_SIZE);
								}
							}
							}
							////////入射シフト量チャープここまで
							if (x == 0){
								intNowPCWwidthOffset = intPCWwidthChirp;
								for (y2 = intPcwWid-1; y2 >= 0; y2--){
									Pnum[x][y2].Y -= intNowPCWwidthOffset;
								}
								intPreviousPCWwidthOffset = intNowPCWwidthOffset;
							}
							else{
								// PCW_WIDTH_CHIRP用微調整
								intNowPCWwidthOffset = intPCWwidthChirp - (int) (intPCWwidthChirp * (x / (double) (input_PitchShiftPcw_Xend-1)));
								//intNowPCWwidthOffset = intPCWwidthChirp - (int) (intPCWwidthChirp * ((x + 1) / (double) (input_PitchShiftPcw_Xend-1)));
								if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
									for (y2 = intPcwWid-1; y2 >= 0; y2--){
										Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
									}
									intPreviousPCWwidthOffset = intNowPCWwidthOffset;
								}
							}
						}
						else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){

							////////出射シフト量チャープここから(田村)
							if(intChirp2ndLsPer != 0){
							//3列目
							if (y == intPcwWid - 3){
								if (intChirp2ndLsPer == 0){
									s_x3 = 0;
								}
								else if( INT_DIV (intSx3Per, intChirp2ndLsPer) * (output_PitchShiftPcw_Xend - x-2)<INT_DIV(SX3,CELL_SIZE)){
									s_x3 = INT_DIV (intSx3Per, intChirp2ndLsPer) * (output_PitchShiftPcw_Xend - x-2); 
								}
								else{
									s_x3 = INT_DIV(SX3,CELL_SIZE);
								}
							}
							//2列目
							if (y == intPcwWid - 2){
								if (intChirp2ndLsPer == 0){
									s_x2 = 0;
								}
								else if( INT_DIV (intSx2Per, intChirp2ndLsPer) * (output_PitchShiftPcw_Xend - x-2)<INT_DIV(SX2,CELL_SIZE)){
									s_x2 = INT_DIV (intSx2Per, intChirp2ndLsPer) * (output_PitchShiftPcw_Xend - x-2); 
								}
								else{
									s_x2 = INT_DIV(SX2,CELL_SIZE);
								}
							}
							}
						//////////出射シフト量チャープここまで

							if (x == output_PitchShiftChirpPcw_Xend){
								//for (y2 = intPcwWid-1; y2 >= 0; y2--){
								//	Pnum[x][y2].Y -= intNowPCWwidthOffset;
								//}
								intPreviousPCWwidthOffset = 0;

								intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (input_PitchShiftPcw_Xend-1)) + 0.9);
								if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
									for (y2 = intPcwWid-1; y2 >= 0; y2--){
										if (y2 % 2 == 0){
											Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
										}
									}
									intPreviousPCWwidthOffset = intNowPCWwidthOffset;
								}
								intPreviousPCWwidthOffset = 0;
							}
							else{
								// PCW_WIDTH_CHIRP用微調整
								//intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend) / (double) (input_PitchShiftPcw_Xend-1)) + 0.8);


								for (y2 = intPcwWid-1; y2 >= 0; y2--){
									if (y2 % 2 == 1){
										intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend) / (double) (input_PitchShiftPcw_Xend-1)) + 0.9);
									}
									else{
										intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (input_PitchShiftPcw_Xend-1)) + 0.9);
									}
									if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
										Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
									}
								}

								intPreviousPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend) / (double) (input_PitchShiftPcw_Xend-1)) + 0.9);
							}
						}				
					}
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//if (x < input_PitchShiftPcw_Xend){
					//	for (y2 = intPcwWid-1; y2 >= 0; y2--){
					//		if (x == 0){
					//				intNowPCWwidthOffset = intPCWwidthChirp;
					//				Pnum[x][y2].Y -= intNowPCWwidthOffset;
					//				intPreviousPCWwidthOffset = intPCWwidthChirp;
					//			}
					//		}
					//		else{
					//			if (PCW_WIDTH_CHIRP != 0){
					//				intNowPCWwidthOffset = intPCWwidthChirp - (int) (intPCWwidthChirp * (x / (double) (input_PitchShiftPcw_Xend-1)));
					//				if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
					//					Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
					//				}
					//			}
					//		}
					//	}
					//}
					//else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){
					//	for (y2 = intPcwWid-1; y2 >= 0; y2--){
					//		if (x == output_PitchShiftChirpPcw_Xend){
					//			if (PCW_WIDTH_CHIRP != 0){
					//				intNowPCWwidthOffset = intPCWwidthChirp;
					//				Pnum[x][y2].Y -= intNowPCWwidthOffset;
					//				intPreviousPCWwidthOffset = intPCWwidthChirp;
					//			}
					//		}
					//		else{
					//			if (PCW_WIDTH_CHIRP != 0){
					//				intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x-output_PitchShiftChirpPcw_Xend) / (double) (input_PitchShiftPcw_Xend-1)));
					//				if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
					//					Pnum[x][y2].Y -= -(-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
					//				}
					//			}

					//		}
					//	}
					//}



					//}

					// 通常PCWあり
					//else{
					//	if (x != 0){
					//		Pnum[x][y].X = Pnum[x-1][y].X + intPitchX; 		// intPitchXだけ+X座標に置く
					//		Pnum[x][y].Y = Pnum[x-1][y].Y; 					// 同じY座標に置く
					//	}

					//	s_x3 = 0; 

					//	if (y == intPcwWid - LSPCW_ROW){

					//		// 入射 通常PCW
					//		if (intNormPcwPer != 0 && x < input_NormPcw_Xend + 1){
					//		}
					//		// 入射 チャープLSPCW
					//		if (x >= input_NormPcw_Xend && x < input_Chirp_Ls_Xend + 1){
					//			if (intChirp3rdLsPer == 0){
					//				s_x3 = 0;
					//			}
					//			else{
					//				s_x3 = INT_DIV (intSx3Per, intChirp3rdLsPer) * (x + 1 - input_NormPcw_Xend); 
					//			}
					//		}
					//		// LSPCW
					//		if (x >= input_Chirp_Ls_Xend && x < Lspcw_Xend){
					//			s_x3 = intSx3Per; 
					//		}
					//		// 出射 チャープLSPCW
					//		if (x >= Lspcw_Xend && x < output_Chirp_Ls_Xend){
					//			if (intChirp3rdLsPer == 0){
					//				s_x3 = 0;
					//			}
					//			else{
					//				s_x3 = INT_DIV (intSx3Per, intChirp3rdLsPer) * (output_Chirp_Ls_Xend - x); 
					//			}
					//		}
					//		// 入射 通常PCW
					//		if (x >= output_Chirp_Ls_Xend && x < output_PCW_Xend){
					//		}
					//	}
					//}

					/******************** 格子シフト量チャープ(2013/7/19) ********************/
					if (PITCH_SHIFT_MAX > PITCH && intNormPcwPer == 0 && intPitchShiftPcwPer + intPitchShiftChirpPcwPer != 0){
						//// 入射 格子定数変化PCW
						//if (x < input_PitchShiftPcw_Xend){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}
						//// 出射 格子定数変化PCW
						//else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}
						// 格子定数変化PCWとの出射接続部
						//if (x <= 0 + 1 || x >= intPcwPer - (1+1) ){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}

						//// 入射 格子定数変化PCW
						//if (x < input_PitchShiftPcw_Xend){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}
						//// 出射 格子定数変化PCW
						//else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}
						//// 格子定数変化PCWとの出射接続部
						//if (x == output_Chirp_Ls_Xend - 1){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}

						//// 入出射対象構造
						//if (x == output_PitchShiftChirpPcw_Xend - 1) {
						//	if (y % 2 == 1){
						//		continue;
						//	}
						//}
					}
					/******************** 格子シフト量チャープ(2013/7/19) ********************/

					if (SX2 == 0 && SX4 == 0)
					mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
					if (SX3 == 0 && SX4 == 0)
					mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
					if (SX2 == 0 && SX3 == 0)
					mcircle(Pnum[x][y].X + s_x4, Pnum[x][y].Y, z, 1);

					// 出射CORE細線導波路との接続部分の調整
					if(x == intPcwPer - 1){
						// 格子定数変化PCWあり
						/******************** 格子シフト量チャープ(2013/7/19) ********************/
						if (intNormPcwPer == 0 && intPitchShiftPcwPer + intPitchShiftChirpPcwPer != 0){						 
							//if (intNormPcwPer == 0 && intPitchShiftPcwPer != 0){
							/******************** 格子シフト量チャープ(2013/7/19) ********************/

							if(y % 2 == 1){
								//								//Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								//#if PCW_Air_Or_SiO2
								//								Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								//#else
								//								Pnum[x][y].X += intPitchX; 		// intPitchXだけ+X座標に置く
								//#endif

								if (Pnum[x][y].X > 0){
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
									if (PITCH_SHIFT_MAX == PITCH){
										Pnum[x][y].X += intPitchX;
									}
									else{
										Pnum[x][y].X += intPitchShiftX; 		// intPitchXだけ+X座標に置く
										//Pnum[x][y].X += (int) ((intPitchShiftX - intPitchX) / (double) (intPitchShiftChirpPcwPer-1)); 		// intPitchXだけ+X座標に置く
									}
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
								}
								else if (Pnum[x-1][y].X > 0){
									Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}
								else if (Pnum[x-2][y].X > 0){
									Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}

								int poo = 0;

								if (PCW_WIDTH_CHIRP != 0){
									// PCW_WIDTH_CHIRP用微調整
									//Pnum[x][y].Y -= 0;
									poo = INT_DIV (PCW_WIDTH_CHIRP, CELL_SIZE);
									Pnum[x][y].Y -= INT_DIV(poo, (intPitchShiftPcwPer-1));
								}


								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								if (PITCH_SHIFT_MAX != PITCH){
									Pnum[x][y].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchShiftY * (intPcwWid-1 - y) - intPCWwidthChirp;
								}
								/******************** 格子シフト量チャープ(2013/7/19) ********************/

								if (SX2 == 0 && SX4 == 0)
									mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
								if (SX3 == 0 && SX4 == 0)
									mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
								if (SX2 == 0 && SX3 == 0)
									mcircle(Pnum[x][y].X + s_x4, Pnum[x][y].Y, z, 1);
							}
						}
						else{
							if(y % 2 == 1){
								//								//Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								//#if PCW_Air_Or_SiO2
								//								Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								//#else
								//								Pnum[x][y].X += intPitchX; 		// intPitchXだけ+X座標に置く
								//#endif
								if (Pnum[x][y].X > 0){
									Pnum[x][y].X += intPitchX; 		// intPitchXだけ+X座標に置く
								}
								else if (Pnum[x-1][y].X > 0){
									Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}
								else if (Pnum[x-2][y].X > 0){
									Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}

								if (PCW_WIDTH_CHIRP != 0){
									// PCW_WIDTH_CHIRP用微調整
									//Pnum[x][y].Y -= 0;
									int poo;
									poo = INT_DIV (PCW_WIDTH_CHIRP, CELL_SIZE);
									Pnum[x][y].Y -= INT_DIV(poo, (intPitchShiftPcwPer-1));
								}

								if (SX2 == 0 && SX4 == 0)
									mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
								if (SX3 == 0 && SX4 == 0)
									mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
								if (SX2 == 0 && SX3 == 0)
									mcircle(Pnum[x][y].X + s_x4, Pnum[x][y].Y, z, 1);
							}
						}
					}

					// 線欠陥から数えた偶数列円孔の端を配置する場合に使用
					//if(y % 2 == 0){
					//	if(x == 0){
					//		mcircle(Pnum[x][y].X - intPitchShiftX, Pnum[x][y].Y, z, 1); 
					//	}
					//	if(x == intPcwPer - 1){
					//		mcircle(Pnum[x][y].X + intPitchShiftX, Pnum[x][y].Y, z, 1); 
					//	}
					//}

				}
				/****************************** LSPCW ******************************/

			}
		}
		// CORE細線導波路の位置を計算
		intWirePer2 = Pnum[intPcwPer-1][intPcwWid-1].X + intRadius;		// 出射COREスラブの開始点
		intWirePer3 = intWirePer2 + intWireLen2;						// 出射COREスラブの終了点

	}
	/****************************** フォトニック結晶 ******************************/



	/****************************** 入出射細線導波路 ******************************/
	int intPcwSislabOffset;

	// 全面スラブになっているので，細線以外の部分を空気に変更
	if (PCW_SiSLAB_OFFSET != 0){
		intPcwSislabOffset = INT_DIV(PCW_SiSLAB_OFFSET, CELL_SIZE);
	}
	else{
		intPcwSislabOffset = 0;
	}

	for (z = zmax_all - intSlabHeigPer; z < (zmax_all + 1); z++){
		for (y = 0; y < ymax_all - intWireWid_2; y++){

			// 入射
			if (PCW_SiSLAB_OFFSET != 0){
				//intPcwSislabOffset;
			}
			//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset; x++){		// 配列の引数に使用するので-1
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++){		// 配列の引数に使用するので-1
			//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
			//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}

			// 出射
			if (PCW_SiSLAB_OFFSET != 0){
				//intPcwSislabOffset--;		//入出射同一構造にするためのおまじない
			}
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++){
			//for (x = intWirePer2 + intPcwSislabOffset+8; x < intWirePer3-8; x++){	//nondoped
			//for (x = intWirePer2 + intPcwSislabOffset-6; x < intWirePer3; x++){		//doped
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}

		//田村
		//入射
		/*for (y = 0; y < ymax_all - intWireWid_2; y++){
			for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -8; x < intWireLen1 - 1 - intPcwSislabOffset - 7; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-1; y++){
			for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -6; x < intWireLen1 - 1 - intPcwSislabOffset - 5; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-2; y++){
			for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -4; x < intWireLen1 - 1 - intPcwSislabOffset - 3; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-4; y++){
			for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -2; x < intWireLen1 - 1 - intPcwSislabOffset - 2; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-6; y++){
			for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -1; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-6; y++){
			for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -1; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		//出射
		for (y = 0; y < ymax_all - intWireWid_2; y++){
			for (x = intWirePer2 + intPcwSislabOffset+6; x < intWirePer2 + intPcwSislabOffset+8; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-1; y++){
			for (x = intWirePer2 + intPcwSislabOffset+4; x < intWirePer2 + intPcwSislabOffset+6; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-2; y++){
			for (x = intWirePer2 + intPcwSislabOffset+2; x < intWirePer2 + intPcwSislabOffset+4; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-4; y++){
			for (x = intWirePer2 + intPcwSislabOffset+1; x < intWirePer2 + intPcwSislabOffset+2; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-6; y++){
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer2 + intPcwSislabOffset+1; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}
		for (y = 0; y < ymax_all - intWireWid_2-6; y++){
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer2 + intPcwSislabOffset+1; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD; 
				ALL_epsilonx[x][y][z] = epsilon2; 
				ALL_epsilony[x][y][z] = epsilon2; 
				ALL_epsilonz[x][y][z] = epsilon2; 
			}
		}*/
		//田村ここまで
	}
	/****************************** 入出射細線導波路 ******************************/



	/****************************** 対称境界部分の誘電率の設定 ******************************/

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

	/****************************** 対称境界部分の誘電率の設定 ******************************/



	/****************************** 各ノードにモデルを分割 ******************************/
	if(irank != IRANK_MAX){
		for(x = 0; x < xmax+1; x++){
			for(y = 0; y < ymax+1; y++){
				for(z = 0; z < zmax+1; z++){
					cell[x][y][z] = ALL_cell[irank*(xmax-1)+x][y][z]; 
					epsilonx[x][y][z] = ALL_epsilonx[irank*(xmax-1)+x][y][z]; 
					epsilony[x][y][z] = ALL_epsilony[irank*(xmax-1)+x][y][z]; 
					epsilonz[x][y][z] = ALL_epsilonz[irank*(xmax-1)+x][y][z]; 
				}
			}
		}
	}
	else{
		for(x = 0; x < xmax+1; x++){
			for(y = 0; y < ymax+1; y++){
				for(z = 0; z < zmax+1; z++){
					cell[x][y][z] = ALL_cell[irank*(xmax)+x][y][z]; 
					epsilonx[x][y][z] = ALL_epsilonx[irank*(xmax)+x][y][z]; 
					epsilony[x][y][z] = ALL_epsilony[irank*(xmax)+x][y][z]; 
					epsilonz[x][y][z] = ALL_epsilonz[irank*(xmax)+x][y][z]; 
				}
			}
		}
	}
	/****************************** 各ノードにモデルを分割 ******************************/


	/****************************** 共通パラメータの設定 ******************************/

	// 励振点，観測面の設定 (XMAXは "のりしろ" 部分を含めていることに注意)
	intExctPortNum = intExctLen / (XMAX - 1);
	intObseInPortNum = intObseLen1 / (XMAX - 1);
	intObseOutPortNum = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) / (XMAX - 1);
	if (NODE % 2 != 0){
		intObseCenPortNum = XMAX_ALL / 2 / (XMAX - 1);	// 奇数個並列計算のとき
	}
	else{
		intObseCenPortNum = XMAX_ALL / 2 / (XMAX - 1) - 1;	// 偶数個並列計算のとき
	}

	intExctLenPart = intExctLen % (XMAX - 1) - 1;		// 配列の引数に使用するので-1
	intObseLenPart1 = intObseLen1 % (XMAX - 1) - INT_DIV(intObseInter, 2) - 1;
	intObseLenPart2 = intObseLenPart1 + intObseInter;
	intObseLenPart3 = intObseLen1 % (XMAX - 1);
	intObseLenPart4 = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) % (XMAX - 1) - INT_DIV(intObseInter, 2);
	intObseLenPart5 = intObseLenPart4 + intObseInter; 
	intObseLenPart6 = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) % (XMAX - 1);
	//intObseLenPart7 = (XMAX_ALL / 2) % (XMAX - 1) - INT_DIV(intObseInter, 2);
	//intObseLenPart7 = (XMAX_ALL / 2) % (XMAX - 1);
	if (NODE % 2 != 0){
		intObseLenPart7 = (XMAX_ALL / 2) % (XMAX - 1) - 10;		// 奇数個並列計算のとき
	}
	else{
		intObseLenPart7 = (XMAX - 1) - 10;		// 偶数個並列計算のとき
	}
	//intObseLenPart8 = intObseLenPart7 + intObseInter; 
	//intObseLenPart9 = (XMAX_ALL / 2) % (XMAX - 1);

	/****************************** 共通パラメータの設定 ******************************/
}



//誘電率の割り当て
void set_epsilon(){

	//誘電率分布の出力(モデルの確認)
	int tag1 = 1; 

#if _FDTD

	/****************************** 計算実行時 ******************************/
	int node; 

	MPI_Status status; 

	//XY平面
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			epsilon_xy[x][y] = epsilonx[x][y][intSlabCen-1]; 
		}
	}
	if(irank == IRANK_MIN){
		for(x = 0; x<xmax; x++){
			for(y = 0; y < ymax+1; y++){
				fprintf (fpallepsilonx, "%e¥t", epsilon_xy[x][y]); 
				fprintf (fpepsilonx, "%e¥t", epsilon_xy[x][y]); 
			}
			fprintf (fpallepsilonx, "¥n"); 
			fprintf (fpepsilonx, "¥n"); 
		}
	}
	if(irank != IRANK_MIN){

		MPI_Send (&epsilon_xy[0][0], (XMAX+1)*(YMAX+1), MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD); 

		for(x = 1; x<xmax; x++){
			for(y = 0; y < ymax+1; y++){
				fprintf (fpepsilonx, "%e¥t", epsilon_xy[x][y]); 
			}
			fprintf (fpepsilonx, "¥n"); 
		}
	}
	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){

			MPI_Recv (&epsilon_xy[0][0], (XMAX+1)*(YMAX+1), MPI_DOUBLE, node, tag1, MPI_COMM_WORLD, &status); 

			if(node == IRANK_MAX){
				for(x = 1; x < xmax-1; x++){
					for(y = 0; y < ymax+1; y++){
						fprintf(fpallepsilonx, "%e¥t", epsilon_xy[x][y]); 
					}
					fprintf(fpallepsilonx, "¥n"); 
				}
			}
			else{
				for(x = 1; x < xmax; x++){
					for(y = 0; y < ymax+1; y++){
						fprintf (fpallepsilonx, "%e¥t", epsilon_xy[x][y]); 
					}
					fprintf (fpallepsilonx, "¥n"); 
				}
			}
		}
	}
	/****************************** 計算実行時 ******************************/
#else

	/****************************** モデル確認時 ******************************/
	char fname[40],dir_name[50];	//ファイル名格納変数	

	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			epsilon_xy[x][y] = epsilonx[x][y][intSlabCen-1]; 
		}
	}
	for(x = 0; x<xmax; x++){
		for(y = 0; y<ymax+1; y++){
			fprintf (fpallepsilonx, "%e¥t", epsilon_xy[x][y]); 
			fprintf (fpepsilonx, "%e¥t", epsilon_xy[x][y]); 
		}
		fprintf (fpallepsilonx, "¥n"); 
		fprintf (fpepsilonx, "¥n"); 
	}

#if _PROGRAM_TEST
#if	_MODEL_ALL_EPSILON
	//for(z = 0; z < zmax+1; z++){
	for(z = 0; z < zmax+1; z++){
		//for(z = intCladHeight1 - 1; z < zmax+1; z++){
		sprintf(fname,"/AllEpsilon_x_%d.txt",z);
		fpAllEpsilon = fopen(strcat(strcpy(dir_name, "FOR_TEST"),fname),"w");
		sprintf(fname,"/Epsilon_x_%d.txt",z);
		fpEpsilon = fopen(strcat(strcpy(dir_name, "FOR_TEST"),fname),"w");

		for(x=0; x < xmax; x++){
			for(y=0; y < ymax; y++){
				fprintf (fpAllEpsilon, "%e¥t", epsilonx[x][y][z]);
				fprintf (fpEpsilon, "%e¥t", epsilonx[x][y][z]);
			}
			fprintf (fpAllEpsilon, "¥n");
			fprintf (fpEpsilon, "¥n");
		}

		fclose(fpAllEpsilon);
		fclose(fpEpsilon);
	}
#endif
#else
	for(z = 0; z < zmax+1; z++){
		sprintf(fname,"/AllEpsilon_x_%d.txt",z);
		fpAllEpsilon = fopen(strcat(strcpy(dir_name, "FOR_TEST"),fname),"w");
		sprintf(fname,"/Epsilon_x_%d.txt",z);
		fpEpsilon = fopen(strcat(strcpy(dir_name, "FOR_TEST"),fname),"w");

		for(x=0; x < xmax; x++){
			for(y=0; y < ymax+1; y++){
				fprintf (fpAllEpsilon, "%e¥t", epsilonx[x][y][z]);
				fprintf (fpEpsilon, "%e¥t", epsilonx[x][y][z]);
			}
			fprintf (fpAllEpsilon, "¥n");
			fprintf (fpEpsilon, "¥n");
		}

		fclose(fpAllEpsilon);
		fclose(fpEpsilon);
	}
#endif

	/****************************** モデル確認時 ******************************/
#endif

	//YZ平面
	for(y = 0; y < ymax+1; y++){
		for(z = 0; z < zmax+1; z++){
			epsilon_yz[y][z] = epsilony[intObseLenPart1][y][z]; 
		}
	}
	for(y = 0; y < ymax+1; y++){
		for(z = 0; z < zmax+1; z++){
			fprintf(fpepsilony, "%e¥t", epsilon_yz[y][z]); 
		}
		fprintf(fpepsilony, "¥n"); 
	}

	//ZX平面 (Y:境界面)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			epsilon_zx[x][z] = epsilonz[x][ymax][z]; 
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			fprintf(fpepsilonz, "%e¥t", epsilon_zx[x][z]); 
		}
		fprintf(fpepsilonz, "¥n"); 
	}
	//ZX平面 (Y:中心)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			epsilon_zx2[x][z] = epsilonz[x][ymax/2][z]; 
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			fprintf(fpepsilonz2, "%e¥t", epsilon_zx2[x][z]); 
		}
		fprintf(fpepsilonz2, "¥n"); 
	}

	////YZ平面
	//for(y = 0; y < ymax+1; y++){
	//	for(z = 0; z <= zmax; z++){
	//		epsilon_yz[y][z] = epsilony[intObseLenPart1][y][z]; 
	//	}
	//}
	//for(y = 0; y < ymax+1; y++){
	//	for(z = 0; z <= zmax; z++){
	//		fprintf(fpepsilony, "%e¥t", epsilon_yz[y][z]); 
	//	}
	//	fprintf(fpepsilony, "¥n"); 
	//}

	////ZX平面 (Y:境界面)
	//for(x = 0; x < xmax; x++){
	//	for(z = 0; z < zmax+1; z++){
	//		epsilon_zx[x][z] = epsilonz[x][ymax][z]; 
	//	}
	//}
	//for(x = 0; x < xmax; x++){
	//	for(z = 0; z < zmax+1; z++){
	//		fprintf(fpepsilonz, "%e¥t", epsilon_zx[x][z]); 
	//	}
	//	fprintf(fpepsilonz, "¥n"); 
	//}
	////ZX平面 (Y:中心)
	//for(x = 0; x < xmax; x++){
	//	for(z = 0; z < zmax+1; z++){
	//		epsilon_zx2[x][z] = epsilonz[x][ymax/2][z]; 
	//	}
	//}
	//for(x = 0; x < xmax; x++){
	//	for(z = 0; z < zmax+1; z++){
	//		fprintf(fpepsilonz2, "%e¥t", epsilon_zx2[x][z]); 
	//	}
	//	fprintf(fpepsilonz2, "¥n"); 
	//}


	// ファイルポインタを閉じる
	if (irank == IRANK_MIN){
		fclose(fpallepsilonx);
	}
	fclose(fpepsilonx); 
	fclose(fpepsilony); 
	fclose(fpepsilonz); 


	//MPI_Barrier(MPI_COMM_WORLD); 			// 一度同期をとる

	//// メモリの開放
	//for(x = 0; x < XMAX+1; x++)	free(epsilon_xy[x]); 
	//for(y = 0; y < YMAX+1; y++)	free(epsilon_yz[y]); 
	//for(x = 0; x < XMAX+1; x++)	free(epsilon_zx[x]); 

	//free(epsilon_xy); 
	//free(epsilon_yz); 
	//free(epsilon_zx); 
}


// 励振関数
void source_func(){
	if(irank == intExctPortNum){
		// 励振点の設定
		int x = intExctLenPart;
		for(int y = ex_y_st; y < ex_y_ed; y++){
			for(int z = ex_z_st; z < ex_z_ed; z++){
#if _EXITATION_FUNC	// CW励振
				//面内正弦分布励振の場合 ←解析空間が偶数セルか奇数セルかで励振が異なるのでその都度注意

				// スラブ厚の半分のセル数:偶数 導波路幅の半分のセル数:偶数
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st - 1)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st - 1)) * sin(omega0*n*dt); 

				// スラブ厚の半分のセル数:奇数 導波路幅の半分のセル数:奇数
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt); 
				Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt); // 01
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st - 1)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt); // 02
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st - 1)) * sin(omega0*n*dt); // 03
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st - 1)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st - 1)) * sin(omega0*n*dt); // 04
#else	// Gaussian励振
				Hz[x][y][z] += 1000 * cos(omega0*(n-Npeak)*dt) * exp(-(SQ(sigma*dt*(n-Npeak))/2)) * cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st));
#endif
			}
		}
	}

	/****************************** 磁界の対称境界条件(4回対称) ******************************/

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


// モデルへの励振点，観測点の記録
void observation_func(){

	//int y, z; 

	if(irank == intObseInPortNum){ //入射

		for(int x = intObseLenPart1; x < intObseLenPart2; x++){
			/****************************** 観測面の修正(2013/8/8) ******************************/
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
					//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
					/****************************** 観測面の修正(2013/8/8) ******************************/
					if((y == YMAX-1) && (z == (intSlabCen-1))){		//活性層断面中央点の時を記録
						cell[x][y][z] = 4; 					//中央点確認用
					}
					cell[x][y][z] += OBSERVATION; 		//面確認用
				}
			}
		}
	}
	if (irank == intExctPortNum){
		int x;
		x = intExctLenPart;
		for(int y = ex_y_st; y <= ex_y_ed-1; y++){		//プラス1しているのはセル数の関係
			for(int z = ex_z_st; z <= ex_z_ed-1; z++){
				cell[x][y][z] += EXITATION; 		//励振面確認用
			}
		}
	}

	if(irank == intObseOutPortNum){ //出射 NODE 2
		for(int x = intObseLenPart4; x < intObseLenPart5; x++){
			/****************************** 観測面の修正(2013/8/8) ******************************/
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
					//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
					/****************************** 観測面の修正(2013/8/8) ******************************/

					if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
						cell[x][y][z] = 4; 					// 中央点確認用
					}
					cell[x][y][z] += OBSERVATION; 			// 面確認用
				}
			}
		}
	}

	if(irank == intObseCenPortNum){ //出射 NODE 2
		//for(int x = intObseLenPart7; x < intObseLenPart8; x++){
		int x = intObseLenPart7;
		/****************************** 観測面の修正(2013/8/8) ******************************/
		for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
			for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
				//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
				//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
				/****************************** 観測面の修正(2013/8/8) ******************************/

				if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
					cell[x][y][z] = 4; 					// 中央点確認用
				}
				cell[x][y][z] += OBSERVATION; 			// 面確認用
			}
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
void absorpt_bound_condition(){

	/****************************** 対称境界条件 ******************************/			

	// 4回対称
	for(z = 0; z < zmax+1; z++){
		Eyn1x00[ymax][z] = Eyn1x00[ymax-1][z];
		Eyn1x01[ymax][z] = Eyn1x01[ymax-1][z];
		Eyn1xm0[ymax][z] = Eyn1xm0[ymax-1][z];
		Eyn1xm1[ymax][z] = Eyn1xm1[ymax-1][z];
	}
	for(x = 0; x < xmax+1; x++){
		Eyn1z00[x][ymax] = Eyn1z00[x][ymax-1];
		Eyn1z01[x][ymax] = Eyn1z01[x][ymax-1];
	}
	for(x = 0; x < xmax; x++){
		Exn1z00[x][ymax+1] = -Exn1z00[x][ymax-1];
		Exn1z01[x][ymax+1] = -Exn1z01[x][ymax-1];
	}
	for(z = 0; z < zmax; z++){
		Ezn1x00[ymax+1][z] = -Ezn1x00[ymax-1][z];
		Ezn1x01[ymax+1][z] = -Ezn1x01[ymax-1][z];
		Ezn1xm0[ymax+1][z] = -Ezn1xm0[ymax-1][z];
		Ezn1xm1[ymax+1][z] = -Ezn1xm1[ymax-1][z];
	}

	// 8回対称
	for(y = 0; y < ymax+1; y++){
		Ezn1x00[y][zmax] = -Ezn1x00[y][zmax-1];
		Ezn1x01[y][zmax] = -Ezn1x01[y][zmax-1];
		Ezn1xm0[y][zmax] = -Ezn1xm0[y][zmax-1];
		Ezn1xm1[y][zmax] = -Ezn1xm1[y][zmax-1];
	}
	for(x = 0; x < xmax+1; x++){
		Ezn1y00[x][zmax] = -Ezn1y00[x][zmax-1];
		Ezn1y01[x][zmax] = -Ezn1y01[x][zmax-1];
	}
	for(x = 0; x < xmax; x++){
		Exn1y00[x][zmax+1] = Exn1y00[x][zmax-1];
		Exn1y01[x][zmax+1] = Exn1y01[x][zmax-1];
	}
	for(y = 0; y <= ymax-1; y++){
		Eyn1x00[y][zmax+1] = Eyn1x00[y][zmax-1];
		Eyn1x01[y][zmax+1] = Eyn1x01[y][zmax-1];
		Eyn1xm0[y][zmax+1] = Eyn1xm0[y][zmax-1];
		Eyn1xm1[y][zmax+1] = Eyn1xm1[y][zmax-1];
	}

    /****************************** 対称境界条件 ******************************/




	/****************************** Murの2次の吸収境界条件(Ex) ******************************/			

	double u1ax1, u2ax1,u3ax1, u4ax1;
	double u1bx1, u2bx1,u3bx1, u4bx1;
	double u2xa1;

	double velo_dt;

	if(irank != IRANK_MAX){	
		for(x = 1; x < xmax; x++){
			for(z = 1; z < zmax+1; z++){
				velo_dt = (C0 / sqrt(epsilonx[x][0][z]/epsilon0) ) * dt;

				u1ax1 = (velo_dt - dy) / (velo_dt + dy); 
				u2ax1 = (2.0 * dy) / (velo_dt + dy);
				u3ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy) );
				u4ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy) ); 

				Ex[x][0][z] = -Exn2y01[x][z]
				+u1ax1 * (Ex[x][1][z] + Exn2y00[x][z])
					+u2ax1 * (Exn1y00[x][z] + Exn1y01[x][z])
					+u3ax1 * (Exn1y00[x+1][z] - 2.0 * Exn1y00[x][z] + Exn1y00[x-1][z] + Exn1y01[x+1][z] - 2.0 * Exn1y01[x][z] + Exn1y01[x-1][z])
					+u4ax1 * (Exn1y00[x][z+1] - 2.0 * Exn1y00[x][z] + Exn1y00[x][z-1] + Exn1y01[x][z+1] - 2.0 * Exn1y01[x][z] + Exn1y01[x][z-1]); 
			}
		}
		for(x = 1; x < xmax; x++){
			for(y = 1; y < ymax+1; y++){
				velo_dt = (C0 / sqrt(epsilonx[x][y][0]/epsilon0) ) * dt;

				u1bx1 = (velo_dt - dz) / (velo_dt + dz);
				u2bx1 = (2.0 * dz) / (velo_dt + dz);
				u3bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) ); 
				u4bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) ); 

				Ex[x][y][0] = -Exn2z01[x][y]
				+u1bx1 * (Ex[x][y][1] + Exn2z00[x][y])
					+u2bx1 * (Exn1z00[x][y] + Exn1z01[x][y])
					+u3bx1 * (Exn1z00[x+1][y] - 2.0 * Exn1z00[x][y] + Exn1z00[x-1][y] + Exn1z01[x+1][y] - 2.0*Exn1z01[x][y] + Exn1z01[x-1][y])
					+u4bx1 * (Exn1z00[x][y+1] - 2.0 * Exn1z00[x][y] + Exn1z00[x][y-1] + Exn1z01[x][y+1] - 2.0*Exn1z01[x][y] + Exn1z01[x][y-1]); 
			}
		}
	}
	else{	
		for(x = 1; x < xmax-1; x++){
			for(z = 1; z < zmax+1; z++){
				velo_dt = (C0 / sqrt(epsilonx[x][0][z]/epsilon0) ) * dt;

				u1ax1 = (velo_dt - dy) / (velo_dt + dy); 
				u2ax1 = (2.0 * dy) /(velo_dt + dy); 
				u3ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy)); 			
				u4ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy)); 

				Ex[x][0][z] = -Exn2y01[x][z]
				+u1ax1 * (Ex[x][1][z] + Exn2y00[x][z])
					+u2ax1 * (Exn1y00[x][z] + Exn1y01[x][z])
					+u3ax1 * (Exn1y00[x+1][z] - 2.0 * Exn1y00[x][z] + Exn1y00[x-1][z] + Exn1y01[x+1][z] - 2.0 * Exn1y01[x][z] + Exn1y01[x-1][z])
					+u4ax1 * (Exn1y00[x][z+1] - 2.0 * Exn1y00[x][z] + Exn1y00[x][z-1] + Exn1y01[x][z+1] - 2.0 * Exn1y01[x][z] + Exn1y01[x][z-1]); 
			}
		}
		for(x = 1; x < xmax-1; x++){
			for(y = 1; y < ymax+1; y++){
				velo_dt = (C0 / sqrt(epsilonx[x][y][0]/epsilon0) ) * dt;

				u1bx1 = (velo_dt - dz) / (velo_dt + dz);
				u2bx1 = (2.0 * dz) / (velo_dt + dz);
				u3bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) ); 
				u4bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) ); 

				Ex[x][y][0] = -Exn2z01[x][y]
				+u1bx1 * (Ex[x][y][1] + Exn2z00[x][y])
					+u2bx1 * (Exn1z00[x][y] + Exn1z01[x][y])
					+u3bx1 * (Exn1z00[x+1][y] - 2.0 * Exn1z00[x][y] + Exn1z00[x-1][y] + Exn1z01[x+1][y] - 2.0 * Exn1z01[x][y] + Exn1z01[x-1][y])
					+u4bx1 * (Exn1z00[x][y+1] - 2.0 * Exn1z00[x][y] + Exn1z00[x][y-1] + Exn1z01[x][y+1] - 2.0 * Exn1z01[x][y] + Exn1z01[x][y-1]); 
			}
		}
	}

	/****************************** Murの2次の吸収境界条件(Ex) ******************************/			

	double u1xa, u1xc;
	double u2xa, u2xc;

	/****************************** Murの1次の吸収境界条件(Ex) ******************************/			

	for(y = 1; y < ymax+1; y++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonx[0][y][0]/epsilon0) ) * dt;
			u2xa = (velo_dt - dz) / (velo_dt + dz); 

			Ex[0][y][0] = Exn1z01[0][y] + u2xa * (Ex[0][y][1] - Exn1z00[0][y]); 
		}

		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilonx[xmax-1][y][0]/epsilon0) ) * dt;
			u2xc = (velo_dt - dz) / (velo_dt + dz); 

			Ex[xmax-1][y][0] = Exn1z01[xmax-1][y] + u2xc * (Ex[xmax-1][y][1] - Exn1z00[xmax-1][y]);
		}
	}

	for(z = 1; z < zmax; z++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonx[0][0][z]/epsilon0) ) * dt;
			u1xa = (velo_dt - dy) / (velo_dt + dy); 

			Ex[0][0][z] = Exn1y01[0][z] + u1xa * (Ex[0][1][z] - Exn1y00[0][z]); 
		}

		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilonx[xmax-1][0][z]/epsilon0) ) * dt;
			u1xc = (velo_dt - dy) / (velo_dt + dy);

			Ex[xmax-1][0][z] = Exn1y01[xmax-1][z] + u1xc * (Ex[xmax-1][1][z] - Exn1y00[xmax-1][z]); 
		}
	}


	// 辺(Murの1次の吸収境界条件) -- y平面とz平面からそれぞれ算出される値の平均値を取る
	if (irank != IRANK_MIN){
		for(x = 0; x < xmax; x++){
			velo_dt = (C0 / sqrt(epsilonx[x][0][0]/epsilon0) ) * dt;

			u2xa1 = (velo_dt - dx) / (velo_dt + dx);

			Ex[x][0][0] = 0.5 * (Exn1z01[x][0] + u2xa1 * (Ex[x][0][1] - Exn1z00[x][0])
				+ Exn1y01[x][0] + u2xa1 * (Ex[x][1][0] - Exn1y00[x][0]) );
		}
	}
	else{
		for(x = 1; x < xmax; x++){
			velo_dt = (C0 / sqrt(epsilonx[x][0][0]/epsilon0) ) * dt;

			u2xa1 = (velo_dt - dx) / (velo_dt + dx);

			Ex[x][0][0] = 0.5 * (Exn1z01[x][0] + u2xa1 * (Ex[x][0][1] - Exn1z00[x][0])
				+ Exn1y01[x][0] + u2xa1 * (Ex[x][1][0] - Exn1y00[x][0]) );
		}
	}

	/****************************** Murの1次の吸収境界条件(Ex) ******************************/

	double u1by1, u2by1, u3by1, u4by1;
	double u1cy1, u2cy1, u3cy1, u4cy1;
	double u1cy2, u2cy2, u3cy2, u4cy2;

	/****************************** Murの2次の吸収境界条件(Ey) ******************************/

	for(x = 1; x < xmax; x++){
		for(y = 1; y < ymax; y++){
			velo_dt = (C0 / sqrt(epsilony[x][y][0]/epsilon0) ) * dt;

			u1by1 = (velo_dt - dz) / (velo_dt + dz); 
			u2by1 = (2.0 * dz) / (velo_dt + dz); 
			u3by1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) ); 
			u4by1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) ); 

			Ey[x][y][0] = -Eyn2z01[x][y]
			+ u1by1 * (Ey[x][y][1] + Eyn2z00[x][y])
				+ u2by1 * (Eyn1z00[x][y] + Eyn1z01[x][y])
				+ u3by1 * (Eyn1z00[x+1][y] - 2.0 * Eyn1z00[x][y] + Eyn1z00[x-1][y] + Eyn1z01[x+1][y] - 2.0 * Eyn1z01[x][y] + Eyn1z01[x-1][y])
				+ u4by1 * (Eyn1z00[x][y+1] - 2.0 * Eyn1z00[x][y] + Eyn1z00[x][y-1] + Eyn1z01[x][y+1] - 2.0 * Eyn1z01[x][y] + Eyn1z01[x][y-1]); 
		}
	}
	if(irank == IRANK_MIN){
		for(y = 1; y < ymax; y++){
			for(z = 1; z < zmax+1; z++){
				velo_dt = (C0 / sqrt(epsilony[0][y][z]/epsilon0) ) * dt;

				u1cy1 = (velo_dt - dx) / (velo_dt + dx); 
				u2cy1 = (2.0 * dx) / (velo_dt + dx); 
				u3cy1 = (dx * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) ); 
				u4cy1 = (dx * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) ); 

				Ey[0][y][z] = -Eyn2x01[y][z]
				+u1cy1 * (Ey[1][y][z] + Eyn2x00[y][z])
					+u2cy1 * (Eyn1x00[y][z] + Eyn1x01[y][z])
					+u3cy1 * (Eyn1x00[y+1][z] - 2.0 * Eyn1x00[y][z] + Eyn1x00[y-1][z] + Eyn1x01[y+1][z] - 2.0 * Eyn1x01[y][z] + Eyn1x01[y-1][z])
					+u4cy1 * (Eyn1x00[y][z+1] - 2.0 * Eyn1x00[y][z] + Eyn1x00[y][z-1] + Eyn1x01[y][z+1] - 2.0 * Eyn1x01[y][z] + Eyn1x01[y][z-1]); 
			}
		}
	}			
	if(irank == IRANK_MAX){
		for(y = 1; y < ymax; y++){
			for(z = 1; z < zmax+1; z++){
				velo_dt = (C0 / sqrt(epsilony[xmax][y][z]/epsilon0) ) * dt;

				u1cy2 = (velo_dt - dx) / (velo_dt + dx); 
				u2cy2 = (2.0 * dx) / (velo_dt + dx); 
				u3cy2 = (dx * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) ); 
				u4cy2 = (dx * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) ); 

				Ey[xmax][y][z] =  -Eyn2xm1[y][z]
				+ u1cy2 * (Ey[xmax-1][y][z] + Eyn2xm0[y][z])
					+ u2cy2 * (Eyn1xm0[y][z] + Eyn1xm1[y][z])
					+ u3cy2 * (Eyn1xm0[y+1][z] - 2.0 * Eyn1xm0[y][z] + Eyn1xm0[y-1][z] + Eyn1xm1[y+1][z] - 2.0 * Eyn1xm1[y][z] + Eyn1xm1[y-1][z])
					+ u4cy2 * (Eyn1xm0[y][z+1] - 2.0 * Eyn1xm0[y][z] + Eyn1xm0[y][z-1] + Eyn1xm1[y][z+1] - 2.0 * Eyn1xm1[y][z] + Eyn1xm1[y][z-1]); 
			}
		}
	}

	/****************************** Murの2次の吸収境界条件(Ey) ******************************/


	double u2ya, u3ya, u3yb;
	double u2ya1;
	double u2yc1;

	/****************************** Murの1次の吸収境界条件(Ey) ******************************/

	for(x = 1; x < xmax; x++){
		velo_dt = (C0 / sqrt(epsilony[x][0][0]/epsilon0) ) * dt;
		u2ya = (velo_dt - dz) / (velo_dt + dz); 
		Ey[x][0][0] = Eyn1z01[x][0] + u2ya * (Ey[x][0][1] - Eyn1z00[x][0]); 
	}

	for(z = 1; z < zmax+1; z++){
		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilony[0][0][z]/epsilon0) ) * dt;

			u3ya = (velo_dt - dx) / (velo_dt + dx); 
			Ey[0][0][z] = Eyn1x01[0][z] + u3ya * (Ey[1][0][z] - Eyn1x00[0][z]); 
		}
		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilony[xmax][0][0]/epsilon0) ) * dt;

			u3yb = (velo_dt - dx) / (velo_dt + dx);
			Ey[xmax][0][z] = Eyn1xm1[0][z] + u3yb * (Ey[xmax-1][0][z] - Eyn1xm0[0][z]);
		}
	}

	// 辺(Murの1次の吸収境界条件) --x平面とz平面からそれぞれ算出される値の平均値を取る
	for(y = 0; y < ymax; y++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonx[0][y][0]/epsilon0) ) * dt;

			u2ya1 = (velo_dt - dz) / (velo_dt + dz); 
			Ey[0][y][0] = 0.5 * (Eyn1z01[0][y] + u2ya1 * (Ey[0][y][1] - Eyn1z00[0][y])
				+ Eyn1x01[y][0] + u2ya1 * (Ey[1][y][0] - Eyn1x00[y][0])); 
		}
		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilony[xmax][y][0]/epsilon0) ) * dt;

			u2yc1 = (velo_dt - dz) / (velo_dt + dz);
			Ey[xmax][y][0] = 0.5*(Eyn1z01[xmax][y] + u2yc1 * (Ey[xmax][y][1] - Eyn1z00[xmax][y])
				+ Eyn1xm1[y][0] + u2yc1 * (Ey[xmax-1][y][0] - Eyn1xm0[y][0]));
		}
	}
	/****************************** Murの1次の吸収境界条件(Ey) ******************************/

	double u1az1, u2az1, u3az1, u4az1;
	double u1cz1, u2cz1, u3cz1, u4cz1;
	double u1cz2, u2cz2, u3cz2, u4cz2;

	/****************************** Murの2次の吸収境界条件(Ez) ******************************/

	for(x = 1; x < xmax; x++){
		for(z = 1; z < zmax; z++){
			velo_dt = (C0 / sqrt(epsilonz[x][0][z] / epsilon0) ) * dt;

			u1az1 = (velo_dt - dy) / (velo_dt + dy); 
			u2az1 = (2.0 * dy) / (velo_dt + dy); 
			u3az1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy) ); 
			u4az1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy) ); 

			Ez[x][0][z] = -Ezn2y01[x][z]
			+ u1az1 * (Ez[x][1][z] + Ezn2y00[x][z])
				+ u2az1 * (Ezn1y00[x][z] + Ezn1y01[x][z])
				+ u3az1 * (Ezn1y00[x+1][z] - 2.0 * Ezn1y00[x][z] + Ezn1y00[x-1][z] + Ezn1y01[x+1][z] - 2.0 * Ezn1y01[x][z] + Ezn1y01[x-1][z])
				+ u4az1 * (Ezn1y00[x][z+1] - 2.0 * Ezn1y00[x][z] + Ezn1y00[x][z-1] + Ezn1y01[x][z+1] - 2.0 * Ezn1y01[x][z] + Ezn1y01[x][z-1]); 
		}
	}

	for(y = 1; y < ymax+1; y++){
		for(z = 1; z < zmax; z++){
			if(irank == IRANK_MIN){
				velo_dt = (C0 / sqrt(epsilonz[0][y][z] / epsilon0) ) * dt;

				u1cz1 = (velo_dt - dx) / (velo_dt + dx); 
				u2cz1 = (2.0 * dx) / (velo_dt + dx); 
				u3cz1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) ); 
				u4cz1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) ); 

				Ez[0][y][z] = -Ezn2x01[y][z]
				+ u1cz1 * (Ezn2x00[y][z] + Ez[1][y][z])
					+ u2cz1 * (Ezn1x00[y][z] + Ezn1x01[y][z])
					+ u3cz1 * (Ezn1x00[y+1][z] - 2.0 * Ezn1x00[y][z] + Ezn1x00[y-1][z] + Ezn1x01[y+1][z] - 2.0 * Ezn1x01[y][z] + Ezn1x01[y-1][z])
					+ u4cz1 * (Ezn1x00[y][z+1] - 2.0 * Ezn1x00[y][z] + Ezn1x00[y][z-1] + Ezn1x01[y][z+1] - 2.0 * Ezn1x01[y][z] + Ezn1x01[y][z-1]);
			}
			if(irank == IRANK_MAX){
				velo_dt = (C0 / sqrt(epsilonz[xmax][y][z] / epsilon0) ) * dt;

				u1cz2 = (velo_dt - dx) / (velo_dt + dx); 
				u2cz2 = (2.0 * dx) / (velo_dt + dx); 
				u3cz2 = (dy * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) ); 
				u4cz2 = (dz * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) ); 

				Ez[xmax][y][z] = -Ezn2xm1[y][z]
				+ u1cz2 * (Ezn2xm0[y][z] + Ez[xmax-1][y][z])
					+ u2cz2 * (Ezn1xm1[y][z] + Ezn1xm0[y][z])
					+ u3cz2 * (Ezn1xm1[y+1][z] - 2.0 * Ezn1xm1[y][z] + Ezn1xm1[y-1][z] + Ezn1xm0[y+1][z] - 2.0 * Ezn1xm0[y][z] + Ezn1xm0[y-1][z])
					+ u4cz2 * (Ezn1xm1[y][z+1] - 2.0 * Ezn1xm1[y][z] + Ezn1xm1[y][z-1] + Ezn1xm0[y][z+1] - 2.0 * Ezn1xm0[y][z] + Ezn1xm0[y][z-1]); 
			}
		}
	}

	/****************************** Murの2次の吸収境界条件(Ez) ******************************/

	double u1za, u3za, u3zb;
	double u1za1, u1zb1;

	/****************************** Murの1次の吸収境界条件(Ez) ******************************/
	for(x = 1; x < xmax; x++){
		velo_dt = (C0 / sqrt(epsilonz[x][0][0] / epsilon0) ) * dt;
		u1za = (velo_dt - dy) / (velo_dt + dy); 

		Ez[x][0][0] = Ezn1y01[x][0] + u1za * (Ez[x][1][0] - Ezn1y00[x][0]); 
	}

	for(y = 1; y < ymax+1; y++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonz[0][y][0] / epsilon0) ) * dt;
			u3za = (velo_dt - dx) / (velo_dt + dx); 

			Ez[0][y][0] = Ezn1x01[y][0] + u3za * (Ez[1][y][0] - Ezn1x00[y][0]); 
		}
		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilonz[xmax][y][0] / epsilon0) ) * dt;
			u3zb = (velo_dt - dx) / (velo_dt + dx);

			Ez[xmax][y][0] = Ezn1xm1[y][0] + u3zb * (Ez[xmax-1][y][0] - Ezn1xm0[y][0]); 
		}
	}

	// 辺(Murの1次の吸収境界条件) --x平面とy平面からそれぞれ算出される値の平均値を取る
	for(z = 0; z < zmax+1; z++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonz[0][0][z] / epsilon0) ) * dt;
			u1za1 = (velo_dt - dy) / (velo_dt + dy);

			Ez[0][0][z] = 0.5 * (Ezn1y01[0][z] + u1za1 * (Ez[0][1][z] - Ezn1y00[0][z])
				+ Ezn1x01[0][z] + u1za1 * (Ez[1][0][z] - Ezn1x00[0][z]) ); 
		}
		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilonz[xmax][0][z] / epsilon0) ) * dt;
			u1zb1 = (velo_dt - dy) / (velo_dt + dy); 

			Ez[xmax][0][z] = 0.5 * (Ezn1y01[xmax][z] + u1zb1 * (Ez[xmax][1][z] - Ezn1y00[xmax][z])
				+ Ezn1xm1[0][z] + u1zb1 * (Ez[xmax-1][0][z] - Ezn1xm0[0][z])); 
		}
	}
	/****************************** Murの1次の吸収境界条件(Ez) ******************************/

}



/*電界の保存*/
void saving_electric_field(){

	// Ex
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax+1; z++){
			Exn2y00[x][z] = Exn1y00[x][z]; 
			Exn1y00[x][z] = Ex[x][0][z]; 
			Exn2y01[x][z] = Exn1y01[x][z]; 
			Exn1y01[x][z] = Ex[x][1][z]; 
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			Exn2z00[x][y] = Exn1z00[x][y]; 
			Exn1z00[x][y] = Ex[x][y][0]; 
			Exn2z01[x][y] = Exn1z01[x][y]; 
			Exn1z01[x][y] = Ex[x][y][1]; 
		}
	}

	// Ey
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Eyn2z00[x][y] = Eyn1z00[x][y]; 
			Eyn1z00[x][y] = Ey[x][y][0]; 
			Eyn2z01[x][y] = Eyn1z01[x][y]; 
			Eyn1z01[x][y] = Ey[x][y][1]; 
		}
	}
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax+1; z++){
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
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax; z++){
			Ezn2y00[x][z] = Ezn1y00[x][z]; 
			Ezn1y00[x][z] = Ez[x][0][z]; 
			Ezn2y01[x][z] = Ezn1y01[x][z]; 
			Ezn1y01[x][z] = Ez[x][1][z]; 
		}
	}
	for(y = 0; y < ymax+1; y++){
		for(z = 0; z < zmax; z++){
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


//モデルの出力
void output_model(){

	//char cell_xy[XMAX][YMAX]; 
	int tag2 = 2; 
	int x, y, z;

#if _FDTD
	/****************************** 計算実行時 ******************************/
	int node; 

	MPI_Status status; 

	z = intSlabCen - 1;

	// XY平面
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			cell_xy[x][y] = cell[x][y][z]; 
		}
	}

	if(irank == IRANK_MIN){
		for(x = 0; x < xmax; x++){
			for(y = 0; y < ymax; y++){
				fprintf(allmodel_xy, "%d¥t", cell_xy[x][y]); 
				fprintf(model_xy, "%d¥t", cell_xy[x][y]); 
			}
			fprintf(allmodel_xy, "¥n"); 
			fprintf(model_xy, "¥n"); 
		}
		fclose(model_xy);
	}

	// それぞれ分割部のモデル
	if(irank != IRANK_MIN){
		MPI_Send(&cell_xy[0][0], (xmax)*(ymax), MPI_INT, 0, tag2, MPI_COMM_WORLD); 
		for(x = 1; x < xmax; x++){
			for(y = 0; y < ymax; y++){
				fprintf(model_xy, "%d¥t", cell_xy[x][y]); 
			}
			fprintf(model_xy, "¥n"); 
		}

		fclose(model_xy); 
	}

	// 全体モデル生成
	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){
			if(node == IRANK_MAX){

				// 最終段は"のりしろ"が無いので，x方向を-1する
				MPI_Recv(&cell_xy[0][0], (xmax-1)*(ymax), MPI_INT, node, tag2, MPI_COMM_WORLD, &status); 

				for(x = 1; x < xmax-1; x++){
					for(y = 0; y < ymax; y++){
						fprintf(allmodel_xy, "%d¥t", cell_xy[x][y]); 
					}
					fprintf(allmodel_xy, "¥n"); 
				}
			}
			else{

				MPI_Recv(&cell_xy[0][0], (xmax)*(ymax), MPI_INT, node, tag2, MPI_COMM_WORLD, &status); 

				for(x = 1; x < xmax; x++){
					for(y = 0; y < ymax; y++){
						fprintf(allmodel_xy, "%d¥t", cell_xy[x][y]); 
					}
					fprintf(allmodel_xy, "¥n"); 
				}
			}
		}
		fclose(allmodel_xy); 
	}


	// YZ平面
	if(irank == intObseInPortNum){ // 入射
		x = intObseLenPart1;
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				cell_yz[y][z] = cell[x][y][z]; 
			}
		}
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				fprintf(allmodel_yz1, "%d¥t", cell_yz[y][z]); 
			}
			fprintf(allmodel_yz1, "¥n");  
		}
		fclose(allmodel_yz1);
	}

	if(irank == intObseOutPortNum){ // 出射
		x = intObseLenPart4;
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				cell_yz[y][z] = cell[x][y][z]; 
			}
		}
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				fprintf(allmodel_yz4, "%d¥t", cell_yz[y][z]); 
			}
			fprintf(allmodel_yz4, "¥n");  
		}
		fclose(allmodel_yz4);
	}

	if(irank == intObseCenPortNum){			// 中央
		x = intObseLenPart7;
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				cell_yz[y][z] = cell[x][y][z]; 
			}
		}
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				fprintf(allmodel_yz7, "%d¥t", cell_yz[y][z]); 
			}
			fprintf(allmodel_yz7, "¥n");  
		}
		fclose(allmodel_yz7);
	}
	/****************************** 計算実行時 ******************************/

#else	
	/****************************** モデル確認時 ******************************/
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			cell_xy[x][y] = cell[x][y][intSlabCen-1]; 
		}
	}

	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			fprintf(allmodel_xy, "%d¥t", cell_xy[x][y]); 
			fprintf(model_xy, "%d¥t", cell_xy[x][y]); 
		}
		fprintf(allmodel_xy, "¥n"); 
		fprintf(model_xy, "¥n"); 
	}
	fclose(model_xy); 

	x = intObseLenPart1;
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax; z++){
			cell_yz[y][z] = cell[x][y][z]; 
		}
	}
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax; z++){
			fprintf(allmodel_yz1, "%d¥t", cell_yz[y][z]); 
		}
		fprintf(allmodel_yz1, "¥n");  
	}
	fclose(allmodel_yz1);

	if(irank == intObseOutPortNum){ // 出射
		x = intObseLenPart4;
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				cell_yz[y][z] = cell[x][y][z]; 
			}
		}
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				fprintf(allmodel_yz4, "%d¥t", cell_yz[y][z]); 
			}
			fprintf(allmodel_yz4, "¥n");  
		}
		fclose(allmodel_yz4);
	}

	if(irank == intObseCenPortNum){			// 中央
		x = intObseLenPart7;
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				cell_yz[y][z] = cell[x][y][z]; 
			}
		}
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				fprintf(allmodel_yz7, "%d¥t", cell_yz[y][z]); 
			}
			fprintf(allmodel_yz7, "¥n");  
		}
		fclose(allmodel_yz7);
	}
	/****************************** モデル確認時 ******************************/
#endif

}


void output_field_write(char *dir_name_def){

	char fname[40], dir_name[50]; 	//ファイル名格納変数	
	int node; 
	int tag3 = 3; 
	int pi1, pj1, pk1; 
	MPI_Status status; 
	FILE *HZ1; 
	//FILE *HZ1, *HZ2; 

	pi1 = x_cen; 
	pj1 = y_cen; 
	pk1 = z_cen; 

	printf("n = %d¥n", n); 

	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			field_xy[x][y] = Hz[x][y][ex_z_ed-1]; 	//全てのノードで電磁界成分を2次元配列に格納する．
		}
	}

	// モデル出力ファイルポインタの初期化
	if(irank == IRANK_MIN){
		sprintf(fname, "/Field_Hz_XY_%d_01.txt", n);
		HZ1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname), "w"); 
		for(x = 0; x < xmax; x++){
			for(y = 0; y < ymax; y++){
				fprintf(HZ1, "%e¥t", field_xy[x][y]); 
			}
			fprintf(HZ1, "¥n"); 
		}
	}

	// モデルをホストに送信
	else{
		if(irank != IRANK_MAX){
			MPI_Send(&field_xy[0][0], (xmax)*(ymax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 		// ノード0以外のノードがノード0に電磁界成分を送る．
		}
		if(irank == IRANK_MAX){
			MPI_Send(&field_xy[0][0], (xmax-1)*(ymax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 	// ノード0以外のノードがノード0に電磁界成分を送る．
		}
	}

	// 受信したモデルから全モデルを作成
	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){		// ノード0がノード1から順にデータを受け取り出力していく．
			if(node == IRANK_MAX){					// ノードisize-1のみ1セル小さく設定しているため条件文で分岐
				MPI_Recv(&field_xy[0][0], (xmax-1)*(ymax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status); 
				for(x = 1; x < xmax-1; x++){
					for(y = 0; y < ymax; y++){
						fprintf(HZ1, "%e¥t", field_xy[x][y]); 
					}
					fprintf(HZ1, "¥n"); 
				}
			}
			else{
				MPI_Recv(&field_xy[0][0], (xmax)*(ymax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status); 
				for(x = 1; x < xmax; x++){
					for(y = 0; y < ymax; y++){
						fprintf(HZ1, "%e¥t", field_xy[x][y]); 
					}
					fprintf(HZ1, "¥n"); 
				}
			}
		}

		// ファイルポインタを閉じる
		fclose(HZ1); 		
	}


	// YZ平面の電界分布の出力
	int x;
	double E_yz;
	FILE *EYZ1, *EYZ2, *EYZ3;
	char fname2[40], fname3[40], fname4[40];

	if(irank == intObseInPortNum){ //入射
		x = intObseLenPart1;
		sprintf(fname2, "/Field_E_YZ_%d_01.txt", n);
		EYZ1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname2), "w"); 

		for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
			for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				fprintf(EYZ1, "%e¥t", E_yz);
			}
			fprintf(EYZ1, "¥n");
		}
		fclose(EYZ1);
	}

	if(irank == intObseOutPortNum){			// 出射
		x = intObseLenPart4;
		sprintf(fname3, "/Field_E_YZ_%d_04.txt", n);
		EYZ2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname3), "w"); 

		for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
			for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				fprintf(EYZ2, "%e¥t", E_yz);
			}
			fprintf(EYZ2, "¥n");
		}
		fclose(EYZ2);
	}

	if(irank == intObseCenPortNum){			// 出射
		x = intObseLenPart7;
		sprintf(fname4, "/Field_E_YZ_%d_07.txt", n);
		EYZ3 = fopen(strcat(strcpy(dir_name, dir_name_def), fname4), "w"); 

		for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
			for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				fprintf(EYZ3, "%e¥t", E_yz);
			}
			fprintf(EYZ3, "¥n");
		}
		fclose(EYZ3);
	}
}

//ファイル出力
void output_field(char *dir_name_def){

	//double field_xy[XMAX][YMAX]; 	// Hz-field のファイル出力 (面垂直方向の磁界成分)

	if(n <= Nmax - Fcut){
		// 動作確認のためのファイル出力
		if(n == Ncheck){
			output_field_write (dir_name_def); 
		}

		// 定期的なファイル出力
		if(n % Ncutfield == 0){
			output_field_write (dir_name_def); 
		}
	}
	if((n >= Nmax - Fcut) && (n <= Nmax)){

		// 安定点でのファイル出力
		if(n % Ncutfield2 == 0){
			output_field_write (dir_name_def); 
		}
	}
}


void calc_poynting_powerHz(){				// 透過スペクトル計算用の観測面中央でのHzの出力

	//入力部分での poynting power (x方向)
	//double Eyin, Ezin, Hyin, Hzin; 		// 各成分の積算保存変数
	//double EyHz, EzHy; 					// 積分値を掛け算を保存する変数とセル中央の値

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
}



void calc_poynting_power(){				// 吉田作成パワー評価プログラム

	//入力部分での poynting power (x方向)
	double Eyin, Ezin, Hyin, Hzin; 		// 各成分の積算保存変数
	double EyHz, EzHy; 					// 積分値を掛け算を保存する変数とセル中央の値

	double pmax_01 = 0;		// 入力パワーの最大値を記録する変数
	double pmax_03 = 0;		// 出力パワーの最大値を記録
	double pmin_01 = 0;		// 出力パワーの最小値を記録
	double pmin_03 = 0;		// 出力パワーの最小値を記録

	int x = 0; 

	if(irank == intObseInPortNum){ //入射

		if(n == Nmax - Tcut){ //ある時刻までの入出力パワーの最大値と最小値を保持する変数
			powermax_in = 0.0; 
			powermin_in = 0.0; 
			powermax_out = 0.0; 
			powermin_out = 0.0; 
		}

		for(x = intObseLenPart1; x < intObseLenPart2; x++){

			// 初期化
			Eyin = 0; 	Ezin = 0; 	Hyin = 0; 	Hzin = 0; 	EyHz = 0; 	EzHy = 0; 		

			/****************************** 観測面の修正(2013/8/8) ******************************/
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
					//	for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
					/****************************** 観測面の修正(2013/8/8) ******************************/


					Eyin = 0.25 * (Ey[x][y][z] + Ey[x+1][y][z] + Ey[x][y][z+1] + Ey[x+1][y][z+1]); 
					Hzin = 0.50 * (Hz[x][y][z] + Hz[x][y][z+1]); 
					Ezin = 0.25 * (Ez[x][y][z] + Ez[x+1][y][z] + Ez[x][y+1][z] + Ez[x+1][y+1][z]); 
					Hyin = 0.50 * (Hy[x][y][z] + Hy[x][y+1][z]); 

					EyHz += Eyin * Hzin; 
					EzHy += Ezin * Hyin;
				}
			}

			//ポインティングパワーをファイルに保存
			if(x == intObseLenPart3){
				fprintf(fppoynt1, "%e¥n", EyHz - EzHy);
			}
			//if(x == intObseLenPartHz1){
			//	fprintf(fppoynt1h, "%e¥n", EyHz - EzHy); 
			//}

			if(n >= Nmax - Tcut){
				if(pmax_01 < EyHz - EzHy){
					pmax_01 = EyHz - EzHy;
				}

				if(pmin_01 > EyHz - EzHy){
					pmin_01 = EyHz - EzHy;
				}
			}
		}

		if(n >= Nmax - Tcut){
			if(powermax_in < pmax_01){
				powermax_in = pmax_01;
			}
			if(powermin_in > pmin_01){
				powermin_in = pmin_01;
			}
		}
		if(n == Nmax){
			fprintf(avpoynt1, "%s", "Input_Power¥t"); 
			fprintf(avpoynt1, "%e¥n", powermax_in); 
			fprintf(avpoynt1, "%s", "Reflection¥t"); 
			fprintf(avpoynt1, "%e¥n", powermin_in); 
			fprintf(avpoynt1, "%s", "Degree_of_Reflection¥t"); 
			fprintf(avpoynt1, "%e¥n", powermin_in/powermax_in); 
		}
	}

	if(irank == intObseOutPortNum){			// 出射

		if(n == Nmax - Tcut){ //ある時刻までの入出力パワーの最大値と最小値を保持する変数
			powermax_in = 0.0; 
			powermin_in = 0.0; 
			powermax_out = 0.0; 
			powermin_out = 0.0; 
		}

		for(x = intObseLenPart4; x < intObseLenPart5; x++){

			//初期化
			Eyin = 0; 	Ezin = 0; 	Hyin = 0; 	Hzin = 0; 	EyHz = 0; 	EzHy = 0; 		

			/****************************** 観測面の修正(2013/8/8) ******************************/
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
					//	for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
					/****************************** 観測面の修正(2013/8/8) ******************************/

					Eyin = 0.25 * (Ey[x][y][z] + Ey[x+1][y][z] + Ey[x][y][z+1] + Ey[x+1][y][z+1]); 
					Hzin = 0.50 * (Hz[x][y][z] + Hz[x][y][z+1]); 
					Ezin = 0.25 * (Ez[x][y][z] + Ez[x+1][y][z] + Ez[x][y+1][z] + Ez[x+1][y+1][z]); 
					Hyin = 0.50 * (Hy[x][y][z] + Hy[x][y+1][z]); 

					//if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
					//	cell[x][y][z] = 4; 					// 中央点確認用
					//}

					EyHz += Eyin * Hzin; 
					EzHy += Ezin * Hyin;					// ポインティングパワーを各セル毎に足す
				}
			}

			if(x == intObseLenPart6){
				fprintf(fppoynt5, "%e¥n", EyHz - EzHy); 		// ポインティングパワーをファイルに保存
			}
			//if(x == intObseLenPartHz5){
			//	fprintf(fppoynt5h, "%e¥n", EyHz - EzHy);
			//}

			if(n >= Nmax - Tcut){
				if(pmax_03 < EyHz - EzHy){
					pmax_03 = EyHz - EzHy; 
				}
				if(pmin_03 > EyHz - EzHy){
					pmin_03 = EyHz - EzHy; 
				}
			}
		}
		if(n >= Nmax - Tcut){
			if(powermax_out < pmax_03){
				powermax_out = pmax_03; 
			}

			if(powermin_out > pmin_03){
				powermin_out = pmin_03; 
			}
		}
		if(n == Nmax){
			fprintf(avpoynt5, "%s", "Transmission¥t"); 
			fprintf(avpoynt5, "%e¥n", powermax_out); 
			fprintf(avpoynt5, "%s", "Reflection¥t"); 
			fprintf(avpoynt5, "%e¥n¥n", powermin_out); 
		}
	}
}



void mcircle(int x_circ, int y_circ, int z_circ, int type){

	double R; 

	//半径セル数の計算
	if(type == 1)	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています
	else if(type == 2)	R = ((dblRadius2*1.0e10)/(dx*1.0e10)); 
	else if(type == 3)	R = ((dblRadius3*1.0e10)/(dx*1.0e10)); 
	else if(type == 5)	R = ((dblRadius5*1.0e10)/(dx*1.0e10)); 
	else if(type == 6)	R = ((dblRadius6*1.0e10)/(dx*1.0e10)); 
	else if(type == 7)	R = ((dblRadius7*1.0e10)/(dx*1.0e10)); 
	else if(type == 8)	R = ((dblRadius8*1.0e10)/(dx*1.0e10)); 
	else			R = ((dblRadius4*1.0e10)/(dx*1.0e10)); 

	rightquartercircle1(x_circ, y_circ, z_circ, type, R); 
	leftquartercircle1(x_circ-1, y_circ, z_circ, type, R); 
	rightquartercircle2(x_circ, y_circ-1, z_circ, type, R); 
	leftquartercircle2(x_circ-1, y_circ-1, z_circ, type, R); 
}


void halfcircle(int x_circ, int y_circ, int z_circ, int type){

	double R; 

	//半径セル数の計算
	if(type == 1)	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています
	else if(type == 2)	R = ((dblRadius2*1.0e10)/(dx*1.0e10)); 
	else if(type == 3)	R = ((dblRadius3*1.0e10)/(dx*1.0e10)); 
	else if(type == 5)	R = ((dblRadius5*1.0e10)/(dx*1.0e10)); 
	else if(type == 6)	R = ((dblRadius6*1.0e10)/(dx*1.0e10)); 
	else if(type == 7)	R = ((dblRadius7*1.0e10)/(dx*1.0e10)); 
	else if(type == 8)	R = ((dblRadius8*1.0e10)/(dx*1.0e10)); 
	else			R = ((dblRadius4*1.0e10)/(dx*1.0e10)); 

	rightquartercircle2(x_circ, y_circ-1, z_circ, type, R); 
	leftquartercircle2(x_circ-1, y_circ-1, z_circ, type, R); 
}


void rightquartercircle1(int x_circ, int y_circ, int z_circ, int type, double R){

	int x, y, Ie, Je; 
	double r; 

	Ie = (int) (x_circ+R-1); 
	Je = (int) (y_circ+R-1); 
	for(x = x_circ; x <= Ie; x++){
		for(y = y_circ; y <= Je; y++){
			r = sqrt(double((x-x_circ+1) * (x-x_circ+1) + (y-y_circ+1) * (y-y_circ+1)) ) - 0.5; 
			if(r <= R){
				if(type == 1 || type == 3|| type == 5 || type == 6 || type == 7 || type == 8){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX; 
					ALL_epsilonx[x][y][z_circ] = epsilon2; 
					ALL_epsilony[x][y][z_circ] = epsilon2; 
					ALL_epsilonz[x][y][z_circ] = epsilon2; 
				}else if(type == 2){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2; 
					ALL_epsilonx[x][y][z_circ] = epsilon2; 
					ALL_epsilony[x][y][z_circ] = epsilon2; 
					ALL_epsilonz[x][y][z_circ] = epsilon2; 
				}else{
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3; 
				}
			}
		}
	}
}


void leftquartercircle1(int x_circ, int y_circ, int z_circ, int type, double R){

	int x, y, Ie, Je; 
	double r; 

	Ie = (int) (x_circ-R+1); 
	Je = (int) (y_circ+R-1); 
	for(x = x_circ; x >= Ie; x--){
		for(y = y_circ; y <= Je; y++){
			r = sqrt(double((x-x_circ-1) * (x-x_circ-1) + (y-y_circ+1) * (y-y_circ+1)) ) - 0.5; 
			if(r <= R){
				if(type == 1 || type == 3|| type == 5 || type == 6 || type == 7 || type == 8){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX; 
					ALL_epsilonx[x][y][z_circ] = epsilon2; 
					ALL_epsilony[x][y][z_circ] = epsilon2; 
					ALL_epsilonz[x][y][z_circ] = epsilon2; 
				}else if(type == 2){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2; 
					ALL_epsilonx[x][y][z_circ] = epsilon2; 
					ALL_epsilony[x][y][z_circ] = epsilon2; 
					ALL_epsilonz[x][y][z_circ] = epsilon2; 
				}else{
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3; 
				}
			}
		}
	}
}


void rightquartercircle2(int x_circ, int y_circ, int z_circ, int type, double R){

	int x, y, Ie, Je; 
	double r; 

	Ie = (int) (x_circ+R-1); 
	Je = (int) (y_circ-R+1); 
	for(x = x_circ; x <= Ie; x++){
		for(y = y_circ; y >= Je; y--){
			r = sqrt(double((x-x_circ+1) * (x-x_circ+1) + (y-y_circ-1) * (y-y_circ-1)) ) - 0.5; 
			if(r <= R){
				if(type == 1 || type == 3|| type == 5 || type == 6 || type == 7 || type == 8){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX; 
					ALL_epsilonx[x][y][z_circ] = epsilon2; 
					ALL_epsilony[x][y][z_circ] = epsilon2; 
					ALL_epsilonz[x][y][z_circ] = epsilon2; 
				}else if(type == 2){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2; 
					ALL_epsilonx[x][y][z_circ] = epsilon2; 
					ALL_epsilony[x][y][z_circ] = epsilon2; 
					ALL_epsilonz[x][y][z_circ] = epsilon2; 
				}else{
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3; 
				}
			}
		}
	}
}

void leftquartercircle2(int x_circ, int y_circ, int z_circ, int type, double R){

	int x, y, Ie, Je; 
	double r; 

	Ie = (int) (x_circ-R+1); 
	Je = (int) (y_circ-R+1); 
	for(x = x_circ; x >= Ie; x--){
		for(y = y_circ; y >= Je; y--){
			r = sqrt(double((x-x_circ-1) * (x-x_circ-1) + (y-y_circ-1) * (y-y_circ-1)) ) - 0.5; 
			if(r <= R){
				if(type == 1 || type == 3 || type == 5 || type == 6 || type == 7 || type == 8){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX; 
					ALL_epsilonx[x][y][z_circ] = epsilon2; 
					ALL_epsilony[x][y][z_circ] = epsilon2; 
					ALL_epsilonz[x][y][z_circ] = epsilon2; 
				}else if(type == 2){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2; 
					ALL_epsilonx[x][y][z_circ] = epsilon2; 
					ALL_epsilony[x][y][z_circ] = epsilon2; 
					ALL_epsilonz[x][y][z_circ] = epsilon2; 
				}else{
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3; 
				}
			}
		}
	}
}


