/*
3 Dimentional Finite Difference Time Domain Calculation with Periodic Boundary Condition
ver. 5.15
programed by IWAI Takeshi(1999-2002)
original PML programed by Super Doctor SAKAI Atsushi(1996-2002)
imaginary part of PML programed by
IWAI Takeshi(1999-2002)

2001/5/1 ver. 1.0 プログラムの立ち上げ(2 Dimensional version)
2001/5/16 ver. 1.1 プログラム修正
2001/5/18 ver. 1.2 プログラム修正(波数ベクトルの修正)
2001/5/21 ver. 1.3 プログラム完成
2001/5/23 ver. 1.4 プログラム修正(計算範囲の修正)
2001/5/23 ver. 2.0 プログラム完成(導波路計算導入)
2001/5/24 ver. 2.1 プログラム修正(フィールド吐き出し成分の追加)
2001/5/28 ver. 2.1 プログラム修正(File Closeの設定追加)
2001/5/29 ver. 2.2 プログラム修正(波数ベクトルの修正)
2001/5/30 ver. 2.3 プログラム修正(波数ベクトルG-Pの追加)

2001/6/12 ver. 1.0 プログラムの立ち上げ(3 Dimensional version)
2001/7/9  ver. 1.1 プログラム完成

2001/7/17 ver. 3.0 x方向PML吸収境界導入プログラムの立ち上げ
2001/8/2  ver. 3.1 バグ取り修正
2001/10/24 ver. 3.2 コンパイル時に警告が出るバグを修正 細かいコメント追加
2001/10/24 ver. 3.3 クラッド層の追加
*/
/*
arranged by Eiichi Mizuta(2003-)

2003/10/1 ver. 4.0b モデルをDBR型フォトニック結晶導波路へ変更
2003/10/7 ver. 4.1b 電界パラメータを配列で定義
2003/10/8 ver. 4.2b PML層を隣接するモデルと同化 モデルの吐き出し，フィールドの吐き出し，電磁界パルスの吐き出しの変更 時間の吐き出しの追加
2003/10/9 ver. 4.3b パルスピーク位置を変更 吸収境界条件のEx-PMLの範囲の修正 電磁界計算のcondition2の範囲の修正
2003/10/10 ver. 4.4b モデルに円孔の電界パラメータを追加
2003/10/15 ver. 4.5b 電磁界計算のefieldの範囲の修正 入射波中心波長の変更
2003/10/17 ver. 4.6b 電界パラメータの定義の修正 屈折率分布の周期境界条件(y方向)の修正
2003/10/18 ver. 4.6 PML層を均一媒質へ戻す
2003/10/18 ver. 4.6a プログラム終了時にフィールドを吐き出すように修正 ルーチン終了後に時間を吐き出すように修正
2003/10/19 ver. 4.7 円孔部分のDBR透過屈折率に空気フィリングファクタを考慮 GaInAsPとInPの定数を入れ替え
2003/10/29 ver. 4.71 励振パルス幅を変更 励振位置を中央のみに変更
2003/11/13 ver. 4.72 円孔のペア数を活性層上下で別に定義 時間吐き出しの修正 活性層の厚さを1セル単位で変えられるように変更
2003/11/14 ver. 4.73 計算ステップ終了間際にフィールドを4/π位相がずれる間隔で吐き出しを行うように変更
2003/12/14 ver. 4.74 円孔付近の膜厚を変える方法を変更
2004/2/8 ver. 4.75 導波路幅を調整できるように変更
2004/2/17 ver. 4.76 imaxの()範囲の変更

2004/5/24 ver. 5.00 3層構造プログラムとしてモデル部の変更
2004/5/25 ver. 5.01 TEモードとTMモードのデータ吐き出しを分離 conditionファイルの吐き出しデータの内容を変更
2004/5/31 ver. 5.02 TEモードとTMモードの成分の間違いを修正
*/

/*
arranged by Daisuke Mori(2005-)

2005/4/6 ver. 5.03 モデルをSiスラブエアブリッジに変更 励振パルス幅を変更 励振を複数点に
2005/4/8 ver. 5.10 電磁界成分の吐き出しをCSV形式テキストに モデルに導波路中央への円孔追加 mcircle()に円孔半径の引数を導入
2005/4/11 ver. 5.11 フーリエ変換のファイル読み込みをxlsからテキストに修正 ikkを外部読み込みするように変更
2005/4/22 ver. 5.12 フーリエ変換のファイル読み込みが失敗していたのを修正
2005/5/16 ver. 5.13 小円孔の掘り残しモデル化に対応
2006/4/5 ver. 5.14 導波路中央の円孔半径の出力を追加(RADIUS2)
2007/1/15 ver. 5.15 導波路脇1列2列の円孔半径の出力を追加(RADIUS3,4)
2007/1/15 ver. 5.16 円孔の描画を変更(前より小さくなる)

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "FFT.c"
extern void rdft(int n, int isgn, double *a);

/************************************ 各種定数の定義 ************************************/
#define PI 3.141592					// 円周率
#define C0 2.997924e8				// 真空中の光速 [m/s]
#define EPSILON0 8.854e-12			// 真空誘電率 [F/m]
#define MU0 (PI*(4.0e-7))			// 真空中の透磁率
#define Square(x) ((x)*(x))

/****************************** プログラム全体に関する宣言 ******************************/
#define FDTD 0						// FDTD計算 0:モデルを吐き出すモード
// 1:ON 0:condition.txtとnn.xlsファイルを掃きだし，計算データが存在する場合はフーリエ変換を行う
#define PC_WAVEGUIDE 1				// Photonic Crystals WaveGuide計算 0:OFF(単位セル) 1:1列欠損導波路
//#define unitcell_number 9			// 単位セルの数(1:単位セル 奇数:1列欠損導波路(x方向単位セル奇数個))

//変更場所(始め)
#define unitcell_number 11	//7		// 単位セルの数(1:単位セル 奇数:1列欠損導波路(x方向単位セル奇数個))
//変更場所(終り)

// PC_WAVEGUIDEとunitcell_numberの設定内容は統一させること
#define K_STARTPOINT 4				// 波数の計算方向開始点の設定 1:G-X 2:G-J 3:X-J 4:G-P(導波路計算の場合)
#define K_ENDPOINT 4				// 波数の計算方向終了点の設定
// このプログラムではx方向に吸収境界を設定しているので，G-Pのみ対応
//#define IKK_STARTPOINT 1			// 波数の最初の値 0から始めない方がいい
//#define IKK_ENDPOINT 18				// 波数の最後の値 IKK_MAXまで計算すると端まで計算したことになる

//変更場所v5.11(始め)
//#define IKK_STARTPOINT 10			// 波数の最初の値 0から始めない方がいい
//#define IKK_ENDPOINT 10				// 波数の最後の値 IKK_MAXまで計算すると端まで計算したことになる
int IKK_STARTPOINT;				// 波数の最初の値 0から始めない方がいい
int IKK_ENDPOINT;					// 波数の最後の値 IKK_MAXまで計算すると端まで計算したことになる
//変更場所(終り)

//#define IKK_MAX 20					// 波数の分割数 バンドを求める際の点数に相当

//変更場所(始め)
#define IKK_MAX 50					// 波数の分割数 バンドを求める際の点数に相当
//変更場所(終り)

#define AIR 1						// 配列nn[][][]における空気の定数
//#define SI 0						// 配列nn[][][]におけるSiの定数

//変更場所(始め)
#define SIO2 2	// 配列nn[][][]におけるSiO2の定数
//変更場所(終わり)
					
//#define ALL_SPACE_MEDIUM SI			// 空間を初期化する媒質

#define CLAD_MEDIUM AIR			// クラッド層媒質
#define CIRCLE_MEDIUM AIR			// 結晶内媒質

//変更場所(始め)
#define InP_MEDIUM 2				// 配列nn[][][]におけるInPの定数
#define Act_MEDIUM 3				// 配列nn[][][]におけるActの定数
//変更場所(終り)

#define REF_INDEX1 3.43				// Siの屈折率(コア)
//#define REF_INDEX2 1.0					// Airクラッドの屈折率
#define REF_INDEX2 1.444				// SiO2クラッドの屈折率

//変更場所(始め)
#define InP_INDEX 3.1700			// InPの屈折率(使用していない値)
//#define Act_INDEX 3.534			// 活性層の屈折率******変更した*******(Ag-As2Se3)
#define Act_INDEX REF_INDEX1		// コアの屈折率
//変更場所(終り)

/********************************* PMLプログラムの宣言 **********************************/
#define L 32						// PML層数
#define M 2.0						// 導電率の分布を与える次数
#define order 16					// sigma_max算出に使う反射係数[R0=10^(-order)]の次数

/***************************** 設計値・励振関数に関する宣言 *****************************/
//#define PITCH  0.28e-6				// 三角配列のピッチ
//#define RADIUS 0.08e-6				// 円孔，円柱の半径
//#define CORE 0.16e-6				// コアの厚さ
//#define CLAD 0.80e-6				// クラッドの厚さ

//変更場所(始め)
#define PITCH   0.399e-6				// 三角配列のピッチ
#define RADIUS  0.110e-6				// 円孔，円柱の半径1列目
#define RADIUSZ 0.110e-6			// 円孔，円柱の半径2列目
#define RADIUS3 0.110e-6			// 円孔，円柱の半径3列目
#define RADIUS4 0.110e-6			// 円孔，円柱の半径4列目
#define RADIUS5 0.110e-6			// 円孔，円柱の半径5列目
#define RADIUS6 0.110e-6			// 円孔，円柱の半径6列目
#define RADIUS7 0.110e-6			// 円孔，円柱の半径7列目
#define RADIUSex1 0.000e-6			// 円孔，円柱の半径ex1列目（基本的に0にすること）
#define RADIUSex2 0.000e-6			// 円孔，円柱の半径ex2列目（基本的に0にすること）
#define RADIUSd 0.095e-6			//異なる円孔の半径       //追加take
#define LENGTH_L 0.420e-6			//追加take


#define RADIUS2 0.000e-6			// 円孔，円柱の半径2
//#define RADIUS3 0.14e-6			// 円孔，円柱の半径3(導波路脇)
//#define RADIUS4 0.14e-6			// 円孔，円柱の半径4(導波路脇2列目)

#define AIR_clad 0.0e-6				// 空気層の厚さ
#define InP_circle_over 1.5e-6		// 活性層上部の円孔深さ
#define InP_circle_down	1.5e-6		// 活性層下部の円孔深さ
#define InP_sub 0.0e-6				// InP基板の厚さ
#define Act 0.210e-6					// 活性層の厚さ
#define Act_remain 0.00e-6			// 活性層の掘り残し
#define Width 0						// 導波路幅の調整(偶数を入力すること)

#define FIRST_Sx 0.000e-6				// 1列目のx方向のシフト量(セルサイズで割り切れること！！）
#define FIRST_Sy 0.000e-6               // 1列目のy方向のシフト量(+で狭める方向)

#define SECOND_Sx 0.000e-6				// 2列目のx方向のシフト量(セルサイズで割り切れること！！）
#define SECOND_Sy 0.000e-6              // 2列目のy方向のシフト量(+で狭める方向)

#define THIRD_Sx 0.000e-6				// 3列目のx方向のシフト量(セルサイズで割り切れること！！）
#define THIRD_Sy 0.000e-6               //　3列目のy方向のシフト量(+で狭める方向)

#define FOURTH_Sx 0.000e-6				// 4列目のx方向のシフト量(セルサイズで割り切れること！！）
#define FOURTH_Sy 0.000e-6             // 4列目のy方向のシフト量(+で狭める方

#define FIFTH_Sx 0.000e-6				// 5列目のx方向のシフト量(セルサイズで割り切れること！！）
#define FIFTH_Sy 0.000e-6               // 5列目のy方向のシフト量(+で狭める方

#define SIXTH_Sx 0.000e-6				// 6列目のx方向のシフト量(セルサイズで割り切れること！！）
#define SIXTH_Sy 0.000e-6               // 6列目のy方向のシフト量(+で狭める方

#define SEVENTH_Sx 0.000e-6		   	    // 7列目のx方向のシフト量(セルサイズで割り切れること！！）
#define SEVENTH_Sy 0.000e-6             // 7列目のy方向のシフト量(+で狭める方

#define Ex1_Sx 0.000e-6		   	    // 7列目のx方向のシフト量(セルサイズで割り切れること！！）
#define Ex1_Sy 0.000e-6             // 7列目のy方向のシフト量(+で狭める方

#define Ex2_Sx 0.000e-6		   	    // 7列目のx方向のシフト量(セルサイズで割り切れること！！）
#define Ex2_Sy 0.000e-6             // 7列目のy方向のシフト量(+で狭める方


//変更場所(終り)

//#define NORMALFREQ 0.30				// 入射波中心規格化周波数
//#define LAMBDA (PITCH/NORMALFREQ)	// 入射波中心波長．これは変更しない

//変更場所(始め)
#define LAMBDA 1.541e-6				// 入射波中心波長．これは変更しない(1.45e-6)
//変更場所(終り)

//#define CELLLAMBDA 0.700e-6			// セルサイズを決める波長，通常は変更しない
#define excite 1					// 励振の種類 1:TE励振 2:TM励振
#define omega0 (2.0*PI*C0/LAMBDA)	// 励振関数の角周波数
//#define sigma  (omega0*0.3/1.5)		// 励振関数のパルス幅を決める定数
//#define PULSEPEAK 2000				// 励振関数パルスのピークとなる時間 NORMALFREQを変えた場合は励振パルスが切れないように注意

//変更場所(始め)
#define sigma  (omega0*0.005)		// 励振関数のパルス幅を決める定数(バンド計算*0.03，モード断面積*0.005)
#define PULSEPEAK 22000			// 励振関数パルスのピークとなる時間 NORMALFREQを変えた場合は励振パルスが切れないように注意(バンド計算1000，モード断面積22000)
//変更場所(終り)

#define nmax 262144  				// 時間ステップの最大値 基本的に2の乗数
// 2^15=32768 2^16=65536 2^17=131072 2^18=262144 2^19=524288 2^20=1048576 2^21=2097152

/********************************** 基本的にいじらない **********************************/
//#define lambda(x) (CELLLAMBDA/(REF_INDEX1*x))	// 波長/(屈折率*x)
//#define dx lambda(10)							// 微小空間分割サイズX方向
//#define dy lambda(10)							// 微小空間分割サイズY方向
//#define dz lambda(10)							// 微小空間分割サイズZ方向

//変更場所(始め)
#define dx 0.021e-6								// 微小空間分割サイズX方向
#define dy 0.021e-6								// 微小空間分割サイズY方向
#define dz 0.021e-6								// 微小空間分割サイズZ方向
//変更場所(終り)

#define dt (dx*0.0185e-7)						// 微小時間分割サイズ(Courant安定条件に注意！)
//#define EPSILON1 (Square(REF_INDEX1)*EPSILON0)	// Siの誘電率

#define EPSILON2 (Square(REF_INDEX2)*EPSILON0)				// クラッドの誘電率

//変更場所(始め)
#define InP_EPSILON (Square(InP_INDEX)*EPSILON0)			// InPの誘電率
#define Act_EPSILON (Square(Act_INDEX)*EPSILON0)			// 活性層の誘電率
//変更場所(終り)

#define cnstM (dt/MU0)							// Yeeルーチン中で使う定数(磁界)
//#define cnstE (dt/EPSILON0)						// Yeeルーチン中で使う定数(真空中の電界)
//#define cnstE1 (dt/EPSILON1)					// Yeeルーチン中で使う定数(Siコア中の電界)
#define cnstE2 (dt/EPSILON2)					// Yeeルーチン中で使う定数(SiO2クラッド中の電界)

//変更場所(始め)
#define Air_E (dt/EPSILON2)						// Yeeルーチン中で使う真空中の電界 //10/27にSIO2に変更(cnstE2と同じ)
#define InP_E (dt/InP_EPSILON)					// InP中の電界
#define Act_E (dt/Act_EPSILON)					// 活性層中の電界
//変更場所(終り)

/******************************** フーリエ変換関連の宣言 ********************************/
#define STEP nmax+1000				// 入力ファイルの値の総数(FDTDで計算したときのステップ数)
#define F_NUMBER 5					// ファイル中の観測点の数
//#define ELE 6						// 成分の数(Ex,Ey,Ez,Hx,Hy,Hz)の順番に計算

//変更場所(始め)
#define ELE 3						// 成分の数(TE励振:Ex,Ey,Hz, TM励振:Ez,Hx,Hy)
//変更場所(終り)

/*********************************** セル数などの宣言 ***********************************/
const int A=(int)((PITCH*1.0e9)/(dx*1.0e9));	// セル長で規格化したprogram用のフォトニック結晶ピッチ
//#define B 24									// A*root(3)の値を入れる(三角格子の場合) 必ず偶数になるようにする
//const int C=(int)((CORE*1.0e9)/(dz*1.0e9));		// セル長で規格化したprogram用のコア厚
//const int D=(int)((CLAD*1.0e9)/(dz*1.0e9));		// セル長で規格化したprogram用のクラッド厚

//変更場所(始め)
#define B 32																// A*root(3)の値を入れる(三角格子の場合) 必ず偶数になるようにする
const int AIR_clad_unit=(int)((AIR_clad*1.0e9)/(dz*1.0e9));					// セル長で規格化したprogram用のAir厚
const int InP_circle_over_unit=(int)((InP_circle_over*1.0e9)/(dz*1.0e9));	// セル長で規格化したprogram用の上部円孔厚
const int InP_circle_down_unit=(int)((InP_circle_down*1.0e9)/(dz*1.0e9));	// セル長で規格化したprogram用の下部円孔厚
const int InP_sub_unit=(int)((InP_sub*1.0e9)/(dz*1.0e9));					// セル長で規格化したprogram用のInP基板厚
const int Act_unit=(int)((Act*1.0e9)/(dz*1.0e9));							// セル長で規格化したprogram用の活性層厚
const int Act_remain_unit=(int)((Act_remain*1.0e9)/(dz*1.0e9));							// セル長で規格化したprogram用の活性層掘り残し
//変更場所(終り)

//#define imax (B*unitcell_number+L+L)			// x方向のセル数

//変更場所(始め)
#define imax (B*unitcell_number+L+L+Width)		// x方向のセル数
//変更場所(始め)

const int jmax=A;								// y方向のセル数
//#define kmax (C*12+L+L)							// z方向のセル数

//変更場所(始め)
#define kmax (AIR_clad_unit+InP_circle_over_unit+InP_circle_down_unit+InP_sub_unit+Act_unit+L+L)	//z方向のセル数
//変更場所(終り)

#define imove 5									// 観測点を定義するパラメータ(解析空間の中心からの距離)
#define jmove 5									// 観測点を定義するパラメータ(解析空間の中心からの距離)

//変更場所(始め)
#define kmove 5									// 観測点を定義するパラメータ(活性層の中心からの距離)
//変更場所(終り)

/******************************* フィールド吐き出しの設定 *******************************/
#define PrintStat 260010								// 電磁界分布を出す時間ステップの始め		
#define PrintStat_f 100
#define PrintEnd 260210							// 電磁界分布を出す時間ステップの終わり(nmax)		
#define PrintEnd_f 260000
//#define PrintNum 1000000						// 電磁界分布を出す時間ステップの間隔 出さないときはnmaxよりも大きい値を指定

//変更場所(始め)
#define PrintNum 10 						// 電磁界分布を出す時間ステップの間隔 出さないときはnmaxよりも大きい値を指定(通常20000)
#define PrintNum_f 100
//変更場所(終り)

/************************************ 各種関数の宣言 ************************************/
void Init_Space(void);							// 解析モデルを設定
void k_vector_phase(void);						// 波数ベクトルと位相情報の設定
void exc_function(void);						// 励振関数の定義
void calc_hfield(void);							// 磁界計算
void calc_h_periodic_boundary_condition(void);	// 磁界計算用周期境界条件
void calc_efield(void);							// 電界計算
void calc_e_periodic_boundary_condition(void);	// 電界計算用周期境界条件
void calc_e_periodic_boundary_condition2(void);	// 電界計算用周期境界条件
void efield_PML(void);							// PML電界計算
void hfield_PML(void);							// PML磁界計算
void circle_arrange(void);						// 円孔の配置
void mcircle(int ,int ,double);						// make circle function
void halfcircle1(int ,int ,double);					// make half circle function 1
void halfcircle2(int ,int ,double);					// make half circle function 2
void halfcircle3(int ,int ,double);					// make half circle function 3
void halfcircle4(int ,int ,double);					// make half circle function 4
void rightquartercircle1(int ,int ,double);
void leftquartercircle1(int ,int ,double);
void rightquartercircle2(int ,int ,double);
void leftquartercircle2(int ,int ,double);
//******変更した*******(ここから)
void msiftcircle(int ,int ,double);
void rightsiftquartercircle1(int ,int ,double);	
void rightsiftquartercircle2(int ,int ,double);	
void leftsiftquartercircle1(int , int , double );
void leftsiftquartercircle2(int , int , double );
void mtangle1(int ,int ,double );				//変更take
void mtangle2(int ,int ,double );				//変更take
//******変更した*******(ここまで)
void nnprint(void);								// 屈折率プリント
void file_open_field(void);						// フィールド吐き出し用ファイルオープン
void field_print(int);							// electric field and magnetic field print
void file_close_field(void);					// フィールド吐き出し用ファイルクローズ
void file_open(void);							// パルス吐き出し用ファイルオープン
void Poyntingprint(int);						// Ex Ey Ez Hx Hy Hz for "poynting calc" print
void file_close(void);							// パルス吐き出し用ファイルクローズ
void initialize_matrix(void);					// 配列の初期化
void fft_main(void);							// フーリエ変換プログラム
void Output_Condition(void);					// 各パラメータの値をconditionファイルに掃き出す

void file_open_field_f(void);						// フィールド吐き出し用ファイルオープン
void field_print_f(int);							// electric field and magnetic field print
void file_close_field_f(void);					// フィールド吐き出し用ファイルクローズ


/****************************** 各種パラメータ・配列の宣言 ******************************/
double kkx,kky,kkz,G1x,G1y,G2x,G2y,Gz;							// 波数ベクトルと逆格子ベクトル
double cos_phx,sin_phx,cos_phy,sin_phy,cos_phxy,sin_phxy;		// 周期境界条件を導入する際の位相情報を含ませた変数

//変更場所(始め)
double cnstE[imax][jmax+1][kmax];									// 解析空間の物質ごとの電界
//変更場所(終り)

double re_ex[imax][jmax+1][kmax+1],re_ey[imax+1][jmax][kmax+1],re_ez[imax+1][jmax+1][kmax];	// 解析空間の電界実数部成分
double im_ex[imax][jmax+1][kmax+1],im_ey[imax+1][jmax][kmax+1],im_ez[imax+1][jmax+1][kmax];	// 解析空間の電界虚数部成分
double re_hx[imax+1][jmax+1][kmax],re_hy[imax][jmax+1][kmax],re_hz[imax][jmax+1][kmax+1];	// 解析空間の磁界実数部成分
double im_hx[imax+1][jmax+1][kmax],im_hy[imax][jmax+1][kmax],im_hz[imax][jmax+1][kmax+1];	// 解析空間の磁界虚数部成分
double re_hx_jmax[imax+1][kmax],im_hx_jmax[imax+1][kmax];		// 周期境界条件用磁界成分
double re_hz_jmax[imax][kmax+1],im_hz_jmax[imax][kmax+1];		// 周期境界条件用磁界成分
static short int nn[imax][jmax+1][kmax];						// 屈折率分布
//FILE *EX;			FILE *EY;			FILE *EZ;
//FILE *HX;			FILE *HY;			FILE *HZ;

//変更場所(始め)
FILE *EX_XY;			FILE *EY_XY;			FILE *EZ_XY;
FILE *EX_XZ_center;		FILE *EY_XZ_center;		FILE *EZ_XZ_center;
FILE *EX_XZ_end;		FILE *EY_XZ_end;		FILE *EZ_XZ_end;
FILE *EX_YZ;			FILE *EY_YZ;			FILE *EZ_YZ;
FILE *HX_XY;			FILE *HY_XY;			FILE *HZ_XY;
FILE *HX_XZ_center;		FILE *HY_XZ_center;		FILE *HZ_XZ_center;
FILE *HX_XZ_end;		FILE *HY_XZ_end;		FILE *HZ_XZ_end;
FILE *HX_YZ;			FILE *HY_YZ;			FILE *HZ_YZ;
//変更場所(終り)

FILE *fPoynt1Ex;	FILE *fPoynt1Ey;	FILE *fPoynt1Ez;
FILE *fPoynt1Hx;	FILE *fPoynt1Hy;	FILE *fPoynt1Hz;
FILE *fp1;
FILE *fp2;
//int n,m;							// n:時間ステップ m:波数方向設定判定子

//変更場所(始め)
int o=16384,n,m;					// o:時間ステップを吐き出すステップ n:時間ステップ m:波数方向設定判定子
//変更場所(終り)

int ikk;							// ikk:波数分割子
int ci,cj,ck;						// 解析空間の中心座標
int i,j,k;

//変更場所(始め)
static int a,b;						// a:活性層より下のセル数 b:SiO2層より下のセル数
//変更場所(終り)

int effective_circle_cell;
double effective_2r_a;

char fname[30];

char fikk[30];

/**************************** PML関連パラメータ・配列の宣言 *****************************/
//#define PML_coefficient cnstE1			// PML内で使用するパラメータ(基本はcnstEでOK 場合によりcnstE1と使い分けること)

//変更場所(始め)
#define PML_coefficient Act_E			// PML内で使用するパラメータ(基本はcnstEでOK 場合によりcnstE1と使い分けること)
//変更場所(終り)

double sigma_max;						// 外壁での導電率
double sigma_x,sigma_y,sigma_z;			// 導電率
double u1,u2;							// PML計算に使う定数
double term1,term2;						// PML計算に使う定数
double re_exz010[L][jmax+1][L],re_exy010[L][jmax+1][L],re_exz210[L][jmax+1][L],re_exy210[L][jmax+1][L],re_exz012[L][jmax+1][L],re_exy012[L][jmax+1][L],re_exz212[L][jmax+1][L],re_exy212[L][jmax+1][L];
double im_exz010[L][jmax+1][L],im_exy010[L][jmax+1][L],im_exz210[L][jmax+1][L],im_exy210[L][jmax+1][L],im_exz012[L][jmax+1][L],im_exy012[L][jmax+1][L],im_exz212[L][jmax+1][L],im_exy212[L][jmax+1][L];
double re_ezx010[L][jmax+1][L],re_ezy010[L][jmax+1][L],re_ezx210[L][jmax+1][L],re_ezy210[L][jmax+1][L],re_ezx012[L][jmax+1][L],re_ezy012[L][jmax+1][L],re_ezx212[L][jmax+1][L],re_ezy212[L][jmax+1][L];
double im_ezx010[L][jmax+1][L],im_ezy010[L][jmax+1][L],im_ezx210[L][jmax+1][L],im_ezy210[L][jmax+1][L],im_ezx012[L][jmax+1][L],im_ezy012[L][jmax+1][L],im_ezx212[L][jmax+1][L],im_ezy212[L][jmax+1][L];
double re_hyx010[L][jmax+1][L],re_hyz010[L][jmax+1][L],re_hyx210[L][jmax+1][L],re_hyz210[L][jmax+1][L],re_hyx012[L][jmax+1][L],re_hyz012[L][jmax+1][L],re_hyx212[L][jmax+1][L],re_hyz212[L][jmax+1][L];
double im_hyx010[L][jmax+1][L],im_hyz010[L][jmax+1][L],im_hyx210[L][jmax+1][L],im_hyz210[L][jmax+1][L],im_hyx012[L][jmax+1][L],im_hyz012[L][jmax+1][L],im_hyx212[L][jmax+1][L],im_hyz212[L][jmax+1][L];
double re_eyx010[L][jmax][L],re_eyz010[L][jmax][L],re_eyx210[L][jmax][L],re_eyz210[L][jmax][L],re_eyx012[L][jmax][L],re_eyz012[L][jmax][L],re_eyx212[L][jmax][L],re_eyz212[L][jmax][L];
double im_eyx010[L][jmax][L],im_eyz010[L][jmax][L],im_eyx210[L][jmax][L],im_eyz210[L][jmax][L],im_eyx012[L][jmax][L],im_eyz012[L][jmax][L],im_eyx212[L][jmax][L],im_eyz212[L][jmax][L];
double re_hxz010[L][jmax][L],re_hxy010[L][jmax][L],re_hxz210[L][jmax][L],re_hxy210[L][jmax][L],re_hxz012[L][jmax][L],re_hxy012[L][jmax][L],re_hxz212[L][jmax][L],re_hxy212[L][jmax][L];
double im_hxz010[L][jmax][L],im_hxy010[L][jmax][L],im_hxz210[L][jmax][L],im_hxy210[L][jmax][L],im_hxz012[L][jmax][L],im_hxy012[L][jmax][L],im_hxz212[L][jmax][L],im_hxy212[L][jmax][L];
double re_hzx010[L][jmax][L],re_hzy010[L][jmax][L],re_hzx210[L][jmax][L],re_hzy210[L][jmax][L],re_hzx012[L][jmax][L],re_hzy012[L][jmax][L],re_hzx212[L][jmax][L],re_hzy212[L][jmax][L];
double im_hzx010[L][jmax][L],im_hzy010[L][jmax][L],im_hzx210[L][jmax][L],im_hzy210[L][jmax][L],im_hzx012[L][jmax][L],im_hzy012[L][jmax][L],im_hzx212[L][jmax][L],im_hzy212[L][jmax][L];
double re_exz110[imax-L][jmax+1][L],re_exy110[imax-L][jmax+1][L],re_exz112[imax-L][jmax+1][L],re_exy112[imax-L][jmax+1][L];
double im_exz110[imax-L][jmax+1][L],im_exy110[imax-L][jmax+1][L],im_exz112[imax-L][jmax+1][L],im_exy112[imax-L][jmax+1][L];
double re_hyx110[imax-L][jmax+1][L],re_hyz110[imax-L][jmax+1][L],re_hyx112[imax-L][jmax+1][L],re_hyz112[imax-L][jmax+1][L];
double im_hyx110[imax-L][jmax+1][L],im_hyz110[imax-L][jmax+1][L],im_hyx112[imax-L][jmax+1][L],im_hyz112[imax-L][jmax+1][L];
double re_eyx011[L][jmax][kmax-L+1],re_eyz011[L][jmax][kmax-L+1],re_eyx211[L][jmax][kmax-L+1],re_eyz211[L][jmax][kmax-L+1];
double im_eyx011[L][jmax][kmax-L+1],im_eyz011[L][jmax][kmax-L+1],im_eyx211[L][jmax][kmax-L+1],im_eyz211[L][jmax][kmax-L+1];
double re_hzx011[L][jmax][kmax-L+1],re_hzy011[L][jmax][kmax-L+1],re_hzx211[L][jmax][kmax-L+1],re_hzy211[L][jmax][kmax-L+1];
double im_hzx011[L][jmax][kmax-L+1],im_hzy011[L][jmax][kmax-L+1],im_hzx211[L][jmax][kmax-L+1],im_hzy211[L][jmax][kmax-L+1];
double re_eyx110[imax-L+1][jmax][L],re_eyz110[imax-L+1][jmax][L],re_eyx112[imax-L+1][jmax][L],re_eyz112[imax-L+1][jmax][L];
double im_eyx110[imax-L+1][jmax][L],im_eyz110[imax-L+1][jmax][L],im_eyx112[imax-L+1][jmax][L],im_eyz112[imax-L+1][jmax][L];
double re_hxz110[imax-L+1][jmax][L],re_hxy110[imax-L+1][jmax][L],re_hxz112[imax-L+1][jmax][L],re_hxy112[imax-L+1][jmax][L];
double im_hxz110[imax-L+1][jmax][L],im_hxy110[imax-L+1][jmax][L],im_hxz112[imax-L+1][jmax][L],im_hxy112[imax-L+1][jmax][L];
double re_ezx011[L][jmax+1][kmax-L],re_ezy011[L][jmax+1][kmax-L],re_ezx211[L][jmax+1][kmax-L],re_ezy211[L][jmax+1][kmax-L];
double im_ezx011[L][jmax+1][kmax-L],im_ezy011[L][jmax+1][kmax-L],im_ezx211[L][jmax+1][kmax-L],im_ezy211[L][jmax+1][kmax-L];
double re_hyx011[L][jmax+1][kmax-L],re_hyz011[L][jmax+1][kmax-L],re_hyx211[L][jmax+1][kmax-L],re_hyz211[L][jmax+1][kmax-L];
double im_hyx011[L][jmax+1][kmax-L],im_hyz011[L][jmax+1][kmax-L],im_hyx211[L][jmax+1][kmax-L],im_hyz211[L][jmax+1][kmax-L];
double re_exz011[L][jmax+1][kmax-L+1],re_exy011[L][jmax+1][kmax-L+1],re_exz211[L][jmax+1][kmax-L+1],re_exy211[L][jmax+1][kmax-L+1];
double im_exz011[L][jmax+1][kmax-L+1],im_exy011[L][jmax+1][kmax-L+1],im_exz211[L][jmax+1][kmax-L+1],im_exy211[L][jmax+1][kmax-L+1];
double re_ezx110[imax-L+1][jmax+1][L],re_ezy110[imax-L+1][jmax+1][L],re_ezx112[imax-L+1][jmax+1][L],re_ezy112[imax-L+1][jmax+1][L];
double im_ezx110[imax-L+1][jmax+1][L],im_ezy110[imax-L+1][jmax+1][L],im_ezx112[imax-L+1][jmax+1][L],im_ezy112[imax-L+1][jmax+1][L];
double re_hxz011[L][jmax][kmax-L],re_hxy011[L][jmax][kmax-L],re_hxz211[L][jmax][kmax-L],re_hxy211[L][jmax][kmax-L];
double im_hxz011[L][jmax][kmax-L],im_hxy011[L][jmax][kmax-L],im_hxz211[L][jmax][kmax-L],im_hxy211[L][jmax][kmax-L];
double re_hzx110[imax-L][jmax][L],re_hzy110[imax-L][jmax][L],re_hzx112[imax-L][jmax][L],re_hzy112[imax-L][jmax][L];
double im_hzx110[imax-L][jmax][L],im_hzy110[imax-L][jmax][L],im_hzx112[imax-L][jmax][L],im_hzy112[imax-L][jmax][L];

/************************************ メインルーチン ************************************/
int main(){
	printf("dx = %e, dy = %e, dz = %e\n",dx,dy,dz);	// セルサイズの表示
	printf("dt = %e\n",dt);							// 微小時間ステップの表示
//	printf("半径のセル数 = %d\n",(int)((RADIUS*5.0e9)/(dx*5.0e9)));	// 半径のセル数の表示
//	printf("ピッチのセル数 = %d\n",A);				// ピッチのセル数の表示
//	printf("コア層のセル数 = %d\n",C);				// コア層のセル数の表示
//	printf("クラッド層のセル数 = %d\n",D);			// クラッド層のセル数の表示

//変更場所v5.11(始め)
	FILE *ikk_file;
	ikk_file = fopen("ikk.txt","r");
	fscanf(ikk_file,"%d",&IKK_STARTPOINT);
	IKK_ENDPOINT = IKK_STARTPOINT;
//変更場所(終り)

//変更場所(始め)
	printf("円孔の半径のセル数 = %d",(int)((RADIUS*5.0e9)/(dx*5.0e9)));					// 円孔の半径のセル数の表示
	printf("		円孔の半径のセル数2 = %d\n",(int)((RADIUS2*5.0e9)/(dx*5.0e9)));		// 円孔の半径のセル数の表示2
	printf("ピッチのセル数 = %d\n",A);													// ピッチのセル数の表示
	printf("導波路幅変化のセル数 = %d\n",Width);	
	printf("1列目シフトのセル数 = %d\n",int(FIRST_Sx/dx));								// ピッチのセル数の表示
	printf("1列目シフトのセル数 = %d\n",int(FIRST_Sy/dy));
	printf("2列目シフトのセル数 = %d\n",int(SECOND_Sx/dx));
	printf("2列目シフトのセル数 = %d\n",int(SECOND_Sy/dy));
	printf("3列目シフトのセル数 = %d\n",int(THIRD_Sx/dx));
	printf("3列目シフトのセル数 = %d\n",int(THIRD_Sy/dy));
	printf("4列目シフトのセル数 = %d\n",int(FOURTH_Sx/dx));
	printf("4列目シフトのセル数 = %d\n",int(FOURTH_Sy/dy));
	printf("活性層上部の円孔深さ = %d",InP_circle_over_unit);							// 活性層上部の円孔深さのセル数表示
	printf("	活性層下部の円孔深さ = %d\n",InP_circle_down_unit);						// 活性層下部の円孔深さのセル数表示
	printf("InP基板のセル数 = %d",InP_sub_unit);										// InP基板のセル数の表示
	printf("		活性層のセル数 = %d",Act_unit);										// 活性層のセル数の表示
	printf("	空気層のセル数 = %d\n",AIR_clad_unit);									// 空気層のセル数の表示
	
	if(excite == 1){
		printf("TE-like mode\n");
	}
	if(excite == 2){
		printf("TM-like mode\n");
	}
//変更場所(終り)

	ci=imax/2;								// 解析空間の中心を設定
	cj=jmax/2;								// 解析空間の中心を設定
//	ck=kmax/2;								// 解析空間の中心を設定
	printf("imax = %d	jmax = %d	kmax = %d\n",imax,jmax,kmax);		// 解析空間座標の表示
//	printf("ci = %d, cj = %d, ck = %d\n",ci,cj,ck);					// 解析空間の中心座標の表示

	sigma_max=-((M+1)*EPSILON0*C0*log(1.0*pow(double(10),double(-order))))/(2*L*dx);	// PMLパラメータ(sigma_max)の計算
//	printf("Number of PML : %d, PML constant M : %.2e, Reflection order : %d, sigma_max : %e\n",L,M,order,sigma_max);	// PMLパラメータの表示

	Init_Space();							// 解析モデルを設定
	Output_Condition();						// 各種パラメータをcondition.txtファイルに出力

#if FDTD==1
	printf("FDTD simulation START\n");

	for(m=K_STARTPOINT;m<=K_ENDPOINT;m++){
		for(ikk=IKK_STARTPOINT;ikk<=IKK_ENDPOINT;ikk++){
			file_open();
			k_vector_phase();
			initialize_matrix();
			printf("m = %d/%d, ikk = %d/%d\n",m,K_ENDPOINT,ikk,IKK_ENDPOINT);

//変更場所(始め)
			FILE *tf;
			tf = fopen("time.txt","w");
			fprintf(tf,"FDTD start: ");
			fprintf(tf,"\n");

			char time[9];
			_strtime( time );
			printf("FDTD start:\n");
			printf("t=0/%d %s\n",nmax,time);
			fprintf(tf,"n = 0 : ");
			fputs(time,tf);
			fprintf(tf,"\n");
//変更場所(終り)

			for(n=1;n<=nmax;n++){
//				char time[9];
//				_strtime( time );
//				if(n%1000==0) printf("t=%d/%d %s\n",n,nmax,time);

				exc_function();

				calc_hfield();
				hfield_PML();
				calc_h_periodic_boundary_condition();
				calc_efield();
				calc_e_periodic_boundary_condition2();
				efield_PML();
				calc_e_periodic_boundary_condition();

				Poyntingprint(1);

//				if(n>=PrintStat && n<=PrintEnd && (n%PrintNum)==0){

//変更場所(始め)
				if(n>=PrintStat && n<=PrintEnd && (n%PrintNum)==0 || n==nmax){
//変更場所(終り)

					file_open_field();
					field_print(n);
					file_close_field();
				}

				if(n>=PrintStat_f && n<=PrintEnd_f && (n%PrintNum_f)==0 || n==nmax){
//変更場所(終り)

					file_open_field_f();
					field_print_f(n);
					file_close_field_f();
				}




//変更場所(始め)
				_strtime( time );
				if(n%10000==0 || n==o){
					printf("t=%d/%d %s\n",n,nmax,time);
					fprintf(tf,"n = %d : ",n);
					fputs(time,tf);
					fprintf(tf,"\n");
					if(n==o){
						o=o*2;
					}
				}
//変更場所(終り)

			}
			file_close();

//変更場所(始め)
			fclose(tf);
//変更場所(終り)

		}
	}
#endif

	Poyntingprint(0);
	nnprint();
	fft_main();

	// FFTファイルの移動
	char OutputFile1[20], OutputFile2[20];
	m = 4;	// 数字直打ちなのでこれ以降にプログラムを追加するときには注意が必要
	sprintf(OutputFile1, "0%d_FFT%s.xls", m, fikk);
	sprintf(OutputFile2, "..\\%s", OutputFile1);
	if( rename( OutputFile1, OutputFile2 ) != 0 )
	
	return 0;
}
/********************************* メインルーチン終わり *********************************/


/************************************* モデルの定義 *************************************/
void Init_Space(){
//	int i,j,k;

//変更場所(始め)
	int i,j,k;
//変更場所(終り)

	for(i=0;i<imax;i++){
		for(j=0;j<=jmax;j++){
//			for(k=0;k<kmax;k++){
//				nn[i][j][k]=ALL_SPACE_MEDIUM;	// 解析空間を初期化
//			}
//			for(k=0;k<(ck-(C/2)-D);k++){
//				nn[i][j][k]=CIRCLE_MEDIUM;		// クラッド層下部を空気層とする
//			}
//			for(k=(ck-(C/2)-D);k<(ck-(C/2));k++){
//				nn[i][j][k]=CLAD_MEDIUM;		// クラッド層を定義
//			}
//			for(k=(ck+(C/2));k<kmax;k++){
//				nn[i][j][k]=CIRCLE_MEDIUM;		// コア層上部を空気層とする
//			}

//変更場所(始め)
			for(k=0;k<kmax;k++){
				nn[i][j][k]=Act_MEDIUM;							// 解析空間の初期化(全パラメータを活性層に)
				cnstE[i][j][k]=Act_E;
			}
			for(k=L;k<(L+InP_sub_unit);k++){					// InP_subの定義
				nn[i][j][k]=InP_MEDIUM;
				cnstE[i][j][k]=InP_E;
			}
//変更場所(始め)
			a=InP_sub_unit+L;
			for(k=a;k<(a+InP_circle_down_unit);k++){			// InP_circle_downの定義
				nn[i][j][k]=CLAD_MEDIUM;
				cnstE[i][j][k]=cnstE2;
			}
			b=InP_circle_down_unit+a;
			for(k=b;k<(b+Act_unit);k++){						// 活性層を定義
				nn[i][j][k]=Act_MEDIUM;
				cnstE[i][j][k]=Act_E;
			}
			b=Act_unit+b;
			for(k=b;k<(b+InP_circle_over_unit);k++){			// InP_circle_overの定義
				nn[i][j][k]=CLAD_MEDIUM;
				cnstE[i][j][k]=cnstE2;
//変更場所(終わり)
			}
			b=InP_circle_over_unit+b;
			for(k=b;k<(b+AIR_clad_unit);k++){					// AIR_cladの定義
				nn[i][j][k]=AIR;
				cnstE[i][j][k]=Air_E;
			}
//変更場所(終り)

		}
	}



//変更場所(始め)
	ck=a+InP_circle_down_unit+Act_unit/2;						//z方向の解析空間の中心を活性層中心に設定
	printf("ci = %d	cj = %d		ck = %d\n",ci,cj,ck);			// 解析空間の中心座標の表示

	printf("Number of PML : %d		PML constant M : %.2e\nReflection order : %d		sigma_max : %e\n",L,M,order,sigma_max);	// PMLパラメータの表示
//変更場所(終り)

	printf("triangular array ON\n");
	circle_arrange();						// 解析モデルに円孔を配置

	for(i=0;i<imax;i++){
		for(k=0;k<kmax;k++){
			nn[i][jmax][k]=nn[i][0][k];		// 屈折率分布の周期境界条件(y方向)

//変更場所(始め)
			cnstE[i][jmax][k]=cnstE[i][0][k];	// 屈折率分布の周期境界条件(y方向)
//変更場所(終り)

		}
	}
}
/*
何故屈折率分布にも周期境界条件を導入するのかというと，
これを導入しないと電磁界計算の周期境界面で実は屈折率定数がいい加減な値が代入されてしまう．
何気にかなり重要．
*/

void circle_arrange(){
	int place;
	for(place=0;place<=(unitcell_number*2);place++){
//******変更した*******(ここから)
//シフト列3列目(=偶数)のときの円孔生成
		if(place == unitcell_number-3){
			msiftcircle(((B*place)/2)+L+THIRD_Sy/dx,A/2+THIRD_Sx/dy,RADIUS3);
			msiftcircle(((B*place)/2)+L+THIRD_Sy/dx,A*3/2+THIRD_Sx/dy,RADIUSd);                  //変更take
		}
		else if(place == unitcell_number+3){
			msiftcircle(((B*place)/2)+L-THIRD_Sy/dx+Width,A/2+THIRD_Sx/dy,RADIUS3);
			msiftcircle(((B*place)/2)+L-THIRD_Sy/dx,A*3/2+THIRD_Sx/dy,RADIUSd);                  //変更take
		}
//シフト列1列目(=偶数)のときの円孔生成
		else if(place == unitcell_number-1){
			msiftcircle(((B*place)/2)+L+FIRST_Sy/dx,A/2+FIRST_Sx/dy,RADIUS);
			msiftcircle(((B*place)/2)+L+FIRST_Sy/dx,A*3/2+FIRST_Sx/dy,RADIUSd);                  //変更take
		}
		else if(place == unitcell_number+1){
			msiftcircle(((B*place)/2)+L-FIRST_Sy/dx+Width,A/2+FIRST_Sx/dy,RADIUS);
			msiftcircle(((B*place)/2)+L-FIRST_Sy/dx+Width,A*3/2+FIRST_Sx/dy,RADIUSd);                  //変更take
		}
//シフト列5列目(=偶数)のときの円孔生成
		else if(place == unitcell_number-5){
			msiftcircle(((B*place)/2)+L+FIFTH_Sy/dx,A/2+FIFTH_Sx/dy,RADIUS5);
			msiftcircle(((B*place)/2)+L+FIFTH_Sy/dx,3*A/2+FIFTH_Sx/dy,RADIUSd);                  //変更take
		}
		else if(place == unitcell_number+5){
			msiftcircle(((B*place)/2)+L-FIFTH_Sy/dx+Width,A/2+FIFTH_Sx/dy,RADIUS5);
			msiftcircle(((B*place)/2)+L-FIFTH_Sy/dx+Width,A*3/2+FIFTH_Sx/dy,RADIUSd);                  //変更take
		}/*
		else if(place == unitcell_number+10){
			mtangle1(((B*place)/2)+L-SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,LENGTH_L);
//			mtangle1(((B*place)/2)+L-SEVENTH_Sy/dx,A*3/2+SEVENTH_Sx/dy,LENGTH_L);
		}
		else if(place == unitcell_number-10){
			mtangle2(((B*place)/2)+L+SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,LENGTH_L);
//			mtangle2(((B*place)/2)+L+SEVENTH_Sy/dx,3*A/2+SEVENTH_Sx/dy,LENGTH_L);
		}*/
//シフト列7列目(=偶数)のときの円孔生成
		else if(place == unitcell_number-7){
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //変更take
		}
		else if(place == unitcell_number+7){
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //変更take
		}
//シフト列9列目(=偶数)のときの円孔生成
		else if(place == unitcell_number-9){
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //変更take
		}
		else if(place == unitcell_number+9){
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //変更take
		}
//シフト列11列目(=偶数)のときの円孔生成
		else if(place == unitcell_number-11){
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //変更take
		}
		else if(place == unitcell_number+11){
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //変更take
		}
//シフト列ex1,ex2列目のときの円孔生成
		else if(place == unitcell_number-4){
			msiftcircle(((B*place)/2)+L+Ex1_Sy/dx,A/2+Ex1_Sx/dy,RADIUSex1);
			msiftcircle(((B*place)/2)+L+Ex2_Sy/dx,A/2+Ex2_Sx/dy,RADIUSex2);
		}
		else if(place == unitcell_number+4){
			msiftcircle(((B*place)/2)+L-Ex1_Sy/dx+Width,A/2+Ex1_Sx/dy,RADIUSex1);
			msiftcircle(((B*place)/2)+L-Ex2_Sy/dx+Width,A/2+Ex2_Sx/dy,RADIUSex2);
		}




#if PC_WAVEGUIDE==0		// 単位セル
		if(place%2==1){
			halfcircle1(B/2,0,RADIUS);
			halfcircle2(B/2,A,RADIUS);
#endif
#if PC_WAVEGUIDE==1		// 1列欠損導波路
		if(place%2==1 && place!=unitcell_number){
//			halfcircle1(((B*place)/2)+L,0,RADIUS);
//			halfcircle2(((B*place)/2)+L,A,RADIUS);


//変更場所(始め)
//******変更した*******(ここから)
//奇数列かつ1,3列目でないとき
				if(place == unitcell_number-2){
					msiftcircle(((B*place)/2)+L+SECOND_Sy/dx,A+SECOND_Sx/dy,RADIUSZ);
					msiftcircle(((B*place)/2)+L+SECOND_Sy/dx,A*2+SECOND_Sx/dy,RADIUSd);                  //変更take
				}
				else if(place == unitcell_number+2){
					msiftcircle(((B*place)/2)+L-SECOND_Sy/dx+Width,A+SECOND_Sx/dy,RADIUSZ);
					msiftcircle(((B*place)/2)+L-SECOND_Sy/dx+Width,2*A+SECOND_Sx/dy,RADIUSd);                  //変更take
				}
				else if(place == unitcell_number-4){
					msiftcircle(((B*place)/2)+L+FOURTH_Sy/dx,A+FOURTH_Sx/dy,RADIUS4);
					msiftcircle(((B*place)/2)+L+FOURTH_Sy/dx,2*A+FOURTH_Sx/dy,RADIUSd);                  //変更take
				}
				else if(place == unitcell_number+4){
					msiftcircle(((B*place)/2)+L-FOURTH_Sy/dx+Width,A+FOURTH_Sx/dy,RADIUS4);
					msiftcircle(((B*place)/2)+L-FOURTH_Sy/dx+Width,2*A+FOURTH_Sx/dy,RADIUSd);                  //変更take
				}
				else if(place == unitcell_number-6){
					msiftcircle(((B*place)/2)+L+SIXTH_Sy/dx,A+SIXTH_Sx/dy,RADIUS6);
					msiftcircle(((B*place)/2)+L+SIXTH_Sy/dx,2*A+SIXTH_Sx/dy,RADIUSd);                  //変更take
				}
				else if(place == unitcell_number+6){
					msiftcircle(((B*place)/2)+L-SIXTH_Sy/dx+Width,A+SIXTH_Sx/dy,RADIUS6);
					msiftcircle(((B*place)/2)+L-SIXTH_Sy/dx+Width,2*A+SIXTH_Sx/dy,RADIUSd);                  //変更take
				}
				else if(place == unitcell_number-8){
					msiftcircle(((B*place)/2)+L+SIXTH_Sy/dx,A+SIXTH_Sx/dy,RADIUS6);
					msiftcircle(((B*place)/2)+L+SIXTH_Sy/dx,2*A+SIXTH_Sx/dy,RADIUSd);                  //変更take
				}
				else if(place == unitcell_number+8){
					msiftcircle(((B*place)/2)+L-SIXTH_Sy/dx+Width,A+SIXTH_Sx/dy,RADIUS6);
					msiftcircle(((B*place)/2)+L-SIXTH_Sy/dx+Width,2*A+SIXTH_Sx/dy,RADIUSd);                  //変更take
				}
				else{
					if(place<unitcell_number){
						halfcircle1(((B*place)/2)+L,0,RADIUSd);
						mcircle(((B*place)/2)+L,A,RADIUS);                 //変更take
						halfcircle2(((B*place)/2)+L,2*A,RADIUSd);                  //変更take
					}
					if(place>unitcell_number){
						halfcircle1(((B*place)/2)+L+Width,0,RADIUSd);
						mcircle(((B*place)/2)+L+Width,A,RADIUS6);                  //変更take
						halfcircle2(((B*place)/2)+L+Width,2*A,RADIUSd);                  //変更take
					}
				}
//******変更した*******(ここまで)

				
//******変更した*******(ここから)
/*
			if(place<unitcell_number-2){
				halfcircle1(((B*place)/2)+L,0,RADIUS);
				halfcircle2(((B*place)/2)+L,A,RADIUS);
			}
			if(place>unitcell_number+2){
				halfcircle1(((B*place)/2)+L+Width,0,RADIUS);
				halfcircle2(((B*place)/2)+L+Width,A,RADIUS);
			}
			if(place==unitcell_number-2){
				halfcircle1(((B*place)/2)+L,0,RADIUS4);
				halfcircle2(((B*place)/2)+L,A,RADIUS4);
			}
			if(place==unitcell_number+2){
				halfcircle1(((B*place)/2)+L+Width,0,RADIUS4);
				halfcircle2(((B*place)/2)+L+Width,A,RADIUS4);
			}
*/
//******変更した*******(ここまで)
//変更場所(終り)

#endif
		}

//変更場所v5.10(始め)
#if PC_WAVEGUIDE==1		// 1列欠損導波路
		if(place==unitcell_number){
			mcircle(((B*place)/2)+L+Width/2,A/2,RADIUS2);		//導波路中央位相シフト円孔
		}
#endif
//変更場所(終り)
	
	}
	printf("\n");
}
/********************************** モデルの定義終わり **********************************/


/*********************************** 円孔の定義と生成 ***********************************/
void mcircle(int x, int y, double r){
	rightquartercircle1(x,y,r);
	leftquartercircle1(x-1,y,r);
	rightquartercircle2(x,y-1,r);
	leftquartercircle2(x-1,y-1,r);
	printf("[%4d %4d]",x,y);
}
//******変更した*******(ここから)
//シフトした円孔がはみ出したとき用
void msiftcircle(int x, int y, double r){
	rightsiftquartercircle1(x,y,r);
	leftsiftquartercircle1(x-1,y,r);
	rightsiftquartercircle2(x,y-1,r);
	leftsiftquartercircle2(x-1,y-1,r);
	printf("[%4d %4d]",x,y);
}

void rightsiftquartercircle1(int x, int y, double RAD){
	int i,j,Ie,Je;
	double R=((RAD*5.0e9)/(dx*5.0e9));	//5.16
	double r;							//5.16
	int circle_cell=0;
//	int R=(int)((RAD*5.0e9)/(dx*5.0e9));
//	int r;
	Ie=(x+R-1);
	Je=(y+R-1);
	for(i=x;i<=Ie;i++){
		for(j=y;j<=Je;j++){
//			for(k=(ck-(C/2));k<(ck+(C/2));k++){

//変更場所(始め)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//変更場所(終り)

//変更場所v5.13(始め)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//変更場所(終り)

				r=sqrt(double((i-x+1)*(i-x+1)+(j-y+1)*(j-y+1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//変更場所(始め)
				if(r<=R){	//5.16
					if(j>2*A-1){                                                           //変更take
						nn[i][j-2*A][k]=CIRCLE_MEDIUM;
						cnstE[i][j-2*A][k]=Air_E;
					}else{
						nn[i][j][k]=CIRCLE_MEDIUM;
						cnstE[i][j][k]=Air_E;
					}
					if (k==a) circle_cell++;
				}
//変更場所(終り)

			}
		}
	}
	circle_cell*=4;
	printf ("\ncircle_cell = %d\t 2r/a = %f\n", circle_cell, sqrt((double)circle_cell/PI)/23.0);
	effective_circle_cell = circle_cell;
	effective_2r_a = 2 * sqrt((double)circle_cell * dx * dx / PI) / PITCH;
}
void rightsiftquartercircle2(int x, int y, double RAD){
	int i,j,Ie,Je;
	double R=((RAD*5.0e9)/(dx*5.0e9));	//5.16
	double r;							//5.16
	int circle_cell=0;
	Ie=(x+R-1);
	Je=(y-R+1);
	for(i=x;i<=Ie;i++){
		for(j=y;j>=Je;j--){
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
				r=sqrt(double((i-x+1)*(i-x+1)+(j-y-1)*(j-y-1)))-0.5;	//5.16

			if(r<=R){	//5.16
					if(j>2*A-1){                                                           //変更take
						nn[i][j-2*A][k]=CIRCLE_MEDIUM;
						cnstE[i][j-2*A][k]=Air_E;
					}else{
						nn[i][j][k]=CIRCLE_MEDIUM;
						cnstE[i][j][k]=Air_E;
					}
					if (k==a) circle_cell++;
				}

			}
		}
	}
	circle_cell*=4;
	printf ("\ncircle_cell = %d\t 2r/a = %f\n", circle_cell, sqrt((double)circle_cell/PI)/23.0);
	effective_circle_cell = circle_cell;
	effective_2r_a = 2 * sqrt((double)circle_cell * dx * dx / PI) / PITCH;
}

void leftsiftquartercircle1(int x, int y, double RAD){
	int i,j,Ie,Je;
	double R=((RAD*5.0e9)/(dx*5.0e9));	//5.16
	double r;							//5.16
//	int R=(int)((RAD*5.0e9)/(dx*5.0e9));
//	int r;
	Ie=(x-R+1);
	Je=(y+R-1);
	for(i=x;i>=Ie;i--){
		for(j=y;j<=Je;j++){
//			for(k=(ck-(C/2));k<(ck+(C/2));k++){

//変更場所(始め)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//変更場所(終り)

//変更場所v5.13(始め)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//変更場所(終り)

				r = sqrt(double((i-x-1)*(i-x-1)+(j-y+1)*(j-y+1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//変更場所(始め)
				if(r<=R){	//5.16
					if(j>2*A-1){                                                           //変更take
						nn[i][j-2*A][k]=CIRCLE_MEDIUM;
						cnstE[i][j-2*A][k]=Air_E;
					}else{
						nn[i][j][k]=CIRCLE_MEDIUM;
						cnstE[i][j][k]=Air_E;
					}
				}
//変更場所(終り)

			}
		}
	}
}

void leftsiftquartercircle2(int x, int y, double RAD){
	int i,j,Ie,Je;
	double R=((RAD*5.0e9)/(dx*5.0e9));	//5.16
	double r;							//5.16
	Ie=(x-R+1);
	Je=(y-R+1);
	for(i=x;i>=Ie;i--){
		for(j=y;j>=Je;j--){
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;

				r = sqrt(double((i-x-1)*(i-x-1)+(j-y-1)*(j-y-1)))-0.5;	//5.16
				if(r<=R){	//5.16
					if(j>2*A-1){                                                           //変更take
						nn[i][j-2*A][k]=CIRCLE_MEDIUM;
						cnstE[i][j-2*A][k]=Air_E;
					}else{
						nn[i][j][k]=CIRCLE_MEDIUM;
						cnstE[i][j][k]=Air_E;
					}

				}

			}
		}
	}
}

void mtangle1(int x, int y, double RAD){
	int i,j,Ie,Je;	
	double R=((RAD*5.0e9)/(dy*5.0e9));					
	Ie=imax;
	Je=(y+R-1);
	for(i=x;i<=Ie;i++){
		for(j=y;j<=Je;j++){
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
				if(j>2*A-1){                                                           //変更take
					nn[i][j-2*A][k]=CIRCLE_MEDIUM;
					cnstE[i][j-2*A][k]=Air_E;
				}else{
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
				}
			}
		}

	}
}

void mtangle2(int x, int y, double RAD){
	int i,j,Ie,Je;	
	double R=((RAD*5.0e9)/(dy*5.0e9));					
	Ie=1;
	Je=(y+R-1);
	for(i=x;i>=Ie;i--){
		for(j=y;j<=Je;j++){
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
				if(j>2*A-1){                                                           //変更take
					nn[i][j-2*A][k]=CIRCLE_MEDIUM;
					cnstE[i][j-2*A][k]=Air_E;
				}else{
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
				}
			}
		}

	}
}




//******変更した*******(ここまで)

void halfcircle1(int x, int y, double r){
	rightquartercircle1(x,y,r);
	leftquartercircle1(x-1,y,r);
	printf("[%4d %4d]",x,y);
}

void halfcircle2(int x, int y, double r){
	rightquartercircle2(x,y-1,r);
	leftquartercircle2(x-1,y-1,r);
	printf("[%4d %4d]",x,y);
}

void halfcircle3(int x, int y, double r){
	rightquartercircle1(x,y,r);
	rightquartercircle2(x,y-1,r);
	printf("[%4d %4d]",x,y);
}

void halfcircle4(int x, int y, double r){
	leftquartercircle1(x-1,y,r);
	leftquartercircle2(x-1,y-1,r);
	printf("[%4d %4d]",x,y);
}

void rightquartercircle1(int x, int y, double RAD){
	int i,j,Ie,Je;
	double R=((RAD*5.0e9)/(dx*5.0e9));	//5.16
	double r;							//5.16
	int circle_cell=0;
//	int R=(int)((RAD*5.0e9)/(dx*5.0e9));
//	int r;
	Ie=(x+R-1);
	Je=(y+R-1);
	for(i=x;i<=Ie;i++){
		for(j=y;j<=Je;j++){
//			for(k=(ck-(C/2));k<(ck+(C/2));k++){

//変更場所(始め)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//変更場所(終り)

//変更場所v5.13(始め)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//変更場所(終り)

				r=sqrt(double((i-x+1)*(i-x+1)+(j-y+1)*(j-y+1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//変更場所(始め)
				if(r<=R){	//5.16
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
					if (k==a) circle_cell++;
				}
//変更場所(終り)

			}
		}
	}
	circle_cell*=4;
	printf ("\ncircle_cell = %d\t 2r/a = %f\n", circle_cell, sqrt((double)circle_cell/PI)/23.0);
}

void leftquartercircle1(int x, int y, double RAD){
	int i,j,Ie,Je;
	double R=((RAD*5.0e9)/(dx*5.0e9));	//5.16
	double r;							//5.16
//	int R=(int)((RAD*5.0e9)/(dx*5.0e9));
//	int r;
	Ie=(x-R+1);
	Je=(y+R-1);
	for(i=x;i>=Ie;i--){
		for(j=y;j<=Je;j++){
//			for(k=(ck-(C/2));k<(ck+(C/2));k++){

//変更場所(始め)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//変更場所(終り)

//変更場所v5.13(始め)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//変更場所(終り)

				r = sqrt(double((i-x-1)*(i-x-1)+(j-y+1)*(j-y+1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//変更場所(始め)
				if(r<=R){	//5.16
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
				}
//変更場所(終り)

			}
		}
	}
}

void rightquartercircle2(int x, int y, double RAD){
	int i,j,Ie,Je;
	double R=((RAD*5.0e9)/(dx*5.0e9));	//5.16
	double r;							//5.16
//	int R=(int)((RAD*5.0e9)/(dx*5.0e9));
//	int r;
	Ie=(x+R-1);
	Je=(y-R+1);
	for(i=x;i<=Ie;i++){
		for(j=y;j>=Je;j--){
//			for(k=(ck-(C/2));k<(ck+(C/2));k++){

//変更場所(始め)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//変更場所(終り)

//変更場所v5.13(始め)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//変更場所(終り)

				r = sqrt(double((i-x+1)*(i-x+1)+(j-y-1)*(j-y-1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//変更場所(始め)
				if(r<=R){	//5.16
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
				}
//変更場所(終り)

			}
		}
	}
}

void leftquartercircle2(int x, int y, double RAD){
	int i,j,Ie,Je;
	double R=((RAD*5.0e9)/(dx*5.0e9));	//5.16
	double r;							//5.16
//	int R=(int)((RAD*5.0e9)/(dx*5.0e9));
//	int r;
	Ie=(x-R+1);
	Je=(y-R+1);
	for(i=x;i>=Ie;i--){
		for(j=y;j>=Je;j--){
//			for(k=(ck-(C/2));k<(ck+(C/2));k++){

//変更場所(始め)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//変更場所(終り)

//変更場所v5.13(始め)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//変更場所(終り)

				r = sqrt(double((i-x-1)*(i-x-1)+(j-y-1)*(j-y-1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//変更場所(始め)
				if(r<=R){	//5.16
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
				}
//変更場所(終り)

			}
		}
	}
}

/*
何故いちいち円孔を4分割して定義しているかというと，
FDTD計算ではそうしないと真円を作ることができないから．
*/
/******************************** 円孔の定義と生成終わり ********************************/


/************************************ 位相情報の計算 ************************************/
void k_vector_phase(){
	double kk0X,kk0J,GG,R1x,R1y,R2x,R2y,Rz,phasex,phasey,phasexy;

	R1x=imax*dx;		R1y=0;
	R2x=0;				R2y=jmax*dy;
	Rz=0;

	GG=2.0*PI/((2.0*(PITCH/2)*sqrt(3.0))/2);

	G1x=2.0*PI/R1x;		G1y=0;
	G2x=0;				G2y=2.0*PI/R2y;
	Gz=0;

	kk0J=G2y/2;			kk0X=G1x/2;

	if(m==1){	// G-X
		kkx=G1x*ikk/IKK_MAX;
		kky=0;
	}
	if(m==2){	// G-J
		kkx=0;
		kky=(2*G2y/3)*ikk/IKK_MAX;
	}
	if(m==3){	// X-J
		kkx=G1x;
		kky=(G1x/sqrt(3.0))*ikk/IKK_MAX;
	}
	if(m==4){	// G-P(導波路計算の場合)
		kkx=0;
		kky=kk0J*ikk/IKK_MAX;
	}
	kkz=0;

	phasex=kkx*R1x+kky*R1y;
	phasey=kkx*R2x+kky*R2y;
	phasexy=kkx*R1x+kky*R2y;

	cos_phx=cos(phasex);
	sin_phx=sin(phasex);
	cos_phy=cos(phasey);
	sin_phy=sin(phasey);
	cos_phxy=cos(phasexy);
	sin_phxy=sin(phasexy);
}
/********************************* 位相情報の計算終わり *********************************/


/*************************************** 励振関数 ***************************************/
void exc_function(){
	double exc_cos,exc_sin;
	exc_cos=cos(omega0*dt*(n-PULSEPEAK))*exp(-((sigma*dt*(n-PULSEPEAK))*(sigma*dt*(n-PULSEPEAK))/2));
	exc_sin=sin(omega0*dt*(n-PULSEPEAK))*exp(-((sigma*dt*(n-PULSEPEAK))*(sigma*dt*(n-PULSEPEAK))/2));

	if(excite==1){
		i=5;j=6;k=3;

//変更場所(始め)
//		i=0;j=0;k=0;
//変更場所(終り)

		re_hz[ci+i][cj+j][ck+k]+=exc_cos;		im_hz[ci+i][cj+j][ck+k]+=exc_sin;
		re_ex[ci+i][cj+j][ck+k]+=exc_cos;		im_ex[ci+i][cj+j][ck+k]+=exc_sin;
		re_ey[ci+i][cj+j][ck+k]+=exc_cos;		im_ey[ci+i][cj+j][ck+k]+=exc_sin;
		i=-2;j=-4;k=-1;
		re_hz[ci+i][cj+j][ck+k]+=exc_cos;		im_hz[ci+i][cj+j][ck+k]+=exc_sin;
		re_ex[ci+i][cj+j][ck+k]+=exc_cos;		im_ex[ci+i][cj+j][ck+k]+=exc_sin;
		re_ey[ci+i][cj+j][ck+k]+=exc_cos;		im_ey[ci+i][cj+j][ck+k]+=exc_sin;
		i=-4;j=3;k=-2;
		re_hz[ci+i][cj+j][ck+k]+=exc_cos;		im_hz[ci+i][cj+j][ck+k]+=exc_sin;
		re_ex[ci+i][cj+j][ck+k]+=exc_cos;		im_ex[ci+i][cj+j][ck+k]+=exc_sin;
		re_ey[ci+i][cj+j][ck+k]+=exc_cos;		im_ey[ci+i][cj+j][ck+k]+=exc_sin;
	}
	if(excite==2){
		i=3;j=3;k=3;
		re_hx[ci+i][cj+j][ck+k]+=exc_cos;		im_hx[ci+i][cj+j][ck+k]+=exc_sin;
		re_hy[ci+i][cj+j][ck+k]+=exc_cos;		im_hy[ci+i][cj+j][ck+k]+=exc_sin;
		re_ez[ci+i][cj+j][ck+k]+=exc_cos;		im_ez[ci+i][cj+j][ck+k]+=exc_sin;
		i=-2;j=-3;k=-5;
		re_hx[ci+i][cj+j][ck+k]+=exc_cos;		im_hx[ci+i][cj+j][ck+k]+=exc_sin;
		re_hy[ci+i][cj+j][ck+k]+=exc_cos;		im_hy[ci+i][cj+j][ck+k]+=exc_sin;
		re_ez[ci+i][cj+j][ck+k]+=exc_cos;		im_ez[ci+i][cj+j][ck+k]+=exc_sin;
		i=1;j=2;k=-2;
		re_hx[ci+i][cj+j][ck+k]+=exc_cos;		im_hx[ci+i][cj+j][ck+k]+=exc_sin;
		re_hy[ci+i][cj+j][ck+k]+=exc_cos;		im_hy[ci+i][cj+j][ck+k]+=exc_sin;
		re_ez[ci+i][cj+j][ck+k]+=exc_cos;		im_ez[ci+i][cj+j][ck+k]+=exc_sin;
	}
}
/*
励振は中心座標(ci,cj,ck)から(i,j,k)ずらした位置で与えている．
つまり，モデルの設定範囲を超えないように注意すること．
*/
/************************************ 励振関数終わり ************************************/


/****************************** 電磁界の計算＆周期境界条件 ******************************/
void calc_hfield(){
	double re_dhx,re_dhy,re_dhz;
	double im_dhx,im_dhy,im_dhz;

	for(i=L;i<=imax-L;i++){
		for(j=0;j<jmax;j++){
			for(k=L;k<kmax-L;k++){
				re_dhx=((re_ey[i][j][k+1]-re_ey[i][j][k])/dz)-((re_ez[i][j+1][k]-re_ez[i][j][k])/dy);
				im_dhx=((im_ey[i][j][k+1]-im_ey[i][j][k])/dz)-((im_ez[i][j+1][k]-im_ez[i][j][k])/dy);
				re_hx[i][j][k]=re_hx[i][j][k]+cnstM*re_dhx;
				im_hx[i][j][k]=im_hx[i][j][k]+cnstM*im_dhx;
			}
		}
	}

	for(i=L;i<imax-L;i++){
		for(j=0;j<=jmax;j++){
			for(k=L;k<kmax-L;k++){
				re_dhy=((re_ez[i+1][j][k]-re_ez[i][j][k])/dx)-((re_ex[i][j][k+1]-re_ex[i][j][k])/dz);
				im_dhy=((im_ez[i+1][j][k]-im_ez[i][j][k])/dx)-((im_ex[i][j][k+1]-im_ex[i][j][k])/dz);
				re_hy[i][j][k]=re_hy[i][j][k]+cnstM*re_dhy;
				im_hy[i][j][k]=im_hy[i][j][k]+cnstM*im_dhy;
			}
		}
	}


	for(i=L;i<imax-L;i++){
		for(j=0;j<jmax;j++){
			for(k=L;k<=kmax-L;k++){
				re_dhz=((re_ex[i][j+1][k]-re_ex[i][j][k])/dy)-((re_ey[i+1][j][k]-re_ey[i][j][k])/dx);
				im_dhz=((im_ex[i][j+1][k]-im_ex[i][j][k])/dy)-((im_ey[i+1][j][k]-im_ey[i][j][k])/dx);
				re_hz[i][j][k]=re_hz[i][j][k]+cnstM*re_dhz;
				im_hz[i][j][k]=im_hz[i][j][k]+cnstM*im_dhz;
			}
		}
	}
}

void calc_h_periodic_boundary_condition(){
	for(i=1;i<imax;i++){
		for(k=1;k<kmax-1;k++){
			re_hx_jmax[i][k]=re_hx[i][0][k]*cos_phy+im_hx[i][0][k]*sin_phy;
			im_hx_jmax[i][k]=im_hx[i][0][k]*cos_phy-re_hx[i][0][k]*sin_phy;
			re_hx[i][jmax][k]=re_hx[i][0][k]*cos_phy+im_hx[i][0][k]*sin_phy;
			im_hx[i][jmax][k]=im_hx[i][0][k]*cos_phy-re_hx[i][0][k]*sin_phy;
		}
	}
	for(i=1;i<imax-1;i++){
		for(k=1;k<kmax;k++){
			re_hz_jmax[i][k]=re_hz[i][0][k]*cos_phy+im_hz[i][0][k]*sin_phy;
			im_hz_jmax[i][k]=im_hz[i][0][k]*cos_phy-re_hz[i][0][k]*sin_phy;
			re_hz[i][jmax][k]=re_hz[i][0][k]*cos_phy+im_hz[i][0][k]*sin_phy;
			im_hz[i][jmax][k]=im_hz[i][0][k]*cos_phy-re_hz[i][0][k]*sin_phy;
		}
	}
}

void calc_efield(){
	double re_dex,re_dey,re_dez;
	double im_dex,im_dey,im_dez;

	for(i=L;i<imax-L;i++){
		for(j=1;j<jmax;j++){
//			for(k=L;k<kmax-L;k++){

//変更場所(始め)
			for(k=L;k<=kmax-L;k++){
//変更場所(終り)

				re_dex=((re_hz[i][j][k]-re_hz[i][j-1][k])/dy)-((re_hy[i][j][k]-re_hy[i][j][k-1])/dz);
				im_dex=((im_hz[i][j][k]-im_hz[i][j-1][k])/dy)-((im_hy[i][j][k]-im_hy[i][j][k-1])/dz);
//				switch(nn[i][j][k]){
//					case AIR:
//						re_ex[i][j][k]=re_ex[i][j][k]+cnstE*re_dex;
//						im_ex[i][j][k]=im_ex[i][j][k]+cnstE*im_dex;
//						break;
//					case SI:
//						re_ex[i][j][k]=re_ex[i][j][k]+cnstE1*re_dex;
//						im_ex[i][j][k]=im_ex[i][j][k]+cnstE1*im_dex;
//						break;
//					case SIO2:
//						re_ex[i][j][k]=re_ex[i][j][k]+cnstE2*re_dex;
//						im_ex[i][j][k]=im_ex[i][j][k]+cnstE2*im_dex;
//						break;
//				}

//変更場所(始め)
				re_ex[i][j][k]=re_ex[i][j][k]+cnstE[i][j][k]*re_dex;
				im_ex[i][j][k]=im_ex[i][j][k]+cnstE[i][j][k]*im_dex;
//変更場所(終り)

			}
		}
	}

//	for(i=L;i<imax-L;i++){

//変更場所(始め)
	for(i=L;i<=imax-L;i++){
//変更場所(終り)

		for(j=0;j<jmax;j++){
//			for(k=L;k<kmax-L;k++){

//変更場所(始め)
			for(k=L;k<=kmax-L;k++){
//変更場所(終り)

				re_dey=((re_hx[i][j][k]-re_hx[i][j][k-1])/dz)-((re_hz[i][j][k]-re_hz[i-1][j][k])/dx);
				im_dey=((im_hx[i][j][k]-im_hx[i][j][k-1])/dz)-((im_hz[i][j][k]-im_hz[i-1][j][k])/dx);
//				switch(nn[i][j][k]){
//					case AIR:
//						re_ey[i][j][k]=re_ey[i][j][k]+cnstE*re_dey;
//						im_ey[i][j][k]=im_ey[i][j][k]+cnstE*im_dey;
//						break;
//					case SI:
//						re_ey[i][j][k]=re_ey[i][j][k]+cnstE1*re_dey;
//						im_ey[i][j][k]=im_ey[i][j][k]+cnstE1*im_dey;
//						break;
//					case SIO2:
//						re_ey[i][j][k]=re_ey[i][j][k]+cnstE2*re_dey;
//						im_ey[i][j][k]=im_ey[i][j][k]+cnstE2*im_dey;
//						break;
//				}

//変更場所(始め)
				re_ey[i][j][k]=re_ey[i][j][k]+cnstE[i][j][k]*re_dey;
				im_ey[i][j][k]=im_ey[i][j][k]+cnstE[i][j][k]*im_dey;
//変更場所(終り)

			}
		}
	}

//	for(i=L;i<imax-L;i++){

//変更場所(始め)
	for(i=L;i<=imax-L;i++){
//変更場所(終り)

		for(j=1;j<jmax;j++){
			for(k=L;k<kmax-L;k++){
				re_dez=((re_hy[i][j][k]-re_hy[i-1][j][k])/dx)-((re_hx[i][j][k]-re_hx[i][j-1][k])/dy);
				im_dez=((im_hy[i][j][k]-im_hy[i-1][j][k])/dx)-((im_hx[i][j][k]-im_hx[i][j-1][k])/dy);
//				switch(nn[i][j][k]){
//					case AIR:
//						re_ez[i][j][k]=re_ez[i][j][k]+cnstE*re_dez;
//						im_ez[i][j][k]=im_ez[i][j][k]+cnstE*im_dez;
//						break;
//					case SI:
//						re_ez[i][j][k]=re_ez[i][j][k]+cnstE1*re_dez;
//						im_ez[i][j][k]=im_ez[i][j][k]+cnstE1*im_dez;
//						break;
//					case SIO2:
//						re_ez[i][j][k]=re_ez[i][j][k]+cnstE2*re_dez;
//						im_ez[i][j][k]=im_ez[i][j][k]+cnstE2*im_dez;
//						break;
//				}

//変更場所(始め)
				re_ez[i][j][k]=re_ez[i][j][k]+cnstE[i][j][k]*re_dez;
				im_ez[i][j][k]=im_ez[i][j][k]+cnstE[i][j][k]*im_dez;
//変更場所(終り)

			}
		}
	}
}

void calc_e_periodic_boundary_condition(){
	for(i=1;i<imax-1;i++){
		for(k=1;k<kmax;k++){
			re_ex[i][0][k]=re_ex[i][jmax][k]*cos_phy-im_ex[i][jmax][k]*sin_phy;
			im_ex[i][0][k]=im_ex[i][jmax][k]*cos_phy+re_ex[i][jmax][k]*sin_phy;
		}
	}
	for(i=1;i<imax;i++){
		for(k=1;k<kmax-1;k++){
			re_ez[i][0][k]=re_ez[i][jmax][k]*cos_phy-im_ez[i][jmax][k]*sin_phy;
			im_ez[i][0][k]=im_ez[i][jmax][k]*cos_phy+re_ez[i][jmax][k]*sin_phy;
		}
	}
}

void calc_e_periodic_boundary_condition2(){
	double re_dex,re_dez;
	double im_dex,im_dez;

	for(i=L;i<imax-L;i++){
//		for(k=L;k<kmax-L;k++){

//変更場所(始め)
		for(k=L;k<=kmax-L;k++){
//変更場所(終り)

			re_dex=((re_hz_jmax[i][k]-re_hz[i][jmax-1][k])/dy)-((re_hy[i][jmax][k]-re_hy[i][jmax][k-1])/dz);
			im_dex=((im_hz_jmax[i][k]-im_hz[i][jmax-1][k])/dy)-((im_hy[i][jmax][k]-im_hy[i][jmax][k-1])/dz);
//			switch(nn[i][jmax][k]){
//				case AIR:
//					re_ex[i][jmax][k]=re_ex[i][jmax][k]+cnstE*re_dex;
//					im_ex[i][jmax][k]=im_ex[i][jmax][k]+cnstE*im_dex;
//					break;
//				case SI:
//					re_ex[i][jmax][k]=re_ex[i][jmax][k]+cnstE1*re_dex;
//					im_ex[i][jmax][k]=im_ex[i][jmax][k]+cnstE1*im_dex;
//					break;
//				case SIO2:
//					re_ex[i][jmax][k]=re_ex[i][jmax][k]+cnstE2*re_dex;
//					im_ex[i][jmax][k]=im_ex[i][jmax][k]+cnstE2*im_dex;
//					break;
//			}

//変更場所(始め)
			re_ex[i][jmax][k]=re_ex[i][jmax][k]+cnstE[i][jmax][k]*re_dex;
			im_ex[i][jmax][k]=im_ex[i][jmax][k]+cnstE[i][jmax][k]*im_dex;
//変更場所(終り)

		}
	}

//	for(i=L;i<imax-L;i++){

//変更場所(始め)
	for(i=L;i<=imax-L;i++){
//変更場所(終り)

		for(k=L;k<kmax-L;k++){
			re_dez=((re_hy[i][jmax][k]-re_hy[i-1][jmax][k])/dx)-((re_hx_jmax[i][k]-re_hx[i][jmax-1][k])/dy);
			im_dez=((im_hy[i][jmax][k]-im_hy[i-1][jmax][k])/dx)-((im_hx_jmax[i][k]-im_hx[i][jmax-1][k])/dy);
//			switch(nn[i][jmax][k]){
//				case AIR:
//					re_ez[i][jmax][k]=re_ez[i][jmax][k]+cnstE*re_dez;
//					im_ez[i][jmax][k]=im_ez[i][jmax][k]+cnstE*im_dez;
//					break;
//				case SI:
//					re_ez[i][jmax][k]=re_ez[i][jmax][k]+cnstE1*re_dez;
//					im_ez[i][jmax][k]=im_ez[i][jmax][k]+cnstE1*im_dez;
//					break;
//				case SIO2:
//					re_ez[i][jmax][k]=re_ez[i][jmax][k]+cnstE2*re_dez;
//					im_ez[i][jmax][k]=im_ez[i][jmax][k]+cnstE2*im_dez;
//					break;
//			}

//変更場所(始め)
			re_ez[i][jmax][k]=re_ez[i][jmax][k]+cnstE[i][jmax][k]*re_dez;
			im_ez[i][jmax][k]=im_ez[i][jmax][k]+cnstE[i][jmax][k]*im_dez;
//変更場所(終り)

		}
	}
}
/*************************** 電磁界の計算＆周期境界条件終わり ***************************/


/************************************* 吸収境界条件 *************************************/
void efield_PML(){

/*********** Ex-PML ***********/
	sigma_x=0.0;sigma_y=0.0;sigma_z=0.0;
	u1=0.0;u2=0.0;term1=0.0;term2=0.0;

	for(i=1;i<L;i++){
		for(j=1;j<=jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(010)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz010[i][j][k]=term1*re_exz010[i][j][k]-term2*(re_hy[i][j][k]-re_hy[i][j][k-1])/dz;
				im_exz010[i][j][k]=term1*im_exz010[i][j][k]-term2*(im_hy[i][j][k]-im_hy[i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exy010[i][j][k]=term1*re_exy010[i][j][k]+term2*(re_hz[i][j][k]-re_hz[i][j-1][k])/dy;
				im_exy010[i][j][k]=term1*im_exy010[i][j][k]+term2*(im_hz[i][j][k]-im_hz[i][j-1][k])/dy;

				//Exz+Exy
				re_ex[i][j][k]=re_exz010[i][j][k]+re_exy010[i][j][k];
				im_ex[i][j][k]=im_exz010[i][j][k]+im_exy010[i][j][k];

			/*(210)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz210[i][j][k]=term1*re_exz210[i][j][k]-term2*(re_hy[imax-1-i][j][k]-re_hy[imax-1-i][j][k-1])/dz;
				im_exz210[i][j][k]=term1*im_exz210[i][j][k]-term2*(im_hy[imax-1-i][j][k]-im_hy[imax-1-i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exy210[i][j][k]=term1*re_exy210[i][j][k]+term2*(re_hz[imax-1-i][j][k]-re_hz[imax-1-i][j-1][k])/dy;
				im_exy210[i][j][k]=term1*im_exy210[i][j][k]+term2*(im_hz[imax-1-i][j][k]-im_hz[imax-1-i][j-1][k])/dy;

				//Exz+Exy
				re_ex[imax-1-i][j][k]=re_exz210[i][j][k]+re_exy210[i][j][k];
				im_ex[imax-1-i][j][k]=im_exz210[i][j][k]+im_exy210[i][j][k];

			/*(012)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz012[i][j][k]=term1*re_exz012[i][j][k]-term2*(re_hy[i][j][kmax-k]-re_hy[i][j][kmax-k-1])/dz;
				im_exz012[i][j][k]=term1*im_exz012[i][j][k]-term2*(im_hy[i][j][kmax-k]-im_hy[i][j][kmax-k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exy012[i][j][k]=term1*re_exy012[i][j][k]+term2*(re_hz[i][j][kmax-k]-re_hz[i][j-1][kmax-k])/dy;
				im_exy012[i][j][k]=term1*im_exy012[i][j][k]+term2*(im_hz[i][j][kmax-k]-im_hz[i][j-1][kmax-k])/dy;

				//Exz+Exy		
				re_ex[i][j][kmax-k]=re_exz012[i][j][k]+re_exy012[i][j][k];
				im_ex[i][j][kmax-k]=im_exz012[i][j][k]+im_exy012[i][j][k];

			/*(212)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz212[i][j][k]=term1*re_exz212[i][j][k]-term2*(re_hy[imax-i-1][j][kmax-k]-re_hy[imax-i-1][j][kmax-k-1])/dz;
				im_exz212[i][j][k]=term1*im_exz212[i][j][k]-term2*(im_hy[imax-i-1][j][kmax-k]-im_hy[imax-i-1][j][kmax-k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exy212[i][j][k]=term1*re_exy212[i][j][k]+term2*(re_hz[imax-i-1][j][kmax-k]-re_hz[imax-i-1][j-1][kmax-k])/dy;
				im_exy212[i][j][k]=term1*im_exy212[i][j][k]+term2*(im_hz[imax-i-1][j][kmax-k]-im_hz[imax-i-1][j-1][kmax-k])/dy;

				//Exz+Exy
				re_ex[imax-i-1][j][kmax-k]=re_exz212[i][j][k]+re_exy212[i][j][k];
				im_ex[imax-i-1][j][kmax-k]=im_exz212[i][j][k]+im_exy212[i][j][k];
			}
		}
	}

	for(i=1;i<L;i++){
//		for(j=0;j<=jmax;j++){

//変更場所(始め)
		for(j=1;j<=jmax;j++){
//変更場所(終り)

			for(k=L;k<=kmax-L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=0.0;

			/*(011)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz011[i][j][k]=term1*re_exz011[i][j][k]-term2*(re_hy[i][j][k]-re_hy[i][j][k-1])/dz;
				im_exz011[i][j][k]=term1*im_exz011[i][j][k]-term2*(im_hy[i][j][k]-im_hy[i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exy011[i][j][k]=term1*re_exy011[i][j][k]+term2*(re_hz[i][j][k]-re_hz[i][j-1][k])/dy;
				im_exy011[i][j][k]=term1*im_exy011[i][j][k]+term2*(im_hz[i][j][k]-im_hz[i][j-1][k])/dy;

				//Exz+Exy
				re_ex[i][j][k]=re_exz011[i][j][k]+re_exy011[i][j][k];
				im_ex[i][j][k]=im_exz011[i][j][k]+im_exy011[i][j][k];
			
			/*(211)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz211[i][j][k]=term1*re_exz211[i][j][k]-term2*(re_hy[imax-1-i][j][k]-re_hy[imax-1-i][j][k-1])/dz;
				im_exz211[i][j][k]=term1*im_exz211[i][j][k]-term2*(im_hy[imax-1-i][j][k]-im_hy[imax-1-i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exy211[i][j][k]=term1*re_exy211[i][j][k]+term2*(re_hz[imax-1-i][j][k]-re_hz[imax-1-i][j-1][k])/dy;
				im_exy211[i][j][k]=term1*im_exy211[i][j][k]+term2*(im_hz[imax-1-i][j][k]-im_hz[imax-1-i][j-1][k])/dy;

				//Exz+Exy
				re_ex[imax-1-i][j][k]=re_exz211[i][j][k]+re_exy211[i][j][k];
				im_ex[imax-1-i][j][k]=im_exz211[i][j][k]+im_exy211[i][j][k];
			}
		}
	}

	for(i=L;i<imax-L;i++){
//		for(j=0;j<=jmax;j++){

//変更場所(始め)
		for(j=1;j<=jmax;j++){
//変更場所(終り)

			for(k=1;k<L;k++){

				sigma_x=0.0;
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(110)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz110[i][j][k]=term1*re_exz110[i][j][k]-term2*(re_hy[i][j][k]-re_hy[i][j][k-1])/dz;
				im_exz110[i][j][k]=term1*im_exz110[i][j][k]-term2*(im_hy[i][j][k]-im_hy[i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exy110[i][j][k]=term1*re_exy110[i][j][k]+term2*(re_hz[i][j][k]-re_hz[i][j-1][k])/dy;
				im_exy110[i][j][k]=term1*im_exy110[i][j][k]+term2*(im_hz[i][j][k]-im_hz[i][j-1][k])/dy;

				//Exz+Exy
				re_ex[i][j][k]=re_exz110[i][j][k]+re_exy110[i][j][k];
				im_ex[i][j][k]=im_exz110[i][j][k]+im_exy110[i][j][k];
			
			/*(112)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz112[i][j][k]=term1*re_exz112[i][j][k]-term2*(re_hy[i][j][kmax-k]-re_hy[i][j][kmax-k-1])/dz;
				im_exz112[i][j][k]=term1*im_exz112[i][j][k]-term2*(im_hy[i][j][kmax-k]-im_hy[i][j][kmax-k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exy112[i][j][k]=term1*re_exy112[i][j][k]+term2*(re_hz[i][j][kmax-k]-re_hz[i][j-1][kmax-k])/dy;
				im_exy112[i][j][k]=term1*im_exy112[i][j][k]+term2*(im_hz[i][j][kmax-k]-im_hz[i][j-1][kmax-k])/dy;

				//Exz+Exy
				re_ex[i][j][kmax-k]=re_exz112[i][j][k]+re_exy112[i][j][k];
				im_ex[i][j][kmax-k]=im_exz112[i][j][k]+im_exy112[i][j][k];
			}
		}
	}

/*********** Ey-PML ***********/ 
	sigma_x=0.0;sigma_y=0.0;sigma_z=0.0;
	u1=0.0;u2=0.0;term1=0.0;term2=0.0;

	for(i=1;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(010)*/
				//Eyx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx010[i][j][k]=term1*re_eyx010[i][j][k]-term2*(re_hz[i][j][k]-re_hz[i-1][j][k])/dx;
				im_eyx010[i][j][k]=term1*im_eyx010[i][j][k]-term2*(im_hz[i][j][k]-im_hz[i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyz010[i][j][k]=term1*re_eyz010[i][j][k]+term2*(re_hx[i][j][k]-re_hx[i][j][k-1])/dz;
				im_eyz010[i][j][k]=term1*im_eyz010[i][j][k]+term2*(im_hx[i][j][k]-im_hx[i][j][k-1])/dz;

				//Eyx+Eyz
				re_ey[i][j][k]=re_eyx010[i][j][k]+re_eyz010[i][j][k];
				im_ey[i][j][k]=im_eyx010[i][j][k]+im_eyz010[i][j][k];

			/*(210)*/
				//Eyx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx210[i][j][k]=term1*re_eyx210[i][j][k]-term2*(re_hz[imax-i][j][k]-re_hz[imax-i-1][j][k])/dx;
				im_eyx210[i][j][k]=term1*im_eyx210[i][j][k]-term2*(im_hz[imax-i][j][k]-im_hz[imax-i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyz210[i][j][k]=term1*re_eyz210[i][j][k]+term2*(re_hx[imax-i][j][k]-re_hx[imax-i][j][k-1])/dz;
				im_eyz210[i][j][k]=term1*im_eyz210[i][j][k]+term2*(im_hx[imax-i][j][k]-im_hx[imax-i][j][k-1])/dz;

				//Eyx+Eyz
				re_ey[imax-i][j][k]=re_eyx210[i][j][k]+re_eyz210[i][j][k];
				im_ey[imax-i][j][k]=im_eyx210[i][j][k]+im_eyz210[i][j][k];

			/*(012)*/
				//Eyx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx012[i][j][k]=term1*re_eyx012[i][j][k]-term2*(re_hz[i][j][kmax-k]-re_hz[i-1][j][kmax-k])/dx;
				im_eyx012[i][j][k]=term1*im_eyx012[i][j][k]-term2*(im_hz[i][j][kmax-k]-im_hz[i-1][j][kmax-k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyz012[i][j][k]=term1*re_eyz012[i][j][k]+term2*(re_hx[i][j][kmax-k]-re_hx[i][j][kmax-k-1])/dz;
				im_eyz012[i][j][k]=term1*im_eyz012[i][j][k]+term2*(im_hx[i][j][kmax-k]-im_hx[i][j][kmax-k-1])/dz;

				//Eyx+Eyz
				re_ey[i][j][kmax-k]=re_eyx012[i][j][k]+re_eyz012[i][j][k];
				im_ey[i][j][kmax-k]=im_eyx012[i][j][k]+im_eyz012[i][j][k];

			/*(212)*/
				//Eyx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx212[i][j][k]=term1*re_eyx212[i][j][k]-term2*(re_hz[imax-i][j][kmax-k]-re_hz[imax-i-1][j][kmax-k])/dx;
				im_eyx212[i][j][k]=term1*im_eyx212[i][j][k]-term2*(im_hz[imax-i][j][kmax-k]-im_hz[imax-i-1][j][kmax-k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyz212[i][j][k]=term1*re_eyz212[i][j][k]+term2*(re_hx[imax-i][j][kmax-k]-re_hx[imax-i][j][kmax-k-1])/dz;
				im_eyz212[i][j][k]=term1*im_eyz212[i][j][k]+term2*(im_hx[imax-i][j][kmax-k]-im_hx[imax-i][j][kmax-k-1])/dz;

				//Eyx+Eyz
				re_ey[imax-i][j][kmax-k]=re_eyx212[i][j][k]+re_eyz212[i][j][k];
				im_ey[imax-i][j][kmax-k]=im_eyx212[i][j][k]+im_eyz212[i][j][k];
			}
		}
	}

	for(i=1;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=L;k<=kmax-L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=0.0;

			/*(011)*/
				//Eyx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx011[i][j][k]=term1*re_eyx011[i][j][k]-term2*(re_hz[i][j][k]-re_hz[i-1][j][k])/dx;
				im_eyx011[i][j][k]=term1*im_eyx011[i][j][k]-term2*(im_hz[i][j][k]-im_hz[i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyz011[i][j][k]=term1*re_eyz011[i][j][k]+term2*(re_hx[i][j][k]-re_hx[i][j][k-1])/dz;
				im_eyz011[i][j][k]=term1*im_eyz011[i][j][k]+term2*(im_hx[i][j][k]-im_hx[i][j][k-1])/dz;

				//Eyx+Eyz
				re_ey[i][j][k]=re_eyx011[i][j][k]+re_eyz011[i][j][k];
				im_ey[i][j][k]=im_eyx011[i][j][k]+im_eyz011[i][j][k];

			/*(211)*/
				//Eyx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx211[i][j][k]=term1*re_eyx211[i][j][k]-term2*(re_hz[imax-i][j][k]-re_hz[imax-i-1][j][k])/dx;
				im_eyx211[i][j][k]=term1*im_eyx211[i][j][k]-term2*(im_hz[imax-i][j][k]-im_hz[imax-i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyz211[i][j][k]=term1*re_eyz211[i][j][k]+term2*(re_hx[imax-i][j][k]-re_hx[imax-i][j][k-1])/dz;
				im_eyz211[i][j][k]=term1*im_eyz211[i][j][k]+term2*(im_hx[imax-i][j][k]-im_hx[imax-i][j][k-1])/dz;

				//Eyx+Eyz
				re_ey[imax-i][j][k]=re_eyx211[i][j][k]+re_eyz211[i][j][k];
				im_ey[imax-i][j][k]=im_eyx211[i][j][k]+im_eyz211[i][j][k];
			}
		}
	}

	for(i=L;i<=imax-L;i++){
		for(j=0;j<jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=0.0;
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(110)*/
				//Eyx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx110[i][j][k]=term1*re_eyx110[i][j][k]-term2*(re_hz[i][j][k]-re_hz[i-1][j][k])/dx;
				im_eyx110[i][j][k]=term1*im_eyx110[i][j][k]-term2*(im_hz[i][j][k]-im_hz[i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyz110[i][j][k]=term1*re_eyz110[i][j][k]+term2*(re_hx[i][j][k]-re_hx[i][j][k-1])/dz;
				im_eyz110[i][j][k]=term1*im_eyz110[i][j][k]+term2*(im_hx[i][j][k]-im_hx[i][j][k-1])/dz;

				//Eyx+Eyz
				re_ey[i][j][k]=re_eyx110[i][j][k]+re_eyz110[i][j][k];
				im_ey[i][j][k]=im_eyx110[i][j][k]+im_eyz110[i][j][k];
			
			/*(112)*/
				//Eyx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx112[i][j][k]=term1*re_eyx112[i][j][k]-term2*(re_hz[i][j][kmax-k]-re_hz[i-1][j][kmax-k])/dx;
				im_eyx112[i][j][k]=term1*im_eyx112[i][j][k]-term2*(im_hz[i][j][kmax-k]-im_hz[i-1][j][kmax-k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyz112[i][j][k]=term1*re_eyz112[i][j][k]+term2*(re_hx[i][j][kmax-k]-re_hx[i][j][kmax-k-1])/dz;
				im_eyz112[i][j][k]=term1*im_eyz112[i][j][k]+term2*(im_hx[i][j][kmax-k]-im_hx[i][j][kmax-k-1])/dz;

				//Eyx+Eyz
				re_ey[i][j][kmax-k]=re_eyx112[i][j][k]+re_eyz112[i][j][k];
				im_ey[i][j][kmax-k]=im_eyx112[i][j][k]+im_eyz112[i][j][k];
			}
		}
	}

/*********** Ez-PML ***********/ 
	sigma_x=0.0;sigma_y=0.0;sigma_z=0.0;
	u1=0.0;u2=0.0;term1=0.0;term2=0.0;

	for(i=1;i<L;i++){
		for(j=1;j<=jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(010)*/
				//Ezx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx010[i][j][k]=term1*re_ezx010[i][j][k]+term2*(re_hy[i][j][k]-re_hy[i-1][j][k])/dx;
				im_ezx010[i][j][k]=term1*im_ezx010[i][j][k]+term2*(im_hy[i][j][k]-im_hy[i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezy010[i][j][k]=term1*re_ezy010[i][j][k]-term2*(re_hx[i][j][k]-re_hx[i][j-1][k])/dy;
				im_ezy010[i][j][k]=term1*im_ezy010[i][j][k]-term2*(im_hx[i][j][k]-im_hx[i][j-1][k])/dy;

				//Ezx+Ezy
				re_ez[i][j][k]=re_ezx010[i][j][k]+re_ezy010[i][j][k];
				im_ez[i][j][k]=im_ezx010[i][j][k]+im_ezy010[i][j][k];

			/*(210)*/
				//Ezx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx210[i][j][k]=term1*re_ezx210[i][j][k]+term2*(re_hy[imax-i][j][k]-re_hy[imax-i-1][j][k])/dx;
				im_ezx210[i][j][k]=term1*im_ezx210[i][j][k]+term2*(im_hy[imax-i][j][k]-im_hy[imax-i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezy210[i][j][k]=term1*re_ezy210[i][j][k]-term2*(re_hx[imax-i][j][k]-re_hx[imax-i][j-1][k])/dy;
				im_ezy210[i][j][k]=term1*im_ezy210[i][j][k]-term2*(im_hx[imax-i][j][k]-im_hx[imax-i][j-1][k])/dy;

				//Ezx+Ezy
				re_ez[imax-i][j][k]=re_ezx210[i][j][k]+re_ezy210[i][j][k];
				im_ez[imax-i][j][k]=im_ezx210[i][j][k]+im_ezy210[i][j][k];

			/*(012)*/
				//Ezx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx012[i][j][k]=term1*re_ezx012[i][j][k]+term2*(re_hy[i][j][kmax-1-k]-re_hy[i-1][j][kmax-1-k])/dx;
				im_ezx012[i][j][k]=term1*im_ezx012[i][j][k]+term2*(im_hy[i][j][kmax-1-k]-im_hy[i-1][j][kmax-1-k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezy012[i][j][k]=term1*re_ezy012[i][j][k]-term2*(re_hx[i][j][kmax-1-k]-re_hx[i][j-1][kmax-1-k])/dy;
				im_ezy012[i][j][k]=term1*im_ezy012[i][j][k]-term2*(im_hx[i][j][kmax-1-k]-im_hx[i][j-1][kmax-1-k])/dy;

				//Ezx+Ezy
				re_ez[i][j][kmax-1-k]=re_ezx012[i][j][k]+re_ezy012[i][j][k];
				im_ez[i][j][kmax-1-k]=im_ezx012[i][j][k]+im_ezy012[i][j][k];

			/*(212)*/
				//Ezx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx212[i][j][k]=term1*re_ezx212[i][j][k]+term2*(re_hy[imax-i][j][kmax-1-k]-re_hy[imax-i-1][j][kmax-1-k])/dx;
				im_ezx212[i][j][k]=term1*im_ezx212[i][j][k]+term2*(im_hy[imax-i][j][kmax-1-k]-im_hy[imax-i-1][j][kmax-1-k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezy212[i][j][k]=term1*re_ezy212[i][j][k]-term2*(re_hx[imax-i][j][kmax-1-k]-re_hx[imax-i][j-1][kmax-1-k])/dy;
				im_ezy212[i][j][k]=term1*im_ezy212[i][j][k]-term2*(im_hx[imax-i][j][kmax-1-k]-im_hx[imax-i][j-1][kmax-1-k])/dy;

				//Ezx+Ezy
				re_ez[imax-i][j][kmax-1-k]=re_ezx212[i][j][k]+re_ezy212[i][j][k];
				im_ez[imax-i][j][kmax-1-k]=im_ezx212[i][j][k]+im_ezy212[i][j][k];
			}
		}
	}

	for(i=1;i<L;i++){
		for(j=1;j<=jmax;j++){
			for(k=L;k<kmax-L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=0.0;

			/*(011)*/
				//Ezx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx011[i][j][k]=term1*re_ezx011[i][j][k]+term2*(re_hy[i][j][k]-re_hy[i-1][j][k])/dx;
				im_ezx011[i][j][k]=term1*im_ezx011[i][j][k]+term2*(im_hy[i][j][k]-im_hy[i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezy011[i][j][k]=term1*re_ezy011[i][j][k]-term2*(re_hx[i][j][k]-re_hx[i][j-1][k])/dy;
				im_ezy011[i][j][k]=term1*im_ezy011[i][j][k]-term2*(im_hx[i][j][k]-im_hx[i][j-1][k])/dy;

				//Ezx+Ezy
				re_ez[i][j][k]=re_ezx011[i][j][k]+re_ezy011[i][j][k];
				im_ez[i][j][k]=im_ezx011[i][j][k]+im_ezy011[i][j][k];

			/*(211)*/
				//Ezx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx211[i][j][k]=term1*re_ezx211[i][j][k]+term2*(re_hy[imax-i][j][k]-re_hy[imax-i-1][j][k])/dx;
				im_ezx211[i][j][k]=term1*im_ezx211[i][j][k]+term2*(im_hy[imax-i][j][k]-im_hy[imax-i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezy211[i][j][k]=term1*re_ezy211[i][j][k]-term2*(re_hx[imax-i][j][k]-re_hx[imax-i][j-1][k])/dy;
				im_ezy211[i][j][k]=term1*im_ezy211[i][j][k]-term2*(im_hx[imax-i][j][k]-im_hx[imax-i][j-1][k])/dy;

				//Ezx+Ezy
				re_ez[imax-i][j][k]=re_ezx211[i][j][k]+re_ezy211[i][j][k];
				im_ez[imax-i][j][k]=im_ezx211[i][j][k]+im_ezy211[i][j][k];
			}
		}
	}

	for(i=L;i<=imax-L;i++){
		for(j=1;j<=jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=0.0;
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(110)*/
				//Ezx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx110[i][j][k]=term1*re_ezx110[i][j][k]+term2*(re_hy[i][j][k]-re_hy[i-1][j][k])/dx;
				im_ezx110[i][j][k]=term1*im_ezx110[i][j][k]+term2*(im_hy[i][j][k]-im_hy[i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezy110[i][j][k]=term1*re_ezy110[i][j][k]-term2*(re_hx[i][j][k]-re_hx[i][j-1][k])/dy;
				im_ezy110[i][j][k]=term1*im_ezy110[i][j][k]-term2*(im_hx[i][j][k]-im_hx[i][j-1][k])/dy;

				//Exz+Exy
				re_ez[i][j][k]=re_ezx110[i][j][k]+re_ezy110[i][j][k];
				im_ez[i][j][k]=im_ezx110[i][j][k]+im_ezy110[i][j][k];
			
			/*(112)*/
				//Ezx
				u1=(sigma_x*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx112[i][j][k]=term1*re_ezx112[i][j][k]+term2*(re_hy[i][j][kmax-1-k]-re_hy[i-1][j][kmax-1-k])/dx;
				im_ezx112[i][j][k]=term1*im_ezx112[i][j][k]+term2*(im_hy[i][j][kmax-1-k]-im_hy[i-1][j][kmax-1-k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//変更場所(終り)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezy112[i][j][k]=term1*re_ezy112[i][j][k]-term2*(re_hx[i][j][kmax-1-k]-re_hx[i][j-1][kmax-1-k])/dy;
				im_ezy112[i][j][k]=term1*im_ezy112[i][j][k]-term2*(im_hx[i][j][kmax-1-k]-im_hx[i][j-1][kmax-1-k])/dy;

				//Ezx+Ezy
				re_ez[i][j][kmax-1-k]=re_ezx112[i][j][k]+re_ezy112[i][j][k];
				im_ez[i][j][kmax-1-k]=im_ezx112[i][j][k]+im_ezy112[i][j][k];
			}
		}
	}
}

void hfield_PML(){

/*********** Hx-PML ***********/ 
	sigma_x=0.0;sigma_y=0.0;sigma_z=0.0;
	u1=0.0;u2=0.0;term1=0.0;term2=0.0;

	for(i=1;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(010)*/
				//Hxz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz010[i][j][k]=term1*re_hxz010[i][j][k]+term2*(re_ey[i][j][k+1]-re_ey[i][j][k])/dz;
				im_hxz010[i][j][k]=term1*im_hxz010[i][j][k]+term2*(im_ey[i][j][k+1]-im_ey[i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxy010[i][j][k]=term1*re_hxy010[i][j][k]-term2*(re_ez[i][j+1][k]-re_ez[i][j][k])/dy;
				im_hxy010[i][j][k]=term1*im_hxy010[i][j][k]-term2*(im_ez[i][j+1][k]-im_ez[i][j][k])/dy;

				//Hxz+Hxy
				re_hx[i][j][k]=re_hxz010[i][j][k]+re_hxy010[i][j][k];
				im_hx[i][j][k]=im_hxz010[i][j][k]+im_hxy010[i][j][k];

			/*(210)*/
				//Hxz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz210[i][j][k]=term1*re_hxz210[i][j][k]+term2*(re_ey[imax-i][j][k+1]-re_ey[imax-i][j][k])/dz;
				im_hxz210[i][j][k]=term1*im_hxz210[i][j][k]+term2*(im_ey[imax-i][j][k+1]-im_ey[imax-i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxy210[i][j][k]=term1*re_hxy210[i][j][k]-term2*(re_ez[imax-i][j+1][k]-re_ez[imax-i][j][k])/dy;
				im_hxy210[i][j][k]=term1*im_hxy210[i][j][k]-term2*(im_ez[imax-i][j+1][k]-im_ez[imax-i][j][k])/dy;

				//Hxz+Hxy		
				re_hx[imax-i][j][k]=re_hxz210[i][j][k]+re_hxy210[i][j][k];
				im_hx[imax-i][j][k]=im_hxz210[i][j][k]+im_hxy210[i][j][k];
				
			/*(012)*/
				//Hxz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz012[i][j][k]=term1*re_hxz012[i][j][k]+term2*(re_ey[i][j][kmax-1-k+1]-re_ey[i][j][kmax-1-k])/dz;
				im_hxz012[i][j][k]=term1*im_hxz012[i][j][k]+term2*(im_ey[i][j][kmax-1-k+1]-im_ey[i][j][kmax-1-k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxy012[i][j][k]=term1*re_hxy012[i][j][k]-term2*(re_ez[i][j+1][kmax-1-k]-re_ez[i][j][kmax-1-k])/dy;
				im_hxy012[i][j][k]=term1*im_hxy012[i][j][k]-term2*(im_ez[i][j+1][kmax-1-k]-im_ez[i][j][kmax-1-k])/dy;

				//Hxz+Hxy
				re_hx[i][j][kmax-1-k]=re_hxz012[i][j][k]+re_hxy012[i][j][k];
				im_hx[i][j][kmax-1-k]=im_hxz012[i][j][k]+im_hxy012[i][j][k];

			/*(212)*/
				//Hxz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz212[i][j][k]=term1*re_hxz212[i][j][k]+term2*(re_ey[imax-i][j][kmax-1-k+1]-re_ey[imax-i][j][kmax-1-k])/dz;
				im_hxz212[i][j][k]=term1*im_hxz212[i][j][k]+term2*(im_ey[imax-i][j][kmax-1-k+1]-im_ey[imax-i][j][kmax-1-k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxy212[i][j][k]=term1*re_hxy212[i][j][k]-term2*(re_ez[imax-i][j+1][kmax-1-k]-re_ez[imax-i][j][kmax-1-k])/dy;
				im_hxy212[i][j][k]=term1*im_hxy212[i][j][k]-term2*(im_ez[imax-i][j+1][kmax-1-k]-im_ez[imax-i][j][kmax-1-k])/dy;

				//Hxz+Hxy
				re_hx[imax-i][j][kmax-1-k]=re_hxz212[i][j][k]+re_hxy212[i][j][k];
				im_hx[imax-i][j][kmax-1-k]=im_hxz212[i][j][k]+im_hxy212[i][j][k];
			}
		}
	}

	for(i=1;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=L;k<kmax-L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=0.0;

			/*(011)*/
				//Hxz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz011[i][j][k]=term1*re_hxz011[i][j][k]+term2*(re_ey[i][j][k+1]-re_ey[i][j][k])/dz;
				im_hxz011[i][j][k]=term1*im_hxz011[i][j][k]+term2*(im_ey[i][j][k+1]-im_ey[i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxy011[i][j][k]=term1*re_hxy011[i][j][k]-term2*(re_ez[i][j+1][k]-re_ez[i][j][k])/dy;
				im_hxy011[i][j][k]=term1*im_hxy011[i][j][k]-term2*(im_ez[i][j+1][k]-im_ez[i][j][k])/dy;

				//Hxz+Hxy
				re_hx[i][j][k]=re_hxz011[i][j][k]+re_hxy011[i][j][k];
				im_hx[i][j][k]=im_hxz011[i][j][k]+im_hxy011[i][j][k];
				
			/*(211)*/
				//Hxz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz211[i][j][k]=term1*re_hxz211[i][j][k]+term2*(re_ey[imax-i][j][k+1]-re_ey[imax-i][j][k])/dz;
				im_hxz211[i][j][k]=term1*im_hxz211[i][j][k]+term2*(im_ey[imax-i][j][k+1]-im_ey[imax-i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxy211[i][j][k]=term1*re_hxy211[i][j][k]-term2*(re_ez[imax-i][j+1][k]-re_ez[imax-i][j][k])/dy;
				im_hxy211[i][j][k]=term1*im_hxy211[i][j][k]-term2*(im_ez[imax-i][j+1][k]-im_ez[imax-i][j][k])/dy;

				//Hxz+Hxy
				re_hx[imax-i][j][k]=re_hxz211[i][j][k]+re_hxy211[i][j][k];
				im_hx[imax-i][j][k]=im_hxz211[i][j][k]+im_hxy211[i][j][k];
			}
		}
	}

	for(i=L;i<=imax-L;i++){
		for(j=0;j<jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=0.0;
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(110)*/
				//Hxz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz110[i][j][k]=term1*re_hxz110[i][j][k]+term2*(re_ey[i][j][k+1]-re_ey[i][j][k])/dz;
				im_hxz110[i][j][k]=term1*im_hxz110[i][j][k]+term2*(im_ey[i][j][k+1]-im_ey[i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxy110[i][j][k]=term1*re_hxy110[i][j][k]-term2*(re_ez[i][j+1][k]-re_ez[i][j][k])/dy;
				im_hxy110[i][j][k]=term1*im_hxy110[i][j][k]-term2*(im_ez[i][j+1][k]-im_ez[i][j][k])/dy;

				//Hxz+Hxy
				re_hx[i][j][k]=re_hxz110[i][j][k]+re_hxy110[i][j][k];
				im_hx[i][j][k]=im_hxz110[i][j][k]+im_hxy110[i][j][k];

			/*(112)*/
				//Hxz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz112[i][j][k]=term1*re_hxz112[i][j][k]+term2*(re_ey[i][j][kmax-1-k+1]-re_ey[i][j][kmax-1-k])/dz;
				im_hxz112[i][j][k]=term1*im_hxz112[i][j][k]+term2*(im_ey[i][j][kmax-1-k+1]-im_ey[i][j][kmax-1-k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxy112[i][j][k]=term1*re_hxy112[i][j][k]-term2*(re_ez[i][j+1][kmax-1-k]-re_ez[i][j][kmax-1-k])/dy;
				im_hxy112[i][j][k]=term1*im_hxy112[i][j][k]-term2*(im_ez[i][j+1][kmax-1-k]-im_ez[i][j][kmax-1-k])/dy;

				//Hxz+Hxy
				re_hx[i][j][kmax-1-k]=re_hxz112[i][j][k]+re_hxy112[i][j][k];
				im_hx[i][j][kmax-1-k]=im_hxz112[i][j][k]+im_hxy112[i][j][k];
			}
		}
	}

/*********** Hy-PML ***********/ 
	sigma_x=0.0;sigma_y=0.0;sigma_z=0.0;
	u1=0.0;u2=0.0;term1=0.0;term2=0.0;

	for(i=1;i<L;i++){
		for(j=0;j<=jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(010)*/
				//Hyx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx010[i][j][k]=term1*re_hyx010[i][j][k]+term2*(re_ez[i+1][j][k]-re_ez[i][j][k])/dx;
				im_hyx010[i][j][k]=term1*im_hyx010[i][j][k]+term2*(im_ez[i+1][j][k]-im_ez[i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyz010[i][j][k]=term1*re_hyz010[i][j][k]-term2*(re_ex[i][j][k+1]-re_ex[i][j][k])/dz;
				im_hyz010[i][j][k]=term1*im_hyz010[i][j][k]-term2*(im_ex[i][j][k+1]-im_ex[i][j][k])/dz;

				//Hyx+Hyz
				re_hy[i][j][k]=re_hyx010[i][j][k]+re_hyz010[i][j][k];
				im_hy[i][j][k]=im_hyx010[i][j][k]+im_hyz010[i][j][k];

			/*(210)*/
				//Hyx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx210[i][j][k]=term1*re_hyx210[i][j][k]+term2*(re_ez[imax-1-i+1][j][k]-re_ez[imax-1-i][j][k])/dx;
				im_hyx210[i][j][k]=term1*im_hyx210[i][j][k]+term2*(im_ez[imax-1-i+1][j][k]-im_ez[imax-1-i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyz210[i][j][k]=term1*re_hyz210[i][j][k]-term2*(re_ex[imax-1-i][j][k+1]-re_ex[imax-1-i][j][k])/dz;
				im_hyz210[i][j][k]=term1*im_hyz210[i][j][k]-term2*(im_ex[imax-1-i][j][k+1]-im_ex[imax-1-i][j][k])/dz;

				//Hyx+Hyz
				re_hy[imax-1-i][j][k]=re_hyx210[i][j][k]+re_hyz210[i][j][k];
				im_hy[imax-1-i][j][k]=im_hyx210[i][j][k]+im_hyz210[i][j][k];

			/*(012)*/
				//Hyx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx012[i][j][k]=term1*re_hyx012[i][j][k]+term2*(re_ez[i+1][j][kmax-1-k]-re_ez[i][j][kmax-1-k])/dx;
				im_hyx012[i][j][k]=term1*im_hyx012[i][j][k]+term2*(im_ez[i+1][j][kmax-1-k]-im_ez[i][j][kmax-1-k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyz012[i][j][k]=term1*re_hyz012[i][j][k]-term2*(re_ex[i][j][kmax-1-k+1]-re_ex[i][j][kmax-1-k])/dz;
				im_hyz012[i][j][k]=term1*im_hyz012[i][j][k]-term2*(im_ex[i][j][kmax-1-k+1]-im_ex[i][j][kmax-1-k])/dz;

				//Hyx+Hyz
				re_hy[i][j][kmax-1-k]=re_hyx012[i][j][k]+re_hyz012[i][j][k];
				im_hy[i][j][kmax-1-k]=im_hyx012[i][j][k]+im_hyz012[i][j][k];

			/*(210)*/
				//Hyx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx212[i][j][k]=term1*re_hyx212[i][j][k]+term2*(re_ez[imax-1-i+1][j][kmax-1-k]-re_ez[imax-1-i][j][kmax-1-k])/dx;
				im_hyx212[i][j][k]=term1*im_hyx212[i][j][k]+term2*(im_ez[imax-1-i+1][j][kmax-1-k]-im_ez[imax-1-i][j][kmax-1-k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyz212[i][j][k]=term1*re_hyz212[i][j][k]-term2*(re_ex[imax-1-i][j][kmax-1-k+1]-re_ex[imax-1-i][j][kmax-1-k])/dz;
				im_hyz212[i][j][k]=term1*im_hyz212[i][j][k]-term2*(im_ex[imax-1-i][j][kmax-1-k+1]-im_ex[imax-1-i][j][kmax-1-k])/dz;

				//Hyx+Hyz
				re_hy[imax-1-i][j][kmax-1-k]=re_hyx212[i][j][k]+re_hyz212[i][j][k];
				im_hy[imax-1-i][j][kmax-1-k]=im_hyx212[i][j][k]+im_hyz212[i][j][k];
			}
		}
	}

	for(i=1;i<L;i++){
		for(j=0;j<=jmax;j++){
			for(k=L;k<kmax-L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=0.0;

			/*(011)*/
				//Hyx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx011[i][j][k]=term1*re_hyx011[i][j][k]+term2*(re_ez[i+1][j][k]-re_ez[i][j][k])/dx;
				im_hyx011[i][j][k]=term1*im_hyx011[i][j][k]+term2*(im_ez[i+1][j][k]-im_ez[i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyz011[i][j][k]=term1*re_hyz011[i][j][k]-term2*(re_ex[i][j][k+1]-re_ex[i][j][k])/dz;
				im_hyz011[i][j][k]=term1*im_hyz011[i][j][k]-term2*(im_ex[i][j][k+1]-im_ex[i][j][k])/dz;

				//Hyx+Hyz
				re_hy[i][j][k]=re_hyx011[i][j][k]+re_hyz011[i][j][k];
				im_hy[i][j][k]=im_hyx011[i][j][k]+im_hyz011[i][j][k];

			/*(211)*/
				//Hyx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx211[i][j][k]=term1*re_hyx211[i][j][k]+term2*(re_ez[imax-1-i+1][j][k]-re_ez[imax-1-i][j][k])/dx;
				im_hyx211[i][j][k]=term1*im_hyx211[i][j][k]+term2*(im_ez[imax-1-i+1][j][k]-im_ez[imax-1-i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyz211[i][j][k]=term1*re_hyz211[i][j][k]-term2*(re_ex[imax-1-i][j][k+1]-re_ex[imax-1-i][j][k])/dz;
				im_hyz211[i][j][k]=term1*im_hyz211[i][j][k]-term2*(im_ex[imax-1-i][j][k+1]-im_ex[imax-1-i][j][k])/dz;

				//Hyx+Hyz
				re_hy[imax-1-i][j][k]=re_hyx211[i][j][k]+re_hyz211[i][j][k];
				im_hy[imax-1-i][j][k]=im_hyx211[i][j][k]+im_hyz211[i][j][k];
			}
		}
	}

	for(i=L;i<imax-L;i++){
		for(j=0;j<=jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=0.0;
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(110)*/
				//Hyx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx110[i][j][k]=term1*re_hyx110[i][j][k]+term2*(re_ez[i+1][j][k]-re_ez[i][j][k])/dx;
				im_hyx110[i][j][k]=term1*im_hyx110[i][j][k]+term2*(im_ez[i+1][j][k]-im_ez[i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyz110[i][j][k]=term1*re_hyz110[i][j][k]-term2*(re_ex[i][j][k+1]-re_ex[i][j][k])/dz;
				im_hyz110[i][j][k]=term1*im_hyz110[i][j][k]-term2*(im_ex[i][j][k+1]-im_ex[i][j][k])/dz;

				//Hyx+Hyz
				re_hy[i][j][k]=re_hyx110[i][j][k]+re_hyz110[i][j][k];
				im_hy[i][j][k]=im_hyx110[i][j][k]+im_hyz110[i][j][k];

			/*(112)*/
				//Hyx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx112[i][j][k]=term1*re_hyx112[i][j][k]+term2*(re_ez[i+1][j][kmax-1-k]-re_ez[i][j][kmax-1-k])/dx;
				im_hyx112[i][j][k]=term1*im_hyx112[i][j][k]+term2*(im_ez[i+1][j][kmax-1-k]-im_ez[i][j][kmax-1-k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyz112[i][j][k]=term1*re_hyz112[i][j][k]-term2*(re_ex[i][j][kmax-1-k+1]-re_ex[i][j][kmax-1-k])/dz;
				im_hyz112[i][j][k]=term1*im_hyz112[i][j][k]-term2*(im_ex[i][j][kmax-1-k+1]-im_ex[i][j][kmax-1-k])/dz;

				//Hyx+Hyz
				re_hy[i][j][kmax-1-k]=re_hyx112[i][j][k]+re_hyz112[i][j][k];
				im_hy[i][j][kmax-1-k]=im_hyx112[i][j][k]+im_hyz112[i][j][k];
			}
		}
	}

/*********** Hz-PML ***********/ 
	sigma_x=0.0;sigma_y=0.0;sigma_z=0.0;
	u1=0.0;u2=0.0;term1=0.0;term2=0.0;

	for(i=1;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(010)*/
				//Hzx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx010[i][j][k]=term1*re_hzx010[i][j][k]-term2*(re_ey[i+1][j][k]-re_ey[i][j][k])/dx;
				im_hzx010[i][j][k]=term1*im_hzx010[i][j][k]-term2*(im_ey[i+1][j][k]-im_ey[i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzy010[i][j][k]=term1*re_hzy010[i][j][k]+term2*(re_ex[i][j+1][k]-re_ex[i][j][k])/dy;
				im_hzy010[i][j][k]=term1*im_hzy010[i][j][k]+term2*(im_ex[i][j+1][k]-im_ex[i][j][k])/dy;

				//Hzx+Hzy
				re_hz[i][j][k]=re_hzx010[i][j][k]+re_hzy010[i][j][k];
				im_hz[i][j][k]=im_hzx010[i][j][k]+im_hzy010[i][j][k];
				
			/*(210)*/
				//Hzx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx210[i][j][k]=term1*re_hzx210[i][j][k]-term2*(re_ey[imax-1-i+1][j][k]-re_ey[imax-1-i][j][k])/dx;
				im_hzx210[i][j][k]=term1*im_hzx210[i][j][k]-term2*(im_ey[imax-1-i+1][j][k]-im_ey[imax-1-i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzy210[i][j][k]=term1*re_hzy210[i][j][k]+term2*(re_ex[imax-1-i][j+1][k]-re_ex[imax-1-i][j][k])/dy;
				im_hzy210[i][j][k]=term1*im_hzy210[i][j][k]+term2*(im_ex[imax-1-i][j+1][k]-im_ex[imax-1-i][j][k])/dy;

				//Hzx+Hzy
				re_hz[imax-1-i][j][k]=re_hzx210[i][j][k]+re_hzy210[i][j][k];
				im_hz[imax-1-i][j][k]=im_hzx210[i][j][k]+im_hzy210[i][j][k];

			/*(012)*/
				//Hzx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx012[i][j][k]=term1*re_hzx012[i][j][k]-term2*(re_ey[i+1][j][kmax-k]-re_ey[i][j][kmax-k])/dx;
				im_hzx012[i][j][k]=term1*im_hzx012[i][j][k]-term2*(im_ey[i+1][j][kmax-k]-im_ey[i][j][kmax-k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzy012[i][j][k]=term1*re_hzy012[i][j][k]+term2*(re_ex[i][j+1][kmax-k]-re_ex[i][j][kmax-k])/dy;
				im_hzy012[i][j][k]=term1*im_hzy012[i][j][k]+term2*(im_ex[i][j+1][kmax-k]-im_ex[i][j][kmax-k])/dy;

				//Hzx+Hzy
				re_hz[i][j][kmax-k]=re_hzx012[i][j][k]+re_hzy012[i][j][k];
				im_hz[i][j][kmax-k]=im_hzx012[i][j][k]+im_hzy012[i][j][k];

			/*(212)*/
				//Hzx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx212[i][j][k]=term1*re_hzx212[i][j][k]-term2*(re_ey[imax-1-i+1][j][kmax-k]-re_ey[imax-1-i][j][kmax-k])/dx;
				im_hzx212[i][j][k]=term1*im_hzx212[i][j][k]-term2*(im_ey[imax-1-i+1][j][kmax-k]-im_ey[imax-1-i][j][kmax-k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzy212[i][j][k]=term1*re_hzy212[i][j][k]+term2*(re_ex[imax-1-i][j+1][kmax-k]-re_ex[imax-1-i][j][kmax-k])/dy;
				im_hzy212[i][j][k]=term1*im_hzy212[i][j][k]+term2*(im_ex[imax-1-i][j+1][kmax-k]-im_ex[imax-1-i][j][kmax-k])/dy;

				//Hzx+Hzy
				re_hz[imax-1-i][j][kmax-k]=re_hzx212[i][j][k]+re_hzy212[i][j][k];
				im_hz[imax-1-i][j][kmax-k]=im_hzx212[i][j][k]+im_hzy212[i][j][k];
			}
		}
	}

	for(i=1;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=L;k<=kmax-L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=0.0;

			/*(011)*/
				//Hzx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx011[i][j][k]=term1*re_hzx011[i][j][k]-term2*(re_ey[i+1][j][k]-re_ey[i][j][k])/dx;
				im_hzx011[i][j][k]=term1*im_hzx011[i][j][k]-term2*(im_ey[i+1][j][k]-im_ey[i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzy011[i][j][k]=term1*re_hzy011[i][j][k]+term2*(re_ex[i][j+1][k]-re_ex[i][j][k])/dy;
				im_hzy011[i][j][k]=term1*im_hzy011[i][j][k]+term2*(im_ex[i][j+1][k]-im_ex[i][j][k])/dy;

				//Hzx+Hzy
				re_hz[i][j][k]=re_hzx011[i][j][k]+re_hzy011[i][j][k];
				im_hz[i][j][k]=im_hzx011[i][j][k]+im_hzy011[i][j][k];

			/*(211)*/
				//Hzx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx211[i][j][k]=term1*re_hzx211[i][j][k]-term2*(re_ey[imax-1-i+1][j][k]-re_ey[imax-1-i][j][k])/dx;
				im_hzx211[i][j][k]=term1*im_hzx211[i][j][k]-term2*(im_ey[imax-1-i+1][j][k]-im_ey[imax-1-i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzy211[i][j][k]=term1*re_hzy211[i][j][k]+term2*(re_ex[imax-1-i][j+1][k]-re_ex[imax-1-i][j][k])/dy;
				im_hzy211[i][j][k]=term1*im_hzy211[i][j][k]+term2*(im_ex[imax-1-i][j+1][k]-im_ex[imax-1-i][j][k])/dy;

				//Hzx+Hzy
				re_hz[imax-1-i][j][k]=re_hzx211[i][j][k]+re_hzy211[i][j][k];
				im_hz[imax-1-i][j][k]=im_hzx211[i][j][k]+im_hzy211[i][j][k];
			}
		}
	}

	for(i=L;i<imax-L;i++){
		for(j=0;j<jmax;j++){
			for(k=1;k<L;k++){

				sigma_x=0.0;
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(110)*/
				//Hzx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx110[i][j][k]=term1*re_hzx110[i][j][k]-term2*(re_ey[i+1][j][k]-re_ey[i][j][k])/dx;
				im_hzx110[i][j][k]=term1*im_hzx110[i][j][k]-term2*(im_ey[i+1][j][k]-im_ey[i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzy110[i][j][k]=term1*re_hzy110[i][j][k]+term2*(re_ex[i][j+1][k]-re_ex[i][j][k])/dy;
				im_hzy110[i][j][k]=term1*im_hzy110[i][j][k]+term2*(im_ex[i][j+1][k]-im_ex[i][j][k])/dy;

				//Hzx+Hzy
				re_hz[i][j][k]=re_hzx110[i][j][k]+re_hzy110[i][j][k];
				im_hz[i][j][k]=im_hzx110[i][j][k]+im_hzy110[i][j][k];

			/*(112)*/
				//Hzx
				u1=(sigma_x*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx112[i][j][k]=term1*re_hzx112[i][j][k]-term2*(re_ey[i+1][j][kmax-k]-re_ey[i][j][kmax-k])/dx;
				im_hzx112[i][j][k]=term1*im_hzx112[i][j][k]-term2*(im_ey[i+1][j][kmax-k]-im_ey[i][j][kmax-k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//変更場所(始め)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//変更場所(終り)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzy112[i][j][k]=term1*re_hzy112[i][j][k]+term2*(re_ex[i][j+1][kmax-k]-re_ex[i][j][kmax-k])/dy;
				im_hzy112[i][j][k]=term1*im_hzy112[i][j][k]+term2*(im_ex[i][j+1][kmax-k]-im_ex[i][j][kmax-k])/dy;

				//Hzx+Hzy
				re_hz[i][j][kmax-k]=re_hzx112[i][j][k]+re_hzy112[i][j][k];
				im_hz[i][j][kmax-k]=im_hzx112[i][j][k]+im_hzy112[i][j][k];
			}
		}
	}
}
/********************************** 吸収境界条件終わり **********************************/


/*************************** 解析モデル＆フィールドの吐き出し ***************************/
void nnprint(){
//	FILE *NN;
//	int i,j;
//	NN=fopen("nn.xls","w");
//	for(i=L;i<imax-L;i++){
//		for(j=0;j<jmax;j++){
//			fprintf(NN,"%d\t",nn[i][j][ck]);
//		}
//    fprintf(NN,"\n");
//	}
//	fclose(NN);

//変更場所(始め)
	FILE *MODEL_XY;
	FILE *MODEL_XZ;
	FILE *MODEL_YZ;
	int i,j,k;
	//MODEL_XY=fopen("model_xy.txt","w");
	//MODEL_XZ=fopen("model_xz.txt","w");
	//MODEL_YZ=fopen("model_yz.txt","w");
	MODEL_XY=fopen("model_xy.xls","w");
	MODEL_XZ=fopen("model_xz.xls","w");
	MODEL_YZ=fopen("model_yz.xls","w");

	for(i=L;i<imax-L;i++){
		for(j=0;j<jmax;j++){
			fprintf(MODEL_XY,"%d\t",nn[i][j][ck]);
		}
    fprintf(MODEL_XY,"\n");
	}

	for(i=L;i<imax-L;i++){
		for(k=L;k<kmax-L;k++){
			fprintf(MODEL_XZ,"%d\t",nn[i][cj][k]);
		}
    fprintf(MODEL_XZ,"\n");
	}

	for(j=0;j<jmax;j++){
		for(k=L;k<kmax-L;k++){
			fprintf(MODEL_YZ,"%d\t",nn[ci][j][k]);
		}
    fprintf(MODEL_YZ,"\n");
	}

	fclose(MODEL_XY);
	fclose(MODEL_XZ);
	fclose(MODEL_YZ);
//変更場所(終り)

}

void file_open_field(){
	if(ikk<10) sprintf(fikk,"%d_0%d_%05d",m,ikk,n); else sprintf(fikk,"%d_%d_%05d",m,ikk,n);

//	sprintf(fname,"%sEx_field.xls",fikk);EX=fopen(fname,"w");
//	sprintf(fname,"%sEy_field.xls",fikk);EY=fopen(fname,"w");
//	sprintf(fname,"%sEz_field.xls",fikk);EZ=fopen(fname,"w");
//	sprintf(fname,"%sHx_field.xls",fikk);HX=fopen(fname,"w");
//	sprintf(fname,"%sHy_field.xls",fikk);HY=fopen(fname,"w");
//	sprintf(fname,"%sHz_field.xls",fikk);HZ=fopen(fname,"w");

//変更場所(始め)
	if(excite == 1){
		sprintf(fname,"%sEx_XY_field.txt",fikk);EX_XY=fopen(fname,"w");
		sprintf(fname,"%sEx_XZ_center_field.txt",fikk);EX_XZ_center=fopen(fname,"w");
		sprintf(fname,"%sEx_XZ_end_field.txt",fikk);EX_XZ_end=fopen(fname,"w");
		sprintf(fname,"%sEx_YZ_field.txt",fikk);EX_YZ=fopen(fname,"w");

		sprintf(fname,"%sEy_XY_field.txt",fikk);EY_XY=fopen(fname,"w");
		sprintf(fname,"%sEy_XZ_center_field.txt",fikk);EY_XZ_center=fopen(fname,"w");
		sprintf(fname,"%sEy_XZ_end_field.txt",fikk);EY_XZ_end=fopen(fname,"w");
		sprintf(fname,"%sEy_YZ_field.txt",fikk);EY_YZ=fopen(fname,"w");

		sprintf(fname,"%sHz_XY_field.txt",fikk);HZ_XY=fopen(fname,"w");
		sprintf(fname,"%sHz_XZ_center_field.txt",fikk);HZ_XZ_center=fopen(fname,"w");
		sprintf(fname,"%sHz_XZ_end_field.txt",fikk);HZ_XZ_end=fopen(fname,"w");
		sprintf(fname,"%sHz_YZ_field.txt",fikk);HZ_YZ=fopen(fname,"w");
	}
	if(excite == 2){
		sprintf(fname,"%sEz_XY_field.txt",fikk);EZ_XY=fopen(fname,"w");
		sprintf(fname,"%sEz_XZ_center_field.txt",fikk);EZ_XZ_center=fopen(fname,"w");
		sprintf(fname,"%sEz_XZ_end_field.txt",fikk);EZ_XZ_end=fopen(fname,"w");
		sprintf(fname,"%sEz_YZ_field.txt",fikk);EZ_YZ=fopen(fname,"w");

		sprintf(fname,"%sHx_XY_field.txt",fikk);HX_XY=fopen(fname,"w");
		sprintf(fname,"%sHx_XZ_center_field.txt",fikk);HX_XZ_center=fopen(fname,"w");
		sprintf(fname,"%sHx_XZ_end_field.txt",fikk);HX_XZ_end=fopen(fname,"w");
		sprintf(fname,"%sHx_YZ_field.txt",fikk);HX_YZ=fopen(fname,"w");

		sprintf(fname,"%sHy_XY_field.txt",fikk);HY_XY=fopen(fname,"w");
		sprintf(fname,"%sHy_XZ_center_field.txt",fikk);HY_XZ_center=fopen(fname,"w");
		sprintf(fname,"%sHy_XZ_end_field.txt",fikk);HY_XZ_end=fopen(fname,"w");
		sprintf(fname,"%sHy_YZ_field.txt",fikk);HY_YZ=fopen(fname,"w");
	}
//変更場所(終り)

}

void field_print(int n){
//	int i,j;

//変更場所(始め)
	int i,j,k;
//変更場所(終り)

	printf("Output Field t = %d\n",n);

	for(i=L;i<imax-L;i++){
		for(j=0;j<jmax;j++){
//			fprintf(EX,"%e\t",re_ex[i][j][ck]);
//			fprintf(EY,"%e\t",re_ey[i][j][ck]);
//			fprintf(EZ,"%e\t",re_ez[i][j][ck]);
//			fprintf(HX,"%e\t",re_hx[i][j][ck]);
//			fprintf(HY,"%e\t",re_hy[i][j][ck]);
//			fprintf(HZ,"%e\t",re_hz[i][j][ck]);

//変更場所(始め)
			if(excite == 1){
				fprintf(EX_XY,"%e\t",re_ex[i][j][ck]);
				fprintf(EY_XY,"%e\t",re_ey[i][j][ck]);
				fprintf(HZ_XY,"%e\t",re_hz[i][j][ck]);
			}
			if(excite == 2){
				fprintf(EZ_XY,"%e\t",re_ez[i][j][ck]);
				fprintf(HX_XY,"%e\t",re_hx[i][j][ck]);
				fprintf(HY_XY,"%e\t",re_hy[i][j][ck]);
			}
//変更場所(終り)

		}
//		fprintf(EX,"\n");
//		fprintf(EY,"\n");
//		fprintf(EZ,"\n");
//		fprintf(HX,"\n");
//		fprintf(HY,"\n");
//		fprintf(HZ,"\n");

//変更場所(始め)
		if(excite == 1){
			fprintf(EX_XY,"\n");
			fprintf(EY_XY,"\n");
			fprintf(HZ_XY,"\n");
		}
		if(excite == 2){
			fprintf(EZ_XY,"\n");
			fprintf(HX_XY,"\n");
			fprintf(HY_XY,"\n");
		}
//変更場所(終り)

	}

//変更場所(始め)
	for(i=L;i<imax-L;i++){
		for(k=L;k<kmax-L;k++){
			if(excite == 1){
				fprintf(EX_XZ_center,"%e\t",re_ex[i][cj][k]);
				fprintf(EY_XZ_center,"%e\t",re_ey[i][cj][k]);
				fprintf(HZ_XZ_center,"%e\t",re_hz[i][cj][k]);

				fprintf(EX_XZ_end,"%e\t",re_ex[i][0][k]);
				fprintf(EY_XZ_end,"%e\t",re_ey[i][0][k]);
				fprintf(HZ_XZ_end,"%e\t",re_hz[i][0][k]);
			}
			if(excite == 2){
				fprintf(EZ_XZ_center,"%e\t",re_ez[i][cj][k]);
				fprintf(HX_XZ_center,"%e\t",re_hx[i][cj][k]);
				fprintf(HY_XZ_center,"%e\t",re_hy[i][cj][k]);

				fprintf(EZ_XZ_end,"%e\t",re_ez[i][0][k]);
				fprintf(HX_XZ_end,"%e\t",re_hx[i][0][k]);
				fprintf(HY_XZ_end,"%e\t",re_hy[i][0][k]);
			}
		}
		if(excite == 1){
			fprintf(EX_XZ_center,"\n");
			fprintf(EY_XZ_center,"\n");
			fprintf(HZ_XZ_center,"\n");

			fprintf(EX_XZ_end,"\n");
			fprintf(EY_XZ_end,"\n");
			fprintf(HZ_XZ_end,"\n");
		}
		if(excite == 2){
			fprintf(EZ_XZ_center,"\n");
			fprintf(HX_XZ_center,"\n");
			fprintf(HY_XZ_center,"\n");

			fprintf(EZ_XZ_end,"\n");
			fprintf(HX_XZ_end,"\n");
			fprintf(HY_XZ_end,"\n");
		}
	}

	for(j=0;j<jmax;j++){
		for(k=L;k<kmax-L;k++){
			if(excite == 1){
				fprintf(EX_YZ,"%e\t",re_ex[ci][j][k]);
				fprintf(EY_YZ,"%e\t",re_ey[ci][j][k]);
				fprintf(HZ_YZ,"%e\t",re_hz[ci][j][k]);
			}
			if(excite == 2){
				fprintf(EZ_YZ,"%e\t",re_ez[ci][j][k]);
				fprintf(HX_YZ,"%e\t",re_hx[ci][j][k]);
				fprintf(HY_YZ,"%e\t",re_hy[ci][j][k]);
			}
		}
		if(excite == 1){
			fprintf(EX_YZ,"\n");
			fprintf(EY_YZ,"\n");
			fprintf(HZ_YZ,"\n");
		}
		if(excite == 2){
			fprintf(EZ_YZ,"\n");
			fprintf(HX_YZ,"\n");
			fprintf(HY_YZ,"\n");
		}
	}
//変更場所(終り)

}


void file_close_field(){
//	fclose(EX);
//	fclose(EY);
//	fclose(EZ);
//	fclose(HX);
//	fclose(HY);
//	fclose(HZ);

//変更場所(始め)
	if(excite == 1){
		fclose(EX_XY);	fclose(EX_XZ_center);	fclose(EX_XZ_end);	fclose(EX_YZ);
		fclose(EY_XY);	fclose(EY_XZ_center);	fclose(EY_XZ_end);	fclose(EY_YZ);
		fclose(HZ_XY);	fclose(HZ_XZ_center);	fclose(HZ_XZ_end);	fclose(HZ_YZ);
	}
	if(excite == 2){
		fclose(EZ_XY);	fclose(EZ_XZ_center);	fclose(EZ_XZ_end);	fclose(EZ_YZ);
		fclose(HX_XY);	fclose(HX_XZ_center);	fclose(HX_XZ_end);	fclose(HX_YZ);
		fclose(HY_XY);	fclose(HY_XZ_center);	fclose(HY_XZ_end);	fclose(HY_YZ);
	}
//変更場所(終り)

}


void file_open_field_f(){
	if(ikk<10) sprintf(fikk,"%d_0%d_%05d",m,ikk,n); else sprintf(fikk,"%d_%d_%05d",m,ikk,n);

//	sprintf(fname,"%sEx_field.xls",fikk);EX=fopen(fname,"w");
//	sprintf(fname,"%sEy_field.xls",fikk);EY=fopen(fname,"w");
//	sprintf(fname,"%sEz_field.xls",fikk);EZ=fopen(fname,"w");
//	sprintf(fname,"%sHx_field.xls",fikk);HX=fopen(fname,"w");
//	sprintf(fname,"%sHy_field.xls",fikk);HY=fopen(fname,"w");
//	sprintf(fname,"%sHz_field.xls",fikk);HZ=fopen(fname,"w");

//変更場所(始め)
	if(excite == 1){
		sprintf(fname,"%sHz_XZ_center_field.txt",fikk);HZ_XZ_center=fopen(fname,"w");
	}
	if(excite == 2){
		sprintf(fname,"%sEz_XZ_center_field.txt",fikk);EZ_XZ_center=fopen(fname,"w");
	}
//変更場所(終り)

}


void field_print_f(int n){

//変更場所(始め)
	int i,j,k;
//変更場所(終り)

	printf("Output Field t = %d\n",n);
    i = ci;
	k = ck;

//変更場所(始め)

			if(excite == 1){
				fprintf(HZ_XZ_center,"%e\t",re_hz[i][cj][k]);
			}
			if(excite == 2){
				fprintf(EZ_XZ_center,"%e\t",re_ez[i][cj][k]);
			}
		if(excite == 1){
			fprintf(HZ_XZ_center,"\n");
		}
		if(excite == 2){
			fprintf(EZ_XZ_center,"\n");
		}
	
}

void file_close_field_f(){

//変更場所(始め)
	if(excite == 1){
			fclose(HZ_XZ_center);
	}
	if(excite == 2){
			fclose(EZ_XZ_center);	
	}
//変更場所(終り)

}



/*
解析モデル，フィールドともz座標はckのみの表示とした．
x-z or y-z方向で出力したい場合は手動で変えること．
*/
/************************ 解析モデル＆フィールドの吐き出し終わり ************************/


/******************************** 電磁界のパルス吐き出し ********************************/
void file_open(){
	if(ikk<10) sprintf(fikk,"%d_0%d",m,ikk); else sprintf(fikk,"%d_%d",m,ikk);

//変更場所v5.03(始め)
	if(excite == 1){
		sprintf(fname,"%sEx.txt",fikk);fPoynt1Ex=fopen(fname,"w");
		sprintf(fname,"%sEy.txt",fikk);fPoynt1Ey=fopen(fname,"w");
		sprintf(fname,"%sHz.txt",fikk);fPoynt1Hz=fopen(fname,"w");
	}
	if(excite == 2){
		sprintf(fname,"%sEz.txt",fikk);fPoynt1Ez=fopen(fname,"w");
		sprintf(fname,"%sHx.txt",fikk);fPoynt1Hx=fopen(fname,"w");
		sprintf(fname,"%sHy.txt",fikk);fPoynt1Hy=fopen(fname,"w");
	}
}
//変更場所(終り)

void Poyntingprint(int flag){
	int X1=ci;		int X2=ci+imove;	int X3=ci-imove;
	int Y1=cj;		int Y2=cj+jmove;	int Y3=cj-jmove;

//変更場所(始め)
	int Z1=ck;		int Z2=ck+kmove;	int Z3=ck-kmove;
//変更場所(終り)

	if(flag==0){
		printf("測定点");
		printf("\n");
//		printf("x=%d y=%d\n",X1,Y1);		// Pin(ci,cj)解析空間の中心
//		printf("x=%d y=%d\n",X2,Y2);		// Pout1(ci+imove,cj+jmove)
//		printf("x=%d y=%d\n",X2,Y3);		// Pout2(ci+imove,cj-jmove)
//		printf("x=%d y=%d\n",X3,Y2);		// Pout3(ci-imove,cj+jmove)
//		printf("x=%d y=%d\n",X3,Y3);		// Pout4(ci-imove,cj-jmove)

//変更場所(始め)
		printf("x=%d	y=%d	z=%d\n",X1,Y1,Z1);		// Pin(ci,cj)解析空間の中心
		printf("x=%d	y=%d	z=%d\n",X2,Y2,Z1);		// Pout(ci+imove,cj+jmove,ck)
		printf("x=%d	y=%d	z=%d\n",X1,Y3,Z1);		// Pout(ci,cj-jmove,ck)
		printf("x=%d	y=%d	z=%d\n",X1,Y2,Z3);		// Pout(ci,cj+jmove,ck-kmove)
		printf("x=%d	y=%d	z=%d\n",X3,Y3,Z2);		// Pout(ci-imove,cj-jmove,ck+kmove)
//変更場所(終り)

	}

	if(flag==1){
//		fprintf(fPoynt1Ex,"%e\t%e\t%e\t%e\t%e\n",re_ex[X1][Y1][ck],re_ex[X2][Y2][ck],re_ex[X2][Y3][ck],re_ex[X3][Y2][ck],re_ex[X3][Y3][ck]);
//		fprintf(fPoynt1Ey,"%e\t%e\t%e\t%e\t%e\n",re_ey[X1][Y1][ck],re_ey[X2][Y2][ck],re_ey[X2][Y3][ck],re_ey[X3][Y2][ck],re_ey[X3][Y3][ck]);
//		fprintf(fPoynt1Ez,"%e\t%e\t%e\t%e\t%e\n",re_ez[X1][Y1][ck],re_ez[X2][Y2][ck],re_ez[X2][Y3][ck],re_ez[X3][Y2][ck],re_ez[X3][Y3][ck]);
//		fprintf(fPoynt1Hx,"%e\t%e\t%e\t%e\t%e\n",re_hx[X1][Y1][ck],re_hx[X2][Y2][ck],re_hx[X2][Y3][ck],re_hx[X3][Y2][ck],re_hx[X3][Y3][ck]);
//		fprintf(fPoynt1Hy,"%e\t%e\t%e\t%e\t%e\n",re_hy[X1][Y1][ck],re_hy[X2][Y2][ck],re_hy[X2][Y3][ck],re_hy[X3][Y2][ck],re_hy[X3][Y3][ck]);
//		fprintf(fPoynt1Hz,"%e\t%e\t%e\t%e\t%e\n",re_hz[X1][Y1][ck],re_hz[X2][Y2][ck],re_hz[X2][Y3][ck],re_hz[X3][Y2][ck],re_hz[X3][Y3][ck]);

//変更場所(始め)
		if(excite == 1){
			fprintf(fPoynt1Ex,"%e,%e,%e,%e,%e\n",re_ex[X1][Y1][Z1],re_ex[X2][Y2][Z1],re_ex[X1][Y3][Z1],re_ex[X1][Y2][Z3],re_ex[X3][Y3][Z2]);
			fprintf(fPoynt1Ey,"%e,%e,%e,%e,%e\n",re_ey[X1][Y1][Z1],re_ey[X2][Y2][Z1],re_ey[X1][Y3][Z1],re_ey[X1][Y2][Z3],re_ey[X3][Y3][Z2]);
			fprintf(fPoynt1Hz,"%e,%e,%e,%e,%e\n",re_hz[X1][Y1][Z1],re_hz[X2][Y2][Z1],re_hz[X1][Y3][Z1],re_hz[X1][Y2][Z3],re_hz[X3][Y3][Z2]);
		}
		if(excite == 2){
			fprintf(fPoynt1Ez,"%e,%e,%e,%e,%e\n",re_ez[X1][Y1][Z1],re_ez[X2][Y2][Z1],re_ez[X1][Y3][Z1],re_ez[X1][Y2][Z3],re_ez[X3][Y3][Z2]);
			fprintf(fPoynt1Hx,"%e,%e,%e,%e,%e\n",re_hx[X1][Y1][Z1],re_hx[X2][Y2][Z1],re_hx[X1][Y3][Z1],re_hx[X1][Y2][Z3],re_hx[X3][Y3][Z2]);
			fprintf(fPoynt1Hy,"%e,%e,%e,%e,%e\n",re_hy[X1][Y1][Z1],re_hy[X2][Y2][Z1],re_hy[X1][Y3][Z1],re_hy[X1][Y2][Z3],re_hy[X3][Y3][Z2]);
		}
//変更場所(終り)
	
	}
}

void file_close(){
	if(excite == 1){
		fclose(fPoynt1Ex);
		fclose(fPoynt1Ey);
		fclose(fPoynt1Hz);
	}
	if(excite == 2){
		fclose(fPoynt1Ez);
		fclose(fPoynt1Hx);
		fclose(fPoynt1Hy);
	}
}
/***************************** 電磁界のパルス吐き出し終わり *****************************/


/************************************* 配列の初期化 *************************************/
void initialize_matrix(){
	for(i=0;i<imax;i++){
		for(j=0;j<=jmax;j++){
			for(k=0;k<=kmax;k++){
				re_ex[i][j][k]=0.0;			im_ex[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<=imax;i++){
		for(j=0;j<jmax;j++){
			for(k=0;k<=kmax;k++){
				re_ey[i][j][k]=0.0;			im_ey[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<=imax;i++){
		for(j=0;j<=jmax;j++){
			for(k=0;k<kmax;k++){
				re_ez[i][j][k]=0.0;			im_ez[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<=imax;i++){
		for(j=0;j<jmax;j++){
			for(k=0;k<kmax;k++){
				re_hx[i][j][k]=0.0;			im_hx[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<imax;i++){
		for(j=0;j<=jmax;j++){
			for(k=0;k<kmax;k++){
				re_hy[i][j][k]=0.0;			im_hy[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<imax;i++){
		for(j=0;j<jmax;j++){
			for(k=0;k<=kmax;k++){
				re_hz[i][j][k]=0.0;			im_hz[i][j][k]=0.0;
			}
		}
	}

	for(i=0;i<=imax;i++){
		for(k=0;k<kmax;k++){
			re_hx_jmax[i][k]=0.0;			im_hx_jmax[i][k]=0.0;
		}
	}
	for(i=0;i<imax;i++){
		for(k=0;k<=kmax;k++){
			re_hz_jmax[i][k]=0.0;			im_hz_jmax[i][k]=0.0;
		}
	}

	for(i=0;i<L;i++){
		for(j=0;j<=jmax;j++){
			for(k=0;k<L;k++){
				re_exz010[i][j][k]=0.0;	re_exy010[i][j][k]=0.0;	re_exz210[i][j][k]=0.0;	re_exy210[i][j][k]=0.0;	re_exz012[i][j][k]=0.0;	re_exy012[i][j][k]=0.0;	re_exz212[i][j][k]=0.0;	re_exy212[i][j][k]=0.0;
				im_exz010[i][j][k]=0.0;	im_exy010[i][j][k]=0.0;	im_exz210[i][j][k]=0.0;	im_exy210[i][j][k]=0.0;	im_exz012[i][j][k]=0.0;	im_exy012[i][j][k]=0.0;	im_exz212[i][j][k]=0.0;	im_exy212[i][j][k]=0.0;
				re_ezx010[i][j][k]=0.0;	re_ezy010[i][j][k]=0.0;	re_ezx210[i][j][k]=0.0;	re_ezy210[i][j][k]=0.0;	re_ezx012[i][j][k]=0.0;	re_ezy012[i][j][k]=0.0;	re_ezx212[i][j][k]=0.0;	re_ezy212[i][j][k]=0.0;
				im_ezx010[i][j][k]=0.0;	im_ezy010[i][j][k]=0.0;	im_ezx210[i][j][k]=0.0;	im_ezy210[i][j][k]=0.0;	im_ezx012[i][j][k]=0.0;	im_ezy012[i][j][k]=0.0;	im_ezx212[i][j][k]=0.0;	im_ezy212[i][j][k]=0.0;
				re_hyx010[i][j][k]=0.0;	re_hyz010[i][j][k]=0.0;	re_hyx210[i][j][k]=0.0;	re_hyz210[i][j][k]=0.0;	re_hyx012[i][j][k]=0.0;	re_hyz012[i][j][k]=0.0;	re_hyx212[i][j][k]=0.0;	re_hyz212[i][j][k]=0.0;
				im_hyx010[i][j][k]=0.0;	im_hyz010[i][j][k]=0.0;	im_hyx210[i][j][k]=0.0;	im_hyz210[i][j][k]=0.0;	im_hyx012[i][j][k]=0.0;	im_hyz012[i][j][k]=0.0;	im_hyx212[i][j][k]=0.0;	im_hyz212[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=0;k<L;k++){
				re_eyx010[i][j][k]=0.0;	re_eyz010[i][j][k]=0.0;	re_eyx210[i][j][k]=0.0;	re_eyz210[i][j][k]=0.0;	re_eyx012[i][j][k]=0.0;	re_eyz012[i][j][k]=0.0;	re_eyx212[i][j][k]=0.0;	re_eyz212[i][j][k]=0.0;
				im_eyx010[i][j][k]=0.0;	im_eyz010[i][j][k]=0.0;	im_eyx210[i][j][k]=0.0;	im_eyz210[i][j][k]=0.0;	im_eyx012[i][j][k]=0.0;	im_eyz012[i][j][k]=0.0;	im_eyx212[i][j][k]=0.0;	im_eyz212[i][j][k]=0.0;
				re_hxz010[i][j][k]=0.0;	re_hxy010[i][j][k]=0.0;	re_hxz210[i][j][k]=0.0;	re_hxy210[i][j][k]=0.0;	re_hxz012[i][j][k]=0.0;	re_hxy012[i][j][k]=0.0;	re_hxz212[i][j][k]=0.0;	re_hxy212[i][j][k]=0.0;
				im_hxz010[i][j][k]=0.0;	im_hxy010[i][j][k]=0.0;	im_hxz210[i][j][k]=0.0;	im_hxy210[i][j][k]=0.0;	im_hxz012[i][j][k]=0.0;	im_hxy012[i][j][k]=0.0;	im_hxz212[i][j][k]=0.0;	im_hxy212[i][j][k]=0.0;
				re_hzx010[i][j][k]=0.0;	re_hzy010[i][j][k]=0.0;	re_hzx210[i][j][k]=0.0;	re_hzy210[i][j][k]=0.0;	re_hzx012[i][j][k]=0.0;	re_hzy012[i][j][k]=0.0;	re_hzx212[i][j][k]=0.0;	re_hzy212[i][j][k]=0.0;
				im_hzx010[i][j][k]=0.0;	im_hzy010[i][j][k]=0.0;	im_hzx210[i][j][k]=0.0;	im_hzy210[i][j][k]=0.0;	im_hzx012[i][j][k]=0.0;	im_hzy012[i][j][k]=0.0;	im_hzx212[i][j][k]=0.0;	im_hzy212[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<imax-L;i++){
		for(j=0;j<=jmax;j++){
			for(k=0;k<L;k++){
				re_exz110[i][j][k]=0.0;	re_exy110[i][j][k]=0.0;	re_exz112[i][j][k]=0.0;	re_exy112[i][j][k]=0.0;
				im_exz110[i][j][k]=0.0;	im_exy110[i][j][k]=0.0;	im_exz112[i][j][k]=0.0;	im_exy112[i][j][k]=0.0;
				re_hyx110[i][j][k]=0.0;	re_hyz110[i][j][k]=0.0;	re_hyx112[i][j][k]=0.0;	re_hyz112[i][j][k]=0.0;
				im_hyx110[i][j][k]=0.0;	im_hyz110[i][j][k]=0.0;	im_hyx112[i][j][k]=0.0;	im_hyz112[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=0;k<=kmax-L;k++){
				re_eyx011[i][j][k]=0.0;	re_eyz011[i][j][k]=0.0;	re_eyx211[i][j][k]=0.0;	re_eyz211[i][j][k]=0.0;
				im_eyx011[i][j][k]=0.0;	im_eyz011[i][j][k]=0.0;	im_eyx211[i][j][k]=0.0;	im_eyz211[i][j][k]=0.0;
				re_hzx011[i][j][k]=0.0;	re_hzy011[i][j][k]=0.0;	re_hzx211[i][j][k]=0.0;	re_hzy211[i][j][k]=0.0;
				im_hzx011[i][j][k]=0.0;	im_hzy011[i][j][k]=0.0;	im_hzx211[i][j][k]=0.0;	im_hzy211[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<=imax-L;i++){
		for(j=0;j<jmax;j++){
			for(k=0;k<L;k++){
				re_eyx110[i][j][k]=0.0;	re_eyz110[i][j][k]=0.0;	re_eyx112[i][j][k]=0.0;	re_eyz112[i][j][k]=0.0;
				im_eyx110[i][j][k]=0.0;	im_eyz110[i][j][k]=0.0;	im_eyx112[i][j][k]=0.0;	im_eyz112[i][j][k]=0.0;
				re_hxz110[i][j][k]=0.0;	re_hxy110[i][j][k]=0.0;	re_hxz112[i][j][k]=0.0;	re_hxy112[i][j][k]=0.0;
				im_hxz110[i][j][k]=0.0;	im_hxy110[i][j][k]=0.0;	im_hxz112[i][j][k]=0.0;	im_hxy112[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<L;i++){
		for(j=0;j<=jmax;j++){
			for(k=0;k<kmax-L;k++){
				re_ezx011[i][j][k]=0.0;	re_ezy011[i][j][k]=0.0;	re_ezx211[i][j][k]=0.0;	re_ezy211[i][j][k]=0.0;
				im_ezx011[i][j][k]=0.0;	im_ezy011[i][j][k]=0.0;	im_ezx211[i][j][k]=0.0;	im_ezy211[i][j][k]=0.0;
				re_hyx011[i][j][k]=0.0;	re_hyz011[i][j][k]=0.0;	re_hyx211[i][j][k]=0.0;	re_hyz211[i][j][k]=0.0;
				im_hyx011[i][j][k]=0.0;	im_hyz011[i][j][k]=0.0;	im_hyx211[i][j][k]=0.0;	im_hyz211[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<L;i++){
		for(j=0;j<=jmax;j++){
			for(k=0;k<=kmax-L;k++){
				re_exz011[i][j][k]=0.0;	re_exy011[i][j][k]=0.0;	re_exz211[i][j][k]=0.0;	re_exy211[i][j][k]=0.0;
				im_exz011[i][j][k]=0.0;	im_exy011[i][j][k]=0.0;	im_exz211[i][j][k]=0.0;	im_exy211[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<=imax-L;i++){
		for(j=0;j<=jmax;j++){
			for(k=0;k<L;k++){
				re_ezx110[i][j][k]=0.0;	re_ezy110[i][j][k]=0.0;	re_ezx112[i][j][k]=0.0;	re_ezy112[i][j][k]=0.0;
				im_ezx110[i][j][k]=0.0;	im_ezy110[i][j][k]=0.0;	im_ezx112[i][j][k]=0.0;	im_ezy112[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<L;i++){
		for(j=0;j<jmax;j++){
			for(k=0;k<kmax-L;k++){
				re_hxz011[i][j][k]=0.0;	re_hxy011[i][j][k]=0.0;	re_hxz211[i][j][k]=0.0;	re_hxy211[i][j][k]=0.0;
				im_hxz011[i][j][k]=0.0;	im_hxy011[i][j][k]=0.0;	im_hxz211[i][j][k]=0.0;	im_hxy211[i][j][k]=0.0;
			}
		}
	}
	for(i=0;i<imax-L;i++){
		for(j=0;j<jmax;j++){
			for(k=0;k<L;k++){
				re_hzx110[i][j][k]=0.0;	re_hzy110[i][j][k]=0.0;	re_hzx112[i][j][k]=0.0;	re_hzy112[i][j][k]=0.0;
				im_hzx110[i][j][k]=0.0;	im_hzy110[i][j][k]=0.0;	im_hzx112[i][j][k]=0.0;	im_hzy112[i][j][k]=0.0;
			}
		}
	}
}
/********************************** 配列の初期化終わり **********************************/


/******************************** フーリエ変換プログラム ********************************/
void fft_main(){
	double *xr;
	static double a[F_NUMBER][STEP];
	int mf,nf,nf2;
	int f_ele,f_num,step;
	int fact,acc,count;
	double ws;
	double *data[ELE+1];
	float val;
	char InputFile[20];
	char OutputFile[20];
//	char *file[6]={"Ex","Ey","Ez","Hx","Hy","Hz"};

//変更場所(始め)
	char *file[3];
	if(excite == 1){
		file[0]="Ex";
		file[1]="Ey";
		file[2]="Hz";
	}
	if(excite == 2){
		file[0]="Ez";
		file[1]="Hx";
		file[2]="Hy";
	}
//変更場所(終り)

	ws=2*PI/dt;

	printf("FFT start\n");
	for(m=K_STARTPOINT;m<=K_ENDPOINT;m++){
		for(ikk=IKK_STARTPOINT;ikk<=IKK_ENDPOINT;ikk++){
			printf("Direction = %d, ikk = %d\n",m,ikk);
			for(f_ele=0;f_ele<ELE;f_ele++){
				printf("現在FFT中 : %s\n",file[f_ele]);

				if(ikk<10){
					sprintf(fikk,"0%d",ikk);
				}
				else{
					sprintf(fikk,"%d",ikk);
				}
//				sprintf(InputFile,"%d_%s%s.xls",m,fikk,file[f_ele]);
				sprintf(InputFile,"%d_%s%s.txt",m,fikk,file[f_ele]);	//v5.11
				if((fp1=fopen(InputFile,"r"))==NULL){
					printf("ファイルが見つかりません．---- %s\n",InputFile);
					exit(1);
				}
				for(step=0,count=0;step<STEP;step++){
					for(f_num=0;f_num<F_NUMBER;f_num++){
						if(fscanf(fp1,"%e,",&val)!=EOF){	//変更場所v5.12
							a[f_num][step]=(double)val;
							count++;
						}
					}
				}
				count/=F_NUMBER;
				for(mf=0,fact=2;fact<=count;mf++){
					fact*=2;
				}
				nf=fact/2;
				nf2=nf/2;
				acc=count-nf;
				if((xr=(double *)calloc(1,sizeof(double)*nf))==NULL){
					printf("xrが確保できません メモリ不足\n");
					exit(1);
				}
				if((data[f_ele+1]=(double *)calloc(1,sizeof(double)*nf2))==NULL){
					printf("data %dが確保できません メモリ不足\n",(f_ele+1));
					exit(1);
				}
				for(i=0;i<F_NUMBER;i++){
					for(j=0;j<nf;j++){
						xr[j]=a[i][j+acc];
					}
					rdft(nf,1,xr);
					for(j=0;j<nf2;j++){
						*(xr+j)=sqrt((*(xr+j*2))*(*(xr+j*2))+(*(xr+j*2+1))*(*(xr+j*2+1)));
					}
					for(j=0;j<nf2;j++){
						*(data[f_ele+1]+j)+=*(xr+j);
					}
				}
				free(xr);
				fclose(fp1);
			}
			if((data[0]=(double *)calloc(1,sizeof(double)*nf2))==NULL){
				printf("data 0が確保できません メモリ不足\n");
				exit(1);
			}
			for(j = 0;j<nf2;j++){
				*(data[0]+j)=PITCH/(2*PI*C0/(j*ws/nf));
			}
			sprintf(OutputFile,"0%d_FFT%s.xls",m,fikk);
			fp2=fopen(OutputFile,"w");
			for(i=0;i<nf2;i++){
				for(j=0;j<=ELE;j++){
					fprintf(fp2,"%e\t",*(data[j]+i));
				}
				fprintf(fp2,"\n");
			}
			for(j=0;j<=ELE;j++){
				free(data[j]);
			}
			fclose(fp2);
		}
	}
	printf("FFT end\n");
}
/***************************** フーリエ変換プログラム終わり *****************************/


/********************************* パラメータの吐き出し *********************************/
void Output_Condition(void){
	FILE *OC;
	OC = fopen("Condition.txt","w");
//	fprintf(OC,"AIR = %d\nSI = %d\nSIO2 = %d\nCIRCLE_MEDIUM = %d\nimax = %d\njmax = %d\nkmax = %d\nci = %d\ncj = %d\nck = %d\nPITCH = %.4e\nRADIUS = %.4e\nCORE = %.4e\nCLAD = %.4e\nNORMALFREQ = %.2f\nCELLLAMBDA = %.4e\nsigma = %.4f\nPULSEPEAK = %d\nREF_INDEX1 = %.4f\nnmax = %d\ndx = %e\ndy = %e\ndz = %e\ndt = %e"
//		,AIR,SI,SIO2,CIRCLE_MEDIUM,imax,jmax,kmax,ci,cj,ck,PITCH,RADIUS,CORE,CLAD,NORMALFREQ,CELLLAMBDA,(sigma/omega0),PULSEPEAK,REF_INDEX1,nmax,dx,dy,dz,dt);

//変更場所(始め)
	fprintf(OC,"AIR = %d\nInP = %d\nAct = %d\nSIO2 = %d\nCIRCLE_MEDIUM = %d\nimax = %d\njmax = %d\nkmax = %d\nci = %d\ncj = %d\nck = %d\nInP_INDEX = %.3f\nAct_INDEX = %.3f\nClad_INDEX = %.3f\nPITCH = %.4e\nB = %d\nRADIUS = %.4e\nRADIUS2 = %.4e\nAir-hole Cells = %d\nEffective 2r/a = %.4f\nAir_clad = %.4e\nInP_over = %.4e\nInP_down = %.4e\nInP_sub = %.4e\nAct = %.4e\nWidth = %d\n3rd_Shift = %.3e\nLAMBDA = %.3e\nsigma = %.4f\nPULSEPEAK = %d\nnmax = %d\ndx = %e\ndy = %e\ndz = %e\ndt = %e\nk vector = %d/%d\n"

		,AIR,InP_MEDIUM,Act_MEDIUM,SIO2,CIRCLE_MEDIUM
		,imax,jmax,kmax,ci,cj,ck
		,InP_INDEX,Act_INDEX,REF_INDEX2
		,PITCH,B,RADIUS,RADIUS2,effective_circle_cell,effective_2r_a,AIR_clad,InP_circle_over,InP_circle_down,InP_sub,Act,Width,THIRD_Sx
		,LAMBDA,(sigma/omega0),PULSEPEAK,nmax
		,dx,dy,dz,dt
		,IKK_ENDPOINT,IKK_MAX);
//変更場所(終り)

	fclose(OC);
}
/****************************** パラメータの吐き出し終わり ******************************/
