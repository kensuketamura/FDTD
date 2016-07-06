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


#define _CRT_SECURE_NO_WARNINGS //	警告を発生させないようにする

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//baba lab
#include <direct.h>

//kuramitsu lab
//#include <sys/stat.h>
//#include <sys/types.h>

#include <stdlib.h>
#include <string>
#include "mpi.h"
#include "grobal_function.cpp"
#include "parameter.h"
#include "module0.h"

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
void output_field(char*); 
void output_field_write(char *); 
void output_model(); 
void mcircle(int, int, int);
void rightquartercircle1(int, int, int, double); 
void leftquartercircle1(int, int, int, double); 
void rightquartercircle2(int, int, int, double); 
void leftquartercircle2(int, int, int, double); 




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
	MPI_Comm_size (MPI_COMM_WORLD, &isize); 
	MPI_Comm_rank (MPI_COMM_WORLD, &irank); 
	MPI_Get_processor_name (processor_name, &namelen); 

	if (isize != ISIZE){
		printf ("MPIで設定した計算機の台数(%d)がプログラム中の値と一致しません．\n終了します\n", ISIZE); 
		return 0; 
	}

	printf ("%d分割並列処理スタート\n", isize); 
	printf ("Process %d on %s\n", irank, processor_name); 

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
			//_strtime(time); 
			//fprintf(fpparameter, "Start Time:\t %s\n", time); 
			s_time = MPI_Wtime(); 
		}

		// 電磁界計算
		for(n = 1 ; n <= Nmax; n++){

			// 時間ステップ数の表示
			if(n % Ncut == 0){
				//_strtime(time); 
				printf("n = %d, \t\t", n); 
				//printf("time = %s\n", time); 
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

			// 一度同期をとる
			MPI_Barrier(MPI_COMM_WORLD); 
#endif
			if(n == 1) {
				output_model(); 		// モデルの出力
				set_epsilon(); 			// 誘電率の割り当て
			}
		}

		if (irank == IRANK_MIN){
			//_strtime(time);
			//fprintf(fpparameter, "End Time:\t %s\n", time); 	/*計算終了時刻の出力*/
			//時刻の出力
			e_time = MPI_Wtime(); 
			printf ("\ntime = %f\n", e_time - s_time); 
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

	//baba lab
	_mkdir(strcpy(dir_name, dir_name_def)); 		// 振り分けできるかテスト
	
	//kuramitsu lab
	//mkdir(strcpy(dir_name, dir_name_def), 0755); 		// 振り分けできるかテスト

	if (irank == IRANK_MIN){
		fpparameter = fopen (strcat(strcpy(dir_name, dir_name_def), "/Parameter.txt"), "w"); 
		//allmodel_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_xy.txt"), "w"); 
		//fpallepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/All_Epsilon_xy.txt"), "w"); 
	}


	model_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xy.txt"), "w"); 		// 振り分けできるかテスト
	model_yz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_yz.txt"), "w"); 
	model_xz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xz.txt"), "w"); 


	fpepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_xy.txt"), "w"); 
	fpepsilony = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz.txt"), "w"); 
	fpepsilonz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx.txt"), "w"); 
	//fpepsilony2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz2.txt"), "w"); 
	//fpepsilonz2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx2.txt"), "w"); 

	
}


/*出力用ファイルを閉じる*/
void file_close(){

	if (irank == IRANK_MIN){
		fclose(fpparameter);  
		//fclose(fpallepsilonx); 
	}

	fclose(model_xy);
	fclose(model_yz); 
	fclose(model_xz); 

	fclose(fpepsilonx); 
	fclose(fpepsilony); 
	fclose(fpepsilonz); 
	//fclose(fpepsilony2); 
	//fclose(fpepsilonz2); 

}


// 計算用パラメータの設定と出力
void parameter(char* dir_name){

	if (irank == IRANK_MIN){

		fprintf(fpparameter, "XMAX_ALL = %d\n", XMAX_ALL); 
		fprintf(fpparameter, "YMAX_ALL = %d\n", YMAX_ALL); 
		fprintf(fpparameter, "ZMAX_ALL = %d\n", ZMAX_ALL); 
		fprintf(fpparameter, "Nodes = %d\n", NODE); 
		fprintf(fpparameter, "Cell Size [nm] = %d\n", CELL_SIZE); 
		fprintf(fpparameter, "Time Step [s] = %e\n", dt); 
		fprintf(fpparameter, "Final Time Step = %d\n", Nmax); 
		fprintf(fpparameter, "Final Time [s] = %e\n", (double) Nmax * dt); 
		fprintf(fpparameter, "\n"); 

		fprintf(fpparameter, "Upper Clad Index = %lf\n", n_clad); 
		fprintf(fpparameter, "Slab Index = %lf\n", n_core); 
		fprintf(fpparameter, "Upper Height [nm] = %d\n", CLAD_HEIGHT1); 
		fprintf(fpparameter, "Slab Height [nm] = %d\n", SLAB_HEIGHT); 
		fprintf(fpparameter, "\n"); 

		fprintf(fpparameter, "Exctation [nm] = %d\n", EXCT_LEN);
		fprintf(fpparameter, "Exctation to Observation [nm] = %d\n", EXCT_OBSE_LEN); 

		fprintf(fpparameter, "\n"); 
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


	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				// セルの目印
				cell[x][y][z] = 0; 
			}	
		}	
	}


	/****************************** Murの吸収境界条件 ******************************/

	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax+1; z++){
			Exn2y00[x][z] = 0.0; 
			Exn1y00[x][z] = 0.0; 
			Exn2y01[x][z] = 0.0; 
			Exn1y01[x][z] = 0.0; 

		}	
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			Exn2z00[x][y] = 0.0; 
			Exn1z00[x][y] = 0.0; 
			Exn2z01[x][y] = 0.0; 
			Exn1z01[x][y] = 0.0; 

		}	
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Eyn2z00[x][y] = 0.0; 
			Eyn1z00[x][y] = 0.0; 
			Eyn2z01[x][y] = 0.0; 
			Eyn1z01[x][y] = 0.0; 

		}	
	}
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax+1; z++){
			Eyn2x00[y][z] = 0.0;
			Eyn1x00[y][z] = 0.0;
			Eyn2x01[y][z] = 0.0;
			Eyn1x01[y][z] = 0.0; 
			Eyn2xm1[y][z] = 0.0;
			Eyn1xm1[y][z] = 0.0;
			Eyn2xm0[y][z] = 0.0;
			Eyn1xm0[y][z] = 0.0; 
		}	
	}
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax; z++){
			Ezn2y00[x][z] = 0.0;
			Ezn1y00[x][z] = 0.0;
			Ezn2y01[x][z] = 0.0;
			Ezn1y01[x][z] = 0.0; 

		}	
	}
	for(y = 0; y <= ymax; y++){
		for(z = 0; z <= zmax-1; z++){
			Ezn2x00[y][z] = 0.0;
			Ezn1x00[y][z] = 0.0;
			Ezn2x01[y][z] = 0.0;
			Ezn1x01[y][z] = 0.0; 
			Ezn2xm1[y][z] = 0.0;
			Ezn1xm1[y][z] = 0.0;
			Ezn2xm0[y][z] = 0.0;
			Ezn1xm0[y][z] = 0.0; 
		}	
	}

	/****************************** Murの吸収境界条件 ******************************/
}


// モデルの設定
void modeling(){

	int n_temp; 		//屈折率の値保存用
	double epsilon_temp; 		//誘電率の値保存用

	/****************************** スラブの形成 ******************************/

	//for(x = 0; x < xmax_all+1; x++){
	//	for(y = 0; y < ymax_all; y++){
	//		for(z = 0; z < zmax_all; z++){		

	//			n_temp = CLAD; 
	//			epsilon_temp = epsilon2; 

	//			if(z < air_hc){			//空気層に設定
	//			}
	//			if(z >= air_hc && z < (air_hc + intCladHeight1) ){			//上部クラッドに設定
	//			}
	//			if(z >= (air_hc + intCladHeight1) && z < (air_hc + intCladHeight1 + intSlabHeigPer)){	// スラブに設定
	//				n_temp = CORE; 
	//				epsilon_temp = epsilon1; 
	//			}

	//			ALL_cell[x][y][z] = n_temp; 
	//			ALL_epsilonx[x][y][z] = epsilon_temp; 
	//			ALL_epsilony[x][y][z] = epsilon_temp; 
	//			ALL_epsilonz[x][y][z] = epsilon_temp; 
	//		}
	//	}
	//}
	/****************************** スラブの形成 ******************************/


	/****************************** 入出射細線導波路 ******************************/
	//int intPcwSislabOffset;

	//// 全面スラブになっているので，細線以外の部分を空気に変更
	//if (PCW_SiSLAB_OFFSET != 0){
	//	intPcwSislabOffset = INT_DIV(PCW_SiSLAB_OFFSET, CELL_SIZE);
	//}
	//else{
	//	intPcwSislabOffset = 0;
	//}

	//for (z = zmax_all - intSlabHeigPer; z < (zmax_all + 1); z++){
	//	for (y = 0; y < ymax_all - intWireWid_2; y++){

	//		// 入射
	//		if (PCW_SiSLAB_OFFSET != 0){
	//		}
	//		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++){		// 配列の引数に使用するので-1
	//			ALL_cell[x][y][z] = CLAD; 
	//			ALL_epsilonx[x][y][z] = epsilon2; 
	//			ALL_epsilony[x][y][z] = epsilon2; 
	//			ALL_epsilonz[x][y][z] = epsilon2; 
	//		}

	//		//// 出射
	//		//if (PCW_SiSLAB_OFFSET != 0){
	//		//}
	//		//for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++){
	//		//	ALL_cell[x][y][z] = CLAD; 
	//		//	ALL_epsilonx[x][y][z] = epsilon2; 
	//		//	ALL_epsilony[x][y][z] = epsilon2; 
	//		//	ALL_epsilonz[x][y][z] = epsilon2; 
	//		//}
	//	}
	//}
	/****************************** 入出射細線導波路 ******************************/



	/****************************** 対称境界部分の誘電率の設定 ******************************/
	//後で考える

	/*for(x = 0; x < xmax_all+1; x++){
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
	}*/

	/****************************** 対称境界部分の誘電率の設定 ******************************/



	/****************************** 各ノードにモデルを分割 ******************************/
	//if(irank != IRANK_MAX){
	//	for(x = 0; x < xmax+1; x++){
	//		for(y = 0; y < ymax+1; y++){
	//			for(z = 0; z < zmax+1; z++){
	//				cell[x][y][z] = ALL_cell[irank*(xmax-1)+x][y][z]; 
	//				epsilonx[x][y][z] = ALL_epsilonx[irank*(xmax-1)+x][y][z]; 
	//				epsilony[x][y][z] = ALL_epsilony[irank*(xmax-1)+x][y][z]; 
	//				epsilonz[x][y][z] = ALL_epsilonz[irank*(xmax-1)+x][y][z]; 
	//			}
	//		}
	//	}
	//}
	//else{
	//	for(x = 0; x < xmax+1; x++){
	//		for(y = 0; y < ymax+1; y++){
	//			for(z = 0; z < zmax+1; z++){
	//				cell[x][y][z] = ALL_cell[irank*(xmax)+x][y][z]; 
	//				epsilonx[x][y][z] = ALL_epsilonx[irank*(xmax)+x][y][z]; 
	//				epsilony[x][y][z] = ALL_epsilony[irank*(xmax)+x][y][z]; 
	//				epsilonz[x][y][z] = ALL_epsilonz[irank*(xmax)+x][y][z]; 
	//			}
	//		}
	//	}
	//}
	/****************************** 各ノードにモデルを分割 ******************************/

	/****************************** 分割したモデル ******************************/
	if(irank != IRANK_MAX){
		for(x = 0; x < xmax+1; x++){
			for(y = 0; y < ymax+1; y++){
				for(z = 0; z < zmax+1; z++){

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
					cell[x][y][z] = n_temp; 
					epsilonx[x][y][z] = epsilon_temp; 
					epsilony[x][y][z] = epsilon_temp; 
					epsilonz[x][y][z] = epsilon_temp; 
				}
			}
		}
	}
	else{
		for(x = 0; x < xmax+1; x++){
			for(y = 0; y < ymax+1; y++){
				for(z = 0; z < zmax+1; z++){

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
					cell[x][y][z] = n_temp; 
					epsilonx[x][y][z] = epsilon_temp; 
					epsilony[x][y][z] = epsilon_temp; 
					epsilonz[x][y][z] = epsilon_temp;  

				}
			}
		}
	}

	//フォトニック結晶

 //   int x2,y2;
	//if(irank == IRANK_MIN){
	//	for(x2 = 0; x2 < Row_x; x2++){		
	//		for(y2 = 0; y2 < Row_y; y2++){
	//			if(y2 % 2 == 0){
	//				y = ymax - Wid - intPitchY*y2;
	//				x = intPitchX * x2 + intPitchX * 0.5;
	//			}else{
	//				x = intPitchX * x2;
	//				y = ymax - Wid - intPitchY*y2;
	//			}
	//			for(z = 0; z < zmax+1; z++){
	//				mcircle(x, y, z);				
	//			}
	//		}
	//	}
	//}


	/****************************** 分割したモデル ******************************/


	/****************************** 共通パラメータの設定 ******************************/

	// 励振点，観測面の設定 (XMAXは "のりしろ" 部分を含めていることに注意)
	intExctPortNum = intExctLen / (XMAX - 1);

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
				//fprintf (fpallepsilonx, "%e\t", epsilon_xy[x][y]); 
				fprintf (fpepsilonx, "%e\t", epsilon_xy[x][y]); 
			}
			//fprintf (fpallepsilonx, "\n"); 
			fprintf (fpepsilonx, "\n"); 
		}
	}
	if(irank != IRANK_MIN){
		for(x = 1; x<xmax; x++){
			for(y = 0; y < ymax+1; y++){
				fprintf (fpepsilonx, "%e\t", epsilon_xy[x][y]); 
			}
			fprintf (fpepsilonx, "\n"); 
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
			//fprintf (fpallepsilonx, "%e\t", epsilon_xy[x][y]); 
			fprintf (fpepsilonx, "%e\t", epsilon_xy[x][y]); 
		}
		fprintf (fpallepsilonx, "\n"); 
		fprintf (fpepsilonx, "\n"); 
	}

	/****************************** モデル確認時 ******************************/
#endif

	//ZX平面 (Y:境界面)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			epsilon_zx[x][z] = epsilonz[x][ymax][z]; 
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			fprintf(fpepsilonz, "%e\t", epsilon_zx[x][z]); 
		}
		fprintf(fpepsilonz, "\n"); 
	}
	//ZX平面 (Y:中心)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			epsilon_zx2[x][z] = epsilonz[x][ymax/2][z]; 
		}
	}

	// ファイルポインタを閉じる
	if (irank == IRANK_MIN){
		//fclose(fpallepsilonx);
	}
	fclose(fpepsilonx); 
	fclose(fpepsilony); 
	fclose(fpepsilonz); 

}


// 励振関数
void source_func(){

	int x, y, z;

	if(irank == intExctPortNum){

		// 励振点の設定
		x = intExctLenPart;

		for(y = ex_y_st; y < ex_y_ed; y++){
			for(z = ex_z_st; z < ex_z_ed; z++){
#if _EXITATION_FUNC	// CW励振


				//面内正弦分布励振の場合 ←解析空間が偶数セルか奇数セルかで励振が異なるのでその都度注意

				// スラブ厚の半分のセル数:偶数 導波路幅の半分のセル数:偶数
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st - 1)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st - 1)) * sin(omega0*n*dt); 

				// スラブ厚の半分のセル数:奇数 導波路幅の半分のセル数:奇数
				Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt); // 01

#else	// Gaussian励振

				//Hz[x][(YMAX+1)/2][intSlabCen] += 1000 * cos(omega0*(n-Npeak)*dt) * exp(-(SQ(sigma*dt*(n-Npeak))/2)); 		
				Hz[x][y][z] += 1000 * cos(omega0*(n-Npeak)*dt) * exp(-(SQ(sigma*dt*(n-Npeak))/2)) * cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)); 
#endif
			}
		}
	}


	/****************************** 磁界の対称境界条件(4回対称) ******************************/

	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax; z++){
			Hx[x][ymax][z] = Hx[x][ymax-1][z];		// 偶関数
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			Hz[x][ymax][z] = Hz[x][ymax-1][z];		// 偶関数
		}
	}
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			Hy[x][y][zmax] = -Hy[x][y][zmax-1];		// 奇関数
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Hx[x][y][zmax] = -Hx[x][y][zmax-1];		// 奇関数
		}
	}

	/****************************** 磁界の対称境界条件(4回対称) ******************************/
}




void calc_efield(){

	double dex, dey, dez;
	double cnstEx, cnstEy, cnstEz;

	// Ex
	for(x = 0; x < xmax; x++){
		for(y = 1; y < ymax+1; y++){		// Exはy軸に対して奇関数
			for(z = 1; z < zmax+1; z++){
				cnstEx = dt / epsilonx[x][y][z]; 
				dex = ( (Hz[x][y][z] - Hz[x][y-1][z]) / dy) - ( (Hy[x][y][z] - Hy[x][y][z-1]) / dz); 
				Ex[x][y][z] = Ex[x][y][z] + cnstEx * dex; 
			}
		}
	}

	// Ey
	for(x = 1; x < xmax; x++){	
		for(y = 0; y < ymax; y++){
			for(z = 1; z < zmax+1; z++){	
				cnstEy = dt / epsilony[x][y][z]; 
				dey = ( (Hx[x][y][z] - Hx[x][y][z-1]) / dz)-( (Hz[x][y][z] - Hz[x-1][y][z]) / dx); 
				Ey[x][y][z] = Ey[x][y][z] + cnstEy * dey; 
			}
		}
	}

	// Ez	
	for(x = 1; x < xmax; x++){	
		for(y = 1; y < ymax+1; y++){		// Ezはy軸に対して奇関数
			for(z = 0; z < zmax; z++){
				cnstEz = dt / epsilonz[x][y][z]; 
				dez = ( (Hy[x][y][z] - Hy[x-1][y][z]) / dx) - ( (Hx[x][y][z] - Hx[x][y-1][z]) / dy); 
				Ez[x][y][z] = Ez[x][y][z] + cnstEz * dez; 
			}
		}
	}


	/****************************** 電界の対称境界条件 ******************************/

	// 境界面で反対称となる電界成分の境界面上の値を0としている
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			Ex[x][ymax][z] = 0.0;		// 奇関数
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax; z++){
			Ez[x][ymax][z] = 0.0;		// 奇関数
		}
	}
	/****************************** 電界の対称境界条件 ******************************/
}



void calc_hfield(){	

	double dhx, dhy, dhz;

	// Hx
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				dhx = ( (Ey[x][y][z+1] - Ey[x][y][z]) / dz) - ( (Ez[x][y+1][z] - Ez[x][y][z]) / dy); 
				Hx[x][y][z] = Hx[x][y][z] + cnstHxyz * dhx; 
			}
		}
	}

	// Hy
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax; z++){
				dhy = ( (Ez[x+1][y][z] - Ez[x][y][z]) / dx) - ( (Ex[x][y][z+1] - Ex[x][y][z]) / dz); 
				Hy[x][y][z] = Hy[x][y][z] + cnstHxyz * dhy; 
			}
		}
	}

	// Hz
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax+1; z++){
				dhz = ((Ex[x][y+1][z] - Ex[x][y][z]) / dy) - ((Ey[x+1][y][z] - Ey[x][y][z]) / dx); 
				Hz[x][y][z] = Hz[x][y][z] + cnstHxyz * dhz; 
			}
		}
	}
}


// Mur2次，1次の吸収境界条件から端面の計算
void absorpt_bound_condition(){

	/****************************** 対称境界条件 ******************************/			

	// 2回対称
	//for(z = 0; z < zmax+1; z++){
	//	Exn1y00[xmax][z] = Exn1y00[xmax-1][z];
	//	Exn1y01[xmax][z] = Exn1y01[xmax-1][z];
	//	Exn1ym0[xmax][z] = Exn1ym0[xmax-1][z];
	//	Exn1ym1[xmax][z] = Exn1ym1[xmax-1][z];
	//}
	//for(y = 0; y < ymax+1; y++){
	//	Exn1z00[xmax][y] = Exn1z00[xmax-1][y];
	//	Exn1z01[xmax][y] = Exn1z01[xmax-1][y];
	//	Exn1zm0[xmax][y] = Exn1zm0[xmax-1][y];
	//	Exn1zm1[xmax][y] = Exn1zm1[xmax-1][y];
	//}
	//for(y = 0; y < ymax; y++){
	//	Eyn1z00[xmax+1][y] = -Eyn1z00[xmax-1][y];
	//	Eyn1z01[xmax+1][y] = -Eyn1z01[xmax-1][y];
	//	Eyn1zm0[xmax+1][y] = -Eyn1zm0[xmax-1][y];
	//	Eyn1zm1[xmax+1][y] = -Eyn1zm1[xmax-1][y];
	//}
	//for(z = 0; z < zmax; z++){
	//	Ezn1y00[xmax+1][z] = -Ezn1y00[xmax-1][z];
	//	Ezn1y01[xmax+1][z] = -Ezn1y01[xmax-1][z];
	//	Ezn1ym0[xmax+1][z] = -Ezn1ym0[xmax-1][z];
	//	Ezn1ym1[xmax+1][z] = -Ezn1ym1[xmax-1][z];
	//}

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
		//Eyn1zm0[x][ymax] = Eyn1zm0[x][ymax-1];
		//Eyn1zm1[x][ymax] = Eyn1zm1[x][ymax-1];
	}
	for(x = 0; x < xmax; x++){
		Exn1z00[x][ymax+1] = -Exn1z00[x][ymax-1];
		Exn1z01[x][ymax+1] = -Exn1z01[x][ymax-1];
		//Exn1zm0[x][ymax+1] = -Exn1zm0[x][ymax-1];
		//Exn1zm1[x][ymax+1] = -Exn1zm1[x][ymax-1];
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
		//Ezn1ym0[x][zmax] = -Ezn1ym0[x][zmax-1];
		//Ezn1ym1[x][zmax] = -Ezn1ym1[x][zmax-1];
	}
	for(x = 0; x < xmax; x++){
		Exn1y00[x][zmax+1] = Exn1y00[x][zmax-1];
		Exn1y01[x][zmax+1] = Exn1y01[x][zmax-1];
		//Exn1ym0[x][zmax+1] = Exn1ym0[x][zmax-1];
		//Exn1ym1[x][zmax+1] = Exn1ym1[x][zmax-1];
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
				//fprintf(allmodel_xy, "%d\t", cell_xy[x][y]); 
				fprintf(model_xy, "%d\t", cell_xy[x][y]); 
			}
	/*		fprintf(allmodel_xy, "\n"); */
			fprintf(model_xy, "\n"); 
		}
		fclose(model_xy);
	}

	// それぞれ分割部のモデル
	if(irank != IRANK_MIN){
		for(x = 1; x < xmax; x++){
			for(y = 0; y < ymax; y++){
				fprintf(model_xy, "%d\t", cell_xy[x][y]); 
			}
			fprintf(model_xy, "\n"); 
		}

		fclose(model_xy); 
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
			//fprintf(allmodel_xy, "%d\t", cell_xy[x][y]); 
			fprintf(model_xy, "%d\t", cell_xy[x][y]); 
		}
		//fprintf(allmodel_xy, "\n"); 
		fprintf(model_xy, "\n"); 
	}
	fclose(model_xy); 

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

	printf("n = %d\n", n); 

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
				fprintf(HZ1, "%e\t", field_xy[x][y]); 
			}
			fprintf(HZ1, "\n"); 
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
						fprintf(HZ1, "%e\t", field_xy[x][y]); 
					}
					fprintf(HZ1, "\n"); 
				}
			}
			else{
				MPI_Recv(&field_xy[0][0], (xmax)*(ymax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status); 
				for(x = 1; x < xmax; x++){
					for(y = 0; y < ymax; y++){
						fprintf(HZ1, "%e\t", field_xy[x][y]); 
					}
					fprintf(HZ1, "\n"); 
				}
			}
		}

		// ファイルポインタを閉じる
		fclose(HZ1); 		
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



void mcircle(int x_circ, int y_circ, int z_circ){

	double R; 

	//半径セル数の計算
	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています

	rightquartercircle1(x_circ, y_circ, z_circ, R); 
	leftquartercircle1(x_circ-1, y_circ, z_circ, R); 
	rightquartercircle2(x_circ, y_circ-1, z_circ, R); 
	leftquartercircle2(x_circ-1, y_circ-1, z_circ, R); 
}


void halfcircle(int x_circ, int y_circ, int z_circ){

	double R; 

	//半径セル数の計算
	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています

	rightquartercircle2(x_circ, y_circ-1, z_circ, R); 
	leftquartercircle2(x_circ-1, y_circ-1, z_circ, R); 
}


void rightquartercircle1(int x_circ, int y_circ, int z_circ, double R){

	int x, y, Ie, Je; 
	double r; 

	Ie = (int) (x_circ+R-1); 
	Je = (int) (y_circ+R-1); 
	for(x = x_circ; x <= Ie; x++){
		for(y = y_circ; y <= Je; y++){
			r = sqrt(double((x-x_circ+1) * (x-x_circ+1) + (y-y_circ+1) * (y-y_circ+1))) - 0.5; 
			if(r <= R){
				cell[x][y][z_circ] = CIRCLE_REF_INDEX; 
				epsilonx[x][y][z_circ] = epsilon2; 
				epsilony[x][y][z_circ] = epsilon2; 
				epsilonz[x][y][z_circ] = epsilon2; 
			}
		}
	}
}


void leftquartercircle1(int x_circ, int y_circ, int z_circ, double R){

	int x, y, Ie, Je; 
	double r; 

	Ie = (int) (x_circ-R+1); 
	Je = (int) (y_circ+R-1); 
	for(x = x_circ; x >= Ie; x--){
		for(y = y_circ; y <= Je; y++){
			r = sqrt(double((x-x_circ-1) * (x-x_circ-1) + (y-y_circ+1) * (y-y_circ+1))) - 0.5; 
			if(r <= R){
				cell[x][y][z_circ] = CIRCLE_REF_INDEX; 
				epsilonx[x][y][z_circ] = epsilon2; 
				epsilony[x][y][z_circ] = epsilon2; 
				epsilonz[x][y][z_circ] = epsilon2; 
			}
		}
	}
}


void rightquartercircle2(int x_circ, int y_circ, int z_circ, double R){

	int x, y, Ie, Je; 
	double r; 

	Ie = (int) (x_circ+R-1); 
	Je = (int) (y_circ-R+1); 
	for(x = x_circ; x <= Ie; x++){
		for(y = y_circ; y >= Je; y--){
			r = sqrt(double((x-x_circ+1) * (x-x_circ+1) + (y-y_circ-1) * (y-y_circ-1))) - 0.5; 
			if(r <= R){
				cell[x][y][z_circ] = CIRCLE_REF_INDEX; 
				epsilonx[x][y][z_circ] = epsilon2; 
				epsilony[x][y][z_circ] = epsilon2; 
				epsilonz[x][y][z_circ] = epsilon2; 

			}
		}
	}
}

void leftquartercircle2(int x_circ, int y_circ, int z_circ, double R){

	int x, y, Ie, Je; 
	double r; 

	Ie = (int) (x_circ-R+1); 
	Je = (int) (y_circ-R+1); 
	for(x = x_circ; x >= Ie; x--){
		for(y = y_circ; y >= Je; y--){
			r = sqrt(double((x-x_circ-1) * (x-x_circ-1) + (y-y_circ-1) * (y-y_circ-1))) - 0.5; 
			if(r <= R){
				cell[x][y][z_circ] = CIRCLE_REF_INDEX; 
				epsilonx[x][y][z_circ] = epsilon2; 
				epsilony[x][y][z_circ] = epsilon2; 
				epsilonz[x][y][z_circ] = epsilon2; 
	
			}
		}
	}
}



