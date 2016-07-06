/*
3����_FDTD�@�ɂ��d���E��� ver. 2.01
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

#define _FDTD 1		// FDTD�v�Z			0 : ���f���|���o��(�v���v���Z�b�T�ŃR���p�C����ύX������)
//										1 : �v�Z���s

#define _BAND_CALCULATION 0			// �v�Z�̎�� �o���h�v�Z
#define _PROPAGATION_CALCULATION 1	// �v�Z�̎�� �`���v�Z

#define _CALCULATION_TYPE _PROPAGATION_CALCULATION	// �v�Z�̎��

#define _EXITATION_FUNC 1	// ��U�֐��̎��		0 : Gaussian 
//													1 : CW


#define _CRT_SECURE_NO_WARNINGS //	�x���𔭐������Ȃ��悤�ɂ���

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

//�T�u���[�e�B��
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

	// MPI�ɂ��ʐM�̊J�n
	MPI_Init (&argc, &argv); 
	MPI_Comm_size (MPI_COMM_WORLD, &isize); 
	MPI_Comm_rank (MPI_COMM_WORLD, &irank); 
	MPI_Get_processor_name (processor_name, &namelen); 

	if (isize != ISIZE){
		printf ("MPI�Őݒ肵���v�Z�@�̑䐔(%d)���v���O�������̒l�ƈ�v���܂���D\n�I�����܂�\n", ISIZE); 
		return 0; 
	}

	printf ("%d�������񏈗��X�^�[�g\n", isize); 
	printf ("Process %d on %s\n", irank, processor_name); 

	// �ׂ̌v�Z�@�̔ԍ��̎w��
	left = irank - 1; 
	if(irank == IRANK_MIN){
		left = MPI_PROC_NULL; 
	}
	right = irank + 1; 
	if(irank == IRANK_MAX){
		right = MPI_PROC_NULL; 
	}

	// dir_name (��U�g��) �̔z�񒷂����J��Ԃ�
	for(int dir_count = 0; dir_count < (sizeof(dir_name) / sizeof(dir_name[0]) ); dir_count++){

		initialize_matrix(); 						// �z��̏�����
		modeling(); 								// ���f���̐ݒ�
		file_open(dir_name[dir_count]); 			// �t�@�C�����J��
		parameter(dir_name[dir_count]); 			// �p�����[�^�̐ݒ�Əo��


		// �v�Z�J�n�����̏o��
		if (irank == IRANK_MIN){
			//_strtime(time); 
			//fprintf(fpparameter, "Start Time:\t %s\n", time); 
			s_time = MPI_Wtime(); 
		}

		// �d���E�v�Z
		for(n = 1 ; n <= Nmax; n++){

			// ���ԃX�e�b�v���̕\��
			if(n % Ncut == 0){
				//_strtime(time); 
				printf("n = %d, \t\t", n); 
				//printf("time = %s\n", time); 
			}

			// ��U�֐��̐ݒ�
			source_func(); 

#if _FDTD

			// ��x�������Ƃ�(�����̓m�[�h�Ԃő��x�ɂ΂������������)
			MPI_Barrier (MPI_COMM_WORLD); 	

			// �d�E�̌v�Z
			calc_efield(); 							

			// �z�����E�����ɂ��[�ʂ̌v�Z
			absorpt_bound_condition();

			// ��x�������Ƃ�
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

			// �d�E�̕ۑ�
			saving_electric_field();

			// ���E�̌v�Z
			calc_hfield();

			MPI_Sendrecv( &Hy[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send, 
				&Hy[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status); 
			MPI_Sendrecv( &Hz[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send, 
				&Hz[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);

			// �t�B�[���h�̏o��
			output_field (dir_name[dir_count]); 

			// �|�C���e�B���O�p���[�v�Z�Əo��

			// ��x�������Ƃ�
			MPI_Barrier(MPI_COMM_WORLD); 
#endif
			if(n == 1) {
				output_model(); 		// ���f���̏o��
				set_epsilon(); 			// �U�d���̊��蓖��
			}
		}

		if (irank == IRANK_MIN){
			//_strtime(time);
			//fprintf(fpparameter, "End Time:\t %s\n", time); 	/*�v�Z�I�������̏o��*/
			//�����̏o��
			e_time = MPI_Wtime(); 
			printf ("\ntime = %f\n", e_time - s_time); 
		}

		file_close(); 			// �t�@�C�������
	}

	//MPI_Finalize(); 			// MPI���I������
#elif _CALCULATION_TYPE == _BAND_CALCULATION

#endif
}



// �o�͗p�t�@�C�����J��
void file_open(char* dir_name_def){
	char dir_name[40]; 

	//baba lab
	_mkdir(strcpy(dir_name, dir_name_def)); 		// �U�蕪���ł��邩�e�X�g
	
	//kuramitsu lab
	//mkdir(strcpy(dir_name, dir_name_def), 0755); 		// �U�蕪���ł��邩�e�X�g

	if (irank == IRANK_MIN){
		fpparameter = fopen (strcat(strcpy(dir_name, dir_name_def), "/Parameter.txt"), "w"); 
		//allmodel_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_xy.txt"), "w"); 
		//fpallepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/All_Epsilon_xy.txt"), "w"); 
	}


	model_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xy.txt"), "w"); 		// �U�蕪���ł��邩�e�X�g
	model_yz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_yz.txt"), "w"); 
	model_xz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xz.txt"), "w"); 


	fpepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_xy.txt"), "w"); 
	fpepsilony = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz.txt"), "w"); 
	fpepsilonz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx.txt"), "w"); 
	//fpepsilony2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz2.txt"), "w"); 
	//fpepsilonz2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx2.txt"), "w"); 

	
}


/*�o�͗p�t�@�C�������*/
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


// �v�Z�p�p�����[�^�̐ݒ�Əo��
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

	// ��U�֐��萔�̐ݒ�
	lambda = atof(dir_name) * 1e-9; 		
	omega0 = 2.0*PI*C0/lambda; 
	sigma = omega0 * delta_omega;
}

/*�z��̏�����*/
void initialize_matrix(){

	//�e�m�[�h�̍��W
	if(irank != IRANK_MAX){				
		xmax = XMAX; 
		ymax = YMAX; 
		zmax = ZMAX; 
	}

	//�Ō�̃m�[�h�����̂肵��s�v�Ȃ̂�x������1�Z��������
	if(irank == IRANK_MAX){				
		xmax = XMAX - 1; 
		ymax = YMAX; 
		zmax = ZMAX; 
	}

	// ��͋�Ԃ̍ő�l
	xmax_all = XMAX_ALL; 
	ymax_all = YMAX_ALL; 
	zmax_all = ZMAX_ALL; 

	// ��͋�Ԃ̒��S���W
	x_cen = xmax/2; 
	y_cen = ymax/2; 
	z_cen = zmax/2; 	

	//���f���̒��S�Ɖ�͋�Ԃ̒��S�͂P�Z��������Ă���̂ŗv����
	x_model_cen = x_cen + 1; 
	y_model_cen = y_cen + 1; 

	int x, y, z;

	for (x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax+1; z++){
				// �d�E
				Ex[x][y][z] = 0.0; 
				Ey[x][y][z] = 0.0;
				Ez[x][y][z] = 0.0; 

				// ���E
				Hx[x][y][z] = 0.0;
				Hy[x][y][z] = 0.0; 
				Hz[x][y][z] = 0.0; 
			}	
		}	
	}


	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax+1; z++){
				// �U�d��(��������K�v���Ȃ��悤�ȁD�D�D)
				epsilonx[x][y][z] = EPSILON0; 
				epsilony[x][y][z] = EPSILON0; 
				epsilonz[x][y][z] = EPSILON0; 
			}	
		}	
	}


	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				// �Z���̖ڈ�
				cell[x][y][z] = 0; 
			}	
		}	
	}


	/****************************** Mur�̋z�����E���� ******************************/

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

	/****************************** Mur�̋z�����E���� ******************************/
}


// ���f���̐ݒ�
void modeling(){

	int n_temp; 		//���ܗ��̒l�ۑ��p
	double epsilon_temp; 		//�U�d���̒l�ۑ��p

	/****************************** �X���u�̌`�� ******************************/

	//for(x = 0; x < xmax_all+1; x++){
	//	for(y = 0; y < ymax_all; y++){
	//		for(z = 0; z < zmax_all; z++){		

	//			n_temp = CLAD; 
	//			epsilon_temp = epsilon2; 

	//			if(z < air_hc){			//��C�w�ɐݒ�
	//			}
	//			if(z >= air_hc && z < (air_hc + intCladHeight1) ){			//�㕔�N���b�h�ɐݒ�
	//			}
	//			if(z >= (air_hc + intCladHeight1) && z < (air_hc + intCladHeight1 + intSlabHeigPer)){	// �X���u�ɐݒ�
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
	/****************************** �X���u�̌`�� ******************************/


	/****************************** ���o�ˍא����g�H ******************************/
	//int intPcwSislabOffset;

	//// �S�ʃX���u�ɂȂ��Ă���̂ŁC�א��ȊO�̕�������C�ɕύX
	//if (PCW_SiSLAB_OFFSET != 0){
	//	intPcwSislabOffset = INT_DIV(PCW_SiSLAB_OFFSET, CELL_SIZE);
	//}
	//else{
	//	intPcwSislabOffset = 0;
	//}

	//for (z = zmax_all - intSlabHeigPer; z < (zmax_all + 1); z++){
	//	for (y = 0; y < ymax_all - intWireWid_2; y++){

	//		// ����
	//		if (PCW_SiSLAB_OFFSET != 0){
	//		}
	//		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++){		// �z��̈����Ɏg�p����̂�-1
	//			ALL_cell[x][y][z] = CLAD; 
	//			ALL_epsilonx[x][y][z] = epsilon2; 
	//			ALL_epsilony[x][y][z] = epsilon2; 
	//			ALL_epsilonz[x][y][z] = epsilon2; 
	//		}

	//		//// �o��
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
	/****************************** ���o�ˍא����g�H ******************************/



	/****************************** �Ώ̋��E�����̗U�d���̐ݒ� ******************************/
	//��ōl����

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

	/****************************** �Ώ̋��E�����̗U�d���̐ݒ� ******************************/



	/****************************** �e�m�[�h�Ƀ��f���𕪊� ******************************/
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
	/****************************** �e�m�[�h�Ƀ��f���𕪊� ******************************/

	/****************************** �����������f�� ******************************/
	if(irank != IRANK_MAX){
		for(x = 0; x < xmax+1; x++){
			for(y = 0; y < ymax+1; y++){
				for(z = 0; z < zmax+1; z++){

					n_temp = CLAD; 
					epsilon_temp = epsilon2; 

					if(z < air_hc){			//��C�w�ɐݒ�
					}
					if(z >= air_hc && z < (air_hc + intCladHeight1) ){			//�㕔�N���b�h�ɐݒ�
					}
					if(z >= (air_hc + intCladHeight1) && z < (air_hc + intCladHeight1 + intSlabHeigPer)){	// �X���u�ɐݒ�
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

					if(z < air_hc){			//��C�w�ɐݒ�
					}
					if(z >= air_hc && z < (air_hc + intCladHeight1) ){			//�㕔�N���b�h�ɐݒ�
					}
					if(z >= (air_hc + intCladHeight1) && z < (air_hc + intCladHeight1 + intSlabHeigPer)){	// �X���u�ɐݒ�
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

	//�t�H�g�j�b�N����

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


	/****************************** �����������f�� ******************************/


	/****************************** ���ʃp�����[�^�̐ݒ� ******************************/

	// ��U�_�C�ϑ��ʂ̐ݒ� (XMAX�� "�̂肵��" �������܂߂Ă��邱�Ƃɒ���)
	intExctPortNum = intExctLen / (XMAX - 1);

	/****************************** ���ʃp�����[�^�̐ݒ� ******************************/
}



//�U�d���̊��蓖��
void set_epsilon(){

	//�U�d�����z�̏o��(���f���̊m�F)
	int tag1 = 1; 

#if _FDTD

	/****************************** �v�Z���s�� ******************************/
	int node; 

	MPI_Status status; 

	//XY����
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

	/****************************** �v�Z���s�� ******************************/
#else

	/****************************** ���f���m�F�� ******************************/
	char fname[40],dir_name[50];	//�t�@�C�����i�[�ϐ�	

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

	/****************************** ���f���m�F�� ******************************/
#endif

	//ZX���� (Y:���E��)
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
	//ZX���� (Y:���S)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			epsilon_zx2[x][z] = epsilonz[x][ymax/2][z]; 
		}
	}

	// �t�@�C���|�C���^�����
	if (irank == IRANK_MIN){
		//fclose(fpallepsilonx);
	}
	fclose(fpepsilonx); 
	fclose(fpepsilony); 
	fclose(fpepsilonz); 

}


// ��U�֐�
void source_func(){

	int x, y, z;

	if(irank == intExctPortNum){

		// ��U�_�̐ݒ�
		x = intExctLenPart;

		for(y = ex_y_st; y < ex_y_ed; y++){
			for(z = ex_z_st; z < ex_z_ed; z++){
#if _EXITATION_FUNC	// CW��U


				//�ʓ��������z��U�̏ꍇ ����͋�Ԃ������Z������Z�����ŗ�U���قȂ�̂ł��̓s�x����

				// �X���u���̔����̃Z����:���� ���g�H���̔����̃Z����:����
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st - 1)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st - 1)) * sin(omega0*n*dt); 

				// �X���u���̔����̃Z����:� ���g�H���̔����̃Z����:�
				Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt); // 01

#else	// Gaussian��U

				//Hz[x][(YMAX+1)/2][intSlabCen] += 1000 * cos(omega0*(n-Npeak)*dt) * exp(-(SQ(sigma*dt*(n-Npeak))/2)); 		
				Hz[x][y][z] += 1000 * cos(omega0*(n-Npeak)*dt) * exp(-(SQ(sigma*dt*(n-Npeak))/2)) * cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)); 
#endif
			}
		}
	}


	/****************************** ���E�̑Ώ̋��E����(4��Ώ�) ******************************/

	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax; z++){
			Hx[x][ymax][z] = Hx[x][ymax-1][z];		// ��֐�
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			Hz[x][ymax][z] = Hz[x][ymax-1][z];		// ��֐�
		}
	}
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			Hy[x][y][zmax] = -Hy[x][y][zmax-1];		// ��֐�
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Hx[x][y][zmax] = -Hx[x][y][zmax-1];		// ��֐�
		}
	}

	/****************************** ���E�̑Ώ̋��E����(4��Ώ�) ******************************/
}




void calc_efield(){

	double dex, dey, dez;
	double cnstEx, cnstEy, cnstEz;

	// Ex
	for(x = 0; x < xmax; x++){
		for(y = 1; y < ymax+1; y++){		// Ex��y���ɑ΂��Ċ�֐�
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
		for(y = 1; y < ymax+1; y++){		// Ez��y���ɑ΂��Ċ�֐�
			for(z = 0; z < zmax; z++){
				cnstEz = dt / epsilonz[x][y][z]; 
				dez = ( (Hy[x][y][z] - Hy[x-1][y][z]) / dx) - ( (Hx[x][y][z] - Hx[x][y-1][z]) / dy); 
				Ez[x][y][z] = Ez[x][y][z] + cnstEz * dez; 
			}
		}
	}


	/****************************** �d�E�̑Ώ̋��E���� ******************************/

	// ���E�ʂŔ��Ώ̂ƂȂ�d�E�����̋��E�ʏ�̒l��0�Ƃ��Ă���
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			Ex[x][ymax][z] = 0.0;		// ��֐�
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax; z++){
			Ez[x][ymax][z] = 0.0;		// ��֐�
		}
	}
	/****************************** �d�E�̑Ώ̋��E���� ******************************/
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


// Mur2���C1���̋z�����E��������[�ʂ̌v�Z
void absorpt_bound_condition(){

	/****************************** �Ώ̋��E���� ******************************/			

	// 2��Ώ�
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

	// 4��Ώ�
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

	// 8��Ώ�
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
	/****************************** �Ώ̋��E���� ******************************/			




	/****************************** Mur��2���̋z�����E����(Ex) ******************************/			

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

	/****************************** Mur��2���̋z�����E����(Ex) ******************************/			

	double u1xa, u1xc;
	double u2xa, u2xc;

	/****************************** Mur��1���̋z�����E����(Ex) ******************************/			

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


	// ��(Mur��1���̋z�����E����) -- y���ʂ�z���ʂ��炻�ꂼ��Z�o�����l�̕��ϒl�����
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

	/****************************** Mur��1���̋z�����E����(Ex) ******************************/

	double u1by1, u2by1, u3by1, u4by1;
	double u1cy1, u2cy1, u3cy1, u4cy1;
	double u1cy2, u2cy2, u3cy2, u4cy2;

	/****************************** Mur��2���̋z�����E����(Ey) ******************************/

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

	/****************************** Mur��2���̋z�����E����(Ey) ******************************/


	double u2ya, u3ya, u3yb;
	double u2ya1;
	double u2yc1;

	/****************************** Mur��1���̋z�����E����(Ey) ******************************/

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

	// ��(Mur��1���̋z�����E����) --x���ʂ�z���ʂ��炻�ꂼ��Z�o�����l�̕��ϒl�����
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
	/****************************** Mur��1���̋z�����E����(Ey) ******************************/

	double u1az1, u2az1, u3az1, u4az1;
	double u1cz1, u2cz1, u3cz1, u4cz1;
	double u1cz2, u2cz2, u3cz2, u4cz2;

	/****************************** Mur��2���̋z�����E����(Ez) ******************************/

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

	/****************************** Mur��2���̋z�����E����(Ez) ******************************/

	double u1za, u3za, u3zb;
	double u1za1, u1zb1;

	/****************************** Mur��1���̋z�����E����(Ez) ******************************/
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

	// ��(Mur��1���̋z�����E����) --x���ʂ�y���ʂ��炻�ꂼ��Z�o�����l�̕��ϒl�����
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
	/****************************** Mur��1���̋z�����E����(Ez) ******************************/

}



/*�d�E�̕ۑ�*/
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


//���f���̏o��
void output_model(){

	//char cell_xy[XMAX][YMAX]; 
	int tag2 = 2; 
	int x, y, z;

#if _FDTD
	/****************************** �v�Z���s�� ******************************/
	int node; 

	MPI_Status status; 

	z = intSlabCen - 1;

	// XY����
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

	// ���ꂼ�ꕪ�����̃��f��
	if(irank != IRANK_MIN){
		for(x = 1; x < xmax; x++){
			for(y = 0; y < ymax; y++){
				fprintf(model_xy, "%d\t", cell_xy[x][y]); 
			}
			fprintf(model_xy, "\n"); 
		}

		fclose(model_xy); 
	}

	/****************************** �v�Z���s�� ******************************/

#else	
	/****************************** ���f���m�F�� ******************************/
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

	/****************************** ���f���m�F�� ******************************/
#endif

}


void output_field_write(char *dir_name_def){

	char fname[40], dir_name[50]; 	//�t�@�C�����i�[�ϐ�	
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
			field_xy[x][y] = Hz[x][y][ex_z_ed-1]; 	//�S�Ẵm�[�h�œd���E������2�����z��Ɋi�[����D
		}
	}

	// ���f���o�̓t�@�C���|�C���^�̏�����
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

	// ���f�����z�X�g�ɑ��M
	else{
		if(irank != IRANK_MAX){
			MPI_Send(&field_xy[0][0], (xmax)*(ymax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 		// �m�[�h0�ȊO�̃m�[�h���m�[�h0�ɓd���E�����𑗂�D
		}
		if(irank == IRANK_MAX){
			MPI_Send(&field_xy[0][0], (xmax-1)*(ymax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 	// �m�[�h0�ȊO�̃m�[�h���m�[�h0�ɓd���E�����𑗂�D
		}
	}

	// ��M�������f������S���f�����쐬
	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){		// �m�[�h0���m�[�h1���珇�Ƀf�[�^���󂯎��o�͂��Ă����D
			if(node == IRANK_MAX){					// �m�[�hisize-1�̂�1�Z���������ݒ肵�Ă��邽�ߏ������ŕ���
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

		// �t�@�C���|�C���^�����
		fclose(HZ1); 		
	}

}

//�t�@�C���o��
void output_field(char *dir_name_def){

	//double field_xy[XMAX][YMAX]; 	// Hz-field �̃t�@�C���o�� (�ʐ��������̎��E����)

	if(n <= Nmax - Fcut){
		// ����m�F�̂��߂̃t�@�C���o��
		if(n == Ncheck){
			output_field_write (dir_name_def); 
		}

		// ����I�ȃt�@�C���o��
		if(n % Ncutfield == 0){
			output_field_write (dir_name_def); 
		}
	}
	if((n >= Nmax - Fcut) && (n <= Nmax)){

		// ����_�ł̃t�@�C���o��
		if(n % Ncutfield2 == 0){
			output_field_write (dir_name_def); 
		}
	}
}



void mcircle(int x_circ, int y_circ, int z_circ){

	double R; 

	//���a�Z�����̌v�Z
	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//�v�Z�덷��h�����߂Ɍ��グ���Ă��܂�

	rightquartercircle1(x_circ, y_circ, z_circ, R); 
	leftquartercircle1(x_circ-1, y_circ, z_circ, R); 
	rightquartercircle2(x_circ, y_circ-1, z_circ, R); 
	leftquartercircle2(x_circ-1, y_circ-1, z_circ, R); 
}


void halfcircle(int x_circ, int y_circ, int z_circ){

	double R; 

	//���a�Z�����̌v�Z
	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//�v�Z�덷��h�����߂Ɍ��グ���Ă��܂�

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



