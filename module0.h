/*****************************************************************************/
// ��̌v�Z�@�ł̉�͗̈�
/*****************************************************************************/
#define XMAX (INT_DIV(XMAX_ALL, NODE) + 1)
#define YMAX YMAX_ALL
#define ZMAX ZMAX_ALL


/*****************************************************************************/
// ��U�֐�
/*****************************************************************************/
//char *dir_name[] = {"1550", "1580"}; 		// ��U�֐��̔g�� [nm]
//char *dir_name[] = {"1550"}; 		// ��U�֐��̔g�� [nm]
//char *dir_name[] = {"1555", "1560", "1570", "1575"}; 		// ��U�֐��̔g�� [nm]
double lambda; 											// ��U�֐��̔g�� [m] (�v���O�������� dir_name ��double�����đ��)
double omega0; 											// ��U�֐��̊p���g�� [s^-1]

// �K�E�V�A���p���X
double sigma; 											// �L���蕝�����肷��萔
//static const double delta_omega = 0.1; 						// ���S���g���ŋK�i�������l�S��
//static const int Npeak = 3000; 								// �s�[�N�X�e�b�v��


/*****************************************************************************/
// Model�|���o���p�̒萔
/*****************************************************************************/
#define CLAD 0									// �N���b�h
#define CORE 1									// �R�A
#define GaInAsP	2								// GaInAsP
#define AIR_GaInAsP	3							// CLAD/GaInAsP
#define EXITATION 30							// ��U�_//���f���ゾ��31
#define OBSERVATION 20							// �ϑ��_//���f���ゾ��21
#define CIRCLE_REF_INDEX	CLAD				//�֐�mcircle�ŏ������ސ���
#define CIRCLE_REF_INDEX2	CLAD		//�֐�mcircle2�ŏ������ސ���
#define CIRCLE_REF_INDEX3	2				//�ꏊ�m�F�p�h�b�g�̐F�w��


/*****************************************************************************/
// �t�H�g�j�b�N�������f���p�錾
/*****************************************************************************/
//���������p�����[�^
struct PNUM {int X; int Y; }; 					//�~���̒��S���W��^����ϐ�.sankaku�Ŏg�p


/*****************************************************************************/
// �O���[�o���ϐ� (MPI�̒ʐM�ł̓O���[�o���ϐ����g�p���Ȃ��ƃG���[��������͗l)
/*****************************************************************************/

#if _CALCULATION_TYPE == _PROPAGATION_CALCULATION

// �d���E (XMAX+1�͕���v�Z�ł�"�̂肵��"���K�v�Ȃ��߁CYMAX+1, ZMAX+1�͑Ώ̋��E�����ŕK�v�Ȃ���)
double Ex[XMAX+1][YMAX+1][ZMAX+1]; 
double Ey[XMAX+1][YMAX+1][ZMAX+1]; 
double Ez[XMAX+1][YMAX+1][ZMAX+1]; 
double Hx[XMAX+1][YMAX+1][ZMAX+1]; 
double Hy[XMAX+1][YMAX+1][ZMAX+1]; 
double Hz[XMAX+1][YMAX+1][ZMAX+1]; 

// �U�d�� (YMAX+1, ZMAX+1�͑Ώ̋��E�����ŕK�v�Ȃ���)
double ALL_epsilonx[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1]; 
double ALL_epsilony[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1]; 
double ALL_epsilonz[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1]; 
double epsilonx[XMAX][YMAX+1][ZMAX+1]; 
double epsilony[XMAX][YMAX+1][ZMAX+1]; 
double epsilonz[XMAX][YMAX+1][ZMAX+1]; 
int ALL_cell[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];		// ���ۂ�ALL_cell[XMAX_ALL+1][YMAX_ALL][ZMAX_ALL]�ŗǂ����C�v���O�����̓s����U�d�̂̔z��Ɠ��������ɂ���
int cell[XMAX][YMAX+1][ZMAX+1];							// ���ۂ�cell[XMAX][YMAX][ZMAX]�ŗǂ����C�v���O�����̓s����U�d�̂̔z��Ɠ��������ɂ���
double epsilon_xy[XMAX+1][YMAX+1], epsilon_yz[YMAX+1][ZMAX+1], epsilon_zx[XMAX+1][ZMAX+1]; 
double epsilon_zx2[XMAX+1][ZMAX+1]; 
int cell_xy[XMAX][YMAX], cell_yz[YMAX][ZMAX]; 
int cell_zx[ZMAX][XMAX];//��16/1/6


/*�z�����E�����K�p�̂Ƃ��Ɏg���d�E�̔z��*/
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

// ���E���z�o�͗p
double field_xy[XMAX][YMAX]; 	// Hz-field �̃t�@�C���o�� (�ʐ��������̎��E����)
double field_yz[YMAX][ZMAX];	// Hz-field �̃t�@�C���o�� (�ʕ����̎��E����)
double field_zx_Hz[ZMAX][XMAX];	// Hy-field �̃t�@�C���o�� (�ʕ����̎��E����)
double field_zx_Hy[ZMAX][XMAX];	// Hy-field �̃t�@�C���o�� (�ʕ����̎��E����)

int x, y, z; 				// ���U���W
int xmax, ymax, zmax; 			// ��ԕ����̍ő�l
int xmax_all, ymax_all, zmax_all; //�����O�̍ő�l
//int n, Nmax; 				//n: ���ԃX�e�b�v�CNmax: ���ԃX�e�b�v�̍ő�l
int n; 					//���ԃX�e�b�v��
//int icut, jcut, kcut; 		//icut, jcut, kcut: �t�B�[���h�o�͂̂Ƃ��Ɏg��
//int PrintStat; 
//int PrintEnd; 
//int source_i1, source_i2; 
//int source_j1, source_j2; 
//int source_k; 
int x_cen, y_cen, z_cen; 			//x_cen, y_cen, z_cen: ��͋�Ԃ̒��S���U���W
int x_model_cen, y_model_cen; 	//x_model_cen, y_model_cen:���f���̒��S���U���W
//char fname[30]; 
int irank, isize; 

// �d���E�v�Z�Ɏg���萔�̐錾
//double dex, dey, dez; 				//�d�E�̌v�Z�Ɏg�������l	dex: Ex, dey: Ey, dez: Ez
//double cnstEx, cnstEy, cnstEz; 	//�d�E�̌v�Z�Ɏg���萔		
//double dhx, dhy, dhz; 				//���E�̌v�Z�Ɏg�������l	dhx: Hx, dhy: Hy, dhz: Hz
static const double cnstHxyz = dt / MU0; 					//���E�̌v�Z�Ɏg���萔		cnstHxyz: Magnetic field calculation


///*��U�p�萔�̐錾*/
//double omega0; 	//omega0: ��U�֐��̊p���g��
//double sigma; 	//�K�E�V�A���p���X�̂��߂̍L���蕝�����肷��萔
//int Npeak; 		//�K�E�V�A���p���X�̃s�[�N�X�e�b�v��

/*���f���ݒ�Ɏg���萔�̐錾 -- �֐�modeling()�ȊO�ł��g���萔������̂ł����Ő錾����D*/
double b; //�f�B�X�N����
int bk; //�f�B�X�N�����̗��U�l

double w; 		//�x�����̔���
int wij; 		//�x�����̔����̗��U�l
double h; 		//�x���̍���
int hk; 			//�x���̍����̗��U�l
//double np; 		//�x���̋��ܗ�
/*
#define CLAD 0
#define AlGaAs3 1
#define AlGaAs8 2
#define GaAs 3
#define AlAs 4
#define AlOxAs 5

double htslab; 		//�X���u�̌����̔���
int ihtslab; 		//�X���u�̌����̔����̗��U�l
double rhole; 		//�~�E�̔��a
int irhole; 			//�~�E���a�̗��U�l
double hphole; 		//�~�E�z��s�b�`�̔���
int ihphole; 		//�~�E�z��s�b�`�̔����̗��U�l
double hpholev; 		//�~�E�z��s�b�`�̔����́�3�{
int ihpholev; 		//�~�E�z��s�b�`�̔����́�3�{�̗��U�l
*/
//���o�̓p���[�̍ő�l�ƍŏ��l���L�^����ϐ�
double powermax_in; 
double powermin_in; 
double powermax_out; 
double powermin_out; 

//double n_core; 		//�����w���ܗ�
//double n_clad; 		//�N���b�h�w���ܗ�
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
/*�ϑ��_*/
int inputi, inputj, inputk; 
int outputi, outputj, outputk; 
double pyn_in, pyn_out; 

/*WG�`���[�v�\��*/
#define WG_chirp_gradient 1
#define WG_chirp_gradient2 0.5
#define WG_chirp_gradient3 0.08 //BOUNDARYLINE == 9�Ŏg�p
#define WG_chirp_gradient4 0.04 //BOUNDARYLINE == 10�Ŏg�p
#define WG_chirp_gradient5 0.015 //BOUNDARYLINE == 14�Ŏg�p
#define WG_chirp_gradient6 0.020 //BOUNDARYLINE == 15�Ŏg�p
#define WG_chirp_gradient7 0.030 //BOUNDARYLINE == 16�Ŏg�p
#define WG_chirp_gradient8 0.020 //BOUNDARYLINE == 17�Ŏg�p
#define WG_chirp_gradient9 0.012 //BOUNDARYLINE == 18�Ŏg�p
#define WG_chirp_gradient10 0.010 //BOUNDARYLINE == 19�Ŏg�p
#define WG_chirp_gradient11 0.008 //BOUNDARYLINE == 20�Ŏg�p
#define WG_chirp_gradient12 0.012 //BOUNDARYLINE == 22�Ŏg�p
#define WG_chirp_gradient13 0.016 //BOUNDARYLINE == 23�Ŏg�p
#define WG_chirp_gradient14 0.020 //BOUNDARYLINE == 24�Ŏg�p
#define WG_chirp_gradient15 0.024 //BOUNDARYLINE == 25�Ŏg�p 1/7
#define WG_chirp_gradient16 0.004 //BOUNDARYLINE == 26�Ŏg�p 1/7
#define WG_chirp_gradient17 0.002 //BOUNDARYLINE == 27�Ŏg�p 1/7


#define WG_chirp_off_in_x  187
#define WG_chirp_off_in_x2  187-16
#define WG_chirp_off_in_x21  291 // 1/4
#define WG_chirp_off_in_x33  316

#define WG_chirp_off_in_x_2  0 //BOUNDARYLINE == 9�Ŏg�p


#define WG_chirp_off_out_x 800+1 //����Ńe�[�p���s���߂��Ȃ��悤�ɐ���
#define WG_chirp_off_out_x2 800+1+16
#define WG_chirp_off_out_x9 788+200 //BOUNDARYLINE == 14�Ŏg�p
#define WG_chirp_off_out_x18 788+200+42 //BOUNDARYLINE == 18�Ŏg�p
#define WG_chirp_off_out_x19 788+200+63 //BOUNDARYLINE == 19�Ŏg�p
#define WG_chirp_off_out_x20 788+200+84 //BOUNDARYLINE == 20�Ŏg�p
#define WG_chirp_off_out_x21 788+200+84 //BOUNDARYLINE == 21�Ŏg�p

#define WG_chirp_off_out_y 75.5+3+40+16
#define WG_chirp_off_out_y9 75.5+3+40+16+1  //BOUNDARYLINE == 9�Ŏg�p
#define WG_chirp_off_out_y15 75.5+3+40+16+1  //BOUNDARYLINE == 15�Ŏg�p
#define WG_chirp_off_out_y16 75.5+3+40+16  //BOUNDARYLINE == 16�Ŏg�p
#define WG_chirp_off_out_y17 75.5+3+40+16+1  //BOUNDARYLINE == 17�Ŏg�p
#define WG_chirp_off_out_y18 75.5+3+40+16+1 //BOUNDARYLINE == 18�Ŏg�p
#define WG_chirp_off_out_y19 75.5+3+40+16+1  //BOUNDARYLINE == 19�Ŏg�p
#define WG_chirp_off_out_y20 75.5+3+40+16+1-10  //BOUNDARYLINE == 20�Ŏg�p
#define WG_chirp_off_out_y21 75.5+3+40+16+1  //BOUNDARYLINE == 21�Ŏg�p



/*�t�@�C���|�C���^�錾*/
//���f��			//Ez				//Hz
FILE *model_xy; 		FILE *fpfez_xy; 		FILE *fpfhz_xy; 
FILE *model_yz; 		FILE *fpfez_yz; 		FILE *fpfhz_yz; 
FILE *model_xz; 		FILE *fpfez_zx; 		FILE *fpfhz_zx;//�� 
FILE *model_xy2; 	FILE *fpfez2_xy; 	FILE *fpfhz2_xy; 
FILE *allmodel_xy; 	 
FILE *allmodel_yz1, *allmodel_yz4, *allmodel_yz7;
FILE *allmodel_zx1, *allmodel_zx4, *allmodel_zx7; //��16/1/6 
FILE *allmodel_zx; //��16/1/6 

//�U�d��
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

FILE *fpparameter; 	//�v�Z�p�����[�^�ۑ��t�@�C��
FILE *fppoynt_para; 		//�|�C���e�B���O�p���[Ex*Hy, Ey*Hx�e�����ۑ��t�@�C��
FILE *fppoynt1; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt2; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt3; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt4; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt5; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt6; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt1h; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt2h; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt3h; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt4h; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt5h; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppoynt6h; 		//�|�C���e�B���O�p���[�ۑ��t�@�C��
FILE *fppowerHz1; 		//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz2; 		//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz3; 		//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz4; 	//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz5; 	//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz6; 	//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz1h; 	//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz2h; 	//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz3h; 	//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz4h; 	//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz5h; 	//�f�ʑS��Hz�ۑ��p
FILE *fppowerHz6h; 	//�f�ʑS��Hz�ۑ��p
FILE *fpHz1; 		// ���ˊϑ���Hz�ۑ��t�@�C��
FILE *fpHz5; 		// �o�ˊϑ���Hz�ۑ��t�@�C��
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
