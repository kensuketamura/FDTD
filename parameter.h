/*****************************************************************************/
// ����v�Z���̌v�Z�@�̑䐔
/*****************************************************************************/
#if _FDTD
#define NODE 1
#else
#define NODE 1
#endif

// ����v�Z���̌v�Z�@�̃����N
#define IRANK_MAX (NODE-1)		// �ő�l
#define IRANK_MIN 0				// �ŏ��l(���̐����̊��蓖�Ă�ꂽ�v�Z�@�Ɍ��ʂ�����)
#define ISIZE NODE				// �T�C�Y

#define TRUE 1
#define FALSE 0
/*****************************************************************************/
// �v�Z�p�����[�^		�P�ʂ�nm
// (�z��̒�`��v���v���Z�b�T�Ŏg�p���邽�߁C�}�N����)
// �e�l�̓Z���T�C�Y�̐����{�ɂȂ��Ă��邱��!!
/*****************************************************************************/


/****************************** 3��ڊi�q�V�t�g�\�� ******************************/
#define CELL_SIZE 21			// �Z���T�C�Y
#define PITCH 399				// PC �i�q�萔
#define PITCH_SHIFT_MAX 399 //480		// �i�q�萔�ω�PCW��PC�i�q�萔�̍ő�l

#define SLAB_HEIGHT 210		// �X���u��
#define CLAD_HEIGHT1 525	// �㕔�N���b�h����
#define CLAD_HEIGHT2 0		// �����N���b�h����
#define AIR_HEIGHT 0		// ��C�w����

#define RADIUS 120			// PC�̕W���~�E���a
//#define SX3 80				// �`��(X)������3��ڊi�q�V�t�g��(SX2,SX4=0�łȂ��Ǝg���Ȃ��v���P!!main��1380�s��)
//#define SX1 0				// �`��(X)������1��ڊi�q�V�t�g��
//#define SX2 0				// �`��(X)������2��ڊi�q�V�t�g��(SX3,SX4=0�łȂ��Ǝg���Ȃ��v���P!!main��1380�s��)
//#define SX4 0				//�`��(X)������4��ڊi�q�V�t�g��(SX2,SX3=0�łȂ��Ǝg���Ȃ��v���P!!main��1380�s��)
//#define SY 0				// ��(Y)�����̓��g�H�S�̊i�q�V�t�g��(�Ώ̋��E���g�p���Ă���̂Ŏ��ۂɂ͂��̔{)

#define EXCT_LEN 840					// ��U�_ (���f���̍��[����̋���)
#define EXCT_OBSE_LEN 975				// ��U�_����ϑ��ʂ̒��S�܂ł̋���
#define OBSE_WIRE_LEN 2540				// �ϑ��ʂ̒��S����א����g�H�[�܂ł̋���
#define OBSE_INTER (PITCH)				// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
#define WIRE_OUTPUT_LEN (EXCT_LEN + 45)	// �o�ˍא����g�H�̒���(���ł��̐����̓Z���T�C�Y��NODE���̐����{�ɂ��邽�߂Ɏg�p)
#define WIRE_OUTPUT_OFFSET 0			// �o�ˍא����g�H�̃X���u�I�[�̒���
#define WIRE_WID_OFFSET 180				// �א����̃I�t�Z�b�g�ʂ̔���(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define PCW_SiSLAB_TERMINATION_LEN 255	// PCW����CORE�X���u�I�[�̒���
#define PCW_SiSLAB_OFFSET 0				// PCW�c��CORE�X���u�̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define PCW_WIDTH_CHIRP 180				// PCW���̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)

// ����ȉ��͎�����	
//#define NORM_PCW_PER 0				// �ʏ�PCW������
//#define CHIRP_3RD_LS_PER 0			// 3��ڊi�q�V�t�g�ʃ`���[�vLSPCW�������D 2~4��ڂ܂őΉ��D���V�t�g��/������<cellsize���ƃ`���[�v���Ȃ�
//#define CHIRP_2ND_LS_PER 0			//���g�H���`���[�v�ƃV�t�g�ʃ`���[�v�𓯎��ɍs���ۂ̃`���[�vLSPCW�������DPITCH_SHIFT_PER��菬�����Ȃ��ƃ_��?�@2,3��ڂ̂ݑΉ��D���V�t�g��/������<cellsize���ƃ`���[�v���Ȃ�
//#define PITCH_SHIFT_PER 5				// �i�q�萔�ω�PCW�̎�����//�C��=>���g�H���`���[�v�̎�����
//#define LSPCW_SHIFT_DESCRETE FALSE			// �i�q�萔�ω�PCW�̂Ƃ��C�͂��߂���LSPCW�ɂ���ꍇ��FALSE�@���Ȃ��ꍇ��TRUE
//#define PITCH_SHIFT_CHIRP_PER 0		// �i�q�萔�ω��`���[�vPCW�̎�����
//#define LSPCW_PER 20				// LSPCW������
//#define PCW_WID 6				// PCW�̗�
// ����ȏ�͎�����

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// �ʏ�PCW��
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// �`���[�vLSPCW��
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW��
#define WIRE_OUTPUT_OFFSET_PER (WIRE_OUTPUT_OFFSET, CELL_SIZE)	// �o�ˍא����g�H�̃X���u�I�[�̒���

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// ���ˊϑ��ʂ̒��S���W (���f���̍��[����̋���)
#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)	// �o�ˊϑ��ʂ̒��S���W (���f���̉E�[����̋���)

#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)			// ���ˍא��� 
#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)		// �o�ˍא���

/****************************** 3��ڊi�q�V�t�g�\�� ******************************/

/*****************************************************************************/
// �S��͗̈�
/*****************************************************************************/
#define XMAX_ALL 200	// SiO2 1340
#define YMAX_ALL 163	//163 PCW_WID:6/163  PCW_WID:8/209  PCW_WID:10/255
#define ZMAX_ALL 42		// CLAD_HEIGHT1:525/42  CLAD_HEIGHT1:750/57  CLAD_HEIGHT1:990/73


/*****************************************************************************/
// �Z���T�C�Y [m]
/*****************************************************************************/
static const double dblCellSize = CELL_SIZE * 1e-9; 
static const double dx = dblCellSize; 
static const double dy = dblCellSize; 
static const double dz = dblCellSize; 
static const double inv_dx = 1/dx; 
static const double inv_dy = 1/dy; 
static const double inv_dz = 1/dz; 

/*****************************************************************************/
// ���ԃX�e�b�v (�N�[�����g�̈�������Ȃǂɒ���)
/*****************************************************************************/
#if _FDTD

/*-------------------- CELL_SIZE:15nm --------------------*/
static const double dt = 2.8e-17; 			// ���ԃX�e�b�v[s]
static const int Nmax = 100; 				// �ŏI���ԃX�e�b�v
/*-------------------- CELL_SIZE:15nm --------------------*/

static const int Ncut = 50000; 				// ���ԃX�e�b�v����\��������Ԋu
static const int Tcut = 1000; 				// �p���[�̕��ς̎Z�o���J�n���鎞�ԃX�e�b�v  (�ŏI�v�Z�X�e�b�v����̍�)
static const int Fcut = 200; 				// �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v����̍�)

#else
static const double dt = 2.8e-17; 			// ���ԃX�e�b�v[s]
static const int Ncut = 5; 				// ���ԃX�e�b�v����\��������Ԋu
static const int Tcut = 30; 			// �G�l���M�[�̕��ς̎Z�o���J�n���鎞�ԃX�e�b�v
static const int Fcut = 30; 			// �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v����̍�)
static const int Nmax = 1; 			// �ŏI���ԃX�e�b�v
#endif

static const int Ncheck = 10; 					// ����m�F�p�̃t�B�[���h���o�͂��鎞�ԃX�e�b�v
static const int Ncutfield = Ncut; 			// �t�B�[���h���o�͂��鎞�ԃX�e�b�v
static const int Ncutfield2 = 5; 			// �����Ԃł̃t�B�[���h���o�͂��鎞�ԃX�e�b�v�Ԋu

/*****************************************************************************/
// ������[MKSA�n]
/*****************************************************************************/
#define PI 3.141592
#define C0 2.997924e8			// �^�󒆂̌��� [m/s]
#define EPSILON0 8.854e-12		// �^�󒆂̗U�d�� [F/m]
#define MU0 (PI*(4.0e-7))		// �^�󒆂̓����� [N/A^2]


/*****************************************************************************/
// �X���u�Ɖ�͗̈�
/*****************************************************************************/
static const double dblSlabHeig = SLAB_HEIGHT * 1.0e-9; 				// �X���u��
static const double dblCladHeight1 = CLAD_HEIGHT1 * 1.0e-9; 				// �㕔�N���b�h�w��
static const int	intSlabHeigPer = INT_DIV (dblSlabHeig/2.0, dblCellSize); 				//�X���u��(�Ώ̋��E�������g�p���Ă���̂ŁC�X���u����1/2)
static const int	intCladHeight1 = INT_DIV (CLAD_HEIGHT1, CELL_SIZE); 				//�㕔�N���b�h��
static const int	air_hc = 0; 					//�N���b�h�㕔��C�w��
static const int	intSlabCen = air_hc + intCladHeight1 + intSlabHeigPer; 			//�����w���S�Z��(+1��[intSlabHeigPer]����Z�����̂Ƃ��ɒ����l�ɂ��邽��)

/*****************************************************************************/
// �t�H�g�j�b�N�������g�H
/*****************************************************************************/
static const double dblPitchCellComp = dblCellSize/2.0; 			// �i�q�萔�̊ۂ�
static const double dblPitch = PITCH * 1e-9; 					// �~�E�i�q�萔

static const double dblRadius = RADIUS * 1e-9; 					// �~�E���a
static const double dblDiamter = 2.0 * dblRadius; 				// �~�E���a

static const int intPitchX = INT_DIV (PITCH, CELL_SIZE); 		// �i�q�萔�̃Z���T�C�Y(X����)
static const int intPitchY = (INT_DIV((PITCH * sqrt(3.0)/2 + 0.5), CELL_SIZE)); 		// �i�q�萔�̃Z���T�C�Y(Y����)	+0.5�͎l�̌ܓ��̂���
static const int intRadius = INT_DIV (RADIUS, CELL_SIZE); 		// �~�E���a

static const int Row_x = 10;
static const int Row_y = 9;
static const int Wid = 10;


/*****************************************************************************/
// �ޗ��̋��ܗ��ƗU�d��
/*****************************************************************************/
static const double nw = 4.5; //�ʎq��ˑw�̋��ܗ�
static const double nsch1 = 3.265; //SCH1�w�̋��ܗ� 3.265:�o���h�M���b�v�g��1100nm
static const double nsch2 = 3.292; //SCH2�w�̋��ܗ� 3.292:�o���h�M���b�v�g��1150nm
static const double nsch3 = 3.32; //SCH3�w�̋��ܗ� 3.32:�o���h�M���b�v�g��1200nm
static const double np = 3.17; //�x���̋��ܗ�

static const double n_core = 3.43; 		// �R�A�̋��ܗ�(CORE)
static const double n_clad = 1.444; 		// �N���b�h���ܗ�(SiO2)

static const double epsilon0 = EPSILON0; 					// �^��̗U�d��
static const double epsilon1 = EPSILON0 * SQ(n_core); 		// �R�A�̗U�d��
static const double epsilon2 = EPSILON0 * SQ(n_clad); 		// �N���b�h�̗U�d��


/*****************************************************************************/
// ��U�֐�
/*****************************************************************************/
#if _FDTD
#if _EXITATION_FUNC
char *dir_name[] = {"1550"}; 		// ��U�֐��̔g�� [nm]
#else
char *dir_name[] = {"1550"}; 		// ��U�֐��̔g�� [nm]
#endif
#else
char *dir_name[] = {"1550"}; 		// ��U�֐��̔g�� [nm]
#endif

static const double delta_omega = 0.05; 						// ���S���g���ŋK�i�������l�S��
static const int Npeak = 500; 								// �s�[�N�X�e�b�v��

// ��U�_�̍��W
static const int ex_y_st = YMAX_ALL - 10; 		// ���g�H�f�ʎn�Z����(��) ����͋�Ԃ̒��ԃZ�����W���瓱�g�H��(����1/2�l�ɂȂ��Ă���)�������Ă���
static const int ex_y_ed = YMAX_ALL; 					// ���g�H�f�ʏI�Z����(��)
static const int ex_z_st = ZMAX_ALL - intSlabHeigPer; 	// ���g�H�f�ʎn�Z����(�c) ����͋�Ԃ̒��ԃZ�����W���瓱�g�H��(��1/2�l)�������Ă���
static const int ex_z_ed = ZMAX_ALL; 			// ���g�H�f�ʏI�Z����(�c)


/*****************************************************************************/
// �ϑ��_�Ɨ�U�_�̐ݒ� (�v�Z�덷��h�����߂Ɍ��グ���Čv�Z)
/*****************************************************************************/

// �����O�̗�U�_�C�ϑ��� (������̓v���O�������ŏ���)
static const int intExctLen = INT_DIV (EXCT_LEN, CELL_SIZE); //��U�_ �����O�̗�U�_

// �|�C���e�B���O�p���[�̍ő�l�C�ŏ��l�̌v�Z���
static const int intObseInter = INT_DIV (OBSE_INTER, CELL_SIZE); 	// �ϑ����

// ������̗�U�_�C�ϑ���
int intExctPortNum;	// ��U�_������v�Z�@�̔ԍ�

int intExctLenPart;	// ��U�_

// �|�C���e�B���O�p���[�̍ő�l�̎Z�o�p

static const int WGlength = (int)((2.00e-6*1e10)/(dx*1e10)); //��`���g�H��(���ˑ��D�Z�����ɕϊ�)
