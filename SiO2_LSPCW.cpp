/*
3 Dimentional Finite Difference Time Domain Calculation with Periodic Boundary Condition
ver. 5.15
programed by IWAI Takeshi(1999-2002)
original PML programed by Super Doctor SAKAI Atsushi(1996-2002)
imaginary part of PML programed by
IWAI Takeshi(1999-2002)

2001/5/1 ver. 1.0 �v���O�����̗����グ(2 Dimensional version)
2001/5/16 ver. 1.1 �v���O�����C��
2001/5/18 ver. 1.2 �v���O�����C��(�g���x�N�g���̏C��)
2001/5/21 ver. 1.3 �v���O��������
2001/5/23 ver. 1.4 �v���O�����C��(�v�Z�͈͂̏C��)
2001/5/23 ver. 2.0 �v���O��������(���g�H�v�Z����)
2001/5/24 ver. 2.1 �v���O�����C��(�t�B�[���h�f���o�������̒ǉ�)
2001/5/28 ver. 2.1 �v���O�����C��(File Close�̐ݒ�ǉ�)
2001/5/29 ver. 2.2 �v���O�����C��(�g���x�N�g���̏C��)
2001/5/30 ver. 2.3 �v���O�����C��(�g���x�N�g��G-P�̒ǉ�)

2001/6/12 ver. 1.0 �v���O�����̗����グ(3 Dimensional version)
2001/7/9  ver. 1.1 �v���O��������

2001/7/17 ver. 3.0 x����PML�z�����E�����v���O�����̗����グ
2001/8/2  ver. 3.1 �o�O���C��
2001/10/24 ver. 3.2 �R���p�C�����Ɍx�����o��o�O���C�� �ׂ����R�����g�ǉ�
2001/10/24 ver. 3.3 �N���b�h�w�̒ǉ�
*/
/*
arranged by Eiichi Mizuta(2003-)

2003/10/1 ver. 4.0b ���f����DBR�^�t�H�g�j�b�N�������g�H�֕ύX
2003/10/7 ver. 4.1b �d�E�p�����[�^��z��Œ�`
2003/10/8 ver. 4.2b PML�w��אڂ��郂�f���Ɠ��� ���f���̓f���o���C�t�B�[���h�̓f���o���C�d���E�p���X�̓f���o���̕ύX ���Ԃ̓f���o���̒ǉ�
2003/10/9 ver. 4.3b �p���X�s�[�N�ʒu��ύX �z�����E������Ex-PML�͈̔͂̏C�� �d���E�v�Z��condition2�͈̔͂̏C��
2003/10/10 ver. 4.4b ���f���ɉ~�E�̓d�E�p�����[�^��ǉ�
2003/10/15 ver. 4.5b �d���E�v�Z��efield�͈̔͂̏C�� ���˔g���S�g���̕ύX
2003/10/17 ver. 4.6b �d�E�p�����[�^�̒�`�̏C�� ���ܗ����z�̎������E����(y����)�̏C��
2003/10/18 ver. 4.6 PML�w���ψ�}���֖߂�
2003/10/18 ver. 4.6a �v���O�����I�����Ƀt�B�[���h��f���o���悤�ɏC�� ���[�`���I����Ɏ��Ԃ�f���o���悤�ɏC��
2003/10/19 ver. 4.7 �~�E������DBR���ߋ��ܗ��ɋ�C�t�B�����O�t�@�N�^���l�� GaInAsP��InP�̒萔�����ւ�
2003/10/29 ver. 4.71 ��U�p���X����ύX ��U�ʒu�𒆉��݂̂ɕύX
2003/11/13 ver. 4.72 �~�E�̃y�A���������w�㉺�ŕʂɒ�` ���ԓf���o���̏C�� �����w�̌�����1�Z���P�ʂŕς�����悤�ɕύX
2003/11/14 ver. 4.73 �v�Z�X�e�b�v�I���ԍۂɃt�B�[���h��4/�Έʑ��������Ԋu�œf���o�����s���悤�ɕύX
2003/12/14 ver. 4.74 �~�E�t�߂̖�����ς�����@��ύX
2004/2/8 ver. 4.75 ���g�H���𒲐��ł���悤�ɕύX
2004/2/17 ver. 4.76 imax��()�͈͂̕ύX

2004/5/24 ver. 5.00 3�w�\���v���O�����Ƃ��ă��f�����̕ύX
2004/5/25 ver. 5.01 TE���[�h��TM���[�h�̃f�[�^�f���o���𕪗� condition�t�@�C���̓f���o���f�[�^�̓��e��ύX
2004/5/31 ver. 5.02 TE���[�h��TM���[�h�̐����̊ԈႢ���C��
*/

/*
arranged by Daisuke Mori(2005-)

2005/4/6 ver. 5.03 ���f����Si�X���u�G�A�u���b�W�ɕύX ��U�p���X����ύX ��U�𕡐��_��
2005/4/8 ver. 5.10 �d���E�����̓f���o����CSV�`���e�L�X�g�� ���f���ɓ��g�H�����ւ̉~�E�ǉ� mcircle()�ɉ~�E���a�̈����𓱓�
2005/4/11 ver. 5.11 �t�[���G�ϊ��̃t�@�C���ǂݍ��݂�xls����e�L�X�g�ɏC�� ikk���O���ǂݍ��݂���悤�ɕύX
2005/4/22 ver. 5.12 �t�[���G�ϊ��̃t�@�C���ǂݍ��݂����s���Ă����̂��C��
2005/5/16 ver. 5.13 ���~�E�̌@��c�����f�����ɑΉ�
2006/4/5 ver. 5.14 ���g�H�����̉~�E���a�̏o�͂�ǉ�(RADIUS2)
2007/1/15 ver. 5.15 ���g�H�e1��2��̉~�E���a�̏o�͂�ǉ�(RADIUS3,4)
2007/1/15 ver. 5.16 �~�E�̕`���ύX(�O��菬�����Ȃ�)

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "FFT.c"
extern void rdft(int n, int isgn, double *a);

/************************************ �e��萔�̒�` ************************************/
#define PI 3.141592					// �~����
#define C0 2.997924e8				// �^�󒆂̌��� [m/s]
#define EPSILON0 8.854e-12			// �^��U�d�� [F/m]
#define MU0 (PI*(4.0e-7))			// �^�󒆂̓�����
#define Square(x) ((x)*(x))

/****************************** �v���O�����S�̂Ɋւ���錾 ******************************/
#define FDTD 0						// FDTD�v�Z 0:���f����f���o�����[�h
// 1:ON 0:condition.txt��nn.xls�t�@�C����|�������C�v�Z�f�[�^�����݂���ꍇ�̓t�[���G�ϊ����s��
#define PC_WAVEGUIDE 1				// Photonic Crystals WaveGuide�v�Z 0:OFF(�P�ʃZ��) 1:1�񌇑����g�H
//#define unitcell_number 9			// �P�ʃZ���̐�(1:�P�ʃZ�� �:1�񌇑����g�H(x�����P�ʃZ�����))

//�ύX�ꏊ(�n��)
#define unitcell_number 11	//7		// �P�ʃZ���̐�(1:�P�ʃZ�� �:1�񌇑����g�H(x�����P�ʃZ�����))
//�ύX�ꏊ(�I��)

// PC_WAVEGUIDE��unitcell_number�̐ݒ���e�͓��ꂳ���邱��
#define K_STARTPOINT 4				// �g���̌v�Z�����J�n�_�̐ݒ� 1:G-X 2:G-J 3:X-J 4:G-P(���g�H�v�Z�̏ꍇ)
#define K_ENDPOINT 4				// �g���̌v�Z�����I���_�̐ݒ�
// ���̃v���O�����ł�x�����ɋz�����E��ݒ肵�Ă���̂ŁCG-P�̂ݑΉ�
//#define IKK_STARTPOINT 1			// �g���̍ŏ��̒l 0����n�߂Ȃ���������
//#define IKK_ENDPOINT 18				// �g���̍Ō�̒l IKK_MAX�܂Ōv�Z����ƒ[�܂Ōv�Z�������ƂɂȂ�

//�ύX�ꏊv5.11(�n��)
//#define IKK_STARTPOINT 10			// �g���̍ŏ��̒l 0����n�߂Ȃ���������
//#define IKK_ENDPOINT 10				// �g���̍Ō�̒l IKK_MAX�܂Ōv�Z����ƒ[�܂Ōv�Z�������ƂɂȂ�
int IKK_STARTPOINT;				// �g���̍ŏ��̒l 0����n�߂Ȃ���������
int IKK_ENDPOINT;					// �g���̍Ō�̒l IKK_MAX�܂Ōv�Z����ƒ[�܂Ōv�Z�������ƂɂȂ�
//�ύX�ꏊ(�I��)

//#define IKK_MAX 20					// �g���̕����� �o���h�����߂�ۂ̓_���ɑ���

//�ύX�ꏊ(�n��)
#define IKK_MAX 50					// �g���̕����� �o���h�����߂�ۂ̓_���ɑ���
//�ύX�ꏊ(�I��)

#define AIR 1						// �z��nn[][][]�ɂ������C�̒萔
//#define SI 0						// �z��nn[][][]�ɂ�����Si�̒萔

//�ύX�ꏊ(�n��)
#define SIO2 2	// �z��nn[][][]�ɂ�����SiO2�̒萔
//�ύX�ꏊ(�I���)
					
//#define ALL_SPACE_MEDIUM SI			// ��Ԃ�����������}��

#define CLAD_MEDIUM AIR			// �N���b�h�w�}��
#define CIRCLE_MEDIUM AIR			// �������}��

//�ύX�ꏊ(�n��)
#define InP_MEDIUM 2				// �z��nn[][][]�ɂ�����InP�̒萔
#define Act_MEDIUM 3				// �z��nn[][][]�ɂ�����Act�̒萔
//�ύX�ꏊ(�I��)

#define REF_INDEX1 3.43				// Si�̋��ܗ�(�R�A)
//#define REF_INDEX2 1.0					// Air�N���b�h�̋��ܗ�
#define REF_INDEX2 1.444				// SiO2�N���b�h�̋��ܗ�

//�ύX�ꏊ(�n��)
#define InP_INDEX 3.1700			// InP�̋��ܗ�(�g�p���Ă��Ȃ��l)
//#define Act_INDEX 3.534			// �����w�̋��ܗ�******�ύX����*******(Ag-As2Se3)
#define Act_INDEX REF_INDEX1		// �R�A�̋��ܗ�
//�ύX�ꏊ(�I��)

/********************************* PML�v���O�����̐錾 **********************************/
#define L 32						// PML�w��
#define M 2.0						// ���d���̕��z��^���鎟��
#define order 16					// sigma_max�Z�o�Ɏg�����ˌW��[R0=10^(-order)]�̎���

/***************************** �݌v�l�E��U�֐��Ɋւ���錾 *****************************/
//#define PITCH  0.28e-6				// �O�p�z��̃s�b�`
//#define RADIUS 0.08e-6				// �~�E�C�~���̔��a
//#define CORE 0.16e-6				// �R�A�̌���
//#define CLAD 0.80e-6				// �N���b�h�̌���

//�ύX�ꏊ(�n��)
#define PITCH   0.399e-6				// �O�p�z��̃s�b�`
#define RADIUS  0.110e-6				// �~�E�C�~���̔��a1���
#define RADIUSZ 0.110e-6			// �~�E�C�~���̔��a2���
#define RADIUS3 0.110e-6			// �~�E�C�~���̔��a3���
#define RADIUS4 0.110e-6			// �~�E�C�~���̔��a4���
#define RADIUS5 0.110e-6			// �~�E�C�~���̔��a5���
#define RADIUS6 0.110e-6			// �~�E�C�~���̔��a6���
#define RADIUS7 0.110e-6			// �~�E�C�~���̔��a7���
#define RADIUSex1 0.000e-6			// �~�E�C�~���̔��aex1��ځi��{�I��0�ɂ��邱�Ɓj
#define RADIUSex2 0.000e-6			// �~�E�C�~���̔��aex2��ځi��{�I��0�ɂ��邱�Ɓj
#define RADIUSd 0.095e-6			//�قȂ�~�E�̔��a       //�ǉ�take
#define LENGTH_L 0.420e-6			//�ǉ�take


#define RADIUS2 0.000e-6			// �~�E�C�~���̔��a2
//#define RADIUS3 0.14e-6			// �~�E�C�~���̔��a3(���g�H�e)
//#define RADIUS4 0.14e-6			// �~�E�C�~���̔��a4(���g�H�e2���)

#define AIR_clad 0.0e-6				// ��C�w�̌���
#define InP_circle_over 1.5e-6		// �����w�㕔�̉~�E�[��
#define InP_circle_down	1.5e-6		// �����w�����̉~�E�[��
#define InP_sub 0.0e-6				// InP��̌���
#define Act 0.210e-6					// �����w�̌���
#define Act_remain 0.00e-6			// �����w�̌@��c��
#define Width 0						// ���g�H���̒���(��������͂��邱��)

#define FIRST_Sx 0.000e-6				// 1��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define FIRST_Sy 0.000e-6               // 1��ڂ�y�����̃V�t�g��(+�ŋ��߂����)

#define SECOND_Sx 0.000e-6				// 2��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define SECOND_Sy 0.000e-6              // 2��ڂ�y�����̃V�t�g��(+�ŋ��߂����)

#define THIRD_Sx 0.000e-6				// 3��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define THIRD_Sy 0.000e-6               //�@3��ڂ�y�����̃V�t�g��(+�ŋ��߂����)

#define FOURTH_Sx 0.000e-6				// 4��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define FOURTH_Sy 0.000e-6             // 4��ڂ�y�����̃V�t�g��(+�ŋ��߂��

#define FIFTH_Sx 0.000e-6				// 5��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define FIFTH_Sy 0.000e-6               // 5��ڂ�y�����̃V�t�g��(+�ŋ��߂��

#define SIXTH_Sx 0.000e-6				// 6��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define SIXTH_Sy 0.000e-6               // 6��ڂ�y�����̃V�t�g��(+�ŋ��߂��

#define SEVENTH_Sx 0.000e-6		   	    // 7��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define SEVENTH_Sy 0.000e-6             // 7��ڂ�y�����̃V�t�g��(+�ŋ��߂��

#define Ex1_Sx 0.000e-6		   	    // 7��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define Ex1_Sy 0.000e-6             // 7��ڂ�y�����̃V�t�g��(+�ŋ��߂��

#define Ex2_Sx 0.000e-6		   	    // 7��ڂ�x�����̃V�t�g��(�Z���T�C�Y�Ŋ���؂�邱�ƁI�I�j
#define Ex2_Sy 0.000e-6             // 7��ڂ�y�����̃V�t�g��(+�ŋ��߂��


//�ύX�ꏊ(�I��)

//#define NORMALFREQ 0.30				// ���˔g���S�K�i�����g��
//#define LAMBDA (PITCH/NORMALFREQ)	// ���˔g���S�g���D����͕ύX���Ȃ�

//�ύX�ꏊ(�n��)
#define LAMBDA 1.541e-6				// ���˔g���S�g���D����͕ύX���Ȃ�(1.45e-6)
//�ύX�ꏊ(�I��)

//#define CELLLAMBDA 0.700e-6			// �Z���T�C�Y�����߂�g���C�ʏ�͕ύX���Ȃ�
#define excite 1					// ��U�̎�� 1:TE��U 2:TM��U
#define omega0 (2.0*PI*C0/LAMBDA)	// ��U�֐��̊p���g��
//#define sigma  (omega0*0.3/1.5)		// ��U�֐��̃p���X�������߂�萔
//#define PULSEPEAK 2000				// ��U�֐��p���X�̃s�[�N�ƂȂ鎞�� NORMALFREQ��ς����ꍇ�͗�U�p���X���؂�Ȃ��悤�ɒ���

//�ύX�ꏊ(�n��)
#define sigma  (omega0*0.005)		// ��U�֐��̃p���X�������߂�萔(�o���h�v�Z*0.03�C���[�h�f�ʐ�*0.005)
#define PULSEPEAK 22000			// ��U�֐��p���X�̃s�[�N�ƂȂ鎞�� NORMALFREQ��ς����ꍇ�͗�U�p���X���؂�Ȃ��悤�ɒ���(�o���h�v�Z1000�C���[�h�f�ʐ�22000)
//�ύX�ꏊ(�I��)

#define nmax 262144  				// ���ԃX�e�b�v�̍ő�l ��{�I��2�̏搔
// 2^15=32768 2^16=65536 2^17=131072 2^18=262144 2^19=524288 2^20=1048576 2^21=2097152

/********************************** ��{�I�ɂ�����Ȃ� **********************************/
//#define lambda(x) (CELLLAMBDA/(REF_INDEX1*x))	// �g��/(���ܗ�*x)
//#define dx lambda(10)							// ������ԕ����T�C�YX����
//#define dy lambda(10)							// ������ԕ����T�C�YY����
//#define dz lambda(10)							// ������ԕ����T�C�YZ����

//�ύX�ꏊ(�n��)
#define dx 0.021e-6								// ������ԕ����T�C�YX����
#define dy 0.021e-6								// ������ԕ����T�C�YY����
#define dz 0.021e-6								// ������ԕ����T�C�YZ����
//�ύX�ꏊ(�I��)

#define dt (dx*0.0185e-7)						// �������ԕ����T�C�Y(Courant��������ɒ��ӁI)
//#define EPSILON1 (Square(REF_INDEX1)*EPSILON0)	// Si�̗U�d��

#define EPSILON2 (Square(REF_INDEX2)*EPSILON0)				// �N���b�h�̗U�d��

//�ύX�ꏊ(�n��)
#define InP_EPSILON (Square(InP_INDEX)*EPSILON0)			// InP�̗U�d��
#define Act_EPSILON (Square(Act_INDEX)*EPSILON0)			// �����w�̗U�d��
//�ύX�ꏊ(�I��)

#define cnstM (dt/MU0)							// Yee���[�`�����Ŏg���萔(���E)
//#define cnstE (dt/EPSILON0)						// Yee���[�`�����Ŏg���萔(�^�󒆂̓d�E)
//#define cnstE1 (dt/EPSILON1)					// Yee���[�`�����Ŏg���萔(Si�R�A���̓d�E)
#define cnstE2 (dt/EPSILON2)					// Yee���[�`�����Ŏg���萔(SiO2�N���b�h���̓d�E)

//�ύX�ꏊ(�n��)
#define Air_E (dt/EPSILON2)						// Yee���[�`�����Ŏg���^�󒆂̓d�E //10/27��SIO2�ɕύX(cnstE2�Ɠ���)
#define InP_E (dt/InP_EPSILON)					// InP���̓d�E
#define Act_E (dt/Act_EPSILON)					// �����w���̓d�E
//�ύX�ꏊ(�I��)

/******************************** �t�[���G�ϊ��֘A�̐錾 ********************************/
#define STEP nmax+1000				// ���̓t�@�C���̒l�̑���(FDTD�Ōv�Z�����Ƃ��̃X�e�b�v��)
#define F_NUMBER 5					// �t�@�C�����̊ϑ��_�̐�
//#define ELE 6						// �����̐�(Ex,Ey,Ez,Hx,Hy,Hz)�̏��ԂɌv�Z

//�ύX�ꏊ(�n��)
#define ELE 3						// �����̐�(TE��U:Ex,Ey,Hz, TM��U:Ez,Hx,Hy)
//�ύX�ꏊ(�I��)

/*********************************** �Z�����Ȃǂ̐錾 ***********************************/
const int A=(int)((PITCH*1.0e9)/(dx*1.0e9));	// �Z�����ŋK�i������program�p�̃t�H�g�j�b�N�����s�b�`
//#define B 24									// A*root(3)�̒l������(�O�p�i�q�̏ꍇ) �K�������ɂȂ�悤�ɂ���
//const int C=(int)((CORE*1.0e9)/(dz*1.0e9));		// �Z�����ŋK�i������program�p�̃R�A��
//const int D=(int)((CLAD*1.0e9)/(dz*1.0e9));		// �Z�����ŋK�i������program�p�̃N���b�h��

//�ύX�ꏊ(�n��)
#define B 32																// A*root(3)�̒l������(�O�p�i�q�̏ꍇ) �K�������ɂȂ�悤�ɂ���
const int AIR_clad_unit=(int)((AIR_clad*1.0e9)/(dz*1.0e9));					// �Z�����ŋK�i������program�p��Air��
const int InP_circle_over_unit=(int)((InP_circle_over*1.0e9)/(dz*1.0e9));	// �Z�����ŋK�i������program�p�̏㕔�~�E��
const int InP_circle_down_unit=(int)((InP_circle_down*1.0e9)/(dz*1.0e9));	// �Z�����ŋK�i������program�p�̉����~�E��
const int InP_sub_unit=(int)((InP_sub*1.0e9)/(dz*1.0e9));					// �Z�����ŋK�i������program�p��InP���
const int Act_unit=(int)((Act*1.0e9)/(dz*1.0e9));							// �Z�����ŋK�i������program�p�̊����w��
const int Act_remain_unit=(int)((Act_remain*1.0e9)/(dz*1.0e9));							// �Z�����ŋK�i������program�p�̊����w�@��c��
//�ύX�ꏊ(�I��)

//#define imax (B*unitcell_number+L+L)			// x�����̃Z����

//�ύX�ꏊ(�n��)
#define imax (B*unitcell_number+L+L+Width)		// x�����̃Z����
//�ύX�ꏊ(�n��)

const int jmax=A;								// y�����̃Z����
//#define kmax (C*12+L+L)							// z�����̃Z����

//�ύX�ꏊ(�n��)
#define kmax (AIR_clad_unit+InP_circle_over_unit+InP_circle_down_unit+InP_sub_unit+Act_unit+L+L)	//z�����̃Z����
//�ύX�ꏊ(�I��)

#define imove 5									// �ϑ��_���`����p�����[�^(��͋�Ԃ̒��S����̋���)
#define jmove 5									// �ϑ��_���`����p�����[�^(��͋�Ԃ̒��S����̋���)

//�ύX�ꏊ(�n��)
#define kmove 5									// �ϑ��_���`����p�����[�^(�����w�̒��S����̋���)
//�ύX�ꏊ(�I��)

/******************************* �t�B�[���h�f���o���̐ݒ� *******************************/
#define PrintStat 260010								// �d���E���z���o�����ԃX�e�b�v�̎n��		
#define PrintStat_f 100
#define PrintEnd 260210							// �d���E���z���o�����ԃX�e�b�v�̏I���(nmax)		
#define PrintEnd_f 260000
//#define PrintNum 1000000						// �d���E���z���o�����ԃX�e�b�v�̊Ԋu �o���Ȃ��Ƃ���nmax�����傫���l���w��

//�ύX�ꏊ(�n��)
#define PrintNum 10 						// �d���E���z���o�����ԃX�e�b�v�̊Ԋu �o���Ȃ��Ƃ���nmax�����傫���l���w��(�ʏ�20000)
#define PrintNum_f 100
//�ύX�ꏊ(�I��)

/************************************ �e��֐��̐錾 ************************************/
void Init_Space(void);							// ��̓��f����ݒ�
void k_vector_phase(void);						// �g���x�N�g���ƈʑ����̐ݒ�
void exc_function(void);						// ��U�֐��̒�`
void calc_hfield(void);							// ���E�v�Z
void calc_h_periodic_boundary_condition(void);	// ���E�v�Z�p�������E����
void calc_efield(void);							// �d�E�v�Z
void calc_e_periodic_boundary_condition(void);	// �d�E�v�Z�p�������E����
void calc_e_periodic_boundary_condition2(void);	// �d�E�v�Z�p�������E����
void efield_PML(void);							// PML�d�E�v�Z
void hfield_PML(void);							// PML���E�v�Z
void circle_arrange(void);						// �~�E�̔z�u
void mcircle(int ,int ,double);						// make circle function
void halfcircle1(int ,int ,double);					// make half circle function 1
void halfcircle2(int ,int ,double);					// make half circle function 2
void halfcircle3(int ,int ,double);					// make half circle function 3
void halfcircle4(int ,int ,double);					// make half circle function 4
void rightquartercircle1(int ,int ,double);
void leftquartercircle1(int ,int ,double);
void rightquartercircle2(int ,int ,double);
void leftquartercircle2(int ,int ,double);
//******�ύX����*******(��������)
void msiftcircle(int ,int ,double);
void rightsiftquartercircle1(int ,int ,double);	
void rightsiftquartercircle2(int ,int ,double);	
void leftsiftquartercircle1(int , int , double );
void leftsiftquartercircle2(int , int , double );
void mtangle1(int ,int ,double );				//�ύXtake
void mtangle2(int ,int ,double );				//�ύXtake
//******�ύX����*******(�����܂�)
void nnprint(void);								// ���ܗ��v�����g
void file_open_field(void);						// �t�B�[���h�f���o���p�t�@�C���I�[�v��
void field_print(int);							// electric field and magnetic field print
void file_close_field(void);					// �t�B�[���h�f���o���p�t�@�C���N���[�Y
void file_open(void);							// �p���X�f���o���p�t�@�C���I�[�v��
void Poyntingprint(int);						// Ex Ey Ez Hx Hy Hz for "poynting calc" print
void file_close(void);							// �p���X�f���o���p�t�@�C���N���[�Y
void initialize_matrix(void);					// �z��̏�����
void fft_main(void);							// �t�[���G�ϊ��v���O����
void Output_Condition(void);					// �e�p�����[�^�̒l��condition�t�@�C���ɑ|���o��

void file_open_field_f(void);						// �t�B�[���h�f���o���p�t�@�C���I�[�v��
void field_print_f(int);							// electric field and magnetic field print
void file_close_field_f(void);					// �t�B�[���h�f���o���p�t�@�C���N���[�Y


/****************************** �e��p�����[�^�E�z��̐錾 ******************************/
double kkx,kky,kkz,G1x,G1y,G2x,G2y,Gz;							// �g���x�N�g���Ƌt�i�q�x�N�g��
double cos_phx,sin_phx,cos_phy,sin_phy,cos_phxy,sin_phxy;		// �������E�����𓱓�����ۂ̈ʑ������܂܂����ϐ�

//�ύX�ꏊ(�n��)
double cnstE[imax][jmax+1][kmax];									// ��͋�Ԃ̕������Ƃ̓d�E
//�ύX�ꏊ(�I��)

double re_ex[imax][jmax+1][kmax+1],re_ey[imax+1][jmax][kmax+1],re_ez[imax+1][jmax+1][kmax];	// ��͋�Ԃ̓d�E����������
double im_ex[imax][jmax+1][kmax+1],im_ey[imax+1][jmax][kmax+1],im_ez[imax+1][jmax+1][kmax];	// ��͋�Ԃ̓d�E����������
double re_hx[imax+1][jmax+1][kmax],re_hy[imax][jmax+1][kmax],re_hz[imax][jmax+1][kmax+1];	// ��͋�Ԃ̎��E����������
double im_hx[imax+1][jmax+1][kmax],im_hy[imax][jmax+1][kmax],im_hz[imax][jmax+1][kmax+1];	// ��͋�Ԃ̎��E����������
double re_hx_jmax[imax+1][kmax],im_hx_jmax[imax+1][kmax];		// �������E�����p���E����
double re_hz_jmax[imax][kmax+1],im_hz_jmax[imax][kmax+1];		// �������E�����p���E����
static short int nn[imax][jmax+1][kmax];						// ���ܗ����z
//FILE *EX;			FILE *EY;			FILE *EZ;
//FILE *HX;			FILE *HY;			FILE *HZ;

//�ύX�ꏊ(�n��)
FILE *EX_XY;			FILE *EY_XY;			FILE *EZ_XY;
FILE *EX_XZ_center;		FILE *EY_XZ_center;		FILE *EZ_XZ_center;
FILE *EX_XZ_end;		FILE *EY_XZ_end;		FILE *EZ_XZ_end;
FILE *EX_YZ;			FILE *EY_YZ;			FILE *EZ_YZ;
FILE *HX_XY;			FILE *HY_XY;			FILE *HZ_XY;
FILE *HX_XZ_center;		FILE *HY_XZ_center;		FILE *HZ_XZ_center;
FILE *HX_XZ_end;		FILE *HY_XZ_end;		FILE *HZ_XZ_end;
FILE *HX_YZ;			FILE *HY_YZ;			FILE *HZ_YZ;
//�ύX�ꏊ(�I��)

FILE *fPoynt1Ex;	FILE *fPoynt1Ey;	FILE *fPoynt1Ez;
FILE *fPoynt1Hx;	FILE *fPoynt1Hy;	FILE *fPoynt1Hz;
FILE *fp1;
FILE *fp2;
//int n,m;							// n:���ԃX�e�b�v m:�g�������ݒ蔻��q

//�ύX�ꏊ(�n��)
int o=16384,n,m;					// o:���ԃX�e�b�v��f���o���X�e�b�v n:���ԃX�e�b�v m:�g�������ݒ蔻��q
//�ύX�ꏊ(�I��)

int ikk;							// ikk:�g�������q
int ci,cj,ck;						// ��͋�Ԃ̒��S���W
int i,j,k;

//�ύX�ꏊ(�n��)
static int a,b;						// a:�����w��艺�̃Z���� b:SiO2�w��艺�̃Z����
//�ύX�ꏊ(�I��)

int effective_circle_cell;
double effective_2r_a;

char fname[30];

char fikk[30];

/**************************** PML�֘A�p�����[�^�E�z��̐錾 *****************************/
//#define PML_coefficient cnstE1			// PML���Ŏg�p����p�����[�^(��{��cnstE��OK �ꍇ�ɂ��cnstE1�Ǝg�������邱��)

//�ύX�ꏊ(�n��)
#define PML_coefficient Act_E			// PML���Ŏg�p����p�����[�^(��{��cnstE��OK �ꍇ�ɂ��cnstE1�Ǝg�������邱��)
//�ύX�ꏊ(�I��)

double sigma_max;						// �O�ǂł̓��d��
double sigma_x,sigma_y,sigma_z;			// ���d��
double u1,u2;							// PML�v�Z�Ɏg���萔
double term1,term2;						// PML�v�Z�Ɏg���萔
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

/************************************ ���C�����[�`�� ************************************/
int main(){
	printf("dx = %e, dy = %e, dz = %e\n",dx,dy,dz);	// �Z���T�C�Y�̕\��
	printf("dt = %e\n",dt);							// �������ԃX�e�b�v�̕\��
//	printf("���a�̃Z���� = %d\n",(int)((RADIUS*5.0e9)/(dx*5.0e9)));	// ���a�̃Z�����̕\��
//	printf("�s�b�`�̃Z���� = %d\n",A);				// �s�b�`�̃Z�����̕\��
//	printf("�R�A�w�̃Z���� = %d\n",C);				// �R�A�w�̃Z�����̕\��
//	printf("�N���b�h�w�̃Z���� = %d\n",D);			// �N���b�h�w�̃Z�����̕\��

//�ύX�ꏊv5.11(�n��)
	FILE *ikk_file;
	ikk_file = fopen("ikk.txt","r");
	fscanf(ikk_file,"%d",&IKK_STARTPOINT);
	IKK_ENDPOINT = IKK_STARTPOINT;
//�ύX�ꏊ(�I��)

//�ύX�ꏊ(�n��)
	printf("�~�E�̔��a�̃Z���� = %d",(int)((RADIUS*5.0e9)/(dx*5.0e9)));					// �~�E�̔��a�̃Z�����̕\��
	printf("		�~�E�̔��a�̃Z����2 = %d\n",(int)((RADIUS2*5.0e9)/(dx*5.0e9)));		// �~�E�̔��a�̃Z�����̕\��2
	printf("�s�b�`�̃Z���� = %d\n",A);													// �s�b�`�̃Z�����̕\��
	printf("���g�H���ω��̃Z���� = %d\n",Width);	
	printf("1��ڃV�t�g�̃Z���� = %d\n",int(FIRST_Sx/dx));								// �s�b�`�̃Z�����̕\��
	printf("1��ڃV�t�g�̃Z���� = %d\n",int(FIRST_Sy/dy));
	printf("2��ڃV�t�g�̃Z���� = %d\n",int(SECOND_Sx/dx));
	printf("2��ڃV�t�g�̃Z���� = %d\n",int(SECOND_Sy/dy));
	printf("3��ڃV�t�g�̃Z���� = %d\n",int(THIRD_Sx/dx));
	printf("3��ڃV�t�g�̃Z���� = %d\n",int(THIRD_Sy/dy));
	printf("4��ڃV�t�g�̃Z���� = %d\n",int(FOURTH_Sx/dx));
	printf("4��ڃV�t�g�̃Z���� = %d\n",int(FOURTH_Sy/dy));
	printf("�����w�㕔�̉~�E�[�� = %d",InP_circle_over_unit);							// �����w�㕔�̉~�E�[���̃Z�����\��
	printf("	�����w�����̉~�E�[�� = %d\n",InP_circle_down_unit);						// �����w�����̉~�E�[���̃Z�����\��
	printf("InP��̃Z���� = %d",InP_sub_unit);										// InP��̃Z�����̕\��
	printf("		�����w�̃Z���� = %d",Act_unit);										// �����w�̃Z�����̕\��
	printf("	��C�w�̃Z���� = %d\n",AIR_clad_unit);									// ��C�w�̃Z�����̕\��
	
	if(excite == 1){
		printf("TE-like mode\n");
	}
	if(excite == 2){
		printf("TM-like mode\n");
	}
//�ύX�ꏊ(�I��)

	ci=imax/2;								// ��͋�Ԃ̒��S��ݒ�
	cj=jmax/2;								// ��͋�Ԃ̒��S��ݒ�
//	ck=kmax/2;								// ��͋�Ԃ̒��S��ݒ�
	printf("imax = %d	jmax = %d	kmax = %d\n",imax,jmax,kmax);		// ��͋�ԍ��W�̕\��
//	printf("ci = %d, cj = %d, ck = %d\n",ci,cj,ck);					// ��͋�Ԃ̒��S���W�̕\��

	sigma_max=-((M+1)*EPSILON0*C0*log(1.0*pow(double(10),double(-order))))/(2*L*dx);	// PML�p�����[�^(sigma_max)�̌v�Z
//	printf("Number of PML : %d, PML constant M : %.2e, Reflection order : %d, sigma_max : %e\n",L,M,order,sigma_max);	// PML�p�����[�^�̕\��

	Init_Space();							// ��̓��f����ݒ�
	Output_Condition();						// �e��p�����[�^��condition.txt�t�@�C���ɏo��

#if FDTD==1
	printf("FDTD simulation START\n");

	for(m=K_STARTPOINT;m<=K_ENDPOINT;m++){
		for(ikk=IKK_STARTPOINT;ikk<=IKK_ENDPOINT;ikk++){
			file_open();
			k_vector_phase();
			initialize_matrix();
			printf("m = %d/%d, ikk = %d/%d\n",m,K_ENDPOINT,ikk,IKK_ENDPOINT);

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
				if(n>=PrintStat && n<=PrintEnd && (n%PrintNum)==0 || n==nmax){
//�ύX�ꏊ(�I��)

					file_open_field();
					field_print(n);
					file_close_field();
				}

				if(n>=PrintStat_f && n<=PrintEnd_f && (n%PrintNum_f)==0 || n==nmax){
//�ύX�ꏊ(�I��)

					file_open_field_f();
					field_print_f(n);
					file_close_field_f();
				}




//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

			}
			file_close();

//�ύX�ꏊ(�n��)
			fclose(tf);
//�ύX�ꏊ(�I��)

		}
	}
#endif

	Poyntingprint(0);
	nnprint();
	fft_main();

	// FFT�t�@�C���̈ړ�
	char OutputFile1[20], OutputFile2[20];
	m = 4;	// �������ł��Ȃ̂ł���ȍ~�Ƀv���O������ǉ�����Ƃ��ɂ͒��ӂ��K�v
	sprintf(OutputFile1, "0%d_FFT%s.xls", m, fikk);
	sprintf(OutputFile2, "..\\%s", OutputFile1);
	if( rename( OutputFile1, OutputFile2 ) != 0 )
	
	return 0;
}
/********************************* ���C�����[�`���I��� *********************************/


/************************************* ���f���̒�` *************************************/
void Init_Space(){
//	int i,j,k;

//�ύX�ꏊ(�n��)
	int i,j,k;
//�ύX�ꏊ(�I��)

	for(i=0;i<imax;i++){
		for(j=0;j<=jmax;j++){
//			for(k=0;k<kmax;k++){
//				nn[i][j][k]=ALL_SPACE_MEDIUM;	// ��͋�Ԃ�������
//			}
//			for(k=0;k<(ck-(C/2)-D);k++){
//				nn[i][j][k]=CIRCLE_MEDIUM;		// �N���b�h�w��������C�w�Ƃ���
//			}
//			for(k=(ck-(C/2)-D);k<(ck-(C/2));k++){
//				nn[i][j][k]=CLAD_MEDIUM;		// �N���b�h�w���`
//			}
//			for(k=(ck+(C/2));k<kmax;k++){
//				nn[i][j][k]=CIRCLE_MEDIUM;		// �R�A�w�㕔����C�w�Ƃ���
//			}

//�ύX�ꏊ(�n��)
			for(k=0;k<kmax;k++){
				nn[i][j][k]=Act_MEDIUM;							// ��͋�Ԃ̏�����(�S�p�����[�^�������w��)
				cnstE[i][j][k]=Act_E;
			}
			for(k=L;k<(L+InP_sub_unit);k++){					// InP_sub�̒�`
				nn[i][j][k]=InP_MEDIUM;
				cnstE[i][j][k]=InP_E;
			}
//�ύX�ꏊ(�n��)
			a=InP_sub_unit+L;
			for(k=a;k<(a+InP_circle_down_unit);k++){			// InP_circle_down�̒�`
				nn[i][j][k]=CLAD_MEDIUM;
				cnstE[i][j][k]=cnstE2;
			}
			b=InP_circle_down_unit+a;
			for(k=b;k<(b+Act_unit);k++){						// �����w���`
				nn[i][j][k]=Act_MEDIUM;
				cnstE[i][j][k]=Act_E;
			}
			b=Act_unit+b;
			for(k=b;k<(b+InP_circle_over_unit);k++){			// InP_circle_over�̒�`
				nn[i][j][k]=CLAD_MEDIUM;
				cnstE[i][j][k]=cnstE2;
//�ύX�ꏊ(�I���)
			}
			b=InP_circle_over_unit+b;
			for(k=b;k<(b+AIR_clad_unit);k++){					// AIR_clad�̒�`
				nn[i][j][k]=AIR;
				cnstE[i][j][k]=Air_E;
			}
//�ύX�ꏊ(�I��)

		}
	}



//�ύX�ꏊ(�n��)
	ck=a+InP_circle_down_unit+Act_unit/2;						//z�����̉�͋�Ԃ̒��S�������w���S�ɐݒ�
	printf("ci = %d	cj = %d		ck = %d\n",ci,cj,ck);			// ��͋�Ԃ̒��S���W�̕\��

	printf("Number of PML : %d		PML constant M : %.2e\nReflection order : %d		sigma_max : %e\n",L,M,order,sigma_max);	// PML�p�����[�^�̕\��
//�ύX�ꏊ(�I��)

	printf("triangular array ON\n");
	circle_arrange();						// ��̓��f���ɉ~�E��z�u

	for(i=0;i<imax;i++){
		for(k=0;k<kmax;k++){
			nn[i][jmax][k]=nn[i][0][k];		// ���ܗ����z�̎������E����(y����)

//�ύX�ꏊ(�n��)
			cnstE[i][jmax][k]=cnstE[i][0][k];	// ���ܗ����z�̎������E����(y����)
//�ύX�ꏊ(�I��)

		}
	}
}
/*
���̋��ܗ����z�ɂ��������E�����𓱓�����̂��Ƃ����ƁC
����𓱓����Ȃ��Ɠd���E�v�Z�̎������E�ʂŎ��͋��ܗ��萔�����������Ȓl���������Ă��܂��D
���C�ɂ��Ȃ�d�v�D
*/

void circle_arrange(){
	int place;
	for(place=0;place<=(unitcell_number*2);place++){
//******�ύX����*******(��������)
//�V�t�g��3���(=����)�̂Ƃ��̉~�E����
		if(place == unitcell_number-3){
			msiftcircle(((B*place)/2)+L+THIRD_Sy/dx,A/2+THIRD_Sx/dy,RADIUS3);
			msiftcircle(((B*place)/2)+L+THIRD_Sy/dx,A*3/2+THIRD_Sx/dy,RADIUSd);                  //�ύXtake
		}
		else if(place == unitcell_number+3){
			msiftcircle(((B*place)/2)+L-THIRD_Sy/dx+Width,A/2+THIRD_Sx/dy,RADIUS3);
			msiftcircle(((B*place)/2)+L-THIRD_Sy/dx,A*3/2+THIRD_Sx/dy,RADIUSd);                  //�ύXtake
		}
//�V�t�g��1���(=����)�̂Ƃ��̉~�E����
		else if(place == unitcell_number-1){
			msiftcircle(((B*place)/2)+L+FIRST_Sy/dx,A/2+FIRST_Sx/dy,RADIUS);
			msiftcircle(((B*place)/2)+L+FIRST_Sy/dx,A*3/2+FIRST_Sx/dy,RADIUSd);                  //�ύXtake
		}
		else if(place == unitcell_number+1){
			msiftcircle(((B*place)/2)+L-FIRST_Sy/dx+Width,A/2+FIRST_Sx/dy,RADIUS);
			msiftcircle(((B*place)/2)+L-FIRST_Sy/dx+Width,A*3/2+FIRST_Sx/dy,RADIUSd);                  //�ύXtake
		}
//�V�t�g��5���(=����)�̂Ƃ��̉~�E����
		else if(place == unitcell_number-5){
			msiftcircle(((B*place)/2)+L+FIFTH_Sy/dx,A/2+FIFTH_Sx/dy,RADIUS5);
			msiftcircle(((B*place)/2)+L+FIFTH_Sy/dx,3*A/2+FIFTH_Sx/dy,RADIUSd);                  //�ύXtake
		}
		else if(place == unitcell_number+5){
			msiftcircle(((B*place)/2)+L-FIFTH_Sy/dx+Width,A/2+FIFTH_Sx/dy,RADIUS5);
			msiftcircle(((B*place)/2)+L-FIFTH_Sy/dx+Width,A*3/2+FIFTH_Sx/dy,RADIUSd);                  //�ύXtake
		}/*
		else if(place == unitcell_number+10){
			mtangle1(((B*place)/2)+L-SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,LENGTH_L);
//			mtangle1(((B*place)/2)+L-SEVENTH_Sy/dx,A*3/2+SEVENTH_Sx/dy,LENGTH_L);
		}
		else if(place == unitcell_number-10){
			mtangle2(((B*place)/2)+L+SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,LENGTH_L);
//			mtangle2(((B*place)/2)+L+SEVENTH_Sy/dx,3*A/2+SEVENTH_Sx/dy,LENGTH_L);
		}*/
//�V�t�g��7���(=����)�̂Ƃ��̉~�E����
		else if(place == unitcell_number-7){
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //�ύXtake
		}
		else if(place == unitcell_number+7){
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //�ύXtake
		}
//�V�t�g��9���(=����)�̂Ƃ��̉~�E����
		else if(place == unitcell_number-9){
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //�ύXtake
		}
		else if(place == unitcell_number+9){
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //�ύXtake
		}
//�V�t�g��11���(=����)�̂Ƃ��̉~�E����
		else if(place == unitcell_number-11){
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L+SEVENTH_Sy/dx,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //�ύXtake
		}
		else if(place == unitcell_number+11){
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A/2+SEVENTH_Sx/dy,RADIUS7);
			msiftcircle(((B*place)/2)+L-SEVENTH_Sy/dx+Width,A*3/2+SEVENTH_Sx/dy,RADIUSd);                  //�ύXtake
		}
//�V�t�g��ex1,ex2��ڂ̂Ƃ��̉~�E����
		else if(place == unitcell_number-4){
			msiftcircle(((B*place)/2)+L+Ex1_Sy/dx,A/2+Ex1_Sx/dy,RADIUSex1);
			msiftcircle(((B*place)/2)+L+Ex2_Sy/dx,A/2+Ex2_Sx/dy,RADIUSex2);
		}
		else if(place == unitcell_number+4){
			msiftcircle(((B*place)/2)+L-Ex1_Sy/dx+Width,A/2+Ex1_Sx/dy,RADIUSex1);
			msiftcircle(((B*place)/2)+L-Ex2_Sy/dx+Width,A/2+Ex2_Sx/dy,RADIUSex2);
		}




#if PC_WAVEGUIDE==0		// �P�ʃZ��
		if(place%2==1){
			halfcircle1(B/2,0,RADIUS);
			halfcircle2(B/2,A,RADIUS);
#endif
#if PC_WAVEGUIDE==1		// 1�񌇑����g�H
		if(place%2==1 && place!=unitcell_number){
//			halfcircle1(((B*place)/2)+L,0,RADIUS);
//			halfcircle2(((B*place)/2)+L,A,RADIUS);


//�ύX�ꏊ(�n��)
//******�ύX����*******(��������)
//��񂩂�1,3��ڂłȂ��Ƃ�
				if(place == unitcell_number-2){
					msiftcircle(((B*place)/2)+L+SECOND_Sy/dx,A+SECOND_Sx/dy,RADIUSZ);
					msiftcircle(((B*place)/2)+L+SECOND_Sy/dx,A*2+SECOND_Sx/dy,RADIUSd);                  //�ύXtake
				}
				else if(place == unitcell_number+2){
					msiftcircle(((B*place)/2)+L-SECOND_Sy/dx+Width,A+SECOND_Sx/dy,RADIUSZ);
					msiftcircle(((B*place)/2)+L-SECOND_Sy/dx+Width,2*A+SECOND_Sx/dy,RADIUSd);                  //�ύXtake
				}
				else if(place == unitcell_number-4){
					msiftcircle(((B*place)/2)+L+FOURTH_Sy/dx,A+FOURTH_Sx/dy,RADIUS4);
					msiftcircle(((B*place)/2)+L+FOURTH_Sy/dx,2*A+FOURTH_Sx/dy,RADIUSd);                  //�ύXtake
				}
				else if(place == unitcell_number+4){
					msiftcircle(((B*place)/2)+L-FOURTH_Sy/dx+Width,A+FOURTH_Sx/dy,RADIUS4);
					msiftcircle(((B*place)/2)+L-FOURTH_Sy/dx+Width,2*A+FOURTH_Sx/dy,RADIUSd);                  //�ύXtake
				}
				else if(place == unitcell_number-6){
					msiftcircle(((B*place)/2)+L+SIXTH_Sy/dx,A+SIXTH_Sx/dy,RADIUS6);
					msiftcircle(((B*place)/2)+L+SIXTH_Sy/dx,2*A+SIXTH_Sx/dy,RADIUSd);                  //�ύXtake
				}
				else if(place == unitcell_number+6){
					msiftcircle(((B*place)/2)+L-SIXTH_Sy/dx+Width,A+SIXTH_Sx/dy,RADIUS6);
					msiftcircle(((B*place)/2)+L-SIXTH_Sy/dx+Width,2*A+SIXTH_Sx/dy,RADIUSd);                  //�ύXtake
				}
				else if(place == unitcell_number-8){
					msiftcircle(((B*place)/2)+L+SIXTH_Sy/dx,A+SIXTH_Sx/dy,RADIUS6);
					msiftcircle(((B*place)/2)+L+SIXTH_Sy/dx,2*A+SIXTH_Sx/dy,RADIUSd);                  //�ύXtake
				}
				else if(place == unitcell_number+8){
					msiftcircle(((B*place)/2)+L-SIXTH_Sy/dx+Width,A+SIXTH_Sx/dy,RADIUS6);
					msiftcircle(((B*place)/2)+L-SIXTH_Sy/dx+Width,2*A+SIXTH_Sx/dy,RADIUSd);                  //�ύXtake
				}
				else{
					if(place<unitcell_number){
						halfcircle1(((B*place)/2)+L,0,RADIUSd);
						mcircle(((B*place)/2)+L,A,RADIUS);                 //�ύXtake
						halfcircle2(((B*place)/2)+L,2*A,RADIUSd);                  //�ύXtake
					}
					if(place>unitcell_number){
						halfcircle1(((B*place)/2)+L+Width,0,RADIUSd);
						mcircle(((B*place)/2)+L+Width,A,RADIUS6);                  //�ύXtake
						halfcircle2(((B*place)/2)+L+Width,2*A,RADIUSd);                  //�ύXtake
					}
				}
//******�ύX����*******(�����܂�)

				
//******�ύX����*******(��������)
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
//******�ύX����*******(�����܂�)
//�ύX�ꏊ(�I��)

#endif
		}

//�ύX�ꏊv5.10(�n��)
#if PC_WAVEGUIDE==1		// 1�񌇑����g�H
		if(place==unitcell_number){
			mcircle(((B*place)/2)+L+Width/2,A/2,RADIUS2);		//���g�H�����ʑ��V�t�g�~�E
		}
#endif
//�ύX�ꏊ(�I��)
	
	}
	printf("\n");
}
/********************************** ���f���̒�`�I��� **********************************/


/*********************************** �~�E�̒�`�Ɛ��� ***********************************/
void mcircle(int x, int y, double r){
	rightquartercircle1(x,y,r);
	leftquartercircle1(x-1,y,r);
	rightquartercircle2(x,y-1,r);
	leftquartercircle2(x-1,y-1,r);
	printf("[%4d %4d]",x,y);
}
//******�ύX����*******(��������)
//�V�t�g�����~�E���͂ݏo�����Ƃ��p
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

//�ύX�ꏊ(�n��)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//�ύX�ꏊ(�I��)

//�ύX�ꏊv5.13(�n��)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//�ύX�ꏊ(�I��)

				r=sqrt(double((i-x+1)*(i-x+1)+(j-y+1)*(j-y+1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//�ύX�ꏊ(�n��)
				if(r<=R){	//5.16
					if(j>2*A-1){                                                           //�ύXtake
						nn[i][j-2*A][k]=CIRCLE_MEDIUM;
						cnstE[i][j-2*A][k]=Air_E;
					}else{
						nn[i][j][k]=CIRCLE_MEDIUM;
						cnstE[i][j][k]=Air_E;
					}
					if (k==a) circle_cell++;
				}
//�ύX�ꏊ(�I��)

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
					if(j>2*A-1){                                                           //�ύXtake
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

//�ύX�ꏊ(�n��)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//�ύX�ꏊ(�I��)

//�ύX�ꏊv5.13(�n��)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//�ύX�ꏊ(�I��)

				r = sqrt(double((i-x-1)*(i-x-1)+(j-y+1)*(j-y+1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//�ύX�ꏊ(�n��)
				if(r<=R){	//5.16
					if(j>2*A-1){                                                           //�ύXtake
						nn[i][j-2*A][k]=CIRCLE_MEDIUM;
						cnstE[i][j-2*A][k]=Air_E;
					}else{
						nn[i][j][k]=CIRCLE_MEDIUM;
						cnstE[i][j][k]=Air_E;
					}
				}
//�ύX�ꏊ(�I��)

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
					if(j>2*A-1){                                                           //�ύXtake
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
				if(j>2*A-1){                                                           //�ύXtake
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
				if(j>2*A-1){                                                           //�ύXtake
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




//******�ύX����*******(�����܂�)

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

//�ύX�ꏊ(�n��)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//�ύX�ꏊ(�I��)

//�ύX�ꏊv5.13(�n��)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//�ύX�ꏊ(�I��)

				r=sqrt(double((i-x+1)*(i-x+1)+(j-y+1)*(j-y+1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//�ύX�ꏊ(�n��)
				if(r<=R){	//5.16
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
					if (k==a) circle_cell++;
				}
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//�ύX�ꏊ(�I��)

//�ύX�ꏊv5.13(�n��)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//�ύX�ꏊ(�I��)

				r = sqrt(double((i-x-1)*(i-x-1)+(j-y+1)*(j-y+1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//�ύX�ꏊ(�n��)
				if(r<=R){	//5.16
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
				}
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//�ύX�ꏊ(�I��)

//�ύX�ꏊv5.13(�n��)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//�ύX�ꏊ(�I��)

				r = sqrt(double((i-x+1)*(i-x+1)+(j-y-1)*(j-y-1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//�ύX�ꏊ(�n��)
				if(r<=R){	//5.16
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
				}
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
			for(k=a;k<(a+InP_circle_down_unit+Act_unit+InP_circle_over_unit);k++){
//�ύX�ꏊ(�I��)

//�ύX�ꏊv5.13(�n��)
				if (RAD ==RADIUS2 && k==(a+Act_unit+InP_circle_over_unit-Act_remain_unit))
					break;
//�ύX�ꏊ(�I��)

				r = sqrt(double((i-x-1)*(i-x-1)+(j-y-1)*(j-y-1)))-0.5;	//5.16
//				r=(int)sqrt((i-x)*(i-x)+(j-y)*(j-y));
//				if(r<R) nn[i][j][k]=CIRCLE_MEDIUM;

//�ύX�ꏊ(�n��)
				if(r<=R){	//5.16
					nn[i][j][k]=CIRCLE_MEDIUM;
					cnstE[i][j][k]=Air_E;
				}
//�ύX�ꏊ(�I��)

			}
		}
	}
}

/*
���̂��������~�E��4�������Ē�`���Ă��邩�Ƃ����ƁC
FDTD�v�Z�ł͂������Ȃ��Ɛ^�~����邱�Ƃ��ł��Ȃ�����D
*/
/******************************** �~�E�̒�`�Ɛ����I��� ********************************/


/************************************ �ʑ����̌v�Z ************************************/
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
	if(m==4){	// G-P(���g�H�v�Z�̏ꍇ)
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
/********************************* �ʑ����̌v�Z�I��� *********************************/


/*************************************** ��U�֐� ***************************************/
void exc_function(){
	double exc_cos,exc_sin;
	exc_cos=cos(omega0*dt*(n-PULSEPEAK))*exp(-((sigma*dt*(n-PULSEPEAK))*(sigma*dt*(n-PULSEPEAK))/2));
	exc_sin=sin(omega0*dt*(n-PULSEPEAK))*exp(-((sigma*dt*(n-PULSEPEAK))*(sigma*dt*(n-PULSEPEAK))/2));

	if(excite==1){
		i=5;j=6;k=3;

//�ύX�ꏊ(�n��)
//		i=0;j=0;k=0;
//�ύX�ꏊ(�I��)

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
��U�͒��S���W(ci,cj,ck)����(i,j,k)���炵���ʒu�ŗ^���Ă���D
�܂�C���f���̐ݒ�͈͂𒴂��Ȃ��悤�ɒ��ӂ��邱�ƁD
*/
/************************************ ��U�֐��I��� ************************************/


/****************************** �d���E�̌v�Z���������E���� ******************************/
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

//�ύX�ꏊ(�n��)
			for(k=L;k<=kmax-L;k++){
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
				re_ex[i][j][k]=re_ex[i][j][k]+cnstE[i][j][k]*re_dex;
				im_ex[i][j][k]=im_ex[i][j][k]+cnstE[i][j][k]*im_dex;
//�ύX�ꏊ(�I��)

			}
		}
	}

//	for(i=L;i<imax-L;i++){

//�ύX�ꏊ(�n��)
	for(i=L;i<=imax-L;i++){
//�ύX�ꏊ(�I��)

		for(j=0;j<jmax;j++){
//			for(k=L;k<kmax-L;k++){

//�ύX�ꏊ(�n��)
			for(k=L;k<=kmax-L;k++){
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
				re_ey[i][j][k]=re_ey[i][j][k]+cnstE[i][j][k]*re_dey;
				im_ey[i][j][k]=im_ey[i][j][k]+cnstE[i][j][k]*im_dey;
//�ύX�ꏊ(�I��)

			}
		}
	}

//	for(i=L;i<imax-L;i++){

//�ύX�ꏊ(�n��)
	for(i=L;i<=imax-L;i++){
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
				re_ez[i][j][k]=re_ez[i][j][k]+cnstE[i][j][k]*re_dez;
				im_ez[i][j][k]=im_ez[i][j][k]+cnstE[i][j][k]*im_dez;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
		for(k=L;k<=kmax-L;k++){
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
			re_ex[i][jmax][k]=re_ex[i][jmax][k]+cnstE[i][jmax][k]*re_dex;
			im_ex[i][jmax][k]=im_ex[i][jmax][k]+cnstE[i][jmax][k]*im_dex;
//�ύX�ꏊ(�I��)

		}
	}

//	for(i=L;i<imax-L;i++){

//�ύX�ꏊ(�n��)
	for(i=L;i<=imax-L;i++){
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
			re_ez[i][jmax][k]=re_ez[i][jmax][k]+cnstE[i][jmax][k]*re_dez;
			im_ez[i][jmax][k]=im_ez[i][jmax][k]+cnstE[i][jmax][k]*im_dez;
//�ύX�ꏊ(�I��)

		}
	}
}
/*************************** �d���E�̌v�Z���������E�����I��� ***************************/


/************************************* �z�����E���� *************************************/
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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz010[i][j][k]=term1*re_exz010[i][j][k]-term2*(re_hy[i][j][k]-re_hy[i][j][k-1])/dz;
				im_exz010[i][j][k]=term1*im_exz010[i][j][k]-term2*(im_hy[i][j][k]-im_hy[i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz210[i][j][k]=term1*re_exz210[i][j][k]-term2*(re_hy[imax-1-i][j][k]-re_hy[imax-1-i][j][k-1])/dz;
				im_exz210[i][j][k]=term1*im_exz210[i][j][k]-term2*(im_hy[imax-1-i][j][k]-im_hy[imax-1-i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz012[i][j][k]=term1*re_exz012[i][j][k]-term2*(re_hy[i][j][kmax-k]-re_hy[i][j][kmax-k-1])/dz;
				im_exz012[i][j][k]=term1*im_exz012[i][j][k]-term2*(im_hy[i][j][kmax-k]-im_hy[i][j][kmax-k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz212[i][j][k]=term1*re_exz212[i][j][k]-term2*(re_hy[imax-i-1][j][kmax-k]-re_hy[imax-i-1][j][kmax-k-1])/dz;
				im_exz212[i][j][k]=term1*im_exz212[i][j][k]-term2*(im_hy[imax-i-1][j][kmax-k]-im_hy[imax-i-1][j][kmax-k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
		for(j=1;j<=jmax;j++){
//�ύX�ꏊ(�I��)

			for(k=L;k<=kmax-L;k++){

				sigma_x=sigma_max*pow((double)(L-i)/L,M);
				sigma_y=0.0;
				sigma_z=0.0;

			/*(011)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz011[i][j][k]=term1*re_exz011[i][j][k]-term2*(re_hy[i][j][k]-re_hy[i][j][k-1])/dz;
				im_exz011[i][j][k]=term1*im_exz011[i][j][k]-term2*(im_hy[i][j][k]-im_hy[i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz211[i][j][k]=term1*re_exz211[i][j][k]-term2*(re_hy[imax-1-i][j][k]-re_hy[imax-1-i][j][k-1])/dz;
				im_exz211[i][j][k]=term1*im_exz211[i][j][k]-term2*(im_hy[imax-1-i][j][k]-im_hy[imax-1-i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
		for(j=1;j<=jmax;j++){
//�ύX�ꏊ(�I��)

			for(k=1;k<L;k++){

				sigma_x=0.0;
				sigma_y=0.0;
				sigma_z=sigma_max*pow((double)(L-k)/L,M);

			/*(110)*/
				//Exz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz110[i][j][k]=term1*re_exz110[i][j][k]-term2*(re_hy[i][j][k]-re_hy[i][j][k-1])/dz;
				im_exz110[i][j][k]=term1*im_exz110[i][j][k]-term2*(im_hy[i][j][k]-im_hy[i][j][k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_exz112[i][j][k]=term1*re_exz112[i][j][k]-term2*(re_hy[i][j][kmax-k]-re_hy[i][j][kmax-k-1])/dz;
				im_exz112[i][j][k]=term1*im_exz112[i][j][k]-term2*(im_hy[i][j][kmax-k]-im_hy[i][j][kmax-k-1])/dz;

				//Exy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx010[i][j][k]=term1*re_eyx010[i][j][k]-term2*(re_hz[i][j][k]-re_hz[i-1][j][k])/dx;
				im_eyx010[i][j][k]=term1*im_eyx010[i][j][k]-term2*(im_hz[i][j][k]-im_hz[i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx210[i][j][k]=term1*re_eyx210[i][j][k]-term2*(re_hz[imax-i][j][k]-re_hz[imax-i-1][j][k])/dx;
				im_eyx210[i][j][k]=term1*im_eyx210[i][j][k]-term2*(im_hz[imax-i][j][k]-im_hz[imax-i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx012[i][j][k]=term1*re_eyx012[i][j][k]-term2*(re_hz[i][j][kmax-k]-re_hz[i-1][j][kmax-k])/dx;
				im_eyx012[i][j][k]=term1*im_eyx012[i][j][k]-term2*(im_hz[i][j][kmax-k]-im_hz[i-1][j][kmax-k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx212[i][j][k]=term1*re_eyx212[i][j][k]-term2*(re_hz[imax-i][j][kmax-k]-re_hz[imax-i-1][j][kmax-k])/dx;
				im_eyx212[i][j][k]=term1*im_eyx212[i][j][k]-term2*(im_hz[imax-i][j][kmax-k]-im_hz[imax-i-1][j][kmax-k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx011[i][j][k]=term1*re_eyx011[i][j][k]-term2*(re_hz[i][j][k]-re_hz[i-1][j][k])/dx;
				im_eyx011[i][j][k]=term1*im_eyx011[i][j][k]-term2*(im_hz[i][j][k]-im_hz[i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx211[i][j][k]=term1*re_eyx211[i][j][k]-term2*(re_hz[imax-i][j][k]-re_hz[imax-i-1][j][k])/dx;
				im_eyx211[i][j][k]=term1*im_eyx211[i][j][k]-term2*(im_hz[imax-i][j][k]-im_hz[imax-i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx110[i][j][k]=term1*re_eyx110[i][j][k]-term2*(re_hz[i][j][k]-re_hz[i-1][j][k])/dx;
				im_eyx110[i][j][k]=term1*im_eyx110[i][j][k]-term2*(im_hz[i][j][k]-im_hz[i-1][j][k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_eyx112[i][j][k]=term1*re_eyx112[i][j][k]-term2*(re_hz[i][j][kmax-k]-re_hz[i-1][j][kmax-k])/dx;
				im_eyx112[i][j][k]=term1*im_eyx112[i][j][k]-term2*(im_hz[i][j][kmax-k]-im_hz[i-1][j][kmax-k])/dx;

				//Eyz
				u1=(sigma_z*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx010[i][j][k]=term1*re_ezx010[i][j][k]+term2*(re_hy[i][j][k]-re_hy[i-1][j][k])/dx;
				im_ezx010[i][j][k]=term1*im_ezx010[i][j][k]+term2*(im_hy[i][j][k]-im_hy[i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx210[i][j][k]=term1*re_ezx210[i][j][k]+term2*(re_hy[imax-i][j][k]-re_hy[imax-i-1][j][k])/dx;
				im_ezx210[i][j][k]=term1*im_ezx210[i][j][k]+term2*(im_hy[imax-i][j][k]-im_hy[imax-i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx012[i][j][k]=term1*re_ezx012[i][j][k]+term2*(re_hy[i][j][kmax-1-k]-re_hy[i-1][j][kmax-1-k])/dx;
				im_ezx012[i][j][k]=term1*im_ezx012[i][j][k]+term2*(im_hy[i][j][kmax-1-k]-im_hy[i-1][j][kmax-1-k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx212[i][j][k]=term1*re_ezx212[i][j][k]+term2*(re_hy[imax-i][j][kmax-1-k]-re_hy[imax-i-1][j][kmax-1-k])/dx;
				im_ezx212[i][j][k]=term1*im_ezx212[i][j][k]+term2*(im_hy[imax-i][j][kmax-1-k]-im_hy[imax-i-1][j][kmax-1-k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx011[i][j][k]=term1*re_ezx011[i][j][k]+term2*(re_hy[i][j][k]-re_hy[i-1][j][k])/dx;
				im_ezx011[i][j][k]=term1*im_ezx011[i][j][k]+term2*(im_hy[i][j][k]-im_hy[i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx211[i][j][k]=term1*re_ezx211[i][j][k]+term2*(re_hy[imax-i][j][k]-re_hy[imax-i-1][j][k])/dx;
				im_ezx211[i][j][k]=term1*im_ezx211[i][j][k]+term2*(im_hy[imax-i][j][k]-im_hy[imax-i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx110[i][j][k]=term1*re_ezx110[i][j][k]+term2*(re_hy[i][j][k]-re_hy[i-1][j][k])/dx;
				im_ezx110[i][j][k]=term1*im_ezx110[i][j][k]+term2*(im_hy[i][j][k]-im_hy[i-1][j][k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_ezx112[i][j][k]=term1*re_ezx112[i][j][k]+term2*(re_hy[i][j][kmax-1-k]-re_hy[i-1][j][kmax-1-k])/dx;
				im_ezx112[i][j][k]=term1*im_ezx112[i][j][k]+term2*(im_hy[i][j][kmax-1-k]-im_hy[i-1][j][kmax-1-k])/dx;

				//Ezy
				u1=(sigma_y*PML_coefficient)/2.0;
				u2=PML_coefficient;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//				u2=cnstE[i][j][k];
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz010[i][j][k]=term1*re_hxz010[i][j][k]+term2*(re_ey[i][j][k+1]-re_ey[i][j][k])/dz;
				im_hxz010[i][j][k]=term1*im_hxz010[i][j][k]+term2*(im_ey[i][j][k+1]-im_ey[i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz210[i][j][k]=term1*re_hxz210[i][j][k]+term2*(re_ey[imax-i][j][k+1]-re_ey[imax-i][j][k])/dz;
				im_hxz210[i][j][k]=term1*im_hxz210[i][j][k]+term2*(im_ey[imax-i][j][k+1]-im_ey[imax-i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz012[i][j][k]=term1*re_hxz012[i][j][k]+term2*(re_ey[i][j][kmax-1-k+1]-re_ey[i][j][kmax-1-k])/dz;
				im_hxz012[i][j][k]=term1*im_hxz012[i][j][k]+term2*(im_ey[i][j][kmax-1-k+1]-im_ey[i][j][kmax-1-k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz212[i][j][k]=term1*re_hxz212[i][j][k]+term2*(re_ey[imax-i][j][kmax-1-k+1]-re_ey[imax-i][j][kmax-1-k])/dz;
				im_hxz212[i][j][k]=term1*im_hxz212[i][j][k]+term2*(im_ey[imax-i][j][kmax-1-k+1]-im_ey[imax-i][j][kmax-1-k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz011[i][j][k]=term1*re_hxz011[i][j][k]+term2*(re_ey[i][j][k+1]-re_ey[i][j][k])/dz;
				im_hxz011[i][j][k]=term1*im_hxz011[i][j][k]+term2*(im_ey[i][j][k+1]-im_ey[i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz211[i][j][k]=term1*re_hxz211[i][j][k]+term2*(re_ey[imax-i][j][k+1]-re_ey[imax-i][j][k])/dz;
				im_hxz211[i][j][k]=term1*im_hxz211[i][j][k]+term2*(im_ey[imax-i][j][k+1]-im_ey[imax-i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz110[i][j][k]=term1*re_hxz110[i][j][k]+term2*(re_ey[i][j][k+1]-re_ey[i][j][k])/dz;
				im_hxz110[i][j][k]=term1*im_hxz110[i][j][k]+term2*(im_ey[i][j][k+1]-im_ey[i][j][k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hxz112[i][j][k]=term1*re_hxz112[i][j][k]+term2*(re_ey[i][j][kmax-1-k+1]-re_ey[i][j][kmax-1-k])/dz;
				im_hxz112[i][j][k]=term1*im_hxz112[i][j][k]+term2*(im_ey[i][j][kmax-1-k+1]-im_ey[i][j][kmax-1-k])/dz;

				//Hxy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx010[i][j][k]=term1*re_hyx010[i][j][k]+term2*(re_ez[i+1][j][k]-re_ez[i][j][k])/dx;
				im_hyx010[i][j][k]=term1*im_hyx010[i][j][k]+term2*(im_ez[i+1][j][k]-im_ez[i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx210[i][j][k]=term1*re_hyx210[i][j][k]+term2*(re_ez[imax-1-i+1][j][k]-re_ez[imax-1-i][j][k])/dx;
				im_hyx210[i][j][k]=term1*im_hyx210[i][j][k]+term2*(im_ez[imax-1-i+1][j][k]-im_ez[imax-1-i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx012[i][j][k]=term1*re_hyx012[i][j][k]+term2*(re_ez[i+1][j][kmax-1-k]-re_ez[i][j][kmax-1-k])/dx;
				im_hyx012[i][j][k]=term1*im_hyx012[i][j][k]+term2*(im_ez[i+1][j][kmax-1-k]-im_ez[i][j][kmax-1-k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx212[i][j][k]=term1*re_hyx212[i][j][k]+term2*(re_ez[imax-1-i+1][j][kmax-1-k]-re_ez[imax-1-i][j][kmax-1-k])/dx;
				im_hyx212[i][j][k]=term1*im_hyx212[i][j][k]+term2*(im_ez[imax-1-i+1][j][kmax-1-k]-im_ez[imax-1-i][j][kmax-1-k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx011[i][j][k]=term1*re_hyx011[i][j][k]+term2*(re_ez[i+1][j][k]-re_ez[i][j][k])/dx;
				im_hyx011[i][j][k]=term1*im_hyx011[i][j][k]+term2*(im_ez[i+1][j][k]-im_ez[i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx211[i][j][k]=term1*re_hyx211[i][j][k]+term2*(re_ez[imax-1-i+1][j][k]-re_ez[imax-1-i][j][k])/dx;
				im_hyx211[i][j][k]=term1*im_hyx211[i][j][k]+term2*(im_ez[imax-1-i+1][j][k]-im_ez[imax-1-i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx110[i][j][k]=term1*re_hyx110[i][j][k]+term2*(re_ez[i+1][j][k]-re_ez[i][j][k])/dx;
				im_hyx110[i][j][k]=term1*im_hyx110[i][j][k]+term2*(im_ez[i+1][j][k]-im_ez[i][j][k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hyx112[i][j][k]=term1*re_hyx112[i][j][k]+term2*(re_ez[i+1][j][kmax-1-k]-re_ez[i][j][kmax-1-k])/dx;
				im_hyx112[i][j][k]=term1*im_hyx112[i][j][k]+term2*(im_ez[i+1][j][kmax-1-k]-im_ez[i][j][kmax-1-k])/dx;

				//Hyz
				u1=(sigma_z*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_z*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx010[i][j][k]=term1*re_hzx010[i][j][k]-term2*(re_ey[i+1][j][k]-re_ey[i][j][k])/dx;
				im_hzx010[i][j][k]=term1*im_hzx010[i][j][k]-term2*(im_ey[i+1][j][k]-im_ey[i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx210[i][j][k]=term1*re_hzx210[i][j][k]-term2*(re_ey[imax-1-i+1][j][k]-re_ey[imax-1-i][j][k])/dx;
				im_hzx210[i][j][k]=term1*im_hzx210[i][j][k]-term2*(im_ey[imax-1-i+1][j][k]-im_ey[imax-1-i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx012[i][j][k]=term1*re_hzx012[i][j][k]-term2*(re_ey[i+1][j][kmax-k]-re_ey[i][j][kmax-k])/dx;
				im_hzx012[i][j][k]=term1*im_hzx012[i][j][k]-term2*(im_ey[i+1][j][kmax-k]-im_ey[i][j][kmax-k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx212[i][j][k]=term1*re_hzx212[i][j][k]-term2*(re_ey[imax-1-i+1][j][kmax-k]-re_ey[imax-1-i][j][kmax-k])/dx;
				im_hzx212[i][j][k]=term1*im_hzx212[i][j][k]-term2*(im_ey[imax-1-i+1][j][kmax-k]-im_ey[imax-1-i][j][kmax-k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx011[i][j][k]=term1*re_hzx011[i][j][k]-term2*(re_ey[i+1][j][k]-re_ey[i][j][k])/dx;
				im_hzx011[i][j][k]=term1*im_hzx011[i][j][k]-term2*(im_ey[i+1][j][k]-im_ey[i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx211[i][j][k]=term1*re_hzx211[i][j][k]-term2*(re_ey[imax-1-i+1][j][k]-re_ey[imax-1-i][j][k])/dx;
				im_hzx211[i][j][k]=term1*im_hzx211[i][j][k]-term2*(im_ey[imax-1-i+1][j][k]-im_ey[imax-1-i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx110[i][j][k]=term1*re_hzx110[i][j][k]-term2*(re_ey[i+1][j][k]-re_ey[i][j][k])/dx;
				im_hzx110[i][j][k]=term1*im_hzx110[i][j][k]-term2*(im_ey[i+1][j][k]-im_ey[i][j][k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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

//�ύX�ꏊ(�n��)
//				u1=(sigma_x*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

				u2=cnstM;
				term1=(1.0-u1)/(1.0+u1);
				term2=u2/(1.0+u1);

				re_hzx112[i][j][k]=term1*re_hzx112[i][j][k]-term2*(re_ey[i+1][j][kmax-k]-re_ey[i][j][kmax-k])/dx;
				im_hzx112[i][j][k]=term1*im_hzx112[i][j][k]-term2*(im_ey[i+1][j][kmax-k]-im_ey[i][j][kmax-k])/dx;

				//Hzy
				u1=(sigma_y*PML_coefficient)/2.0;

//�ύX�ꏊ(�n��)
//				u1=(sigma_y*cnstE[i][j][k])/2.0;
//�ύX�ꏊ(�I��)

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
/********************************** �z�����E�����I��� **********************************/


/*************************** ��̓��f�����t�B�[���h�̓f���o�� ***************************/
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

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

}

void file_open_field(){
	if(ikk<10) sprintf(fikk,"%d_0%d_%05d",m,ikk,n); else sprintf(fikk,"%d_%d_%05d",m,ikk,n);

//	sprintf(fname,"%sEx_field.xls",fikk);EX=fopen(fname,"w");
//	sprintf(fname,"%sEy_field.xls",fikk);EY=fopen(fname,"w");
//	sprintf(fname,"%sEz_field.xls",fikk);EZ=fopen(fname,"w");
//	sprintf(fname,"%sHx_field.xls",fikk);HX=fopen(fname,"w");
//	sprintf(fname,"%sHy_field.xls",fikk);HY=fopen(fname,"w");
//	sprintf(fname,"%sHz_field.xls",fikk);HZ=fopen(fname,"w");

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

}

void field_print(int n){
//	int i,j;

//�ύX�ꏊ(�n��)
	int i,j,k;
//�ύX�ꏊ(�I��)

	printf("Output Field t = %d\n",n);

	for(i=L;i<imax-L;i++){
		for(j=0;j<jmax;j++){
//			fprintf(EX,"%e\t",re_ex[i][j][ck]);
//			fprintf(EY,"%e\t",re_ey[i][j][ck]);
//			fprintf(EZ,"%e\t",re_ez[i][j][ck]);
//			fprintf(HX,"%e\t",re_hx[i][j][ck]);
//			fprintf(HY,"%e\t",re_hy[i][j][ck]);
//			fprintf(HZ,"%e\t",re_hz[i][j][ck]);

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

		}
//		fprintf(EX,"\n");
//		fprintf(EY,"\n");
//		fprintf(EZ,"\n");
//		fprintf(HX,"\n");
//		fprintf(HY,"\n");
//		fprintf(HZ,"\n");

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

	}

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

}


void file_close_field(){
//	fclose(EX);
//	fclose(EY);
//	fclose(EZ);
//	fclose(HX);
//	fclose(HY);
//	fclose(HZ);

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

}


void file_open_field_f(){
	if(ikk<10) sprintf(fikk,"%d_0%d_%05d",m,ikk,n); else sprintf(fikk,"%d_%d_%05d",m,ikk,n);

//	sprintf(fname,"%sEx_field.xls",fikk);EX=fopen(fname,"w");
//	sprintf(fname,"%sEy_field.xls",fikk);EY=fopen(fname,"w");
//	sprintf(fname,"%sEz_field.xls",fikk);EZ=fopen(fname,"w");
//	sprintf(fname,"%sHx_field.xls",fikk);HX=fopen(fname,"w");
//	sprintf(fname,"%sHy_field.xls",fikk);HY=fopen(fname,"w");
//	sprintf(fname,"%sHz_field.xls",fikk);HZ=fopen(fname,"w");

//�ύX�ꏊ(�n��)
	if(excite == 1){
		sprintf(fname,"%sHz_XZ_center_field.txt",fikk);HZ_XZ_center=fopen(fname,"w");
	}
	if(excite == 2){
		sprintf(fname,"%sEz_XZ_center_field.txt",fikk);EZ_XZ_center=fopen(fname,"w");
	}
//�ύX�ꏊ(�I��)

}


void field_print_f(int n){

//�ύX�ꏊ(�n��)
	int i,j,k;
//�ύX�ꏊ(�I��)

	printf("Output Field t = %d\n",n);
    i = ci;
	k = ck;

//�ύX�ꏊ(�n��)

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

//�ύX�ꏊ(�n��)
	if(excite == 1){
			fclose(HZ_XZ_center);
	}
	if(excite == 2){
			fclose(EZ_XZ_center);	
	}
//�ύX�ꏊ(�I��)

}



/*
��̓��f���C�t�B�[���h�Ƃ�z���W��ck�݂̂̕\���Ƃ����D
x-z or y-z�����ŏo�͂������ꍇ�͎蓮�ŕς��邱�ƁD
*/
/************************ ��̓��f�����t�B�[���h�̓f���o���I��� ************************/


/******************************** �d���E�̃p���X�f���o�� ********************************/
void file_open(){
	if(ikk<10) sprintf(fikk,"%d_0%d",m,ikk); else sprintf(fikk,"%d_%d",m,ikk);

//�ύX�ꏊv5.03(�n��)
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
//�ύX�ꏊ(�I��)

void Poyntingprint(int flag){
	int X1=ci;		int X2=ci+imove;	int X3=ci-imove;
	int Y1=cj;		int Y2=cj+jmove;	int Y3=cj-jmove;

//�ύX�ꏊ(�n��)
	int Z1=ck;		int Z2=ck+kmove;	int Z3=ck-kmove;
//�ύX�ꏊ(�I��)

	if(flag==0){
		printf("����_");
		printf("\n");
//		printf("x=%d y=%d\n",X1,Y1);		// Pin(ci,cj)��͋�Ԃ̒��S
//		printf("x=%d y=%d\n",X2,Y2);		// Pout1(ci+imove,cj+jmove)
//		printf("x=%d y=%d\n",X2,Y3);		// Pout2(ci+imove,cj-jmove)
//		printf("x=%d y=%d\n",X3,Y2);		// Pout3(ci-imove,cj+jmove)
//		printf("x=%d y=%d\n",X3,Y3);		// Pout4(ci-imove,cj-jmove)

//�ύX�ꏊ(�n��)
		printf("x=%d	y=%d	z=%d\n",X1,Y1,Z1);		// Pin(ci,cj)��͋�Ԃ̒��S
		printf("x=%d	y=%d	z=%d\n",X2,Y2,Z1);		// Pout(ci+imove,cj+jmove,ck)
		printf("x=%d	y=%d	z=%d\n",X1,Y3,Z1);		// Pout(ci,cj-jmove,ck)
		printf("x=%d	y=%d	z=%d\n",X1,Y2,Z3);		// Pout(ci,cj+jmove,ck-kmove)
		printf("x=%d	y=%d	z=%d\n",X3,Y3,Z2);		// Pout(ci-imove,cj-jmove,ck+kmove)
//�ύX�ꏊ(�I��)

	}

	if(flag==1){
//		fprintf(fPoynt1Ex,"%e\t%e\t%e\t%e\t%e\n",re_ex[X1][Y1][ck],re_ex[X2][Y2][ck],re_ex[X2][Y3][ck],re_ex[X3][Y2][ck],re_ex[X3][Y3][ck]);
//		fprintf(fPoynt1Ey,"%e\t%e\t%e\t%e\t%e\n",re_ey[X1][Y1][ck],re_ey[X2][Y2][ck],re_ey[X2][Y3][ck],re_ey[X3][Y2][ck],re_ey[X3][Y3][ck]);
//		fprintf(fPoynt1Ez,"%e\t%e\t%e\t%e\t%e\n",re_ez[X1][Y1][ck],re_ez[X2][Y2][ck],re_ez[X2][Y3][ck],re_ez[X3][Y2][ck],re_ez[X3][Y3][ck]);
//		fprintf(fPoynt1Hx,"%e\t%e\t%e\t%e\t%e\n",re_hx[X1][Y1][ck],re_hx[X2][Y2][ck],re_hx[X2][Y3][ck],re_hx[X3][Y2][ck],re_hx[X3][Y3][ck]);
//		fprintf(fPoynt1Hy,"%e\t%e\t%e\t%e\t%e\n",re_hy[X1][Y1][ck],re_hy[X2][Y2][ck],re_hy[X2][Y3][ck],re_hy[X3][Y2][ck],re_hy[X3][Y3][ck]);
//		fprintf(fPoynt1Hz,"%e\t%e\t%e\t%e\t%e\n",re_hz[X1][Y1][ck],re_hz[X2][Y2][ck],re_hz[X2][Y3][ck],re_hz[X3][Y2][ck],re_hz[X3][Y3][ck]);

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)
	
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
/***************************** �d���E�̃p���X�f���o���I��� *****************************/


/************************************* �z��̏����� *************************************/
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
/********************************** �z��̏������I��� **********************************/


/******************************** �t�[���G�ϊ��v���O���� ********************************/
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

//�ύX�ꏊ(�n��)
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
//�ύX�ꏊ(�I��)

	ws=2*PI/dt;

	printf("FFT start\n");
	for(m=K_STARTPOINT;m<=K_ENDPOINT;m++){
		for(ikk=IKK_STARTPOINT;ikk<=IKK_ENDPOINT;ikk++){
			printf("Direction = %d, ikk = %d\n",m,ikk);
			for(f_ele=0;f_ele<ELE;f_ele++){
				printf("����FFT�� : %s\n",file[f_ele]);

				if(ikk<10){
					sprintf(fikk,"0%d",ikk);
				}
				else{
					sprintf(fikk,"%d",ikk);
				}
//				sprintf(InputFile,"%d_%s%s.xls",m,fikk,file[f_ele]);
				sprintf(InputFile,"%d_%s%s.txt",m,fikk,file[f_ele]);	//v5.11
				if((fp1=fopen(InputFile,"r"))==NULL){
					printf("�t�@�C����������܂���D---- %s\n",InputFile);
					exit(1);
				}
				for(step=0,count=0;step<STEP;step++){
					for(f_num=0;f_num<F_NUMBER;f_num++){
						if(fscanf(fp1,"%e,",&val)!=EOF){	//�ύX�ꏊv5.12
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
					printf("xr���m�ۂł��܂��� �������s��\n");
					exit(1);
				}
				if((data[f_ele+1]=(double *)calloc(1,sizeof(double)*nf2))==NULL){
					printf("data %d���m�ۂł��܂��� �������s��\n",(f_ele+1));
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
				printf("data 0���m�ۂł��܂��� �������s��\n");
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
/***************************** �t�[���G�ϊ��v���O�����I��� *****************************/


/********************************* �p�����[�^�̓f���o�� *********************************/
void Output_Condition(void){
	FILE *OC;
	OC = fopen("Condition.txt","w");
//	fprintf(OC,"AIR = %d\nSI = %d\nSIO2 = %d\nCIRCLE_MEDIUM = %d\nimax = %d\njmax = %d\nkmax = %d\nci = %d\ncj = %d\nck = %d\nPITCH = %.4e\nRADIUS = %.4e\nCORE = %.4e\nCLAD = %.4e\nNORMALFREQ = %.2f\nCELLLAMBDA = %.4e\nsigma = %.4f\nPULSEPEAK = %d\nREF_INDEX1 = %.4f\nnmax = %d\ndx = %e\ndy = %e\ndz = %e\ndt = %e"
//		,AIR,SI,SIO2,CIRCLE_MEDIUM,imax,jmax,kmax,ci,cj,ck,PITCH,RADIUS,CORE,CLAD,NORMALFREQ,CELLLAMBDA,(sigma/omega0),PULSEPEAK,REF_INDEX1,nmax,dx,dy,dz,dt);

//�ύX�ꏊ(�n��)
	fprintf(OC,"AIR = %d\nInP = %d\nAct = %d\nSIO2 = %d\nCIRCLE_MEDIUM = %d\nimax = %d\njmax = %d\nkmax = %d\nci = %d\ncj = %d\nck = %d\nInP_INDEX = %.3f\nAct_INDEX = %.3f\nClad_INDEX = %.3f\nPITCH = %.4e\nB = %d\nRADIUS = %.4e\nRADIUS2 = %.4e\nAir-hole Cells = %d\nEffective 2r/a = %.4f\nAir_clad = %.4e\nInP_over = %.4e\nInP_down = %.4e\nInP_sub = %.4e\nAct = %.4e\nWidth = %d\n3rd_Shift = %.3e\nLAMBDA = %.3e\nsigma = %.4f\nPULSEPEAK = %d\nnmax = %d\ndx = %e\ndy = %e\ndz = %e\ndt = %e\nk vector = %d/%d\n"

		,AIR,InP_MEDIUM,Act_MEDIUM,SIO2,CIRCLE_MEDIUM
		,imax,jmax,kmax,ci,cj,ck
		,InP_INDEX,Act_INDEX,REF_INDEX2
		,PITCH,B,RADIUS,RADIUS2,effective_circle_cell,effective_2r_a,AIR_clad,InP_circle_over,InP_circle_down,InP_sub,Act,Width,THIRD_Sx
		,LAMBDA,(sigma/omega0),PULSEPEAK,nmax
		,dx,dy,dz,dt
		,IKK_ENDPOINT,IKK_MAX);
//�ύX�ꏊ(�I��)

	fclose(OC);
}
/****************************** �p�����[�^�̓f���o���I��� ******************************/
