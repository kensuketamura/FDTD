#define INT_COMP 1e-2								// int�^���m�̏��Z�̂Ƃ��̌덷�␳�p�W��
#define CARRY_DIV_ORD 1e+10								// ���グ���Z�ł̌��グ�̐�
//#define INT_DIV(x, y) ((x * CARRY_DIV_ORD) / (y * CARRY_DIV_ORD) + INT_COMP)		// ���グ���Z
#define INT_DIV(x, y) ((int) (((double) x / (double) y) + INT_COMP))		// ���グ���Z
#define SQ(x) (x * x)									// 2��
