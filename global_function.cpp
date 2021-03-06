#define INT_COMP 1e-2								// int型同士の除算のときの誤差補正用係数
#define CARRY_DIV_ORD 1e+10								// 桁上げ除算での桁上げの数
//#define INT_DIV(x, y) ((x * CARRY_DIV_ORD) / (y * CARRY_DIV_ORD) + INT_COMP)		// 桁上げ除算
#define INT_DIV(x, y) ((int) (((double) x / (double) y) + INT_COMP))		// 桁上げ除算
#define SQ(x) (x * x)									// 2乗
