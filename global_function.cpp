#define INT_COMP 1e-2								// intŒ^“¯m‚ÌœZ‚Ì‚Æ‚«‚ÌŒë·•â³—pŒW”
#define CARRY_DIV_ORD 1e+10								// Œ…ã‚°œZ‚Å‚ÌŒ…ã‚°‚Ì”
//#define INT_DIV(x, y) ((x * CARRY_DIV_ORD) / (y * CARRY_DIV_ORD) + INT_COMP)		// Œ…ã‚°œZ
#define INT_DIV(x, y) ((int) (((double) x / (double) y) + INT_COMP))		// Œ…ã‚°œZ
#define SQ(x) (x * x)									// 2æ
