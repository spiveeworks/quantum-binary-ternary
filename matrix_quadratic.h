#pragma once

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define range(i, n) for (size_t i = 0; i < (n); i++)

int gcd(int x, int y) {
	while (y != 0) {
		int rem = x % y;
		x = y;
		y = rem;
	}
	if (x < 0) {
		x = -x;
	}
	return x;
}

// stupidly slow when x is big, but intended use case is around int_sqrt(36)=6
unsigned int_sqrt(unsigned x) {
	unsigned i = 1;
	while (i*i <= x) {
		i++;
	}
	return i-1;
}
unsigned int_quadpart(unsigned x) {
	unsigned i = 2;
	unsigned result = 1;
	while (i*i <= x) {
		if (x % (i*i) == 0) {
			result *= i;
			x /= i*i;
		} else {
			i++;
		}
	}
	return result;
}

// x represents the real number x.wh * sqrt(x.qu)
typedef struct radnum {
	int whole;
	unsigned qu;
} radnum;

radnum radnum_reduce(radnum x) {
	if (x.whole == 0 || x.qu == 0) {
		return (radnum){0,1};
	}
	unsigned q = int_quadpart(x.qu);
	return (radnum){x.whole*(int)q, x.qu/(q*q)};
}

radnum radnum_add(radnum x, radnum y) {
	if (x.whole == 0) {
		return y;
	}
	if (y.whole == 0) {
		return x;
	}
	if (x.qu != y.qu) {
		printf("Incompatible radicals: %u and %u\n", x.qu, y.qu);
		exit(1);
	}
	return (radnum){x.whole+y.whole, x.qu};
}

radnum radnum_mul(radnum x, radnum y) {
	return (radnum){x.whole*y.whole, x.qu*y.qu};
}

void radnum_print(radnum x) {
	if (x.whole < 0) {
		printf("-");
		x.whole = -x.whole;
	}
	if (x.whole != 1 || x.qu == 1) {
		printf("%d", x.whole);
	}
	if (x.qu != 1) {
		printf("sqrt(%d)", x.qu);
	}
}

typedef struct num {
	radnum re;
	radnum im;
	unsigned de;
} num;

num num_reduce(num x) {
	x.re = radnum_reduce(x.re);
	x.im = radnum_reduce(x.im);
	int d = gcd(gcd(x.re.whole, x.im.whole), (int)x.de);
	x.re.whole /= (int)d;
	x.im.whole /= (int)d;
	x.de /= d;
	return x;
}

void num_print(num x) {
	//printf("{%d, %u, %d, %u, %u}", x.re.whole, x.re.qu, x.im.whole, x.im.qu, x.de);
	bool paren = x.re.whole != 0 && x.im.whole != 0 && x.de != 1;
	if (paren) {
		printf("(");
	}
	if (x.re.whole == 0 && x.im.whole == 0) {
		printf("0");
	} else {
		if (x.re.whole != 0) {
			radnum_print(x.re);
		}
		if (x.re.whole != 0 && x.im.whole > 0) {
			printf("+");
		}
		if (x.im.whole == 1 && x.im.qu == 1) {
			printf("i");
		} else if (x.im.whole == -1 && x.im.qu == 1) {
			printf("-i");
		} else if (x.im.whole != 0) {
			radnum_print(x.im);
			printf("i");
		}
	}
	if (paren) {
		printf(")");
	}
	if (x.de != 1) {
		printf("/%u", x.de);
	}
}

num num_add(num x, num y) {
	x.re.whole *= (int)y.de;
	x.im.whole *= (int)y.de;
	y.re.whole *= (int)x.de;
	y.im.whole *= (int)x.de;
	return (num){radnum_add(x.re, y.re), radnum_add(x.im, y.im), x.de*y.de};
}

num num_mul(num x, num y) {
	num result;
	// x.de *= y.de;
	// y.de *= x.de;
	radnum xryr = radnum_reduce(radnum_mul(x.re, y.re));
	radnum xiyr = radnum_reduce(radnum_mul(x.im, y.re));
	radnum xryi = radnum_reduce(radnum_mul(x.re, y.im));
	radnum xiyi = radnum_reduce(radnum_mul(x.im, y.im));
	xiyi.whole = -xiyi.whole;
	result.re = radnum_add(xryr, xiyi);
	result.im = radnum_add(xryi, xiyr);
	result.de = x.de * y.de;
	return result;
}

num num_mod_sq(num x) {
	int mod_sq = x.re.whole * x.re.whole * (int)x.re.qu
		+ x.im.whole * x.im.whole * (int)x.im.qu;
	return (num){{mod_sq, 1}, {0, 1}, x.de*x.de};
}

num num_inv(num x) {
	int mod_sq = x.re.whole * x.re.whole * (int)x.re.qu
		+ x.im.whole * x.im.whole * (int)x.im.qu;
	x.re.whole *= (int)x.de;
	x.im.whole *= -(int)x.de;
	x.de = mod_sq;
	return x;
}

num num_inv_sqrt(num x) {
	if (x.im.whole != 0 && x.im.qu != 0) {
		printf("Tried to take invsqrt of non-real number\n");
		exit(1);
	}
	if (x.re.qu != 1) {
		printf("Tried to take invsqrt of non-rational number\n");
		exit(1);
	}
	int sgn = 1;
	int nu = (int)x.de;
	int de = x.re.whole;
	if (de < 0) {
		sgn = -1;
		de = -de;
	}
	return (num){{sgn, nu * de}, {0, 1}, de};
}

/* My matrices are just pointers to num
 * (assumed to be pointers to an array of size n*n)
 * C allows statically declared arrays, and we would like `num x[4][4];` to
 * declare a usable matrix, and `x[i][j]` to extract an appropriate term
 * this means `x[i]` will yield a `num[4]` which by convention ought to be a
 * row vector, and so the equivalent index from a flat array will be
 * `x[i*n + j]`
 *
 * we could use a newtype like `typedef struct {num *mp;} mat;`, possibly add a
 * size field as well, but this does not save much trouble, since we aren't yet
 * using vectors, and usually share size between multiple objects.
 *
 * I will probably come to regret this if I continue to use this code.
 */

// performs matrix multiplication
void mat_mul(int n, num *out, num *x, num *y) {
	range(i, n) {
		range (k, n) {
			num sum = {{0, 1}, {0, 1}, 1};
			range (j, n) {
				sum = num_add(sum, num_mul(x[i*n+j], y[j*n+k]));
			}
			out[i*n+k] = sum;
		}
	}
}

void mat_reduce(int n, num *x) {
	range (i, n) {
		range (j, n) {
			x[i*n+j] = num_reduce(x[i*n+j]);
		}
	}
}

void mat_print(int n, num *x, char*sepcol, char*seprow) {
	range (i, n) {
		range (j, n) {
			num_print(x[i*n+j]);
			printf("%s", sepcol);
		}
		printf("%s", seprow);
	}
}

