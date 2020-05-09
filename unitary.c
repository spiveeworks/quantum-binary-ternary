#include <stdio.h>
#include <stdlib.h>

#define range(i, n) for (size_t i = 0; i < (n); i++)

int gcd(int x, int y) {
	while (y != 0) {
		int rem = x % y;
		x = y;
		y = rem;
	}
	return x;
}

// stupidly slow when x is big, but intended use case is around int_sqrt(36)=6
int int_sqrt(int x) {
	int i = 1;
	while (i*i <= x) {
		i++;
	}
	return i-1;
}
int int_quadpart(int x) {
	int i = 2;
	int result = 1;
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

// x represents the complex number (x.re + i*x.im)/sqrt(x.qu)
typedef struct num {
	int re;
	int im;
	int qu;
} num;

num num_mul(num x, num y) {
	num result;
	result.re = x.re * y.re - x.im * y.im;
	result.im = x.re * y.im + x.im * y.re;
	result.qu = x.qu * y.qu;
	return result;
}

num num_add(num x, num y) {
	num result;
	if (x.re == 0 && x.im == 0) {
		return y;
	}
	if (y.re == 0 && y.im == 0) {
		return x;
	}
	int dqu = gcd(x.qu, y.qu);
	int xmodqu = x.qu / dqu;
	int ymodqu = y.qu / dqu;
	int xmod = int_sqrt(xmodqu);
	int ymod = int_sqrt(ymodqu);
	if (xmod*xmod != xmodqu || ymod * ymod != ymodqu) {
		printf("Incompatible radicals: %d and %d\n", x.qu, y.qu);
		exit(1);
	}
	result.re = x.re * ymod + y.re * xmod;
	result.im = x.im * ymod + y.im * xmod;
	result.qu = x.qu * ymodqu;
	return result;
}

// these (as with all num_* functions) are convenience functions, only use them
// when they are MORE succinct than writing the arithmetic out by hand
num num_neg(num x) {
	return (num){-x.re,-x.im,x.qu};
}
num num_conj(num x) {
	return (num){x.re,-x.im, x.qu};
}

// have to rationalize denominator and realize numerator
// sqrt(q)/(x+iy) = q(x-iy)/sqrt(q)(x^2+y^2)
num num_inv(num x) {
	num result;
	result.re = x.re * x.qu;
	result.im = -x.im * x.qu;
	int mod = x.re * x.re + x.im * x.im;
	result.qu = x.qu * mod * mod;
	return result;
}

// assuming a number is real and rational, take the inverse square root
// (nr/dr)^-0.5 = d/sqrt(n*d)
num num_invsqrt(num x) {
	if (x.im != 0) {
		printf("Tried to take invsqrt of non-real number\n");
		exit(1);
	}
	int d = int_sqrt(x.qu);
	if (d*d != x.qu) {
		printf("Tried to take invsqrt of irrational number\n");
		exit(1);
	}
	int r = gcd(x.re, d);
	if (r < 0) {
		r = -r;
	}
	int n = x.re / r;
	d = d / r;
	return (num){d, 0, n*d};
}

num num_mod_sq(num x) {
	return (num){x.re*x.re + x.im*x.im, 0, x.qu*x.qu};
}

num num_mod(num x) {
	int sq = x.re*x.re + x.im*x.im;
	return (num){sq, 0, sq*x.qu};
}

/*
int num_as_int(num x) {
	if (x.qu == 0) {
		printf("Tried to interpret non-finite num as integer\n");
		exit(1);
	}
	if (x.im != 0) {
		printf("Tried to interpret non-real num as integer\n");
		exit(1);
	}
	if (x.re == 0) {
		return 0;
	}
	int d = int_sqrt(x.qu);
	if (d*d != x.qu) {
		printf("Tried to interpret irrational num as integer\n");
		exit(1);
	}
	int result = x.re / d;
	if (result * d != x.re) {
		printf("Tried to interpret fractional num as integer\n");
		exit(1);
	}
	return result;
}
*/

num num_reduce(num x) {
	int d = int_quadpart(x.qu);
	d = gcd(d, x.re);
	d = gcd(d, x.im);
	if (d < 0) {
		d = -d;
	}
	return (num){x.re/d, x.im/d, x.qu/(d*d)};
}

void num_print(num x) {
	if (x.re == 0 && x.im == 0) {
		printf("0");
	} else {
		if (x.re != 0 && x.im != 0 && x.qu != 1) {
			printf("(");
		}
		if (x.re != 0) {
			printf("%d", x.re);
		}
		if (x.re != 0 && x.im > 0) {
			printf("+");
		}
		if (x.im == 1) {
			printf("i");
		} else if (x.im == -1) {
			printf("-i");
		} else if (x.im != 0) {
			printf("%di", x.im);
		}
		if (x.re != 0 && x.im != 0 && x.qu != 1) {
			printf(")");
		}
		if (x.qu != 1) {
			printf("/");
			int d = int_quadpart(x.qu);
			int qu = x.qu / (d*d);
			if (d != 1) {
				printf("%d", d);
			}
			if (qu != 1) {
				printf("sqrt(%d)", qu);
			}
		}
	}
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
			num sum = {0, 0, 1};
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

// generates a two-level matrix to eliminate an entry from a given matrix
void eliminator(int n, int i, int j, num *out, num *x) {
	num a = x[j*n+j];
	num b = x[i*n+j];
	num modsq = num_add(num_mod_sq(a), num_mod_sq(b));
	num mod = num_invsqrt(modsq);
	// printf("|"); num_print(a); printf("|^2 + |"); num_print(b); printf("|^2 = "); num_print(modsq); printf("\n");
	// printf("1/sqrt("); num_print(modsq); printf(") = "); num_print(mod); printf("\n");
	a = num_mul(a, mod);
	b = num_mul(b, mod);
	range (k, n) {
		range (l, n) {
			out[k*n+l] = (num){0,0,1};
		}
	}
	range (k, n) {
		out[k*n+k] = (num){1, 0, 1};
	}
	out[j*n+j] = num_conj(a);
	out[i*n+i] = num_neg(a);
	out[i*n+j] = b;
	out[j*n+i] = num_conj(b);
}

int main() {
	num fourier[4][4];
	range(i, 4) {
		range (j, 4) {
			switch ((i*j)%4) {
				case 0:
					fourier[i][j] = (num) {1, 0, 4};
					break;
				case 1:
					fourier[i][j] = (num) {0, 1, 4};
					break;
				case 2:
					fourier[i][j] = (num) {-1, 0, 4};
					break;
				case 3:
					fourier[i][j] = (num) {0, -1, 4};
			}
		}
	}
	printf("fourier: \n");
	mat_print(4, &fourier[0][0], ", ", "\n");

	const int REMOVE = (4*(4-1))/2-1;
	num u[REMOVE][4][4];
	num fred[REMOVE][4][4];
	int j = 0;
	int i = j+1;
	range(ui, REMOVE) {
		num *fprev = ui > 0 ? &fred[ui-1][0][0] : &fourier[0][0];
		num *ucurr = &u[ui][0][0];
		num *fcurr = &fred[ui][0][0];

		eliminator(4, i, j, ucurr, fprev);
		mat_reduce(4, ucurr);
		printf("\nu[%lu]: \n", ui);
		mat_print(4, ucurr, ", ", "\n");

		mat_mul(4, fcurr, ucurr, fprev);
		mat_reduce(4, fcurr);
		printf("\nfred[%lu]: \n", ui);
		mat_print(4, fcurr, ", ", "\n");

		if (i == 3) {
			j++;
			i = j+1;
		} else {
			i++;
		}
	}
}

