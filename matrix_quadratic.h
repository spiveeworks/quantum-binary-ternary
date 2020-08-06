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

// this type represents the field Q[sqrt(-1), sqrt(2), sqrt(3)]
// @Performance change all of the num_ procedures to use pointers? it got big
typedef struct num {
	int c[2][2][2]; // coefficients
	unsigned de; // denominator
} num;

const int MONOMIAL_SQUARES[2][2][2] = {
	{{1,3},{2,6}},
	{{-1,-3},{-2,-6}},
};

num num_from_int(int nu) {
	num out = {};
	out.c[0][0][0] = nu;
	out.de = 1;
	return out;
}

num num_reduce(num x) {
	int d = x.de;
	range(i, 2) {
		range(j, 2) {
			range(k, 2) {
				d = gcd(d, x.c[i][j][k]);
			}
		}
	}
	x.de /= d;
	range(i, 2) {
		range(j, 2) {
			range(k, 2) {
				x.c[i][j][k] /= d;
			}
		}
	}
	return x;
}

void num_print(num x) {
	int terms = 0;
	range(i, 2) {
		range(j, 2) {
			range(k, 2) {
				if (x.c[i][j][k] != 0) {
					terms += 1;
				}
			}
		}
	}
	//printf("{%d, %u, %d, %u, %u}", x.re.whole, x.re.qu, x.im.whole, x.im.qu, x.de);
	if (terms == 0) {
		printf("0");
	} else if (terms >= 2 && x.de != 1) {
		printf("(");
	}
	bool first_term = true;
	range(i, 2) {
		range(j, 2) {
			range(k, 2) {
				if (x.c[i][j][k] != 0) {
					if (x.c[i][j][k] > 0 && !first_term) {
						printf("+");
					}
					if (x.c[i][j][k] != 1 || (i == 0 && j == 0 && k == 0)) {
						printf("%d", x.c[i][j][k]);
					}
					if (j == 1 || k == 1) {
						printf("sqrt(%d)", MONOMIAL_SQUARES[0][j][k]);
					}
					if (i == 1) {
						printf("i");
					}
					first_term = false;
				}
			}
		}
	}
	if (terms >= 2 && x.de != 1) {
		printf(")");
		printf("/%u", x.de);
	}
}

num num_add(num x, num y) {
	num out;
	out.de = x.de * y.de;
	range(i, 2) {
		range(j, 2) {
			range(k, 2) {
				out.c[i][j][k] = x.c[i][j][k] * y.de + y.c[i][j][k] * x.de;
			}
		}
	}
	return out;
}

num num_mul(num x, num y) {
	int out_c[3][3][3];
	range(i1, 2) {
		range(j1, 2) {
			range(k1, 2) {
				range(i2, 2) {
					range(j2, 2) {
						range(k2, 2) {
							out_c[i1+i2][j1+j2][k1+k2] =
								x.c[i1][j1][k1] * y.c[i2][j2][k2];
						}
					}
				}
			}
		}
	}
	range(i, 3) {
		range(j, 3) {
			out_c[i][j][0] += out_c[i][j][2] * MONOMIAL_SQUARES[0][0][1];
			out_c[i][j][2] = 0;
			out_c[i][0][j] += out_c[i][2][j] * MONOMIAL_SQUARES[0][1][0];
			out_c[i][2][j] = 0;
			out_c[0][i][j] += out_c[2][i][j] * MONOMIAL_SQUARES[1][0][0];
			out_c[2][i][j] = 0;
		}
	}
	num out;
	range(i, 2) {
		range(j, 2) {
			range(k, 2) {
				out.c[i][j][k] = out_c[i][j][k];
			}
		}
	}
	out.de = x.de * y.de;
	return out;
}

num num_conj_i(num x) {
	range(j, 2) {
		range(k, 2) {
			x.c[1][j][k] = -x.c[1][j][k];
		}
	}
	return x;
}

num num_conj_2(num x) {
	range(i, 2) {
		range(k, 2) {
			x.c[i][1][k] = -x.c[i][1][k];
		}
	}
	return x;
}

num num_conj_3(num x) {
	range(i, 2) {
		range(j, 2) {
			x.c[i][j][1] = -x.c[i][j][1];
		}
	}
	return x;
}

num num_mod_sq(num x) {
	return num_mul(x, num_conj_i(x));
}

num num_inv(num x) {
	num out = {};
	out.c[0][0][0] = x.de;
	out.de = 1;
	x.de = 1;

	num conj_i = num_conj_i(x);
	x = num_mul(x, conj_i);
	out = num_mul(out, conj_i);

	num conj_2 = num_conj_2(x);
	x = num_mul(x, conj_2);
	out = num_mul(out, conj_2);

	num conj_3 = num_conj_3(x);
	x = num_mul(x, conj_3);
	out = num_mul(out, conj_3);

	out.de = x.c[0][0][0];
	return out;
}

num num_inv_sqrt(num x) {
	range(i, 2) {
		range(j, 2) {
			range(k, 2) {
				bool zero = x.c[i][j][k] == 0;
				bool is_const = i == 0 && j == 0 && k == 0;
				if (!zero && !is_const) {
					printf("Tried to take invsqrt of non-rational number\n");
					exit(1);
				}
			}
		}
	}
	int x_nu = x.c[0][0][0];

	int i = 0;
	if (x_nu < 0) {
		x_nu = -x_nu;
		i = 1;
	}

	int d = gcd(x_nu, x.de);
	x_nu /= d;
	x.de /= d;

	num out = {};
	out.de = int_quadpart(x_nu);
	x_nu /= out.de;
	int out_nu = int_quadpart(x.de);
	x.de /= out_nu;

	int j = 0;
	if (x_nu % 2 == 0) {
		j = 1;
		out.de *= 2;
		x_nu /= 2;
	}
	if (x.de % 2 == 0) {
		j = 1;
		x.de /= 2;
	}

	int k = 0;
	if (x_nu % 3 == 0) {
		k = 1;
		out.de *= 3;
		x_nu /= 3;
	}
	if (x.de % 3 == 0) {
		k = 1;
		x.de /= 3;
	}

	if (x.de != 1) {
		printf("Tried to take invsqrt of multiple of %d/%d\n", x_nu, x.de);
		exit(1);
	}
	if (x_nu != 1) {
		printf("Tried to take invsqrt of multiple of %d\n", x_nu);
		exit(1);
	}

	out.c[i][j][k] = out_nu;

	return out;
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

void mat_zero(int n, num *out) {
	range(i, n) {
		range(j, n) {
			out[i*n+j] = num_from_int(0);
		}
	}
}

// performs matrix multiplication
void mat_mul(int n, num *out, num *x, num *y) {
	range(i, n) {
		range (k, n) {
			num sum = num_from_int(0);
			range (j, n) {
				sum = num_add(sum, num_mul(x[i*n+j], y[j*n+k]));
			}
			out[i*n+k] = sum;
		}
	}
}

// takes the Kronecker product of an m x m matrix and an n x n matrix
// out must be an mn x mn matrix
void mat_kronecker(int m, int n, num *out, num *x, num *y) {
	range(i1, m) {
		range(i2, n) {
			range(j1, m) {
				range(j2, n) {
					out[(i1*n+i2)*m*n+(j1*n+j2)] =
						num_mul(x[i1*m+j1], y[i2*n+j2]);
				}
			}
		}
	}
}

void mat_ident(int n, num *out) {
	mat_zero(n, out);
	range(i, n) {
		out[i*n+i].c[0][0][0] = 1;
	}
}

// initializes the permutation matrix that maps basis vector i to i+1
void mat_shift(int n, num *out) {
	mat_zero(n, out);
	range(i, n) {
		out[((i+1)%n)*n+i].c[0][0][0] = 1;
	}
}

void mat_clock(int n, num *out, num *phase, size_t phase_len) {
	mat_zero(n, out);
	range(i, n) {
		out[i*n + i] = phase[i%phase_len];
	}
}

void mat_fourier(int n, num *out, num *phase, size_t phase_len) {
	num scale = num_inv_sqrt(num_from_int(n));
	range(i, n) {
		range(j, n) {
			num coord = phase[i*j%phase_len];
			out[i*n+j] = num_reduce(num_mul(scale, coord));
		}
	}
}

// generates a diagonal matrix which along with the previous 3 functions give a
// basis for the Clifford quotient group
void mat_cliff_diag(int n, num *out, num *phase, size_t phase_len) {
	mat_zero(n, out);
	range(i, n) {
		out[i*n+i] = phase[i*(n-i)%phase_len];
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

