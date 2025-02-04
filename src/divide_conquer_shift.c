#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>

char x = 'x';



// because we split in half, we know that f1 is factorizable by (x+a)^len/2
// so shifted_f = f1(x+a) * (x+a)^len/2 + f0(x+a)
// we can split the polynomials until we find f1 and f0 constant, so f1(x+a)=f1
void divide_conquer(fmpz_poly_t g, const fmpz_poly_t f, fmpz_poly_t *precomputed, slong k)
{
    if (k == 0)
    {
        fmpz_poly_set(g, f);
        return;
    }

    // f is size 2**(k) so the half is 2**(k-1)
    slong half = 1 << (k - 1);
    fmpz_poly_t f0, f1;

    fmpz_poly_init(f0);
    fmpz_poly_init(f1);

    // f0 = first half coefs
    fmpz_poly_set_trunc(f0, f, half);

    // f1 = last half coefs
    fmpz_poly_shift_right(f1, f, half);
    fmpz_poly_set_trunc(f1, f1, half);

    divide_conquer(f0, f0, precomputed, k - 1);
    divide_conquer(f1, f1, precomputed, k - 1);

    fmpz_poly_mul(f1, precomputed[k - 1], f1);
    fmpz_poly_add(g, f0, f1);

    fmpz_poly_clear(f0);
    fmpz_poly_clear(f1);
}

void poly_shift(fmpz_poly_t g, fmpz_poly_t poly, slong a)
{

    /// precomputation

    slong len = fmpz_poly_length(poly);
    printf("\n len=%ld\n", len);

    // Check if len is a power of two
    if ((len & (len - 1)) != 0)
    {
        flint_printf("Error: Polynomial length must be a power of two.\n");
        flint_abort();
    }

    slong m = 0;
    slong temp_len = len;
    // shift until temp_len = 0 to compute 2^m=len
    while (temp_len >>= 1)
        m++;

    printf("\n m=%ld\n", m);

    fmpz_poly_t *precomputed = flint_malloc(m * sizeof(fmpz_poly_t));
    for (slong i = 0; i < m; i++)
    {
        fmpz_poly_init(precomputed[i]);
    }

    fmpz_poly_set_coeff_si(precomputed[0], 1, 1);
    fmpz_poly_set_coeff_si(precomputed[0], 0, a);

    // square (x+a) m times by using the previous one
    for (slong i = 1; i < m; i++)
    {
        fmpz_poly_sqr(precomputed[i], precomputed[i - 1]);
        printf("\n");
        fmpz_poly_print_pretty(precomputed[i], &x);
        printf("\n");
    }

    // g = result of poly(x+a)
    divide_conquer(g, poly, precomputed, m);

    for (slong i = 0; i < m; i++)
        fmpz_poly_clear(precomputed[i]);
    flint_free(precomputed);
}

int main()
{

    // polynomials
    fmpz_poly_t P;
    fmpz_poly_init(P);

    // used for the randomness
    flint_rand_t state;

    // number of coefficients of the polynomial
    slong len;

    // number of bits per coefficient
    flint_bitcnt_t bits;

    len = 8;
    bits = 10;
    flint_randinit(state);
    // bug if time(NULL) used twice
    flint_randseed(state, time(NULL), 0);

    fmpz_poly_randtest(P, state, len, bits);
    while (P->length != len)
    {
        fmpz_poly_randtest(P, state, len, bits);
    }

    printf("\npoly length=%ld\n", P->length);

    printf("\nPolynôme aléatoire :\n");
    fmpz_poly_print(P);
    printf("\n");
    fmpz_poly_print_pretty(P, &x);
    printf("\n\n");

    slong a = 50;

    fmpz_t a_flint;
    fmpz_init_set_ui(a_flint, a);

    fmpz_poly_t result;
    fmpz_poly_init(result);

    fmpz_poly_t result_flint;
    fmpz_poly_init(result_flint);

    poly_shift(result, P, a);

    printf("\nResult :\n");
    fmpz_poly_print(result);
    printf("\n");
    fmpz_poly_print_pretty(result, &x);
    printf("\n\n");

    fmpz_poly_taylor_shift_horner(result_flint, P, a_flint);
    printf("\nResult :\n");
    fmpz_poly_print(result_flint);
    printf("\n");
    fmpz_poly_print_pretty(result_flint, &x);
    printf("\n\n");

    fmpz_poly_clear(result);
    fmpz_poly_clear(result_flint);
    fmpz_poly_clear(P);
    return 0;
}
