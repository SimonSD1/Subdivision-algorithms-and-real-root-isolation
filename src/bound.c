#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>


// bound the abs of the complex roots
void Cauchy_bound(fmpq_t bound, const fmpz_poly_t poly)
{
    fmpz_t max_coeff;
    fmpz_init(max_coeff);

    // find max of |a_i|
    for (slong i = 0; i < poly->length ; i++)
    {
        // cmpabs = compare abs
        if (fmpz_cmpabs(poly->coeffs + i, max_coeff) > 0)
        {
            fmpz_set(max_coeff, poly->coeffs + i);
        }
    }

    //  1 + max(|a_i|) / |a_n|
    fmpq_set_fmpz(bound, max_coeff);
    fmpq_div_fmpz(bound, bound, poly->coeffs + poly->length - 1);
    fmpq_add_si(bound, bound, 1);
    fmpq_abs(bound, bound);

    fmpz_clear(max_coeff);
}

// presented in https://arxiv.org/pdf/2008.11039 and https://faculty.e-ce.uth.gr/akritas/phd_thesis_vigklas.pdf
// bound positive roots
void local_max_bound_implementation(fmpz_t ub, fmpz_poly_t poly) {
    slong n = fmpz_poly_degree(poly);
    if (n < 1) {
        fmpz_zero(ub);
        return;
    }

    slong j = n;
    slong t = 1;

    fmpz_t coef_i, coef_j;
    fmpz_init(coef_i);
    fmpz_init(coef_j);

    // Check if the leading coefficient is negative; if so, negate the polynomial
    fmpz_poly_get_coeff_fmpz(coef_i, poly, n);
    if (fmpz_sgn(coef_i) < 0) {
        fmpz_poly_neg(poly, poly);
    }

    fmpz_t tempub, old_tempub, r, r_pow_k, K;
    fmpz_init(tempub);

    fmpz_zero(ub);

    for (slong i = n ; i >= 0; i--) {
        fmpz_poly_get_coeff_fmpz(coef_i, poly, i);
        fmpz_poly_get_coeff_fmpz(coef_j, poly, j);

        if (fmpz_sgn(coef_i) < 0) {
            // Compute target = 2^t * (-coef_i)
            fmpz_neg(tempub, coef_i);
            fmpz_mul_2exp(tempub, tempub, t);

            // ceil division
            fmpz_cdiv_q(tempub,tempub,coef_j);

            slong k = j - i;

            // tempub = integer part of k-root of tempub so we add one
            fmpz_root(tempub, tempub, k);
            fmpz_add_si(tempub,tempub,1);
            

            // Update ub if K is larger
            if (fmpz_cmp(tempub, ub) > 0) {
                fmpz_set(ub, tempub);
                t++;
            }
        } else {
            // Update j and t if current coefficient is larger than coef_j
            if (fmpz_cmp(coef_i, coef_j) > 0) {
                j = i;
                t = 1;
            }
        }
    }

    fmpz_clear(coef_i);
    fmpz_clear(coef_j);
    fmpz_clear(tempub);
}

int main()
{

    char x= 'x';

    // polynomials
    fmpz_poly_t P;

    // used for the randomness
    flint_rand_t state;

    // number of coefficients of the polynomial
    slong len;

    // number of bits per coefficient
    flint_bitcnt_t bits;

    len = 5;
    bits = 10;
    flint_randinit(state);
    // bug if time(NULL) used twice
    flint_randseed(state, time(NULL), 0);

    fmpz_poly_randtest(P, state, len, bits);

    printf("\npoly length=%ld\n",P->length);

    printf("\nPolynôme aléatoire :\n");
    fmpz_poly_print(P);
    printf("\n");
    fmpz_poly_print_pretty(P,&x);
    printf("\n\n");

    fmpz_t bound_flint;
    fmpz_init(bound_flint);

    fmpz_poly_bound_roots(bound_flint, P);

    printf("bound_flint= : ");
    fmpz_print(bound_flint);
    printf("\n");

    // fmpq_t bound;
    fmpq_t bound_complex;
    fmpq_init(bound_complex);
    Cauchy_bound(bound_complex, P);

    printf("bound= : ");
    fmpq_print(bound_complex);
    printf("\n");

    fmpz_t bound_pos;
    local_max_bound_implementation(bound_pos,P);
    printf("bound_pos= : ");
    fmpz_print(bound_pos);
    printf("\n");

    return 0;
}
