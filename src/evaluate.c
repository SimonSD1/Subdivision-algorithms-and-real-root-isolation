#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>

void evaluate_0(fmpz_t result, fmpz_poly_t poly)
{
    fmpz_poly_get_coeff_fmpz(result, poly, 0);
}

void evaluate_1(fmpz_t result, fmpz_poly_t poly)
{
    fmpz_init(result);
    fmpz_t temp;
    fmpz_init(temp);

    for (int i = 0; i < poly->length; i++)
    {
        fmpz_poly_get_coeff_fmpz(temp, poly, i);
        fmpz_add(result, result, temp);
    }
    fmpz_clear(temp);
}

void evaluate_half(fmpz_t result, fmpz_poly_t poly)
{
    fmpz_init(result);
    fmpz_t temp;
    fmpz_init(temp);

    for (slong i = 0; i < poly->length; i++)
    {
        fmpz_poly_get_coeff_fmpz(temp, poly, i);
        fmpz_fdiv_q_2exp(temp,temp, i);
        fmpz_add(result, result, temp);
    }

    fmpz_clear(temp);
}

int main(int argc, char const *argv[])
{

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

    printf("Polynôme aléatoire :\n");
    fmpz_poly_print(P);
    printf("\n");

    fmpz_t result;
    fmpz_init(result);
    fmpq_t resultq;
    fmpq_init(resultq);

    fmpq_t half;
    fmpq_init(half);
    fmpq_set_si(half,1,2);

    evaluate_0(result, P);
    printf("P(0)= ");
    fmpz_print(result);
    printf("\n");

    evaluate_1(result, P);
    printf("P(1)= ");
    fmpz_print(result);
    printf("\n");

    evaluate_half(result, P);
    printf("P(1/2)= ");
    fmpz_print(result);
    printf("\n");
    fmpz_poly_evaluate_fmpq(resultq,P,half);
    fmpq_print(resultq);
    printf("\n");

    fmpz_clear(result);

    return 0;
}
