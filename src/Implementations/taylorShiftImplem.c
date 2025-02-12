#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/taylorShift_implem.h"

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


// the one from flint uses horner for polynomial with degree less than 50
// and parallelism
void poly_shift(fmpz_poly_t g, fmpz_poly_t poly, fmpz_t a)
{

    /// precomputation

    slong len = fmpz_poly_length(poly);

    printf("\n len = %ld \n",len);

    fmpz_print(a);
    printf("\n");

    slong m;
    fmpz_t len_fmpz;
    fmpz_init_set_si(len_fmpz,len);

    m=fmpz_clog_ui(len_fmpz,2);


    fmpz_poly_t *precomputed = flint_malloc(m * sizeof(fmpz_poly_t));
    if(!precomputed){
        printf("error allocating\n");
    }
    for (slong i = 0; i < m; i++)
    {
        fmpz_poly_init(precomputed[i]);
    }

    fmpz_poly_set_coeff_si(precomputed[0], 1, 1);
    fmpz_poly_set_coeff_fmpz(precomputed[0], 0, a);

    // square (x+a) m times by using the previous one
    for (slong i = 1; i < m; i++)
    {
        fmpz_poly_sqr(precomputed[i], precomputed[i - 1]);
    }

    // g = result of poly(x+a)
    divide_conquer(g, poly, precomputed, m);

    for (slong i = 0; i < m; i++)
        fmpz_poly_clear(precomputed[i]);
    flint_free(precomputed);
}

void naiveShift(fmpz_poly_t result, fmpz_poly_t poly, fmpz_t a){
    fmpz_poly_init(result);
    slong length=poly->length;

    fmpz_poly_t power;
    fmpz_poly_init(power);
    fmpz_poly_set_coeff_si(power, 0,1);

    fmpz_poly_t x_plus_a;
    fmpz_poly_init(x_plus_a);
    fmpz_poly_set_coeff_si(x_plus_a, 1, 1);
    fmpz_poly_set_coeff_fmpz(x_plus_a, 0, a);
    
    fmpz_t coeff;
    fmpz_init(coeff);

    fmpz_poly_t temp;
    fmpz_poly_init(temp);
    
    for(slong i=0; i<length; i++){
        fmpz_poly_get_coeff_fmpz(coeff,poly,i);

        // poly_i * (x+a)^i
        fmpz_poly_scalar_mul_fmpz(temp,power,coeff);

        fmpz_poly_add(result,result,temp);

        fmpz_poly_mul(power,power,x_plus_a);
    }
}

