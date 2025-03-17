#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>
#include "../HeaderFiles/taylorShift_implem.h"


fmpz_poly_t *global_precomputed = NULL;
slong global_precomputed_size = 0;

// because we split in half, we know that f1 is factorizable by (x+a)^len/2
// so shifted_f = f1(x+a) * (x+a)^len/2 + f0(x+a)
// we can split the polynomials until we find f1 and f0 constant, so f1(x+a)=f1
void divide_conquer_plus_one(fmpz_poly_t g, const fmpz_poly_t f, fmpz_poly_t *precomputed, slong k, fmpz_t a, slong cut)
{

    if(f->length<cut){
        fmpz_poly_taylor_shift_horner(g,f,a);
        return;
    }

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

    divide_conquer_plus_one(f0, f0, precomputed, k - 1,a,cut);
    divide_conquer_plus_one(f1, f1, precomputed, k - 1,a,cut);

    fmpz_poly_mul(f1, precomputed[k - 1], f1);
    fmpz_poly_add(g, f0, f1);

    fmpz_poly_clear(f0);
    fmpz_poly_clear(f1);
}


// the one from flint uses horner for polynomial with degree less than 50
// and parallelism
void poly_shift_plus_one(fmpz_poly_t g, fmpz_poly_t poly, fmpz_t a, slong cut)
{
    slong len = fmpz_poly_length(poly);
    //printf("len=%d\n", (int) len);

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
        // square the previous one in the i-th
        fmpz_poly_sqr(precomputed[i], precomputed[i - 1]);
    }

    // g = result of poly(x+a)
    divide_conquer_plus_one(g, poly, precomputed, m,a, cut);

    for (slong i = 0; i < m; i++)
        fmpz_poly_clear(precomputed[i]);
    flint_free(precomputed);
}





void poly_shift_plus_one_Precomputed(fmpz_poly_t g, fmpz_poly_t poly, slong cut)
{
    if(poly->length<cut){
        fmpz_t shift;
        fmpz_init(shift);
        fmpz_set_si(shift, 1);
        fmpz_poly_taylor_shift_horner(g,poly,shift);
        return;
    }
    //printf("len=%d\n", (int) poly->length);
    fmpz_t shift;
    fmpz_init(shift);
    fmpz_set_si(shift, 1);
    slong len = fmpz_poly_length(poly);


    slong m;
    fmpz_t len_fmpz;
    fmpz_init_set_si(len_fmpz,len);

    m=fmpz_clog_ui(len_fmpz,2);

    if (m == 0)
    {
        fmpz_poly_set(g, poly);
        return;
    }


    FILE* filePoly;
    char pathFile[80];
    fmpz_poly_t precomputed;
    fmpz_poly_init(precomputed);


    sprintf(pathFile, "../DATA/Precomputation_DivConq/precomputation%d.txt", (int) (m-1));
    filePoly = fopen(pathFile, "r");

    if (filePoly == NULL) {
        printf("The file is not opened. The program will "
            "now exit.\nm-1 = %d\n", (int) (m-1));
        exit(0);
    }

    if(!fmpz_poly_fread(filePoly, precomputed))
        printf("Could not read the FILE");

    fclose(filePoly);
    // f is size 2**(k) so the half is 2**(k-1)
    slong half = 1 << (m - 1);
    fmpz_poly_t f0, f1;

    fmpz_poly_init(f0);
    fmpz_poly_init(f1);

    // f0 = first half coefs
    fmpz_poly_set_trunc(f0, poly, half);

    // f1 = last half coefs
    fmpz_poly_shift_right(f1, poly, half);
    fmpz_poly_set_trunc(f1, f1, half);

    // g = result of poly(x+a)
    poly_shift_plus_one_Precomputed(f0,f0,cut);
    poly_shift_plus_one_Precomputed(f1,f1,cut);

    fmpz_poly_mul(f1, precomputed, f1);
    fmpz_poly_add(g, f0, f1);

    fmpz_poly_clear(precomputed);
    fmpz_clear(len_fmpz);
    fmpz_clear(shift);
}






void load_precomputed_polynomials(slong max_m) {
    global_precomputed_size = max_m;
    global_precomputed = flint_malloc(max_m * sizeof(fmpz_poly_t));
    
    for (slong i = 0; i < max_m; i++) {
        fmpz_poly_init(global_precomputed[i]);
        char pathFile[80];
        sprintf(pathFile, "../DATA/Precomputation_DivConq/precomputation%d.txt", (int)i);
        FILE *file = fopen(pathFile, "r");
        if (!file || !fmpz_poly_fread(file, global_precomputed[i])) {
            printf("Failed to load precomputation%d.txt\n", (int)i);
            exit(EXIT_FAILURE);
        }
        fclose(file);
    }
}

void poly_shift_plus_one_Precomputed2(fmpz_poly_t g, const fmpz_poly_t poly, slong cut, slong current_m) {
    if (poly->length < cut) {
        fmpz_t shift;
        fmpz_init(shift);
        fmpz_set_si(shift, 1);
        fmpz_poly_taylor_shift_horner(g, poly, shift);
        fmpz_clear(shift);
        return;
    }

    if (current_m == 0) {
        fmpz_poly_set(g, poly);
        return;
    }

    slong half = 1 << (current_m - 1);
    fmpz_poly_t f0, f1;
    fmpz_poly_init(f0);
    fmpz_poly_init(f1);

    fmpz_poly_set_trunc(f0, poly, half);
    fmpz_poly_shift_right(f1, poly, half);
    fmpz_poly_set_trunc(f1, f1, half);

    poly_shift_plus_one_Precomputed2(f0, f0, cut, current_m - 1);
    poly_shift_plus_one_Precomputed2(f1, f1, cut, current_m - 1);

    fmpz_poly_mul(f1, global_precomputed[current_m - 1], f1);
    fmpz_poly_add(g, f0, f1);

    fmpz_poly_clear(f0);
    fmpz_poly_clear(f1);
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

