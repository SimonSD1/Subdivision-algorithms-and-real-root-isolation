#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_vec.h>
#include <flint/ulong_extras.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>
#include "../HeaderFiles/taylorShift_implem.h"
#include  <flint/thread_pool.h>
#include  <flint/thread_support.h>


fmpz **global_precomputed = NULL;
slong *global_precomputed_len = NULL;
slong global_precomputed_size = 0;



void load_precomputed_polynomials(slong max_m) {
    global_precomputed_size = max_m;
    
    // Allocate array of fmpz* pointers and length storage
    global_precomputed = flint_malloc(max_m * sizeof(fmpz *));
    global_precomputed_len = flint_malloc(max_m * sizeof(slong));
    
    for (slong i = 0; i < max_m; i++) {
        fmpz_poly_t tmp;
        fmpz_poly_init(tmp);
        
        char pathFile[80];
        sprintf(pathFile, "DATA/Precomputation_DivConq/precomputation%d.txt", (int)i);
        FILE *file = fopen(pathFile, "r");
        if (!file || !fmpz_poly_fread(file, tmp)) {
            printf("Failed to load precomputation%d.txt\n", (int)i);
            exit(EXIT_FAILURE);
        }
        fclose(file);
        
        // Store coefficients as raw vector
        global_precomputed_len[i] = tmp->length;
        global_precomputed[i] = _fmpz_vec_init(tmp->length);
        _fmpz_vec_set(global_precomputed[i], tmp->coeffs, tmp->length);
        
        fmpz_poly_clear(tmp);
    }
}

void free_global_precomputed() {
    if (global_precomputed) {
        for (slong i = 0; i < global_precomputed_size; i++) {
            _fmpz_vec_clear(global_precomputed[i], global_precomputed_len[i]);
        }
        flint_free(global_precomputed);
        flint_free(global_precomputed_len);
    }
    global_precomputed = NULL;
    global_precomputed_len = NULL;
    global_precomputed_size = 0;
}




void divide_conquer_recursive(fmpz * poly, slong len) {
    fmpz *tmp, *tmp2;
    slong len1, len2;
    slong bits, cutoff;
    fmpz_t shift;

    if (len < 50) {
        fmpz_init_set_si(shift, 1);
        _fmpz_poly_taylor_shift_horner(poly, shift, len);
        fmpz_clear(shift);
        return;
    }

    bits = _fmpz_vec_max_bits(poly, len);
    bits = FLINT_ABS(bits);
    cutoff = 100 + 10 * n_sqrt(FLINT_MAX(bits - FLINT_BITS, 0));

    cutoff = FLINT_MIN(cutoff, 1000);
    if (len < cutoff) {
        fmpz_init_set_si(shift, 1);
        _fmpz_poly_taylor_shift_horner(poly, shift, len);
        fmpz_clear(shift);
        return;
    }

    fmpz_t len_fmpz;
    fmpz_init_set_si(len_fmpz,len);
    slong m = fmpz_clog_ui(len_fmpz,2);
    int precomp;
    //pour éviter déséquilibre
    if (1 << (m - 2) > (len - (1 << (m - 1)))) {
        len1 = 1 << (m - 2);
        precomp = 2;
    }
    else {
        len1 = 1 << (m - 1);
        precomp = 1;
    }
    len2 = len - len1;

    divide_conquer_recursive(poly, len1);
    divide_conquer_recursive(poly + len1, len2);


    tmp = _fmpz_vec_init(len1 + 1);
    tmp2 = _fmpz_vec_init(len);

    slong idx = m - precomp;
    _fmpz_vec_set(tmp, global_precomputed[idx], global_precomputed_len[idx]);


    _fmpz_poly_mul(tmp2, tmp, len1 + 1, poly + len1, len2);

    _fmpz_vec_add(poly, poly, tmp2, len1);
    _fmpz_vec_set(poly + len1, tmp2 + len1, len2);

    _fmpz_vec_clear(tmp, len1 + 1);
    _fmpz_vec_clear(tmp2, len);
}

void poly_shift_plus_one_Precomputed(fmpz_poly_t result, const fmpz_poly_t poly) {
    if (poly != result)
        fmpz_poly_set(result, poly);

    divide_conquer_recursive(result->coeffs, result->length);
}







void divide_conquer_recursive2(fmpz * poly, slong len) {
    fmpz *tmp, *tmp2;
    slong len1, len2;
    slong bits, cutoff;
    fmpz_t shift;

    if (len < 50) {
        fmpz_init_set_si(shift, 1);
        _fmpz_poly_taylor_shift_horner(poly, shift, len);
        fmpz_clear(shift);
        return;
    }

    bits = _fmpz_vec_max_bits(poly, len);
    bits = FLINT_ABS(bits);
    cutoff = 100 + 10 * n_sqrt(FLINT_MAX(bits - FLINT_BITS, 0));

    cutoff = FLINT_MIN(cutoff, 1000);
    if (len < cutoff) {
        fmpz_init_set_si(shift, 1);
        _fmpz_poly_taylor_shift_horner(poly, shift, len);
        fmpz_clear(shift);
        return;
    }

    fmpz_t len_fmpz;
    fmpz_init_set_si(len_fmpz,len);
    slong m = fmpz_clog_ui(len_fmpz,2);
    //pour éviter déséquilibre
    if (1 << (m - 2) > (len - (1 << (m - 1)))) {
        len1 = 1 << (m - 2);
    }
    else {
        len1 = 1 << (m - 1);
    }
    len2 = len - len1;

    divide_conquer_recursive2(poly, len1);
    divide_conquer_recursive2(poly + len1, len2);


    tmp = _fmpz_vec_init(len1 + 1);
    tmp2 = _fmpz_vec_init(len);

    
    /* Now we generate (x+c)^len1 using binomial expansion. It's redundant
    to do this in all branches of the tree, but since it's just O(d),
    it's going to be cheap compared to the actual multiplications
    anyway. */
    fmpz_one(tmp);
    for (slong k = 1; k <= len1; k++)
    {
        if (k > len1 - k)
        {
            fmpz_set(tmp + k, tmp + len1 - k);
        }
        else
        {
           fmpz_mul_ui(tmp + k, tmp + k - 1, len1 + 1 - k);
           fmpz_divexact_ui(tmp + k, tmp + k, k);
        }
    }


    _fmpz_poly_mul(tmp2, tmp, len1 + 1, poly + len1, len2);

    _fmpz_vec_add(poly, poly, tmp2, len1);
    _fmpz_vec_set(poly + len1, tmp2 + len1, len2);

    _fmpz_vec_clear(tmp, len1 + 1);
    _fmpz_vec_clear(tmp2, len);
}

void poly_shift_plus_one_Non_Precomputed(fmpz_poly_t result, const fmpz_poly_t poly) {
    if (poly != result)
        fmpz_poly_set(result, poly);

    divide_conquer_recursive2(result->coeffs, result->length);
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

