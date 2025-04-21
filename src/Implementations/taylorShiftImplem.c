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
        fmpz *tmp = _fmpz_vec_init(i + 1);
        
        /* Now we generate (x+c)^len1 using binomial expansion. It's redundant
        to do this in all branches of the tree, but since it's just O(d),
        it's going to be cheap compared to the actual multiplications
        anyway. */
        fmpz_one(tmp);
        for (slong k = 1; k <= i; k++)
        {
            if (k > i - k)
            {
                fmpz_set(tmp + k, tmp + i - k);
            }
            else
            {
            fmpz_mul_ui(tmp + k, tmp + k - 1, i + 1 - k);
            fmpz_divexact_ui(tmp + k, tmp + k, k);
            }
        }
        
        // Store coefficients as raw vector
        global_precomputed_len[i] = i+1;
        global_precomputed[i] = _fmpz_vec_init(global_precomputed_len[i]);
        _fmpz_vec_set(global_precomputed[i], tmp, global_precomputed_len[i]);
        
        _fmpz_vec_clear(tmp, i + 1);
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

    len1 = len/2;
    len2 = len - len1;

    divide_conquer_recursive(poly, len1);
    divide_conquer_recursive(poly + len1, len2);


    tmp = _fmpz_vec_init(len1 + 1);
    tmp2 = _fmpz_vec_init(len);

    _fmpz_vec_set(tmp, global_precomputed[len1], global_precomputed_len[len1]);


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


double log_two(slong x)
{
    if (x <= 0)
    {
        printf("Error: log2 undefined for x <= 0\n");
        return -1.0;
    }

    double result = 0.0;
    while (x > 1)
    {
        x /= 2;
        result += 1.0;
    }
    return result;
}


void fmpz_compute_binom_poly(fmpz_poly_t binom_poly, slong n)
{
    fmpz_poly_one(binom_poly);
    for (slong k = 1; k <= n; k++)
    {
        if (k > n - k)
        {
            fmpz_set(binom_poly->coeffs + k, binom_poly->coeffs + n - k);
        }
        else
        {
            fmpz_mul_ui(binom_poly->coeffs + k, binom_poly->coeffs + k - 1, n + 1 - k);
            fmpz_divexact_ui(binom_poly->coeffs + k, binom_poly->coeffs + k, k);
        }
    }
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


slong compute_power_array(fmpz_poly_t** power_array, fmpz_poly_t poly, slong threshold){
  
    slong len = poly->length;
    slong block_len = len;

    // Reduce block length until it is <= threshold
    while (block_len > threshold)
        block_len /= 2;

    // Find the nearest power of 2
    int l = log_two(len / block_len);

    if(l==0){
        return;
    }

    (*power_array) = malloc(l * sizeof(fmpz_poly_t));


   
    for (int i = 0; i < l; i++)
        fmpz_poly_init2((*power_array)[i], (block_len + 1) * (1 << i));

    fmpz_compute_binom_poly((*power_array)[0], block_len);

    
    for (int i = 1; i < l; i++)
    {
        _fmpz_poly_mul(
            (*power_array)[i]->coeffs,
            (*power_array)[i - 1]->coeffs, (block_len + 1) * (1 << (i - 1)),
            (*power_array)[i - 1]->coeffs, (block_len + 1) * (1 << (i - 1)));
    }

    return l;
}


void iterative_taylor_shift_precompute(fmpz_poly_t result, const fmpz_poly_t poly, slong threshold, fmpz_poly_t* power_array)
{
    slong len = poly->length;
    slong block_len = len;

    // Reduce block length until it is <= threshold
    while (block_len > threshold)
        block_len /= 2;

    // Find the nearest power of 2
    slong l = log_two(len / block_len);
    slong nblocks = 1 << l;
    slong last_block_len = len - (nblocks - 1) * block_len;

    fmpz_t one;
    fmpz_init_set_ui(one, 1);

    // Base case: if polynomial is small enough, do direct Taylor shift
    if (len <= threshold)
    {
        fmpz_poly_taylor_shift_horner(result, poly, one);
        fmpz_clear(one);
        return ;
    }

    fmpz_poly_t temp_poly, mul_res;
    fmpz_poly_init2(temp_poly, len);
    fmpz_poly_init2(mul_res, len);

    temp_poly->length = len;
    _fmpz_vec_set(temp_poly->coeffs, poly->coeffs, len);

    // Apply Taylor shift (x+1) to each block individually
    for (slong i = 0; i < nblocks - 1; i++)
        _fmpz_poly_taylor_shift_horner(temp_poly->coeffs + (i * block_len), one, block_len);
    _fmpz_poly_taylor_shift_horner(temp_poly->coeffs + (nblocks - 1) * block_len, one, last_block_len);

    slong current_len = block_len;
    slong current_nblocks = nblocks;
    slong current_last_block_len = last_block_len;
    int iter = 0;

    while (iter < l)
    {
        fmpz *binom_coeffs = power_array[iter]->coeffs;
        slong n_pairs = current_nblocks / 2;

        for (slong i = 0; i < n_pairs - 1; i++)
        {
            fmpz *p0 = temp_poly->coeffs + (2 * i * current_len);
            fmpz *p1 = temp_poly->coeffs + ((2 * i + 1) * current_len);
            fmpz *res = temp_poly->coeffs + (i * 2 * current_len);

            _fmpz_poly_mul(mul_res->coeffs, p1, current_len, binom_coeffs, current_len + 1);
            _fmpz_poly_add(res, p0, current_len, mul_res->coeffs, 2 * current_len);
        }

        // Last pair might be unequal size
        fmpz *p0 = temp_poly->coeffs + (2 * (n_pairs - 1) * current_len);
        fmpz *p1 = temp_poly->coeffs + ((2 * (n_pairs - 1) + 1) * current_len);
        fmpz *res = temp_poly->coeffs + ((n_pairs - 1) * 2 * current_len);

        _fmpz_poly_mul(mul_res->coeffs, p1, current_last_block_len, binom_coeffs, current_len + 1);
        _fmpz_poly_add(res, p0, current_len, mul_res->coeffs, current_len + current_last_block_len);

        current_last_block_len += current_len;
        current_len *= 2;
        current_nblocks /= 2;
        iter++;
    }

    temp_poly->length = len;
    fmpz_poly_swap(result, temp_poly);

    fmpz_clear(one);
    fmpz_poly_clear(mul_res);
    fmpz_poly_clear(temp_poly);
}
