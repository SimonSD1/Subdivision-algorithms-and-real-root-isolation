#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/bound.h"
#include "../HeaderFiles/functionsForTests.h"

#include <string.h>

int descartes_rule(const fmpz_poly_t poly) {
    int sign_changes = 0;
    int last_sign = 0;
    slong len = fmpz_poly_length(poly);

    fmpz_t coeff;
    fmpz_init(coeff);
    
    for (slong i = 0; i < len; i++) {
        fmpz_poly_get_coeff_fmpz(coeff,poly, i);
        if (!fmpz_is_zero(coeff) ){
            int sign = fmpz_sgn(coeff);
            if (last_sign != 0 && sign != last_sign) {
                sign_changes++;
            }
            last_sign = sign;
        }
    }
    return sign_changes;
}