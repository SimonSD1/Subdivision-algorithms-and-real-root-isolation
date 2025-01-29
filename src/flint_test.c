#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>

//https://flintlib.org/doc/fmpz_poly.html

int main() {
    fmpz_poly_t poly1, poly2, result;

    fmpz_poly_init(poly1);
    fmpz_poly_init(poly2);
    fmpz_poly_init(result);


    // on cree F et G en indiquant le nombre de coef puis les coefs
    // coef faible a gauche
    fmpz_poly_set_str(poly1, "3  1 2 3"); 
    fmpz_poly_set_str(poly2, "2  1 4");   

    printf("F : ");
    fmpz_poly_print(poly1);
    printf("\n");

    printf("G : ");
    fmpz_poly_print(poly2);
    printf("\n");
    fmpz_poly_mul(result, poly1, poly2);

    printf("Résultat F*G: ");
    fmpz_poly_print(result);
    printf("\n");

    //long int
    fmpz_t res,x;
    fmpz_init(res);
    fmpz_set_ui(x, 5);
    fmpz_poly_evaluate_fmpz(res,poly1,x);
    printf("resultat de l'evaluation:\n");
    fmpz_print(res);
    printf("\n");
    

    fmpz_poly_clear(poly1);
    fmpz_poly_clear(poly2);
    fmpz_poly_clear(result);

    return 0;
}
