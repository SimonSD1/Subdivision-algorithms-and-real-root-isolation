#include <flint/flint.h>
#include <flint/fmpz_poly.h>

int main() {
    fmpz_poly_t poly1, poly2, result;

    // Initialisation des polynômes
    fmpz_poly_init(poly1);
    fmpz_poly_init(poly2);
    fmpz_poly_init(result);

    // Définir des coefficients pour poly1 et poly2
    fmpz_poly_set_str(poly1, "3  1 2 3"); // polynôme 3x^2 + 2x + 1
    fmpz_poly_set_str(poly2, "2  1 4");   // polynôme 4x + 1

    // Produit des deux polynômes
    fmpz_poly_mul(result, poly1, poly2);

    // Affichage du résultat
    printf("Résultat : ");
    fmpz_poly_print(result);
    printf("\n");

    // Libération de la mémoire
    fmpz_poly_clear(poly1);
    fmpz_poly_clear(poly2);
    fmpz_poly_clear(result);

    return 0;
}
