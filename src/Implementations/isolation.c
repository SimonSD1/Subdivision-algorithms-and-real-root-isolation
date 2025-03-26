#include "../HeaderFiles/isolation.h"

// Voici ce que je vous propose d'implémenter :
// Une fonction qui prend en entrée :
//* Un polynôme pol en une variable
//* c un entier flint (fmpz_t)
//* k un entier machine (int32_t)
//* une tableau de solutions codées sous la forme (a, p) -> une solution est dans
//(a / 2^p, (a+1) / 2^p) => ici il faut définir les structures de données
// adéquates.
//* un pointeur sur un entier qui encode le nombre de racines trouvées.
//
// Cette fonction va isoler les racines de pol sur l'intervalle
//(c / 2^k, (c+1) / 2^k)
//
// MAIS on suppose que les changements de variables ont été faits sur pol
// pour qu'on soit ramené à regarder les racines dans l'intervalle (0, 1)
//
//- étape 1 : est-ce que 0 est racine -> à implémenter et mettre à jour ce qui
// doit être mis à jour en fonction du résultat (penser à diviser le polynôme par
// x)
//
//- étape 2 : est-ce que 1/2 est racine ?
//
//- étape 3 : calculer le nombre de variations de signes des coefficients de pol
//
//    (a) Si ce nombre est 0 que pouvez-vous en déduire ?
//    (b) Si ce nombre n'est pas 0 : appliquer la règle des signes de descartes à P
//    tel que P = Q(X+1) où Q est le numérateur de pol( 1/ x ) -> récupérer le nbre
//    de variations de signes (on l'appelle nb) => ça revient à appliquer Descartes
//    sur pol( 1 / (x+1))
//
//    => faire une fonction "naive" mais aussi penser à ce qu'on a dit sur la
//    croissance des coefficients quand on fait le Taylor shift
//
//    => qu'est-ce que nb dit sur le nbre de racines ?
//
//    -> distinguer le cas nb = 0
//    -> distinguer le cas nb = 1
//    -> distinguer le cas nb >= 2
//       => dans ce cas, il y a une bisection et il faut rappeler l'algo sur les
//       intervalles (0, 1/2) puis (1/2, 1) mais attention, dans les appels
//       récursifs il faut bien mettre à jour c et k.

void div_by_x(fmpz_poly_t pol)
{
    fmpz_poly_shift_right(pol, pol, 1);
}

int fmpz_poly_is_half_root(fmpz_poly_t pol)
{

    fmpz_t numerator;
    fmpz_init(numerator);
    fmpz_zero(numerator);

    fmpz_t temp;
    fmpz_init(temp);

    slong d = pol->length - 1;

    for (slong i = d; i >= 0; i--)
    {
        fmpz_poly_get_coeff_fmpz(temp, pol, i);
        fmpz_mul_2exp(temp, temp, d - i);
        fmpz_add(numerator, numerator, temp);
    }

    fmpz_clear(numerator);
    fmpz_clear(temp);

    if (fmpz_is_zero(numerator))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void var_change_to_inf(fmpz_poly_t result, fmpz_poly_t pol){
    fmpz_poly_init(result);


}

// Cette fonction va isoler les racines de pol sur l'intervalle
//(c / 2^k, (c+1) / 2^k)
void isolation(fmpz_poly_t pol, fmpz_t c, __int32_t k, solution *solutions, int *nb_sol)
{
    //- étape 1 : est-ce que 0 est racine -> à implémenter et mettre à jour ce qui
    // doit être mis à jour en fonction du résultat (penser à diviser le polynôme par
    // x)

    fmpz_t eval;
    fmpz_init(eval);

    evaluate_0(eval, pol);

    if (fmpz_is_zero(eval))
    {
        // on divise par x
        div_by_x(pol);

        // ajouter la racine a la liste des solutions
        *solutions[*nb_sol].c = c;
        *solutions[*nb_sol].k = k;
        *solutions[*nb_sol].h = 0;
        *nb_sol++;
    }

    // étape 3 : calculer le nombre de variations de signes des coefficients de pol sur [0;+inf]
    // on utilise un changement de variable
    int sign_changes = descartes_rule(pol);

    printf("\nsigne_changes = %d \n", sign_changes);

    if (sign_changes == 0)
    {
        // rien dans cet interval
        return;
    }

    else if (sign_changes == 1)
    {
        // il y a exactement une solution dans cet interval

        *solutions[*nb_sol].c = c;
        *solutions[*nb_sol].k = k;
        *nb_sol++;
        return;
    }

    else
    {
        // il y a au plus 2 solution dans cet interval, fait la bisection

        // on appel a gauche sur p(x/2)

        

        isolation()


        // on appel a droite sur p((x+1)/2)
    }
}
