#include <stdlib.h>
#include <math.h>

#include "../Matrix/src/matrix.h"


typedef struct {
    double bar_size;
    struct {
        double initial, final;
    } frontier_condition;
    double (*initial_condition) (double position);
} EDPHeatEquation_ST;

/* ========================================================================== */

double _initialEqualtion(double x) {
    if (x >= 0.0 && x <= 0.5) return 2.0 * x;
    else return 2.0 * (1.0 - x);
}

/* ========================================================================== */

static double _calculateElement (
    double r, double elem1, double elem2, double elem3
) {
    return r * elem1 + (1.0 - 2.0 * r) * elem2 + r * elem3;
}

/* ========================================================================== */

MTXMatrix edpGetEquationGraph (
    EDPHeatEquation_ST equation,
    double t_final, unsigned t_part, unsigned x_part
) {
    double k = t_final / (double) t_part;
    double h = equation.bar_size / (double) x_part;
    double i, j, elem1, elem2, elem3, r = k / (h * h);
    MTXMatrix points = mtxNew(t_part, x_part);

    for (i = 1; i < x_part - 1; i ++)
        mtxSetElement(points, 0, i, equation.initial_condition(i * h));
    for (i = 0; i < t_part; i ++) {
        mtxSetElement(points, i, 0, equation.frontier_condition.initial);
        mtxSetElement(
            points, i, x_part - 1, equation.frontier_condition.final
        );
    }

    for (i = 1; i < t_part; i ++)
        for (j = 1; j < x_part - 1; j ++) {
            elem1 = mtxGetElement(points, i - 1, j - 1);
            elem2 = mtxGetElement(points, i - 1, j);
            elem3 = mtxGetElement(points, i - 1, j + 1);
            mtxSetElement(
                points, i, j, _calculateElement(r, elem1, elem2, elem3)
            );
        }

    return points;
}

int main (void) {
    EDPHeatEquation_ST eq = { 1.0, { 0.0, 0.0 }, _initialEqualtion };
    MTXMatrix result = NULL;

    result = edpGetEquationGraph(eq, 0.003, 4, 10);

    mtxPrint(result);

    return 0;
}

