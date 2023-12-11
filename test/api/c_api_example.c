#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "cpx.h"

#define check(x, ...)                                                          \
  _Generic((x), int : check_int, double : check_double)(x, __VA_ARGS__)

static inline bool check_int(int actual, int expected, const char *msg) {
  if (expected == actual) {
    return true;
  }
  fprintf(stderr, "FAIL: %s: expected %d, got %d\n", msg, expected, actual);
  return false;
}

static inline bool check_double(double actual, double expected, double tol,
                                const char *msg) {
  if (fabs(expected - actual) < tol) {
    return true;
  }
  fprintf(stderr, "FAIL: %s: expected %g, got %g\n", msg, expected, actual);
  return false;
}

int testFirst() {

    double energies[6];
    cpx_calculation_type calc = NULL;

    calc = cpx_newCalculation();

    cpx_loadparam("xtb","water",calc);
    cpx_loadsolvent("water",calc);
    cpx_loadsolute("hexadecane",calc);
    cpx_calculate(calc,"crs","water", -51.635659985185,0.3,298.15,1E-4);
    cpx_getenergies(calc,energies);
    cpx_deleteCalculation(calc);

    double sum1 = energies[0] + energies[1] + energies[2] + energies[3] + energies[4] + energies[5];

    if (!check(sum1, -1.718436e-03, 1.0e-5, "Energy for Hexadecane in Water does not match"))
        return 1;

    return 0;
}

int testSecond() {

    double energies[6];
    cpx_calculation_type calc = NULL;

    calc = cpx_newCalculation();

    char crsparam[200];
    strcpy(crsparam,getenv("CPXHOME"));
    strcat(crsparam,"/DB/xtb/crs.param_h2o");
    char smdparam[200];
    strcpy(smdparam,getenv("CPXHOME"));
    strcat(smdparam,"/DB/xtb/smd_h2o");

    cpx_readparam(crsparam,smdparam,calc);
    cpx_loadsolvent("water",calc);
    cpx_loadsolute("hexadecane",calc);
    cpx_calculate(calc,"crs","water", -51.635659985185,0.3,298.15,1E-4);
    cpx_getenergies(calc,energies);
    cpx_deleteCalculation(calc);

    double sum1 = energies[0] + energies[1] + energies[2] + energies[3] + energies[4] + energies[5];

    if (!check(sum1, -1.718436e-03, 1.0e-5, "Energy for Hexadecane in Water does not match"))
        return 1;

    return 0;
}

int testThird() {

    double energies[6];
    cpx_calculation_type calc = NULL;

    calc = cpx_newCalculation();

    cpx_loadparam("xtb","octanol",calc);
    cpx_readCosmoFile(calc,"Solvent","DB/xtb/octanol.cosmo");
    cpx_readCosmoFile(calc,"Solute","DB/xtb/hexadecane.cosmo");
    cpx_calculate(calc,"crs","octanol", -51.635659985185,0.3,298.15,1E-4);
    cpx_getenergies(calc,energies);
    cpx_deleteCalculation(calc);

    double sum1 = energies[0] + energies[1] + energies[2] + energies[3] + energies[4] + energies[5];

    printf("%e\n",sum1);
    if (!check(sum1, -1.532934E-2, 1.0e-5, "Energy for Hexadecane in Octanol does not match"))
        return 1;

    return 0;
}
int main(int argc, char **argv) {
  int stat = 0;
  stat += testFirst();
  stat += testSecond();
  stat += testThird();
  return stat > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
