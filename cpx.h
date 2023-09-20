#pragma once

#define CPX_API_ENTRY
#define CPX_API_CALL
#define CPX_API_SUFFIX__VERSION_1_0_0

#ifdef __cplusplus
extern "C" {
#else
#include <stdbool.h>
#endif

typedef struct _cpx_calculation_type *cpx_calculation_type;

extern CPX_API_ENTRY cpx_calculation_type CPX_API_CALL
cpx_newCalculation(void) CPX_API_SUFFIX__VERSION_1_0_0;

extern CPX_API_ENTRY void CPX_API_CALL
cpx_loadparam(char* /* Method identifier */,
              char* /* Solvent name */,
              cpx_calculation_type /* Calculation object */
              ) CPX_API_SUFFIX__VERSION_1_0_0;

extern CPX_API_ENTRY void CPX_API_CALL
cpx_readparam(char* /* CRS Parameter File */,
              char* /* SMD Parameter File */,
              cpx_calculation_type /* Calculation object */
              ) CPX_API_SUFFIX__VERSION_1_0_0;

extern CPX_API_ENTRY void CPX_API_CALL
cpx_loadsolvent(char* /* Solvent name */,
                cpx_calculation_type /* Calculation object */
                ) CPX_API_SUFFIX__VERSION_1_0_0;

extern CPX_API_ENTRY void CPX_API_CALL
cpx_loadsolute(char* /* Solute name */,
               cpx_calculation_type /* Calculation object */
               ) CPX_API_SUFFIX__VERSION_1_0_0;

extern CPX_API_ENTRY void CPX_API_CALL
cpx_calculate(cpx_calculation_type /* Calculation object */,
              char* /* Method identifier */,
              char* /* Solvent name */,
              double /* Gas Phase Energy*/,
              double /* Probe radius */,
              double /* Temperature */,
              double /* Convergence criterium */
              ) CPX_API_SUFFIX__VERSION_1_0_0;

extern CPX_API_ENTRY void CPX_API_CALL
cpx_readCosmoFile(cpx_calculation_type /* Calculation object */,
                  char* /* Solvent or Solute?*/,
                  char* /* Cosmo file name */
                  ) CPX_API_SUFFIX__VERSION_1_0_0;

extern CPX_API_ENTRY void CPX_API_CALL
cpx_getenergies(cpx_calculation_type /* Calculation object */,
                double* /* Energy array (is, cc, res, smd, ss, shift) */
                ) CPX_API_SUFFIX__VERSION_1_0_0;

extern CPX_API_ENTRY void CPX_API_CALL
cpx_deleteCalculation(cpx_calculation_type /* Calculation object */
                      ) CPX_API_SUFFIX__VERSION_1_0_0;