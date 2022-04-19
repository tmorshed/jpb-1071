#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _ampa_reg(void);
extern void _axnode75_reg(void);
extern void _gabaa_reg(void);
extern void _GPi_cacum_reg(void);
extern void _GPi_HCN_reg(void);
extern void _GPi_HVA_reg(void);
extern void _GPi_KDRf_reg(void);
extern void _GPi_KDRs_reg(void);
extern void _gpi_reg(void);
extern void _GPi_myions_reg(void);
extern void _GPi_NaL_reg(void);
extern void _GPi_Na_reg(void);
extern void _GPi_SK_reg(void);
extern void _parak75_reg(void);
extern void _STh_reg(void);
extern void _STN_Cacum_reg(void);
extern void _STN_CaT_reg(void);
extern void _STN_HVA_reg(void);
extern void _STN_Ih_reg(void);
extern void _STN_KDR_reg(void);
extern void _STN_Kv31_reg(void);
extern void _STN_myions_reg(void);
extern void _STN_NaL_reg(void);
extern void _STN_Na_reg(void);
extern void _STN_sKCa_reg(void);
extern void _tmgsyn_reg(void);
extern void _train_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"Mechanisms//ampa.mod\"");
    fprintf(stderr," \"Mechanisms//axnode75.mod\"");
    fprintf(stderr," \"Mechanisms//gabaa.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_cacum.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_HCN.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_HVA.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_KDRf.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_KDRs.mod\"");
    fprintf(stderr," \"Mechanisms//gpi.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_myions.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_NaL.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_Na.mod\"");
    fprintf(stderr," \"Mechanisms//GPi_SK.mod\"");
    fprintf(stderr," \"Mechanisms//parak75.mod\"");
    fprintf(stderr," \"Mechanisms//STh.mod\"");
    fprintf(stderr," \"Mechanisms//STN_Cacum.mod\"");
    fprintf(stderr," \"Mechanisms//STN_CaT.mod\"");
    fprintf(stderr," \"Mechanisms//STN_HVA.mod\"");
    fprintf(stderr," \"Mechanisms//STN_Ih.mod\"");
    fprintf(stderr," \"Mechanisms//STN_KDR.mod\"");
    fprintf(stderr," \"Mechanisms//STN_Kv31.mod\"");
    fprintf(stderr," \"Mechanisms//STN_myions.mod\"");
    fprintf(stderr," \"Mechanisms//STN_NaL.mod\"");
    fprintf(stderr," \"Mechanisms//STN_Na.mod\"");
    fprintf(stderr," \"Mechanisms//STN_sKCa.mod\"");
    fprintf(stderr," \"Mechanisms//tmgsyn.mod\"");
    fprintf(stderr," \"Mechanisms//train.mod\"");
    fprintf(stderr, "\n");
  }
  _ampa_reg();
  _axnode75_reg();
  _gabaa_reg();
  _GPi_cacum_reg();
  _GPi_HCN_reg();
  _GPi_HVA_reg();
  _GPi_KDRf_reg();
  _GPi_KDRs_reg();
  _gpi_reg();
  _GPi_myions_reg();
  _GPi_NaL_reg();
  _GPi_Na_reg();
  _GPi_SK_reg();
  _parak75_reg();
  _STh_reg();
  _STN_Cacum_reg();
  _STN_CaT_reg();
  _STN_HVA_reg();
  _STN_Ih_reg();
  _STN_KDR_reg();
  _STN_Kv31_reg();
  _STN_myions_reg();
  _STN_NaL_reg();
  _STN_Na_reg();
  _STN_sKCa_reg();
  _tmgsyn_reg();
  _train_reg();
}

#if defined(__cplusplus)
}
#endif
