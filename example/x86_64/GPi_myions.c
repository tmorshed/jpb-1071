/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__GPi_myions
#define _nrn_initial _nrn_initial__GPi_myions
#define nrn_cur _nrn_cur__GPi_myions
#define _nrn_current _nrn_current__GPi_myions
#define nrn_jacob _nrn_jacob__GPi_myions
#define nrn_state _nrn_state__GPi_myions
#define _net_receive _net_receive__GPi_myions 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define nai0 _p[0]
#define nao0 _p[1]
#define ki0 _p[2]
#define ko0 _p[3]
#define nai _p[4]
#define Dnai _p[5]
#define nao _p[6]
#define Dnao _p[7]
#define ki _p[8]
#define Dki _p[9]
#define ko _p[10]
#define Dko _p[11]
#define v _p[12]
#define _ion_ko	*_ppvar[0]._pval
#define _ion_ki	*_ppvar[1]._pval
#define _style_k	*((int*)_ppvar[2]._pvoid)
#define _ion_nao	*_ppvar[3]._pval
#define _ion_nai	*_ppvar[4]._pval
#define _style_na	*((int*)_ppvar[5]._pvoid)
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_GPi_myions", _hoc_setdata,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "nai0_GPi_myions", "mM",
 "nao0_GPi_myions", "mM",
 "ki0_GPi_myions", "mM",
 "ko0_GPi_myions", "mM",
 0,0
};
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GPi_myions",
 "nai0_GPi_myions",
 "nao0_GPi_myions",
 "ki0_GPi_myions",
 "ko0_GPi_myions",
 0,
 0,
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	nai0 = 10;
 	nao0 = 150;
 	ki0 = 140;
 	ko0 = 6.24;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_check_conc_write(_prop, prop_ion, 0);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0]._pval = &prop_ion->param[2]; /* ko */
 	_ppvar[1]._pval = &prop_ion->param[1]; /* ki */
 	_ppvar[2]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for k */
 prop_ion = need_memb(_na_sym);
 nrn_check_conc_write(_prop, prop_ion, 0);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[3]._pval = &prop_ion->param[2]; /* nao */
 	_ppvar[4]._pval = &prop_ion->param[1]; /* nai */
 	_ppvar[5]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for na */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _GPi_myions_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("na", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,(void*)0, (void*)0, (void*)0, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "#k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "#na_ion");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GPi_myions /home/taha/git/jpb-1071/example/Mechanisms/GPi_myions.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "default extra/intra-cellular ion concentrations for GPe";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 2);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 1);
   nrn_update_ion_pointer(_na_sym, _ppvar, 3, 2);
   nrn_update_ion_pointer(_na_sym, _ppvar, 4, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
 {
   nai = nai0 ;
   nao = nao0 ;
   ki = ki0 ;
   ko = ko0 ;
   }

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ko = _ion_ko;
  ki = _ion_ki;
  nao = _ion_nao;
  nai = _ion_nai;
 initmodel(_p, _ppvar, _thread, _nt);
  _ion_ko = ko;
  _ion_ki = ki;
  nrn_wrote_conc(_k_sym, (&(_ion_ki)) - 1, _style_k);
  _ion_nao = nao;
  _ion_nai = nai;
  nrn_wrote_conc(_na_sym, (&(_ion_nai)) - 1, _style_na);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ko = _ion_ko;
  ki = _ion_ki;
  nao = _ion_nao;
  nai = _ion_nai;
  _ion_ko = ko;
  _ion_ki = ki;
  _ion_nao = nao;
  _ion_nai = nai;
}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/taha/git/jpb-1071/example/Mechanisms/GPi_myions.mod";
static const char* nmodl_file_text = 
  "TITLE default extra/intra-cellular ion concentrations for GPe\n"
  "\n"
  "COMMENT\n"
  " Used to set: nai, nao, ki, ko\n"
  "   so that NEURON's default values are not used\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX GPi_myions\n"
  "    USEION k WRITE ko, ki\n"
  "    USEION na WRITE nao, nai\n"
  "	RANGE nai0, nao0, ki0, ko0\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (molar)	= (1/liter)\n"
  "    (mM)	= (millimolar)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	nai0 = 10	(mM)\n"
  "	nao0 = 150	(mM)\n"
  "	ki0  = 140	(mM)\n"
  "	ko0  = 6.24	(mM)\n"
  "}\n"
  "	\n"
  "STATE {\n"
  "    nai (mM)\n"
  "    nao (mM)\n"
  "    ki (mM)\n"
  "    ko (mM)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    nai = nai0\n"
  "    nao = nao0\n"
  "    ki = ki0\n"
  "    ko = ko0\n"
  "}\n"
  ;
#endif
