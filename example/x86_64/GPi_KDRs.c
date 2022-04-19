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
 
#define nrn_init _nrn_init__GPi_KDRs
#define _nrn_initial _nrn_initial__GPi_KDRs
#define nrn_cur _nrn_cur__GPi_KDRs
#define _nrn_current _nrn_current__GPi_KDRs
#define nrn_jacob _nrn_jacob__GPi_KDRs
#define nrn_state _nrn_state__GPi_KDRs
#define _net_receive _net_receive__GPi_KDRs 
#define _f_settables _f_settables__GPi_KDRs 
#define settables settables__GPi_KDRs 
#define states states__GPi_KDRs 
 
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
#define gk _p[0]
#define pslowd _p[1]
#define pslowi _p[2]
#define ek _p[3]
#define ki _p[4]
#define Dpslowd _p[5]
#define Dpslowi _p[6]
#define ik _p[7]
#define pinf_sd _p[8]
#define ptau_sd _p[9]
#define pinf_si _p[10]
#define ptau_si _p[11]
#define v _p[12]
#define _g _p[13]
#define _ion_ki	*_ppvar[0]._pval
#define _ion_ek	*_ppvar[1]._pval
#define _ion_ik	*_ppvar[2]._pval
#define _ion_dikdv	*_ppvar[3]._pval
 
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
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_settables(void);
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
 "setdata_GPi_KDRs", _hoc_setdata,
 "settables_GPi_KDRs", _hoc_settables,
 0, 0
};
 
static void _check_settables(double*, Datum*, Datum*, NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, int _type) {
   _check_settables(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define usetable usetable_GPi_KDRs
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_GPi_KDRs", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gk_GPi_KDRs", "mho/cm2",
 0,0
};
 static double delta_t = 1;
 static double pslowi0 = 0;
 static double pslowd0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_GPi_KDRs", &usetable_GPi_KDRs,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GPi_KDRs",
 "gk_GPi_KDRs",
 0,
 0,
 "pslowd_GPi_KDRs",
 "pslowi_GPi_KDRs",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gk = 0.003;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 1, 1);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* ki */
 	_ppvar[1]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _GPi_KDRs_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GPi_KDRs /home/taha/git/jpb-1071/example/Mechanisms/GPi_KDRs.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_pinf_sd;
 static double *_t_ptau_sd;
 static double *_t_pinf_si;
 static double *_t_ptau_si;
static int _reset;
static char *modelname = "potassium delayed rectifier Kv2.1 membrane channel for GPi neuron model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_settables(_threadargsprotocomma_ double);
static int settables(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_settables(_threadargsprotocomma_ double _lv);
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   settables ( _threadargscomma_ v ) ;
   Dpslowd = ( pinf_sd - pslowd ) / ptau_sd ;
   Dpslowi = ( pinf_si - pslowi ) / ptau_si ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 settables ( _threadargscomma_ v ) ;
 Dpslowd = Dpslowd  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ptau_sd )) ;
 Dpslowi = Dpslowi  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ptau_si )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   settables ( _threadargscomma_ v ) ;
    pslowd = pslowd + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ptau_sd)))*(- ( ( ( pinf_sd ) ) / ptau_sd ) / ( ( ( ( - 1.0 ) ) ) / ptau_sd ) - pslowd) ;
    pslowi = pslowi + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ptau_si)))*(- ( ( ( pinf_si ) ) / ptau_si ) / ( ( ( ( - 1.0 ) ) ) / ptau_si ) - pslowi) ;
   }
  return 0;
}
 static double _mfac_settables, _tmin_settables;
  static void _check_settables(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_settables =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_settables)/400.; _mfac_settables = 1./_dx;
   for (_i=0, _x=_tmin_settables; _i < 401; _x += _dx, _i++) {
    _f_settables(_p, _ppvar, _thread, _nt, _x);
    _t_pinf_sd[_i] = pinf_sd;
    _t_ptau_sd[_i] = ptau_sd;
    _t_pinf_si[_i] = pinf_si;
    _t_ptau_si[_i] = ptau_si;
   }
   _sav_celsius = celsius;
  }
 }

 static int settables(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv) { 
#if 0
_check_settables(_p, _ppvar, _thread, _nt);
#endif
 _n_settables(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_settables(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_settables(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_settables * (_lv - _tmin_settables);
 if (isnan(_xi)) {
  pinf_sd = _xi;
  ptau_sd = _xi;
  pinf_si = _xi;
  ptau_si = _xi;
  return;
 }
 if (_xi <= 0.) {
 pinf_sd = _t_pinf_sd[0];
 ptau_sd = _t_ptau_sd[0];
 pinf_si = _t_pinf_si[0];
 ptau_si = _t_ptau_si[0];
 return; }
 if (_xi >= 400.) {
 pinf_sd = _t_pinf_sd[400];
 ptau_sd = _t_ptau_sd[400];
 pinf_si = _t_pinf_si[400];
 ptau_si = _t_ptau_si[400];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 pinf_sd = _t_pinf_sd[_i] + _theta*(_t_pinf_sd[_i+1] - _t_pinf_sd[_i]);
 ptau_sd = _t_ptau_sd[_i] + _theta*(_t_ptau_sd[_i+1] - _t_ptau_sd[_i]);
 pinf_si = _t_pinf_si[_i] + _theta*(_t_pinf_si[_i+1] - _t_pinf_si[_i]);
 ptau_si = _t_ptau_si[_i] + _theta*(_t_ptau_si[_i+1] - _t_ptau_si[_i]);
 }

 
static int  _f_settables ( _threadargsprotocomma_ double _lv ) {
   double _lrate_k ;
 _lrate_k = pow( 4.0 , ( ( 36.0 - 23.0 ) / 10.0 ) ) ;
   pinf_sd = 1.0 / ( 1.0 + exp ( - ( _lv + 21.98 ) / 8.245 ) ) ;
   ptau_sd = _lrate_k * 5.0 / ( exp ( ( _lv + 30.0 ) / 14.0 ) + exp ( - ( _lv + 30.0 ) / 20.0 ) ) + 5.0 ;
   pinf_si = 0.74 / ( 1.0 + exp ( ( _lv + 25.0 ) / 12.0 ) ) + 0.26 ;
   ptau_si = _lrate_k * 500.0 / ( exp ( ( _lv - 50.0 ) / 20.0 ) + exp ( ( 50.0 - _lv ) / 20.0 ) ) + 800.0 ;
    return 0; }
 
static void _hoc_settables(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 settables ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ki = _ion_ki;
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ki = _ion_ki;
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  pslowi = pslowi0;
  pslowd = pslowd0;
 {
   settables ( _threadargscomma_ v ) ;
   pslowd = pinf_sd ;
   pslowi = pinf_si ;
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

#if 0
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
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
  ki = _ion_ki;
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ik = ( gk ) * pslowd * pslowi * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
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
  ki = _ion_ki;
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
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
  ki = _ion_ki;
  ek = _ion_ek;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(pslowd) - _p;  _dlist1[0] = &(Dpslowd) - _p;
 _slist1[1] = &(pslowi) - _p;  _dlist1[1] = &(Dpslowi) - _p;
   _t_pinf_sd = makevector(401*sizeof(double));
   _t_ptau_sd = makevector(401*sizeof(double));
   _t_pinf_si = makevector(401*sizeof(double));
   _t_ptau_si = makevector(401*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/taha/git/jpb-1071/example/Mechanisms/GPi_KDRs.mod";
static const char* nmodl_file_text = 
  "TITLE potassium delayed rectifier Kv2.1 membrane channel for GPi neuron model\n"
  "\n"
  "COMMENT\n"
  " Potassium Kv2.1 membrane channel for GPi exhibiting *slow* deactivating\n"
  " components.  Based on derived kinetics from Baranauskas 1999.  Their\n"
  " primary experiments were performed at 35degC.\n"
  "\n"
  " Baranauskas 1999 -- Say in paper: Vh=-18mV, Vc=7mV\n"
  "    Based on my own fitting of their data: Vh=-21.98mV, Vc=8.245mV \n"
  " Hernandez 1999 -- tau at -40mV is 28.80ms\n"
  "    Scaled the Gillies2006 curve to match this point\n"
  "\n"
  " Disregard Q10 given the proximity of the experiments to 36degC\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "UNITS {\n"
  "    (mV) = (millivolt)\n"
  "    (mA) = (milliamp)\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX GPi_KDRs\n"
  "    USEION k READ ki,ek WRITE ik\n"
  "    RANGE gk\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    v             (mV)\n"
  "    dt            (ms)\n"
  "    gk = 0.0030   (mho/cm2)  : Baranauskas 1999\n"
  "    ek\n"
  "    ki\n"
  "    celsius\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    pslowd\n"
  "    pslowi\n"
  "}\n"
  "\n"
  "ASSIGNED { \n"
  "    ik (mA/cm2)\n"
  "    pinf_sd\n"
  "    ptau_sd (ms)\n"
  "    pinf_si\n"
  "    ptau_si (ms)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE states METHOD cnexp\n"
  "    ik = (gk)*pslowd*pslowi*(v-ek)\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  "    settables(v)\n"
  "    pslowd = pinf_sd\n"
  "    pslowi = pinf_si\n"
  "}\n"
  "\n"
  "DERIVATIVE states {  \n"
  "    settables(v)\n"
  "    pslowd' = (pinf_sd-pslowd)/ptau_sd\n"
  "    pslowi' = (pinf_si-pslowi)/ptau_si\n"
  "}\n"
  "\n"
  "PROCEDURE settables(v) {\n"
  "	\n"
  "    LOCAL rate_k\n"
  "    TABLE pinf_sd, ptau_sd, pinf_si, ptau_si DEPEND celsius FROM -100 TO 100 WITH 400\n"
  "\n"
  "    rate_k =  4^((36-23)/10)\n"
  "    pinf_sd = 1/(1+exp(-(v+21.98)/8.245))                        : slow deactivating component (parameters from Baranauskus1999)\n"
  "    :pinf_sd = 1/(1+exp((17.5-v)/10))                              : slow deactivating component (from Mohapatra/Held modelDB)\n"
  "    ptau_sd = rate_k*5/(exp((v+30)/14) + exp(-(v+30)/20)) + 5     : slow deactivating component (from Mohapatra/Held modelDB)\n"
  "\n"
  "    pinf_si = 0.74/(1 + exp((v+25)/12)) + 0.26                    : slow inactivating component (from Mohapatra/Held modelDB)\n"
  "    ptau_si = rate_k*500/(exp((v-50)/20) + exp((50-v)/20))+800    : slow inactivating component (from Mohapatra/Held modelDB)\n"
  "\n"
  "}\n"
  "\n"
  "UNITSON\n"
  "\n"
  ;
#endif
