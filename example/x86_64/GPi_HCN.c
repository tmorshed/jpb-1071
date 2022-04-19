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
 
#define nrn_init _nrn_init__GPi_Ih
#define _nrn_initial _nrn_initial__GPi_Ih
#define nrn_cur _nrn_cur__GPi_Ih
#define _nrn_current _nrn_current__GPi_Ih
#define nrn_jacob _nrn_jacob__GPi_Ih
#define nrn_state _nrn_state__GPi_Ih
#define _net_receive _net_receive__GPi_Ih 
#define _f_setinf _f_setinf__GPi_Ih 
#define integrate integrate__GPi_Ih 
#define setinf setinf__GPi_Ih 
 
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
#define gh _p[0]
#define ih _p[1]
#define f _p[2]
#define Df _p[3]
#define finf _p[4]
#define ftau _p[5]
#define v _p[6]
#define _g _p[7]
 
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
 static void _hoc_setinf(void);
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
 "setdata_GPi_Ih", _hoc_setdata,
 "setinf_GPi_Ih", _hoc_setinf,
 0, 0
};
 
static void _check_setinf(double*, Datum*, Datum*, NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, int _type) {
   _check_setinf(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define eih eih_GPi_Ih
 double eih = -56.11;
#define usetable usetable_GPi_Ih
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_GPi_Ih", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "eih_GPi_Ih", "mV",
 "gh_GPi_Ih", "mho/cm2",
 "ih_GPi_Ih", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double f0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "eih_GPi_Ih", &eih_GPi_Ih,
 "usetable_GPi_Ih", &usetable_GPi_Ih,
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
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GPi_Ih",
 "gh_GPi_Ih",
 0,
 "ih_GPi_Ih",
 0,
 "f_GPi_Ih",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 8, _prop);
 	/*initialize range parameters*/
 	gh = 0.001;
 	_prop->param = _p;
 	_prop->param_size = 8;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _GPi_HCN_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 8, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GPi_Ih /home/taha/git/jpb-1071/example/Mechanisms/GPi_HCN.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_finf;
 static double *_t_ftau;
static int _reset;
static char *modelname = "Potassium Ih channel for GPi neuron model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_setinf(_threadargsprotocomma_ double);
static int setinf(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_setinf(_threadargsprotocomma_ double _lv);
 static int _slist1[1], _dlist1[1];
 static int integrate(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   setinf ( _threadargscomma_ v ) ;
   Df = ( finf - f ) / ftau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 setinf ( _threadargscomma_ v ) ;
 Df = Df  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ftau )) ;
  return 0;
}
 /*END CVODE*/
 static int integrate (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   setinf ( _threadargscomma_ v ) ;
    f = f + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ftau)))*(- ( ( ( finf ) ) / ftau ) / ( ( ( ( - 1.0 ) ) ) / ftau ) - f) ;
   }
  return 0;
}
 static double _mfac_setinf, _tmin_setinf;
  static void _check_setinf(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_setinf =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_setinf)/400.; _mfac_setinf = 1./_dx;
   for (_i=0, _x=_tmin_setinf; _i < 401; _x += _dx, _i++) {
    _f_setinf(_p, _ppvar, _thread, _nt, _x);
    _t_finf[_i] = finf;
    _t_ftau[_i] = ftau;
   }
   _sav_celsius = celsius;
  }
 }

 static int setinf(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv) { 
#if 0
_check_setinf(_p, _ppvar, _thread, _nt);
#endif
 _n_setinf(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_setinf(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_setinf(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_setinf * (_lv - _tmin_setinf);
 if (isnan(_xi)) {
  finf = _xi;
  ftau = _xi;
  return;
 }
 if (_xi <= 0.) {
 finf = _t_finf[0];
 ftau = _t_ftau[0];
 return; }
 if (_xi >= 400.) {
 finf = _t_finf[400];
 ftau = _t_ftau[400];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 finf = _t_finf[_i] + _theta*(_t_finf[_i+1] - _t_finf[_i]);
 ftau = _t_ftau[_i] + _theta*(_t_ftau[_i+1] - _t_ftau[_i]);
 }

 
static int  _f_setinf ( _threadargsprotocomma_ double _lv ) {
   double _lvhalf , _lsfactor ;
 _lvhalf = - 75.0 ;
   _lsfactor = 5.5 ;
   finf = 1.0 / ( 1.0 + exp ( ( _lv - _lvhalf ) / _lsfactor ) ) ;
   ftau = ( 1.0 / ( exp ( - 14.59 - 0.086 * _lv ) + exp ( - 1.87 + 0.0701 * _lv ) ) ) ;
    return 0; }
 
static void _hoc_setinf(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_setinf(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 setinf ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
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
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  f = f0;
 {
   setinf ( _threadargscomma_ v ) ;
   f = finf ;
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
 _check_setinf(_p, _ppvar, _thread, _nt);
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ih = gh * f * ( v - eih ) ;
   }
 _current += ih;

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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
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
 {   integrate(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(f) - _p;  _dlist1[0] = &(Df) - _p;
   _t_finf = makevector(401*sizeof(double));
   _t_ftau = makevector(401*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/taha/git/jpb-1071/example/Mechanisms/GPi_HCN.mod";
static const char* nmodl_file_text = 
  "TITLE Potassium Ih channel for GPi neuron model\n"
  "\n"
  "COMMENT\n"
  "\n"
  " HCN1 and HCN2 channels are present in the EP (see Fig4 of Chan2004)\n"
  " HCN2 immunoreactivity appeared to be a bit stronger, but this may be\n"
  " due to the methods used in the study.  Ih modeled from an LGN relay neuron:\n"
  "    -Huguenard & McCormick 1992 (recordings at 35.5degC -- disregard Q10)\n"
  "    -McCormick & Pape 1990 (recordings at 36.1degC -- disregard Q10)\n"
  " Code amended from Gillies2006\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  "    (mV) = (millivolt)\n"
  "    (mA) = (milliamp)\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX GPi_Ih\n"
  "    NONSPECIFIC_CURRENT ih\n"
  "    RANGE gh\n"
  "    GLOBAL eih\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    v (mV)\n"
  "    dt (ms)\n"
  "    gh    = 0.001 (mho/cm2)\n"
  "    eih   = -56.11 (mV)\n"
  "    celsius\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    f\n"
  "}\n"
  "\n"
  "ASSIGNED { \n"
  "    ih (mA/cm2)\n"
  "    finf  \n"
  "    ftau (ms)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE integrate METHOD cnexp\n"
  "    ih = gh*f*(v-eih)\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  "    setinf(v)\n"
  "    f = finf\n"
  "}\n"
  "\n"
  "DERIVATIVE integrate {\n"
  "    setinf(v)\n"
  "    f' = (finf - f)/ftau\n"
  "}\n"
  "\n"
  "PROCEDURE setinf(v) {\n"
  "    \n"
  "    LOCAL vhalf, sfactor\n"
  "    TABLE finf, ftau DEPEND celsius FROM -100 TO 100 WITH 400\n"
  "\n"
  "    vhalf = -75     : membrane potential when Gh is half-activated (mV)\n"
  "    sfactor = 5.5   : slope factor which determines the steepness of the fitted curve\n"
  "    finf = 1.0/(1+exp((v-vhalf)/sfactor))\n"
  "    ftau = (1.0/(exp(-14.59-0.086*v)+exp(-1.87+0.0701*v)))\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
