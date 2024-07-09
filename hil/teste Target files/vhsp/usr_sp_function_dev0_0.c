// generated using template: cop_main.template---------------------------------------------
/******************************************************************************************
**
**  Module Name: cop_main.c
**  NOTE: Automatically generated file. DO NOT MODIFY!
**  Description:
**            Main file
**
******************************************************************************************/
// generated using template: arm/custom_include.template-----------------------------------


#ifdef __cplusplus
#include <limits>

extern "C" {
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <complex.h>

// x86 libraries:
#include "../include/sp_functions_dev0.h"


#ifdef __cplusplus
}
#endif



// ----------------------------------------------------------------------------------------                // generated using template:generic_macros.template-----------------------------------------
/*********************** Macros (Inline Functions) Definitions ***************************/

// ----------------------------------------------------------------------------------------

#ifndef MAX
#define MAX(value, limit) (((value) > (limit)) ? (value) : (limit))
#endif
#ifndef MIN
#define MIN(value, limit) (((value) < (limit)) ? (value) : (limit))
#endif

// generated using template: VirtualHIL/custom_defines.template----------------------------

typedef unsigned char X_UnInt8;
typedef char X_Int8;
typedef signed short X_Int16;
typedef unsigned short X_UnInt16;
typedef int X_Int32;
typedef unsigned int X_UnInt32;
typedef unsigned int uint;
typedef double real;

// ----------------------------------------------------------------------------------------
// generated using template: custom_consts.template----------------------------------------

// arithmetic constants
#define C_SQRT_2                    1.4142135623730950488016887242097f
#define C_SQRT_3                    1.7320508075688772935274463415059f
#define C_PI                        3.1415926535897932384626433832795f
#define C_E                         2.7182818284590452353602874713527f
#define C_2PI                       6.283185307179586476925286766559f

//@cmp.def.start
//component defines






float  _single_phase_meter1_meter__Ts;
unsigned int   _single_phase_meter1_meter__flag_init;
float  _single_phase_meter1_meter__Ts2 , _single_phase_meter1_meter__ws , _single_phase_meter1_meter__ws_n1 , _single_phase_meter1_meter__Va , _single_phase_meter1_meter__Va_n1 , _single_phase_meter1_meter__Va_n2 , _single_phase_meter1_meter__Vb , _single_phase_meter1_meter__Vb_n1 , _single_phase_meter1_meter__Vb_n2;
float  _single_phase_meter1_meter__win , _single_phase_meter1_meter__win_n1 , _single_phase_meter1_meter__what , _single_phase_meter1_meter__what2 , _single_phase_meter1_meter__norm;
float  _single_phase_meter1_meter__GAMA , _single_phase_meter1_meter__KE , _single_phase_meter1_meter__Vgrid_n1 , _single_phase_meter1_meter__Vgrid_n2 , _single_phase_meter1_meter__k1 , _single_phase_meter1_meter__k2 , _single_phase_meter1_meter__k3 , _single_phase_meter1_meter__k4 , _single_phase_meter1_meter__k5;
float  _single_phase_meter1_meter__wg , _single_phase_meter1_meter__TWO_PI , _single_phase_meter1_meter__ct1 , _single_phase_meter1_meter__ct2;
float  _single_phase_meter1_meter__MAX_FREQ , _single_phase_meter1_meter__MIN_FREQ;
float  _single_phase_meter1_meter__Tsig;
unsigned int   _single_phase_meter1_meter__max_k;
unsigned int   _single_phase_meter1_meter__kMIN;
unsigned int   _single_phase_meter1_meter__kLIM;
float  _single_phase_meter1_meter__temp;
unsigned int   _single_phase_meter1_meter__cnt_cycles;
float  _single_phase_meter1_meter__nc;
float  _single_phase_meter1_meter__vk_0 , _single_phase_meter1_meter__vk_1;
float  _single_phase_meter1_meter__vt , _single_phase_meter1_meter__t0v , _single_phase_meter1_meter__tfv;
float  _single_phase_meter1_meter__sumV2;
float  _single_phase_meter1_meter__v_rms;
unsigned int   _single_phase_meter1_meter__kv;
float  _single_phase_meter1_meter__Va_k1;
float  _single_phase_meter1_meter__sumS , _single_phase_meter1_meter__Pw , _single_phase_meter1_meter__Pact;
float  _single_phase_meter1_meter__ik_0 , _single_phase_meter1_meter__ik_1;
float  _single_phase_meter1_meter__sumI2;
float  _single_phase_meter1_meter__i_rms;
float  _single_phase_meter1_meter__Qint , _single_phase_meter1_meter__Qt;












//@cmp.def.end


//-----------------------------------------------------------------------------------------
// generated using template: common_variables.template-------------------------------------
// true global variables



// const variables

//@cmp.var.start
// variables
double _single_phase_meter1_grid_frequency__out = 60.0;
double _single_phase_meter1__ig_ia2__out;
double _single_phase_meter1__vg_va2__out;
X_UnInt32 _single_phase_meter1_ncycles__out = 1;
double _single_phase_meter1_meter__Igrid;
double _single_phase_meter1_meter__Vgrid;
double _single_phase_meter1_meter__f_grid;
X_UnInt32 _single_phase_meter1_meter__n_cycles;

double _single_phase_meter1_meter__Irms;
double _single_phase_meter1_meter__Po;
double _single_phase_meter1_meter__Qo;
double _single_phase_meter1_meter__So;
double _single_phase_meter1_meter__Valpha;
double _single_phase_meter1_meter__Vbeta;
double _single_phase_meter1_meter__Vrms;
double _single_phase_meter1_meter__f;
double _single_phase_meter1_meter__pf;
double _single_phase_meter1_out_bus__out[11];
//@cmp.var.end

//@cmp.svar.start
// state variables

//@cmp.svar.end

//
// Tunable parameters
//
static struct Tunable_params {
} __attribute__((__packed__)) tunable_params;

void *tunable_params_dev0_cpu0_ptr = &tunable_params;

// Dll function pointers
#if defined(_WIN64)
#else
// Define handles for loading dlls
#endif








// generated using template: virtual_hil/custom_functions.template---------------------------------
void ReInit_user_sp_cpu0_dev0() {
#if DEBUG_MODE
    printf("\n\rReInitTimer");
#endif
    //@cmp.init.block.start
    HIL_OutAO(0x4002, 0.0f);
    HIL_OutAO(0x400a, 0.0f);
    {
        _single_phase_meter1_meter__TWO_PI = 2.0 * M_PI ;
        _single_phase_meter1_meter__KE = 1.414 ;
        _single_phase_meter1_meter__GAMA = 50.0 ;
        _single_phase_meter1_meter__Ts = 1e-05 ;
        _single_phase_meter1_meter__Ts2 = _single_phase_meter1_meter__Ts * _single_phase_meter1_meter__Ts ;
        _single_phase_meter1_meter__ws = 0.0 ;
        _single_phase_meter1_meter__ws_n1 = 0.0 ;
        _single_phase_meter1_meter__Va = 0.0 ;
        _single_phase_meter1_meter__Va_n1 = 0.0 ;
        _single_phase_meter1_meter__Va_n2 = 0.0 ;
        _single_phase_meter1_meter__Vb = 0.0 ;
        _single_phase_meter1_meter__Vb_n1 = 0.0 ;
        _single_phase_meter1_meter__Vb_n2 = 0.0 ;
        _single_phase_meter1_meter__win = 0.0 ;
        _single_phase_meter1_meter__win_n1 = 0.0 ;
        _single_phase_meter1_meter__Vgrid_n1 = 0.0 ;
        _single_phase_meter1_meter__Vgrid_n2 = 0.0 ;
        _single_phase_meter1_meter__ct1 = 0.5 * _single_phase_meter1_meter__Ts ;
        _single_phase_meter1_meter__ct2 = 2.0 * _single_phase_meter1_meter__Ts * _single_phase_meter1_meter__KE ;
        _single_phase_meter1_meter__MAX_FREQ = 90.0 * _single_phase_meter1_meter__TWO_PI ;
        _single_phase_meter1_meter__MIN_FREQ = 20.0 * _single_phase_meter1_meter__TWO_PI ;
        _single_phase_meter1_meter__flag_init = 15u ;
        if ( _single_phase_meter1_meter__Ts < 25e-6 ) _single_phase_meter1_meter__kMIN = 50 ;
        else if ( _single_phase_meter1_meter__Ts < 60e-6 ) _single_phase_meter1_meter__kMIN = 35 ;
        else if ( _single_phase_meter1_meter__Ts < 110e-6 ) _single_phase_meter1_meter__kMIN = 20 ;
        else _single_phase_meter1_meter__kMIN = 10 ;
        _single_phase_meter1_meter__cnt_cycles = 0 ;
        _single_phase_meter1_meter__max_k = 1200 ;
        _single_phase_meter1_meter__kLIM = _single_phase_meter1_meter__max_k - 10 ;
        _single_phase_meter1_meter__kv = 0 ;
        _single_phase_meter1_meter__vk_1 = 1.0 ;
        _single_phase_meter1_meter__vk_0 = _single_phase_meter1_meter__vk_1 ;
        _single_phase_meter1_meter__t0v = _single_phase_meter1_meter__Ts / 2.0 ;
        _single_phase_meter1_meter__sumV2 = 0.0 ;
        _single_phase_meter1_meter__v_rms = 0 ;
        _single_phase_meter1_meter__Va_k1 = 1 ;
        _single_phase_meter1_meter__ik_1 = 0.002 ;
        _single_phase_meter1_meter__ik_0 = _single_phase_meter1_meter__ik_1 ;
        _single_phase_meter1_meter__sumI2 = 0.0 ;
        _single_phase_meter1_meter__i_rms = 0 ;
        _single_phase_meter1_meter__sumS = 0.0 ;
        _single_phase_meter1_meter__Pw = 0.0 ;
        _single_phase_meter1_meter__Pact = 0.0 ;
        _single_phase_meter1_meter__Qint = 0.0 ;
        _single_phase_meter1_meter__Qt = 0.0 ;
    }
    HIL_OutAO(0x4000, 0.0f);
    HIL_OutAO(0x4001, 0.0f);
    HIL_OutAO(0x4003, 0.0f);
    HIL_OutAO(0x4004, 0.0f);
    HIL_OutAO(0x4005, 0.0f);
    HIL_OutAO(0x4006, 0.0f);
    HIL_OutAO(0x4007, 0.0f);
    HIL_OutAO(0x4008, 0.0f);
    HIL_OutAO(0x4009, 0.0f);
    //@cmp.init.block.end
}


// Dll function pointers and dll reload function
#if defined(_WIN64)
// Define method for reloading dll functions
void ReloadDllFunctions_user_sp_cpu0_dev0(void) {
    // Load each library and setup function pointers
}

void FreeDllFunctions_user_sp_cpu0_dev0(void) {
}

#else
// Define method for reloading dll functions
void ReloadDllFunctions_user_sp_cpu0_dev0(void) {
    // Load each library and setup function pointers
}

void FreeDllFunctions_user_sp_cpu0_dev0(void) {
}
#endif

void load_fmi_libraries_user_sp_cpu0_dev0(void) {
#if defined(_WIN64)
#else
#endif
}


void ReInit_sp_scope_user_sp_cpu0_dev0() {
    // initialise SP Scope buffer pointer
}
// generated using template: common_timer_counter_handler.template-------------------------

/*****************************************************************************************/
/**
* This function is the handler which performs processing for the timer counter.
* It is called from an interrupt context such that the amount of processing
* performed should be minimized.  It is called when the timer counter expires
* if interrupts are enabled.
*
*
* @param    None
*
* @return   None
*
* @note     None
*
*****************************************************************************************/

void TimerCounterHandler_0_user_sp_cpu0_dev0() {
#if DEBUG_MODE
    printf("\n\rTimerCounterHandler_0");
#endif
    //////////////////////////////////////////////////////////////////////////
    // Set tunable parameters
    //////////////////////////////////////////////////////////////////////////
    // Generated from the component: Single-phase Meter1.Grid frequency
    // Generated from the component: Single-phase Meter1.ncycles
//////////////////////////////////////////////////////////////////////////
    // Output block
    //////////////////////////////////////////////////////////////////////////
    //@cmp.out.block.start
    // Generated from the component: Single-phase Meter1._Ig.Ia2
    _single_phase_meter1__ig_ia2__out = (HIL_InFloat(0xc80000 + 0x6));
    // Generated from the component: Single-phase Meter1._Vg.Va2
    _single_phase_meter1__vg_va2__out = (HIL_InFloat(0xc80000 + 0x4));
    // Generated from the component: Single-phase Meter1.I_t
    HIL_OutAO(0x4002, (float)_single_phase_meter1__ig_ia2__out);
    // Generated from the component: Single-phase Meter1.V_t
    HIL_OutAO(0x400a, (float)_single_phase_meter1__vg_va2__out);
    // Generated from the component: Single-phase Meter1.Meter
    _single_phase_meter1_meter__Igrid = _single_phase_meter1__ig_ia2__out;
    _single_phase_meter1_meter__Vgrid = _single_phase_meter1__vg_va2__out;
    _single_phase_meter1_meter__f_grid = _single_phase_meter1_grid_frequency__out;
    _single_phase_meter1_meter__n_cycles = _single_phase_meter1_ncycles__out;
    {
        if ( _single_phase_meter1_meter__flag_init )     {
            _single_phase_meter1_meter__flag_init = 0u ;
            if ( _single_phase_meter1_meter__f_grid > 80.0 )         {
                _single_phase_meter1_meter__wg = 80.0 * _single_phase_meter1_meter__TWO_PI ;
            }
            else if ( _single_phase_meter1_meter__f_grid < 30.0 )         {
                _single_phase_meter1_meter__wg = 30.0 * _single_phase_meter1_meter__TWO_PI ;
            }
            else         {
                _single_phase_meter1_meter__wg = _single_phase_meter1_meter__TWO_PI * _single_phase_meter1_meter__f_grid ;
            }
            _single_phase_meter1_meter__what = _single_phase_meter1_meter__wg ;
            _single_phase_meter1_meter__what2 = _single_phase_meter1_meter__what * _single_phase_meter1_meter__what ;
            _single_phase_meter1_meter__Vrms = 0.0 ;
            _single_phase_meter1_meter__Irms = 0.0 ;
            _single_phase_meter1_meter__Po = 0.0 ;
            _single_phase_meter1_meter__Qo = 0.0 ;
            _single_phase_meter1_meter__So = 0.0 ;
            _single_phase_meter1_meter__pf = 0.0 ;
            _single_phase_meter1_meter__Valpha = 0.0 ;
            _single_phase_meter1_meter__Vbeta = 0.0 ;
        }
        else     {
            _single_phase_meter1_meter__k1 = _single_phase_meter1_meter__ct2 * _single_phase_meter1_meter__what ;
            _single_phase_meter1_meter__k2 = _single_phase_meter1_meter__Ts2 * _single_phase_meter1_meter__what2 + _single_phase_meter1_meter__ct2 * _single_phase_meter1_meter__what + 4 ;
            _single_phase_meter1_meter__k3 = 2 * _single_phase_meter1_meter__Ts2 * _single_phase_meter1_meter__what2 - 8 ;
            _single_phase_meter1_meter__k4 = _single_phase_meter1_meter__Ts2 * _single_phase_meter1_meter__what2 - _single_phase_meter1_meter__ct2 * _single_phase_meter1_meter__what + 4 ;
            _single_phase_meter1_meter__k5 = _single_phase_meter1_meter__Ts2 * _single_phase_meter1_meter__KE * _single_phase_meter1_meter__what2 ;
            _single_phase_meter1_meter__Va = ( _single_phase_meter1_meter__k1 * ( _single_phase_meter1_meter__Vgrid - _single_phase_meter1_meter__Vgrid_n2 ) - _single_phase_meter1_meter__k3 * _single_phase_meter1_meter__Va_n1 - _single_phase_meter1_meter__k4 * _single_phase_meter1_meter__Va_n2 ) / _single_phase_meter1_meter__k2 ;
            _single_phase_meter1_meter__Vb = ( _single_phase_meter1_meter__k5 * ( _single_phase_meter1_meter__Vgrid + 2 * _single_phase_meter1_meter__Vgrid_n1 + _single_phase_meter1_meter__Vgrid_n2 ) - _single_phase_meter1_meter__k3 * _single_phase_meter1_meter__Vb_n1 - _single_phase_meter1_meter__k4 * _single_phase_meter1_meter__Vb_n2 ) / _single_phase_meter1_meter__k2 ;
            _single_phase_meter1_meter__norm = _single_phase_meter1_meter__Va * _single_phase_meter1_meter__Va + _single_phase_meter1_meter__Vb * _single_phase_meter1_meter__Vb ;
            if ( _single_phase_meter1_meter__norm < 0.01 ) _single_phase_meter1_meter__norm = 0.01 ;
            _single_phase_meter1_meter__win = ( - _single_phase_meter1_meter__GAMA * _single_phase_meter1_meter__KE * _single_phase_meter1_meter__what * _single_phase_meter1_meter__Vb ) / _single_phase_meter1_meter__norm * ( _single_phase_meter1_meter__Vgrid - _single_phase_meter1_meter__Va ) ;
            _single_phase_meter1_meter__ws = _single_phase_meter1_meter__ws_n1 + _single_phase_meter1_meter__ct1 * ( _single_phase_meter1_meter__win + _single_phase_meter1_meter__win_n1 ) ;
            _single_phase_meter1_meter__what = _single_phase_meter1_meter__ws + _single_phase_meter1_meter__wg ;
            if ( _single_phase_meter1_meter__what > _single_phase_meter1_meter__MAX_FREQ )         {
                _single_phase_meter1_meter__what = _single_phase_meter1_meter__MAX_FREQ ;
                _single_phase_meter1_meter__ws = _single_phase_meter1_meter__MAX_FREQ - _single_phase_meter1_meter__wg ;
            }
            else if ( _single_phase_meter1_meter__what < _single_phase_meter1_meter__MIN_FREQ )         {
                _single_phase_meter1_meter__what = _single_phase_meter1_meter__MIN_FREQ ;
                _single_phase_meter1_meter__ws = _single_phase_meter1_meter__MIN_FREQ - _single_phase_meter1_meter__wg ;
            }
            _single_phase_meter1_meter__what2 = _single_phase_meter1_meter__what * _single_phase_meter1_meter__what ;
            _single_phase_meter1_meter__ws_n1 = _single_phase_meter1_meter__ws ;
            _single_phase_meter1_meter__Va_n2 = _single_phase_meter1_meter__Va_n1 ;
            _single_phase_meter1_meter__Va_n1 = _single_phase_meter1_meter__Va ;
            _single_phase_meter1_meter__Vb_n2 = _single_phase_meter1_meter__Vb_n1 ;
            _single_phase_meter1_meter__Vb_n1 = _single_phase_meter1_meter__Vb ;
            _single_phase_meter1_meter__win_n1 = _single_phase_meter1_meter__win ;
            _single_phase_meter1_meter__Vgrid_n2 = _single_phase_meter1_meter__Vgrid_n1 ;
            _single_phase_meter1_meter__Vgrid_n1 = _single_phase_meter1_meter__Vgrid ;
            if ( _single_phase_meter1_meter__what >= 30.0 * _single_phase_meter1_meter__TWO_PI && _single_phase_meter1_meter__what <= 80.0 * _single_phase_meter1_meter__TWO_PI )         {
                _single_phase_meter1_meter__vk_0 = _single_phase_meter1_meter__Vgrid ;
                _single_phase_meter1_meter__ik_0 = _single_phase_meter1_meter__Igrid ;
                if ( ( ( ( _single_phase_meter1_meter__Va >= 0.0 && _single_phase_meter1_meter__Va_k1 < 0.0 ) || ( _single_phase_meter1_meter__Va <= 0.0 && _single_phase_meter1_meter__Va_k1 > 0.0 ) ) && _single_phase_meter1_meter__kv > _single_phase_meter1_meter__kMIN ) || _single_phase_meter1_meter__kv > _single_phase_meter1_meter__kLIM )             {
                    _single_phase_meter1_meter__vt = fabs ( _single_phase_meter1_meter__Va ) + fabs ( _single_phase_meter1_meter__Va_k1 ) ;
                    _single_phase_meter1_meter__Tsig = _single_phase_meter1_meter__kv ;
                    if ( _single_phase_meter1_meter__vt > 1e-9 )                 {
                        _single_phase_meter1_meter__tfv = fabs ( _single_phase_meter1_meter__Va_k1 ) / _single_phase_meter1_meter__vt ;
                    }
                    else                 {
                        _single_phase_meter1_meter__tfv = 0 ;
                    }
                    _single_phase_meter1_meter__sumV2 += _single_phase_meter1_meter__vk_1 * _single_phase_meter1_meter__vk_1 * _single_phase_meter1_meter__tfv ;
                    _single_phase_meter1_meter__sumV2 /= _single_phase_meter1_meter__Tsig ;
                    _single_phase_meter1_meter__v_rms += sqrt ( _single_phase_meter1_meter__sumV2 ) ;
                    _single_phase_meter1_meter__sumI2 += _single_phase_meter1_meter__ik_1 * _single_phase_meter1_meter__ik_1 * _single_phase_meter1_meter__tfv ;
                    _single_phase_meter1_meter__sumI2 /= _single_phase_meter1_meter__Tsig ;
                    _single_phase_meter1_meter__i_rms += sqrt ( _single_phase_meter1_meter__sumI2 ) ;
                    _single_phase_meter1_meter__Pw += _single_phase_meter1_meter__vk_1 * _single_phase_meter1_meter__ik_1 * _single_phase_meter1_meter__tfv ;
                    _single_phase_meter1_meter__Pact += _single_phase_meter1_meter__Pw / _single_phase_meter1_meter__Tsig ;
                    _single_phase_meter1_meter__Qint += _single_phase_meter1_meter__Vb * _single_phase_meter1_meter__ik_1 * _single_phase_meter1_meter__tfv ;
                    _single_phase_meter1_meter__Qt += _single_phase_meter1_meter__Qint / _single_phase_meter1_meter__Tsig ;
                    if ( _single_phase_meter1_meter__vt > 1e-9 )                 {
                        _single_phase_meter1_meter__t0v = fabs ( _single_phase_meter1_meter__Va ) / _single_phase_meter1_meter__vt ;
                    }
                    else                 {
                        _single_phase_meter1_meter__t0v = 0 ;
                    }
                    _single_phase_meter1_meter__kv = 0 ;
                    _single_phase_meter1_meter__sumV2 = _single_phase_meter1_meter__vk_0 * _single_phase_meter1_meter__vk_0 * _single_phase_meter1_meter__t0v ;
                    _single_phase_meter1_meter__Pw = _single_phase_meter1_meter__vk_0 * _single_phase_meter1_meter__ik_0 * _single_phase_meter1_meter__t0v ;
                    _single_phase_meter1_meter__Qint = _single_phase_meter1_meter__Vb * _single_phase_meter1_meter__ik_0 * _single_phase_meter1_meter__t0v ;
                    _single_phase_meter1_meter__sumI2 = _single_phase_meter1_meter__ik_0 * _single_phase_meter1_meter__ik_0 * _single_phase_meter1_meter__t0v ;
                    _single_phase_meter1_meter__cnt_cycles ++ ;
                }
                else             {
                    _single_phase_meter1_meter__sumV2 += _single_phase_meter1_meter__vk_0 * _single_phase_meter1_meter__vk_0 ;
                    _single_phase_meter1_meter__sumI2 += _single_phase_meter1_meter__ik_0 * _single_phase_meter1_meter__ik_0 ;
                    _single_phase_meter1_meter__Pw += _single_phase_meter1_meter__vk_0 * _single_phase_meter1_meter__ik_0 ;
                    _single_phase_meter1_meter__Qint += _single_phase_meter1_meter__Vb * _single_phase_meter1_meter__ik_0 ;
                }
                if ( _single_phase_meter1_meter__kv < _single_phase_meter1_meter__max_k ) _single_phase_meter1_meter__kv ++ ;
                _single_phase_meter1_meter__vk_1 = _single_phase_meter1_meter__vk_0 ;
                _single_phase_meter1_meter__Va_k1 = _single_phase_meter1_meter__Va ;
                _single_phase_meter1_meter__ik_1 = _single_phase_meter1_meter__ik_0 ;
            }
            else         {
                _single_phase_meter1_meter__Vrms = 0.0 ;
                _single_phase_meter1_meter__Irms = 0.0 ;
                _single_phase_meter1_meter__Po = 0.0 ;
                _single_phase_meter1_meter__Qo = 0.0 ;
                _single_phase_meter1_meter__So = 0.0 ;
                _single_phase_meter1_meter__pf = 0.0 ;
                _single_phase_meter1_meter__v_rms = 0.0 ;
                _single_phase_meter1_meter__i_rms = 0.0 ;
                _single_phase_meter1_meter__Pact = 0.0 ;
                _single_phase_meter1_meter__Qt = 0.0 ;
                _single_phase_meter1_meter__kv = 0u ;
                _single_phase_meter1_meter__cnt_cycles = 0u ;
                _single_phase_meter1_meter__sumV2 = 0.0 ;
                _single_phase_meter1_meter__sumI2 = 0.0 ;
                _single_phase_meter1_meter__Pw = 0.0 ;
                _single_phase_meter1_meter__Qint = 0.0 ;
            }
        }
        if ( _single_phase_meter1_meter__cnt_cycles >= ( _single_phase_meter1_meter__n_cycles * 2 ) )     {
            _single_phase_meter1_meter__nc = _single_phase_meter1_meter__cnt_cycles ;
            _single_phase_meter1_meter__cnt_cycles = 0 ;
            _single_phase_meter1_meter__Vrms = _single_phase_meter1_meter__v_rms / _single_phase_meter1_meter__nc ;
            _single_phase_meter1_meter__Irms = _single_phase_meter1_meter__i_rms / _single_phase_meter1_meter__nc ;
            _single_phase_meter1_meter__So = _single_phase_meter1_meter__Vrms * _single_phase_meter1_meter__Irms ;
            _single_phase_meter1_meter__Po = _single_phase_meter1_meter__Pact / _single_phase_meter1_meter__nc ;
            _single_phase_meter1_meter__Qo = _single_phase_meter1_meter__Qt / _single_phase_meter1_meter__nc ;
            if ( _single_phase_meter1_meter__So ) _single_phase_meter1_meter__pf = _single_phase_meter1_meter__Po / _single_phase_meter1_meter__So ;
            else _single_phase_meter1_meter__pf = 0 ;
            if ( _single_phase_meter1_meter__Qo < 0 )         {
                _single_phase_meter1_meter__pf = - 1.0 * _single_phase_meter1_meter__pf ;
            }
            _single_phase_meter1_meter__v_rms = 0.0 ;
            _single_phase_meter1_meter__i_rms = 0.0 ;
            _single_phase_meter1_meter__Pact = 0.0 ;
            _single_phase_meter1_meter__Qt = 0.0 ;
        }
        _single_phase_meter1_meter__Valpha = _single_phase_meter1_meter__Va ;
        _single_phase_meter1_meter__Vbeta = - 1.0 * _single_phase_meter1_meter__Vb ;
        _single_phase_meter1_meter__f = _single_phase_meter1_meter__what / _single_phase_meter1_meter__TWO_PI ;
    }
    // Generated from the component: Single-phase Meter1.Freq
    HIL_OutAO(0x4000, (float)_single_phase_meter1_meter__f);
    // Generated from the component: Single-phase Meter1.I_RMS
    HIL_OutAO(0x4001, (float)_single_phase_meter1_meter__Irms);
    // Generated from the component: Single-phase Meter1.POWER_P
    HIL_OutAO(0x4003, (float)_single_phase_meter1_meter__Po);
    // Generated from the component: Single-phase Meter1.POWER_PF
    HIL_OutAO(0x4004, (float)_single_phase_meter1_meter__pf);
    // Generated from the component: Single-phase Meter1.POWER_Q
    HIL_OutAO(0x4005, (float)_single_phase_meter1_meter__Qo);
    // Generated from the component: Single-phase Meter1.POWER_S
    HIL_OutAO(0x4006, (float)_single_phase_meter1_meter__So);
    // Generated from the component: Single-phase Meter1.V_RMS
    HIL_OutAO(0x4007, (float)_single_phase_meter1_meter__Vrms);
    // Generated from the component: Single-phase Meter1.V_alpha
    HIL_OutAO(0x4008, (float)_single_phase_meter1_meter__Valpha);
    // Generated from the component: Single-phase Meter1.V_beta
    HIL_OutAO(0x4009, (float)_single_phase_meter1_meter__Vbeta);
    // Generated from the component: Single-phase Meter1.out_bus
    _single_phase_meter1_out_bus__out[0] = _single_phase_meter1__vg_va2__out;
    _single_phase_meter1_out_bus__out[1] = _single_phase_meter1__ig_ia2__out;
    _single_phase_meter1_out_bus__out[2] = _single_phase_meter1_meter__f;
    _single_phase_meter1_out_bus__out[3] = _single_phase_meter1_meter__Vrms;
    _single_phase_meter1_out_bus__out[4] = _single_phase_meter1_meter__Irms;
    _single_phase_meter1_out_bus__out[5] = _single_phase_meter1_meter__Po;
    _single_phase_meter1_out_bus__out[6] = _single_phase_meter1_meter__Qo;
    _single_phase_meter1_out_bus__out[7] = _single_phase_meter1_meter__So;
    _single_phase_meter1_out_bus__out[8] = _single_phase_meter1_meter__pf;
    _single_phase_meter1_out_bus__out[9] = _single_phase_meter1_meter__Valpha;
    _single_phase_meter1_out_bus__out[10] = _single_phase_meter1_meter__Vbeta;
    // Generated from the component: Single-phase Meter1.meas out
//@cmp.out.block.end
    //////////////////////////////////////////////////////////////////////////
    // Update block
    //////////////////////////////////////////////////////////////////////////
    //@cmp.update.block.start
    // Generated from the component: Single-phase Meter1.Meter
    //@cmp.update.block.end
}
// ----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------