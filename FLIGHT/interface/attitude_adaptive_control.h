/*******************************************************************************
 * @Copylift (c) 2020, Jiawang Mou, Inc.
 * @
 * @ pluginTemplate : [description]
 * @
 * @ filename : attitude_adaptive_control.h
 * @ author   : Jiawang Mou(moujiawang@sjtu.edu.cn)
 * @ create   : 2020/11/11 	 10:46:25
 ******************************************************************************/

#ifndef __ATTITUDE_ADAPTIVE_CONTROL_H__
#define __ATTITUDE_ADAPTIVE_CONTROL_H__

////////////////////////////////////////////////////////////////////////////////
// Headers
//
#include "stabilizer_types.h"

////////////////////////////////////////////////////////////////////////////////
// Typedefs & Constants
//
typedef struct {
    float error[3];      //< error
    float prevError[3];  //< previous error
    float errorderiv[3]; //< derivative
    float eGain;         //< error gain in Sa
    float Sa[3];         //< Sa
    float Ka[3][3];      //< positive diagonal gain matrix of Sa
    float T0_e[3];       //< initial torque estimate
    float alpha[6];      //< J_e and T0_e
    float J_e[3];        //< inertia matrix estimate
    float outputLimit;   //< total PID output limit, absolute value. '0' means no limit.
    float AGain[6][6];   //< gains of adaptive member
    float dt;            //< delta-time dt
    float out;           //< out
} ACobject;
////////////////////////////////////////////////////////////////////////////////
// Classes
//

////////////////////////////////////////////////////////////////////////////////
// Functions
//
void attitudeAdadptiveControl(Axis3f gyro, attitude_t* actualAngle, attitude_t* desiredAngle, control_t* output);
void attitudeAdadptiveControlInit(const float dt);
#endif //__ATTITUDE_ADAPTIVE_CONTROL_H__
////////////////////////////////// EOF /////////////////////////////////////////
