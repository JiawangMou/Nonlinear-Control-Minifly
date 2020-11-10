#ifndef __ATTITUDE_ADAPTIVE_CONTROL_H
#define __ATTITUDE_ADAPTIVE_CONTROL_H
#include "stabilizer_types.h"

typedef struct
{
	float error[3];        //< error
	float prevError[3];    //< previous error
	float errorderiv[3];   //< derivative
	float kp;           //< proportional gain
	float ki;           //< integral gain
	float kd;           //< derivative gain
	float outP;         //< proportional output (debugging)
	float outI;         //< integral output (debugging)
	float outD;         //< derivative output (debugging)
	float iLimit;       //< integral limit
	float outputLimit;  //< total PID output limit, absolute value. '0' means no limit.
	float dt;           //< delta-time dt
	float out;			//< out
} ACobject;

void attitudeAdadptiveControl(attitude_t *actualAngle,attitude_t *desiredAngle,control_t *output);

#endif