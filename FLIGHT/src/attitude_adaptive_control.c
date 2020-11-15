#include "attitude_adaptive_control.h"
#include "config.h"
#include "math.h"
#include "maths.h"
#include "matrix.h"



ACobject AcAttitude;

static inline int16_t ControlOutLimit(float in)
{
	if (in > INT16_MAX)
		return INT16_MAX;
	else if (in < -INT16_MAX)
		return -INT16_MAX;
	else
		return (int16_t)in;
}
void attitudeAdadptiveControl(Axis3f gyro, attitude_t* actualAngle, attitude_t* desiredAngle, control_t* output)
{
    fp_angles_t delta;
    float       temp[3] = {0};
    float       matrix[3][3]         = { 0 }; //当前实际欧拉角对应的旋转矩阵
    float       desired_matrix[3][3] = { 0 }; //预达到目标欧拉角对应的旋转矩阵
	
	gyro.x = gyro.x * DEG2RAD;
    gyro.y = gyro.y * DEG2RAD;
    gyro.z = gyro.z * DEG2RAD; 

    delta.angles.roll  = actualAngle->roll * DEG2RAD;
    delta.angles.pitch = actualAngle->pitch * DEG2RAD;
    delta.angles.yaw   = 0;
    buildRotationMatrix(&delta, matrix);

    delta.angles.roll  = desiredAngle->roll * DEG2RAD;
    delta.angles.pitch = desiredAngle->pitch * DEG2RAD;
    delta.angles.yaw   = 0;
    buildRotationMatrix(&delta, desired_matrix);

    float z_desired[3] = { desired_matrix[0][2], desired_matrix[1][2], desired_matrix[2][2] };
    //    float z[3]         = { matrix[0][2], matrix[1][2], matrix[2][2] };
    float x[3] = { matrix[0][0], matrix[1][0], matrix[2][0] };
    float y[3] = { matrix[0][1], matrix[1][1], matrix[2][1] };

    MulMatrixDD(y, z_desired, 1, 3, 1, &AcAttitude.error[0]);
    MulMatrixDD(x, z_desired, 1, 3, 1, &AcAttitude.error[1]);
    AcAttitude.error[1] = -AcAttitude.error[1];
    // AcAttitude.error[2] = 0;
    MulMatrixDD(&AcAttitude.eGain, AcAttitude.error, 3, 3, 1, temp);
    // e_gain * error
    for (int i = 0; i < 3; i++)
        temp[i] = AcAttitude.error[i] * AcAttitude.eGain;
    // Sa = w + e_gain * error;
   
    MatrixADD(temp, gyro.axis, 3, 1, AcAttitude.Sa);
    // e_deriv
    for (int i = 0; i < 3; i++) {
        AcAttitude.errorderiv[i] = (AcAttitude.error[i] - AcAttitude.prevError[i]) / AcAttitude.dt;
        AcAttitude.prevError[i]  = AcAttitude.error[i];
		AcAttitude.error[i] = 0;
    }
    /*
    // J*w
    float J_w[3];
    MulMatrixDD(&acobject->J_e, &gyro, 3, 3, 1, &J_w);
    // eGain*e
    float eGain_error[3];
    MulMatrixDD(&acobject->eGain, &acobject->error, 3, 3, 1, &eGain_error);
    // eGain*e_deriv
    float eGain_errorderiv[3];
    MulMatrixDD(&acobject->eGain, &acobject->errorderiv, 3, 3, 1, &eGain_errorderiv);
    // J_eGain_e_deriv
    float J_e_eGain_errorderiv[3];
    MulMatrixDD(&acobject->J_e, &eGain_errorderiv, 3, 3, 1, &J_e_eGain_errorderiv);
    // eGain*e × J*w
    float eGain_error_cross_J_w[3];
    vector3_crossproduct(&eGain_error,&J_w,&eGain_error_cross_J_w);
    // Ka*Sa
    float Ka_Sa[3];
    MulMatrixDD(&acobject->Ka, &acobject->Sa, 3, 3, 1, &Ka_Sa);
    float output_temp[3];
    MatrixADD(&Ka_Sa, &eGain_error_cross_J_w, 3, 1, &output_temp);
    MatrixADD(&output_temp, &J_e_eGain_errorderiv, 3, 1, &output_temp);
    */
    float Y[3][6] = { 0 };
    Y[0][0]       = -AcAttitude.eGain * AcAttitude.errorderiv[0];
    Y[0][1]       = AcAttitude.eGain * AcAttitude.errorderiv[2] * gyro.y;
    Y[0][2]       = -AcAttitude.eGain * AcAttitude.errorderiv[1] * gyro.z;
    Y[1][0]       = -AcAttitude.eGain * AcAttitude.errorderiv[2] * gyro.x;
    Y[1][1]       = -AcAttitude.eGain * AcAttitude.errorderiv[1];
    Y[1][2]       = AcAttitude.eGain * AcAttitude.errorderiv[0] * gyro.z;
    Y[2][0]       = AcAttitude.eGain * AcAttitude.errorderiv[1] * gyro.x;
    Y[2][1]       = -AcAttitude.eGain * AcAttitude.errorderiv[0] * gyro.y;
    Y[2][2]       = -AcAttitude.eGain * AcAttitude.errorderiv[2];
    Y[0][3]       = 1;
    Y[1][4]       = 1;
    Y[2][5]       = 1;
    for (int i = 0; i < 3; i++) {
        AcAttitude.alpha[i]     = AcAttitude.J_e[i];
        AcAttitude.alpha[i + 3] = AcAttitude.T0_e[i];
    }
    float Y_alpha[3] = { 0 };
    MulMatrixDD(&Y[0][0], AcAttitude.alpha, 3, 6, 1, Y_alpha);
    float Ka_Sa[3] = {0};
    MulMatrixDD(&AcAttitude.Ka[0][0], AcAttitude.Sa, 3, 3, 1, Ka_Sa);
    output->roll  = ControlOutLimit(-Ka_Sa[0] + Y_alpha[0]);
    output->pitch = ControlOutLimit(-Ka_Sa[1] + Y_alpha[1]);
//    output->yaw  = -Ka_Sa[2] + Y_alpha[2];

    float YT[6][3] = { 0 };
    TransMatrixD(&Y[0][0], 3, 6, &YT[0][0]);
    float YT_Sa[6] = { 0 };
    MulMatrixDD(&YT[0][0], AcAttitude.Sa, 6, 3, 1, YT_Sa);
    float Aderiv[6] = { 0 };
    MulMatrixDD(&AcAttitude.AGain[0][0], YT_Sa, 6, 6, 1, Aderiv);

    for (int i = 0; i < 6; i++)
        AcAttitude.alpha[i] -= Aderiv[i] *  AcAttitude.dt;
}


void attitudeAdadptiveControlInit(const float dt)
{
    for (int i = 0; i < 3; i++) {
        AcAttitude.error[i]      = 0;
        AcAttitude.errorderiv[i] = 0;
        AcAttitude.prevError[i]  = 0;
        AcAttitude.Sa[i]         = 0;
        AcAttitude.T0_e[i]       = 0;
        AcAttitude.J_e[i]        = 0;
        for (int j = 0; j < 3; j++) {
            AcAttitude.Ka[i][j] = 0;
        }
    }

    AcAttitude.eGain = 8;

    AcAttitude.Ka[0][0] = 4000;
    AcAttitude.Ka[1][1] = 4000;
    AcAttitude.Ka[2][2] = 4000;

    AcAttitude.T0_e[0] = 0;
    AcAttitude.T0_e[1] = 0;
    AcAttitude.T0_e[2] = 0;

    AcAttitude.J_e[0] = 0.0211f;
    AcAttitude.J_e[1] = 0.0211f;
    AcAttitude.J_e[2] = 0.0366F;
    
    for(int i=0; i<6; i++)
    {
        for(int j=0; j<6; j++)
            AcAttitude.AGain[i][j] = 0;
        AcAttitude.AGain[i][i] = 1;
    }
    AcAttitude.dt = dt;
}

void attitudeAdadptiveControlReset(void)
{
    for (int i = 0; i < 3; i++) {
        AcAttitude.error[i]      = 0;
        AcAttitude.errorderiv[i] = 0;
        AcAttitude.prevError[i]  = 0;
        AcAttitude.Sa[i]         = 0;
    }

    AcAttitude.T0_e[0] = 1;
    AcAttitude.T0_e[1] = 1;
    AcAttitude.T0_e[2] = 1;

    AcAttitude.J_e[0] = 1;
    AcAttitude.J_e[1] = 1;
    AcAttitude.J_e[2] = 1;
}
