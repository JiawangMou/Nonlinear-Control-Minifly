#include "attitude_adaptive_control.h"
#include "config.h"
#include "maths.h"
#include "matrix.h"

#ifdef ADAPTIVE_CONTROL
ACobject AcAttitude;
#endif

void attitudeAdadptiveControl(
    ACobject* acobject, Axis3f gyro, attitude_t* actualAngle, attitude_t* desiredAngle, control_t* output)
{
    fp_angles_t delta;
    float       temp[3];
    float       matrix[3][3]         = { 0 }; //当前实际欧拉角对应的旋转矩阵
    float       desired_matrix[3][3] = { 0 }; //预达到目标欧拉角对应的旋转矩阵

    delta.angles.roll  = actualAngle->roll;
    delta.angles.pitch = actualAngle->pitch;
    delta.angles.yaw   = 0;
    buildRotationMatrix(&delta, matrix);

    delta.angles.roll  = desiredAngle->roll;
    delta.angles.pitch = desiredAngle->pitch;
    delta.angles.yaw   = 0;
    buildRotationMatrix(&delta, desired_matrix);

    float z_desired[3] = { desired_matrix[1][3], desired_matrix[2][3], desired_matrix[3][3] };
    float z[3]         = { matrix[0][2], matrix[1][2], matrix[2][2] };
    float x[3]         = { matrix[0][0], matrix[1][0], matrix[2][0] };
    float y[3]         = { matrix[0][1], matrix[1][1], matrix[2][1] };

    MulMatrixDD(&y, &z_desired, 1, 3, 1, &acobject->error[0]);
    MulMatrixDD(&x, &z_desired, 1, 3, 1, &acobject->error[1]);
    // acobject->error[2] = 0;
    MulMatrixDD(&acobject->eGain, &acobject->error, 3, 3, 1, &temp);
    // Sa = w + e_gain * error;
    MatrixADD(&temp, &gyro, 3, 1, &acobject->Sa);
    // e_deriv
    for (int i = 0; i < 3; i++) {
        acobject->errorderiv[i] = (acobject->error[i] - acobject->prevError[i]) / acobject->dt;
        acobject->prevError[i]  = acobject->error[i];
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
    Y[0][0]       = -acobject->eGain * acobject->errorderiv[0];
    Y[0][1]       =  acobject->eGain * acobject->errorderiv[2] * gyro[1];
    Y[0][2]       = -acobject->eGain * acobject->errorderiv[1] * gyro[2];
    Y[1][0]       = -acobject->eGain * acobject->errorderiv[3] * gyro[0];
    Y[1][1]       = -acobject->eGain * acobject->errorderiv[1];
    Y[1][2]       =  acobject->eGain * acobject->errorderiv[0] * gyro[2];
    Y[2][0]       =  acobject->eGain * acobject->errorderiv[1] * gyro[0];
    Y[2][1]       = -acobject->eGain * acobject->errorderiv[0] * gyro[1];
    Y[2][2]       = -acobject->eGain * acobject->errorderiv[2];

    for (int i = 0; i < 3; i++) {
        acobject->alpha[i]     = acobject->J_e[i];
        acobject->alpha[i + 3] = acobject->T0_e[i];
    }
    float Y_alpha[3] = {0};
    MulMatrixDD(&Y, &alpha, 3, 6, 1, &Y_alpha);
    float Ka_Sa[3];
    MulMatrixDD(&acobject->Ka, &acobject->Sa, 3, 3, 1, &Ka_Sa);
    output->roll  = -Ka_Sa[0] + Y_alpha[0];
    output->pitch = -Ka_Sa[1] + Y_alpha[1];
    output->yaw   = -Ka_Sa[2] + Y_alpha[2];

    float YT[6][3] = {0};
    TransMatrixD(&Y,3,6,&YT);
    float YT_Sa[6] = {0};
    MulMatrixDD(&YT,&acobject->Sa,6,3,1,&YT_Sa);
    float Aderiv[6] = {0};
    MulMatrixDD(&acobject->AGain,&YT_Sa,6,6,1,&Aderiv);

    for(int i=0; i<6; i++)
        acobject->alpha[i] -= Aderiv[i];

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
            AcAttitude.eGain[i][j] = 0;
            AcAttitude.Ka[i][j]    = 0;
        }
    }

    AcAttitude.eGain = 1;


    AcAttitude.Ka[0][0] = 3;
    AcAttitude.Ka[1][1] = 3;
    AcAttitude.Ka[2][2] = 3;

    AcAttitude.T0_e[0] = 1;
    AcAttitude.T0_e[1] = 1;
    AcAttitude.T0_e[2] = 1;

    AcAttitude.J_e[0] = 1;
    AcAttitude.J_e[1] = 1;
    AcAttitude.J_e[2] = 1;

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