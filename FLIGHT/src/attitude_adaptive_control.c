#include "attitude_adaptive_control.h"
#include "maths.h"
#include "matrix.h"

void attitudeAdadptiveControl(ACobject *acobject, attitude_t *actualAngle,attitude_t *desiredAngle,control_t *output)
{
    fp_angles_t delta;
    float matrix[3][3] = {0};
    float desired_matrix[3][3] = {0};

    delta.angles.roll = actualAngle->roll;
    delta.angles.pitch = actualAngle->pitch;
    delta.angles.yaw = 0;    
    buildRotationMatrix(&delta, matrix);

    delta.angles.roll = desiredAngle->roll;
    delta.angles.pitch = desiredAngle->pitch;
    delta.angles.yaw = 0;    
    buildRotationMatrix(&delta, desired_matrix);

    float z_desired[3] = {desired_matrix[1][3],desired_matrix[2][3],desired_matrix[3][3]};
    float z[3] = {matrix[0][2],matrix[1][2],matrix[2][2]};
    float x[3] = {matrix[0][0],matrix[1][0],matrix[2][0]};
    float y[3] = {matrix[0][1],matrix[1][1],matrix[2][1]};    
    
    MulMatrixDD(y, z_desired, 1,3,1,&acobject->error[0]);
    MulMatrixDD(x, z_desired, 1,3,1,&acobject->error[1]);

    


}
