
% roll = 0;
% pitch = 0;
% roll_d = 1;
% pitch_d = 1;
% yawRate = [0;0;0];

S = 'XYZ';
% coder.extrinsic('angle2dcm');
roll
pitch
roll_d
pitch_d
AngRATE
a_e

ModelParam_uavJxx =0.0211;%x轴转动惯量（单位： kg.m^2）
ModelParam_uavJyy =0.0211;%y轴转动惯量（单位： kg.m^2）
ModelParam_uavJzz =0.0366;%z轴转动惯量（单位： kg.m^2）
ModelParam_uavJ_T0 = [ModelParam_uavJxx;ModelParam_uavJyy;ModelParam_uavJzz;0;0;0];
R = zeros(3,3);
R1 = zeros(3,3);
x = zeros(3,1);
y = zeros(3,1);

R = angle2dcm(roll,pitch,0,S) 
R1= angle2dcm(roll_d,pitch_d,0,S)
x = [R(1,1);R(2,1);R(3,1)];
y = [R(1,2);R(2,2);R(3,2)];

z_d = [R1(1,3);R1(2,3);R1(3,3)];

%%m = [1,0,-sin(pitch);0,cos(roll),cos(roll)*sin(pitch);0,-sin(roll),cos(pitch)*cos(roll)];
e = [-dot(y,z_d);dot(x,z_d);0]
edot = diff(e)

Y = zeros(3,6);
T = [1,0,0,0,0,0;
     0,1,0,0,0,0;
     0,0,1,0,0,0;
     0,0,0,1,0,0;
     0,0,0,0,1,0;
     0,0,0,0,0,1];
Sa = T*e + AngRATE;

Y(1,1) = -T*edot(1,1);
Y(1,2) =  T*edot(3,1)*AngRATE(2,1);
Y(1,3) = -T*edot(2,1)*AngRATE(3,1);
Y(2,1) = -T*edot(3,1)*AngRATE(1,1);
Y(2,2) = -T*edot(2,1);
Y(2,3) =  T*edot(1,1)*AngRATE(3,1);
Y(3,1) =  T*edot(2,1)*AngRATE(1,1);
Y(3,2) = -T*edot(1,1)*AngRATE(2,1);
Y(3,3) = -T*edot(3,1);
Y(1,4) = 1;
Y(2,5) = 1;
Y(3,6) = 1;
Y
Ka = 3;
Torque_sum = -Ka*Sa + Y*a_e

adot = -T2*Y'*Sa;