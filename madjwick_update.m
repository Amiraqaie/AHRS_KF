function q = madgwick_update(q, gyro, acc, mag, beta, dt)
% q     : quaternion [q0 q1 q2 q3]'  (scalar-first)
% gyro  : rad/s [gx gy gz]'
% acc   : m/s^2 [ax ay az]'
% mag   : magnetometer [mx my mz]'
% beta  : algorithm gain
% dt    : timestep

q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

gx = gyro(1); gy = gyro(2); gz = gyro(3);
ax = acc(1);  ay = acc(2);  az = acc(3);
mx = mag(1);  my = mag(2);  mz = mag(3);

%% Normalize accelerometer
if norm([ax ay az]) == 0
    return
end
acc = [ax ay az] / norm([ax ay az]);
ax = acc(1); ay = acc(2); az = acc(3);

%% Normalize magnetometer
if norm([mx my mz]) == 0
    return
end
mag = [mx my mz] / norm([mx my mz]);
mx = mag(1); my = mag(2); mz = mag(3);

%% Reference direction of Earth's magnetic field
hx = 2*mx*(0.5 - q2^2 - q3^2) + 2*my*(q1*q2 - q0*q3) + 2*mz*(q1*q3 + q0*q2);
hy = 2*mx*(q1*q2 + q0*q3) + 2*my*(0.5 - q1^2 - q3^2) + 2*mz*(q2*q3 - q0*q1);
bz = 2*mx*(q1*q3 - q0*q2) + 2*my*(q2*q3 + q0*q1) + 2*mz*(0.5 - q1^2 - q2^2);
bx = sqrt(hx^2 + hy^2);

%% Gradient descent correction step
F = [
    2*(q1*q3 - q0*q2) - ax
    2*(q0*q1 + q2*q3) - ay
    2*(0.5 - q1^2 - q2^2) - az
    2*bx*(0.5 - q2^2 - q3^2) + 2*bz*(q1*q3 - q0*q2) - mx
    2*bx*(q1*q2 - q0*q3) + 2*bz*(q0*q1 + q2*q3) - my
    2*bx*(q0*q2 + q1*q3) + 2*bz*(0.5 - q1^2 - q2^2) - mz
];

J = [
    -2*q2,        2*q3,       -2*q0,        2*q1
     2*q1,        2*q0,        2*q3,        2*q2
     0,          -4*q1,       -4*q2,        0
    -2*bz*q2,     2*bz*q3,    -4*bx*q2-2*bz*q0, -4*bx*q3+2*bz*q1
    -2*bx*q3+2*bz*q1, 2*bx*q2+2*bz*q0, 2*bx*q1+2*bz*q3, -2*bx*q0+2*bz*q2
     2*bx*q2,     2*bx*q3-4*bz*q1, 2*bx*q0-4*bz*q2, 2*bx*q1
];

step = J' * F;
step = step / norm(step);

%% Quaternion derivative
q_dot = 0.5 * [
    -q1*gx - q2*gy - q3*gz
     q0*gx + q2*gz - q3*gy
     q0*gy - q1*gz + q3*gx
     q0*gz + q1*gy - q2*gx
] - beta * step;

%% Integrate
q = q + q_dot * dt;
q = q / norm(q);
end
