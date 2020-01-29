function [My_tot, My_vib, My_rot] = getRoVibTransDipoleMomentY(v1, v2, j1, j2, m1, m2)

j3 = j2;
m3 = m2;
j2 = 1;
m2_p1 = 1;
m2_m1 = -1;

i = sqrt(complex(-1));

j123 =    [   j1     j2  j3 ];
m000 =    [    0      0   0 ];
m123_p1 = [  -m1  m2_p1  m3 ];
m123_m1 = [  -m1  m2_m1  m3 ];

% Wigner3j(j123, m000)*Wigner3j(j123, m123_m1) 
% Wigner3j(j123, m000)*Wigner3j(j123, m123_p1)

if abs(j1 - j3) == 1 && abs(m1-m3) == 1 && abs(m1) <= j1 && abs(m3) <= j3
    My_rot = i*(-1)^m1*sqrt(2*pi/3)*sqrt(((2*j1+1)*(2*j2+1)*(2*j3+1))/(4*pi)) ... 
        * (Wigner3j(j123, m000)*(Wigner3j(j123, m123_m1) + Wigner3j(j123, m123_p1)));
    
    if v2 - v1 == 1
        My_vib = sqrt(v2);
    elseif v2 - v1 == -1
        My_vib = sqrt(v1);
    else
        My_vib = 0;
    end
else
    
    My_rot = 0;
    if v2 - v1 == 1
        My_vib = sqrt(v2);
    elseif v2 - v1 == -1
        My_vib = sqrt(v1);
    else
        My_vib = 0;
    end
    
end

% My_rot = i*(-1)^m1*sqrt(2*pi/3)*sqrt(((2*j1+1)*(2*j2+1)*(2*j3+1))/(4*pi)) ... 
%     * (Wigner3j(j123, m000)*(Wigner3j(j123, m123_m1) + Wigner3j(j123, m123_p1)));


My_tot = My_vib * My_rot;