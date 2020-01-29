function [Mz_tot, Mz_vib, Mz_rot] = getRoVibTransDipoleMomentZ(v1, v2, j1, j2, m1, m2)

j3 = j2;
m3 = m2;
j2 = 1;
m2 = 0;

j123 = [  j1 j2  j3 ];
m123 = [ -m1 m2  m3 ];
m000 = [   0  0   0 ];

% Wig1 = Wigner3j(j123, m000)
% Wig2 = Wigner3j(j123, m123)

% if j1 == j3
%     Mz_rot = 0;
%     if v2 - v1 == 1
%         Mz_vib = sqrt(v2);
%     elseif v2 - v1 == -1
%         Mz_vib = sqrt(v1);
%     else
%         Mz_vib = 0;
%     end
if abs(j3 - j1) == 1 && m1 == m3 && abs(m1) <= j1 && abs(m3) <= j3
    
    Mz_rot = (-1)^m1 * sqrt(4*pi/3)*sqrt(((2*j1+1)*(2*j2+1)*(2*j3+1))/(4*pi))*Wigner3j(j123, m000)*Wigner3j(j123, m123);
    if v2 - v1 == 1
        Mz_vib = sqrt(v2);
    elseif v2 - v1 == -1
        Mz_vib = sqrt(v1);
    else
        Mz_vib = 0;
    end    
    
else
    
    Mz_rot = 0;
    if v2 - v1 == 1
        Mz_vib = sqrt(v2);
    elseif v2 - v1 == -1
        Mz_vib = sqrt(v1);
    else
        Mz_vib = 0;
    end
end

% Mz_rot = (-1)^m1 * sqrt(4*pi/3)*sqrt(((2*j1+1)*(2*j2+1)*(2*j3+1))/(4*pi))*Wigner3j(j123, m000)*Wigner3j(j123, m123);
% 
% if v2 - v1 == 1
%     Mz_vib = sqrt(v2);
% 
% elseif v2 - v1 == -1
%     Mz_vib = sqrt(v1);
%     
% else
%     Mz_vib = 0;
% end

Mz_tot = Mz_vib * Mz_rot;
    