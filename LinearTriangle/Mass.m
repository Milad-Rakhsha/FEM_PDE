function [ M ] = Mass(eta,zeta)
%This function evaluates the integrand of the mass matrix integral at
% a given point eta and zeta. This is going to be used in the guass
% quadrate method later.
N_=N();
m11=N_{1}(eta,zeta)*N_{1}(eta,zeta);
m12=N_{1}(eta,zeta)*N_{2}(eta,zeta);
m13=N_{1}(eta,zeta)*N_{3}(eta,zeta);
m22=N_{2}(eta,zeta)*N_{2}(eta,zeta);
m23=N_{2}(eta,zeta)*N_{3}(eta,zeta);
m33=N_{3}(eta,zeta)*N_{3}(eta,zeta);

M=[m11, m12, m13;
   m12, m22, m23;
   m13, m23, m33];


end

