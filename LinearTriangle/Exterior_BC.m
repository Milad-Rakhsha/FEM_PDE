function [ r ] = Exterior_BC(eta)
%This function evaluates the integrand of the mass matrix integral at
% a given point eta and zeta. This is going to be used in the guass
% quadrate method later.
N_=N();
r1=N_{1}(eta,0);
r2=N_{2}(eta,0);
r=[r1;
   r2];
end

