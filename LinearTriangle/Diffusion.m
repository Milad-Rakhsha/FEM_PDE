function [ D ] = Diffusion(eta,zeta,invJ)
%This function evaluates the integrand of the mass matrix integral at
% a given point eta and zeta. This is going to be used in the guass
% quadrate method later.
Nx=N_x();

N1_eta_zeta=[
    Nx{1}(eta,zeta);
    Nx{2}(eta,zeta)];

N2_eta_zeta=[
    Nx{3}(eta,zeta);
    Nx{4}(eta,zeta)];

N3_eta_zeta=[
    Nx{5}(eta,zeta);
    Nx{6}(eta,zeta)];


d11=(invJ*N1_eta_zeta)'*(invJ*N1_eta_zeta);
d12=(invJ*N1_eta_zeta)'*(invJ*N2_eta_zeta);
d13=(invJ*N1_eta_zeta)'*(invJ*N3_eta_zeta);
d22=(invJ*N2_eta_zeta)'*(invJ*N2_eta_zeta);
d23=(invJ*N2_eta_zeta)'*(invJ*N3_eta_zeta);
d33=(invJ*N3_eta_zeta)'*(invJ*N3_eta_zeta);

D= [d11, d12, d13;
    d12, d22, d23;
    d13, d23, d33];

end

