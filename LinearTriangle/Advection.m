function [ Cx,Cy ] = Advection(eta,zeta,invJ)
%This function evaluates the integrand of the mass matrix integral at
% a given point eta and zeta. This is going to be used in the guass
% quadrate method later.
Nx=N_x();
N_=N();
N1=N_{1}(eta,zeta);
N2=N_{2}(eta,zeta);
N3=N_{3}(eta,zeta);


N1_eta_zeta=[
    Nx{1}(eta,zeta);
    Nx{2}(eta,zeta)];

N2_eta_zeta=[
    Nx{3}(eta,zeta);
    Nx{4}(eta,zeta)];

N3_eta_zeta=[
    Nx{5}(eta,zeta);
    Nx{6}(eta,zeta)];

c11= N1*(invJ*N1_eta_zeta);
c12= N1*(invJ*N2_eta_zeta);
c13= N1*(invJ*N3_eta_zeta);
c21= N2*(invJ*N1_eta_zeta);
c22= N2*(invJ*N2_eta_zeta);
c23= N2*(invJ*N3_eta_zeta);
c31= N3*(invJ*N1_eta_zeta);
c32= N3*(invJ*N2_eta_zeta);
c33= N3*(invJ*N3_eta_zeta);

Cx=[c11(1), c12(1), c13(1);
    c21(1), c22(1), c23(1);
    c31(1), c32(1), c33(1)];

Cy=[c11(2), c12(2), c13(2);
    c21(2), c22(2), c23(2);
    c31(2), c32(2), c33(2)];


end