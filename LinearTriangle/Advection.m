function [ C ] = Advection(eta,zeta,invJ,V)
%This function evaluates the integrand of the mass matrix integral at
% a given point eta and zeta. This is going to be used in the guass
% quadrate method later.
Nx=N_x();
N_=N();
N1=N_{1}(eta,zeta);
N2=N_{2}(eta,zeta);
N3=N_{3}(eta,zeta);
V1=V(1,:)';
V2=V(2,:)';
V3=V(3,:)';

N1_eta_zeta=[
    Nx{1}(eta,zeta);
    Nx{2}(eta,zeta)];

N2_eta_zeta=[
    Nx{3}(eta,zeta);
    Nx{4}(eta,zeta)];

N3_eta_zeta=[
    Nx{5}(eta,zeta);
    Nx{6}(eta,zeta)];

c11= N1*(invJ*N1_eta_zeta)'*V1;
c12= N1*(invJ*N2_eta_zeta)'*V2;
c13= N1*(invJ*N3_eta_zeta)'*V3;
c21= N2*(invJ*N1_eta_zeta)'*V1;
c22= N2*(invJ*N2_eta_zeta)'*V2;
c23= N2*(invJ*N3_eta_zeta)'*V3;
c31= N3*(invJ*N1_eta_zeta)'*V1;
c32= N3*(invJ*N2_eta_zeta)'*V2;
c33= N3*(invJ*N3_eta_zeta)'*V3;

C= [c11, c12, c13;
    c21, c22, c23;
    c31, c32, c33];

end