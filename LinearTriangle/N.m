function out = N()
N1=@(eta,zeta) 1-eta-zeta;
N2=@(eta,zeta) eta;
N3=@(eta,zeta) zeta;
out = {N1, N2, N3};
end

