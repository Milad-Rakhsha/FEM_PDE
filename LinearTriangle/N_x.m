function out = N_x()
% Note that this cannot return a vector values for the purpose of 
N11=@(eta,zeta) -1;
N12=@(eta,zeta) -1;

N21=@(eta,zeta) 1;
N22=@(eta,zeta) 0;

N31=@(eta,zeta) 0;
N32=@(eta,zeta) 1;
out = {N11, N12, N21, N22, N31, N32};
end
