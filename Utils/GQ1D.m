function [ int_F ] = GQ1D(fun,a,b,Np)
%GQ2D integrates the fun(x,y) over the range of x:[a,b] y:[c,d] 
% using Np quadrature points

if(Np==1)
    GQP=[0];
    coeff=[2];
elseif(Np==2)
    GQP=[-1/sqrt(3)  +1/sqrt(3) ];
    coeff=[1    1];
elseif(Np==3)
    GQP=[-sqrt(3/5)  0 +sqrt(3/5)];
    coeff=[5/9 8/9 5/9];
elseif(Np==4)
    GQP=[-0.861136 -0.339981 +0.339981 +0.861136];
    coeff=[0.347855 0.652145 0.652145 0.347855];
elseif(Np==5)    
    GQP=[-0.90618 -0.538469  0 +0.538469 +0.90618];
    coeff=[0.236927 0.478629 0.568889 0.478629 0.236927];
else
    disp(['please implement N=', num2str(Np),'GQ2D.m'])
end
       
    %calculate the function value at
    F=fun(GQP);
    int_F=(b-a)/2*sum(coeff.*F);
       
end

