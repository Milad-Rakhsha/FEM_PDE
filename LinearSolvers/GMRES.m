function [ Sol ] = GMRES(A,b,M,x0,maxIter,residual,restart)

m=restart;
V=zeros(length(x0),m);
H=zeros(m+1,m);

Rtn=zeros(m+1,m);
for iter=1:maxIter
    r0=M\(b-A*x0);
    beta=norm(r0,2);
    V(:,1)=r0/beta;
    for j=1:m
        wj=M\(A*V(:,j));
       for i=1:j
           H(i,j)=wj'*V(:,i);
           wj=wj-H(i,j)*V(:,i);
       end
       H(j+1,j)=norm(wj,2);
       V(:,j+1)=wj/H(j+1,j); 
       Hn=H(1:j+1,1:j);
       Rn=Rtn(1:j,1:j);
    end
    e1=[1.0;zeros(m,1)];
    y = lsqlin(H,beta*e1);
%     length(y)
%     size(V)
    xm=x0+V(:,1:m)*y;
    
    res=norm(M\(b-A*xm),2);
    if (res<residual)
        break;
    end
    x0=xm;
    disp(['norm=',num2str(res)]);
end
Sol=xm;
end

