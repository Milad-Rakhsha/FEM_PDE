clear
clc
addpath('./distmesh');
N_GQ=4;

Xh=1;
fd=@(p) ddiff(drectangle(p,-1,Xh,-1,1),dcircle(p,0,0,0.5));
fh=@(p) 0.07+0.02*dcircle(p,0,0,0.5);
[p,t]=distmesh2d(fd,fh,0.1,[-1,-1;Xh,1.2],[-1,-1;-1,Xh;1,-1;1,Xh]);
rs=sqrt(p(:,1).^2+p(:,2).^2);
Dirichlet_Nodes=find(abs(rs-0.5)<0.001);
Neumann_Nodes=find((abs(p(:,1))>0.99) | (abs(p(:,2))>0.99));
% Neumann_Nodes=find(abs(rs-0.5)<0.001);
% Dirichlet_Nodes=find((abs(p(:,1))>0.99) | (abs(p(:,2))>0.99));
mu=1;
source =1;
v_i=[2,2];

% fd=@(p) sqrt(sum(p.^2,2))-1;
% fh=@(p) 0.04-0.00*dcircle(p,0,0,0.5);
% [p,t]=distmesh2d(fd,fh,0.15,[-1,-1;1,1],[-1,0]);
% Dirichlet_Nodes = unique(boundedges(p,t));
% mu=1;
% source =1;
% v_i=[0,0];


% Dirichlet_Nodes
Ne=size(t,1)
Np=size(p,1)
V=repmat(v_i,Np,1);

addpath('./LinearTriangle');

N=N();
Nx=N_x();
K=zeros(Np,Np);
M=zeros(Np,Np);
r=zeros(Np,1);
[W,Points]=GQ_Points(N_GQ)
[p1,p2]=meshgrid(Points,Points);
gq_W=W'*W;
for i=1:Ne
    node=t(i,:);
    Num_nodes=length(node);
    Mi=zeros(Num_nodes,Num_nodes);
    Ci=zeros(Num_nodes,Num_nodes);
    Di=zeros(Num_nodes,Num_nodes);
    Ri=zeros(Num_nodes,1);
    for j=1:N_GQ*N_GQ
        eta=p1(j);
        zeta=p2(j);
        weight=gq_W(j);
        J=Jacobian(eta,zeta,p(node,:));
        detJ=det(J);
        invJ=inv(J);
        
        Mi=Mi+weight*Mass(eta,zeta)*detJ;
        Di=Di+weight*Diffusion(eta,zeta,invJ)*detJ;
        Ci=Ci+weight*Advection(eta,zeta,invJ,V(node,:))*detJ;
        Ri=Ri+weight*Source(eta,zeta)*detJ*source;
        
        for k1=1:length(node)
            for k2=1:length(node)
                K(node(k1),node(k2))=K(node(k1),node(k2))+Ci(k1,k2)+mu*Di(k1,k2);
                M(node(k1),node(k2))=M(node(k1),node(k2))+Mi(k1,k2);
            end
            r(node(k1))=r(node(k1))+Ri(k1);
        end
        
    end
end

% Apply Dirichlet BCs:
K(Dirichlet_Nodes,:)=0;
K(:,Dirichlet_Nodes)=0;
K(Dirichlet_Nodes,Dirichlet_Nodes) = eye(length(Dirichlet_Nodes),length(Dirichlet_Nodes));
r(Dirichlet_Nodes)=0;

u=K\r;
figure;
trisurf(t,p(:,1),p(:,2),0*p(:,1),u,'edgecolor','k','facecolor','interp');
view(2); axis([-1 1 -1 1]); axis equal; colorbar;colormap(jet(256))
