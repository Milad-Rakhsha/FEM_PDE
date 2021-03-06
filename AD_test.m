clear
clc
addpath('./distmesh');
addpath('./LinearSolvers');
disp('Meshing...');

N_GQ=2;

Xh=1;
fd=@(p) ddiff(drectangle(p,-Xh,Xh,-1,1),dcircle(p,0,0,0.5));
fh=@(p) 0.2+0.00*dcircle(p,0,0,0.5);
[p,t]=distmesh2d(fd,fh,0.1,[-1.1*Xh,-1.2;1.1*Xh,1.2],[-Xh,-1;-Xh,1;Xh,-1;Xh,1]);
rs=sqrt(p(:,1).^2+p(:,2).^2);
circle_Nodes=find(abs(rs-0.5)<0.001);
LR_nodes=find((abs(p(:,1))>0.99) | (abs(p(:,2))>0.99));
Neumann_Nodes=circle_Nodes;
Dirichlet_Nodes=LR_nodes;
Dirichlet_Values=0;
Neumann_Value=0.0;
mu=1.0;
source =1;
v_i=[1,0];

% fd=@(p) sqrt(sum(p.^2,2))-1;
% fh=@(p) 0.04-0.00*dcircle(p,0,0,0.5);
% [p,t]=distmesh2d(fd,fh,0.15,[-1,-1;1,1],[-1,0]);
% Dirichlet_Nodes = unique(boundedges(p,t));
% mu=1;
% source =1;
% v_i=[0,0];


% Dirichlet_Nodes
Ne=size(t,1);
Np=size(p,1);
V=repmat(v_i,Np,1);

addpath('./LinearTriangle');

N=N();
Nx=N_x();
K=zeros(Np,Np);
M=zeros(Np,Np);
Cx=zeros(Np,Np);
Cy=zeros(Np,Np);
r=zeros(Np,1);
[W,Points]=GQ_Points(N_GQ);
[p1,p2]=meshgrid(Points,Points);
gq_W=W'*W;
disp('Building System matrices...');
for i=1:Ne
    node=t(i,:);
    Num_nodes=length(node);
    Mi=zeros(Num_nodes,Num_nodes);
    Cix=zeros(Num_nodes,Num_nodes);
    Ciy=zeros(Num_nodes,Num_nodes);
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
        [cx,cy]=Advection(eta,zeta,invJ);
        Cix=Cix+weight*cx*detJ;
        Ciy=Ciy+weight*cy*detJ;
        Ri=Ri+weight*Source(eta,zeta)*detJ*source;
        
        for k1=1:length(node)
            for k2=1:length(node)
                K(node(k1),node(k2))=K(node(k1),node(k2))+...
                    (V(k2,1)*Cix(k1,k2)+V(k2,2)*Ciy(k1,k2))+...
                    mu*Di(k1,k2);
                Cx(node(k1),node(k2))=Cx(node(k1),node(k2))+Cix(k1,k2);
                Cy(node(k1),node(k2))=Cy(node(k1),node(k2))+Ciy(k1,k2);                M(node(k1),node(k2))=M(node(k1),node(k2))+Mi(k1,k2);
            end
            r(node(k1))=r(node(k1))+Ri(k1);
        end
        
    end
end
disp('Imposing Neumann BCs...');
%%% Element exterior
for i=1:Ne
    node=t(i,:);
    Num_nodes=length(node);
    Ri=zeros(2,1);%since this is a linear element 2 nodes per edge
    Num_edge=Num_nodes;
    
    if (size(intersect(node,Neumann_Nodes),1)<2)
        continue;
    end
    
    %%%loop over edges of this element
    for edge=1:Num_edge
        node1=edge;
        node2=edge+1;
        if (edge==Num_edge); node2=1; end
        %%% check to see if this edge is on the boundary
        if(any(Neumann_Nodes(:) == node(node1)) && any(Neumann_Nodes(:) == node(node2)) )
            for j=1:N_GQ
                eta=Points(j);
                weight=W(j);
                L=norm(p(node(node1),:)-p(node(node2),:),2);
                detJ=L/2;
                Ri=Ri-Neumann_Value*weight*Exterior_BC(eta)*detJ;
            end
            r(node(node1))=r(node(node1))+Ri(1);
            r(node(node2))=r(node(node2))+Ri(2);
        end
    end
    
end
disp('Imposing Dirichlet BCs...');
% Apply Dirichlet BCs:
K(Dirichlet_Nodes,:)=0;
K(:,Dirichlet_Nodes)=0;
K(Dirichlet_Nodes,Dirichlet_Nodes) = eye(length(Dirichlet_Nodes),length(Dirichlet_Nodes));
r(Dirichlet_Nodes)=Dirichlet_Values;

disp('Solving...');
u=K\r;
% u=GMRES(K,r,eye(Np),zeros(Np,1),5,1e-8,100);

figure;
trisurf(t,p(:,1),p(:,2),0*p(:,1),u,'edgecolor','k','facecolor','interp');
view(2); axis([-Xh Xh -1 1]); axis equal; colorbar;colormap(jet(256))
