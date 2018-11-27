clear
clc
addpath('./distmesh');
addpath('./LinearSolvers');

U_inlet=1;
N_GQ=2;
disp('Meshing...');

Xh=3;
Radius=1;
fd=@(p) ddiff(drectangle(p,-Xh,Xh,-Xh,Xh),dcircle(p,0,0,Radius));
fh=@(p) 0.02+0.0*dcircle(p,0,0,0.5);
[p,t]=distmesh2d(fd,fh,0.2,[-1.1*Xh,-1.1*Xh;1.1*Xh,1.1*Xh],[-Xh,-Xh;-Xh,Xh;Xh,-Xh;Xh,Xh]);

BC_nodes = unique(boundedges(p,t));
BC_points=p(BC_nodes,:);
R_nodes=BC_nodes(find(BC_points(:,1)>0.99*Xh));
L_nodes=BC_nodes(find(BC_points(:,1)<-0.99*Xh));
TB_nodes=BC_nodes(find(BC_points(:,2)>0.99*Xh));

rs=sqrt(BC_points(:,1).^2+BC_points(:,2).^2);
circle_nodes=BC_nodes(find(rs<Radius*1.1));
circle_points=p(circle_nodes,:);


Dirichlet_Nodes=R_nodes;
Inlet_nodes=L_nodes;

mu=0.001;
source =0;
v_i=[0,0];

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
%%% Element Interior
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
                Cy(node(k1),node(k2))=Cy(node(k1),node(k2))+Ciy(k1,k2);
                M(node(k1),node(k2))=M(node(k1),node(k2))+Mi(k1,k2);
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
    Num_edge=Num_nodes;
    
    if (size(intersect(node,Inlet_nodes),1)<2)
        continue;
    end
    
    %%%loop over edges of this element
    for edge=1:Num_edge
        Ri=zeros(2,1);%since this is a linear element 2 nodes per edge
        node1=edge;
        node2=edge+1;
        if (edge==Num_edge); node2=1; end
        %%% check to see if this edge is on the boundary
        if(any(Inlet_nodes(:) == node(node1)) && any(Inlet_nodes(:) == node(node2)) )
            for j=1:N_GQ
                eta=Points(j);
                weight=W(j);
                L=norm(p(node(node1),:)-p(node(node2),:),2);
                detJ=L/2;
                Ri=Ri-mu*U_inlet*weight*Exterior_BC(eta)*detJ;
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
r(Dirichlet_Nodes)=0.0;


disp('Solving...');
u=K\r;
% u=GMRES(K,r,eye(Np),zeros(Np,1),5,1e-8,100);




figure;
trisurf(t,p(:,1),p(:,2),0*p(:,1),u,'edgecolor','k','facecolor','interp');
view(2); axis([-Xh Xh -Xh Xh]); axis equal; colorbar;colormap(jet(256))
title('phi')
shading interp


vx=inv(M)*Cx*u;
vy=inv(M)*Cy*u;
v=sqrt(vx.^2+vy.^2);

pressure=-1/2*1*(v.^2-U_inlet^2);


figure;
trisurf(t,p(:,1),p(:,2),0*p(:,1),vx,'edgecolor','k','facecolor','interp');
view(2); axis([-Xh Xh -Xh Xh]); axis equal; colorbar;colormap(jet(256))
shading interp
title('Ux')

figure;
trisurf(t,p(:,1),p(:,2),0*p(:,1),vy,'facecolor','interp');
view(2); axis([-Xh Xh -Xh Xh]); axis equal; colorbar;colormap(jet(256))
shading interp
title('Uy')


figure;
trisurf(t,p(:,1),p(:,2),0*p(:,1),v,'facecolor','interp');
view(2); axis([-Xh Xh -Xh Xh]); axis equal; colorbar;colormap(jet(256))
shading interp
title('|U|')


figure;
trisurf(t,p(:,1),p(:,2),0*p(:,1),pressure,'facecolor','interp');
view(2); axis([-Xh Xh -Xh Xh]); axis equal; colorbar;colormap(jet(256))
shading interp
title('Pressure')

% [x,y]=meshgrid(p(:,1),p(:,2));
% [u,v]=meshgrid(vx,vy);
% 
% starty = p(1:200:end,2);
% startx = p(1:200:end,1);
% streamline(x,y,u,v,startx,starty)

