%FEM code to solve the transient 3D Heat equation
%10.4167*dT/dt = [(d^2)T/dx^2 +(d^2)T/dy^2 +(d^2)T/dz^2]
%boundary conditions
%T(0,y,z,t)=25, T(50,y,z,t)=40, dT/dy(x,0,z,t)=0 ,dT/dy(x,6,t)=50(T-Tf),
%dT/dz(x,y,0,t)=0, dT/dz(x,y,1,t)=50(T-Tf)
%Initial condition: T(x,y,0)=25 over the domain.
%using hexahedral elements and forward difference method
%Variable descriptions
%K=element matrix for spatial term (T,xx +T,yy)
%M=element matrix for time-dependent term (T,t)
%F - element forcing vector
%KK -system matrix of K
%MM- system matrix of M
%FF- system forcing vector
%FN - effective system vector
%KN - effective system matrix
%fsol - Solution vector
%sol - time hiostory solution of selected nodes
%gcoord - Coordinate values of each node
%nodes -nodal conectivity of each element
%index - a vector containing system dofs associated with each element
%bcdof-  a vector containing dofs associated with boundary conditions
%bcval - a vector containing boundary condition values associated with  the dofs in bcdof 
%K1- element matrix due to Cauchy-type flux
%F1 - element vector due to flux boundary condition
%index1 - index for nodal dofs with flux


%input data for control parameters
Nel = 10;                            %number of elements
NNel = 8;                            %number of nodes per element
Ndof = 1;                            %number of dofs per node
Nnode = 44;                          %total number of nodes in the system
sdof=Nnode*Ndof;                     %total system dofs
delta_t=1;                         %time step for transient analysis
in_time= 0;                          %initial time
end_time = 1200;                       %termination time
Ntime = fix((end_time-in_time)/delta_t);    %number of time increment
k=10.4167;                              %thermal coefficient
bd1=50;
bd2=50;
Lx=50;
x_div=5;
Ly=6;
y_div=6;
Lz=1;
z_div=1;
nf=10;                                 %number of element boundaries with flux
NNels=4;                             %number of nodes per side of each element
h=500;                                 %convection coefficient
f=25;                                 %ambient temperature

%input data for nodal coordinate values

[nodes,connectivity]=NodalValue3D10elements(Lx,Ly,Lz,x_div,y_div,z_div);

 %input data for boundary conditions
 bcdof(1)=1;                             %first node is constrained
 bcval(1)=bd1;                            %whose described value is 25
 bcdof(2)=11;                             %11th node is constrained
 bcval(2)=bd2;                            %whose described value is 50
 bcdof(3)=12;                             %12th node is constrained
 bcval(3)=bd1;                            %whose described value is 25
bcdof(4)=22;                             %22th node is constrained
 bcval(4)=bd2;                            %whose described value is 50
bcdof(5)=23;                             %23rd node is constrained
 bcval(5)=bd1;                            %whose described value is 25
bcdof(6)=33;                             %33rd node is constrained
 bcval(6)=bd2;                            %whose described value is 50
 bcdof(6)=34;                             %34thd node is constrained
 bcval(6)=bd1;                            %whose described value is 50
 bcdof(6)=44;                             %44th node is constrained
 bcval(6)=bd2;                            %whose described value is 50

 %input for flux boundary conditions
 %nflx(i,j) i- element no: j- four side nodes
 nflx(1,1)=23; nflx(1,2)=24;nflx(1,3)=35;nflx(1,4)=34;
 nflx(2,1)=24; nflx(2,2)=25;nflx(2,3)=36;nflx(2,4)=35;
 nflx(3,1)=25; nflx(3,2)=26;nflx(3,3)=37;nflx(3,4)=36;
 nflx(4,1)=26; nflx(4,2)=27;nflx(4,3)=38;nflx(4,4)=37;
 nflx(5,1)=27; nflx(5,2)=28;nflx(5,3)=39;nflx(5,4)=38;
 nflx(6,1)=28; nflx(6,2)=29;nflx(6,3)=40;nflx(6,4)=39;
 nflx(7,1)=29; nflx(7,2)=30;nflx(7,3)=41;nflx(7,4)=40;
 nflx(8,1)=30; nflx(8,2)=31;nflx(8,3)=42;nflx(8,4)=41;
 nflx(9,1)=31; nflx(9,2)=32;nflx(9,3)=43;nflx(9,4)=42;
 nflx(10,1)=32;nflx(10,2)=33;nflx(10,3)=44;nflx(10,4)=43;

%initialization of matrices and vectors
FF=zeros(sdof,1);                        %initialization of system vector
FN=zeros(sdof,1);                        %initialization of effective system vector
fsol=zeros(sdof,1);                      %solution vector
sol=zeros(1,Ntime+1);                    %vector containing time history solution
KK =zeros(sdof,sdof);                    %initialization of system matrix
MM=zeros(sdof,sdof);                     %initialization of system matrix
KN =zeros(sdof,sdof);                    %effective system matrix
index =zeros(NNel*Ndof,1);               %initialization of index vector


%computation of element matrices,vectors and their assembly
for iel =1:Nel
    nd(1) = connectivity(iel,1);
    nd(2) =connectivity(iel,2);
    nd(3) =connectivity(iel,3);
    nd(4) =connectivity(iel,4);
    nd(5) =connectivity(iel,5);
    nd(6) =connectivity(iel,6);
    nd(7) =connectivity(iel,7);
    nd(8) =connectivity(iel,8);
    x1=nodes(nd(1),1);
    y1=nodes(nd(1),2);
    z1=nodes(nd(1),3);
    x2=nodes(nd(2),1);
    y2=nodes(nd(2),2);
    z2=nodes(nd(2),3);
    x3=nodes(nd(3),1);
    y3=nodes(nd(3),2);
    z3=nodes(nd(3),3);
    x4=nodes(nd(4),1);
    y4=nodes(nd(4),2);
    z4=nodes(nd(4),3);
    x5=nodes(nd(5),1);
    y5=nodes(nd(5),2);
    z5=nodes(nd(5),3);
    x6=nodes(nd(6),1);
    y6=nodes(nd(6),2);
    z6=nodes(nd(6),3);
    x7=nodes(nd(7),1);
    y7=nodes(nd(7),2);
    z7=nodes(nd(7),3);
    x8=nodes(nd(8),1);
    y8=nodes(nd(8),2);
    z8=nodes(nd(8),3);
    a=(x2-x1)/2;
    b=(y5-y1)/2;
    c=(z4-z1)/2;
    

    index =Systemdofs(nd,NNel,Ndof);
    K =HexahedralElementStiffnessMatrix1(a,b,c);
    KK=Assembly1(KK,K,index);
    M=k*HexahedralElementCapacitanceMatrix1(a,b,c);
    MM=Assembly2(MM,M,index);
    
end

%additional computation due to flux boundary condition
for ifx=1:nf
    nds(1)=nflx(ifx,1);
    nds(2)=nflx(ifx,2);
    nds(3)=nflx(ifx,3);
    nds(4)=nflx(ifx,4);
    x1=nodes(nds(1),1);
    y1=nodes(nds(1),2);
    x2=nodes(nds(2),1);
    y2=nodes(nds(2),2);
    eleng=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    
    index1=Fluxelementdof(nds,NNels,Ndof);
    K1=h*CauchyFlux(eleng);
    F1=h*f*fef1l(eleng);
    [KK,FF]=Assembly3(KK,FF,K1,F1,index1);
    
end


% %loop for time integration
for in = 1:sdof
    fsol(in)=25;
end
sol(1)=fsol(20);
KN=2*MM+delta_t*KK;
FN=delta_t*FF+((2*MM)-(delta_t*KK))*fsol;
    [KN,FN]=ApplyBC(KN,FN,bcdof,bcval);
    fsol=KN\FN;
for it=1:Ntime
    
    sol(it+1)=fsol(20);
end

% %plot the solution at node 33
time =0:delta_t:Ntime*delta_t;
plot(time,sol);
xlabel('Time');
ylabel('TemperatureSolution')
    