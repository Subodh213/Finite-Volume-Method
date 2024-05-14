tStart = cputime;
%  Lid Driven Cavity Parameters
nx = 32;
ny = 32;
nu = 0.01; %  viscosity
lx = 1;
ly = 1;
dx = lx/nx;
dy = ly/ny;
%% Boundary conditions
uT = 1.0;
uB = 0.0;
vL = 0.0;
vR = 0.0;
% Reynolds no.
Re = uT*lx/nu;
% dt based on linear advection-diffusion 
dt_max = min([dx/uT , 0.25*dx*dx/nu]);
%% CFL condition
% dt = dt_max;
dt = 6.25E-3;%Added 
CFL_no = uT*dt/dx;


%% Staggereg Grid 

u = zeros(ny+2,nx+2); %preallocating
v = zeros(ny+2,nx+2); %preallocating

%setting boundary conditions
% Set bcs on u velocity grid
    % left wall
    u(2:end-1,2) = 0.0;
    % right wall
    u(2:end-1,end) = 0.0;
    % top wall -ghost cell
    u(1,2:end-1) = 2.0*uT - u(2,2:end-1);
    % bottom wall- ghost cell
    u(end,2:end-1) = 2.0*uB - u(end-1,2:end-1);

% Set bcs on v velocity grid
    % top wall
    v(1,2:end-1) = 0.0;
    % bottom wall
    v(end-1,2:end-1) = 0.0;
    % left wall-ghost cell
    v(1:end-1,1) = 2.0*vL - v(1:end-1,2);
    % right wall-ghost cell
    v(2:end-1,end) = 2.0*vR - v(2:end-1,end-1);

%%
%For calculating L2 error norm
u_old = u;
v_old = v;

delta_u = zeros(ny+2,nx+2);
delta_v = zeros(ny+2,nx+2);

l2_norm_u = [];
l2_norm_v = [];

error_u = 1;
error_v = 1;
%%

%preallocating intermidate velocity (before pressure correction)
u_tilde = zeros(ny+2,nx+2);
v_tilde = zeros(ny+2,nx+2);

%preallocating convection and diffusion for saving n-1 time step values
%for Adam Bashforth Method
temp_Xcon = zeros(ny+2,nx+2);
temp_Xdiff = zeros(ny+2,nx+2);
temp_Ycon = zeros(ny+2,nx+2);
temp_Ydiff = zeros(ny+2,nx+2);

%preallocating pressure
p = zeros(ny+2,nx+2); 

%preallocating convection and diffusion term
con = zeros(1,1);
diff = zeros(1,1);

%% Assemble pressure coefficients using RHS of Poission Equation
%∇2 P = 
Ae = 1/dx^2*ones(ny+2,nx+2);
Aw = 1/dx^2*ones(ny+2,nx+2);
An = 1/dy^2*ones(ny+2,nx+2);
As = 1/dy^2*ones(ny+2,nx+2);
% set the left wall coefficients
Aw(2:end-1,2) = 0.0;
% set the right wall coefficients
Ae(2:end-1,end-1) = 0.0;
% set the top wall coefficients
An(2,2:end-1) = 0.0;
% set the bottom wall coefficients
As(end-1,2:end-1) = 0.0;
Ap = -(Aw + Ae + An + As);

%% Time integration
n = 1;
while error_u > 1e-8 && error_v > 1e-8

    
%% Assemble x-mom first for u_tilde -  interior points 
    for i = 3:nx+1
        for j = 2:ny+1
            ue = 0.5*(u(j,i)+u(j,i+1));
            uw = 0.5*(u(j,i)+u(j,i-1));
            un = 0.5*(u(j,i)+u(j+1,i));
            us = 0.5*(u(j,i)+u(j-1,i));
            vn = 0.5*(v(j+1,i-1)+v(j+1,i));
            vs = 0.5*(v(j,i)+v(j,i-1));

            
            con(j,i) = -(ue*ue - uw*uw)/dx - (un*vn - us*vs)/dy; %Convection
            diff(j,i) = (1/Re)*((u(j,i+1) - 2*u(j,i) + u(j,i-1))/dx^2 + ...
                        (u(j+1,i) - 2*u(j,i) + u(j-1,i))/dy^2) ; %Diffusion
            %Hi = Convection + Diffusion
            % u_tilde  ̃= u_i^n+ ∫H_i dt 
            %∫H_i dt= Δt*[(3/2 * (H_i^n)) - (1/2 * H_i^(n-1))]  Adam Bashforth
            
            u_tilde(j,i) = u(j,i) + dt*(((3/2)*(con(j,i) + diff(j,i))) - ((1/2)*(temp_Xcon(j,i) + temp_Xdiff(j,i)))); %  Adams Bashforth

            temp_Xcon(j,i) = con(j,i); %saving convection term for (n-1) time 
            temp_Xdiff(j,i) = diff(j,i);%saving convection term for (n-1) time
        end
    end

    %% Assemble y-mom first for v_tilde - do interior points only
    for i = 2:nx+1
        for j = 2:ny
            ue = 0.5*(u(j,i+1) + u(j-1,i+1));
            uw = 0.5*(u(j-1,i) + u(j,i));
            ve = 0.5*(v(j,i) + v(j,i+1));
            vw = 0.5*(v(j,i) + v(j,i-1));
            vn = 0.5*(v(j,i) + v(j+1,i));
            vs = 0.5*(v(j,i) + v(j-1,i));

            con(j,i) = -(ve*ue - vw*uw)/dx - (vn*vn - vs*vs)/dy;
            diff(j,i) = (1/Re)*((v(j,i+1) - 2*v(j,i) + v(j,i-1))/dx^2 + ...
                        (v(j+1,i) - 2*v(j,i) + v(j-1,i))/dy^2) ;

            %Hi = Convection + Diffusion
            % u_tilde  ̃= u_i^n+ ∫H_i dt 
            %∫H_i dt= Δt*[(3/2 * (H_i^n)) - (1/2 * H_i^(n-1))]  Adam Bashforth

            v_tilde(j,i) = v(j,i) + dt/2*(3*(con(j,i) + diff(j,i)) - (temp_Ycon(j,i) + temp_Ydiff(j,i))); % 2nd order Adams Bashforth

            temp_Ycon(j,i) = con(j,i); %saving convection term for (n-1) time 
            temp_Ydiff(j,i) = diff(j,i); %saving convection term for (n-1) time
        end
    end

    %Pressure Correction
    %D_i*(G_i  P^(n+1))=  1/∆t  D_i* (u tilde_i ) ̃
    %% rhs: prhs = 1/∆t  D_i* (u tilde_i ) =  1/dt * (du_tilde/dx + dv_tilde/dy)
    divut = zeros(ny+2,nx+2);
    divut(2:end-1,2:end-1) = (u_tilde(2:end-1,3:end) - u_tilde(2:end-1,2:end-1))/dx + (v_tilde(2:end-1,2:end-1) - v_tilde(1:end-2,2:end-1))/dy;
    prhs = divut/dt;
    
    %% Solve pressure poisson equation
    p = sor_solver(p,Ae,Aw,An,As,Ap,prhs,nx,ny);
    
    % Final time advance - do interior only
    % u = ut - dt*dpdx
    u(2:end-1,3:end-1) = u_tilde(2:end-1,3:end-1) - dt*(p(2:end-1,3:end-1)-p(2:end-1,2:end-2))/dx; % (n+1)th time step
    v(2:end-2,2:end-1) = v_tilde(2:end-2,2:end-1) - dt*(p(3:end-1,2:end-1)-p(2:end-2,2:end-1))/dy; % (n+1)th time step
    delta_u = u - u_old;
    delta_v = v - v_old;


%Calculating and Appending L2_error_norms
    l2_norm_u = [l2_norm_u,sqrt(sum(delta_u.^2,"all")/(nx*ny))];
    l2_norm_v = [l2_norm_v,sqrt(sum(delta_v.^2,"all")/(nx*ny))];
    error_u = l2_norm_u(end);
    error_v = l2_norm_v(end);
    fprintf('\t%d \t%f \t%f\n',n,l2_norm_u,l2_norm_v);
    n = n+1;

    %Update U_old for new iteration
    u_old = u;
    v_old = v;
end
time = dt*linspace(2,n,n-1);
figure(1)
% xlim([0,80])
% ylim([1e-8,1])
% subplot(2,1,1)
title('Re = 100  on a 32*32 finite volume grid')
semilogy(time,l2_norm_u,'^');
xlabel('non-dimensional time units')
ylabel('average l2 norm of u vel')
hold on

tEnd = cputime - tStart;
%%
%Load Ghia&Ghia data
T = readtable('GhiaData.csv');
x_coord = T.x_coord;  
y_coord = T.y_coord;  
u_ghia_Re100 = T.u_ghia_Re100;
v_ghia_Re100 = T.v_ghia_Re100;

%%

figure(2)
title('u velocity along vertical line through geometric centre of cavity')
plot(y_coord,u_ghia_Re100)
hold on
% % plot(linspace(0,1,66),u(end:-1:1,34),'--')
% % hold on
% plot(linspace(0,1,130),u(end:-1:1,66),'o')
% % hold on
 plot(linspace(0,1,34),u(end:-1:1,18),'d')

figure(3)
plot(x_coord,v_ghia_Re100)
hold on
% plot(linspace(0,1,66),v(33,end:-1:1),'--')
% hold on
% plot(linspace(0,1,130),v(65,end:-1:1),'o')
% % hold on
plot(linspace(0,1,34),v(17,end:-1:1),'d')
% %% n * n finite volumes
% y_grid_pt = round(1 + (1-y_coord)*(nx+1)); % at x = 0.5
% x_grid_pt = round(1 + x_coord*(nx+1)); % at y = 0.5
% 
% % Book-keeping of u and v velocities at above grid points
% u_midway = zeros(1,1);
% v_midway = zeros(1,1);
% for i = 1:length(y_grid_pt) 
%     u_midway(i,1) = u(y_grid_pt(i),nx/2+2);
% end
% for i = 1:length(x_grid_pt)
%     v_midway(i,1) = v(nx/2+1,x_grid_pt(i));
% end
% figure(2)
% plot(y_coord,u_ghia_Re100,y_coord,u_midway(:,1),'o',y_coord,u_midway(:,2),'s',y_coord,u_midway(:,3),'^')
% xlabel('y coordinate at x = 0.5')
% ylabel('u velocity')
% legend('Ghia paper Table I','32*32 grid,\delta t = 3.125e-3','64*64 grid,\delta t = 5e-4','128*128 grid,\delta t = 1.25e-4');
% figure(3)
% plot(x_coord,v_ghia_Re100,x_coord,v_midway(:,1),'o',x_coord,v_midway(:,2),'s',x_coord,v_midway(:,3),'^')
% legend('Ghia paper Table II','32*32 grid,\delta t = 3.125e-3','64*64 grid,\delta t = 5e-4','128*128 grid,\delta t = 1.25e-4');

%% Save data
% save('G32Dt1.mat', 'time', 'l2_norm_u', 'u', 'v', "x_coord","y_coord","u_ghia_Re100","v_ghia_Re100");
save('G32Dt2d.mat', 'time', 'l2_norm_u', 'u', 'v', "x_coord","y_coord","u_ghia_Re100","v_ghia_Re100","tEnd");
% save('G32Dt3.mat', 'time', 'l2_norm_u', 'u', 'v', "x_coord","y_coord","u_ghia_Re100","v_ghia_Re100");
%% SOR solver function
function phi = sor_solver(phi_0,Ae,Aw,An,As,Ap,rhs_vec,nx,ny)
   w = 2/(1+sin(pi/(nx+1)));
   phi = phi_0;
   r = zeros(nx,ny);
   residual = 1;
%    error = 1;
   tol = 1e-5;
%    residue_history = [];
%    error_history = [];
%    phi = zeros(nx+2,ny+2);
   while residual > tol
       for i = 2:nx+1
           for j = 2:ny+1
               phi(j,i) = w*((rhs_vec(j,i) - Ae(j,i)*phi(j,i+1) - Aw(j,i)*phi(j,i-1) - ...
                   An(j,i)*phi(j+1,i) - As(j,i)*phi(j-1,i))/Ap(j,i)) + (1-w)*phi(j,i);
           end
       end
       residual = 0;
       for i = 2:nx+1
           for j = 2:ny+1
               r(j,i) = rhs_vec(j,i) - Ae(j,i)*phi(j,i+1) - Aw(j,i)*phi(j,i-1) - ...
                   An(j,i)*phi(j+1,i) - As(j,i)*phi(j-1,i) - Ap(j,i)*phi(j,i);
               residual = residual + r(j,i)*r(j,i);
           end
       end
       residual = sqrt(residual/(nx*ny));
%        error = sqrt(sum((phi - phi_0).^2,"all")/(nx*ny));
%        residue_history = [residue_history,residual];
%        error_history = [error_history,sqrt(sum(error,'all')/M^2)];
%        disp(residual);
   end
end