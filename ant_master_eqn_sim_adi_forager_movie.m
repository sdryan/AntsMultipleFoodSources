clear
clc


L = 5; %Length of domain in 1D
dx = .1; % spatial step (assume square domain)
Nx = L/dx; % # of Spatial Steps in Each Dimension
dt = .01; % Time Step
Time = 20.0; % Total Time
Track = 0; %Assign Numbers to output images

% % Parameters
% alpha = .1; %diffusion of random foragers
% sigma = .5; %translation of foragers via pheromone
% omega_pq = 100.0; % transition prob from forager to returner
% omega_qp = 100.0; % transition prob from return to forager
% nu = .1; % advection of returners toward home
% D = .01; % diffusion of pheromone in the environment
% gamma = .01; % evaporation/degradation of pheromone
% A = 1.0; % pheromone deposited at food source


% Parameters
alpha = .1; %diffusion of random foragers
sigma = 100.5; %translation of foragers via pheromone
omega_pq = 100.0; % transition prob from forager to returner
omega_qp = 100.0; % transition prob from return to forager
nu = .15; % advection of returners toward home
D = .2; % diffusion of pheromone in the environment
gamma = .001; % evaporation/degradation of pheromone
A = .001; % pheromone deposited at food source

 count1 = 0;

CFL = alpha*dt/(dx*dx)
r = .5*alpha*dt/(dx*dx);
rs = sigma*dt/(8.0*dx*dx);
ru = .5*D*dt/(dx*dx);

xf(1) = 1.0; % x-coor of food source
xf(2) = 2.5; % y-coor of food source
xn(1) = 4.0; % x-coor of nest
xn(2) = 2.5; % y-coor of nest
% Setup x and y
 for i = 1:(Nx+1)
     for j = 1:(Nx+1)
         x1(i,j) = (i-1)*dx;
         y1(i,j) = (j-1)*dx;
         v_homex(i,j) = (xn(1)-x1(i,j))/sqrt((xn(1)-x1(i,j))^2+(xn(2)-y1(i,j))^2+.00000001); % x-coor Home Pointing Vector
         v_homey(i,j) = (xn(2)-y1(i,j))/sqrt((xn(1)-x1(i,j))^2+(xn(2)-y1(i,j))^2+.00000001); % y-coor Home Pointing Vector
     end
 end

% Initial Conditions ICs
for i = 1:(Nx+1)
    for j = 1:(Nx+1)
        p(i,j,1) = .1*exp(-sqrt((x1(i,j)-xn(1))^2+(y1(i,j)-xn(2))^2));%1/(L*L);
        q(i,j,1) = 0.0;%.1*exp(-sqrt((x1(i,j)-xf(1))^2+(y1(i,j)-xf(2))^2));%1/(L*L);
        u(i,j,1) = 0.0;
        p(i,j,2) = .1*exp(-sqrt((x1(i,j)-xn(1))^2+(y1(i,j)-xn(2))^2));%1/(L*L);
        q(i,j,2) = 0.0;%.1*exp(-sqrt((x1(i,j)-xf(1))^2+(y1(i,j)-xf(2))^2));%1/(L*L);
        u(i,j,2) = 0.0;
    end
end

Time1 = 0.0;

for k = 1:Time/dt
    Time1 = Time1 + dt;

    % Compute Next Timestep for Intermediate Points

    % Stage 1: n+1/2, implicit in x

    %p
        for j = 2:Nx
            for ii = 1:Nx-1
                a_p(ii) = 1.0+2.0*r;
                b_p(ii) = -r-rs*(u(ii+2,j,1)-u(ii,j,1));% b is left (i-1) -D*dt/(dx*dx)
                c_p(ii) = -r+rs*(u(ii+2,j,1)-u(ii,j,1));% c is right (i+1) -D*dt/(dx*dx)-c*dt/dx
                f_p(ii) = p(ii+1,j,1) + r*(p(ii+1,j+1,1)-2.0*p(ii+1,j,1)+p(ii+1,j-1,1))-sigma*dt*p(ii+1,j,1)/(2.0*dx*dx)*(u(ii+2,j,1)+u(ii,j,1)+u(ii+1,j+1,1)+u(ii+1,j-1,1)-4.0*u(ii+1,j,1))-rs*(u(ii+1,j+1,1)-u(ii+1,j-1,1))*(p(ii+1,j+1,1)-p(ii+1,j-1,1))-dt*omega_pq*p(ii+1,j,1)/(2.0)*delta1(x1(ii+1,j),y1(ii+1,j),xf(1),xf(2),dx)+dt*omega_qp*q(ii+1,j,1)/(2.0)*delta1(x1(ii+1,j),y1(ii+1,j),xn(1),xn(2),dx);%f is RHS 
            end
            f_p(1) = f_p(1)-0*b_p(1)*p(1,j,1);
            f_p(Nx-1) = f_p(Nx-1)-0*c_p(Nx-1)*p(Nx+1,j,1);
            p(2:Nx,j,2) = tridiag(a_p,b_p,c_p,f_p);
        end

    %q
    for i = 2:Nx
        for j = 2:Nx
            q(i,j,2) = q(i,j,1) - (nu*dt)/(4.0*dx)*(v_homex(i+1,j)*q(i+1,j,1)-v_homex(i-1,j)*q(i-1,j,1)+v_homey(i,j+1)*q(i,j+1,1)-v_homey(i,j-1)*q(i,j-1,1))+omega_pq*(dt/2.0)*delta1(x1(i,j),y1(i,j),xf(1),xf(2),dx)*p(i,j,1)-omega_qp*(dt/2.0)*delta1(x1(i,j),y1(i,j),xn(1),xn(2),dx)*q(i,j,1);
        end
    end
  
    %u
        for j = 2:Nx
            for ii = 1:Nx-1
                a_u(ii) = 1.0+2.0*ru;
                b_u(ii) = -ru;% b is left (i-1) -D*dt/(dx*dx)
                c_u(ii) = -ru;% c is right (i+1) -D*dt/(dx*dx)-c*dt/dx
                f_u(ii) = u(ii+1,j,1) + ru*(u(ii+1,j+1,1)-2.0*u(ii+1,j,1)+u(ii+1,j-1,1))-dt*gamma/(2.0)*u(ii+1,j,1)+dt*A/(2.0)*q(ii+1,j,1)*exp(-((x1(ii+1,j)-xf(1))^2+(y1(ii+1,j)-xf(2))^2));%f is RHS 
            end
            f_u(1) = f_u(1)-b_u(1)*u(1,j,1);
            f_u(Nx-1) = f_u(Nx-1)-c_u(Nx-1)*u(Nx+1,j,1);
            u(2:Nx,j,2) = tridiag(a_u,b_u,c_u,f_u);
        end

    % for i = 2:Nx
    %     for j = 2:Nx
    %         p(i,j,2) = p(i,j,1) + (alpha*dt)/(dx*dx)*(p(i+1,j,1)+p(i-1,j,1)+p(i,j+1,1)+p(i,j-1,1)-4.0*p(i,j,1))-(sigma*dt)/(dx*dx)*p(i,j,1)*(u(i+1,j,1)+u(i-1,j,1)+u(i,j+1,1)+u(i,j-1,1)-4.0*u(i,j,1))-(sigma*dt)/(4.0*dx*dx)*((u(i+1,j,1)-u(i-1,j,1))*(p(i+1,j,1)-p(i-1,j,1))+(u(i,j+1,1)-u(i,j-1,1))*(p(i,j+1,1)-p(i,j-1,1)))-omega_pq*dt*delta1(x1(i,j),y1(i,j),xf(1),xf(2),dx)*p(i,j,1)+omega_qp*dt*delta1(x1(i,j),y1(i,j),xn(1),xn(2),dx)*q(i,j,1);
    %         q(i,j,2) = q(i,j,1) - (nu*dt)/(2.0*dx)*(q(i,j,1)*(v_homex(i+1,j,1)-v_homex(i-1,j,1)+v_homey(i,j+1,1)-v_homey(i,j-1,1))+v_homex(i,j,1)*(q(i+1,j,1)-q(i-1,j,1))+v_homey(i,j,1)*(q(i,j+1,1)-q(i,j-1,1)))+omega_pq*dt*delta1(x1(i,j),y1(i,j),xf(1),xf(2),dx)*p(i,j,1)-omega_qp*dt*delta1(x1(i,j),y1(i,j),xn(1),xn(2),dx)*q(i,j,1);
    %         u(i,j,2) = u(i,j,1) + (D*dt)/(dx*dx)*(u(i+1,j,1)+u(i-1,j,1)+u(i,j+1,1)+u(i,j-1,1)-4.0*u(i,j,1))-dt*gamma*u(i,j,1)+dt*A*q(i,j,1)*exp(-((x1(i,j)-xf(1))^2+(y1(i,j)-xf(2))^2));
    %     end
    % end

        % BCs
    for i = 1:Nx+1
        p(i,1,2) = p(i,2,2);
        q(i,1,2) = q(i,2,2);
        u(i,1,2) = 0.0;
        p(i,Nx+1,2) = p(i,Nx,2);
        q(i,Nx+1,2) = q(i,Nx,2);
        u(i,Nx+1,2) = 0.0;       
        p(1,i,2) = p(2,i,2);
        q(1,i,2) = q(2,i,2);
        u(1,i,2) = 0.0;
        p(Nx+1,i,2) = p(Nx,i,2);
        q(Nx+1,i,2) = q(Nx,i,2);
        u(Nx+1,i,2) = 0.0; 
    end

    % Corners
     p(1,1,2) = p(2,2,2);
     q(1,1,2) = q(2,2,2);
     u(1,1,2) = u(2,2,2);
     p(Nx+1,1,2) = p(Nx,2,2);
     q(Nx+1,1,2) = q(Nx,2,2);
     u(Nx+1,1,2) = u(Nx,2,2);
     p(1,Nx+1,2) = p(2,Nx,2);
     q(1,Nx+1,2) = q(2,Nx,2);
     u(1,Nx+1,2) = u(2,Nx,2);
     p(Nx+1,Nx+1,2) = p(Nx,Nx,2);
     q(Nx+1,Nx+1,2) = q(Nx,Nx,2);
     u(Nx+1,Nx+1,2) = u(Nx,Nx,2);

% Enforce Prob Denisty Integrates to 1
      total_p = 0.0;
      for i = 2:Nx
          for j = 2:Nx
              total_p = total_p + (p(i,j,2)+q(i,j,2))*dx*dx/(16.0);
          end
      end


    % Update p, q, u
    for i = 1:(Nx+1)
        for j = 1:(Nx+1)
            %enforce positivity
            if (p(i,j,2) < 0)
                p(i,j,2) = -p(i,j,2);
            end
            if (q(i,j,2) < 0)
                q(i,j,2) = -q(i,j,2);
            end
            p(i,j,1) = p(i,j,2)/total_p;
            q(i,j,1) = q(i,j,2)/total_p;
            u(i,j,1) = u(i,j,2);
        end
    end

  % Stage 2: n+1, implicit in y

    %p
        for j = 2:Nx
            for ii = 1:Nx-1
                a_p(ii) = 1.0+2.0*r;
                b_p(ii) = -r-rs*(u(j,ii+2,1)-u(j,ii,1));% b is left (i-1) -D*dt/(dx*dx)
                c_p(ii) = -r+rs*(u(j,ii+2,1)-u(j,ii,1));% c is right (i+1) -D*dt/(dx*dx)-c*dt/dx
                f_p(ii) = p(j,ii+1,1) + r*(p(j+1,ii+1,1)-2.0*p(j,ii+1,1)+p(j-1,ii+1,1))-sigma*dt*p(j,ii+1,1)/(2.0*dx*dx)*(u(j,ii+2,1)+u(j,ii,1)+u(j+1,ii+1,1)+u(j-1,ii+1,1)-4.0*u(j,ii+1,1))-rs*(u(j+1,ii+1,1)-u(j-1,ii+1,1))*(p(j+1,ii+1,1)-p(j-1,ii+1,1))-dt*omega_pq*p(j,ii+1,1)/(2.0)*delta1(x1(j,ii+1),y1(j,ii+1),xf(1),xf(2),dx)+dt*omega_qp*q(j,ii+1,1)/(2.0)*delta1(x1(j,ii+1),y1(j,ii+1),xn(1),xn(2),dx);%f is RHS 
            end
            f_p(1) = f_p(1)-b_p(1)*p(j,1,1);
            f_p(Nx-1) = f_p(Nx-1)-c_p(Nx-1)*p(j,Nx+1,1);
            p(j,2:Nx,2) = tridiag(a_p,b_p,c_p,f_p);
        end

      %q
    for i = 2:Nx
        for j = 2:Nx
            q(i,j,2) = q(i,j,1) - (nu*dt)/(4.0*dx)*(v_homex(i+1,j)*q(i+1,j,1)-v_homex(i-1,j)*q(i-1,j,1)+v_homey(i,j+1)*q(i,j+1,1)-v_homey(i,j-1)*q(i,j-1,1))+omega_pq*(dt/2.0)*delta1(x1(i,j),y1(i,j),xf(1),xf(2),dx)*p(i,j,1)-omega_qp*(dt/2.0)*delta1(x1(i,j),y1(i,j),xn(1),xn(2),dx)*q(i,j,1);
        end
    end

        %u
        for j = 2:Nx
            for ii = 1:Nx-1
                a_u(ii) = 1.0+2.0*ru;
                b_u(ii) = -ru;% b is left (i-1) -D*dt/(dx*dx)
                c_u(ii) = -ru;% c is right (i+1) -D*dt/(dx*dx)-c*dt/dx
                f_u(ii) = u(j,ii+1,1) + ru*(u(j+1,ii+1,1)-2.0*u(j,ii+1,1)+u(j-1,ii+1,1))-dt*gamma/(2.0)*u(j,ii+1,1)+dt*A/(2.0)*q(j,ii+1,1)*exp(-((x1(j,ii+1)-xf(1))^2+(y1(j,ii+1)-xf(2))^2));%f is RHS 
            end
            f_u(1) = f_u(1)-b_u(1)*u(j,1,1);
            f_u(Nx-1) = f_u(Nx-1)-c_u(Nx-1)*u(j,Nx+1,1);
            u(j,2:Nx,2) = tridiag(a_u,b_u,c_u,f_u);
        end

        % BCs
    for i = 1:Nx+1
        p(i,1,2) = p(i,2,2);
        q(i,1,2) = q(i,2,2);
        u(i,1,2) = 0.0;
        p(i,Nx+1,2) = p(i,Nx,2);
        q(i,Nx+1,2) = q(i,Nx,2);
        u(i,Nx+1,2) = 0.0;       
        p(1,i,2) = p(2,i,2);
        q(1,i,2) = q(2,i,2);
        u(1,i,2) = 0.0;
        p(Nx+1,i,2) = p(Nx,i,2);
        q(Nx+1,i,2) = q(Nx,i,2);
        u(Nx+1,i,2) = 0.0; 
    end

    % Corners
     p(1,1,2) = p(2,2,2);
     q(1,1,2) = q(2,2,2);
     u(1,1,2) = u(2,2,2);
     p(Nx+1,1,2) = p(Nx,2,2);
     q(Nx+1,1,2) = q(Nx,2,2);
     u(Nx+1,1,2) = u(Nx,2,2);
     p(1,Nx+1,2) = p(2,Nx,2);
     q(1,Nx+1,2) = q(2,Nx,2);
     u(1,Nx+1,2) = u(2,Nx,2);
     p(Nx+1,Nx+1,2) = p(Nx,Nx,2);
     q(Nx+1,Nx+1,2) = q(Nx,Nx,2);
     u(Nx+1,Nx+1,2) = u(Nx,Nx,2);

% Enforce Prob Denisty Integrates to 1
      total_p = 0.0;
      for i = 2:Nx
          for j = 2:Nx
              total_p = total_p + (p(i,j,2)+q(i,j,2))*dx*dx/(16.0);
          end
      end
    % Update p, q, u
    for i = 1:(Nx+1)
        for j = 1:(Nx+1)
            %enforce positivity
            if (p(i,j,2) < 0)
                p(i,j,2) = -p(i,j,2);
            end
            if (q(i,j,2) < 0)
                q(i,j,2) = -q(i,j,2);
            end
            p(i,j,1) = p(i,j,2)/total_p;
            q(i,j,1) = q(i,j,2)/total_p;
            u(i,j,1) = u(i,j,2);
        end
    end

    count1 = count1 + 1;
    if (count1 > 9)
        count1 = 0;
        Time1
        clf()
        hold on
subplot(2,2,1)

%plot(xn(1),xn(2),'ko','MarkerSize',3)
imagesc(p(:,:,2))
colorbar
title('Forager Density p')
ax = gca; 
ax.FontSize = 16; 

subplot(2,2,2)

imagesc(q(:,:,2))
colorbar
title('Returner Density q')
ax = gca; 
ax.FontSize = 16; 

subplot(2,2,3)

imagesc(u(:,:,2))
colorbar 
title('Pheromone Density u')
ax = gca; 
ax.FontSize = 16; 

drawnow

       Track = Track + 1;
out=sprintf('MasterEqn_mov_%3.3d.png',Track);
filename=sprintf('MasterEqn_mov_%3.3d.png',Track);
saveas(gcf,filename)
    end
end


function [d1] = delta1(a,b,c,d,dx)

if (abs(a-c) <= dx*1.0  && abs(b-d) <= dx*1.0)
    d1 = 1.0;
else
    d1 = 0.0;
end


end

%Tridiagonal Solver With Thomas Algorithm
function y = tridiag( a, b, c, f )

%  Solve the  n x n  tridiagonal system for y:
%
%  [ a(1)  c(1)                                  ] [  y(1)  ]   [  f(1)  ]
%  [ b(2)  a(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
%  [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        b(n-1) a(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
%  [                                 b(n)  a(n)  ] [  y(n)  ]   [  f(n)  ]
%
%  f must be a vector (row or column) of length n
%  a, b, c must be vectors of length n (note that b(1) and c(n) are not used)

% some additional information is at the end of the file

n = length(f);
v = zeros(n,1);   
y = v;
w = a(1);
y(1) = f(1)/w;
for i=2:n
    v(i-1) = c(i-1)/w;
    w = a(i) - b(i)*v(i-1);
    y(i) = ( f(i) - b(i)*y(i-1) )/w;
end
for j=n-1:-1:1
   y(j) = y(j) - v(j)*y(j+1);
end
end


