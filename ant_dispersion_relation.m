clear;

%here we will write code that generates a dispersion curve relating the
%eigenvalue of the linearized operator with the frequency of the
%perturbation.

%In this script, let x be the largest eigenvalue of the spectrum, and let z
%represent the first coordinate of the wavevector k. Fix the second
%coordinate of k to some value. This is ok because the dispersion relation
%will be symmetric.

%diffusion of ants
alpha = 10;

%sensitivity of ants to pheromone
sigma = 100;

%lattice size
M = 5;
N = 5;

%transition propensities
Opq = 100;
Oqp = 1;

%second coordinate of wavevector
k = 1;

%speed of returning ants
for nu = [0.1 1 5 10]
%nu = 1;

%food location
xf = [4,3];

%nest location
y = [M/2,N/2];

%unit vector connecting food to nest
u = (y-xf)/(norm(y-xf));

%term emerging from averging equations
E = (erf(xf(1)) - erf(xf(1)-N))*(erf(xf(2))-erf(xf(2)- M));

%amplitude of pheromone release
%for A = [5,10,15,20,25]
A = 15;

%pheromone diffusion coefficient
D = 25;

%evaporation rate
%for gamma = [0.1 1 10 50 100]
gamma = 1;

%setting up matrix
syms x z
B(1,1) = x + alpha*(k^2 + z^2) + Opq/(M*N);
B(1,2) = -Oqp/(M*N);
B(1,3) = -(sigma/(M*N))*(k^2 + z^2);
B(2,1) = -Opq/(M*N);
B(2,2) = x - 1i*nu*(u*[z;k]) + Oqp/(M*N);
B(2,3) = 0;
B(3,1) = 0;
B(3,2) = -(A*pi*E)/(4*M*N);
B(3,3) = x + gamma + D*(k^2 + z^2);

%determinant
C = det(B);
Z = 0.001:0.001:2;
lz = length(Z);
eigen = zeros(lz,1);

for j = 1:lz
    z = Z(j);
    C1 = subs(C);
    F = real(vpasolve(C1==0));
    F1 = max(F);
    eigen(j) = F1;
end

%plotting figure
figure(1)
plot(Z,eigen,'LineWidth',4)
h = legend('\nu = 0.1','\nu = 1','\nu = 5','\nu = 10');
set(h,'box','off')
set(gca,'fontsize',20)
xlabel('k')
ylabel('Re[\lambda]')
hold on
end

