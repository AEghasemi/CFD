%--------------------------------------------------------------------------
%                 High order  WENO-Explicit Runge Kutta schemes           -
%                                                                         -
%             By: Amirehsan Ghasemi   Email: aghasemi@vols.utk.edu        -
%--------------------------------------------------------------------------

function main
clc
clear all

length_grid = 5;
x  = linspace(0, length_grid, 400);
dx = x(2)-x(1);
z  = length(x);

r = 2;               % WENO Order ==> (2*r-1) choose from r = 1 or r= 2 or r= 3
order_erk = 6;       % Explicit Runge Kutta Order 
epsilon   = 10^-6;   % 10^-5 <= epsilon <= 10^-7
p = r;               % same as r

%%% Try 1: Use function
% u_n = zeros(z,1);    % u_n at(t = 0)
% for i = 1:z
%     u_n(i) = sin(2*pi*x(i))+2;
% end
% BCS at the beginning u_n(1)
% BCS at the end       u_n(end)

%%% Try 2: Use function
u_n = create_u_n(z, 3);

% -------------------------  Time Steps
dt = 1.e-3;
nn = 1.e6;

% -------------------------  Time Marching

for i = 1:nn
    u = erk(u_n,order_erk,dt,dx,r,epsilon,p);
    u_n = u;
    u_n(1) = 3;          % BCS at the beginning: u_n(1)
    u_n(z) = u(z-1);     % BCS at the end      : u_n(end)
    figure(1)            % Plot
    plot(x(2:z-1),u(2:z-1))
    xlabel('x')
    ylabel('u')
    axis([0 length_grid 0 5 ])     
end


end


%##########################################################################
%                       Explicit Runge-Kutta
%##########################################################################

function u = erk(u_n, order_erk, dt, dx, r, epsilon, p)

z = length(u_n);                       % Size(u_n)

s = stages_erk(order_erk);             % Finding Number of Stages
[a,b,~] = coeffs_erk(s,order_erk);     % Call Coefficients of ERK                        
k = zeros(z,1,s);                      % Preallocate

for i = 1:s
    tmp = 0;
    for  j = 1:(i-1)
        tmp = tmp + a(i,j)*k(:,:,j);
    end
    k(:,:,i) = (-1)*get_res(u_n+(dt*tmp),dx,r,epsilon,p);
end

tmp = 0;
for i = 1:s
    tmp = tmp + b(i)*k(:,:,i);
end

u = u_n+(dt*tmp);
end


function out = get_res(u_k,dx,r,epsilon,p)
f = flux(u_k);
out = weno_diff(f,r,epsilon,p,dx);
end


function out = flux(u_k)
z   = length(u_k);
out = zeros(z,1);
for  i = 1:z
    out(i) = (u_k(i))^2/(2);
end
end


%##########################################################################
%                    Compute d(f)/dx Using WENO Scheme
%##########################################################################

function dfdx = weno_diff(f,r,epsilon,p,dx)

z = length(f);

rj = zeros(z,1);             % Preallocate
rj(:,:)   = r;               % Fixing rj = r  
rj(1)     = 1;               % at j = 1
rj(2)     = 2;               % at j = 2
rj(end-1) = 2;               % at j = end-1 
rj(end)   = 1;               % at j = end

fph  = zeros(z-1,1);         % Preallocate 
dfdx = zeros(z,1);           % Preallocate  dfdx(1) = 0 , dfdx(end) = 0


for j = 1:(z-1)
    r = rj(j);
    tmp = 0;
    for k = 0:(r-1)
        for L = 0:(r-1)
            tmp = tmp + (Omegak(f,k,r,j,epsilon,p)*a_rkL(r,k+1,L+1)*f(j+k-r+1+L));
        end
    end
    fph(j) = tmp;
end


dfdx(1) = 0;
dfdx(z) = 0;
for j = 2:(z-1)
    dfdx(j) = (fph(j)-fph(j-1)) / dx;
end

end


%##########################################################################
%                           Compute Omega k
%##########################################################################

function out = Omegak(f,k,r,j,epsilon,p)

tmp = 0;
for kk = 0:(r-1)
    tmp = tmp + alphak(f,kk,r,j,epsilon,p); 
end

out = alphak(f,k,r,j,epsilon,p)/tmp;
end


%##########################################################################
%                           Compute Alpha k 
%##########################################################################

function out = alphak(f,k,r,j,epsilon,p)

out = C_rk(r,k+1)/(epsilon+ISk2(f,k,r,j))^p;

end


%##########################################################################
%               New Smoothness Measurement: Compute ISk2 
%##########################################################################

function out = ISk2(f,k,r,j)


if     ( r == 1 || r == 2 )
    out = ISk(f,k,r,j);
  
elseif ( r == 3 )
    if    ( k == 0 )
        out = (13/12)*(f(j-2)-2*f(j-1)+f(j))^2 + (1/4)*(f(j-2)-4*f(j-1)+3*f(j))^2;
    elseif( k == 1 )
        out = (13/12)*(f(j-1)-2*f(j)+f(j+1))^2 + (1/4)*(f(j-1)-f(j+1))^2;
    elseif( k == 2 )
        out = (13/12)*(f(j)-2*f(j+1)+f(j+2))^2 + (1/4)*(3*f(j)-4*f(j+1)+f(j+2))^2;    
    end   
end    
  
end

%..........................................................................
%                               Compute ISk 
%..........................................................................

function out = ISk(f,k,r,j)
out = 0;

for L = 1:(r-1)
    for i = 1:(r-1)
        out = out + (ddf(f,(j+k+i-r),L))^2 /(r-L);
    end
end

end


%..........................................................................
%                    Compute ddf =>  inner f[j+k+i-r,L]
%..........................................................................

function out = ddf(f,idx1,idx2)

if ( idx2 == 0 )
    out = f(idx1);
else
    out = ddf(f,(idx1+1),(idx2-1)) - ddf(f,idx1,(idx2-1));
end

end


%##########################################################################
%##########################################################################
%                           WENO Coefficients
%##########################################################################
%##########################################################################

function out = C_rk(r,k)

C = zeros(3,3);     % Preallocate
C(1,1) = 1;
C(1,2) = 0;
C(1,3) = 0;

C(2,1) = (1/3);
C(2,2) = (2/3);
C(2,3) = (0);

C(3,1) = (1/10);
C(3,2) = (6/10);
C(3,3) = (3/10);

out = C(r,k);
end


function out = a_rkL(r,k,L)

if     ( r == 1 )
    a(1,1) = 1;
    
elseif ( r == 2 )
    a = zeros(r,r);
    a(1,1) = (-1/2);
    a(1,2) = (3/2);
    
    a(2,1) = (1/2);
    a(2,2) = (1/2);
    
elseif ( r == 3 )
    a = zeros(r,r);
    a(1,1) = (1/3);
    a(1,2) = (-7/6);
    a(1,3) = (11/6);
    
    a(2,1) = (-1/6);
    a(2,2) = (5/6);
    a(2,3) = (1/3);
    
    a(3,1) = (1/3);
    a(3,2) = (5/6);
    a(3,3) = (-1/6);
end

out = a(k,L);
end
 

%##########################################################################
%##########################################################################
%   Finding Number of Stages (s) and Coefficients: [aij] ; [bi] ; [ci] 
%##########################################################################
%##########################################################################

% Finding Number of Statges According To Order
function stage = stages_erk(order)
if     ( order == 1 )
    stage = 1;
elseif ( order == 2 )
    stage = 2;
elseif ( order == 3 )
    stage = 3;
elseif ( order == 4 )
    stage = 4;
elseif ( order == 5 )
    stage = 6;
elseif ( order == 6 )
    stage = 7;
elseif ( order == 7 )
    stage = 9;
elseif ( order == 8 )
    stage = 11;
end

end

%**********************************************
% Finding Coefficients: [aij] ; [bi] ; [ci]
function [a, b, c] = coeffs_erk(s, order)

% [aij] known as Runge–Kutta Matrix
% [bi]  known as Weights
% [ci]  known as Nodes
%**********************************************
% Finding Coefficient
if     ( order == 1)          %(forward) Euler method
    a = zeros(s,s);           % Preallocate
    
    b = zeros(1,s);           % Preallocate
    b(1,1) = 1;
    
    c = zeros(s,1);           % Preallocate
    
elseif ( order == 2)          % Midpoint Method
    a = zeros(s,s);           % Preallocate
    a(2,1) = (1/2);
    
    b = zeros(1,s);           % Preallocate
    b(1,2) = (1);
    
    c = zeros(s,1);           % Preallocate
    c(2,1) = (1/2);
                              
elseif ( order == 3)          % Kutta's Method
    a = zeros(s,s);           % Preallocate
    a(2,1) = (1/2);
  
    b = zeros(1,s);           % Preallocate
    b(1,1) = (1/6);
    b(1,2) = (2/3);
    b(1,3) = (1/6);
    
    c = zeros(s,1);           % Preallocate
    c(2,1) = (1/2);
    c(3,1) = (1);
                               
elseif ( order == 4)          % Original RK Method
    a = zeros(s,s);           % Preallocate
    a(2,1) = (1/2);
    
    a(3,2) = (1/2);
    
    a(4,3) = (1);
    
    b = zeros(1,s);           % Preallocate
    b(1,1) = (1/6);
    b(1,2) = (1/3);
    b(1,3) = (1/3);
    b(1,4) = (1/6);
    
    c = zeros(s,1);           % Preallocate
    c(2,1) = (1/2);
    c(3,1) = (1/2);
    c(4,1) = (1);

elseif ( order == 5)          % Nystrom Method
    a = zeros(s,s);           % Preallocate
    a(2,1) = (1/3);
    
    a(3,1) = (4/25);
    a(3,2) = (6/25);
    
    a(4,1) = (1/4);
    a(4,2) = (-3);
    a(4,3) = (15/4);
    
    a(5,1) = (2/27);
    a(5,2) = (10/9);
    a(5,3) = (-50/81);
    a(5,4) = (8/81);
    
    a(6,1) = (2/25);
    a(6,2) = (12/25);
    a(6,3) = (2/15);
    a(6,4) = (8/75);
    
    b = zeros(1,s);           % Preallocate
    b(1,1) = (23/192);
    b(1,3) = (125/192);
    b(1,5) = (-27/64);
    b(1,6) = (125/192);
    
    c = zeros(s,1);           % Preallocate
    c(2,1) = (1/3);
    c(3,1) = (2/5);
    c(4,1) = (1);    
    c(5,1) = (2/3); 
    c(6,1) = (4/5);
    
elseif ( order == 6 )         % Butcher's Book
    a = zeros(s,s);           % Preallocate
    a(2,1) = (2/5);
        
    a(3,2) = (4/5);
        
    a(4,1) = (169/1458);
    a(4,2) = (110/729);
    a(4,3) = (-65/1458);
        
    a(5,1) = (-44/675);
    a(5,2) = (-88/135);
    a(5,3) = (76/351);
    a(5,4) = (336/325);
        
    a(6,1) = (21/106);
    a(6,3) = (-105/689);
    a(6,4) = (-324/689);
    a(6,5) = (45/106);
        
    a(7,1) = (-2517/4864);
    a(7,2) = (-55/38);
    a(7,3) = (10615/31616);
    a(7,4) = (567/7904);
    a(7,5) = (7245/4864);
    a(7,6) = (2597/2432);
        
    b = zeros(1,s);           % Preallocate
    b(1,3) = (1375/4992);
    b(1,4) = (6561/20384);
    b(1,5) = (3375/12544);
    b(1,6) = (53/768);
    b(1,7) = (19/294);
        
    c = zeros(s,1);           % Preallocate
    c(2,1) = (2/5);
    c(3,1) = (4/5);
    c(4,1) = (2/9);
    c(5,1) = (8/15);
    c(7,1) = (1);   
    
elseif ( order == 7 )         % Butcher's Book
    a = zeros(s,s);           % Preallocate
    a(2,1) = (1/6);
        
    a(3,2) = (1/3);
        
    a(4,1) = (1/8);
    a(4,3) = (3/8);
        
    a(5,1) = (148/1331);
    a(5,3) = (150/1331);
    a(5,4) = (-56/1331);
        
    a(6,1) = (-404/243);
    a(6,3) = (-170/27);
    a(6,4) = (4024/1701);
    a(6,5) = (10648/1701);
        
    a(7,1) = (2466/2401);
    a(7,3) = (1242/343);
    a(7,4) = (-19176/16807);
    a(7,5) = (-51909/16807);
    a(7,6) = (1053/2401);
        
    a(8,1) = (5/154);
    a(8,4) = (96/539);
    a(8,5) = (-1815/20384);
    a(8,6) = (-405/2464);
    a(8,7) = (49/1144);
        
    a(9,1) = (-113/32);
    a(9,3) = (-195/22);
    a(9,4) = (32/7);
    a(9,5) = (29403/3584);
    a(9,6) = (-729/512);
    a(9,7) = (1029/1408);
    a(9,8) = (21/16);
        
    b = zeros(1,s);           % Preallocate
    b(1,4) = (32/105);
    b(1,5) = (1771561/6289920);
    b(1,6) = (243/2560);
    b(1,7) = (16807/74880);
    b(1,8) = (77/1440);
    b(1,9) = (11/270);
        
    c = zeros(s,1);           % Preallocate
    c(2,1) = (1/6);
    c(3,1) = (1/3);
    c(4,1) = (1/2);
    c(5,1) = (2/11);
    c(6,1) = (2/3);
    c(7,1) = (6/7);
    c(9,1) = (1);

elseif ( order == 8 )         % Cooper-Verner Method
    a = zeros(s,s);           % Preallocate
    a(2,1) = (1/6);
        
    a(3,1) = (1/4);
    a(3,2) = (1/4);
        
    a(4,1) = (1/7);
    a(4,2) = ((-7-3*sqrt(21))/98);
    a(4,3) = ((21+5*sqrt(21))/49);
        
    a(5,1) = ((11+sqrt(21))/84);
    a(5,3) = ((18+4*sqrt(21))/63);
    a(5,4) = ((21-sqrt(21))/252);
        
    a(6,1) = ((5+sqrt(21))/48);
    a(6,3) = ((9+sqrt(21))/36);
    a(6,4) = ((-231+14*sqrt(21))/360);
    a(6,5) = ((63-7*sqrt(21))/80);
        
    a(7,1) = ((10-sqrt(21))/42);
    a(7,3) = ((-432+92*sqrt(21))/315);
    a(7,4) = ((633-145*sqrt(21))/90);
    a(7,5) = ((-504+115*sqrt(21))/70);
    a(7,6) = ((63-13*sqrt(21))/35);
        
    a(8,1) = (1/14);
    a(8,5) = ((14-3*sqrt(21))/126);
    a(8,6) = ((13-3*sqrt(21))/63);
    a(8,7) = (1/9);
        
    a(9,1) = (1/32);
    a(9,5) = ((91-21*sqrt(21))/576);
    a(9,6) = (11/72);
    a(9,7) = ((-385-75*sqrt(21))/1152);
    a(9,8) = ((63+13*sqrt(21))/128);
        
    a(10,1) = (1/14);
    a(10,5) = (1/9);
    a(10,6) = ((-733-147*sqrt(21))/2205);
    a(10,7) = ((515+111*sqrt(21))/504);
    a(10,8) = ((-51-11*sqrt(21))/56);
    a(10,9) = ((132+28*sqrt(21))/245);
        
    a(11,5) = ((-42+7*sqrt(21))/18);
    a(11,6) = ((-18+28*sqrt(21))/45);
    a(11,7) = ((-273-53*sqrt(21))/72);
    a(11,8) = ((301+53*sqrt(21))/72);
    a(11,9) = ((28-28*sqrt(21))/45);
    a(11,10)= ((49-7*sqrt(21))/18);
        
        
    b = zeros(1,s);           % Preallocate
    b(1,1) = (1/20);
    b(1,8) = (49/180);
    b(1,9) = (16/45);
    b(1,10)= (49/180);
    b(1,11)= (1/20);
        
    c = zeros(s,1);           % Preallocate
    c(2,1) = (1/2);
    c(3,1) = (1/2);
    c(4,1) = ((7+sqrt(21))/14);
    c(5,1) = ((7+sqrt(21))/14);
    c(6,1) = (1/2);
    c(7,1) = ((7-sqrt(21))/14);
    c(8,1) = ((7-sqrt(21))/14);
    c(9,1) = (1/2);
    c(10,1)= ((7+sqrt(21))/14);
    c(11,1)= (1);
end

end


%##########################################################################
%                            u for try 2
%##########################################################################

function out = create_u_n(n, shift)

out = [(sin(2*pi*linspace(0,1,floor(n/4))) + shift), shift*ones(1,floor(3*n/4))]';

end



