%function steady_state()

%
% Defining the parameters
%

nx = 1000;
dt = 1;
a = 0;
b = 20;
dx = (b-a)/(nx);
x = linspace(a, b-dx, nx);
diff = 1/(6*dx^2);
perf = (2/3)*(10^(-3));

%
% Initial condition
%

u = zeros(nx,1);

%
%  Creating the matrix for the diffusion term
%

main_second_der = -2 * ones(nx,1);
off_secon_der = ones(nx-1,1);
Diff_matrix1 = diag(main_second_der) + diag(off_secon_der, -1) + diag(off_secon_der,1);



main_first_der = zeros(nx,1);
off_up_first_der = ones(nx-1,1);
off_down_first_der = ones(nx-1,1);
for i = 2:(nx-2)
    off_up_first_der(i+1) = off_up_first_der(i+1)/i;
end
for i = 1:(nx-1)
    off_down_first_der(i) = - off_down_first_der(i)/i;
end
Diff_matrix2= diag(main_first_der) + diag(off_down_first_der, -1) + diag(off_up_first_der,1);

Diff3 = Diff_matrix1 + Diff_matrix2;

%
%  Adding in the Neumann condition at the origin
%

Diff3(1,1) = -6;
Diff3(1, 2) = 6;

Diff3 = Diff3* diff;
Diff3 = sparse(Diff3);

%
% Adding in the perfusion term
%

I = speye(nx);
Perf_matrix = perf * speye(nx);

%
% Adding in the source term Q
%

Source = sin(x)./x;
Source = besselj(1,x)./x;
Source(1) = 1;
 news =  5* Source.^2;

news = transpose(news);


%
% Timestepping using the Crank-Nicolson method
%

figure
hold on

trapint = u;

for i = 1:60
    uold = u;
    tt(i)=i*dt;
    AA = (I - (dt/2)*Diff3 + (dt/2)*Perf_matrix);
    
      if tt(i) < 30 
      new_step = (I + (dt/2)*Diff3 - (dt/2)*Perf_matrix)*u + dt*news;
      else
      new_step = (I + (dt/2)*Diff3 - (dt/2)*Perf_matrix)*u; 
      end
     
    u = AA\new_step;
    uu=u;
    xx=x;
    xx(nx+1)=b;
    uu(nx+1)=0;
    %Computing the integral for the damage using trapezoidal rule
    A = 226.8;
    E = 6.28*(10^5);
    R = 8.314;
    %Function to evaluate
    fold = exp(226.8 - E./(R*(uold+310)));
    f    = exp(226.8 - E./(R*(u+310)));

    trapint = trapint+0.5*(fold+f)*dt
    hold on
    
    if tt(i)< 30
    %subplot(1,2,1)
    %plot(xx,uu,'r')
    %xlabel('x')
    %ylabel('Temperature profile')
    %leg1= legend('heating up')
    %else
    %subplot(1,2,2)
    %plot(xx,uu,'b')
    %xlabel('x')
    %ylabel('Temperature profile')
    %leg2= legend('cooling down')
    end
    ux(i)=uu(1);
    %subplot(1,3,2)
    %plot(x, trapint)
    %xlabel('x')
    %ylabel('Damage')
    
end
     hold off 
     %subplot(1,3,3)
     plot(tt,ux, 'r', 'linewidth', 2)
     xlabel('time')
     ylabel('Temperature profile')   