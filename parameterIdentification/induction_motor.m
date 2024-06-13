function  [sys, x0] = induction_motor(t, x, u, flag, params)
% equations
% \Gamma_{m}=K_{t} i_{a}
% J_{m}\ddot{\theta}_{m}+b\dot{\theta}_{m}=\Gamma_{m}-\Gamma_{c}\\
% L_{a}\frac{di_{a}}{dt}+R_{a}i_{a}=v_{a}-K_{e}\dot{\theta}_{m}

% States
%x1 = \omega  \\
%x2 = \phi_{r \alpha} \\
%x3= \phi_{r \beta} \\
%x4 = i_{s \alpha} \\
%x4 = i_{s \beta} \\

% State equations
% dx1 = \dot{\omega} = x(2);
% dx2 = (Kt*ia-\Gamma_{c}-b*\dot{\theta}_{m})/Jtot=(ke*x(3)-u(2)-b*x(2))/Jtot;
% dx3 = (v_{a}-K_{e}\dot{\theta}_{m}-R_{a}i_{a})/La=(u(1)-ke*x(2)-R*x(3))/L;
% params_real=[R;L;ke;b;Jtot;umax];

R_s = params(1) ;
L_s = params(2) ;
R_r = params(3) ;
L_r = params(4) ;
p = params(5) ; 
L_m = params(6) ;
J = params(7) ;
fv = params(8) ;

alpha_1 = (p*L_m)/(J*L_r) ;
alpha_2 = 1/J ;
alpha_3 = fv/J ;

a = (R_r * L_m)/L_r ;
b = R_r / L_r ;

sigma = 1 - (L_m*L_m)/(L_s*L_r) ;

gama_1 = (R_s)/(sigma*L_s) + (R_r*L_m*L_m)/(sigma*L_s*L_r*L_r) ;
gama_2 = ( R_r * L_m )/(sigma*L_s*L_r*L_r) ;
gama_3 = (p*L_m)/(sigma*L_s*L_r) ;
gama_4 = 1/(sigma*L_s) ;

Cl = 0 ;

%% Initial condition
% IC.thetar;IC.omegar;IC.ir;
%states_ic=[IC.omega ; IC.phi_ra ; IC.phi_rb ; IC.i_sa ; IC.i_sb ];
states_ic=[ 0 0 0 0 0 ];

if flag == 1 % If flag = 1, return state derivatives, xDot
    sys(1,1) = alpha_1*(x(2)*x(5) - x(3)*x(4)) - alpha_2*u(1) + alpha_3*x(1) ;
    sys(1,2) = a * x(4) - b * x(2) - p*x(1)*x(3) ; 
    sys(1,3) = a * x(5) - b * x(3) + p*x(1)*x(2) ;
    sys(1,4) = -gama_1*x(4) + gama_2*x(2)+gama_3*x(1)*x(3)+gama_4*u(2) ;
    sys(1,5) = -gama_1*x(5) + gama_2*x(3)-gama_3*x(1)*x(2)+gama_4*u(3) ;
    
elseif flag == 3 % If flag = 3, return system outputs, y
    sys(1,1) = x(1); 
	sys(2,1) = x(2);
    sys(3,1) = x(3); 
	sys(4,1) = x(4);
    sys(5,1) = x(5);

elseif flag == 0 % If flag = 0, return initial condition data, sizes and x0
   %sys = [nb d'états continus; nb d'états discrets; nb de sorties; nb d'entrées;0; 0]; 
    sys = [5; 0; 5; 3; 0; 0];  
    x0 = [states_ic(1); states_ic(2); states_ic(3); states_ic(4); states_ic(5) ];

else 
    % If flag is anything else, no need to return anything
    % since this is a continuous system
    sys = [];
end
