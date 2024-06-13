R_s = 1.633 ; % Ohms
L_s = 0.142 ; % H
R_r = 0.930 ; % Ohms
L_r = 0.076 ; % H
p = 2 ;       % poles   
L_m = 0.099 ; % H
J = 0.029 ;   % moment of inertia
fv = 0.0038 ; % friction

ics = [0 0 0 0 0 ] ;

parameters = [R_s L_s R_r L_r p L_m J fv ]

a = (R_r * L_m)/L_r ;
b = R_r / L_r ;

sigma = 1 - (L_m*L_m)/(L_s*L_r) ;

gamma_1 = (R_s)/(sigma*L_s) + (R_r*L_m*L_m)/(sigma*L_s*L_r*L_r) ;
gamma_2 = ( R_r * L_m )/(sigma*L_s*L_r*L_r) ;
gamma_3 = (p*L_m)/(sigma*L_s*L_r) ;
gamma_4 = 1/(sigma*L_s) ;

alpha_1 = (p*L_m)/(J*L_r) ;
alpha_2 = 1/J ;
alpha_3 = fv/J ;

p = 2 ;
w0 = 152 ; % rad/s

AMPLITUDE = 100 ;
FREQUENCY = 2*pi*50 ;
