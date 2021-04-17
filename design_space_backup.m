clc;
clear all;
% ALL UNITS IN SI %

% Input Parameters to optimize:
% W_L - Wing loading
% weight - aicraft weight (kg)
% b - wing span (m)
% V_cruise - Cruise velocity (m/s)

global rho 
rho = 1.225; %kg/m^3

%design space
W_L = 3:.25:6; %kg/m^2
weight = 0.25:.05:2; %kg
b = .5:.25:5; %m
v_cruise = 20:5:50; %m/s

[W_L,weight,b,v_cruise] = ndgrid(W_L,weight,b,v_cruise);

[S_ref, c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(W_L, weight, b, v_cruise);
[mass_total, mass_bat, mass_motor, mass_empty, L_D] = calc_weight(W_L, weight, b, v_cruise);

%check if mass converges
mass_check = abs(weight-mass_total);
valid = find(mass_check < 0.01); %index of elements that are valid

%throw out non-valid values
W_L2 = W_L(valid);
weight2 = weight(valid);
b2 = b(valid);
v_cruise2 = v_cruise(valid);

S_ref2 = S_ref(valid);
c2 = c(valid);

mass_total2 = mass_total(valid);
mass_bat2 = mass_bat(valid);
mass_motor2 = mass_motor(valid);
mass_empty2 = mass_empty(valid);
L_D2 = L_D(valid);


[V_cruise_op, index] = max(v_cruise2);

%display optimum

W_L_op = W_L2(index);
weight_op = weight2(index)
b_op = b2(index);
v_cruise_op = v_cruise2(index);

S_ref_op = S_ref2(index);
c_op = c2(index);
AR = b_op^2/S_ref_op;

mass_total_op = mass_total2(index)
mass_bat_op = mass_bat2(index)
mass_motor_op = mass_motor2(index)
mass_empty_op = mass_empty2(index)
L_D_op = L_D2(index)

[mass_total, mass_bat, mass_motor, mass_empty, L_D] = calc_weight(W_L_op, weight_op, b_op, v_cruise_op)


function [S_ref, c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(W_L, weight, b, v_air)
S_ref = weight./W_L;
%gives conventional tail dimensions based on wing geometry

%guess some typical values for tail parameters
static_margin = 0.05;
ht_vol_cf = 0.40;
vt_vol_cf = 0.03;
AR_ht = 4;
AR_vt = 1.5; 

c = S_ref./b; %mean chord estimate
l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
l_t = 0.6.*l_fus; %ballpark for length of 1/4 chord to 1/4 tail chord. Should be optimized for weight!

s_htail = ht_vol_cf.*S_ref.*c./l_t;
c_htail = sqrt(s_htail./AR_ht);

s_vtail = vt_vol_cf.*S_ref.*b./l_t;
c_vtail = sqrt(s_vtail./AR_vt);

end 

function [Cd,Cdi,Cd0] = calc_Cd(W_L, weight, b, v_air)
global rho
[S_ref, c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(W_L, weight, b, v_air);

AR = b.^2./S_ref;
CL = 9.81.*W_L./(0.5.*rho.*v_air.^2);
e = .9; %this obviously varies with speed
Cdi = CL.^2./(pi.*AR.*e); %lift induced drag

%estimate fuselage area
l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
r_fus = l_fus./8./2; %assume fuselage fineness ratio of 8

%we model the fuselage as cone nose and cylinder fuselage
%assume nose is 1/5 length of entire fuselage
l_nose = 1/5.*l_fus;

Cdf_nose = 2./sqrt(3).*calc_Cf(l_nose,v_air).*(pi.*r_fus.*sqrt(r_fus.^2+l_nose.^2))./S_ref;
Cdf_fuse = calc_Cf(l_fus,v_air).*(2.*pi.*r_fus.*l_fus)./S_ref;
Cdf_wing = calc_Cf(c,v_air).*(2.*S_ref)./S_ref;
Cdf_tail = calc_Cf(c_htail,v_air).*(2*s_htail+2*s_vtail)./S_ref; 

Cd0 = 1.25.*(Cdf_nose + Cdf_fuse + Cdf_wing + Cdf_tail);

Cd = Cd0 + Cdi;

    function CF = calc_Cf(l,v_air)
        %flat plate assumption for Cf
        nu = 15.52e-6; %m^2/s
        Re = v_air.*l./nu;
        
        if Re > 5e5
            CF = 0.455./(log10(Re).^2.58);
        else
            CF = 1.328./sqrt(Re);
        end
    end
end

function [mass_bat, L_D] = battery_weight(W_L, weight, b, v_air)
    S_ref = weight./W_L;
    [Cd,Cdi,Cd0] = calc_Cd(W_L, weight, b, v_air);
    AR = b.^2./S_ref;
    e = .99;
    K = 1./(pi.*AR.*e);
    L_D = 1./(2.*sqrt(Cd0.*K));
   
    eta_sys = .374; %battery + propellor + motor effiency estimate
    range = 42195/2; %meters
    H_batt = 155*3600; %specific energy density (Joules/kg)
    g = 9.807; %m/s^2
  
    mass_bat = range.*g.*weight./(L_D.*eta_sys.*H_batt);
    
end

function [mass_motor] = motor_weight(W_L, weight, b, v_air)
    global rho
    S_ref = weight./W_L;
    
    %calculates thrust required for cruise at max L/D
    [Cd,Cdi,Cd0] = calc_Cd(W_L, weight, b, v_air);
    Power_cruise = 1/2.*rho.*v_air.^3.*S_ref.*(2.*Cd0); %watts - At max L/D, Cd0 = Cdi
    Power_max = Power_cruise./0.40; %assume motor runs most efficient at 40% max power
    
    motor_p_density = 5.101/0.001; %watts/kg of a motor 
    
    mass_motor =  1./motor_p_density.*Power_max;
end

function [mass_empty] = empty_weight(W_L, weight, b, v_air)
    S_ref = weight./W_L;
    [s_htail, c_htail, s_vtail, c_vtail] = size_plane(W_L, weight, b, v_air);

    %estimate fuselage area
    l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
    r_fus = l_fus./8./2; %assume fuselage fineness ratio of 8
    
    %weight estimations based off walrus
    m_wing = S_ref.*1.31346914; %kg/m^2 x s_wing
    m_htail = s_htail.*1.695180676; %kg/m^2 x s_htail
    m_vtail = s_vtail.*0.792369773; %kg/m^2 x s_vtail 
    m_fuse = (pi.*r_fus.^2).*l_fus.*29.43741259; %kg/m^3 x vol_fuselage 
    
    mass_empty = m_wing + m_htail + m_vtail + m_fuse;
end

function [mass_total, mass_bat, mass_motor, mass_empty, L_D] = calc_weight(W_L, weight, b, v_air)
    [mass_bat, L_D] = battery_weight(W_L, weight, b, v_air);
    [mass_motor] = motor_weight(W_L, weight, b, v_air);
    [mass_empty] = empty_weight(W_L, weight, b, v_air);
    
    mass_total = mass_bat + mass_motor + mass_empty;
end


