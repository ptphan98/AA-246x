clc;
close all;
clear all;
% ALL UNITS IN SI %
% Input Parameters to optimize:
% W_L - Wing loading
% weight - aicraft weight (kg)
% b - wing span (m)
% V_cruise - Cruise velocity (m/s)

global rho 
rho = 1.225; %kg/m^3

%attempt an optimization with constraints
x_0 = [5, 2, 2, 30];
A = [-1 0 0 0;
     1 0 0 0;
     0 -1 0 0;
     0 0 -1 0;
     0 0 0  1];
b = [0; 9; 0; 0; 50]; 

x = fmincon( @(x) optimize( x(1),x(2),x(3),x(4)), x_0, A,b );

W_L_op = x(1)
weight_op = x(2)
b_op = x(3)
v_cruise_op = x(4)
AR = b_op^2/(weight_op/W_L_op)

[mass_total, mass_bat, mass_motor, mass_empty, L_D] = calc_weight(W_L_op, weight_op, b_op, v_cruise_op);
[S_ref, c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(W_L_op, weight_op, b_op, v_cruise_op);
battery_frac = mass_bat/mass_total
motor_frac = mass_motor/mass_total
structural_frac = mass_empty/mass_total
plot_aircraft(W_L_op, weight_op, b_op, v_cruise_op)

%% plot aircraft 
function [] = plot_aircraft(W_L, weight, b, v_cruise)
    boom_dia = 0.0154; %m - visual purpose only
    [S_ref, c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(W_L, weight, b, v_cruise);
    
    %estimate fuselage area
    l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
    r_fus = l_fus./8./2; %assume fuselage fineness ratio of 8

    %we model the fuselage as cone nose and cylinder fuselage
    %assume nose is 1/5 length of entire fuselage
    l_nose = 1/5.*l_fus;
    tail_end = l_t+.75*c_htail;
    
    pgon = polyshape([tail_end tail_end tail_end-l_fus  tail_end-l_fus],[-r_fus r_fus r_fus -r_fus]);
    
    figure;
    hold on
    rectangle('Position',[.24*c 0-boom_dia/2 l_t boom_dia],'EdgeColor','k','LineWidth',1,'FaceColor',[.3 .3 .3]); %boom
    rectangle('Position',[(c/4 + l_t-.25*c_htail) (0-(s_htail/c_htail)/2) c_htail s_htail/c_htail],'EdgeColor','r','LineWidth',3); %horizontal tail
    rectangle('Position',[ (c/4 + l_t-.25*c_vtail) 0 c_vtail s_vtail/c_vtail],'EdgeColor','b','LineWidth',3); %vertical tail
    rectangle('Position',[0 0-b/2 c b],'EdgeColor','k','LineWidth',3); %wing
    plot(pgon)
    axis equal
end

%%
function score = optimize(W_L, weight, b, v_air)    
    %define useful parameters
    S_ref = weight./W_L;
    c = S_ref./b;
    AR = b.^2./S_ref;
    
    [mass_total, mass_bat, mass_motor, mass_empty, L_D] = calc_weight(W_L, weight, b, v_air);
    
    %constraint - add constrains as necessary
    cost1 = abs(weight-mass_total); %make sure mass converges
    cost2 = 10000;
    if AR > 3.5
        cost2  = 0; %force an aspect ratio
    end
    cost3 = 1000;
    if L_D < 15
        cost3 = 0; %set max limit for L/D
    end
   
    
    score = -v_air + cost1.*10 + cost2 + cost3; %maximize velocity given constraints
end


%%
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
l_t = 0.75.*l_fus; %ballpark for length of 1/4 chord to 1/4 tail chord. Should be optimized for weight!

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


