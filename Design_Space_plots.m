clc;
close all;
clear all;
% ALL UNITS IN SI %
% Input Parameters to optimize:
% W_L - Wing loading (Kg/m^2)
% weight - aicraft weight (kg)
% b - wing span (m)
% V_cruise - Cruise velocity (m/s)

global rho 
rho = 1.225; %kg/m^3

%attempt an optimization with constraints
W_L_lo = 0; %Kg/m^2
Weight_lo = 0; %Kg
Span_lo = 0.1:0.05:1; %m
V_cruise_lo = 0; %m/s

W_L_up = 10; %Kg/m^2
Weight_up = 10; %Kg
Span_up = +Inf; %m
V_cruise_up = +Inf; %m/s

W_L_op = zeros(1,length(Span_lo));
weight_op = zeros(1,length(Span_lo));
b_op = zeros(1,length(Span_lo));
v_cruise_op = zeros(1,length(Span_lo));

for i = 1:length(Span_lo)

    x_0 = [5, 2, 2, 30];
    A = [];
    b = [];
    Aeq = [0 0 1 0];
    beq = [Span_lo(i)];
    lb = [W_L_lo,Weight_lo,0,V_cruise_lo];
    ub = [W_L_up,Weight_up,Span_up,V_cruise_up];

    x = fmincon( @(x) optimize( x(1),x(2),x(3),x(4)), x_0,A,b,Aeq,beq,lb,ub, @(x) nonl_const( x(1),x(2),x(3),x(4) ) );

    W_L_op(i) = x(1);
    weight_op(i) = x(2);
    b_op(i) = x(3);
    v_cruise_op(i) = x(4);
    
    i
end

save('design_space_plot_data')
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
    rectangle('Position',[(c/4 + l_t-.25*c_htail) (0-(s_htail/c_htail)/2) c_htail s_htail/c_htail],'EdgeColor','r','LineWidth',3); %horizontal tail
    rectangle('Position',[ (c/4 + l_t-.25*c_vtail) 0 c_vtail s_vtail/c_vtail],'EdgeColor','b','LineWidth',3); %vertical tail
    rectangle('Position',[0 0-b/2 c b],'EdgeColor','k','LineWidth',3); %wing
    plot(pgon)
    axis equal
end

%% optimizing functions

function score = optimize(W_L, weight, b, v_air)    
    %define useful parameters
    S_ref = weight./W_L;
    c = S_ref./b;
    AR = b.^2./S_ref;
    
    [mass_total, mass_bat, mass_motor, mass_empty] = calc_weight(W_L, weight, b, v_air);
    [CL, Cd, Cdi , Cd0, L_D, v_ideal] = calc_aero(W_L, weight, b, v_air);

    score = -v_air; %maximize velocity given constraints
end

% nonlinear constraints
function [c, ceq] = nonl_const(W_L, weight, b, v_air)
    [mass_total, mass_bat, mass_motor, mass_empty] = calc_weight(W_L, weight, b, v_air);
    [CL, Cd, Cdi , Cd0, L_D, v_ideal] = calc_aero(W_L, weight, b, v_air);
    [time,x_p,y_p,y_min] = calc_takeoff(W_L, weight, b, v_air);
    [bend_stress,deflection_span] = calc_beam(W_L, weight, b, v_air);
    
    S_ref = weight./W_L;
    AR = b.^2./S_ref;
    
    %constrain aspect ratio
    AR_lower = 3.5;
    AR_upper = 15;

    max_stress = 3.5*10^9; %failure stress of carbon fiber - assume 3.5 GPa;
    
    c = [- AR + AR_lower, AR - AR_upper, -y_min, bend_stress - max_stress];
    ceq = abs(weight-mass_total);
end

%% aircraft sizing functions

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

function [CL, Cd, Cdi , Cd0, L_D, v_ideal] = calc_aero(W_L, weight, b, v_air)
global rho
[S_ref, c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(W_L, weight, b, v_air);

AR = b.^2./S_ref;
CL = 9.81.*W_L./(0.5.*rho.*v_air.^2);
e = .9; %this obviously varies with speed
K = 1./(pi.*AR.*e);
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

L_D = CL/Cd;

%Find the optimum cruise speed, v_ideal, that maximizes L/D for this
%configuration
Cl_opt = sqrt(Cd0/K);
v_ideal = sqrt( (2*9.81*weight)/(rho*Cl_opt*S_ref) );


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

function [mass_bat] = battery_weight(W_L, weight, b, v_air)
    [CL, Cd, Cdi , Cd0, L_D, v_ideal] = calc_aero(W_L, weight, b, v_air);
   
    eta_sys = .374; %battery + propellor + motor effiency estimate
    range = 42195/2; %meters
    H_batt = 155*3600; %specific energy density (Joules/kg)
    g = 9.807; %m/s^2
  
    mass_bat = range.*g.*weight./(L_D.*eta_sys.*H_batt);
end

function [mass_motor,T_W,Thrust,Power_max] = motor_weight(W_L, weight, b, v_air)
    global rho
    S_ref = weight./W_L;
    
    %calculates thrust required for cruise at v_air
    [CL, Cd, Cdi , Cd0, L_D, v_ideal] = calc_aero(W_L, weight, b, v_air);
    Thrust = 1/2.*rho.*v_air.^2.*S_ref.*(Cd); %Newtons
    Power_cruise = 1/2.*rho.*v_air.^3.*S_ref.*(Cd); %watts
    Power_max = Power_cruise./0.40; %assume motor runs most efficient at 40% max power
    
    motor_p_density = 5.101/0.001; %watts/kg of a motor 
    T_W = Thrust*2/(weight*9.81); %static thrust to weight. assume static is twice cruise thrust 
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
    m_fuse = (pi.*r_fus.^2).*l_fus.*29.43741259; %kg/m^2 x s_fuselage 
    
    mass_empty = m_wing + m_htail + m_vtail + m_fuse;
end

function [mass_total, mass_bat, mass_motor, mass_empty, T_W] = calc_weight(W_L, weight, b, v_air)
    [mass_bat] = battery_weight(W_L, weight, b, v_air);
    [mass_motor,T_W] = motor_weight(W_L, weight, b, v_air);
    [mass_empty] = empty_weight(W_L, weight, b, v_air);
    
    mass_total = mass_bat + mass_motor + mass_empty;
end

function [bend_stress,deflection_span] = calc_beam(W_L, weight, b, v_air)
%define important airfoil properties
S_ref = weight./W_L;
c = S_ref./b;
t_c = 0.08; %thickness to chord ratio of airfoil
g_load = 1.5; %how much load do we expect?

%assume elliptical lift distribution
l_0 = (4*weight*9.81*g_load)/(pi*b); %make sure load in newtons
lift = @(x) l_0.*sqrt(1-(x./(b/2)).^2);

x = 0:0.001:b/2;
y = lift(x);
y_flip = flip(y); %cumtrapz integrates in wrong direction

shear = cumtrapz(x,y_flip); %shear is integral load distribution
moment = cumtrapz(x,shear); %moment is integral shear

%orient shear and moment in right directions
shear = flip(shear);
moment = flip(moment);

%plots - helpful for double checking
% figure; plot(x,y); title('Lift Distribution'); 
% xlabel('Span (m)'); ylabel('Load per unit span (N/m)');
% figure; plot(x,shear); title('Shear');
% xlabel('Span (m)'); ylabel('Shear (N)');
% figure; plot(x,moment); title('Moment')
% xlabel('Span (m)'); ylabel('Moment (N*m)');

%find max_stress
r_spar = 0.75*t_c*c/2; %assume spar diameter is 0.75 thickness of airfoil 
I_c = pi/4*r_spar^4;   %2nd moment of inertia
bend_stress = max(moment)*r_spar/I_c;

%find deflection
E = 228*10^9; %modulus of carbon fiber - 228 GPa

d2v_dx2 = moment/(E*I_c);
rotation = cumtrapz(x,d2v_dx2);
deflection = cumtrapz(x,rotation);

deflection_span = deflection(end)/b; %deflection to span ratio
end

function [time,x_p,y_p,y_min] = calc_takeoff(W_L, weight, b, v_air)
    global rho
    g = 9.81; %gravity
    S_ref = weight./W_L;
    V_throw = 7; %m/s - Velocity we can throw airplane
    height_throw = 1.8; %m - height aircraft is thrown from 
    Cl_max = 1.2; %assume a Cl_max at takeoff
    
    [CL, Cd, Cdi , Cd0, L_D, v_ideal] = calc_aero(W_L, weight, b, v_air);
    [mass_motor,T_W,Thrust,Power_max] = motor_weight(W_L, weight, b, v_air);
    
    %solve takeoff trajectory
    tspan = [0, 5];
    x0 = [0, V_throw, height_throw, 0];
    [time,x] = ode45(@(t,x) odefun(t,x), tspan, x0);
    
    x_p = x(:,1);
    y_p = x(:,3);
    y_min = min(y_p); %make sure the aircraft doesn't hit the ground

    function dxdt = odefun(t, x)
        dxdt = zeros(4,1);
        %X direction of takeoff
        dxdt(1) = x(2); %velocity
        dxdt(2) = Power_max/(weight*x(2)) - 1/(2*weight)*rho*x(2)^2*Cd*S_ref; %acceleration 

        %Y direction of takeoff
        dxdt(3) = x(4); %velocity
        dxdt(4) = 1/(2*weight)*rho*x(2)^2*Cl_max*S_ref - g; %acceleration 
    end
end


