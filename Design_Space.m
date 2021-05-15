% Bulk of this code was written by Philip with input from Jean and Karthik 
% Three functions (empty_weight, calc_beam, and openfile) were written by Yamaan
% with modifications by Philip to get integrated into the code.

clc;
close all;
clear all;
% ALL UNITS IN SI %
% Input Parameters to optimize:
% weight - aicraft weight (kg)
% b - wing span (m)
% V_cruise - Cruise velocity (m/s)

global cst 
cst = struct();
cst.rho = 1.225; %kg/m^3
cst.g = 9.807; %m/s^2
cst.CL_max = 1.2;
cst.V_stall = 7; %m/s
cst.W_L = 1/2 * cst.rho * cst.V_stall^2 * cst.CL_max / cst.g; % wing loading is sized by stall speed
cst.spar_ratio = .7; %percent spar of max airfoil thickness
cst.l_fus = 0.46; %from the CAD
cst.r_fus = 0.040; %fits our lipo comfortably
[~,cst.xin,cst.yin]=openfile('ag0308.dat'); %WING AIRFOIL, EDIT
[~,cst.xtin,cst.ytin]=openfile('naca0008.dat'); %TAIL AIRFOIL, EDIT

% Lower and upper bounds
Weight_lo = 0; %Kg
Span_lo = 0.1; %m
V_cruise_lo = 0; %m/s
l_t_lo = 0; %NOTE: l_t is quarter chord of wing to quarter chord tail

Weight_up = 2.04117; %Kg
Span_up = 1.5; %m
V_cruise_up = +Inf; %m/s
l_t_up = +Inf; %m

x_0 = [1,1.5,30,.75];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [Weight_lo,Span_lo,V_cruise_lo,l_t_lo];
ub = [Weight_up,Span_up,V_cruise_up,l_t_up];

options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluation', 10000);
%options = optimoptions('fmincon','MaxFunctionEvaluation', 10000);
x = fmincon( @(x) optimize( x(1),x(2),x(3),x(4)), x_0,A,b,Aeq,beq,lb,ub, @(x) nonl_const( x(1),x(2),x(3),x(4)),options);

%show optimized aircraft 
weight = x(1);
b = x(2);
v_cruise = x(3);
l_t = x(4);
AR = b^2/(weight/cst.W_L);
S_ref = weight / cst.W_L;
chord = S_ref/b;
time = 21100/v_cruise/60

[CL, Cd, Cdi , Cd0, L_D, v_ideal, Drag] = calc_aero(weight, b, v_cruise, l_t);
[sigma_max, deflection_span] = calc_beam(S_ref, weight, b);
[mass_empty, m_fixed, m_wing, m_htail, m_vtail, m_boom] = empty_weight(S_ref, b, l_t);
[mass_bat] = battery_weight(weight, L_D);
[mass_motor,T_W,Power_max] = motor_weight(weight, Drag, v_cruise);
mass_total = mass_bat + mass_motor + mass_empty;
[c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b, l_t);
CG = calc_CG(S_ref, b, l_t, m_fixed, m_wing, m_htail, m_vtail, m_boom, mass_motor, mass_bat);
b_htail = s_htail/c_htail;
b_vtail = s_vtail/c_vtail;

%output useful info

b
chord
S_ref
b_htail
c_htail
s_htail
b_vtail
c_vtail
s_vtail
l_t
l_boom = l_t - (c*.75) + .18*cst.l_fus

battery_frac = mass_bat/mass_total
motor_frac = mass_motor/mass_total
emptyweight_frac = mass_empty/mass_total

plot_aircraft(weight, S_ref, b, v_cruise, l_t)

%% plot aircraft - written by Philip
function [] = plot_aircraft(weight, S_ref, b, v_air, l_t)
    global cst
    [c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b, l_t);
    
    %zero is defined as quarter chord
    %estimate fuselage area
    l_fus = cst.l_fus; 
    r_fus = cst.r_fus; 
    r_boom = 0.0127; %half an inch
    l_boom = l_t - (c*.75) + .18*l_fus;
    
    %we model the fuselage as cone nose and cylinder fuselage
    tail_quart = l_t;
    
    %plots boom
    pgon = polyshape([tail_quart tail_quart tail_quart-l_boom  tail_quart-l_boom],[-r_boom r_boom r_boom -r_boom]);
    
    %plots fuselage
    pgon2 = polyshape([tail_quart-l_t+c*.75 tail_quart-l_t+c*.75 tail_quart-l_t-l_fus+c*.75  tail_quart-l_t-l_fus+c*.75],[-r_fus r_fus r_fus -r_fus]);
    
    figure;
    hold on
    rectangle('Position',[(tail_quart-.25*c_htail) (0-(s_htail/c_htail)/2) c_htail s_htail/c_htail],'EdgeColor','r','LineWidth',3); %horizontal tail
    rectangle('Position',[ (tail_quart-.25*c_vtail) 0 c_vtail s_vtail/c_vtail],'EdgeColor','b','LineWidth',3); %vertical tail
    rectangle('Position',[-c/4 0-b/2 c b],'EdgeColor','k','LineWidth',3); %wing
    pg = plot(pgon);
    pg.FaceColor = 'k';
    pg2 = plot(pgon2);
    
    %plot CG
    [mass_empty, m_fixed, m_wing, m_htail, m_vtail, m_boom] = empty_weight(S_ref, b, l_t);
    [CL, Cd, Cdi , Cd0, L_D, v_ideal, Drag] = calc_aero(weight, b, v_air, l_t);
    [mass_bat] = battery_weight(weight, L_D);
    [mass_motor,T_W,Power_max] = motor_weight(weight, Drag, v_air);
    [CG,x_check] = calc_CG(S_ref, b, l_t, m_fixed, m_wing, m_htail, m_vtail, m_boom, mass_motor, mass_bat);
    
    plot(CG,0,'*k')
    %plot(x_check,zeros(1,length(x_check)),'*k')
    axis equal
end

%% optimizing functions - written by Philip

function score = optimize(weight, b, v_air, l_t)    
    score = -v_air; %maximize velocity and minimize weight
end

% nonlinear constraints - written by Philip
function [c, ceq] = nonl_const(weight, b, v_air,l_t)
    global cst
    S_ref = weight / cst.W_L; 
    c = S_ref/b;
    [mass_empty, m_fixed, m_wing, m_htail, m_vtail, m_boom] = empty_weight(S_ref, b, l_t);
    [CL, Cd, Cdi , Cd0, L_D, v_ideal, Drag] = calc_aero(weight, b, v_air, l_t);
    [sigma_max, deflection_span] = calc_beam(S_ref, weight, b);
    [mass_bat] = battery_weight(weight, L_D);
    [mass_motor,T_W,Power_max] = motor_weight(weight, Drag, v_air);
    CG = calc_CG(S_ref, b, l_t, m_fixed, m_wing, m_htail, m_vtail, m_boom, mass_motor, mass_bat);
    mass_total = mass_bat + mass_motor + mass_empty;
    l_boom = l_t - (c*.75) + .18*cst.l_fus; 
    
    %Wing Structure Constraints
    max_stress = 1.8*10^9; %failure stress of Carbon fiber
    max_deflection_span = .10; %tip should not deflection more than 10% span.
    
    %constrain aspect ratio
    AR = b.^2./S_ref;
    AR_lower = 0;
    AR_upper = 20;

    %c = [- AR + AR_lower, AR - AR_upper, sigma_max - max_stress,deflection_span-max_deflection_span];
    c = [sigma_max - max_stress,deflection_span-max_deflection_span,mass_total - weight, mass_bat-0.218, l_boom-0.9144];
    ceq = [];
end

%% aircraft sizing functions

% Written by Philip
function [c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b, l_t)
%gives conventional tail dimensions based on wing geometry
%guess some typical values for tail parameters
static_margin = 0.05;
ht_vol_cf = 0.40;
vt_vol_cf = 0.03;
AR_ht = 4;
AR_vt = 1.5; 

c = S_ref./b; %mean chord estimate
l_t = l_t; %ballpark for length of 1/4 chord to 1/4 tail chord. Should be optimized for weight!

s_htail = ht_vol_cf.*S_ref.*c./l_t;
c_htail = sqrt(s_htail./AR_ht);
s_vtail = vt_vol_cf.*S_ref.*b./l_t;
c_vtail = sqrt(s_vtail./AR_vt);
end 

% Written by Philip
function [CL, Cd, Cdi , Cd0, L_D, v_ideal, Drag] = calc_aero(weight, b, v_air,l_t)
global cst
S_ref = weight / cst.W_L;
[c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b, l_t);

AR = b.^2./S_ref;
CL = 9.81.* cst.W_L./(0.5.*cst.rho.*v_air.^2);
e = .9; %this obviously varies with speed
K = 1./(pi.*AR.*e);
Cdi = CL.^2./(pi.*AR.*e); %lift induced drag

%estimate fuselage area
l_fus = cst.l_fus; %measured on CAD
r_fus = cst.r_fus; 

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
v_ideal = sqrt( (2*9.81*weight)/(cst.rho*Cl_opt*S_ref) );
Drag = 1/2 * cst.rho * v_air^2 * S_ref * Cd; %Newtons

    function CF = calc_Cf(l,v_air)
        %flat plate assumption for Cf
        nu = 15.52e-6; %m^2/s
        Re = v_air.*l./nu;
        CF = 0.455./(log10(Re).^2.58); %assume airflow always turbulent
    end
end

% Written by Philip
function [mass_bat] = battery_weight(weight, L_D)
    global cst
    eta_sys = .33; %battery + propellor + motor effiency estimate + ESC
    safety_factor = 1.1;
    range = 42195/2*safety_factor; %meters
    H_batt = 155*3600; %specific energy density (Joules/kg)
    mass_bat = range.*cst.g.*weight./(L_D.*eta_sys.*H_batt);
end

% Written by Philip
function [mass_motor,T_W,Power_max] = motor_weight(weight, Thrust, v_air)
    Power_cruise = Thrust * v_air; %watts
    
    %assume motor power is 2x mechanical power
    %assume motor runs most efficient at 40% max power
    Power_max = Power_cruise./0.40; 
    
    motor_p_density = 5.101/0.001; %watts/kg of a motor 
    T_W = Thrust*2/(weight*9.81); %static thrust to weight. assume static is twice cruise thrust 
    mass_motor =  1./motor_p_density.*Power_max *2; % times 2 for prop efficiency
end

% Written by Yamaan, modified by Philip
function [mass_empty, m_fixed, m_wing, m_htail, m_vtail, m_boom] = empty_weight(S_ref, b, l_t)
    global cst
    [c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b, l_t);
    
    %%%% Material Properties
    rho_cf = 2000; %kg/m^3, density of carbon fiber 
    
    rho_spar = .0644; %kg/m - linear density of tailboom
    rho_tail_boom = .055; %kg/m - linear density of tailboom
    
    rho_g_f = 0.047468; %kg/m^2 - 1.4 oz/yd^2 fiberglass - fuselage
    rho_g_t = 0.023734; %kg/m^2 - .7 oz/yd^2 fiberglass - tails
    fiber_resin_ratio = 1; %assume fiber and resin weight equal
    
    rho_f = 28.03231; %kg/m^3 - Dow blue foam
    
    %%%% CALCULATE FUSELAGE WEIGHT - Currently assumes a cylinder fuselage
    %%%% - Update 5/12/2021: Fuselage weight is now estimated using CAD. This section isn't used.    
    
    %assume fuselage wall thickness is 1 cm of foam
    %To account for fuselage structure, assume fuselage covered in 1.4 oz
    %fiberglass
    
%     rho_g10 = 1800; %kg/m^3
%     n = 4; %number bulkheads;
%     m_bulkheads = n*(pi*cst.r_fus^2)*0.0015875* rho_g10;
%     
%     l_fus = cst.l_fus; %ballpark for fuselage length.
%     r_fus = cst.r_fus; %assume fuselage fineness ratio of 8
%     wall_thick = 0.01; %assume fuselage wall thickness is .5 cm of foam
%     vol_fuse = ((pi.*r_fus^2)-(pi.*(r_fus-wall_thick)^2))*l_fus; 
%     sa_fuse = (2*pi*r_fus)*l_fus;
%     
%     m_fuse_glass = sa_fuse* rho_g_f;
%     m_fuse_glass = m_fuse_glass + m_fuse_glass/fiber_resin_ratio; %account for expoxy weight
%     m_fuse =  vol_fuse*rho_f + m_fuse_glass + m_bulkheads;
    
    %%%% Calculate tailboom weight
    
    m_boom = (l_t - (c*.75)+.18*cst.l_fus)*rho_tail_boom; %Assume the boom only extends 18% into fuselage. From the CAD
    
    %%%% WING CALCS, CURRENTLY ASSUMING RECTANGULAR SECTIONS. IF TAPERED
    %%%% DESIRED, JUST ADD MORE INPUTS TO FUNCTION AND EDIT
    %%%% cr/ct
    
    c = S_ref/b;
    cr = c;
    ct = c;
    
    %Calculate Wing Spar weight
    
    m_w_spar = rho_spar*b; % new calculation from the CAD. Wing spar has been chosen online
    
    xvec = linspace(0,b/2,10000); %m, x positions along wing starting from root to b/2
    cvec = ((b/2 - xvec)*cr + xvec*ct)*2/b; %m, chord at each x position along wing starting from root to b/2

    %%% USE AIRFOIL READER TO CALCULATE WING VOLUME, THEN MASS
    
    %plot(xin,yin);
    a = polyarea(cst.xin,cst.yin); %m^2/m (chord)
    vol = 2*trapz(xvec,a*cvec.^2); %m^3, total wing volume
    m_wing = vol * rho_f +  m_w_spar; %kg, wing mass (assumes solid x section)
    
    
    %%%% TAIL CALCS, CURRENTLY ASSUMING RECTANGULAR SECTIONS. IF TAPERED
    %%%% DESIRED, JUST ADD MORE INPUTS TO FUNCTION AND EDIT
    %%%% crt/ctt/crtv/cttv
    
    bt = s_htail/c_htail; %horizontal tail span
    crt = c_htail; %horizontal root chord
    ctt = c_htail; %horizontal tip chord
    
    btv = s_vtail/c_vtail; %vert tail span
    crtv = c_vtail; %vert root chord
    cttv = c_vtail; %vert tip chord
    
    %TAIL WEIGHT, HORIZONTAL
    xtvec = linspace(0,bt/2,10000); %m, x positions along tail starting from root to bt/2
    ctvec = ((bt/2 - xtvec)*crt + xtvec*ctt)*2/bt; %m, chord at each x position along tail starting from root to b/2

    %plot(xtin,ytin);
    at = polyarea(cst.xtin,cst.ytin); %m^2/m (chord)
    volt = 2*trapz(xtvec,at*ctvec.^2); %m^3, total tail volume
    
    %Calculate horizontal tail structure - assume covered in .7 oz
    %fiberglass
    sa_ht = 2*s_htail;
    m_ht_glass = sa_ht*rho_g_t;
    m_ht_glass = m_ht_glass + m_ht_glass/fiber_resin_ratio;

    %TAIL WEIGHT, VERTICAL
    xtvvec = linspace(0,btv,10000); %m, x positions along tail starting from root to btv
    ctvvec = ((btv - xtvvec)*crtv + xtvvec*cttv)/btv; %m, chord at each x position along tail starting from root to b/2
    
    atv = polyarea(cst.xtin,cst.ytin); %m^2/m (chord)
    voltv = trapz(xtvvec,atv*ctvvec.^2); %m^3, total tail volume
    
    %Calculate vertical tail structure - assume covered in .7 oz
    %fiberglass
    sa_vt = 2*s_vtail;
    m_vt_glass = sa_vt*rho_g_t;
    m_vt_glass = m_vt_glass + m_vt_glass/fiber_resin_ratio;

    m_htail = volt * rho_f + m_ht_glass; %kg, tail horizontal mass (assumes solid x section)
    m_vtail = voltv * rho_f + m_vt_glass; %kg, tail vertical mass (assumes solid x section)
    
    %mass of fuselage + electronics
    m_fixed = 0.425;  %fuselage weight estimated from CAD + payload

    mass_empty = m_boom + m_wing + m_htail + m_vtail + m_fixed;

end

function [CG, x_check] = calc_CG(S_ref, b, l_t, m_fixed, m_wing, m_htail, m_vtail, m_boom, mass_motor, mass_bat)
    global cst
    [c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b, l_t);
    
    % calculate the CG. Assume the zero is at wing leading edge. Postive
    % going towards tail
    x_fixed = ( c*.75 - cst.l_fus/2); %Assume fuselage ends at wing trailing edge. Assume CG fuselage at center.
    w_fixed = m_fixed*x_fixed;
    
    x_wing = c*.25; %Assume CG wing at .5 chord
    w_wing = m_wing*x_wing;
    
    x_htail = l_t + c_htail*.25; %Assume CG tail at .5 chord
    w_htail = m_htail*x_htail;
    
    x_vtail = l_t + c_vtail*.25; %Assume CG tail at .5 chord
    w_vtail = m_vtail*x_vtail;
    
    l_boom = l_t - (c*.75) + .18*cst.l_fus; %assume boom goes .18 into fuselage
    x_boom = l_t - l_boom/2; 
    w_boom = m_boom*x_boom;
    
    x_motor = (c*.75 - cst.l_fus); %assume motor attached to front of fuselage
    w_motor = mass_motor*x_motor;
    
    x_bat = (c*.75 - .385); %location of battery based off CAD
    w_bat = mass_bat*x_bat;
    
    x_check = [x_fixed x_wing x_htail x_vtail x_boom x_motor x_bat];
    mass_total = m_fixed  + m_wing + m_htail + m_vtail + m_boom + mass_motor + mass_bat;
    
    CG = (w_fixed + w_wing + w_htail + w_vtail + w_boom + w_motor + w_bat)/mass_total;
end

% Written by Yamaan, modified by Philip
function [sigma_max, deflection_span] = calc_beam(S_ref, weight, b)
    global cst
    
    %Spar Dimensions
    % https://www.clearwatercomposites.com/product/9-16-x-0-633-carbon-fiber-tube/
    r_o = 0.0160782/2; %m - 0.633 in OD
    r_i = 0.0143002/2; %m - 0.563 in ID
    
    if r_i < 0
        r_i = 0; %make sure r_i isn't negative
    end
    
    %Spar properties
    E = 70*10^9; %MPa taken from online carbon fiber spar
    I = pi*0.25* (r_o^4 - r_i^4); %m^4, derived from spar geometry
    J = pi*0.25* (r_o^4 - r_i^4); %m^4, derived from spar geometry
    ymax = r_o; %maximum spar distance from neutral axis

    %xlim([0 b/2]);
    %ylim([0 cr]);
    
    g_load = 3; %how much load do we expect?
    xvec = linspace(0,b/2,10000); %m, x positions along wing starting from root to b/2
    %wvec = cvec * L_S; %N/m, rectangular lift distribution, do not use
    w0 = (4*weight*g_load*cst.g)/(pi*b);
    wvec = w0 * sqrt(1 - (xvec/(0.5*b)).^2);%N/m, elliptical lift distribution
    wvec_flip = flip(wvec);
    
    Vvec = cumtrapz(xvec,wvec_flip); %N
    Mvec = cumtrapz(xvec,Vvec); %Nm
    
    Vvec = flip(Vvec);
    Mvec = flip(Mvec);
    
    u_primevec = cumtrapz(xvec,Mvec);
    uvec = cumtrapz(xvec,u_primevec)/(E*I); %m
    umax = max(uvec);
    
    deflection_span = umax/b; %deflection to span ratio
    
    V_max = max(Vvec);
    M_max = max(Mvec);
    %T_spar = L * d; neglecting for now in calculation as negligible

    sigma_max = M_max * ymax/I; %* 10^-6; %MPa, max compressive/tensile at root
    %tau_max =T_spar * r/J * 10^-6; %MPa, max shear due to torsion

    
    plots = 0; %toggle for debugging, need root and tip chords
    
    if plots == 1
        %chord as a fcn of span
        %figure(1)
        %fplot(@(x) ((b/2 - x)*cr + x*ct)*2/b, [0 b/2])
        %figure(2)
        %fplot(@(x) ((b/2 - x)*cr + x*ct)*2/b * L_S, [0 b/2])

        figure(3)
        subplot(3,1,1)
        plot(xvec,Vvec);
        xlabel('x (m)')
        ylabel('V(x) (N)')
        subplot(3,1,2)
        plot(xvec,Mvec);
        xlabel('x (m)')
        ylabel('M(x) (Nm)')
        subplot(3,1,3)
        plot(xvec,uvec);
        xlabel('x (m)')
        ylabel('u(x) (m)')
    end
end

% Written by Yamaan
function [n,xin,yin] = openfile(file)
    % Read airfoil profile and place in xin, yin vectors
    fid = fopen(file);
    n = fscanf(fid,'%10f',1); %length of airfoil vectors
    data = fscanf(fid,'%10f %10f \n', [2,n]);
    fclose(fid);
    xin = data(1,:);
    big=max(xin);
    xin=xin./big;
    yin = data(2,:);
    yin=yin./big;
end