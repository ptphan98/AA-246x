global cst 
cst = struct();
cst.rho = 1.225; %kg/m^3
cst.g = 9.807; %m/s^2
cst.CL_max = 1.2;
cst.V_stall = 7; %m/s
cst.W_L = 1/2 * cst.rho * cst.V_stall^2 * cst.CL_max / cst.g; % wing loading is sized by stall speed
cst.spar_ratio = .5; %percent spar of max airfoil thickness

[sigma_max, deflection_span] = calc_beam(.4, 1.5, 1.3)
function [sigma_max, deflection_span] = calc_beam(S_ref, weight, b)
    global cst
    
    %Spar Dimensions
    t_c = .08; %assume a thickness to chord for the wing
    c = S_ref/b; %wing chord - m
    r_o = cst.spar_ratio*t_c*c/2; %m - estimate a thickness for the spar. 
    r_i = r_o - 0.0015875; %m - assume a constant wall thickness of 1/16 inch - from Dragonplate's website
    
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
    
    g_load = 1; %how much load do we expect?
    xvec = linspace(0,b/2,10000); %m, x positions along wing starting from root to b/2
    %wvec = cvec * L_S; %N/m, rectangular lift distribution, do not use
    w0 = (4*weight*g_load*cst.g)/(pi*b);
    wvec = w0 * sqrt(1 - (xvec/(0.5*b)).^2);%N/m, elliptical lift distribution
    wvec_flip = flip(wvec);
    
    Vvec = cumtrapz(xvec,wvec_flip); %N
    Mvec = cumtrapz(xvec,Vvec); %Nm
    
    Vvec = flip(Vvec);
    Mvec = flip(Mvec);
    
    figure;
    plot(xvec,wvec)
    title('lift')
    figure;
    plot(xvec,Vvec)
    title('shear')
    figure;
    plot(xvec,Mvec)
    title('moment')
    
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