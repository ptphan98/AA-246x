
weight = 1.5;
b = 1.3;
v_air = 20;

%% start of function
S_ref = .4;
c = S_ref./b;

%assume elliptical lift distribution
l_0 = (4*weight*9.81)/(pi*b); %make sure load in newtons
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
figure; plot(x,y); title('Lift Distribution'); 
xlabel('Span (m)'); ylabel('Load per unit span (N/m)');
figure; plot(x,shear); title('Shear');
xlabel('Span (m)'); ylabel('Shear (N)');
figure; plot(x,moment); title('Moment')
xlabel('Span (m)'); ylabel('Moment (N*m)');

%find max_stress
t_c = 0.08; %thickness to chord ratio of airfoil
r_spar = t_c*c/2; %assume spar is maximum thickness of airfoil 
I_c = pi/4*r_spar^4;   %2nd moment of inertia
stress = max(moment)*r_spar/I_c;

%find deflection
E = 228*10^9; %modulus of carbon fiber - 228 GPa

d2v_dx2 = moment/(E*I_c);
rotation = cumtrapz(x,d2v_dx2);
deflection = cumtrapz(x,rotation);

deflection_span = deflection(end)/b; %deflection to span ratio




