%% Determine the aperture size D to achieve a required directivity for a given frequency


% For the horn antenna on HFSS the directivity is 13.4dB so q =5 from the Aperure _Efficiency_Analysis paper Fig 3(a)
clc
clear all

f = 10*10^9;
c = 3*10^8;
lamda =1000*(c/f);
k_0 = 2*pi/lamda;
%% ----------------Given an aperture compute the max Dir------------------- run the first section initially
%Set the maximum dimension of the suqare apperture in mm.
D=90;
%Compute the maximum directivity that can be achieved with the predefined
%aperture D
Dir = 10*log10(4*pi* D^2/(lamda^2))
disp(['To achieve a directivity of ',num2str(Dir),'dB, the size of the rectangular aperture should be ',num2str(D),' mm'])

%% ----------------Assuming a maximum Dir find the required aperture size-------------------run the first section initially
%Assume a directivity Dir
Dir=10^(29/10);
%Aperture area A - all units in mm and Dir in dB




D=sqrt(A);

disp(['To achieve a directivity of ',num2str(Dir),'dB, the size of the rectangular aperure should be ',num2str(D),' mm'])
%% Formulas
%ns is the spillover efficiency
%ni is the illumination efficiency
%n_ap=ni*ns
%G_max = D_max*n_ap;

% ns=(2*qf+1)/(2*pi)
% 
% s=sqrt(((x-x0)^2)+((y-y0)^2));
% r=sqrt(x^2+y^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y));
% r0 = sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0));
% feed_pattern = (r0^2+r^2-s^2)/(2*r0*r);

%% Effects of offset angle for a given aperure size D, feed distance R0, feed pattern q,element pattern qe

D = 165; %apperture size
R0 = 117; %Set phase center(feed distance)
F_D = R0/D; %F/D ratio
q=5; % corresponds to feed pattern of 13 dB
qe=1; % corresponds to element pattern of 7.78dB
% theta_0=deg2rad();
x0=0;
y0=0;
%--------------------------------------------------------------------------------------------------------------------
%Calculate the spillover efficiecny for different value of theta_0 which
%correspond to the feed offset

l=1; % l is used to save the integration result in a vector for different values of theta_0
for i = 0:40
theta_0=deg2rad(i);

%spillover numerator equation (10),(12)
spillover_fun=@(x,y)(R0./(sqrt(x.^2+y.^2+(R0.^2).*(sec(theta_0).^2)+(2*R0.*tan(theta_0).*y))).^3).*((sqrt(x0^2+y0^2+(R0^2).*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0)).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2*sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y))))).^(2*q);
In= integral2(spillover_fun,-D/2,D/2,-D/2,D/2);


%spillover denominator equation (7),(9)
ns(l)=((2*q+1)/(2*pi))*In;
l=l+1;
end

%---------------------------------------------------------------------------------------------------------------------
%Calculate the illumination efficiecny for different value of theta_0 which
%correspond to the feed offset
m=1;
for j = 0:40
theta_0=deg2rad(j);
%illumination nominator equation (19),(20)
illum_fun_n = @(x,y) ((R0^qe)./(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^(1+qe)).*(((sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2.*(sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))))).^q;
n_in = integral2(illum_fun_n,-D/2,D/2,-D/2,D/2);

%illumination denominator equation (19),(20)
illum_fun_d = @(x,y) ((R0^2*qe)./(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^(2+2*qe)).*(((sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2.*(sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))))).^(2*q);
d_in = integral2(illum_fun_d,-D/2,D/2,-D/2,D/2);


%total illumination equation (19)
%ni(m) = 4*n_in^2/(pi*(D^2)*d_in)-0.169;  %% have to be fixed
ni(m) = n_in^2/((D^2)*d_in);
m=m+1;

end


% 40 is the maximum theta offset that is investigated
plot(0:40,ns)
hold on
plot(0:40,ni)
hold on
plot(0:40,ni.*ns)

title("Efficiency of Reflectarray for different feed positions")

legend("Spillover efficiency","Illumination efficiency","Apperture efficiency")

grid on

%% Effects of F/D for a given feed offset theta_0,give aperure size D,feed pattern q,element pattern qe

% For the horn antenna on HFSS the directivity is 13.4dB so q =5 from the Aperure _Efficiency_Analysis paper Fig 3(a)
clc
clear all
close

D = 165; %apperture size
R0_init = D/2; %Set phase center(feed distance)
q=4; % corresponds to feed pattern of 13 dB
qe=1; % corresponds to element pattern of 7.78dB
theta_0=deg2rad(0);%the feed offset can be set with theta_0
x0=0;
y0=0;


%-------------------------------------------------------------------------------------------------------------------
%Calculate the spillover efficiecny for different value of feed distance R0 starting from R0_init=D/2 up to D which
%correspond to the F/D from 0.5 to 1



l=1; % l is used to save the integration result in a vector for different values of theta_0
for R0 = R0_init:D


%spillover numerator equation (10),(12)
spillover_fun=@(x,y)(R0./(sqrt(x.^2+y.^2+(R0.^2).*(sec(theta_0).^2)+(2*R0.*tan(theta_0).*y))).^3).*((sqrt(x0^2+y0^2+(R0^2).*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0)).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2*sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y))))).^(2*q);
In= integral2(spillover_fun,-D/2,D/2,-D/2,D/2);


%spillover denominator equation (7),(9)
ns(l)=((2*q+1)/(2*pi))*In;
l=l+1;
end

%----------------------------------------------------------------------------------------------------------------------
%Calculate the illumination efficiecny for different value of feed distance R0 starting from R0_init=D/2 up to D which
%correspond to the F/D from 0.5 to 1
m=1;
for R0 = R0_init:D

%illumination nominator equation (19),(20)
illum_fun_n = @(x,y) ((R0^qe)./(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^(1+qe)).*(((sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2.*(sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))))).^q;
n_in = integral2(illum_fun_n,-D/2,D/2,-D/2,D/2);

%illumination denominator equation (19),(20)
illum_fun_d = @(x,y) ((R0^2*qe)./(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^(2+2*qe)).*(((sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2.*(sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))))).^(2*q);
d_in = integral2(illum_fun_d,-D/2,D/2,-D/2,D/2);


%total illumination equation (19)
%ni(m) = 4*n_in^2/(pi*(D^2)*d_in)-0.169;  %% have to be fixed
ni(m) = n_in^2/((D^2)*d_in);
m=m+1;

end



plot(R0_init:D,ns)
hold on
plot(R0_init:D,ni)
hold on
plot(R0_init:D,ni.*ns)

title("Efficiency of Reflectarray for different F/D")

legend("Spillover efficiency","Illumination efficiency","Apperture efficiency")

ylabel("Efficiency %")
xlabel("Focal Distance F(mm)")
grid on

n_ap = ni.*ns;
disp(['Maximum Aperture efficiency is ',num2str(100*max(n_ap)),'%'])

%find the maximum and the index
[pos] = find(ismember(n_ap, max(n_ap(:))));
%use the index in the possible feed distance to address the correct feed distance
ss = R0_init:D;
disp(['With Feed distance of ',num2str(ss(pos)),'mm'])

grid on

%% Effects of feed pattern q for a given feed offset theta_0,give aperure size D,feed distance R_0,element pattern qe

clc
clear all



D = 165; %apperture size
R0 = 117; %Set phase center(feed distance)
F_D = R0/D; %F/D ratio
qe=1; % corresponds to element pattern of 7.78dB
theta_0=deg2rad(0);%the feed offset can be set with theta_0
x0=0;
y0=0;
%----------------------------------------------------------------------------------------------------------------------
%Calculate the spillover efficiecny for different value of feed pattern q from 2 to 10
l=1;
for q = 2:10

spillover_fun=@(x,y)(R0./(sqrt(x.^2+y.^2+(R0.^2).*(sec(theta_0).^2)+(2*R0.*tan(theta_0).*y))).^3).*((sqrt(x0^2+y0^2+(R0^2).*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0)).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2*sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y))))).^(2*q);
In= integral2(spillover_fun,-D/2,D/2,-D/2,D/2);

ns(l)=((2*q+1)/(2*pi))*In;
l=l+1;
end
%----------------------------------------------------------------------------------------------------------------------
%Calculate the illumination efficiecny for different value of feed pattern q from 2 to 10
m=1;
for q = 2:10
%illumination nominator equation (19),(20)
illum_fun_n = @(x,y) ((R0^qe)./(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^(1+qe)).*(((sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2.*(sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))))).^q;
n_in = integral2(illum_fun_n,-D/2,D/2,-D/2,D/2);

%illumination denominator equation (19),(20)
illum_fun_d = @(x,y) ((R0^2*qe)./(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^(2+2*qe)).*(((sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).^2+(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))).^2-(sqrt(((x-x0).^2)+((y-y0).^2))).^2)./(2.*(sqrt(x0^2+y0^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0)*y0))).*(sqrt(x.^2+y.^2+(R0^2)*(sec(theta_0)^2)+(2*R0*tan(theta_0).*y))))).^(2*q);
d_in = integral2(illum_fun_d,-D/2,D/2,-D/2,D/2);


%total illumination equation (19)
%ni(m) = 4*n_in^2/(pi*(D^2)*d_in)-0.169;  %% have to be fixed
ni(m) = n_in^2/((D^2)*d_in);
m=m+1;

end


n_ap = ni.*ns;
plot(2:10,ni.*ns)
% hold on
% plot(2:10,ni)
% hold on
% plot(2:10,ns)

title("Apperture Efficiency of Reflectarray for different feeds")
ylabel("Efficiency %")
xlabel("Feed pattern q")
grid on

