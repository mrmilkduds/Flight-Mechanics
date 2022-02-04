%This matlab program was written by Joshua Romans and James Stewart and finished on 3/16
%for submission to Dr. Hubner of the University of Alabama for grading in
%the course AEM368. 

%The purpose of this code is to supply the user with meaningful
%flight characteristics given a set of parameters to describe an aircraft.

%Loading Data from specs.txt and ensuring bounds are met------------------
close all;
data_array=load('specs1.txt'); % Order of Arrays should be order in which parameters are given in specs.txt 

if data_array(1) <= 30000 && data_array(1) >= 0 %ft
    altitude=data_array(1); 
elseif data_array(1) > 30000 
    data_array(1) = 30000;
else
    data_array(1) = 0;
end


if (1<= data_array(2)) && (data_array(2)<= 7) %ft
    chordRoot=data_array(2); 
elseif data_array(2) > 7
    chordRoot = 7;
else
    chordRoot = 1;
end


if (15 <= data_array(3)) && (data_array(3)<= 50) %ft
    span=data_array(3);
elseif data_array(3) > 50
    span = 50;
else
    span = 15;
end

    
if (0 <= data_array(4)) && (data_array(4) <=1)    %[]
    taperRatio=data_array(4);
elseif data_array(4) > 1
    taperRatio = 1;
else
    taperRatio = 0;
end

if (0 <= data_array(5)) && (data_array(5) <= 10) %degree (not used)
    sweep = data_array(5);
elseif data_array(5) > 10
    sweep = 10;
else
    sweep = 0;
end

if (1 <= data_array(6)) && (data_array(6) <=2) %[]
    clMax=data_array(6);
elseif data_array(6) > 2
    clMax = 2;
else
    clMax = 1;
end

if (0.015 <= data_array(7)) && (data_array(7) <= 0.04) %[]
    cdo=data_array(7);
elseif data_array(7) > 0.04
    cdo = 0.04;
else
    cdo = 0.015;
end

if (500 <= data_array(8)) && (data_array(8) <= 5000) %lb
    weightGross=data_array(8);
elseif data_array(8) > 5000
    weightGross = 5000;
else
    weightGross = 500;
end

if (20 <= data_array(9)) && (data_array(9) <= 40) %percentage
    fuelCap=data_array(9);
elseif data_array(9) > 40
    fuelCap = 40;
else
    fuelCap = 20;
end

if (.35 <= data_array(10)) &&(data_array(10) <= .6) %1/s
    BSFC=data_array(10);
elseif data_array(10) > .6
    BSFC = .6;
else
    BSFC = .35;
end

if (50 <= data_array(11)) && (data_array(11) <=500) %hp
    powerShaft=data_array(11); 
elseif data_array(11) > 50
    powerShaft = 500;
else
    powerShaft = 50;
end

if (.5 <=data_array(12))&&(data_array(12) <= .9) %[]
    effProp=data_array(12);
elseif data_array(12) > .9
    effProp = .9;
else
    effProp = .6;
end

%-------------------------------------------------------------


%Item 1-------------------------------------------------------

%SSL values to calculate values at altitude
go=32.2; 
R=1716.6; 
rhoSSL=.002377; 
tempSSL=59+459.67; 
pSSL=2116.2;
lapseRate=-.00356616; 

%Calculation at altitude and unit conversions
tempFahr=59+(lapseRate*altitude);   

tempAlt=tempFahr + 459.67; 
pAlt=pSSL*(tempSSL/tempAlt)^(go/(lapseRate*R)); %Pressure Trend equation 
rhoAlt=pAlt/(R*tempAlt); %Ideal Gas Law

%Item 2----------------------------------------------------

planformA=(((taperRatio*chordRoot)+(chordRoot))/2)*span; %area of two trapezoids

aspectRatio=(span)^2/planformA;

%Item 3---------------------------------------------------------------

%Tables from assignment sheet 
AspectRatioT= [20 16 12 8 4 2];
TaperRatioT=[1 ; .8; .6; .4; .2; 0];

Table2=[ .82, 0.882, 0.907, 0.937, 0.972, 0.990;
        .891, 0.909, 0.929, 0.952, 0.979, 0.993;
       .924, 0.938, 0.953, 0.969, 0.987, 0.995;
       .950, 0.960, 0.970, 0.980, 0.992, 0.997;
       .937, 0.942, 0.950, 0.962, 0.980, 0.991;
       .775, 0.783, 0.797, 0.820, 0.865, 0.913;];


e_span=interp2(AspectRatioT,TaperRatioT,Table2,aspectRatio,taperRatio); %2d Linear Interpolation of tables
e_oswald=.7;

%Item 4----------------------------------------------------------

k=(1/(pi*aspectRatio*e_oswald));

MaxL_D=((1/(k*cdo))^.5)/2; %CHECKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

%Item 5----------------------------------------------------------------

%Calculations of variables needed for Re,V@maxL/D
MAC= chordRoot *(2/3) * (( 1 + taperRatio + taperRatio^2 )/( 1 + taperRatio )) ; %M
nu_0=3.62*10^(-7); 
nu=nu_0*((tempAlt/tempSSL)^(3/2))*((tempSSL+198.72)/(tempAlt+198.72)); %At Altitude
 
VmaxL_D=((k/cdo)^(1/4))*((2*weightGross)/(rhoAlt*planformA))^(1/2); %CHECKKKKKKKKKKKKKKKKKKKK
Re=(rhoAlt*MAC*VmaxL_D)/nu;

%Item 6-------------------------------------------------------

%Assigning various weight variables for next few items
FuelCapacity=(fuelCap/100)*weightGross;
weightHalf=weightGross-(.5*FuelCapacity); 
W1=weightGross-FuelCapacity;
W0=weightGross;

WingLoading=weightHalf/planformA;

%Item 7----------------------------------------------------------------------

Vstall=((2*weightHalf)/(rhoAlt*planformA*clMax))^(1/2); 

%Item 8-------------------------------------------------------------------

Cpp=BSFC/(3600*550) ; %Unit Conversion to ft^-1           
Cl_Cd_RMAX=((pi*aspectRatio*e_oswald)/(4*cdo))^(1/2);

%Unit Conversion ft to mi
Range_max=(effProp/Cpp)*(Cl_Cd_RMAX)*log(W0/W1); %in ft
Range_miles=Range_max/5280; %mi

%Item 9----------------------------------------------------------------

%Cl^3/2/cd max
Cl_Cd3_2=(((3*pi*aspectRatio*e_oswald)^(3/4))/(4*cdo^(1/4)));

%Unit Conversion sec to hr
Epp_sec=(effProp/Cpp)*Cl_Cd3_2*((2*rhoAlt*planformA)^(1/2))*(((W1)^(-1/2))-((W0)^(-1/2)));
Epp_hr=Epp_sec/3600;



%Item 10 -----------------------------------------------------

%Unit and variable handling 
powerAvailable = powerShaft*effProp; %in HP
PowerA_Alt=powerAvailable*(rhoAlt/rhoSSL)*550; %now ftlb/s

%Coefficient Calculations for polynomial solution of Vmax
Coeff1=.5*rhoAlt*planformA*cdo;
Coeff2=(2*weightHalf^2)/(rhoAlt*planformA*pi*aspectRatio*e_oswald);
Coeff3=PowerA_Alt;

%Solving for real roots of polynomial
syms f(x);
f(x)=(Coeff1*x^4)+Coeff2-(Coeff3*x);
sol=vpasolve(f);
Realroot=real(sol);
a=(1.4*R*tempAlt)^(1/2);

Vmax=max(Realroot);
Mach=Vmax/a;

%Item 11,12,13 -------------------------------------------------------

%Table Calculates Pa, Pr, R/Cmax, Pex, and Absolute Ceiling by iterating h

%Some table handling stuff
hstep = 10;
ceilingbelow3000 = false;
for h=1:hstep:30000
    
    %Similar to item 1, alt replaced with h to iterate on
    temphFahr=59+(lapseRate*h);    %Temp Calculation, based on altitude in ft
    temph=temphFahr + 459.67; %Temp at alt in degrees Rankine
    ph=pSSL*(tempSSL/temph)^(go/(lapseRate*R)); %Pressure Trend equation psf
    rhoh=ph/(R*temph); %Density Calc, from Ideal Gas Law slug/ft^3
    powerAvailable_h = powerAvailable*(rhoh/rhoSSL)*550;
    
    %Table value assignment
    Table(h,1) = h; %altitude
    Table(h,2) = powerAvailable_h; %Power Available 
    Table(h,3) = (((2*(weightGross*weightGross*weightGross))/(rhoh*planformA))^(1/2))*(1/Cl_Cd3_2); %Power Required
    Table(h,4) = (effProp*powerAvailable_h/weightGross)-(0.8776*((weightGross/(planformA*rhoh*cdo))^(1/2))*(1/((Cl_Cd_RMAX)^(3/2)))); %R/CMax
    Table(h,5) = Table(h,2) - Table(h,3);
    Table(h,6) = (hstep/Table(h,4));
    
    %Find Absolute Ceiling
    if Table(h,5) < 0
        ceilingAbsolute = Table(h,1);
        ceilingbelow3000 = true;
        break
    end
end

%max of all RCmax should find RCmax across altitudes.
rcMax = max(Table(1:ceilingAbsolute,4));

%Basically integrating(read:summing) change in time between each height
%iteration
ttc= (sum( Table( 1:10000, 6 ) ))/60; %in mins

%Item14 ------------------------------------------------------------------

%Unit conversion ft to mi
rangeGlide = altitude*(Cl_Cd_RMAX);
rangeGlide_miles = rangeGlide/5280;
    
%Item 15 -----------------------------------------------------------------

%More table handling things
vMin = Vstall-5;
vMax = Vmax;
jMax = 100;
vStep = ((vMax-vMin)/jMax);

%probably the worst way to due this, put V and Pr in same table and plot
%the columns.
for j=1:jMax
    V = vMin + j*vStep;
    cl = (2*weightGross)/(rhoAlt*V*V*planformA);
    cd = cdo + ((cl*cl)/(aspectRatio*e_oswald*pi));
    cl_cd = (cl/cd);
    powerRequired(j,1) = V;
    powerRequired(j,2) = ((weightGross/cl_cd)*(((2*weightGross)/(rhoAlt*planformA*cl))^(1/2)))/550;
 
end

%figure
vmaxr=((k/cdo)^(.25))*((2*weightHalf)/(rhoAlt*planformA))^(.5);
vmaxe=((k/(3*cdo))^(.25))*((2*weightHalf)/(rhoAlt*planformA))^(.5);

 clR = (2*weightGross)/(rhoAlt*vmaxr*vmaxr*planformA);
    cdR = cdo + ((clR*clR)/(aspectRatio*e_oswald*pi));
    clR_cdR = (clR/cdR);

clE = (2*weightGross)/(rhoAlt*vmaxe*vmaxe*planformA);
    cdE = cdo + ((clE*clE)/(aspectRatio*e_oswald*pi));
    clE_cdE = (clE/cdE);
    
 clStall=(2*weightGross)/(rhoAlt*Vstall*Vstall*planformA);
    cdStall = cdo + ((clStall*clStall)/(aspectRatio*e_oswald*pi));
    clS_cdS = (clStall/cdStall);
powerReqStall=((weightGross/clS_cdS)*(((2*weightGross)/(rhoAlt*planformA*clStall))^(1/2)))/550;    
powerReqR=((weightGross/clR_cdR)*(((2*weightGross)/(rhoAlt*planformA*clR))^(1/2)))/550; 
powerReqE=((weightGross/clE_cdE)*(((2*weightGross)/(rhoAlt*planformA*clE))^(1/2)))/550;
    
%plot(powerRequired(j,1), powerRequired(j,2))
v=powerRequired(:,1);
PR=powerRequired(:,2);
plot(v,PR)
title('Power Required vs. Speed')
xlabel('Speed (ft/s)')
ylabel('Power (hp)') 
hold on;
plot(vmaxr,powerReqR,'bo')
hold on;
plot(vmaxe,powerReqE,'o')
hold on;
plot(Vstall,powerReqStall,'go')
hold on;
legend('Power Required','Max Range', 'Max Endurance / Rate of Climb', 'Stall') 
%OUTPUTS-------------------------------------------------------------------
fprintf('Temperature: %.4g [degF]\n',tempFahr)
fprintf('Pressure: %.4g [psf]\n',pAlt)
fprintf('Density: %.4g [sl/ft^3]\n',rhoAlt)
fprintf('Aspect Ratio: %.4g\n',aspectRatio)
fprintf('Spanwise Efficiency Factor: %.4g \n',e_span)
fprintf('Oswald Efficiency Factor: %.4g \n',e_oswald)
fprintf('Reynolds Number: %.4g \n',Re)
fprintf('Flight Speed at Max L/D: %.4g [ft/s] \n',VmaxL_D)
fprintf('Wing Loading at 50%% Fuel Consumption: %.4g [lb/ft^2] \n',WingLoading)
fprintf('Stall Speed at 50%% Fuel Consumption: %.4g [ft/s] \n',Vstall)
fprintf('Maximum Endurance: %.4g [hr] \n',Epp_hr)
fprintf('Corresponding Speed at 50%% Fuel Consumption: %.4g [ft/s]\n',tempAlt) %CHECK
fprintf('Maximum Range: %.4g [mi] \n',Range_miles)
fprintf('Corresponding Speed at 50%% Fuel Consumption: %.4g [ft/s]\n',tempAlt) %CHECK
fprintf('Maximum Speed: %.4g [ft/s] \n',Vmax)
fprintf('Corresponding Mach Number: %.4g \n',Mach)
if ceilingbelow3000 == false
    fprintf('Absolute Flight Ceiling at MTOW: Greater than 30000 [ft] \n',ceilingAbsolute)
else
    fprintf('Absolute Flight Ceiling at MTOW: %.4g [ft] \n',ceilingAbsolute)
end

fprintf('Maximum Rate of Climb at MTOW: %.4g [ft/s]\n',rcMax) 
fprintf('Time to Climb to 10000 ft from Sea Level: %.4g [min] \n',ttc) 
fprintf('Maximum Glide Range: %.4g [mi]\n',rangeGlide_miles) 



                                    

