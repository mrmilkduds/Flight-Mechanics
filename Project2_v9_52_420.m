%This matlab program was written by Joshua Romans, James Stewart, Raza
%Bajwa, and Emily Bell and finished on 4/19
%for submission to Dr. Hubner of the University of Alabama for grading in
%the course AEM368. 

%The purpose of this code is to supply the user with meaningful
%flight characteristics given a set of parameters to describe an aircraft.

%Loads all data from the three text files
clear all;
close all;
data_array1=load('specs.txt');% Order of Arrays should be order in which parameters are given in specs.txt 
data_array2=load('wing.txt');% Order of Arrays should be order in which parameters are given in wing.txt 
data_array3=load('tail.txt');% Order of Arrays should be order in which parameters are given in tail.txt 

%Checks each data set and assigns all variables
data_array1 = check1(data_array1); 
    altitude=data_array1(1);
    Mach=data_array1(2);
    cdofuse=data_array1(3);
    
data_array2 = check2(data_array2);
    wingarea=data_array2(1);
    ARwing=data_array2(2);
    taperwing=data_array2(3); 
    sweepwingc_4=data_array2(4);
    thickr_wing=data_array2(5);
    xcg_cr_wing=data_array2(6);
    zcg_cr_wing=data_array2(7);
    cm_ac_wb=data_array2(8);

data_array3 = check3(data_array3);  
    tailareapct=data_array3(1);
    ARtail=data_array3(2);
    tapertail=data_array3(3);
    sweeptail=data_array3(4);
    thickr_tail=data_array3(5);
    lt_b_wing=data_array3(6);
    zt_cr_wing=data_array3(7);
    tailangle=data_array3(8);
    dele_dela=data_array3(9);

%--------------------------------------------------------------

%Item 1-------------------------------------------------------

%SSL values to calculate values at altitude
go=32.2; %ft/s^2
R=1716.6; %
rhoSSL=.00237717; %slug/ft3
tempSSL=59+459.67; %degrees Rankine
pSSL=2116.2; %psf
lapseRate= -.00356616; %degrees Rankine/ft

%Calculation at altitude and unit conversions
tempFahr=59+(lapseRate*altitude);   

tempAlt=tempFahr + 459.67;
pAlt=pSSL*(tempSSL/tempAlt)^(go/(lapseRate*R)); %Pressure Trend equation 
rhoAlt=pAlt/(R*tempAlt); %Ideal Gas Law

sweepwingLE=sweepwingc_4+((4/ARwing)*((0-.25)*((1-taperwing)/(1+taperwing))));

%Item 2------------------------------------------------------

%Value of a to get M*a=v
a=(1.4*R*tempAlt)^.5; %ft/s
v=Mach*a; %ft/s

%Item 3------------------------------------------------------

muSSL=3.62*10^-7; %psf
muAlt=muSSL*((tempAlt/tempSSL)^(1.5))*((tempSSL+198.72)/(tempAlt+198.72));
b=(wingarea*ARwing)^.5;
cr=(2*wingarea)/(b*(taperwing+1));
wingMac= (2/3)*cr*((1+taperwing+taperwing^2)/(1+taperwing));
reynolds_wing=(rhoAlt*wingMac*v)/(muAlt);

%Item4-------------------------------------------------------

wingmac_rootchord=wingMac/cr;

%Item5-------------------------------------------------------
y_bar_bar=(b*((2*taperwing)+1))/(6*(taperwing+1));
x_ac=(wingMac/4)+(y_bar_bar*tand(sweepwingLE));
x_ac_cr=x_ac/cr;
h_ac=x_ac/wingMac;

%Item6-------------------------------------------------------
%wing lift curve slope
%utilizes equation for lift curve slope from Lec 2.4.2

AspectRatioT= [ 16 12 8 4 2];
TaperRatioT=[1 ; .8; .6; .4; .2; 0];
Table1= [0.2956 0.2511 0.1952 0.1195 0.0682;
         0.2252 0.1908 0.1478 0.0901 0.0513;
         0.1527 0.1288 0.0991 0.0598 0.0338;
         0.0933 0.0789 0.0610 0.0369 0.0208;
         0.0939 0.0846 0.0715 0.0500 0.0318;
         0.3236 0.3091 0.2850 0.2343 0.1778;];

tauw=interp2(AspectRatioT,TaperRatioT,Table1,ARwing,taperwing); 

awing = (2*pi())/(1+((2*pi())/pi()*ARwing)*(1+tauw)); %in 1/rad
awing = awing*(pi()/180); %in 1/deg

%Item7-------------------------------------------------------
%tail lift curve slope
%same procedure as item 6 but with the tail values instead

taut=interp2(AspectRatioT,TaperRatioT,Table1,ARtail,tapertail);

atail = (2*pi())/(1+((2*pi())/pi()*ARtail)*(1+taut)); %in 1/rad
atail = atail*(pi()/180); %in 1/deg

%Item8-------------------------------------------------------
%Tail volume ratio with mac

tailarea=(tailareapct/100)*wingarea;
cr_tail=(2*tailarea)/(b*(tapertail+1));
tailMac=(2/3)*cr_tail*((1+tapertail+tapertail^2)/(1+tapertail));

lt=(.5*b)*(lt_b_wing);

Vh=((lt)*(tailarea))/(wingarea*wingMac); 


%LE of Tail
btail=(tailarea*ARtail)^.5;
y_bar_bar2=(b*((2*tapertail)+1))/(6*(tapertail+1));
x_ac_tail=(tailMac/4)+(y_bar_bar2*tand(sweeptail));

xcg_tail = cr_tail*xcg_cr_wing;
hcg_tail = xcg_tail/tailMac;  
cg_to_le_tail=lt-x_ac_tail;
LETail=xcg_tail+cg_to_le_tail;


%Item9-------------------------------------------------------
%Trim AoA
%Finding intersect of Cm v AoA plot from Lec 7.2
%treating as linear equation y = mx + b. Where y = zero because of trimmed
%condition. Equations found in example 7.5 from the book

xcg = cr*xcg_cr_wing; %hcg  
hcg=xcg/wingMac;
h_ac=x_ac/wingMac;
cm0=cm_ac_wb + Vh*atail*tailangle; %b value

cmcg_alpha = awing*((hcg - h_ac) - (Vh*(atail/awing)*(1-dele_dela))); %slope
alpha_trim_check = -cm0/cmcg_alpha; %treating alpha as x and rearanging to solve.

alpha2=(cm_ac_wb+(atail*Vh*tailangle))/(awing*(hcg-h_ac)-(atail*Vh*(1-dele_dela)));
alpha_trim=-alpha2;


%Item10-------------------------------------------------------
%static margin
Smargin=(h_ac-hcg)+((atail/awing)*Vh*(1-dele_dela));
zcg = zcg_cr_wing*cr;  %Used later in Cmcg




%Item11-------------------------------------------------------
%Lift Coefficient at trim AoA

clw_trim = awing*alpha_trim; %lift coefficient is equivalent to lift curve slope times AoA of the wingcl
clt_trim = (atail*alpha_trim*(1-dele_dela)) -atail*tailangle;
cl_trim = clw_trim + clt_trim;%*(tailarea/wingarea)

%Item12-------------------------------------------------------
%Weight and Wing loading at trim AoA
%when trimmed W=L

W = cl_trim*.5*rhoAlt*(v^2)*wingarea; %equation for lift
W_loading = W/wingarea;


%Table for e-oswald and e-spanwise

Table2=[0.882, 0.907, 0.937, 0.972, 0.990;
        0.909, 0.929, 0.952, 0.979, 0.993;
        0.938, 0.953, 0.969, 0.987, 0.995;
        0.960, 0.970, 0.980, 0.992, 0.997;
        0.942, 0.950, 0.962, 0.980, 0.991;
        0.783, 0.797, 0.820, 0.865, 0.913;];


e_span_wing=interp2(AspectRatioT,TaperRatioT,Table2,ARwing,taperwing);
e_os_wing=.75*e_span_wing;

e_span_tail=interp2(AspectRatioT,TaperRatioT,Table2,ARtail,tapertail);
e_os_tail=.75*e_span_tail;



%Print all variables with proper units

fprintf('Temperature: %.4g [degF]\n',tempFahr)
fprintf('Pressure: %.4g [psf]\n',pAlt)
fprintf('Density: %.4g [sl/ft3]\n',rhoAlt)
fprintf('Flight Speed: %.4g [ft/s]\n',v)
fprintf('Reynolds #: %.4g \n',reynolds_wing)
fprintf('Wing mac / Root Chord: %.4g \n',wingmac_rootchord)
fprintf('Root Chord LE and Wing AC Ratio: %.4g \n',x_ac_cr)
fprintf('Wing Lift Curve Slope: %.4g [1/deg] \n',awing)
fprintf('Tail Lift Curve Slope: %.4g [1/deg] \n',atail)
fprintf('Tail Volume Ratio: %.4g \n',Vh)
fprintf('Trim Angle of Attack: %.4g [deg]\n',alpha_trim)
fprintf('Static Margin: %.4g [ft]\n',Smargin)
fprintf('Lift Coeff at Trim: %.4g \n',cl_trim)
fprintf('Aircraft Weight: %.4g [lbf]\n',W)
fprintf('Wing Loading: %.4g \n',W_loading)



hold off
%top down plot
figure('name','Wing and Tail to scale')
xtable=[(0),((b/2)*tand(sweepwingLE)), ((b/2)*tand(sweepwingLE)+(taperwing*cr)),(cr)];
ytable=[0,(-b/2),(-b/2), (0)];
plot(xtable, ytable,'b')
hold on

xtable2=[(0),((b/2)*tand(sweepwingLE)), ((b/2)*tand(sweepwingLE)+(taperwing*cr)),(cr)];
ytable2=[0,(b/2),(b/2), (0)];
plot(xtable2, ytable2,'b')
hold on

xtable3=[(LETail),((LETail+(btail/2)*tand(sweeptail))), (LETail+(btail/2)*tand(sweeptail)+(tapertail*cr_tail)),(LETail+cr_tail)];
ytable3=[0,(-btail/2),(-btail/2), (0)];
plot(xtable3,ytable3,'b')

hold on

plot(x_ac,0,'go')
plot(xcg,0,'ro')
hold on
xtable4=[(LETail),((LETail+(btail/2)*tand(sweeptail))), (LETail+(btail/2)*tand(sweeptail)+(tapertail*cr_tail)),(LETail+cr_tail)];
ytable4=[0,(btail/2),(btail/2), (0)];
plot(xtable4,ytable4,'b')

xle=[0, (LETail+cr_tail)];
yle=[0,0];
plot(xle,yle,'b')

%CM_CG
alphatable=[-6:10];
figure('name','cmcg vs AoA');
cm_cg= cm_ac_wb+(atail*Vh*tailangle)+alphatable*((awing*(hcg-h_ac))-(atail*Vh*(1-dele_dela)));
plot(alphatable,cm_cg)

%lift drag polar plot------------------------------------------------------

bt=(tailarea*ARtail)^.5;
cr=(2*tailarea)/(bt*(tapertail+1));
tailMac= (2/3)*cr*((1+tapertail+tapertail^2)/(1+tapertail));
reynolds_tail=(rhoAlt*tailMac*v)/(muAlt);

%CFwing=(.455/(log(reynolds_wing^2.58))-(8700/reynolds_wing));
%DF=CFwing*rhoAlt*v*v*wingarea*wingMac;
%Cwingform_1=1;
%Ctailform_1=1;

kwing=(5.46*(thickr_wing)^2)+(1.55-sind(sweepwingLE))*(thickr_wing)+1; %Form factor correction
ktail=(5.46*(thickr_tail)^2)+(1.55-sind(sweeptail))*(thickr_tail)+1;

AoA = -6;
ldpolar = zeros(5,17);
q = .5*rhoAlt*(v^2);
for i=1:17
clw = awing*AoA;
clt = (atail*AoA*(1-dele_dela)) -atail*tailangle;
cl = clw + clt;%*(tailarea/wingarea);

inducedwing = (clw^2)/(pi()*ARwing*e_os_wing);
inducedtail = (clt^2)/(pi()*ARtail*e_os_tail); 

if reynolds_wing > 500000
    frictionwing = (.455/((log(reynolds_wing))^2.58)) - (1700/reynolds_wing);
   else
    frictionwing = 1.33/(sqrt(reynolds_wing));
end

if reynolds_tail > 500000
    frictiontail = (.455/((log(reynolds_tail))^2.58));
   else
    frictiontail = 1.33/(sqrt(reynolds_tail));
end

Dfuse = cdofuse*q*wingarea;
DInducedWing = inducedwing*q*wingarea;
DInducedTail = inducedtail*q*tailarea;
DFrictionWing = frictionwing*kwing*q*2*wingarea;
DFrictionTail = frictionwing*ktail*q*2*wingarea;


Dtotal = Dfuse + DInducedWing + DInducedTail + DFrictionWing + DFrictionTail;
cd = Dtotal/q/wingarea;
%cd = cdofuse + inducedwing + inducedtail + frictionwing*kwing + frictionwing*ktail;
cl_cd = cl/cd;

Lw = clw*q*wingarea;
Lt = clt*q*tailarea;
Dw = Dfuse + DInducedWing + DFrictionWing;
%Dw = (cdofuse +inducedwing + frictionwing*kwing)*q*wingarea;




cl_cd = cl/cd;

Lw = clw*q*wingarea;
Lt = clt*q*tailarea;
Dw = (cdofuse +inducedwing + frictionwing*kwing)*q*wingarea;
Macw = cm_ac_wb*q*wingarea;
Mcgw = Macw + (Lw + Dw*AoA)*(xcg - x_ac)*wingMac + (Lw*AoA - Dw)*zcg;
Mcgt = -lt/Lt;

Mcg = Mcgw + Mcgt;
cmcg = Mcg/(q*wingarea*wingMac);
ldpolar(:,i) = [AoA,cl,cd,cl_cd,cmcg]; 
AoA = AoA+1;
end
figure('name', 'Lift Drag Polar')
plot(ldpolar(2,:),ldpolar(3,:))

%Cl/Cd v AoA plot----------------------------------------------------------

inducedwing_trim = (clw_trim^2)/(pi()*ARwing*e_os_wing);
inducedtail_trim = (clt_trim^2)/(pi()*ARtail*e_os_tail); 
cd_trim = cdofuse + inducedwing_trim + inducedtail_trim + frictionwing*kwing + frictionwing*ktail;
cl_cd_trim = cl_trim/cd_trim;

figure('name', 'Cl/Cd vs AoA')
plot(ldpolar(1,:),ldpolar(4,:))
hold on 
plot(alpha_trim, cl_cd_trim, 'go')
hold off

figure('name', 'Cmg (?) vs AoA')
plot(ldpolar(1,:),ldpolar(5,:))

%Set of functions that will check all paramaters 
function data_array1 = check1(data_array1) %Function to check specs bounds

min = [0;.05;.01];
max = [30000;.4;.05]; %arrays containing the minimum and max values of each spec

for i=1:3 %loop that runs through and checks all 12 specs
    if data_array1(i) < min(i)     %checks to see if spec is below minimum value
        data_array1(i) = min (i);   %if below corrects spec to be the minimum value
    elseif data_array1(i) > max(i) %checks to see if spec is above maximum value
        data_array1(i) = max(i);    %if above corrects spec to be the maximum value
    end
end
end 
function data_array2 = check2(data_array2) %Function to check wing bounds

min = [25;4;.2;0;.05;-2;-1;-.2];
max = [250;16;1;20;.15;3;1;0]; %arrays containing the minimum and max values of each spec

for i=1:8 %loop that runs through and checks all 12 specs
    if data_array2(i) < min(i)     %checks to see if spec is below minimum value
        data_array2(i) = min (i);   %if below corrects spec to be the minimum value
    elseif data_array2(i) > max(i) %checks to see if spec is above maximum value
        data_array2(i) = max(i);    %if above corrects spec to be the maximum value
    end
end
end 
function data_array3 = check3(data_array3) %Function to check tail bounds

min = [10;3;.4;0;.05;.5;.5;-5;0];
max = [40;3;1;10;.15;2;-1;5;1]; %arrays containing the minimum and max values of each spec

for i=1:9 %loop that runs through and checks all 12 specs
    if data_array3(i) < min(i)     %checks to see if spec is below minimum value
        data_array3(i) = min (i);   %if below corrects spec to be the minimum value
    elseif data_array3(i) > max(i) %checks to see if spec is above maximum value
        data_array3(i) = max(i);    %if above corrects spec to be the maximum value
    end
end
end 