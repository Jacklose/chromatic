%% focusing_sphere.m
% Copyright C. Stubbs Feb 20 2015
% computes focus accommodation for spherical lens

% all dimensions are in units of lens radius. Keep object at least
% one focal length away to get a real image. Typical octopus lens has
% f/1.2, so focal length of 2.4 lens radii. 

% actual 50% points in opsin response are basically at 450 and 550 nm.
% Peak of opsin is 500 nm, and focal length there is 11.87 mm

% Our fit to J&S gives focal shift at 550 of zero, by their convention.
% shift at 500 nm is -0.1289 mm, shift at 450 nm is -0.2936 mm
% so centering this at 500 gives shift to 550 of +0.1289 and to 450 it's
% -0.1647 mm. Converting these to a percentage of 11.87 we get
% redshift=0.1289/11.87 = 1.09 percent redshfit and 
% blueshift=0.1647/11.87= 1.39 percent blueshift

%% initialize
clear all
f=2.4 ; % focal length in lens radius units, middle of opsin response
blueshift = 1.39 % percent focal shift for blue wavelength. Here we use 450. Keep this a positive number
redshift = 1.09 % percent focus shift for red wavelength. Here we use 550
Obj=[8:0.2:10000]; % object distances in lens radius units
lensradius=5; % lens radius in mm, lens diameter is obviously twice this
blurtolerance=0.05 % this is faction of chromatic blur we'll accept as definition of depth of field. Ranges from 0 to 1. 

trigterm=cos(asin(1./Obj));

inverseI=(1./trigterm).*(1/f-1./Obj);

I=1./inverseI; % midband location of image


%% now compute same thing for opsin band edges, presumed to be 50% points
fr=(1+redshift*1E-2)*f % focal length in lens radius units, red light, 3% shift in focus
fb=(1-blueshift*1E-2)*f; % focal lenght in lens radius units, blue light, for 3% shift in focus across the band

trigterm=cos(asin(1./Obj));

inverseIr=(1./trigterm).*(1/fr-1./Obj);
inverseIb=(1./trigterm).*(1/fb-1./Obj);

Ir=1./inverseIr; % red edge image position
Ib=1./inverseIb % blue edge image position

%% compute focal shift as a function of position
dI=Ir-Ib; % difference in image position between red and blue


%% compute depth of field in comparison to chromatic focal shift. 
% first, take the numerical derivative of the focal shift vs. object
% distance. We'll do this for the blue light

dObj=Obj(2)-Obj(1); % find step size in O

diff_I=diff(I); % take numerical derivative of bset focus setting with respect to object distance
% this always returns a difference vector that is 1 shorter than the
% original input array, so tack one on the end to make it same length

diff_I(end+1)=diff_I(end) % it's pretty flat out at large O values so just replicate

diff_I=diff_I-0.0000001 % add on a small value to avoid numerical divergences at zero. Note it's negative so offset must be <0

% account for step size in O
dIdO=diff_I./dObj;

% Now let's pick a definition of depth of field to be the delta-O spacing
% that produces a focal shift that is half the chromatic shift
% we want dO*dIdO=chromatic_shift/blurtolerance. This is the unconfused span. 

unconfused_dO=abs((dI*blurtolerance)./dIdO);

unconfused_fraction=unconfused_dO./Obj;

%% pick out a few ranges of interest
range_in_meters1=0.2;
range_in_meters2=1.0;
range_in_meters3=5.0;
range1=range_in_meters1/(lensradius*1E-3); % find range in units of lens radius, which is in mm
range2=range_in_meters2/(lensradius*1E-3);
range3=range_in_meters3/(lensradius*1E-3);

% use min function to find array elements that come closest to desired ones
[value, index1]=min(abs(Obj-range1));
[value,index2]=min(abs(Obj-range2));
[value,index3]=min(abs(Obj-range3));

% convert to units of meters
coarseObj=lensradius*1E-3*[Obj(index1) Obj(index2) Obj(index3)]; % this is range on which we're focused
coarse_dO=lensradius*1E-3*[unconfused_dO(index1) unconfused_dO(index2) unconfused_dO(index3)]; % this is dz value that induces blur at specified level


%% plots

% chromatic focus vs. range, 3 wavelengths

figure(1)
semilogx(Obj,I,'g','LineWidth',3)
axis([4 1000 1.5 5])
grid on
hold on
semilogx(Obj,Ir,'r','LineWidth',3)
semilogx(Obj,Ib,'b','LineWidth',3)
xlabel('O/R, log scale')
ylabel('I/R, linear scale')
%text(9.88,3,'No Real Image Formed','HorizontalAlignment','center',... 
%	'BackgroundColor',[.7 .9 .7],'FontSize',11)
hold off
shg

%% chromatic focus vs. range, linear plot with real distances
figure(10)
Obj_meters=Obj*lensradius*1E-3;
Image_mm=I*lensradius;
Image_r_mm=Ir*lensradius;
Image_b_mm=Ib*lensradius;

plot(Obj_meters,Image_r_mm,'r',Obj_meters,Image_mm,'g',...
    Obj_meters,Image_b_mm,'b','LineWidth',3)
set(gca,'FontSize',20)
axis([0.2 4 11.6 13])
legend('550 nm','500 nm','450 nm')
grid on

xlabel('Distance to Object, meters')
ylabel('Accommodation, millimeters')
%text(9.88,3,'No Real Image Formed','HorizontalAlignment','center',... 
%	'BackgroundColor',[.7 .9 .7],'FontSize',11)

shg
% 
% figure(2)
% semilogx(Obj,dI,'k','LineWidth',3)
% xlabel('O/R, log scale')
% ylabel('focus shift, in units of lens radius')
% grid on
% shg

%% plot of depth of field vs. range
figure(3)
semilogx(Obj, unconfused_dO,'k','LineWidth',3)
grid on
hold on
% here is unity line:
semilogx(Obj,Obj,'b','LineWidth',3)
axis([8 1000 0 500])
ylabel('depth of field in lens radius units')
shg
hold off


figure(4)
errorbar(log10(Obj),0*Obj, unconfused_dO,'k')
grid on
shg

figure(5)
plot(Obj,unconfused_dO,'b')
grid on
shg

figure(6)
semilogx(Obj,dIdO,'k','LineWidth',3)
grid on
ylabel('dI/dO')
shg


figure(7)
errorbar(log10(Obj),0*Obj, unconfused_fraction,'k')
grid on
title('chromatically unconfused fraction of range to object')
shg

figure(8)
errorbar([1 2 3], coarseObj, coarse_dO,'k.')
ylabel('meters')
title(['Range and depth of field for designated fraction of chromatic image degradation'])
% grid on

