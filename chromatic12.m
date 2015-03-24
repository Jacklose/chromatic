% chromatic12.m, simulates chromatic aberration for Cephalopods
% Released under GPL v2 license by Christopher Stubbs 

% version 12, March 22 2015

%%%%%%%%%%%%% 
% This program computes chromatic aberrations of an image, given
% these inputs:
%  1) pupil structure, as a function of (r,theta).
%  2) solar photon spectrum
%  3) array of reflectance spectra in R,G,B for test pattern image
%  4) opsin response function
%  5) lens diameter and focal length, for chromatic blur computation
%  6) depth-dependent attentuation of sunlight in seawater.
%  7) depth in ocean for this computation, for attentuation calc.
%  8) blue-yellow intensity ratio, called colorbalance.
%  9) measured fraction change in focal length vs. wavelength. 
% 10) vector MTFspacing that contains best-focus wavelengths to be used

% User can specify the pupil mask, and two colors used in the test
% pattern.

% This version loops through focus settings and creates blurred images at
% each focus setting. Uses partial-annular and on-axis circular pupils.

% output is in a few different forms. The MTF values computed at the
% accommodation settings appear in MTF.dat, along with color information as RGB triplets
% and the pupil outer radius (to distinguish pupil masks).

% MTF uses the (max-min)/average definition, so the biggest possible value is 2. 
 
% Another potentially interesting output is a panel of images at the first 
% three accommodation settings. It's useful to pick these to be 450,500 and
% 550 nm in order to span the opsin peak. Note that although the program
% stores a .JPG file it does not look the same as the figure window. I've
% found it more effective to save the Figure as a PDF file, by hand. 

% other spectral plots of interest can also be generated from the data. 

% Finally, the relevant results are stored in a .mat file for future import
% and plotting. Filename for this should be changed each time the program
% is run since it is overwritten. 

% All output files should appear in the current working directory, which
% should be added to the MATLAB path. 

% the program was not written for maximal computational efficiency! It
% takes about 7 minutes to run through three accommodation settings on a
% high end laptop. Be patient. 

%%%%%%%%%%%%%% the fine print %%%%%%%%%%%%%
% The reflectance spectra should have 1 nm spacing, from 350 nm to 900 nm. 
% in the example in the program some real-world data are splined onto this
% wavelength spacing. 

% The other pieces of input information (or computed here if not available)
% are the opsin photon response times the illimination spectrum, again
% spanning 350 to 900 nm, and the pupil shape. Three pupil options are provided, but 
% of course one could always add more if desired. 

% spatial arrays are in focal plane coordinates, with units of 5 microns per
% pixel. This amounts to one pixel per rabdome. 

% In general, spectra are all column vectors, with one row per wavelength.
% note good reference for various reflection spectra is http://speclib.jpl.nasa.gov/search-1

% Image arithmetic in MATLAB:
% 2-d arrays of numbers can be mapped into any color range for displaying. 
% for this, need a 2-d array of numbers, and a color map.

% Matlab handles RGB images as MxNx3 arrays, with RGB values that range from 0 to
% 1. 

% Two-dimensional matrices can be *displayed* with some color map, in
% two ways: 1) image(): direct mapping of intensity onto colormap (cmap) array, or
% 2) imagesc(): scaled so that max and min intensity values are scaled onto color map.
% 
% we can convert from the 3-plane image representation and 2-d intensity
% arrays with rgb2gray() and mat2gray(A,[amin amax]). The first of these takes an RGB
% images and produces an intensity image with entries between 0 and 1. The
% latter takes a 2-d array and map it onto 0 to 1, with optional values amin and amax. 

% for data type double, intensity is scaled from 0 to 1. 
% for data type uint8 and unit16, intensity values range 0 to 255 or 64K. duh.

% we can make rgb images with ind2rgb(Array,map) where map is the colormap
% but that's tough to do with the custom color map used below. It was
% easier to capture the image from the display with getframe()


%% intitialize. Edit these parameters as needed. 
clear all
subplot(1,1,1)

depth = 3; % water depth in meters
lambda=350:900; % array of discrete wavelengths, one nm spacing. 
lambda=lambda'; % make it a column vector
MTFspacing=450:20:650; % array of wavelengths at which images are made and MTF is computed. Stay within 450 and 650

focallength=12; % focal length in mm. We assume f/1.2 beam at full aperture and so a 10 mm diameter lens

diagnosticplots=false; % set to true for verbose plotting

% full pupil:
% % % 
% pupilouterradius=4 % outer radius of pupil, in mm
% pupilinnerradius=0 % inner radius of pupil annulus, in mm
% pupilangulargap=0 % pupil gap full angluar extent in degrees
% % % % 
% % % small pupil, 2 mm diameter
% pupilouterradius=1;
% pupilinnerradius=0;
% pupilangulargap=0;
% % % % %%%

% this is a small pupil area of A=pi*r^2 = pi square mm.

% % % annular pupil:
%%% make dr appropriate to match area of small pupil. Area of annulus is 
%%% A = 2*pi*[(360-angulargap)/360]*rinner*dr which we want to equal pi
%% in order to match small pupil above. So dr=[360/2*(360-gap)*rinner]

%%% uncomment these lines for annular pupil
pupilinnerradius=3.0; % inner radius of pupil annulus, in mm
pupilangulargap=180; % pupil gap full angluar extent in degrees
dr=(360/(2*(360-pupilangulargap)*pupilinnerradius));  % chosen to match 
pupilouterradius=pupilinnerradius+dr; % outer radius of pupil, in mm
%%% 

testrows=200; % sets size of test pattern
testcols=1000;
chirp_factor=25000; % this determines the bar pattern frequency change. Default is 25000. smaller numbers give higher freq 

% Now pick the two colors for the test pattern. Color 1 is left-most, then it 
% transitions to color2. These are RGB values as floating point numbers
% between 0 and 1. Scale each of them to have intensity of one.
red=[1 0 0]; 
green=[0 1 0];
blue=[0 0 1];
yellow=[1 1 0]; % didn't normalize this since need full intensity of green
black=[0 0 0];
white=[1 1 1];
grey=0.95*[1 1 1];

colorbalance=1.05; % this sets the ratio of the blue to the yellow signal. 

 % pick colors here. Color 1 is  left-most bar
 % also designate a file save name (overwrites existing) for this image
 % trio
 
color1=black;
color2=blue; % be sure to change file names below as appropriate
MTFlabel='Annular Pupil Black/Blue';
savename='BlkBlu_annular.mat';
figname='BlkBlu_annular.jpg'; % name of high quality JPG file that is made

% wipe out existing files so we over-write them. 
delete(figname);
delete(savename);
delete(['MTF' filename]);
delete(['accomMTF' filename]);

MTFfile='MTF.dat'; % keep this one, gets appended to

%% initialization for publication quality figures
close all;
lw=3;
msz=3;
width=5;
height=5;

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultAxesLineWidth',3);

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
close all % wipe out figure frame

%% produce solar spectrum
% read in typical solar photon spectrum, from 
% http://www.pvlighthouse.com.au/resources/optics/spectrum%20library/spectrum%20library.aspx
% there is a function that reads in the relevant data file, which has 1 nm
% spacing from 350 to 900 nm. Note this is a photon spectrum. This was
% converted to a relative photon spectrum by multiplying by lambda.

filename='solar_photons.dat';
solarspectrum1=readsolarspectrum(filename);
% convert to column vector, to conform to our standard, and normalize it
solarspectrum1=(solarspectrum1)/max(solarspectrum1);

% plot solar spectrum
% figure(1)
% plot(lambda, solarspectrum,'b.')
% grid on
% title('solar photon spectrum')
% shg

% now account for attenuation in seawater. Data from Morel & Maritonrena(2001)
% units of attenuation are absorption coeff in m^-1

attenuationlengthdata=[
301.703	0.098
315.619	0.081
322.819	0.069
334.373	0.055
353.546	0.039
369.399	0.031
393.589	0.024
422.558	0.024
451.842	0.024
484.446	0.023
496.430	0.029
507.439	0.034
516.631	0.047
521.839	0.049
531.880	0.053
548.697	0.064
567.216	0.080
576.892	0.093
585.546	0.117
594.737	0.175
601.937	0.215
611.091	0.264
628.022	0.316
654.574	0.361
671.875	0.400
687.842	0.457
704.718	0.525
719.347	0.580
720.000 0.580
730.000 0.580
800     0.580
900     0.580
];

z0=spline(attenuationlengthdata(:,1),attenuationlengthdata(:,2),lambda);

% compute the transmission using depth in m
T=exp(-depth.*z0);
% make attenuated solar spectrum
solarspectrum=T.*solarspectrum1;

% normalize to unity at peak
solarspectrum=(solarspectrum)/max(solarspectrum);

% show attenutation of solar spectrum
if (diagnosticplots)
    figure()
    plotyy(lambda,solarspectrum1,lambda,solarspectrum)
end
%% empirical opsin response function, from ALS email Feb 22 2015
opsinsamples=[
300.487	0.133
315.327	0.188
325.523	0.226
334.229	0.252
341.907	0.262
347.582	0.265
363.513	0.256
375.775	0.232
383.964	0.220
398.437	0.211
406.168	0.218
414.548	0.246
423.776	0.305
431.511	0.356
438.551	0.424
445.771	0.506
451.777	0.580
460.755	0.687
469.524	0.801
480.237	0.907
488.942	0.964
491.589	0.989
496.104	1.001
502.110	1.001
511.790	0.969
516.369	0.936
520.954	0.889
525.597	0.825
530.628	0.747
538.120	0.629
543.732	0.535
550.830	0.428
558.961	0.314
566.250	0.219
576.504	0.127
583.149	0.082
591.860	0.048
602.115	0.023
612.886	0.009
623.663	0.006
645.600	4.058e-4
673.398	4.058e-4
694.941	5.999e-4
720.551	0.000
752.418	0.000
776.543	0.000
798.149	0.000
];
% spline this onto our default wavelength spacing of 1 nm.
opsinresponse=spline(opsinsamples(:,1), opsinsamples(:,2),lambda);

if (diagnosticplots)
    figure()
    plot(lambda,opsinresponse,'rv')
    grid on
    hold on
    plot(lambda,solarspectrum,'kv')
    title('opsin response vs. wavelength in nm')
    shg
end

%% compute focus shifts relative to 550 nm, using octopus data from Jaegger and Sands paper
% note this is an extrapolation for lambda < 450 nm and >700nm
focusshift=-5.4676e-05*lambda.^2 + 0.0794*lambda -27.1047; % percent focus shift, this is an array
focusshiftmm=1E-2*focusshift*focallength; % convert to mm of focus shift
% 
% figure(30)
% plot(lambda,focusshiftmm,'k.')
% grid on
% shg

% compute accommodation settings relative to 550 nm for the specfic
% best-focus wavelengths we are investigating. 

accommodation=1E-2*(-5.4676e-05*MTFspacing.^2 + 0.0794*MTFspacing -27.1047)*focallength; % this is an array

%% spectral reflectance functions, for BGR channels

fish=[
    0.3248    0.0431    0.0781
    0.3344    0.0435    0.0820
    0.3442    0.0439    0.0857
    0.3541    0.0443    0.0893
    0.3641    0.0447    0.0928
    0.3742    0.0451    0.0961
    0.3844    0.0456    0.0994
    0.3947    0.0460    0.1025
    0.4050    0.0465    0.1055
    0.4154    0.0470    0.1084
    0.4258    0.0476    0.1111
    0.4362    0.0481    0.1138
    0.4467    0.0487    0.1163
    0.4571    0.0492    0.1187
    0.4674    0.0498    0.1211
    0.4778    0.0504    0.1233
    0.4880    0.0511    0.1254
    0.4983    0.0517    0.1274
    0.5084    0.0523    0.1293
    0.5184    0.0530    0.1310
    0.5283    0.0537    0.1327
    0.5381    0.0543    0.1343
    0.5477    0.0550    0.1358
    0.5573    0.0557    0.1372
    0.5667    0.0564    0.1385
    0.5759    0.0571    0.1397
    0.5851    0.0579    0.1408
    0.5943    0.0586    0.1418
    0.6033    0.0593    0.1428
    0.6123    0.0601    0.1436
    0.6213    0.0608    0.1444
    0.6302    0.0616    0.1450
    0.6391    0.0623    0.1456
    0.6481    0.0631    0.1461
    0.6570    0.0638    0.1465
    0.6660    0.0646    0.1469
    0.6750    0.0654    0.1471
    0.6841    0.0661    0.1473
    0.6933    0.0669    0.1474
    0.7025    0.0677    0.1474
    0.7118    0.0684    0.1474
    0.7213    0.0692    0.1473
    0.7308    0.0699    0.1471
    0.7405    0.0707    0.1468
    0.7504    0.0714    0.1465
    0.7604    0.0722    0.1461
    0.7706    0.0729    0.1457
    0.7809    0.0737    0.1451
    0.7913    0.0744    0.1446
    0.8018    0.0751    0.1439
    0.8121    0.0758    0.1432
    0.8223    0.0765    0.1425
    0.8324    0.0772    0.1416
    0.8421    0.0779    0.1408
    0.8515    0.0786    0.1398
    0.8605    0.0793    0.1389
    0.8690    0.0799    0.1378
    0.8769    0.0806    0.1368
    0.8842    0.0812    0.1356
    0.8909    0.0819    0.1345
    0.8967    0.0825    0.1332
    0.9018    0.0831    0.1320
    0.9061    0.0836    0.1307
    0.9097    0.0842    0.1293
    0.9127    0.0848    0.1279
    0.9150    0.0853    0.1265
    0.9167    0.0858    0.1250
    0.9179    0.0863    0.1235
    0.9185    0.0868    0.1220
    0.9186    0.0873    0.1204
    0.9183    0.0877    0.1188
    0.9176    0.0881    0.1172
    0.9165    0.0885    0.1156
    0.9151    0.0889    0.1139
    0.9134    0.0893    0.1122
    0.9114    0.0896    0.1104
    0.9092    0.0899    0.1087
    0.9068    0.0902    0.1069
    0.9043    0.0905    0.1051
    0.9016    0.0907    0.1033
    0.8989    0.0909    0.1014
    0.8962    0.0911    0.0996
    0.8934    0.0913    0.0977
    0.8907    0.0914    0.0958
    0.8881    0.0915    0.0939
    0.8855    0.0916    0.0920
    0.8832    0.0916    0.0901
    0.8810    0.0916    0.0882
    0.8790    0.0916    0.0863
    0.8773    0.0916    0.0843
    0.8759    0.0915    0.0824
    0.8747    0.0914    0.0805
    0.8738    0.0912    0.0785
    0.8732    0.0910    0.0766
    0.8728    0.0908    0.0747
    0.8727    0.0906    0.0727
    0.8728    0.0903    0.0708
    0.8732    0.0901    0.0689
    0.8739    0.0898    0.0670
    0.8748    0.0895    0.0651
    0.8760    0.0893    0.0632
    0.8774    0.0891    0.0613
    0.8791    0.0889    0.0594
    0.8810    0.0887    0.0575
    0.8832    0.0887    0.0557
    0.8857    0.0886    0.0539
    0.8884    0.0887    0.0521
    0.8914    0.0888    0.0503
    0.8946    0.0890    0.0485
    0.8981    0.0893    0.0468
    0.9018    0.0897    0.0451
    0.9058    0.0902    0.0434
    0.9100    0.0909    0.0417
    0.9145    0.0916    0.0401
    0.9193    0.0926    0.0385
    0.9242    0.0936    0.0369
    0.9295    0.0948    0.0354
    0.9350    0.0962    0.0339
    0.9407    0.0978    0.0324
    0.9467    0.0995    0.0309
    0.9529    0.1015    0.0295
    0.9588    0.1036    0.0282
    0.9633    0.1060    0.0269
    0.9659    0.1086    0.0256
    0.9668    0.1114    0.0244
    0.9660    0.1144    0.0232
    0.9639    0.1177    0.0220
    0.9603    0.1212    0.0210
    0.9556    0.1251    0.0199
    0.9499    0.1291    0.0189
    0.9432    0.1335    0.0180
    0.9357    0.1382    0.0171
    0.9277    0.1431    0.0163
    0.9191    0.1484    0.0155
    0.9101    0.1539    0.0148
    0.9009    0.1598    0.0141
    0.8917    0.1661    0.0135
    0.8824    0.1727    0.0130
    0.8734    0.1796    0.0125
    0.8647    0.1869    0.0121
    0.8565    0.1946    0.0118
    0.8489    0.2027    0.0115
    0.8420    0.2114    0.0113
    0.8359    0.2206    0.0112
    0.8305    0.2304    0.0111
    0.8256    0.2409    0.0111
    0.8210    0.2521    0.0112
    0.8167    0.2641    0.0113
    0.8124    0.2770    0.0116
    0.8081    0.2907    0.0119
    0.8034    0.3054    0.0123
    0.7984    0.3211    0.0128
    0.7928    0.3378    0.0133
    0.7864    0.3557    0.0140
    0.7792    0.3747    0.0147
    0.7709    0.3949    0.0155
    0.7615    0.4165    0.0164
    0.7506    0.4393    0.0174
    0.7383    0.4636    0.0185
    0.7243    0.4889    0.0197
    0.7085    0.5149    0.0210
    0.6907    0.5413    0.0223
    0.6708    0.5677    0.0238
    0.6486    0.5936    0.0254
    0.6240    0.6188    0.0270
    0.5971    0.6427    0.0288
    0.5689    0.6651    0.0307
    0.5404    0.6860    0.0326
    0.5123    0.7056    0.0347
    0.4857    0.7238    0.0369
    0.4611    0.7408    0.0392
    0.4385    0.7567    0.0416
    0.4178    0.7715    0.0441
    0.3990    0.7854    0.0467
    0.3820    0.7983    0.0494
    0.3667    0.8105    0.0523
    0.3530    0.8220    0.0553
    0.3409    0.8328    0.0583
    0.3303    0.8430    0.0615
    0.3211    0.8527    0.0649
    0.3132    0.8618    0.0683
    0.3066    0.8703    0.0718
    0.3012    0.8784    0.0755
    0.2970    0.8859    0.0792
    0.2937    0.8930    0.0831
    0.2915    0.8997    0.0871
    0.2901    0.9060    0.0912
    0.2896    0.9118    0.0954
    0.2898    0.9173    0.0997
    0.2907    0.9225    0.1041
    0.2921    0.9273    0.1087
    0.2941    0.9318    0.1133
    0.2966    0.9361    0.1180
    0.2994    0.9401    0.1229
    0.3025    0.9438    0.1279
    0.3059    0.9474    0.1329
    0.3094    0.9508    0.1381
    0.3130    0.9540    0.1434
    0.3166    0.9571    0.1487
    0.3201    0.9600    0.1542
    0.3235    0.9627    0.1598
    0.3267    0.9653    0.1655
    0.3296    0.9678    0.1713
    0.3321    0.9701    0.1772
    0.3342    0.9723    0.1832
    0.3357    0.9743    0.1892
    0.3367    0.9762    0.1954
    0.3370    0.9779    0.2017
    0.3366    0.9796    0.2081
    0.3355    0.9811    0.2147
    0.3337    0.9824    0.2215
    0.3313    0.9837    0.2285
    0.3283    0.9848    0.2359
    0.3247    0.9859    0.2436
    0.3206    0.9868    0.2518
    0.3160    0.9876    0.2604
    0.3109    0.9883    0.2695
    0.3053    0.9889    0.2792
    0.2994    0.9894    0.2896
    0.2931    0.9898    0.3006
    0.2864    0.9901    0.3123
    0.2794    0.9903    0.3248
    0.2722    0.9904    0.3382
    0.2646    0.9904    0.3524
    0.2569    0.9904    0.3675
    0.2490    0.9903    0.3837
    0.2410    0.9901    0.4006
    0.2328    0.9898    0.4183
    0.2245    0.9895    0.4364
    0.2162    0.9891    0.4549
    0.2079    0.9886    0.4736
    0.1995    0.9881    0.4923
    0.1912    0.9875    0.5108
    0.1830    0.9868    0.5291
    0.1749    0.9861    0.5471
    0.1669    0.9854    0.5650
    0.1591    0.9846    0.5827
    0.1515    0.9838    0.6003
    0.1442    0.9829    0.6176
    0.1371    0.9820    0.6349
    0.1303    0.9810    0.6519
    0.1238    0.9801    0.6688
    0.1177    0.9791    0.6856
    0.1119    0.9780    0.7021
    0.1065    0.9770    0.7183
    0.1014    0.9759    0.7344
    0.0965    0.9748    0.7501
    0.0920    0.9737    0.7656
    0.0877    0.9726    0.7807
    0.0837    0.9715    0.7956
    0.0800    0.9704    0.8100
    0.0765    0.9693    0.8241
    0.0732    0.9682    0.8379
    0.0701    0.9671    0.8511
    0.0673    0.9660    0.8640
    0.0647    0.9649    0.8763
    0.0622    0.9638    0.8881
    0.0600    0.9627    0.8994
    0.0580    0.9617    0.9102
    0.0561    0.9607    0.9203
    0.0544    0.9597    0.9298
    0.0529    0.9587    0.9387
    0.0515    0.9578    0.9469
    0.0503    0.9569    0.9543
    0.0493    0.9561    0.9611
    0.0484    0.9552    0.9673
    0.0476    0.9544    0.9728
    0.0470    0.9537    0.9777
    0.0465    0.9530    0.9820
    0.0461    0.9523    0.9858
    0.0458    0.9516    0.9890
    0.0456    0.9510    0.9917
    0.0455    0.9504    0.9940
    0.0455    0.9498    0.9958
    0.0456    0.9493    0.9972
    0.0458    0.9487    0.9982
    0.0460    0.9483    0.9989
    0.0464    0.9478    0.9992
    0.0467    0.9474    0.9992
    0.0472    0.9470    0.9989
    0.0476    0.9466    0.9984
    0.0481    0.9462    0.9976
    0.0487    0.9459    0.9967
    0.0493    0.9456    0.9955
    0.0499    0.9453    0.9942
    0.0505    0.9451    0.9928
    0.0511    0.9449    0.9913
    0.0517    0.9447    0.9897
    0.0523    0.9445    0.9881
    0.0530    0.9443    0.9865
    0.0535    0.9442    0.9849
    0.0541    0.9441    0.9833
    0.0547    0.9440    0.9817
    0.0552    0.9439    0.9802
    0.0556    0.9438    0.9787
    0.0560    0.9438    0.9773
    0.0564    0.9438    0.9759
    0.0568    0.9438    0.9746
    0.0571    0.9438    0.9732
    0.0574    0.9438    0.9720
    0.0576    0.9439    0.9707
    0.0578    0.9440    0.9695
    0.0580    0.9440    0.9683
    0.0581    0.9441    0.9672
    0.0583    0.9443    0.9661
    0.0583    0.9444    0.9651
    0.0584    0.9445    0.9641
    0.0584    0.9447    0.9631
    0.0584    0.9449    0.9622
    0.0584    0.9451    0.9613
    0.0584    0.9453    0.9604
    0.0583    0.9455    0.9596
    0.0582    0.9457    0.9589
    0.0581    0.9459    0.9581
    0.0580    0.9462    0.9574
    0.0578    0.9464    0.9568
    0.0576    0.9467    0.9562
    0.0575    0.9469    0.9556
    0.0573    0.9472    0.9550
    0.0570    0.9475    0.9545
    0.0568    0.9478    0.9541
    0.0566    0.9481    0.9536
    0.0563    0.9484    0.9532
    0.0560    0.9487    0.9528
    0.0558    0.9491    0.9525
    0.0555    0.9494    0.9522
    0.0552    0.9497    0.9519
    0.0549    0.9500    0.9516
    0.0546    0.9504    0.9514
    0.0543    0.9507    0.9511
    0.0540    0.9511    0.9509
    0.0536    0.9514    0.9508
    0.0533    0.9518    0.9506
    0.0530    0.9521    0.9505
    0.0527    0.9524    0.9504
    0.0524    0.9528    0.9502
    0.0520    0.9531    0.9502
    0.0517    0.9535    0.9501
    0.0514    0.9538    0.9500
    0.0511    0.9542    0.9500
    0.0508    0.9545    0.9500
    0.0505    0.9549    0.9499
    0.0502    0.9552    0.9499
    0.0500    0.9556    0.9499
    0.0497    0.9559    0.9499
    0.0494    0.9562    0.9499
    0.0492    0.9565    0.9499
    0.0490    0.9568    0.9499
    0.0488    0.9572    0.9499
    0.0486    0.9575    0.9500
    0.0484    0.9578    0.9500
    0.0482    0.9581    0.9500
    0.0480    0.9583    0.9500
    0.0479    0.9586    0.9500
    0.0478    0.9589    0.9500
    0.0476    0.9591    0.9501
    0.0475    0.9594    0.9501
    0.0474    0.9596    0.9501
    0.0474    0.9598    0.9501
    0.0473    0.9601    0.9501
    0.0472    0.9603    0.9501
    0.0472    0.9605    0.9501
    0.0471    0.9606    0.9501
    0.0471    0.9608    0.9501
    0.0471    0.9610    0.9501
    0.0471    0.9611    0.9501
    0.0471    0.9612    0.9501
    0.0471    0.9613    0.9501
    0.0471    0.9614    0.9501
    0.0472    0.9615    0.9501
    0.0472    0.9616    0.9500
    0.0473    0.9616    0.9500
    0.0473    0.9616    0.9500
    0.0474    0.9617    0.9500
    0.0475    0.9616    0.9500
    0.0475    0.9616    0.9500
    0.0476    0.9616    0.9500
    0.0477    0.9615    0.9500
    0.0478    0.9614    0.9500
    0.0480    0.9613    0.9500
    0.0481    0.9612    0.9500
    0.0482    0.9610    0.9500
    0.0483    0.9609    0.9500
    0.0485    0.9607    0.9500
    0.0486    0.9605    0.9499
    0.0488    0.9602    0.9499
    0.0490    0.9600    0.9499
    0.0491    0.9597    0.9499
    0.0493    0.9594    0.9499
    0.0495    0.9590    0.9499
    0.0497    0.9587    0.9499
    0.0498    0.9583    0.9499
    0.0500    0.9579    0.9499
    0.0502    0.9574    0.9499
    0.0504    0.9569    0.9499
    0.0506    0.9564    0.9499
    0.0508    0.9559    0.9499
    0.0510    0.9554    0.9500
    0.0513    0.9548    0.9500
    0.0515    0.9542    0.9500
    0.0517    0.9535    0.9500
    0.0519    0.9528    0.9500
    0.0521    0.9521    0.9500
    0.0524    0.9514    0.9500
    0.0526    0.9506    0.9501
    0.0528    0.9498    0.9501
    0.0530    0.9490    0.9501
    0.0533    0.9481    0.9501
    0.0535    0.9472    0.9502
    0.0537    0.9462    0.9502
    0.0540    0.9453    0.9502
    0.0542    0.9442    0.9503
    0.0544    0.9432    0.9503
    0.0547    0.9421    0.9504
    0.0549    0.9410    0.9504
    0.0551    0.9398    0.9505
    0.0554    0.9386    0.9505
    0.0556    0.9374    0.9506
    0.0558    0.9361    0.9506
    0.0560    0.9348    0.9507
    0.0563    0.9334    0.9508
    0.0565    0.9320    0.9508
    0.0567    0.9306    0.9509
    0.0569    0.9291    0.9510
    0.0571    0.9276    0.9511
    0.0574    0.9260    0.9511
    0.0576    0.9244    0.9512
    0.0578    0.9228    0.9513
    0.0580    0.9211    0.9514
    0.0582    0.9193    0.9515
    0.0584    0.9175    0.9516
    0.0586    0.9157    0.9517
    0.0587    0.9138    0.9518
    0.0589    0.9119    0.9519
    0.0591    0.9099    0.9521
    0.0593    0.9079    0.9522
    0.0594    0.9059    0.9523
    0.0596    0.9037    0.9524
    0.0598    0.9016    0.9526
    0.0599    0.8994    0.9527
    0.0601    0.8971    0.9529
    0.0602    0.8948    0.9530
    0.0603    0.8924    0.9532
    0.0604    0.8900    0.9533
    0.0606    0.8876    0.9535
    0.0607    0.8851    0.9537
    0.0608    0.8825    0.9539
    0.0609    0.8799    0.9541
    0.0610    0.8772    0.9542
    0.0610    0.8745    0.9544
    0.0611    0.8717    0.9546
    0.0612    0.8689    0.9548
    0.0612    0.8660    0.9551
    0.0613    0.8630    0.9553
    0.0613    0.8600    0.9555
    0.0613    0.8570    0.9557
    0.0613    0.8538    0.9560
    0.0613    0.8507    0.9562
    0.0613    0.8474    0.9565
    0.0613    0.8441    0.9567
    0.0613    0.8408    0.9570
    0.0612    0.8374    0.9572
    0.0612    0.8339    0.9575
    0.0611    0.8304    0.9578
    0.0611    0.8268    0.9581
    0.0610    0.8231    0.9584
    0.0609    0.8194    0.9587
    0.0608    0.8156    0.9590
    0.0607    0.8118    0.9593
    0.0605    0.8079    0.9596
    0.0604    0.8039    0.9599
    0.0602    0.7999    0.9603
    0.0601    0.7958    0.9606
    0.0599    0.7916    0.9610
    0.0597    0.7874    0.9613
    0.0595    0.7831    0.9617
    0.0592    0.7788    0.9621
    0.0590    0.7743    0.9625
    0.0587    0.7698    0.9628
    0.0585    0.7653    0.9632
    0.0582    0.7606    0.9637
    0.0579    0.7559    0.9641
    0.0576    0.7512    0.9645
    0.0572    0.7463    0.9649
    0.0569    0.7414    0.9654
    0.0565    0.7364    0.9658
    0.0562    0.7314    0.9663
    0.0558    0.7262    0.9667
    0.0553    0.7210    0.9672
    0.0549    0.7158    0.9677
    0.0545    0.7104    0.9682
    0.0540    0.7050    0.9687
    0.0535    0.6995    0.9692
    0.0530    0.6940    0.9697
    0.0525    0.6883    0.9702
    0.0520    0.6826    0.9707
    0.0514    0.6768    0.9713
    0.0508    0.6709    0.9718
    0.0502    0.6650    0.9724
    0.0496    0.6589    0.9730
    0.0490    0.6528    0.9735
    0.0483    0.6467    0.9741
    0.0477    0.6404    0.9747
    0.0470    0.6341    0.9753
    0.0463    0.6276    0.9760
    0.0455    0.6211    0.9766
    0.0448    0.6146    0.9772
    0.0440    0.6079    0.9779
    0.0432    0.6011    0.9785
    0.0424    0.5943    0.9792
    0.0416    0.5874    0.9799
    0.0407    0.5804    0.9806
    0.0398    0.5733    0.9813
    0.0389    0.5662    0.9820
    0.0380    0.5589    0.9827
    0.0370    0.5516    0.9834
    0.0361    0.5442    0.9842
    0.0351    0.5367    0.9849
    0.0341    0.5291    0.9857
    0.0330    0.5214    0.9865
    0.0320    0.5136    0.9872
    0.0309    0.5058    0.9880
    0.0297    0.4978    0.9888
    0.0286    0.4898    0.9897
    0.0274    0.4817    0.9905
    0.0263    0.4735    0.9913
    0.0250    0.4652    0.9922
    0.0238    0.4568    0.9931
    0.0225    0.4483    0.9939
    0.0212    0.4397    0.9948
    0.0199    0.4310    0.9957
    0.0186    0.4223    0.9966
    0.0172    0.4134    0.9975
    0.0158    0.4045    0.9985
    0.0144    0.3954    0.9994
    0.0129    0.3863    1.0004
    0.0115    0.3770    1.0014
    0.0100    0.3677    1.0023
    0.0084    0.3583    1.0033
    0.0069    0.3488    1.0044
    0.0053    0.3392    1.0054
    0.0037    0.3294    1.0064
    0.0020    0.3196    1.0075
    0.0003    0.3097    1.0085
    0.0014    0.2997    1.0096
    0.0031    0.2896    1.0107
    0.0049    0.2794    1.0118
    0.0067    0.2691    1.0129
    0.0085    0.2586    1.0140
    0.0103    0.2481    1.0151
    0.0122    0.2375    1.0163];

% write these into individual vectors, with some sensible normalization
blue_reflectance=fish(:,1); % intensity fixed to get similar integrals in blue and yellow
green_reflectance=fish(:,2);
red_reflectance=fish(:,3);

bluesignal=blue_reflectance.*opsinresponse.*solarspectrum;
yellowsignal=(green_reflectance+red_reflectance).*opsinresponse.*solarspectrum;
ratio=sum(bluesignal)/sum(yellowsignal);
blue_reflectance=colorbalance*blue_reflectance/ratio; % adjust to get desired imbalance between blue and yellow


%% create MTF test pattern. Size is determined by testrows and testcols, set in initialization block at top. 
% this needs to be large enough to accommodate the convolutions later on
% we will create an array that has ones for color1 and zeros to designate color2
testimage=ones(testrows,testcols); % 200 x 1000 initial array
xpattern=1:testcols;

% one dimensional vector with test pattern transition points
% this is a periodic square wave with changing period that depends on x
barpattern=0.5+0.5*square(2*pi*(xpattern.^2/chirp_factor)); % makes 1's and 0's for bar pattern

% now extend this onto the 2-d array as 1's and 0's by stepping through
% columns
for rowcounter=1:testrows
    testimage(rowcounter,:)=testimage(rowcounter,:).*barpattern;
end


%% Create RGB color panes for the test pattern image, using the color1 and
% color2 choices made above. The superposition is made by taking the
% testimage mask and its complement, designated with "~"

% these are two dimensional arrays

redimage=color1(1)*testimage+color2(1)*~testimage;
greenimage=color1(2)*testimage+color2(2)*~testimage;
blueimage=color1(3)*testimage+color2(3)*~testimage;

testRGB=zeros(testrows,testcols,3);

for rowcounter=1:testrows
    for colcounter=1:testcols
        testRGB(rowcounter,colcounter,1)=redimage(rowcounter,colcounter);
        testRGB(rowcounter,colcounter,2)=greenimage(rowcounter,colcounter);
        testRGB(rowcounter,colcounter,3)=blueimage(rowcounter,colcounter);
    end
end

%%
figure()
imshow(testRGB)
shg
% make green one same as red one, for chromatic blur calculation below
% greenpixellogical=redpixellogical;
% greencheckpixels=redcheckpixels;



% % now capture the image that was displayed and write it to a png file
% frame=getframe()
% imwrite(frame.cdata,'testimage.png')
% 

%% make intensity image and find pixels with RGB content

BWimage=redimage+greenimage+blueimage;

if (diagnosticplots)
    figure()
    imshow(redimage)
    title('red channel of original image')
    figure()
    imshow(greenimage)
    title('green channel of original image')
    figure()
    imshow(blueimage)
    title('blue channel of original image')
    shg
end

%% now let's make a stack of images, one at each of a set of discrete wavelengths

% initialize some arrays
redpart=zeros(testrows,testcols,length(lambda));
greenpart=redpart;
bluepart=redpart;
flux=redpart;

% this next object is an array of color images, each of which is NxMx3 for
% RGB
rgb_by_lambda=zeros(testrows,testcols,3,length(lambda));

% now run through each pixel at each wavelength and add up flux
% contributions to make an unblurred intensity ar
for rowcounter=1:testrows
    rowcounter
    for colcounter=1:testcols
        for wavelengthcounter=1:length(lambda)
           redpart(rowcounter,colcounter,wavelengthcounter)=redimage(rowcounter,colcounter).*red_reflectance(wavelengthcounter);
           greenpart(rowcounter,colcounter,wavelengthcounter)=greenimage(rowcounter,colcounter).*green_reflectance(wavelengthcounter);
           bluepart(rowcounter,colcounter,wavelengthcounter)=blueimage(rowcounter,colcounter).*blue_reflectance(wavelengthcounter);
           flux(rowcounter,colcounter,wavelengthcounter)=(redpart(rowcounter,colcounter,wavelengthcounter)+bluepart(rowcounter,colcounter,wavelengthcounter)+...
               greenpart(rowcounter,colcounter,wavelengthcounter))...
               .*solarspectrum(wavelengthcounter).*opsinresponse(wavelengthcounter); % modify according to the illumination   
           rgb_by_lambda(rowcounter,colcounter,1,wavelengthcounter)=redpart(rowcounter,colcounter,wavelengthcounter).*solarspectrum(wavelengthcounter).*opsinresponse(wavelengthcounter);
           rgb_by_lambda(rowcounter,colcounter,2,wavelengthcounter)=greenpart(rowcounter,colcounter,wavelengthcounter).*solarspectrum(wavelengthcounter).*opsinresponse(wavelengthcounter);
           rgb_by_lambda(rowcounter,colcounter,3,wavelengthcounter)=bluepart(rowcounter,colcounter,wavelengthcounter).*solarspectrum(wavelengthcounter).*opsinresponse(wavelengthcounter);
        end
    end
end
imagestack=flux; % just renaming this for consistency with what follows.
intensity=sum(flux,3); % sum over wavelengths, check to be sure test pattern is flat 
intensity=0.99*intensity./(max(max(intensity))); % renormalize for plots
%% line cut across the summed image
if (diagnosticplots)
    figure()
    plot(intensity(100,:))
    shg
end

%% scroll through the image stack, wavelength by wavelength
if (diagnosticplots)
    for foo=1:300
    imshow(rgb_by_lambda(:,:,:,foo))
    shg
    lambda(foo)
    pause(0.05)
    end
end


%% set up pupil mask array. At the plane of the lens, the pupil has a size of 
% 10 mm x 10 mm. As it gets closer to the best-focus it shrinks in size to 
% S=10*df/12 mm, in proportion to the defocus distance df. 
% we can use the interp2(Xorig,Yorig,Z,Xnew,Ynew,'cubic') function to rebin
% the pupil onto the closest 5 micron pixel grid so that convolutions 
% will work OK. We are never more than 0.8mm = 800 microns out of focus so
% the blur kernel is at most 660 microns ~ 130 pixels in extent. 

% Note that the pupil goes through a parity flip as we pass through focus.

% need to be sure we don't close off the annuluar pupil with pixels that
% are too coarse. If we assume minimum pupil opening is 1 mm then 
% we need to be at least 5 microns * 12mm/1mm = 60 microns away to have the
% opening be one pixel in extent. 

pupilmask=zeros(2000,2000);  % with a 5 micron spacing this is 10 mm x 10mm
[pupilX, pupilY]=meshgrid(1:2000,1:2000); % make 2-d coordinate arrays

% convert inner and outer pupil radii to pixel units
r1=pupilinnerradius/5E-3; % pupil units were mm, now in 5 micron pixels at lens plane
r2=pupilouterradius/5E-3;
for maskrow=1:2000
    maskrow % print this out to monitor progress
    for maskcol=1:2000
        masky=maskrow-2000/2; 
        maskx=maskcol-2000/2;
        rpix=sqrt((maskx).^2+(masky).^2);  % distance from center of pupil
        thetapix=atan2d(masky,maskx)+90 ; % rotate so theta=0 is up, angles increase CCW. this is a number between -180 and 180
        if((rpix<=r2)&&(rpix>=r1)&&(abs(thetapix)>=pupilangulargap/2))
            pupilmask(maskrow,maskcol)=1;
        end % end of conditional mask criteria
    end % end of column loop
end % end of mask row loop
%%
if (diagnosticplots)
    figure()
    imshow(pupilmask)
    truesize([250 250])
    shg
end

%% set up an outermost loop over different focus settings

focuscounter=0;
varyfocusstack2=zeros(length(MTFspacing),testrows,testcols); % initialize array. First argument here must be number of focus settings in loop below
blurredimage=zeros(testrows,testcols,300);

for bestfocuslambda=min(MTFspacing):(MTFspacing(end)-MTFspacing(end-1)):max(MTFspacing) % Be sure to go to previous line and initialize the focus stack to number of settings. 
                                 %change lines 1058, 1165
                                 % 1118, 1137 and 1247 and 1291
    
focuscounter=focuscounter+1; % increment this counter, it's initialized so that it's 1 the first time through

%% compute offset from best-focus wavelength. 
% Find  the index that corresponds to best focus at this wavelength
bestindex=floor(bestfocuslambda)-349; % this is an array

% now compute magnitude of focus offset from that best-focus wavelength. 
% sign convention is with positive focusshift being further from lens
% this quantity is negative for wavelengths that focus closer to the lens
% than the best-focus rays 
relativefocusshift=focusshiftmm-focusshiftmm(bestindex); % this is a vector with a value for each wavelength
    
%% now create a stack of blurred images of what the creature would see, and display them, by looping through
% wavelengths, blurring appropriately, and accounting for the opsin
% response. 

for counter=1:300 % this steps through wavelengths, 350 to 650 nm
    
      % now compute pupil screen size in pixels
      
        pupilarraysize=ceil((abs(relativefocusshift(counter))/12)*2000); % rounds upwards to integer pixel value
        
        % we need to re-grid the pupil onto this smaller array of 5 micron pixels
        stride=ceil(2000/pupilarraysize);
        xregrid=1:stride:2000; % new array should be pupilarraysize 
        [newY, newX]=meshgrid(xregrid,xregrid); % it's a square so this is ok

        % interpolate onto new pupil array. For some reason need to rotate
        % it too
        newpupil=interp2(pupilX,pupilY,pupilmask,newX,newY,'cubic');
        newpupil=rot90(newpupil);
        % flip it if this wavelength focuses closer than bestfocus wavelength 
        if(relativefocusshift(counter)>0)
            newpupil=flipud(newpupil); % flip up-down for other side of focus
            newpupil=fliplr(newpupil);
        end
        
        extent=size(newpupil);
        
        if(extent(1)<=3) % is it a really tiny pupil plane? 
            newpupil=zeros(3);   % if so, hardwire to 3x3 kernel with single central pixel=1
            newpupil(2,2)=1;
        end
        
        midpoint=floor(extent(1)/2); % find midpoint of new pupil array
        
        % now normalize to sum of values in this convolution kernel
        integral=sum(sum(newpupil));
        % trap on case where pupil hole has vanished, put it back
            if (integral<1)
                newpupil(midpoint,midpoint)=1;
                integral=1;
            end
            
        blurkernel=newpupil;
        blurkernel=newpupil./integral;  % normalize to a total of one
        sum(sum(blurkernel)); % check it
        
        if(diagnosticplots)
                figure()
                imshow(blurkernel)
                shg
        end
    
    % fuzz out the image using this kernel, produce image of same size and
    % replicate edges as needed
    % this makes a 3-d data structure where x and y are pixels and z is
    % wavelength
    %
    blurredimage(:,:,counter)=imfilter(imagestack(:,:,counter),blurkernel,'same','symmetric');
    % look at it
    if (diagnosticplots)
        figure()
        imshow(blurredimage(:,:,counter))
        shg
        title(lambda(counter))
        pause(0.2)
    end
    counter % spit out loop counter for diagnostic purposes
end % end of wavelength loop

%% Let's look at the final resulting image, added up over all wavelengths
figure()
finalimage=sum(blurredimage,3);
finalimage=finalimage/max(max(finalimage)); % renormalize it to display properly
imshow(finalimage)
labelstring=num2str(bestfocuslambda);
label2=strcat(labelstring,'nm focus, chromatically blurred image, summed over opsin response ');
title(label2)

%% add this image to the stack of varying-focus images
varyfocusstack2(focuscounter,:,:)=finalimage;

end
%%  display the stack of focus-scan images
% show one image to get window setting correct
% also collapse this to a simple 2-d array. The squeeze() function 
% gets rid of pesky additional dimensions that are only 1 deep.

for showcounter=1:length(MTFspacing)  % set this equal to number of focus settings
subplot(1,1,1)
imshow(squeeze(varyfocusstack2(showcounter,:,:)))
shg
pause(1)
% go ahead and compute MTF's while we're in this loop
MTFmax(showcounter)=max(max(squeeze(varyfocusstack2(showcounter,100,750:950))));
MTFmin(showcounter)=min(min(squeeze(varyfocusstack2(showcounter,100,750:950))));
end

MTF=2*(MTFmax-MTFmin)./(MTFmax+MTFmin); %peak-to-peak difference divided by the average.  
% This is a vector with dimension equal to number of focus settings 

%% plot cross section as a line cut
[a, b, c]=size(varyfocusstack2);
midframe=ceil(a/2);

x0=1;
y0=1;
figure()
set(gca,'Units','normalized', 'FontUnits','points', 'FontWeight','normal', 'FontSize',9, 'FontName','Times New Roman')
axes('FontName','Times New Roman') 
    subplot(3,1,1)
    plotline1=squeeze(varyfocusstack2(1,100,:));
    plot(plotline1,'b')
 %   title('focus at 400 nm')
    axis([1 950 0 1])   % trim off the edge that has artifacts
    grid on
    subplot(3,1,2)
    plotline2=squeeze(varyfocusstack2(midframe,100,:));
    plot(plotline2,'k')
    axis([1 950 0 1])
    grid on
  %  title('focus at 450 nm')
    subplot(3,1,3)
    plotline3=squeeze(varyfocusstack2(end,100,:));
    plot(plotline3,'r')
    axis([1 950 0 1])
   % title('focus at 500')
    grid on
    shg

%% load in data set
% load('BlkBlu_full.mat');
%% panel plot. Plots the first, last and central focus settings
% positions are [left bottom width height]

h1=figure();

first=midframe-2;
second=midframe;
third=midframe+2;

position1=[-0.02, 0.8, 1, 0.1]; % top panel -  color image
position2=[-0.02, 0.65, 1, 0.1]; % image 1
position3=[0.1, 0.55, 0.8, 0.1]; % linecut 1
position4=[-0.02, 0.40, 1, 0.1]; % image 2
position5=[0.1, 0.30, 0.8, 0.1]; % linecut 2
position6=[-0.02, 0.15, 1, 0.1]; % image 3
position7=[0.1, 0.05, 0.8, 0.1]; % linecut 3


% trim off right edge to suppress convolution artifacts
imageOrig=testRGB(:,1:950,:);  % color image, no blurring
imagef1=squeeze(varyfocusstack2(first,:,1:950)); % focus 450 nm
imagef2=squeeze(varyfocusstack2(second,:,1:950)); % focus 500 nm
imagef3=squeeze(varyfocusstack2(third,:,1:950)); % focus 550 nm

linecutf1=squeeze(varyfocusstack2(first,100,1:950));
linecutf2=squeeze(varyfocusstack2(second,100,1:950));
linecutf3=squeeze(varyfocusstack2(third,100,1:950));

xscale=1:950; % convert to units of microns for x axis
xscale=xscale/5;

subplot('Position',position1); 
imshow(imageOrig); 
% title(pupiltype,'FontName','Times','FontSize',30)

subplot('Position',position2);
imshow(imagef1,'InitialMagnification','fit','Border','tight')

subplot('Position',position3)
plot(linecutf1,'b'); 
set(gca,'xtick',[],'ytick',[])
%text(20,0.4,'\bf 450 nm','FontName','Times','FontSize',16)
%text(20,0.2,' \bf focus','FontName','Times','FontSize',16)
ylim([0 1.1])

subplot('Position',position4);
imshow(imagef2,'Border','tight')

subplot('Position',position5)
plot(linecutf2,'k'); 
set(gca,'xtick',[],'ytick',[])
%text(20,0.4,'\bf 500 nm','FontName','Times','FontSize',16)
%text(20,0.2,' \bf focus','FontName','Times','FontSize',16)
ylim([0 1.1])

subplot('Position',position6);
imshow(imagef3,'Border','tight')

% since we don't show x axes above this is the only one with units of
% microns
subplot('Position',position7)
plot(xscale,linecutf3,'r'); 
set(gca,'xtick',[],'ytick',[])
%text(4,0.4,'\bf 550 nm','FontName','Times','FontSize',16)
%text(4,0.2,' \bf focus','FontName','Times','FontSize',16)
ylim([0 1.1])

% subplot('Position',position8)
% imshow(zeros(10))

truesize([200,1000])
shg

frame=getframe(h1);
% % now capture the image that was displayed and write it to a JPG file
imwrite(frame.cdata,figname,'Quality',100) % tried 'Mode','lossless' but it fails for some reason. 

%print(figname,'-djpeg','-cmyk')

%% 
% spectral plots
if (diagnosticplots)
    figure()
    plot(lambda,blue_reflectance,'b.')
    grid on
    hold on
    plot(lambda,green_reflectance,'g.')
    plot(lambda,red_reflectance,'r.')
    plot(lambda,opsinresponse.*solarspectrum,'kv')
    shg
    hold off

    bluesignal=blue_reflectance.*opsinresponse.*solarspectrum;
    yellowsignal=(green_reflectance+red_reflectance).*opsinresponse.*solarspectrum;
    ratio=sum(bluesignal)/sum(yellowsignal); % normalize to the same intensity
    yellowsignal=yellowsignal*1.05*ratio; % introduce 10 percent difference

    figure(16)
    plot(lambda,bluesignal,'b')
    hold on
    plot(lambda,yellowsignal,'r')
    grid on
    hold off
    xlabel('wavelength, nm','FontName','Helvetica','FontSize',16)
    ylabel('incident photon irradiance','FontName','Helvetica','FontSize',16)
    shg
end

%% plot MTF vs accommodation
%% MTF plot vs accommodation
h3=figure();
plot(accommodation,MTF,'k*')
hold on
plot(accommodation,MTF,'b')
xlabel('accommodation, mm','FontName','HELVETICA','FontSize',16)
ylabel('MTF','FontName','Helvetica','FontSize',16)
shg
axis([-0.35 0.2 0 2.1]);
text(-0.2,2,MTFlabel,'FontName','HELVETICA','FontSize',20) % for some reason this does not appear in saved JPEG!
hold off

% this gets changed a few
accomMTFfigname=['accomMTF' figname];

frame=getframe(h3);
% % now capture the image that was displayed and write it to a pdf file

imwrite(frame.cdata,accomMTFfigname,'Quality',100)


%% store the resulting data structures

save(savename,'varyfocusstack2','blurredimage','MTF')

% append MTF to the data file

outvec=[MTF color1 color2 pupilouterradius];
dlmwrite(MTFfile,outvec,'-append');

%% some spectral plots

subplot(5,1,1)
plot(lambda, solarspectrum)
text(650,0.8,'solar illumination','FontName','Helvetica','FontSize',16)
ylim([0 1.1])
title(MTFlabel,'FontName','Helvetica','FontSize',16)

shg

subplot(5,1,2)
plot(lambda,opsinresponse)
text(650,0.8,'opsin response','FontName','Helvetica','FontSize',16)
ylim([0 1.1])
shg

subplot(5,1,3)
plot(lambda,blue_reflectance,'b')
ylim([0 1.1])
hold on
plot(lambda,green_reflectance,'g')
plot(lambda,red_reflectance,'r')
shg
ylim([0 1.1])
hold off
text(650,0.8,'reflectances','FontName','Helvetica','FontSize',16)

% these next ones are solar*opsin*reflectance
subplot(5,1,4)
plot(lambda,bluesignal,'b')
hold on
plot(lambda,yellowsignal,'g')
text(650,0.8,'photon signal','FontName','Helvetica','FontSize',16)
if (diagnosticplots)
    plot(MTFspacing,MTF,'k*')
end
ylim([0 1.1])
shg
hold off

subplot(5,1,5)
plot(accommodation,MTF,'k*')
hold on
plot(accommodation,MTF,'k')
shg
ylim([0 1.1])
text(-0.2,1.5,'MTF vs. accommodation, mm','FontName','Helvetica','FontSize',16)
axis([-0.7 0.8 0 2.1])
hold off


%% now loop through them all, with a pause in between
% focussetting=450:50:550;
% while(1)
% for focusdispcounter=1:3
%     labelstring=num2str(focussetting(focusdispcounter));
%     dispframe=squeeze(varyfocusstack2(focusdispcounter,:,:));
%     imshow(dispframe)
%     title(labelstring)
%     shg
%     pause(0.2)
% end
% end
% 
% 





    

