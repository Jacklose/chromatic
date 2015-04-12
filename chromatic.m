% chromatic.m, simulates chromatic aberration for Cephalopods
% Released under GPL v2 license by Christopher Stubbs and Alexander Stubbs
% astubbs@berkeley.edu 

% version 19, April 7, 2015. Includes figure plots for paper

%%%%%%%%%%%%% 
% This program computes chromatic aberrations of an image, given
% these inputs:
%  1) pupil structure, as a function of (r,theta).
%  2) solar photon spectrum
%  3) array of three reflectance spectra for test pattern image
%  4) opsin response function vs. wavelength
%  5) lens diameter and focal length, for chromatic blur computation
%  6) depth-dependent attentuation of sunlight in seawater.
%  7) depth in ocean for this computation, for attentuation calc.
%  8) blue-yellow intensity ratio, called colorbalance, for blue/yellow checkerboard.
%  9) measured fractional change in focal length vs. wavelength. 
% 10) vector MTFspacing that contains best-focus wavelengths to be used

% User can specify the pupil mask, and two colors used in the test
% pattern.

% This version loops through focus settings and creates blurred images at
% each focus setting. Uses partial-annular and on-axis circular pupils.

% There is an interactive stage where the threshold for MTF computation can
% be adjusted. Best to have sound turned on to hear the "beep" at the
% prompt for input!

% output is in a few different forms. The MTF values computed at the
% accommodation settings are appended to MTF.dat, along with color information as RGB triplets
% and the pupil outer radius (to distinguish pupil masks). There is also a
% run-specific spectral data file that contains MTF and input spectrum.

% MTF is image quality factor. We determine this by looking for pixels that exceed some
% multiplicative factor above the mean. 
 
% Another potentially interesting output is a panel of images at the first 
% three accommodation settings. It's useful to pick these to be 450,500 and
% 550 nm in order to span the opsin peak. Note that although the program
% stores a .JPG file it does not look the same as the figure window.

% other spectral plots of interest can also be generated from the data. 

% Finally, the relevant results are stored in a .mat file for future import
% and plotting. Filename for this should be changed each time the program
% is run since it is overwritten. 

% All output files should appear in the current working directory, which
% should be added to the MATLAB path. 

% the program was not written for maximal computational efficiency! It
% takes about 15 minutes to run through three accommodation settings on a
% high end laptop. Be patient. 

% Another feature is the ability to create point spread function images, by
% setting the PSFimage parameter to 'true'. This produces a test image with
% a single non-zero pixel, at pixel(100,500), the center of the image.

% There are a number of other logical switches that can be used to turn
% things on and off. These include
% diagnosticplots - if set true a number of diagnostic plots are produced. 
% PSFimage - if set true, the test image is replaced with single point
% GaussianSpectrum- if set true, get Gaussian narrowband spectrum

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

% Image arithmetic in MATLAB:
% 2-d arrays of numbers can be mapped into any color range for displaying. 
% for this, need a 2-d array of numbers, and a color map.

% Matlab handles RGB images as MxNx3 arrays, with RGB values that range from 0 to
% 1. M is row variable and N is column

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
close all

%%%%%%%%%%%%%%% these parameters don't change very frequently %%%%%%%%
depth = 3; % water depth in meters. default is 3
lambda=350:900; % array of discrete wavelengths, one nm spacing. default is 350:900 and this is essential!
lambda=lambda'; % make it a column vector
MTFspacing=450:5:650; % array of wavelengths at which images are made and MTF is computed. Stay within 450 and 650. default is 450:5:650.
PSFimage=false % make this false for normal operation. Set to true if test image of single point is desired
focallength=12; % focal length in mm. We assume f/1.2 beam at full aperture and so a 10 mm diameter lens. default is 12
lambdastride=5; % spacing between wavelengths for which the image stack is convolved and summed. Default is 5 nm for broadband sources
diagnosticplots=false; % set to true for verbose plotting
GaussianSpectrum=false; % creates spiky spectrum, default is false
if (GaussianSpectrum)
    lambdastride=1;
    MTFspacing=450:1:650
end
testrows=200; % sets size of test pattern. default is 200 and 1000
testcols=1000;

% this next parameter determines the bar pattern frequency change.  Smaller numbers give higher freq
% this is a balance between getting a wide dynamic range in frequency but numbers below 10,000 have major 
% aliasing for annular pupil. 
chirp_factor=13000;  % 13000 seems a good compromise choice for this number. 

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

colorbalance=1.05; % this sets the ratio of the blue to the yellow signal. Set to 99 for no correction
MTFcut=1.5; % threshold for fractional MTF computation. Reduce this number if there is a sharp spectral feature
            % superimposed on other light. This is the multiplicative
            % factor above the mean light level that we consider to be
            % "resolved" bar target. For example for dual Gaussian spectra
            % of equal intensity, the MTF will never exceed 1.5. For broad
            % continuum reflectance spectra, a number like 1.8 seems a good
            % choice. Making this too low lets aliasing artifacts seep
            % through. Here are some good values:
            % single Gaussian line in spectrum: 1.7
            % broad continuum reflectance spectra: 1.9
            % dual Gaussian sources: 1.4. 
            % NOTE- be sure to reduce lambdastride above to 1 and
            % MTFspacing interval to 2 if using narrowband sources. 

%%%% Frequently edited parameters begin here: pupil and color choices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin pupil selection section %%%%%%%%%%%
% comment out all but one of these choices. 

% pupil='Small';
% pupilswitch=1;

% pupil='Full';
% pupilswitch=2;

pupil='Annular'
pupilswitch=3;

% pupil='ReallySmall';
% pupilswitch=4;

%%%%%%%%%%%%% end pupil selection section %%%%%%%%%%%%

 % pick colors here. Color 1 is  left-most bar
 % also designate a file save name (overwrites existing) for this image
 % trio. Use color1=black as default
 
 %%
color1=blue;
color2=yellow; % be sure to change file name below as appropriate. color2 is used for spectral plots
MTFlabel=['BluYel' pupil]; % no spaces, no slashes! Pupil type is appended automatically. Use color info. 


% this does some housekeeping if GaussianSpectrum is set. Leave it alone
if (GaussianSpectrum) 
    color1=black
    color2=blue
end

%%
%%%%%%%%% end of section where choices are edited. %%%%%%%%%%%%

% housekeeping for files and names
%%
savename=[MTFlabel '.mat']; % name of workspace file that is saved
figname=[MTFlabel '.jpg']; % name of high quality JPG file that is made
datafilename=[MTFlabel '.dat']; % name of spectral data file
%%
% wipe out existing files so we over-write them. 
delete(figname);
delete(savename);
delete(datafilename);
 % wipe out previous image versions too
delete(['AccomMTF' figname]); 
delete(['SpecAccomMTF' figname]);

MTFfile='MTF.dat'; % keep this one, gets appended to

warning('off','all') % suppress warnings to command line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% define pupil parameters %%%%%
if (pupilswitch==1)
    % % small pupil, 1 mm diameter
    pupilouterradius=0.5;
    pupilinnerradius=0;
    pupilangulargap=0;
end

if (pupilswitch==2)
% full pupil, 8 m diameter
    pupilouterradius=4 % outer radius of pupil, in mm
    pupilinnerradius=0 % inner radius of pupil annulus, in mm
    pupilangulargap=0 % pupil gap full angluar extent in degrees
end

if (pupilswitch==3)
    % we make this same as small pupil area of A=pi*r^2 = pi square mm.

    % % % annular pupil:
    %%% make dr appropriate to match area of small pupil. Area of annulus is 
    %%% A = 2*pi*[(360-angulargap)/360]*rinner*dr which we want to equal pi
    %%% in order to match small pupil above. So dr=[360/2*(360-gap)*rinner]
    %%% hardwire dr to 0.33
        pupilinnerradius=3.0; % inner radius of pupil annulus, in mm
        pupilangulargap=180; % pupil gap full angluar extent in degrees
        dr=0.3333;  % chosen to match 
        pupilouterradius=pupilinnerradius+dr; % outer radius of pupil, in mm
    %%% 
end

if (pupilswitch==4)
   % really small 
    pupilouterradius=0.2;
    pupilinnerradius=0;
    pupilangulargap=0;
end

%%%%% end of definition of pupil parameters %%%%%%%%%%%

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

%% produce solar photon spectrum
% read in typical solar photon spectrum, from 
% http://www.pvlighthouse.com.au/resources/optics/spectrum%20library/spectrum%20library.aspx
% Array below hhas 1 nm
% spacing from 350 to 900 nm. Note this is a photon spectrum. This was
% converted to a relative photon spectrum by multiplying by lambda.
solarspectrum1= [
0.1839
0.1976
0.2051
0.2076
0.2063
0.2025
0.1975
0.1926
0.1889
0.1877
0.1904
0.1977
0.2085
0.2214
0.2347
0.2470
0.2570
0.2642
0.2684
0.2693
0.2667
0.2607
0.2532
0.2465
0.2427
0.2442
0.2521
0.2638
0.2756
0.2837
0.2846
0.2759
0.2613
0.2457
0.2343
0.2321
0.2422
0.2604
0.2804
0.2962
0.3015
0.2925
0.2747
0.2562
0.2448
0.2486
0.2728
0.3123
0.3592
0.4058
0.4442
0.4686
0.4807
0.4842
0.4829
0.4804
0.4798
0.4814
0.4849
0.4897
0.4957
0.5024
0.5091
0.5152
0.5202
0.5233
0.5242
0.5237
0.5231
0.5233
0.5257
0.5306
0.5363
0.5404
0.5403
0.5337
0.5192
0.5006
0.4825
0.4698
0.4673
0.4781
0.4990
0.5249
0.5509
0.5719
0.5845
0.5901
0.5917
0.5922
0.5945
0.6009
0.6112
0.6243
0.6393
0.6552
0.6712
0.6863
0.6998
0.7106
0.7181
0.7217
0.7224
0.7216
0.7207
0.7209
0.7234
0.7276
0.7326
0.7376
0.7417
0.7443
0.7455
0.7459
0.7461
0.7464
0.7473
0.7487
0.7506
0.7527
0.7549
0.7572
0.7600
0.7635
0.7683
0.7748
0.7829
0.7913
0.7982
0.8018
0.8003
0.7928
0.7812
0.7681
0.7562
0.7483
0.7464
0.7498
0.7570
0.7666
0.7774
0.7879
0.7974
0.8050
0.8102
0.8120
0.8102
0.8058
0.8003
0.7951
0.7917
0.7911
0.7931
0.7972
0.8027
0.8090
0.8156
0.8215
0.8257
0.8274
0.8255
0.8196
0.8109
0.8013
0.7926
0.7866
0.7846
0.7859
0.7894
0.7940
0.7986
0.8023
0.8056
0.8092
0.8139
0.8204
0.8291
0.8390
0.8488
0.8570
0.8623
0.8637
0.8619
0.8580
0.8531
0.8482
0.8444
0.8419
0.8409
0.8415
0.8438
0.8480
0.8533
0.8590
0.8646
0.8692
0.8722
0.8741
0.8752
0.8761
0.8771
0.8787
0.8804
0.8818
0.8824
0.8817
0.8794
0.8762
0.8729
0.8703
0.8693
0.8704
0.8731
0.8763
0.8793
0.8811
0.8812
0.8799
0.8782
0.8768
0.8765
0.8780
0.8807
0.8839
0.8869
0.8890
0.8897
0.8899
0.8905
0.8926
0.8973
0.9050
0.9139
0.9217
0.9258
0.9239
0.9145
0.9002
0.8844
0.8705
0.8620
0.8614
0.8670
0.8762
0.8865
0.8952
0.9004
0.9027
0.9034
0.9038
0.9052
0.9086
0.9134
0.9190
0.9245
0.9291
0.9322
0.9337
0.9338
0.9325
0.9298
0.9261
0.9220
0.9186
0.9166
0.9170
0.9204
0.9256
0.9312
0.9358
0.9380
0.9367
0.9325
0.9266
0.9200
0.9138
0.9088
0.9055
0.9041
0.9049
0.9080
0.9136
0.9209
0.9291
0.9371
0.9442
0.9496
0.9534
0.9559
0.9575
0.9583
0.9586
0.9586
0.9581
0.9571
0.9557
0.9538
0.9511
0.9471
0.9416
0.9342
0.9248
0.9148
0.9061
0.9003
0.8992
0.9040
0.9133
0.9250
0.9372
0.9479
0.9554
0.9603
0.9635
0.9658
0.9681
0.9712
0.9748
0.9783
0.9812
0.9831
0.9835
0.9829
0.9819
0.9812
0.9813
0.9826
0.9844
0.9855
0.9849
0.9815
0.9746
0.9642
0.9507
0.9343
0.9156
0.8951
0.8750
0.8582
0.8470
0.8443
0.8515
0.8660
0.8844
0.9030
0.9181
0.9272
0.9311
0.9318
0.9313
0.9313
0.9333
0.9371
0.9420
0.9471
0.9518
0.9555
0.9581
0.9599
0.9610
0.9615
0.9614
0.9591
0.9529
0.9409
0.9216
0.8940
0.8618
0.8297
0.8021
0.7839
0.7782
0.7822
0.7919
0.8032
0.8119
0.8152
0.8149
0.8141
0.8159
0.8236
0.8390
0.8597
0.8823
0.9033
0.9191
0.9274
0.9297
0.9289
0.9275
0.9285
0.9335
0.9414
0.9499
0.9569
0.9602
0.9586
0.9543
0.9507
0.9509
0.9581
0.9737
0.9906
1.0000
0.9930
0.9606
0.8975
0.8125
0.7179
0.6261
0.5494
0.4979
0.4728
0.4731
0.4976
0.5454
0.6140
0.6948
0.7781
0.8540
0.9126
0.9468
0.9604
0.9598
0.9514
0.9418
0.9362
0.9346
0.9357
0.9384
0.9414
0.9437
0.9449
0.9447
0.9429
0.9393
0.9339
0.9271
0.9197
0.9124
0.9060
0.9010
0.8975
0.8953
0.8943
0.8944
0.8955
0.8971
0.8987
0.8998
0.9001
0.8990
0.8974
0.8957
0.8949
0.8957
0.8982
0.9005
0.9002
0.8948
0.8819
0.8601
0.8325
0.8034
0.7769
0.7573
0.7476
0.7463
0.7505
0.7575
0.7644
0.7690
0.7713
0.7719
0.7713
0.7701
0.7688
0.7683
0.7692
0.7725
0.7788
0.7887
0.8011
0.8148
0.8285
0.8408
0.8508
0.8586
0.8646
0.8692
0.8728
0.8756
0.8776
0.8785
0.8781
0.8762
0.8726
0.8676
0.8613
0.8542
0.8464
0.8385
0.8316
0.8269
0.8260
0.8299
0.8394
0.8523
0.8655
0.8760
0.8811
0.8785
0.8705
0.8599
0.8497
0.8429
0.8416
0.8448
0.8505
0.8568
0.8618
0.8640
0.8637
0.8618
0.8591
0.8565
0.8546
0.8532
0.8522
0.8512
0.8499
0.8482
0.8465
0.8453
0.8452
0.8468
0.8502
0.8540
0.8562
0.8549
0.8483
0.8350
0.8158
0.7921
0.7651
0.7361
0.7066
0.6778
0.6511
0.6277
0.6090  
];

 solarspectrum1=(solarspectrum1)/max(solarspectrum1); % normalize it

  % plot solar spectrum
% figure(1)
% plot(lambda, solarspectrum,'b.')
% grid on
% title('solar photon spectrum')
% shg

% now account for optical attenuation in seawater. Data from Morel & Maritonrena(2001)
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

% convert this into our default spectral spacing of 350:900
z0=spline(attenuationlengthdata(:,1),attenuationlengthdata(:,2),lambda);

% compute the transmission using depth in meters. This is an array since z0
% has spacing of lambda
T=exp(-depth.*z0);
% make attenuated solar spectrum
solarspectrum=T.*solarspectrum1;

% normalize to unity at peak
solarspectrum=(solarspectrum)/max(solarspectrum);

% show initial and attenuted solar photon spectrum
if (diagnosticplots)
    figure()
    plotyy(lambda,solarspectrum1,lambda,solarspectrum)
end
%% empirical opsin response function, from ALS email Feb 22 2015
% reference is Queensland PhD thesis on visual systems
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
focusshiftmm=1E-2*focusshift*focallength; % convert to mm of focus shift. this is also an array

if (diagnosticplots)
    figure()
    plot(lambda,focusshiftmm,'k.')
    grid on
    shg
end
% compute accommodation settings in mm relative to 550 nm for the specfic
% best-focus wavelengths we are investigating. 

accommodation=1E-2*(-5.4676e-05*MTFspacing.^2 + 0.0794*MTFspacing -27.1047)*focallength; % this is an array

%% spectral reflectance functions, for BGR channels
% these are typical marine creature hyperspectral data
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
blue_reflectance=fish(:,1); 
green_reflectance=fish(:,2);
red_reflectance=fish(:,3);

%%% over-ride the reflectance with Gaussians if logic switch is set
if(GaussianSpectrum)
    % % % Multicolor Gaussians 
    % % Gaussians at two distinct wavelengths. index 151 is 500 nm. index 101 is 
    % % 450 nm and index 201 is 550 nm. Let's use Gaussian width of sigma=5 nm.

    blue_reflectance=0.001+0*fish(:,1);
    green_reflectance=0.001+0*fish(:,1);
    red_reflectance=0.001+0*fish(:,1);
    % use 500 and 555 nm.
    center1=151;
    center2=206;
    sigma=2;
    blue_reflectance=0.001+exp(-((lambda-lambda(center1)).^2)/(2*sigma^2))+...
        exp(-((lambda-lambda(center2)).^2)/(2*sigma^2));
    red_reflectance=0.001+0.*lambda;
    green_reflectance=red_reflectance;
    colorbalance=99; % this avoids rescaling of blue reflectance below since yellow has zero flux
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Gaussian spectrum section 


% % % % 
bluesignal=blue_reflectance.*opsinresponse.*solarspectrum;
greensignal=green_reflectance.*opsinresponse.*solarspectrum;
redsignal=red_reflectance.*opsinresponse.*solarspectrum;

%%%%% balance the fluxes to desired ratio in blue and yellow. 
yellowsignal=(green_reflectance+red_reflectance).*opsinresponse.*solarspectrum;
ratio=sum(bluesignal)/sum(yellowsignal);
if (colorbalance~=99) % conditional on not being set to 99, which defeats this
    blue_reflectance=colorbalance*blue_reflectance/ratio; % adjust to get desired imbalance between blue and yellow
    bluesignal=blue_reflectance.*opsinresponse.*solarspectrum;
end

% extract a data vector of input spectrum at MTF wavelength settings.
% Assumes color 1 is black and color 2 is of interest. format is RGB for
% color array indices
fullsignal=color2(1)*redsignal+color2(2)*greensignal+color2(3)*bluesignal;
spectrum=zeros(length(MTFspacing),1);
for counter=1:length(MTFspacing)
    [lambdafoo, lambdaindex]=min(abs(lambda-MTFspacing(counter))); % finds indices at MTFspacing values
    spectrum(counter)=fullsignal(lambdaindex); % assigns spectrum array the appropriate value for this MTFspacing index
    bluesignalMTFspacing(counter)=bluesignal(lambdaindex); % used later for MTF plots
    yellowsignalMTFspacing(counter)=yellowsignal(lambdaindex); % used later for MTF plots
end

spectrum=spectrum./max(spectrum); % normalize it

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

%% PSF image
%%%% single-point test image for PSF determination. skipped for normal operation. 
if (PSFimage)
    testimage=testimage.*0;
    testimage(100,500)=1;
end
%%%%%
%% Create RGB color panes for the test pattern image, using the color1 and
% color2 choices made above. The superposition is made by taking the
% testimage mask and its complement, designated with "~"

% these are two dimensional arrays. Gets the color content from both colors

redimage=color1(1)*testimage+color2(1)*~testimage;
greenimage=color1(2)*testimage+color2(2)*~testimage;
blueimage=color1(3)*testimage+color2(3)*~testimage;

if (PSFimage)
   redimage=testimage;
   greenimage=testimage;
   blueimage=testimage;
end

testRGB=zeros(testrows,testcols,3); % initialize RGB image

for rowcounter=1:testrows
    for colcounter=1:testcols
        testRGB(rowcounter,colcounter,1)=redimage(rowcounter,colcounter);
        testRGB(rowcounter,colcounter,2)=greenimage(rowcounter,colcounter);
        testRGB(rowcounter,colcounter,3)=blueimage(rowcounter,colcounter);
    end
end

%%
figure(1)
imshow(testRGB)
shg

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
% this is the unblurred image hypercube

% initialize some arrays
redpart=zeros(testrows,testcols,length(lambda));
greenpart=redpart;
bluepart=redpart;
flux=redpart;

% this next object is an array of color images, each of which is NxMx3 for
% RGB. So it's a hyperspectral cube of RGB images
rgb_by_lambda=zeros(testrows,testcols,3,length(lambda));

% now run through each pixel at each wavelength and add up flux
% contributions to make an unblurred intensity array
for rowcounter=1:testrows
    creating_test_images_at_row=rowcounter
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

%% scroll through the unblurred image stack, wavelength by wavelength
if (diagnosticplots)
    figure()
    for foo=1:300
    imshow(rgb_by_lambda(:,:,:,foo))
    shg
    lambda(foo)
    pause(0.05)
    end
end


%% set up pupil mask array. At the plane of the lens, the pupil has a size of 
% 10 mm x 10 mm so we start with a 2000 x 2000 array. As it gets closer to the best-focus it shrinks in size to 
% a side length of S=10*df/12 mm, in proportion to the defocus distance df. 
% we can use the interp2(Xorig,Yorig,Z,Xnew,Ynew,'cubic') function to rebin
% the pupil onto the closest 5 micron pixel grid so that convolutions 
% will work OK. We are never more than 0.8mm = 800 microns out of focus so
% the blur kernel is at most 660 microns ~ 130 pixels in extent. 

% Since we require an alignmnet between test image pixels and the pupil
% image, there is some granularity here that ends up generating step
% changes in modulation transfer function. But unless we do something much
% fancier with re-pixelization we'll just live with it. 

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
    making_pupil_mask_at_row=maskrow % print this out to monitor progress
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

%% set up an outermost loop over different focus (i.e. accommodation) settings. 

focuscounter=0;
got_one=false;

% This next array is the hyperspectral image cube after blurring at each
% wavelength, and then summed over wavelength. 
varyfocusstack2=zeros(length(MTFspacing),testrows,testcols); % initialize array. First argument here must be number of focus settings in loop below
blurredimage=zeros(testrows,testcols,300); % this is the temporary image stack made for each focus setting.  

for bestfocuslambda=min(MTFspacing):(MTFspacing(end)-MTFspacing(end-1)):max(MTFspacing) % this loops over focus settings in MTFspacing

focuscounter=focuscounter+1; % increment this counter, it's initialized so that it's 1 the first time through

%% compute offset from best-focus wavelength. 
% Compute the index in wavelength vector that is closest to this best-focus
% wavelength. The floor() function rounds to an integer
bestindex=floor(bestfocuslambda)-349; % this is the wavelength index that corresponds to best-focus lambda this time through loop. 350 nm is index = 1 

% now compute magnitude of focus offset from that of best-focus wavelength. 
% sign convention is with positive focusshift being further from lens
% this quantity is negative for wavelengths that focus closer to the lens
% than the best-focus rays. We want the units to be mm
relativefocusshift=focusshiftmm-focusshiftmm(bestindex); % this is a vector with a value for each wavelength, after subtracting off the accommodation setting
                                                         % setting for this
                                                         % wavelength
                                                          
    
%% now create a stack of blurred images of what the creature would see, and display them, by looping through
% wavelengths, blurring appropriately, and accounting for the opsin
% response and attenuated solar spectrum and reflectance 

  for counter=101:lambdastride:301 % this steps through wavelengths. Truncated range of wavelengths to conform to 
                     % Jagger and Sands range of focus shifts.  Use
                      % counter=101:lambdastride:301 for 450 to 650. For PSF images use
                      % 1:5:301. Be sure to reduce value of lambdastride
                      % for narrow spectral sources. For broadband spectra, lambdastride=5
                      % is OK. 

    
      % now compute pupil screen size in pixels for each wavelength
      
        pupilarraysize=ceil((abs(relativefocusshift(counter))/focallength)*2000); % rounds upwards to integer pixel value This is a scalar
        
        % we need to re-grid the pupil onto this smaller array of 5 micron pixels
        stride=ceil(2000/pupilarraysize); % remapping conversion factor, rounded to an integer
        xregrid=1:stride:2000; % new array should be pupilarraysize 
        [newY, newX]=meshgrid(xregrid,xregrid); % Make new pupil mask coords. it's a square so this is ok

        % interpolate onto new pupil array. For some reason need to rotate
        % it too
        newpupil=interp2(pupilX,pupilY,pupilmask,newX,newY,'cubic'); % this interpolates old pupil onto new grid
        newpupil=rot90(newpupil);
        % flip it in parity if this wavelength focuses closer than bestfocus wavelength 
        if(relativefocusshift(counter)>0)
            newpupil=flipud(newpupil); % flip up-down for other side of focus
            newpupil=fliplr(newpupil);
        end
        
        extent=size(newpupil); % look for cases where it's really become small
        
        if(extent(1)<=3) % is it a really tiny pupil plane? 
            newpupil=zeros(3);   % if so, hardwire to 3x3 kernel with single central pixel=1
            newpupil(2,2)=1;
        end
        
        midpoint=floor(extent(1)/2); % find midpoint of new pupil array
        
        % now normalize the sum of values in this convolution kernel
        integral=sum(sum(newpupil));
        % trap on case where pupil hole has vanished, put it back
            if (integral<1)
                newpupil(midpoint,midpoint)=1; % make central pixel = 1
                integral=1;
            end
            
        blurkernel=newpupil; % done to keep a convention used below. Legacy stuff, sorry. 
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
    % wavelength. This is the blurred hypercube for this focus setting
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
    % send a progress report to the command line
    status=['at wavelength ' num2str(lambda(counter)) ' for focus setting ' num2str(bestfocuslambda)] % spit out loop counters for diagnostic purposes
end % end of wavelength loop

%% Let's look at the final resulting image, added up over all wavelengths

figure(2)
finalimage=sum(blurredimage,3); % this adds up the blurred hypercube along the wavelength direction
finalimage=finalimage*(0.5*testrows*testcols/sum(sum(finalimage))); % they should all have same total number of photons. This enforces flux conservation
finalimagedisp=finalimage/max(max(finalimage)); % renormalize it to display properly. We don't use this for computations
imshow(finalimagedisp)
labelstring=num2str(bestfocuslambda);
label2=strcat(labelstring,' nm focus, chromatically blurred image ');
title(label2, 'FontName','HELVETICA','FontSize',16)
shg

%% add this image to the stack of varying-focus images
varyfocusstack2(focuscounter,:,:)=finalimage; % take this summed, blurred image and insert it into the image array vs. MTFsetting

end % end of focus setting loop

%%  display the stack of focus-scan images, and interactively set MTFcut

% show one image to get window setting correct
% also collapse this to a simple 2-d array. The squeeze() function 
% gets rid of pesky additional dimensions that are only 1 deep.

% the logic of this interactive loop is a bit weird but it works. 

newMTF=MTFcut; % initialize so we prompt for new value first time through

while(newMTF) % loops through iterations of choosing MTFcut value, as long as not zero
        figure(20)
        for showcounter=1:length(MTFspacing)  % set this equal to number of focus settings
                subplot(2,1,1) % upper plot is the image
                imshow(squeeze(varyfocusstack2(showcounter,:,:)))
                labelstring=num2str(MTFspacing(showcounter));
                label3=[MTFlabel, ' ', labelstring,' nm accommodation setting'];
                title(label3,'FontName','Helvetica','FontSize',16)
                truesize([200,1000])

                subplot(2,1,2) % lower plot is the line cut, with MTF threshold shown too
                linecut=squeeze(varyfocusstack2(showcounter,100,:));
                normlinecut=linecut/mean(linecut); % normalize to the mean value
                plot(1:testcols,normlinecut,'r')
                hold on
                plot(1:testcols,[MTFcut+0*testcols],'g'); % draws a line across plot
                hold off
                ylim([0 2.3])
                shg
                pause(0.2)

                %%% compute various image quality factors MTF1-MTFn. 

                % normalize to excursions above/below the mean flux level
                relativeFlux=squeeze(varyfocusstack2(showcounter,100,:));
                relativeFlux=relativeFlux/mean(relativeFlux); % normalize
                % 
                % plot(relativeFlux,'k.')
                % shg

                % MTF1 is number of pixels above some adjustable threshold
                % there are 500 potential array elements that might be above the mean
                MTF1(showcounter)=sum(relativeFlux > MTFcut)/500; % Adds up number of pixels that exceed the threshold, in the single-row line cut. 
                                                                 % this is an array with a value at each MTFspacing value
                % compute line transfer function, slope of transition
                % dark bar on the left has center at 42, and bright bar has center at column 98
                % so let's extract the section from 0 to 98. 

                % compute 10% to 90% dx
                linesegment=squeeze(varyfocusstack2(showcounter,100,42:98));
                linesegment=linesegment/max(linesegment);

                [foo1 index90]=min(abs(linesegment-0.9*(linesegment(end)-linesegment(1))));
                [foo2 index10]=min(abs(linesegment-0.1*(linesegment(end)-linesegment(1))));
                dx1090(showcounter)=index90-index10;
                xvec=index10:index90;
                xvec=xvec';
                yvec=linesegment(index10:index90);
                linfitparams=polyfit(xvec,yvec,1);
                slope(showcounter)=linfitparams(1); % this is slope

                % MTF3 is slope of first color1-color2 transition
                MTF3(showcounter)=slope(showcounter);

                % compute std deviation of pixel histogram
                % MTF4 is simply the standard deviation of the central line
                % cut, scaled up by a factor of two as an arbitrary
                % normalization
                MTF4(showcounter)=2*std(squeeze(varyfocusstack2(showcounter,100,:)))

       
        
                % compute peak intensity at widest bar
                MTF5(showcounter)=max(varyfocusstack2(showcounter,100,80:110));
        
         end % this is end of MTF showcounter loop
         
        % here is where we pick one MTF over the rest
        MTF=MTF4;
        
  %      MTF=MTF./max(MTF);  % renormalize to fractional MTF 
   


    %% MTF plot vs accommodation
    h3=figure(30);
    
    plot(accommodation,MTF,'b', accommodation, MTF, 'b*')
    hold on
   
    set(gca,'FontSize',20)
    set(gca,'box','off');  % turns off upper tick marks for wavelength that follows below
 %   plot(accommodation,spectrum,'k.')
    xlabel('Accommodation, mm (blue)','FontName','HELVETICA','FontSize',20)
    ylabel('MTF (blue), Irradiance (red)','FontName','Helvetica','FontSize',20)
    axis([-0.29 0.17 0 1.1]);
    % plottext=[MTFlabel ' with threshold ' num2str(MTFcut)]
    plottext=[MTFlabel]
   % text(-0.2,1.05,plottext,'FontName','HELVETICA','FontSize',20) % add
   % legend if desired. 
   
   
    
    % add wavelength axis on top and plot the spectra
    ax1_pos = get(gca,'Position'); % position of first axes
    ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
    set(gca,'FontSize',20);
    set(gca,'YTick',[]); % turn off right hand y axis

    axis([450 600 0 1.1]); 
    hold on
   

    if (((sum(color1==blue)==3) || (sum(color2==blue)==3)))
        fooblue=color2
        plot(MTFspacing,bluesignalMTFspacing,'r','Parent',ax2)
        plot(MTFspacing,bluesignalMTFspacing,'rv','Parent',ax2)
        got_one=true;
    end
    if (((sum(color1==yellow)==3) || (sum(color2==yellow)==3)))
        fooyell=color2
        plot(MTFspacing,yellowsignalMTFspacing,'r','Parent',ax2)
        plot(MTFspacing,yellowsignalMTFspacing,'rv','Parent',ax2)
        got_one=true;
    end
    if(got_one==false)
        plot(MTFspacing,spectrum,'r','Parent',ax2)
        plot(MTFspacing,spectrum,'rv','Parent',ax2)
    end
    text(500,1.07,'Wavelength, nm (red)','Parent',ax2,'FontSize',20)
  
    hold off
    hold off
    
    shg
    %%
    
    beep
    beep
    beep
    newMTF=input('New MTF value? (0 to keep this setting)')
    
    if (newMTF~=0)
        MTFcut=newMTF;
    end
    
    got_one=false;
  

end % end of while loop for MTFcut setting. ends if newMTF=0

%% histogram plot
if(diagnosticplots)
figure(25)
for wavecounter=1:length(MTFspacing)
hist(squeeze(varyfocusstack2(wavecounter,100,:)),30);
%    h = findobj(gcf,'type','patch');
%    set(h,'facecolor','','linewidth',2);
axis([0 1.3 0 250])
shg
titlestring=num2str(MTFspacing(wavecounter));
% set(gca,'FontSize',16)
title(titlestring)
pause(0.5)
end
end
%% now capture the image that was displayed and write it to a pdf file
frame=getframe(h3);
accomMTFfigname=['AccomMTF' figname];
imwrite(frame.cdata,accomMTFfigname,'Quality',100)



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

%% This is how you load in an existing data set, for future reference. 
% load('BlkBlu_full.mat');
%% panel plot. Plots the first, last and central focus settings
% positions are [left bottom width height]

h1=figure();

first=1;
second=midframe;
third=length(MTFspacing);

% these are all carefully tweaked to get alignment. Don't change!
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


%% store the resulting data structures

save(savename,'varyfocusstack2','blurredimage','MTF')

% make a file with spectrum and this MTF spectrum

outvec1=[MTFspacing' MTF' spectrum];
dlmwrite(datafilename,outvec1);

% append MTF to the data file

outvec=[MTF color1 color2 pupilouterradius];
dlmwrite(MTFfile,outvec,'-append');

%% some more spectral plots
h4=figure();

subplot(5,1,1)
plot(lambda, solarspectrum)
set(gca,'FontSize',16)
text(650,0.8,['solar illumination ',num2str(depth),' m'],'FontName','Helvetica','FontSize',16)
ylim([0 1.1])
title('Black/Yellow,Annular Pupil','FontName','Helvetica','FontSize',16)

shg

subplot(5,1,2)
plot(lambda,opsinresponse)
set(gca,'FontSize',16)
text(650,0.8,'opsin response','FontName','Helvetica','FontSize',16)
ylim([0 1.1])
shg

subplot(5,1,3)
plot(lambda,blue_reflectance,'b')
set(gca,'FontSize',16)
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
set(gca,'FontSize',16)
hold on
plot(lambda,yellowsignal,'g')
text(650,0.8,'photon signal','FontName','Helvetica','FontSize',16)
if (diagnosticplots)
    plot(MTFspacing,MTF4,'k*')
end
ylim([0 1.1])
shg
hold off

subplot(5,1,5)
plot(MTFspacing,MTF,'k*')
hold on
plot(MTFspacing,MTF,'k')
set(gca,'FontSize',16)
% plot(MTFspacing,spectrum,'r')
shg
ylim([0 1.1])
thislabel=['MTF threshold =', num2str(MTFcut)]; 
text(600,0.8,'MTF vs. accommodation','FontName','Helvetica','FontSize',16)
% text(650, 0.5, thislabel,'FontName','Helvetica','FontSize',16)
%plot(MTFspacing,MTF3,'b')
plot(MTFspacing,MTF,'r')

axis([300 900 0 1])
% set(gca,'xtick',[])

hold off

%% 
frame=getframe(h4);
% % now capture the image that was displayed and write it to a pdf file

imwrite(frame.cdata,['Spec' accomMTFfigname],'Quality',100)

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
%% final plots for paper
h5=figure('Position', [165 39 473 767]);

first=5;
second=midframe;
third=length(MTFspacing);

% these are all carefully tweaked to get alignment. Don't change!
% left-bottom corner, then width, height scaled 0 to 1
position1=[-0.02, 0.55, 1, 0.1]; % top panel -  color image
position2=[-0.02, 0.42, 1, 0.1]; % image 1
position3=[0.1, 0.32, 0.8, 0.1]; % linecut 1
position4=[-0.02, 0.22, 1, 0.1]; % image 2
position5=[0.1, 0.12, 0.8, 0.1]; % linecut 2
%position6=[0.1, 0.02, 0.8, 0.34]; % MTF image

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
plot(linecutf2,'r'); 
set(gca,'xtick',[],'ytick',[])
%text(20,0.4,'\bf 500 nm','FontName','Times','FontSize',16)
%text(20,0.2,' \bf focus','FontName','Times','FontSize',16)
ylim([0 1.1])

% subplot('Position',position6);
% 
% plot(accommodation,MTF,'b')
%     hold on
%     plot(accommodation,spectrum,'r')
%     %legend('MTF','spectrum','FontName','HELVETICA','FontSize',16)
%     plot(accommodation,MTF,'k*')
%     plot(accommodation,spectrum,'k.')
%     xlabel('accommodation, mm','FontName','HELVETICA','FontSize',16)
%     ylabel('MTF','FontName','Helvetica','FontSize',16)
%     axis([-0.35 0.2 -0.1 1.3]);
%     % plottext=[MTFlabel ' with threshold ' num2str(MTFcut)]
%     plottext=[MTFlabel]
%     text(-0.1,1.1,plottext,'FontName','HELVETICA','FontSize',20) % for some reason this does not appear in saved JPEG!
%     hold off
%     shg
%     
truesize([200,1000])
shg

frame=getframe(h5);
% % now capture the image that was displayed and write it to a JPG file
imwrite(frame.cdata,figname,'Quality',100) % tried 'Mode','lossless' but it fails for some reason. 







    

