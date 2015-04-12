% PSF.m, simulates PSF for chromatic aberration for Cephalopods
% Released under GPL v2 license by Christopher Stubbs and Alexander Stubbs
% astubbs@berkeley.edu 

% version 2, April 12 29 2015

%%%%%%%%%%%%% 
% This program computes PSF for the chromatic aberrations of an image, given
% these inputs:
%  1) pupil structure. Currently only axially symmetrical circular pupils
%  2) solar photon spectrum
%  3) array of three reflectance spectra for test pattern image
%  4) opsin response function vs. wavelength
%  5) lens diameter and focal length, for chromatic blur computation
%  6) depth-dependent attentuation of sunlight in seawater.
%  7) depth in ocean for this computation, for attentuation calc.
%  8) blue-yellow intensity ratio, called colorbalance, for blue/yellow checkerboard.
%  9) measured fractional change in focal length vs. wavelength. 
% 10) vector MTFspacing that contains best-focus wavelengths to be used

% User can specify the pupil mask. We assume an illumination spectrum of
% attenuated sunlight, i.e. white reflectance from the object. For this
% computation we do extrapolate the Jagger and Sands chromatic aberration
% changes in focal length down to 350 nm, since we need to include all the
% light within the opsin response function

% This version loops through focus settings and creates blurred images at
% each focus setting. Uses partial-annular and on-axis circular pupils.

% other spectral plots of interest can also be generated from the data. 

% All output files should appear in the current working directory, which
% should be added to the MATLAB path. 

% There are a number of other logical switches that can be used to turn
% things on and off. These include
% diagnosticplots - if set true a number of diagnostic plots are produced.
% half - makes half-annulus pupil
% diffraction - turns diffraction on and off

%% intitialize. Edit these parameters as needed. 
clear all
subplot(1,1,1)
close all

%%%%%%%%%%%%%%% these parameters don't change very frequently %%%%%%%%
depth = 3; % water depth in meters. default is 3
lambda=350:900; % array of discrete wavelengths, one nm spacing. default is 350:900 and this is essential!
lambda=lambda'; % make it a column vector
focallength=12; % focal length in mm. We assume f/1.2 beam at full aperture and so a 10 mm diameter lens. default is 12
diagnosticplots=false; % set to true for verbose plotting
testrows=800; % sets size of test pattern. default is 800 x 800. spatial scale is 1 micron
testcols=800;

%%%% Frequently edited parameters begin here: pupil and focus choices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin selection section. INTEGERS ONLY since these are used in file names %%%%%%%%%%%
bestfocuslambda=480; % best-focus wavelength in nm. Opsin peak is 500 **INTEGERS ONLY**
pupilouterradius=4; % pick 4 for full aperture or annular, 1 for small central  **INTEGERS ONLY**. units are mm
pupilinnerradius=3; % pick zero for central pupil, 3 for annular. make sure it's smaller than outer radius above!
half=true; % set this to true for half-pupil, to false for axisymmetric pupil. 
diffraction=false; % switches diffraction on and off. Only really matters for small pupil (~< 2 mm outer) 
%%%%%%%%%%%%% end pupil selection section %%%%%%%%%%%%

PSFlabel=['PSF_' num2str(bestfocuslambda) '_' num2str(pupilouterradius) '_' num2str(pupilinnerradius)]; % no spaces, no slashes! Pupil type is appended automatically. Use color info. 

%%%%%%%%% end of section where choices are edited. %%%%%%%%%%%%

% housekeeping for files and names

savename=[PSFlabel '.mat']; % name of workspace file that is saved
figname=[PSFlabel '.jpg']; % name of high quality JPG file that is made
datafilename=[PSFlabel '.dat']; % name of FWHM data file. wavelength and pupil size and FWHM go in there

% wipe out existing files so we over-write them. 
delete(figname);
delete(savename);

 warning('off','all') % suppress warnings to command line

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
accommodation=1E-2*(-5.4676e-05*lambda.^2 + 0.0794*lambda -27.1047)*focallength; % this is an array


% % % % this is illumination spectrum on 1 nm spacing, 350 to 900 nm  
fullspectrum=opsinresponse.*solarspectrum;
fullspectrum=fullspectrum/max(fullspectrum); % normalize it

if (diagnosticplots)
    figure(1)
    plot(lambda, fullspectrum, 'k')
    shg

    figure(2)
    plot(lambda,accommodation,'b.')
    shg
end

%% create PSF test pattern. Size is determined by testrows and testcols, set in initialization block at top. 

focuscounter=0;


for bestfocuslambda=350:50:650 % keep for loop
    
blurredimage=zeros(testrows,testcols,301); % this is the temporary image stack made for each wavelength.  
blurredimage2=blurredimage;
scaledimage=blurredimage;
    focuscounter=focuscounter+1; % this is an index for the focus loop
%% compute offset from best-focus wavelength. 
% Compute the index in wavelength vector that is closest to this best-focus
% wavelength. The floor() function rounds to an integer
bestindex=floor(bestfocuslambda)-349; % this is the wavelength index that corresponds to best-focus lambda this time through loop. 350 nm is index = 1 

% now compute magnitude of focus offset from that of best-focus wavelength. 
% sign convention is with positive focusshift being further from lens
% this quantity is negative for wavelengths that focus closer to the lens
% than the best-focus rays. We want the units to be mm
relativefocusshift=focusshiftmm-focusshiftmm(bestindex); % this is a vector with a value for each wavelength, after subtracting off the accommodation setting
                                                         % setting for this best-focus wavelength
                                                          
    
%% now create a an out-of-focus pupil image at each wavelength
midx=floor(testcols/2); % find midpoint of image
midy=midx;
integral=zeros(301,1);


for wavecounter=1:301 % this steps through wavelengths from 350 to 650.   Use
                      % counter=101:301 for 450 to 650. 

    
      % now compute pupil image for each wavelength. Pupil image size is
      % smaller by a factor of pupilscaling
      
        pupilscaling(wavecounter)=abs(relativefocusshift(wavecounter))/focallength; %
        % compute new pupil size, taking into account it's given in mm
        % but pixels are microns
        pupilrad1(wavecounter)=1000*pupilouterradius*pupilscaling(wavecounter);
        pupilrad2(wavecounter)=1000*pupilinnerradius*pupilscaling(wavecounter);
       
        
        for xcounter=1:testcols
            for ycounter=1:testrows
                pixrad=sqrt((xcounter-midx).^2+(ycounter-midy)^2);
                if ((pixrad<=pupilrad1(wavecounter)) && (pixrad>=pupilrad2(wavecounter)))
                    blurredimage(ycounter,xcounter,wavecounter)=1;
                    integral(wavecounter)=integral(wavecounter)+1;
                    if ((half) && (relativefocusshift(wavecounter)<=0) && (ycounter>=midy)) % shorter wavelength than optimal
                        blurredimage(ycounter,xcounter,wavecounter)=0; % wipe out top half of pupil image. note y is row direction
                         integral(wavecounter)=integral(wavecounter)-1; % decrement the counter used to track pixels
                    end
                    if ((half) && (relativefocusshift(wavecounter)>0) && (ycounter<midy)) % longer wavelength than optimal
                        blurredimage(ycounter,xcounter,wavecounter)=0; % wipe out top bottom portion of pupil image
                        integral(wavecounter)=integral(wavecounter)-1;
                         % lowerkill=ycounter
                    end % end of half-pupil wipeout 
                end % end of pixel radius loop
            end % end of ycounter loop
        end % end of x counter loop
        
%         % trap on condition where it's a point smaller than a pixel, add
%         % some flux in the center
        if (integral(wavecounter)==0)
            blurredimage(midx-2:midx+2,midy-2:midy+2,wavecounter)=1;
            integral(wavecounter)=sum(sum(blurredimage(:,:,wavecounter)));
        end
        
        % this next part is a bit clunky but was useful for debugging
        % purposes
        
        % normalize the pupil values to give same flux at each wavelength, before multiplying by spectrum. 
        blurredimage2(:,:,wavecounter)=blurredimage(:,:,wavecounter)/(integral(wavecounter)+1);
        
        if(diffraction) % not really used much- unverified. 
            % figure out diffraction-limited Gaussian and blur the image by
            % that amount. Assume it's set by outer pupil diameter and
            % symmetrical, for now. Focal length and pupil diameters are in mm
            % and lambda in nm. 
            % Use FWHM=f-number times lambda, and then sigma=FWHM/2.35.
            % compute diffraction FWHM in microns
            diff_fwhm=1E6*(focallength/pupilouterradius)*(lambda(wavecounter)*1E-9);
            % convert to a sigma
            sigma=diff_fwhm/2.35
            % make a kernel size that scales with this
            ksize=floor(3*sigma);
            % create the kernel
            kernel=fspecial('gaussian',[ksize ksize],sigma);
            % filter the image
            blurredimage2=imfilter(blurredimage2,kernel,'same');
       
        end
        
        % scale to flux at each wavelength, for this focus setting
        scaledimage(:,:,wavecounter)=blurredimage2(:,:,wavecounter).*fullspectrum(wavecounter);
        
        % verify that the sum of this reproduces the spectrum:
        % plot(squeeze(sum(sum(scaledimage))))
        
    % send a progress report to the command line
   status=['at wavelength ' num2str(lambda(wavecounter))] % spit out loop counters for diagnostic purposes
end % end of wavelength loop


%% Let's look at the final resulting image, added up over all wavelengths

finalimage=sum(scaledimage,3); % this adds up the blurred hypercube along the wavelength direction
finalimage=finalimage/(max(max(finalimage))); % normalize to unity at peak
% compute FWHM by finding place where intensity is 0.5
[diff FWHMpix]=min(min(abs(0.5-finalimage)));
delta=midx-FWHMpix;
FWHM=2*delta

% h2=figure(2); 
% imshow(finalimage)
% labelstring=num2str(bestfocuslambda);
% label2=[' PSF at ' labelstring ' nm focus'];
% title(label2, 'FontName','HELVETICA','FontSize',16)
% shg

h3=figure(3);
surf(log(finalimage+0.001))
colormap default
shading interp
shg

dx=1:testcols;
dx=dx-midx;

h4=figure(4);
plot(dx,finalimage(:,midy),'b')
xlabel('r, microns', 'FontName','HELVETICA','FontSize',16)
ylabel('PSF intensity', 'FontName','HELVETICA','FontSize',16)
focuslabel=['Best-focus set to ' num2str(bestfocuslambda) ' nm' ] ;
pupillabel=['Pupil diameter of ' num2str(2*pupilouterradius) ' mm'];
FWHMlabel=['FWHM= ' num2str(FWHM) ' microns'];
text(-140,1.15,focuslabel,'FontName','HELVETICA','FontSize',16)
text(-140,1.1,pupillabel,'FontName','HELVETICA','FontSize',16)
text(-140,1.05,FWHMlabel,'FontName','HELVETICA','FontSize',16)


axis([-150 150 0 1.2]);
shg


% now capture the radial plot image that was displayed and write it to a jpg file
frame=getframe(h4);
imwrite(frame.cdata,figname,'Quality',100)



%% compute centroid and encircled energy. This has tricky index maneuvers...
x=1:testcols;
y=1:testrows;

phivy=sum(finalimage,2); % collapse final image along columns to get y dependence 
phivx=sum(finalimage,1); % collapse final image along rows to get x dependence

xcenter=sum(phivx.*x)/(sum(phivx));
ycenter=sum(phivy'.*y)/(sum(phivy)); % need to rotate to a single row

% compute a 2-d array that has distances from this centroid

radarray=zeros(testrows,testcols);

for xcounter=1:testcols
    for ycounter=1:testrows
        radarray(ycounter,xcounter)=sqrt((xcounter-xcenter).^2+(ycounter-ycenter).^2);
    end
end

% now convert this into a 1-d vector, reads along each row and makes a
% single column vector
radvec=radarray(:);

% this next step returns both a sorted array, and an index array that gives
% sorting order
[sortedradvec, sortindex]=sort(radvec);

% make a 1-d column vector out of intensity array as well, in same sequence as radvec 

PSFvec=finalimage(:);

% and now (cool, eh?) sort this PSF vector by same index array that sorted radius vector, 
% since it should have same initial order. This is a radial profile of
% flux, about centroid.

sortedPSFvec=PSFvec(sortindex);

% make a plot
% figure(5)
% plot(sortedradvec,sortedPSFvec)
% shg

%% radial integral of enclosed flux

EE=cumsum(sortedPSFvec);
EE=EE/max(EE); % normalize this

EEarray(focuscounter,:)=EE'; % make an array with a line per focus setting

figure(6)
plot(sortedradvec,EE)
axis([0 250 0 1.1])

[EE80, index80]=min(abs(0.8-EE)); % find the entry in EE array closest to 0.8 = 80 percent, and its index
enclosed80=sortedradvec(index80) % this should be the radius that enclosed 80% of flux

Enclosed80values(focuscounter)=enclosed80; % make an array of these enclosed energy values

pause(1) % pause to look at EE number on command line

%% store the resulting data structures

save(savename,'finalimage')

% append to file with the results

outvec1=[bestfocuslambda pupilouterradius pupilinnerradius enclosed80 FWHM xcenter ycenter half diffraction];
dlmwrite('PSF.dat',outvec1,'-append');

outvec2=EE(1:250);
outvec2=outvec2' ; % make it a row vector


dlmwrite('EE.dat',outvec2,'-append'); % write out the first 250 elements of the integrated energy profile

end % end of bestfocus loop

%% make multiple focus enclosed energy plot

h10=figure(10);
arrayindex=[1 2 3 4]

    plot(sortedradvec,EEarray(arrayindex(1),:),'k',sortedradvec,EEarray(arrayindex(2),:),'g',sortedradvec,EEarray(arrayindex(3),:),'b',...
        sortedradvec,EEarray(arrayindex(4),:),'r');
    set(gca,'FontSize',20)
    legend('350','400','450','500')
    xlim([0 250])
    hold on
    plot(sortedradvec, 0.8+0*radvec,'c:');
    hold off
xlabel('r, microns')
ylabel('EE')
% title('Encircled Energy vs. Radius, Small Pupil')
shg

%% print enclosed energy numbers to command line
% convert to Gaussian-equivalent FWHM with FWHM = R80*1.3
FWHM=Enclosed80values*1.3;
accommodations=350:50:650;
disp('Accommodation wavelength and FWHM (microns)')

outvec=[accommodations', FWHM']