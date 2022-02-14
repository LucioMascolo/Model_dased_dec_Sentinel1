clear all; 
close all; 
clc
%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Decompose the Sentinel-1 Stokes vector with the method published in the article:
%% L. Mascol, S.R. Cloude, J.M. Lopez-Sanchez, "Model-based decomposition of dual-pol SAR data: application to Sentinel-1", 
%% Trans. and Geosci. Remote Sens., in print.
%%
%% Inputs:
%% - Elements of the wave coherency matrix C2 (either from HH-HV or VH-VV data)
%
%% Outputs:
%% - mv, mp, alpha and delta
%% - RGB, alpha and HSV images
%% 
%% Authors:
%% - Lucio Mascolo, Institute for Computer Research (IUII), University of Alicante
%% - Shane R. Cloude, Applied Electromagnetics Consultants (AELC)
%% - Juan Manuel Lopez-Sanchez, Institute for Computer Research (IUII), University of Alicante
%% Date: 11/02/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  SET INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - dual-pol mode (either HH-HV or VH-VV):
pol_mode   = 'VH-VV';   

%% - path to the directory containing the elements of the C2 matrix:
dir_c2     = '/Users/luciomascolo/Desktop/test/C2/';

%% - image size (nline x nsample)
nline      = 9582;
nsample    = 8983;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXECUTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% - load the C2 elements

fid              = fopen (strcat(dir_c2,'C11.bin'), 'r');
C11              = fread(fid,[nsample,nline],'single')';
fclose(fid);

fid              = fopen (strcat(dir_c2,'C12_real.bin'), 'r');
C12_real         = fread(fid,[nsample,nline],'single')';
fclose(fid);

fid              = fopen (strcat(dir_c2,'C12_imag.bin'), 'r');
C12_imag         = fread(fid,[nsample,nline],'single')';
fclose(fid);

fid              = fopen (strcat(dir_c2,'C22.bin'), 'r');
C22              = fread(fid,[nsample,nline],'single')';
fclose(fid);

%% - build the Stokes vector from C2
s1               = C11 + C22;
s2               = C11 - C22;
s3               = 2*C12_real;
s4               = 2*C12_imag;
    
%% - apply the decomposition
if strcmp(pol_mode,'HH-HV')==1
    F = 1/2;
else
    F = - 1/2;
end

%% - solve the quadratic (volume term)
a                = 0.75;
b                = -2*(s1-(F*s2));
c                = s1.*s1 - s2.*s2 -s3.*s3 - s4.*s4;
delta1           = b.^2;
delta2           = 4*a*c;
delta            = delta1-delta2;

mv1              = -b + sqrt(delta);
mv1              = mv1/(2*a);
mv2              = (-b - sqrt(delta));
mv2              = mv2/(2*a);

%% - check which solution satisfes s1 > mv
ind1             = find(s1>0); % exclude pixels outside the imaged scene
flag(1)          = isempty(find(s1(ind1)>mv1(ind1), 1));
flag(2)          = isempty(find(s1(ind1)>mv2(ind1), 1));
    
%% solution:
if flag(1) == 0
        mv = mv1;
else
        mv = mv2;
end

%% - obtain mp (polarized term)
mp        = s1-mv;

%% - obtain alpha and delta for the polarized term
if strcmp(pol_mode,'HH-HV')==1
    alpha                = .5*acos((s2-F.*mv)./mp); 
else
    alpha                = .5*acos(-(s2-F.*mv)./mp);
end

delta       = angle(s3+s4*1i);

%% - show RGB images: 
mpdb        = 10*log10(mp);
mvdb        = 10*log10(mv);
ratdb       = mpdb-mvdb;

rgb(:,:,1)  = mat2gray(mpdb,[-20,-2]);
rgb(:,:,2)  = mat2gray(mvdb,[-20,-2]);
rgb(:,:,3)  = mat2gray(ratdb,[-10,10]);

figure();
imagesc(rgb)
title([pol_mode,' RGB'])
axis off

%% - show alpha image:
figure();
imagesc(alpha*180/pi)
title('alpha angle (degrees)') 
axis off
colormap(jet)
colorbar

%% - show HSV image: 
%% Hue (delta):
pscale        = (delta+pi)/(2*pi);% scale the phase in the range 0 to 1
hsv(:,:,1)    = pscale;
%% Saturation (cross-polarized coherence)
ro            = sqrt(C12_real.*C12_real+C12_imag.*C12_imag);
ro            = ro./sqrt(C11.*C22);
hsv(:,:,2)    = ro;
%% Value (span):
spandB        = 10*log10(C11+C22);
hsv(:,:,3)    = mat2gray(mpdb,[-25,5]); % scale between 0 and 1
hsv_final     = hsv2rgb(hsv); % final hsv: depending on the image size, it might take a while
figure()
imagesc(hsv_final)
axis off
title([pol_mode,' HSV']) 


    