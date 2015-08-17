clear all
I=imread ('barbara.jpg');
sz = size(I);
pcheck = size(sz);
if pcheck(1,2) == 3
    grays = rgb2gray(I);
else
    grays = I;
end
grays = imresize(grays,[192 192]);
size_of_image = size(grays);        %gray scale image
figure,imshow(grays),title('Original Image')

nvar = .001;                                      %variance of noise
ngrays = imnoise(grays,'gaussian',0,nvar);        %white gaussian noise
figure,imshow(ngrays),title('noised Image')
Numlev = 2;         %level of decomposition of wavelet domain
[coeffn,bookem] = wavedec2(ngrays,Numlev,'db4');
[coeff,bookem] = wavedec2(grays,Numlev,'db4');
nsd = sqrt(nvar);       % standard deviation of noise
%maxlim = 2*(nsd*sqrt(2*log()))    %threshold value's limit by visu shrink

levs(1) = bookem(1,1)*bookem(1,2);  %size of final Approximation matrix 
for i=1:Numlev
levs(i+1) = bookem(i+1,1)*bookem(i+1,2);      %size of every level
end

%%%%%%%%%%%%%%%%%%%%%%% Noisy coefficients %%%%%%%%%%%%%%%%%%%%%

ncoefflev1 = [];                 %for noisy coefficient of level 1
ncoefflev1(1,:) = coeffn(1,1:levs(1));
n = 1+levs(1);
n1 = 2*levs(2);
for i=2:4
    ncoefflev1(i,:) = coeffn(1,n:n1);
    n = n+levs(2);
    n1 = n1+levs(2);
end

ncoefflev2 = [];                 %for noisy Coefficients of level 2
n = 1+levs(1)+3*levs(2);
n1 = levs(1)+3*levs(2)+levs(3);
for i=1:3
    ncoefflev2(i,:) = coeffn(1,n:n1);
    n = n+levs(3);
    n1 = n1+levs(3);
end

%%%%%%%%%%%%%%%%%%%%%%% Noise free coefficients %%%%%%%%%%%%%%%%%%%%%

coefflev1 = [];                 %for simple coefficient of level 1
coefflev1(1,:) = coeff(1,1:levs(1));
n = 1+levs(1);
n1 = 2*levs(2);
for i=2:4
    coefflev1(i,:) = coeff(1,n:n1);
    n = n+levs(2);
    n1 = n1+levs(2);
end

coefflev2 = [];                 %for simple Coefficients of level 2
n = 1+levs(1)+3*levs(2);
n1 = levs(1)+3*levs(2)+levs(3);
for i=1:3
    coefflev2(i,:) = coeff(1,n:n1);
    n = n+levs(3);
    n1 = n1+levs(3);
end

%%%%%%%%%%%%%%%%%%%%%%% Modified coefficients %%%%%%%%%%%%%%%%%%%%%

 n=15;          %number of Particles
 nv = 4;        %number of variables
 lim = [0 1000;0 1;0.1 10;0.1 4]; %lower and upper bound of variables
 vcf = 2;       %velocity clamping factor
 cc = 3;        %cognitive constant
 sc = 2;        %social constant
 miniw = .4;    %Min Inertia weight
 maxiw = .9;    %Max Inertia weight
 num = 20;  

mcoefflev1 = [];
n = 15;
for i=1:4
    tempnf = coefflev1(i,:);
    tempn = ncoefflev1(i,:);
    [mcoefflev1(i,:),nestlev1(:,:,i)] = Particle_Swarm_Optimizationsc4(n,nv,lim,@add,'min',vcf,cc,sc,miniw,maxiw,num,tempnf,tempn);
end

mcoefflev2 = [];
for i=1:3
    tempnf = coefflev2(i,:);
    tempn = ncoefflev2(i,:);
    [mcoefflev2(i,:),nestlev2(:,:,i)] = Particle_Swarm_Optimizationsc4(n,nv,lim,@add,'min',vcf,cc,sc,miniw,maxiw,num,tempnf,tempn);
end

mcoeffmat = [];         %reconstruction of original signal
mcoeffmat = cat(2,mcoefflev1(1,:),mcoefflev1(2,:),mcoefflev1(3,:),mcoefflev1(4,:),mcoefflev2(1,:),mcoefflev2(2,:),mcoefflev2(3,:));
mimage = waverec2(mcoeffmat,bookem,'db4');
main_image = uint8(mimage);

figure,imshow(main_image),title('Enhanced Image')

modlev1 = [];
modlev2 = [];
%%%%%%%%%%%%%%%%%%%%%%% PSNR & MSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qsize = size(nestlev1,1);
for i=1:qsize
    modlev1(1,:) = thresholdfuncsc4(nestlev1(i,:,1),ncoefflev1(1,:));
    modlev1(2,:) = thresholdfuncsc4(nestlev1(i,:,2),ncoefflev1(2,:));
    modlev1(3,:) = thresholdfuncsc4(nestlev1(i,:,3),ncoefflev1(3,:));
    modlev1(4,:) = thresholdfuncsc4(nestlev1(i,:,4),ncoefflev1(4,:));
    modlev2(1,:) = thresholdfuncsc4(nestlev2(i,:,1),ncoefflev2(1,:));
    modlev2(2,:) = thresholdfuncsc4(nestlev2(i,:,2),ncoefflev2(2,:));
    modlev2(3,:) = thresholdfuncsc4(nestlev2(i,:,3),ncoefflev2(3,:));
    
    modmat = cat(2,modlev1(1,:),modlev1(2,:),modlev1(3,:),modlev1(4,:),modlev2(1,:),modlev2(2,:),modlev2(3,:));
    modimage = waverec2(modmat,bookem,'db4');
    main_image = uint8(modimage);
    
    D = abs(uint8(main_image) - uint8(grays)).^2;
    pmse(1,i) = sum(D(:))/numel(main_image);
    ppsnr(1,i) = 10*log10(255*255/pmse(1,i));
end


mean_Optimized = mean2(main_image)
var_optimzed = std2(main_image)

D = abs(uint8(main_image) - uint8(grays)).^2;
mse = sum(D(:))/numel(main_image)
psnr = 10*log10(255*255/mse)

% %%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%
% summse = sum(cat(1,mselev1,mselev2));
% sumpsnr = sum(cat(1,psnrlev1,psnrlev2));
t = 1:1:qsize;

figure,plot(t,pmse);
xlabel('Number of Iterations');
ylabel('Mean Square Error');
title('MSE vs Iterations')

figure,plot(t,ppsnr);
xlabel('Number of Iterations');
ylabel('Peak Signal To Noise Ratio');
title('PSNR vs Iterations')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%main_image_gray = rgb2gray(main_image);
%grays = rgb2gray(grays);
 K = [0.6 0.6];
 [mssim, ssim_map] = ssim_index(grays,main_image,K);
 mssim
 
 [FSIM, FSIMc] = FeatureSIM(grays,main_image);
 FSIM


% coeffmat2lv1 = [];
% coeffmat2lv1(1,:)  = coeff(1,1:levs(1));  %row vector of final approx. matrix
% n=1+levs(2);
% n1 = n+levs(2);
% for i=1:3
%     coeffmat2lv1(i+1,:)  = coeff(1,n:n1);
%     n = n+levs(2);
%     n1 = n+levs(2);
% end
%  
% n = 5*levs(2);
% n1 = 5*levs(2) + levs(3);
% coeffmatlv2 =[];
% for j=5:(3*Numlev)+1 
%   coeffmatlv2(j,:)  = coeff(1,n:n1);
%   n1 = n1+levs(3);
% end
% 

