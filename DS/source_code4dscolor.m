clear all
I=imread ('peppers_color.jpg');
r = I(:,:,1);
g = I(:,:,2);
b = I(:,:,3);
rr = imresize(r,[192 192]);
gr = imresize(g,[192 192]);
br = imresize(b,[192 192]);
Im(:,:,1) = rr;
Im(:,:,2) = gr;
Im(:,:,3) = br;
%rgrays = imresize(rgrays,[192 192]);
% ggrays = imresize(ggrays,[192 192]);
% bgrays = imresize(bgrays,[192 192]);
size_of_image = size(rr);        %gray scale image
figure,imshow(Im),title('Original Image')

nvar = .001;                                      %variance of noise
ngrays = imnoise(Im,'gaussian',0,nvar);  %white gaussian noise
figure,imshow(ngrays),title('noised Image')
rngrays  = ngrays(:,:,1);
gngrays  = ngrays(:,:,2);
bngrays  = ngrays(:,:,3);

Numlev = 2;         %level of decomposition of wavelet domain
[rcoeffn,rbookem] = wavedec2(rngrays,Numlev,'db4');     %Red Pixels
[rcoeff,rbookem] = wavedec2(rr,Numlev,'db4');

[gcoeffn,gbookem] = wavedec2(gngrays,Numlev,'db4');     %Green pixels
[gcoeff,gbookem] = wavedec2(gr,Numlev,'db4');

[bcoeffn,bbookem] = wavedec2(bngrays,Numlev,'db4');     %Blue Pixels
[bcoeff,bbookem] = wavedec2(br,Numlev,'db4'); 

nsd = sqrt(nvar);       % standard deviation of noise
%maxlim = 2*(nsd*sqrt(2*log()))    %threshold value's limit by visu shrink

levs(1) = rbookem(1,1)*rbookem(1,2);  %size of final Approximation matrix 
for i=1:Numlev
levs(i+1) = rbookem(i+1,1)*rbookem(i+1,2);      %size of every level
end

%%%%%%%%%%%%%%%%%%%%%%% RED Noisy coefficients %%%%%%%%%%%%%%%%%%%%%

rncoefflev1 = [];                 %for noisy coefficient of level 1
rncoefflev1(1,:) = rcoeffn(1,1:levs(1));
n = 1+levs(1);
n1 = 2*levs(2);
for i=2:4
    rncoefflev1(i,:) = rcoeffn(1,n:n1);
    n = n+levs(2);
    n1 = n1+levs(2);
end

rncoefflev2 = [];                 %for noisy Coefficients of level 2
n = 1+levs(1)+3*levs(2);
n1 = levs(1)+3*levs(2)+levs(3);
for i=1:3
    rncoefflev2(i,:) = rcoeffn(1,n:n1);
    n = n+levs(3);
    n1 = n1+levs(3);
end



%%%%%%%%%%%%%%%%%%%%%%% RED Noise free coefficients %%%%%%%%%%%%%%%%%%%%%

rcoefflev1 = [];                 %for simple coefficient of level 1
rcoefflev1(1,:) = rcoeff(1,1:levs(1));
n = 1+levs(1);
n1 = 2*levs(2);
for i=2:4
    rcoefflev1(i,:) = rcoeff(1,n:n1);
    n = n+levs(2);
    n1 = n1+levs(2);
end

rcoefflev2 = [];                 %for simple Coefficients of level 2
n = 1+levs(1)+3*levs(2);
n1 = levs(1)+3*levs(2)+levs(3);
for i=1:3
    rcoefflev2(i,:) = rcoeff(1,n:n1);
    n = n+levs(3);
    n1 = n1+levs(3);
end

%%%%%%%%%%%%%%%%%%%%%%% RED Modified coefficients %%%%%%%%%%%%%%%%%%%%%

rmcoefflev1 = [];

% ds(method,fnc,mydata,popsize,dim,low,up,maxcycle)
% method
% 1: Bijective DSA
% 2: Surjective DSA
% 3: Elitist DSA (strategy 1)
% 4: Elitist DSA (strategy 2)

n = 15;
lower_limit = [0,0,.1,.1];          
upper_limit = [1000,1,20,4];
method = 1 ;
dimensions = 4;
maxIter = 100;

for i=1:4
    tempnf = rcoefflev1(i,:);
    tempn = rncoefflev1(i,:);
    [rmcoefflev1(i,:),rmselev1(i,:),rpsnrlev1(i,:)] = dssc4(method,@add,[],n,dimensions,lower_limit,upper_limit,maxIter,tempnf,tempn,nsd);
end

rmcoefflev2 = [];
for i=1:3
    tempnf = rcoefflev2(i,:);
    tempn = rncoefflev2(i,:);
    [rmcoefflev2(i,:),rmselev2(i,:),rpsnrlev2(i,:)] = dssc4(method,@add,[],n,dimensions,lower_limit,upper_limit,maxIter,tempnf,tempn,nsd);
end

rmcoeffmat = [];         %reconstruction of original signal
rmcoeffmat = cat(2,rmcoefflev1(1,:),rmcoefflev1(2,:),rmcoefflev1(3,:),rmcoefflev1(4,:),rmcoefflev2(1,:),rmcoefflev2(2,:),rmcoefflev2(3,:));
rimage = waverec2(rmcoeffmat,rbookem,'db4');


%%%%%%%%%%%%%%%%%%%%%%% Green Noisy coefficients %%%%%%%%%%%%%%%%%%%%%

gncoefflev1 = [];                 %for noisy coefficient of level 1
gncoefflev1(1,:) = gcoeffn(1,1:levs(1));
n = 1+levs(1);
n1 = 2*levs(2);
for i=2:4
    gncoefflev1(i,:) = gcoeffn(1,n:n1);
    n = n+levs(2);
    n1 = n1+levs(2);
end

gncoefflev2 = [];                 %for noisy Coefficients of level 2
n = 1+levs(1)+3*levs(2);
n1 = levs(1)+3*levs(2)+levs(3);
for i=1:3
    gncoefflev2(i,:) = gcoeffn(1,n:n1);
    n = n+levs(3);
    n1 = n1+levs(3);
end

%%%%%%%%%%%%%%%%%%%%%%% Green Noise free coefficients %%%%%%%%%%%%%%%%%%%%%

gcoefflev1 = [];                 %for simple coefficient of level 1
gcoefflev1(1,:) = gcoeff(1,1:levs(1));
n = 1+levs(1);
n1 = 2*levs(2);
for i=2:4
    gcoefflev1(i,:) = gcoeff(1,n:n1);
    n = n+levs(2);
    n1 = n1+levs(2);
end

gcoefflev2 = [];                 %for simple Coefficients of level 2
n = 1+levs(1)+3*levs(2);
n1 = levs(1)+3*levs(2)+levs(3);
for i=1:3
    gcoefflev2(i,:) = gcoeff(1,n:n1);
    n = n+levs(3);
    n1 = n1+levs(3);
end

%%%%%%%%%%%%%%%%%%%%%%% Green Modified coefficients %%%%%%%%%%%%%%%%%%%%%

gmcoefflev1 = [];
n = 15;
for i=1:4
    tempnf = gcoefflev1(i,:);
    tempn = gncoefflev1(i,:);
    [gmcoefflev1(i,:),gmselev1(i,:),gpsnrlev1(i,:)] = dssc4(method,@add,[],n,dimensions,lower_limit,upper_limit,maxIter,tempnf,tempn,nsd);
end

gmcoefflev2 = [];
for i=1:3
    tempnf = gcoefflev2(i,:);
    tempn = gncoefflev2(i,:);
    [gmcoefflev2(i,:),gmselev2(i,:),gpsnrlev2(i,:)] = dssc4(method,@add,[],n,dimensions,lower_limit,upper_limit,maxIter,tempnf,tempn,nsd);
end

gmcoeffmat = [];         %reconstruction of original signal
gmcoeffmat = cat(2,gmcoefflev1(1,:),gmcoefflev1(2,:),gmcoefflev1(3,:),gmcoefflev1(4,:),gmcoefflev2(1,:),gmcoefflev2(2,:),gmcoefflev2(3,:));
gimage = waverec2(gmcoeffmat,gbookem,'db4');


%%%%%%%%%%%%%%%%%%%%%%% BLUE Noisy coefficients %%%%%%%%%%%%%%%%%%%%%

bncoefflev1 = [];                 %for noisy coefficient of level 1
bncoefflev1(1,:) = bcoeffn(1,1:levs(1));
n = 1+levs(1);
n1 = 2*levs(2);
for i=2:4
    bncoefflev1(i,:) = bcoeffn(1,n:n1);
    n = n+levs(2);
    n1 = n1+levs(2);
end

bncoefflev2 = [];                 %for noisy Coefficients of level 2
n = 1+levs(1)+3*levs(2);
n1 = levs(1)+3*levs(2)+levs(3);
for i=1:3
    bncoefflev2(i,:) = bcoeffn(1,n:n1);
    n = n+levs(3);
    n1 = n1+levs(3);
end

%%%%%%%%%%%%%%%%%%%%%%% BLUE Noise free coefficients %%%%%%%%%%%%%%%%%%%%%

bcoefflev1 = [];                 %for simple coefficient of level 1
bcoefflev1(1,:) = bcoeff(1,1:levs(1));
n = 1+levs(1);
n1 = 2*levs(2);
for i=2:4
    bcoefflev1(i,:) = bcoeff(1,n:n1);
    n = n+levs(2);
    n1 = n1+levs(2);
end

bcoefflev2 = [];                 %for simple Coefficients of level 2
n = 1+levs(1)+3*levs(2);
n1 = levs(1)+3*levs(2)+levs(3);
for i=1:3
    bcoefflev2(i,:) = bcoeff(1,n:n1);
    n = n+levs(3);
    n1 = n1+levs(3);
end

%%%%%%%%%%%%%%%%%%%%%%% BLUE Modified coefficients %%%%%%%%%%%%%%%%%%%%%

bmcoefflev1 = [];
n = 15;
for i=1:4
    tempnf = bcoefflev1(i,:);
    tempn = bncoefflev1(i,:);
    [bmcoefflev1(i,:),bmselev1(i,:),bpsnrlev1(i,:)] = dssc4(method,@add,[],n,dimensions,lower_limit,upper_limit,maxIter,tempnf,tempn,nsd);
end

bmcoefflev2 = [];
for i=1:3
    tempnf = bcoefflev2(i,:);
    tempn = bncoefflev2(i,:);
    [bmcoefflev2(i,:),bmselev2(i,:),bpsnrlev2(i,:)] = dssc4(method,@add,[],n,dimensions,lower_limit,upper_limit,maxIter,tempnf,tempn,nsd);
end

bmcoeffmat = [];         %reconstruction of original signal
bmcoeffmat = cat(2,bmcoefflev1(1,:),bmcoefflev1(2,:),bmcoefflev1(3,:),bmcoefflev1(4,:),bmcoefflev2(1,:),bmcoefflev2(2,:),bmcoefflev2(3,:));
bimage = waverec2(bmcoeffmat,bbookem,'db4');

main_image = [];
main_image(:,:,1) = rimage;
main_image(:,:,2) = gimage;
main_image(:,:,3) = bimage;
main_image = uint8(main_image);
figure,imshow(main_image),title('Enhanced Image')

mean_Optimized = mean2(main_image)
var_optimzed = std2(main_image)

D = abs(uint8(main_image) - uint8(Im)).^2;
mse = sum(D(:))/numel(main_image)
psnr = 10*log10(255*255/mse)

%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%
summse = sum(cat(1,rmselev1,rmselev2,gmselev1,gmselev2,bmselev1,bmselev2));
sumpsnr = sum(cat(1,rpsnrlev1,rpsnrlev2,gpsnrlev1,gpsnrlev2,bpsnrlev1,bpsnrlev2));
tsize = size(summse);
t = 1:1:tsize(2);

figure,plot(t,summse);
xlabel('Number of Iterations');
ylabel('Mean Square Error');
title('MSE vs Iterations (Wavelet Domain)')

main_image_gray = rgb2gray(main_image);
grays = rgb2gray(Im);
 K = [0.6 0.6];
 [mssim, ssim_map] = ssim_index(grays,main_image_gray,K);
 mssim
 
 [FSIM, FSIMc] = FeatureSIM(grays,main_image_gray);
 FSIM

