clear all
I=imread ('house.jpg');
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

nvar = .004;                                      %variance of noise
ngrays = imnoise(grays,'gaussian',0,nvar);        %white gaussian noise
figure,imshow(ngrays),title('noised Image')
Numlev = 2;         %level of decomposition of wavelet domain
smatrix = 4;        %%% size of small matrix
num_decomposition = (size_of_image(1)*size_of_image(2))/(smatrix.^2); 
noisedc_image = ngrays;
original_image = grays;
ngrays = [];
grays = [];
for i=1:(size_of_image(1)/smatrix)          %%%% smatrix*smatrix of noisy image
    for j=1:(size_of_image(2)/smatrix)
     ngrays(:,:,((i-1)*(size_of_image(1)/smatrix))+j) = noisedc_image(((i-1)*smatrix)+1:i*smatrix,((j-1)*smatrix)+1:j*smatrix);
    end
end

for i=1:(size_of_image(1)/smatrix)          %%%% smatrix*smatrix of original image
    for j=1:(size_of_image(2)/smatrix)
     grays(:,:,((i-1)*(size_of_image(1)/smatrix))+j) = original_image(((i-1)*smatrix)+1:i*smatrix,((j-1)*smatrix)+1:j*smatrix);
    end
end

for i1=1:num_decomposition
[coeffn,bookem] = wavedec2(ngrays(:,:,i1),Numlev,'db4');
[coeff,bookem] = wavedec2(grays(:,:,i1),Numlev,'db4');
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
mcoefflev1 = [];

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
maxIter = 10;

for i=1:4
    tempnf = coefflev1(i,:);
    tempn = ncoefflev1(i,:);
    [mcoefflev1(i,:),mselev1(i,:),psnrlev1(i,:)] = dssc6(method,@add,[],n,dimensions,lower_limit,upper_limit,maxIter,tempnf,tempn,nsd);
end

mcoefflev2 = [];
for i=1:3
    tempnf = coefflev2(i,:);
    tempn = ncoefflev2(i,:);
    [mcoefflev2(i,:),mselev2(i,:),psnrlev2(i,:)] = dssc6(method,@add,[],n,dimensions,lower_limit,upper_limit,maxIter,tempnf,tempn,nsd);
end

mcoeffmat = [];         %reconstruction of original signal
mcoeffmat = cat(2,mcoefflev1(1,:),mcoefflev1(2,:),mcoefflev1(3,:),mcoefflev1(4,:),mcoefflev2(1,:),mcoefflev2(2,:),mcoefflev2(3,:));
mimage(:,:,i1) = waverec2(mcoeffmat,bookem,'db4');
i1
end

tmpimage = [];
for i=1:(size_of_image(1)/smatrix)          %%%% smatrix*smatrix of noisy image
    for j=1:(size_of_image(2)/smatrix)
     tmpimage(((i-1)*smatrix)+1:i*smatrix,((j-1)*smatrix)+1:j*smatrix) = mimage(:,:,((i-1)*(size_of_image(1)/smatrix))+j);
    end
end
mimage = tmpimage;
grays = original_image;
main_image = uint8(mimage);

figure,imshow(main_image),title('Enhanced Image')

mean_Optimized = mean2(main_image)
var_optimzed = std2(main_image)

D = abs(uint8(main_image) - uint8(grays)).^2;
mse = sum(D(:))/numel(main_image)
psnr = 10*log10(255*255/mse)

% %%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%
% summse = sum(cat(1,mselev1,mselev2));
% sumpsnr = sum(cat(1,psnrlev1,psnrlev2));
% tsize = size(summse);
% t = 1:1:tsize(2);
% 
% figure,plot(t,summse);
% xlabel('Number of Iterations');
% ylabel('Mean Square Error');
% title('MSE vs Iterations (Wavelet Domain)')
% 
% % figure,plot(t,sumpsnr);
% % xlabel('Number of Iterations');
% % ylabel('Peak Signal To Noise Ratio');
% % title('PSNR vs Iterations (Wavelet Domain)')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

