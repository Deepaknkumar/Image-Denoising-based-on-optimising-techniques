clear all
I=imread ('peppers_color.jpg');
r = I(:,:,1);
g = I(:,:,2);
b = I(:,:,3);
imagesize = 192;
rr = imresize(r,[imagesize imagesize]);
gr = imresize(g,[imagesize imagesize]);
br = imresize(b,[imagesize imagesize]);
Im(:,:,1) = rr;
Im(:,:,2) = gr;
Im(:,:,3) = br;
mrimage = [];
mgimage = [];
mbimage = [];

size_of_image = size(rr);        %gray scale image
figure,imshow(Im),title('Original Image')

nvar = .004;                  %variance of noise
ngrays = imnoise(Im,'gaussian',0,nvar);  %white gaussian noise
figure,imshow(ngrays),title('noised Image')
rngrays  = ngrays(:,:,1);
gngrays  = ngrays(:,:,2);
bngrays  = ngrays(:,:,3);

Numlev = 2;         %level of decomposition of wavelet domain

smatrix = 4;        %%% size of small matrix
numelem =  smatrix^2;  %%% Number of elements in a window 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RED PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rcoeffn,rbookem] = wavedec2(rngrays,Numlev,'db4');     %Red Pixels
[rcoeff,rbookem] = wavedec2(rr,Numlev,'db4');

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

rnestlev1 = [];
rmcoefflev1 = [];
tempcoeff1 = [];
n = 15;
rnumdec1 = ceil(levs(2)/(numelem));
rnumzero1 = (rnumdec1*numelem) - levs(2);
for i=1:4
    rcoefflev1(i,levs(2)+1:rnumdec1*numelem) = zeros(1,rnumzero1);
    rncoefflev1(i,levs(2)+1:rnumdec1*numelem) = zeros(1,rnumzero1);
    for j=1:rnumdec1
        tempnf = rcoefflev1(i,((j-1)*numelem)+1:j*numelem);
        tempn = rncoefflev1(i,((j-1)*numelem)+1:j*numelem);
        [rmcoefflev1(i,((j-1)*numelem)+1:j*numelem),rnestlev1(:,:,i)] = cuckoosc6(n,tempnf,tempn,nsd);
    end
    tempcoeff1(i,:) = rmcoefflev1(i,1:levs(2));
    fprintf('Level 2 Red coefficients matrix number %d processed\n',i)
end
rmcoefflev1 = tempcoeff1;

rnestlev2 = [];
rmcoefflev2 = [];
tempcoeff2 = [];
rnumdec2 = ceil(levs(3)/(numelem));
rnumzero2 = (rnumdec2*numelem) - levs(3);
for i=1:3
    rcoefflev2(i,levs(3)+1:rnumdec2*numelem) = zeros(1,rnumzero2);
    rncoefflev2(i,levs(3)+1:rnumdec2*numelem) = zeros(1,rnumzero2);
    for j=1:rnumdec2
        tempnf = rcoefflev2(i,((j-1)*numelem)+1:j*numelem);
        tempn = rncoefflev2(i,((j-1)*numelem)+1:j*numelem);
        [rmcoefflev2(i,((j-1)*numelem)+1:j*numelem),rnestlev2(:,:,i)] = cuckoosc6(n,tempnf,tempn,nsd);
    end
    tempcoeff2(i,:) = rmcoefflev2(i,1:levs(3));
    fprintf('Level 1 Red coefficients matrix number %d processed\n',i)
end
rmcoefflev2 = tempcoeff2;

rmcoeffmat = [];         %reconstruction of original signal
rmcoeffmat = cat(2,rmcoefflev1(1,:),rmcoefflev1(2,:),rmcoefflev1(3,:),rmcoefflev1(4,:),rmcoefflev2(1,:),rmcoefflev2(2,:),rmcoefflev2(3,:));
rimage = waverec2(rmcoeffmat,rbookem,'db4');



%%%%%%%%%%%%%%%%%%%%%%%% Green Part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[gcoeffn,gbookem] = wavedec2(gngrays,Numlev,'db4');     %Green pixels
[gcoeff,gbookem] = wavedec2(gr,Numlev,'db4');

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
gnestlev1 = [];
gmcoefflev1 = [];
tempcoeff1 = [];
n = 15;
gnumdec1 = ceil(levs(2)/(numelem));
gnumzero1 = (gnumdec1*numelem) - levs(2);
for i=1:4
    gcoefflev1(i,levs(2)+1:gnumdec1*numelem) = zeros(1,gnumzero1);
    gncoefflev1(i,levs(2)+1:gnumdec1*numelem) = zeros(1,gnumzero1);
    for j=1:gnumdec1
        tempnf = gcoefflev1(i,((j-1)*numelem)+1:j*numelem);
        tempn = gncoefflev1(i,((j-1)*numelem)+1:j*numelem);
        [gmcoefflev1(i,((j-1)*numelem)+1:j*numelem),gnestlev1(:,:,i)] = cuckoosc6(n,tempnf,tempn,nsd);
    end
    tempcoeff1(i,:) = gmcoefflev1(i,1:levs(2));
    fprintf('Level 2 Green coefficients matrix number %d processed\n',i)
end
gmcoefflev1 = tempcoeff1;

gnestlev2 = [];
gmcoefflev2 = [];
tempcoeff2 = [];
gnumdec2 = ceil(levs(3)/(numelem));
gnumzero2 = (gnumdec2*numelem) - levs(3);
for i=1:3
    gcoefflev2(i,levs(3)+1:gnumdec2*numelem) = zeros(1,gnumzero2);
    gncoefflev2(i,levs(3)+1:gnumdec2*numelem) = zeros(1,gnumzero2);
    for j=1:gnumdec2
        tempnf = gcoefflev2(i,((j-1)*numelem)+1:j*numelem);
        tempn = gncoefflev2(i,((j-1)*numelem)+1:j*numelem);
        [gmcoefflev2(i,((j-1)*numelem)+1:j*numelem),gnestlev2(:,:,i)] = cuckoosc6(n,tempnf,tempn,nsd);
    end
    tempcoeff2(i,:) = gmcoefflev2(i,1:levs(3));
    fprintf('Level 1 Green coefficients matrix number %d processed\n',i)
end
gmcoefflev2 = tempcoeff2;

gmcoeffmat = [];         %reconstruction of original signal
gmcoeffmat = cat(2,gmcoefflev1(1,:),gmcoefflev1(2,:),gmcoefflev1(3,:),gmcoefflev1(4,:),gmcoefflev2(1,:),gmcoefflev2(2,:),gmcoefflev2(3,:));
gimage = waverec2(gmcoeffmat,gbookem,'db4');



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Blue Part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[bcoeffn,bbookem] = wavedec2(bngrays,Numlev,'db4');     %Blue Pixels
[bcoeff,bbookem] = wavedec2(br,Numlev,'db4'); 

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
bnestlev1 = [];
bmcoefflev1 = [];
tempcoeff1 = [];
n = 15;
bnumdec1 = ceil(levs(2)/(numelem));
bnumzero1 = (bnumdec1*numelem) - levs(2);
for i=1:4
    bcoefflev1(i,levs(2)+1:bnumdec1*numelem) = zeros(1,bnumzero1);
    bncoefflev1(i,levs(2)+1:bnumdec1*numelem) = zeros(1,bnumzero1);
    for j=1:bnumdec1
        tempnf = bcoefflev1(i,((j-1)*numelem)+1:j*numelem);
        tempn = bncoefflev1(i,((j-1)*numelem)+1:j*numelem);
        [bmcoefflev1(i,((j-1)*numelem)+1:j*numelem),bnestlev1(:,:,i)] = cuckoosc6(n,tempnf,tempn,nsd);
    end
    tempcoeff1(i,:) = bmcoefflev1(i,1:levs(2));
    fprintf('Level 2 Blue coefficients matrix number %d processed\n',i)
end
bmcoefflev1 = tempcoeff1;

bnestlev2 = [];
bmcoefflev2 = [];
tempcoeff2 = [];
bnumdec2 = ceil(levs(3)/(numelem));
bnumzero2 = (bnumdec2*numelem) - levs(3);
for i=1:3
    bcoefflev2(i,levs(3)+1:bnumdec2*numelem) = zeros(1,bnumzero2);
    bncoefflev2(i,levs(3)+1:bnumdec2*numelem) = zeros(1,bnumzero2);
    for j=1:bnumdec2
        tempnf = bcoefflev2(i,((j-1)*numelem)+1:j*numelem);
        tempn = bncoefflev2(i,((j-1)*numelem)+1:j*numelem);
        [bmcoefflev2(i,((j-1)*numelem)+1:j*numelem),bnestlev2(:,:,i)] = cuckoosc6(n,tempnf,tempn,nsd);
    end
    tempcoeff2(i,:) = bmcoefflev2(i,1:levs(3));
    fprintf('Level 1 coefficients matrix number %d processed\n',i)
end
bmcoefflev2 = tempcoeff2;

bmcoeffmat = [];         %reconstruction of original signal
bmcoeffmat = cat(2,bmcoefflev1(1,:),bmcoefflev1(2,:),bmcoefflev1(3,:),bmcoefflev1(4,:),bmcoefflev2(1,:),bmcoefflev2(2,:),bmcoefflev2(3,:));
bimage = waverec2(bmcoeffmat,bbookem,'db4');

main_image = [];
main_image(:,:,1) = rimage;
main_image(:,:,2) = gimage;
main_image(:,:,3) = bimage;
main_image = uint8(main_image);
figure,imshow(main_image),title('Enhanced Image')

% rmodlev1 = [];
% rmodlev2 = [];
% gmodlev1 = [];
% gmodlev2 = [];
% bmodlev1 = [];
% bmodlev2 = [];

% %%%%%%%%%%%%%%%%%%%%%%% PSNR & MSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qsize = size(rnestlev1,1);
% for i=1:qsize
%                         %%% RED part %%%%
%     rmodlev1(1,:) = thresholdfuncsc4(rnestlev1(i,:,1),rncoefflev1(1,:));
%     rmodlev1(2,:) = thresholdfuncsc4(rnestlev1(i,:,2),rncoefflev1(2,:));
%     rmodlev1(3,:) = thresholdfuncsc4(rnestlev1(i,:,3),rncoefflev1(3,:));
%     rmodlev1(4,:) = thresholdfuncsc4(rnestlev1(i,:,4),rncoefflev1(4,:));
%     rmodlev2(1,:) = thresholdfuncsc4(rnestlev2(i,:,1),rncoefflev2(1,:));
%     rmodlev2(2,:) = thresholdfuncsc4(rnestlev2(i,:,2),rncoefflev2(2,:));
%     rmodlev2(3,:) = thresholdfuncsc4(rnestlev2(i,:,3),rncoefflev2(3,:));
%     
%     rmodmat = cat(2,rmodlev1(1,:),rmodlev1(2,:),rmodlev1(3,:),rmodlev1(4,:),rmodlev2(1,:),rmodlev2(2,:),rmodlev2(3,:));
%     rmodimage = waverec2(rmodmat,rbookem,'db4');
%     
%                     %%% Green part %%%
%     gmodlev1(1,:) = thresholdfuncsc4(gnestlev1(i,:,1),gncoefflev1(1,:));
%     gmodlev1(2,:) = thresholdfuncsc4(gnestlev1(i,:,2),gncoefflev1(2,:));
%     gmodlev1(3,:) = thresholdfuncsc4(gnestlev1(i,:,3),gncoefflev1(3,:));
%     gmodlev1(4,:) = thresholdfuncsc4(gnestlev1(i,:,4),gncoefflev1(4,:));
%     gmodlev2(1,:) = thresholdfuncsc4(gnestlev2(i,:,1),gncoefflev2(1,:));
%     gmodlev2(2,:) = thresholdfuncsc4(gnestlev2(i,:,2),gncoefflev2(2,:));
%     gmodlev2(3,:) = thresholdfuncsc4(gnestlev2(i,:,3),gncoefflev2(3,:));
%     
%     gmodmat = cat(2,gmodlev1(1,:),gmodlev1(2,:),gmodlev1(3,:),gmodlev1(4,:),gmodlev2(1,:),gmodlev2(2,:),gmodlev2(3,:));
%     gmodimage = waverec2(gmodmat,gbookem,'db4');
%     
%                     %%%% Blue Part %%%%
%     bmodlev1(1,:) = thresholdfuncsc4(bnestlev1(i,:,1),bncoefflev1(1,:));
%     bmodlev1(2,:) = thresholdfuncsc4(bnestlev1(i,:,2),bncoefflev1(2,:));
%     bmodlev1(3,:) = thresholdfuncsc4(bnestlev1(i,:,3),bncoefflev1(3,:));
%     bmodlev1(4,:) = thresholdfuncsc4(bnestlev1(i,:,4),bncoefflev1(4,:));
%     bmodlev2(1,:) = thresholdfuncsc4(bnestlev2(i,:,1),bncoefflev2(1,:));
%     bmodlev2(2,:) = thresholdfuncsc4(bnestlev2(i,:,2),bncoefflev2(2,:));
%     bmodlev2(3,:) = thresholdfuncsc4(bnestlev2(i,:,3),bncoefflev2(3,:));
%     
%     bmodmat = cat(2,bmodlev1(1,:),bmodlev1(2,:),bmodlev1(3,:),bmodlev1(4,:),bmodlev2(1,:),bmodlev2(2,:),bmodlev2(3,:));
%     bmodimage = waverec2(bmodmat,bbookem,'db4');
%     
%                 %%% Main Image %%%%
%     main_image(:,:,1) = rmodimage;
%     main_image(:,:,2) = gmodimage;
%     main_image(:,:,3) = bmodimage;
%     main_image = uint8(main_image);
%     
%     D = abs(uint8(main_image) - uint8(Im)).^2;
%     pmse(1,i) = sum(D(:))/numel(main_image);
%     ppsnr(1,i) = 10*log10(255*255/pmse(1,i));
% end

mean_Optimized = mean2(main_image)
var_optimzed = std2(main_image)

D = abs(uint8(main_image) - uint8(Im)).^2;
mse = sum(D(:))/numel(main_image)
psnr = 10*log10(255*255/mse)

% %%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%
% % summse = sum(cat(1,rmselev1,rmselev2,gmselev1,gmselev2,bmselev1,bmselev2));
% % sumpsnr = sum(cat(1,rpsnrlev1,rpsnrlev2,gpsnrlev1,gpsnrlev2,bpsnrlev1,bpsnrlev2));
% %tsize = size(summse);
% t = 1:1:qsize;
% 
% figure,plot(t,pmse);
% xlabel('Number of Iterations');
% ylabel('Mean Square Error');
% title('MSE vs Iterations')
% 
% figure,plot(t,ppsnr);
% xlabel('Number of Iterations');
% ylabel('Peak Signal To Noise Ratio');
% title('PSNR vs Iterations')
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

main_image_gray = rgb2gray(main_image);
grays = rgb2gray(Im);
 K = [0.6 0.6];
 [mssim, ssim_map] = ssim_index(grays,main_image_gray,K);
 mssim
 
 [FSIM, FSIMc] = FeatureSIM(grays,main_image_gray);
 FSIM

