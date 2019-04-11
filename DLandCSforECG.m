

clear all;clc;close all;
cd 'D:\Backup Drivers\Local Disk D\about cs\matlab\20xx papers\source'
itr=10;
SNR_DCTT=0;SNR_MODD=0; SNR_KSVDD=0;
for i= 1:itr
data0= load('104m.mat');%contains    val: [2x21600 double]
train = data0.val(1,1:400*128);
test = data0.val(1,400*128+1:425*128);

l=400;%Number of trained signal
n=128;%length of signal
Data = zeros(n,l);


for i=1:l
    TrainMat(:,i) = train((i-1)*n+1:i*n);
    %plot(Data(:,i));pause(1)
end

%-------------dictionary learning
param.K = 2*n;% number of atom in dictionary
param.L = 1;param.numIteration = 10;
param.errorFlag = 0;param.preserveDCAtom =0;
param.displayProgress = 0;param.InitializationMethod = 'DataElements';
param.TrueDictionary = randn(n,2*n);%param.InitializationMethod = 'DataElements';
%param.InitializationMethod = 'GivenMatrix';
iniMat = randn(n,param.K);
for i =1: param.K
    iniMat(:,i) = iniMat(:,i)/norm(iniMat(:,i));%normalizing columns of matrix
end
param.initialDictionary = iniMat;
[DicMod, outputMod] = MOD(TrainMat,param);
%save('Sparsifying_ECG_128_256','DicMod');
[DicKSVD,X] = KSVD(TrainMat,param);

%-------------------------------
%-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

x=test';     

%------------------------measurement matrix

N=128;M=N/4;
%A=GenerateMatrix(M,N);
%A=randn(M,N);
A=ones(M,N);A=binornd(A,.5);A=A-.5;A=1/sqrt(M)*A;
dict_DCT = wmpdictionary(N,'LstCpt',{'dct'});
dict_MOD = DicMod;
dict_KSVD = DicKSVD;

%dict1 = wmpdictionary(N1,'lstcpt',{{'Haar',5}});
%load Sparsifying_ECG_512_1024;dict1 = DicMod;

A1_DCT=A*dict_DCT;
A1_MOD=A*dict_MOD;
A1_KSVD=A*dict_KSVD;


%--------------------Sl0 Parameters
for i=1:1
    sigma_off = 0.001;
   A_pinv_DCT = pinv(A1_DCT); A_pinv_MOD = pinv(A1_MOD);
   A_pinv_KSVD = pinv(A1_KSVD);mu_0 = 2;sigma_decrease_factor = 0.5;L = 3;
    %true_s = sparseSigGen4plusNoise(9,floor(27/4),sigma_off);
    if sigma_off>0
        sigma_min = sigma_off*4;
    else
        sigma_min = 0.00001;
    end
end
%------------------------New method Kroneckered

for i=1:length(test)/N
        j=i;
        y=A*x((i-1)*N+1:N*i,1);
        
        xp_DCT = SL0(A1_DCT, y, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv_DCT);
        xp_MOD = SL0(A1_MOD, y, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv_MOD);
        xp_KSVD = SL0(A1_KSVD, y, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv_KSVD);

        zm_DCT=dict_DCT*xp_DCT; 
        zz_DCT(N*j-(N-1):N*j)=zm_DCT(:);
        
        zm_MOD=dict_MOD*xp_MOD; 
        zz_MOD(N*j-(N-1):N*j)=zm_MOD(:);
        
        zm_KSVD=dict_KSVD*xp_KSVD; 
        zz_KSVD(N*j-(N-1):N*j)=zm_KSVD(:);
end

err_DCT = zz_DCT-test;SNR_DCT = 20*log10(norm(test)/norm(err_DCT));
err_MOD = zz_MOD-test;SNR_MOD = 20*log10(norm(test)/norm(err_MOD));
err_KSVD = zz_KSVD-test;SNR_KSVD = 20*log10(norm(test)/norm(err_KSVD));

    
% plot(test);hold on;plot(zz_DCT,'r'); 
% figure;plot(test);hold on;plot(zz_MOD,'r');
% figure;plot(test);hold on;plot(zz_KSVD,'y');

SNR_DCTT=SNR_DCTT+SNR_DCT;
SNR_MODD=SNR_MODD+SNR_MOD;
SNR_KSVDD=SNR_KSVDD+SNR_KSVD;

end
[SNR_DCTT SNR_MODD SNR_KSVDD]/itr



