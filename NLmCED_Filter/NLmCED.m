function DenoisedData = NLmCED(InputData,iter,rho,alpha,num)
% NLmCED_filter is an algorithm for reduction of Rician noise from 3D grayscale images.
% Function Interface:
%     DenoisedData = NLmCED_filter(InputData,iter,num)
% Input Arguments
%     InputData: Noisy data 
%     iter: Number of iterations
%     num : is used to distingush between sythetic data (==1) or invivo(==2)
% Output: 
%     DenoisedData: Final denoised 3D image 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 20xx-20xx   University of Al Manar Tunisia, ENIT and University of Turin (Unito), Italy
% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% AUTHORS:
% Feriel Romdhane - ferielromdhane@yahoo.fr
% Dario Longo - dario.longo@unito.it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[Y,X,Z]=size(InputData);
u = InputData;

%% Parameters 
% the following parameters can be modified 
radii = [num;num;num] ;  % gaussian windows
if (exist('num') ~= 1)
     radii = [1;1;1] ;  
end
rs= [1;1;1]; % search window
beta = 0.08;  
lambda_c = 1;

%% Create the gaussian windows for each direction:
gx = gausswin( 2*radii(1) + 1 ); gx = gx./sum(gx);
gy = gausswin( 2*radii(2) + 1 ); gy = gy./sum(gy);
gz = gausswin( 2*radii(3) + 1 ); gz = gz./sum(gz);

mu = My3DConv( u, gx, gy, gz );
out1 = zeros(Y,X,Z); 

for i=1:iter

%% Constants which are needed with CED eigenmode
[sig] = RicianSTD_NLMCED(double(u)); 
h=beta*sig;
%% create the gradiant 
usigma=imgaussian(u,sig,4*sig);
Gx=derivatives(usigma,'x');
Gy=derivatives(usigma,'y');
Gz=derivatives(usigma,'z');

%% calcutale the structure tensor 
[Jxx, Jxy, Jxz, Jyy, Jyz, Jzz]=StructureTensor3D(Gx,Gy,Gz, rho);
%% get the eigens and vectors values
[mu3,mu2,mu1,v3x,v3y,v3z,v2x,v2y,v2z,v1x,v1y,v1z]= EigenVectors3D(Jxx, Jxy, Jxz, Jyy, Jyz, Jzz);
%% 3D, Coherence Enhancing diffusion 
%% Calculate the parameter K 
K =(mu1-mu2).^2 + (mu1-mu3).^2 + (mu2-mu3).^2;
%% Edge indicator function Cedge
Cp = (mu1-mu2)./(mu1+ mu3); 
Cl = (mu2-mu3)./(mu2+ mu3); 
Cedge = Cl.*(1-Cp);

%% 3D, Proposed Coherence Enhancing diffusion with Tanh function 
lambda1 = alpha;
lambda2 = abs(tanh(Cedge./K));
lambda2 ((mu1<1e-15)&(mu1>-1e-15))=0;
lambda2((mu2<1e-15)&(mu2>-1e-15))=0;
lambda2 ((mu3<1e-15)&(mu3>-1e-15))=0;

lambda3 = alpha + (1.0-alpha)*exp(-0.6931*lambda_c^2./K);
lambda3 ((mu1<1e-15)&(mu1>-1e-15))=0;
lambda3((mu2<1e-15)&(mu2>-1e-15))=0;
lambda3 ((mu3<1e-15)&(mu3>-1e-15))=0;
 
%% components of the diffusion tensor D               
Dxx = lambda1.*v1x.^2   + lambda2.*v2x.^2   + lambda3.*v3x.^2;
Dyy = lambda1.*v1y.^2   + lambda2.*v2y.^2   + lambda3.*v3y.^2;
Dzz = lambda1.*v1z.^2   + lambda2.*v2z.^2   + lambda3.*v3z.^2;

Dxy = lambda1.*v1x.*v1y + lambda2.*v2x.*v2y + lambda3.*v3x.*v3y;
Dxz = lambda1.*v1x.*v1z + lambda2.*v2x.*v2z + lambda3.*v3x.*v3z;
Dyz = lambda1.*v1y.*v1z + lambda2.*v2y.*v2z + lambda3.*v3y.*v3z;

%% In literatue a,b,c,d,e and f are used as variables
a=Dxx; b=Dyy; c=Dzz; d=Dxy; e=Dxz; f=Dyz;               
%% Calculate positive and negative indices
[N1,N2,N3] = size(u);
N1=single(N1);
N2=single(N2);
N3=single(N3);

px = [2:N1,N1];
nx = [1,1:N1-1];

py = [2:N2,N2];
ny = [1,1:N2-1];
pz = [2:N3,N3];
nz = [1,1:N3-1];

% The Stencil Weights
A2 = 0.5*(-f(:,:,nz)-f(:,py,:));
A4 = 0.5*( e(:,:,nz)+e(nx,:,:));
A5 = ( c(:,:,nz)+c);
A6 = 0.5*(-e(:,:,nz)-e(px,:,:));
A8 = 0.5*( f(:,:,nz)+f(:,ny,:));

B1 = 0.5*(-d(nx,:,:) -d(:,py,:));
B2 = (b(:,py,:)+b);
B3 = 0.5*(d(px,:,:)+ d(:,py,:));
B4 = (a(nx,:,:)+a);
B5 = - (a(nx,:,:) + 2*a + a(px,:,:)) ...
      -(b(:,ny,:) + 2*b + b(:,py,:)) ...
      -(c(:,:,nz) + 2*c + c(:,:,pz));
B6 = (a(px,:,:)+a);
B7 = 0.5*(d(nx,:,:)+d(:,ny,:));
B8 = (b(:,ny,:)+b);
B9 = 0.5*(-d(px,:,:)-d(:,ny,:));

C2 = 0.5*(f(:,:,pz) + f(:,py,:));
C4 = 0.5*(-e(:,:,pz)-e(nx,:,:));
C5 = (c(:,:,pz)+c);
C6 = 0.5*(e(:,:,pz)+e(px,:,:));
C8 = 0.5*(-f(:,:,pz)-f(:,ny,:));

     I= A2.*(u(: ,py,nz )-u) + ...
        A4.*(u(nx ,: ,nz)-u) + ...
        A5.*(u(: ,: ,nz)-u)  + ... 
        A6.*(u(px,: ,nz)-u)  + ... 
        A8.*(u(: ,ny,nz)-u)  + ...
        B1.*(u(nx,py,: )-u)  + ...
        B2.*(u(: ,py,: )-u)  + ...
        B3.*(u(px,py,: )-u)  + ...
        B4.*(u(nx,: ,: )-u)  + ...
        B5.*(u(:,: ,: )-u)   + ...
        B6.*(u(px,: ,: )-u)  + ...
        B7.*(u(nx,ny,: )-u)  + ...
        B8.*(u(: ,ny,: )-u)  + ...
        B9.*(u(px,ny,: )-u)  + ...
        C2.*(u(: ,py,pz)-u)  + ...
        C4.*(u(nx,: ,pz)-u)  + ...
        C5.*(u(: ,: ,pz)-u)  + ...
        C6.*(u(px,: ,pz)-u)  + ...
        C8.*(u(: ,ny,pz)-u);
    

%% Loop along the pixels: 
for x=1:X
    for y=1:Y
        for z=1:Z 
                    %%% We are filtering the pixel (x,y,z). First, create a
                    %%% neighborhood around this pixel checking for out-of-bound
                    %%% indices:
                    mx = max( x-rs(1), 1 );
                    MX = min( x+rs(1), X );
                    my = max( y-rs(2), 1 );
                    MY = min( y+rs(2), Y );
                    mz = max( z-rs(3), 1 );
                    MZ = min( z+rs(3), Z );
                    %%% Keep the center values:
                    mu0 = mu(y,x,z);
                    I0 = I(y,x,z);
                    
                    %%% Get the values of the pixels in the whole search neihborhood:
                    vals  = u(my:MY,mx:MX,mz:MZ);
                    %%%  Get the mean values and gradients of the pixels in the whole
                    %%%  search neighborhood:
                    mui   = mu(my:MY,mx:MX,mz:MZ);
                    Ii   = I(my:MY,mx:MX,mz:MZ);
                    
                    % compute the distance
                    dists =  (mui-mu0).*(mui-mu0)+ (Ii -I0).*(Ii -I0);
                    % Normalize the distances:
                    dists = dists.*dists./(h*h);
                    %%Compute the weights:
                    wis   = exp(-dists);
                    %Compute the normalization factor:
                    NORM  = sum(wis(:));
                    % Filter the pixel;
                    pixel = sum(wis(:).*vals(:));
                    % Normalize the pixel
                    if  NORM > 0
                        out1(y,x,z) = pixel/NORM;
                    else
                        out1(y,x,z) = u(y,x,z);
                    end
                end
           
        end
    end

u=out1; 
end
DenoisedData=u;