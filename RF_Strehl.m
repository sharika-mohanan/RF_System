function [bta,DR_NId] = RF_Strehl(NA1, NA2, RI1, RI2, F1, F2_range, Mag1, Mag2, Depz, lambda, PS, IS)

%
% NA1, NA2: numerical aperture of imaging objective O1 and reference 
% objective O2
% RI1, RI2: Refractive index of objective O1 and O2
% F1: focal length of tube lens L1 (um)
% F2_range: range of focal lengths for L1 to give different beta values
% (um)
% Mag1, Mag2: Magnification of objectives
% Depz: distance from focal plane (um)
% lambda: emission wavelength (um)
% PS: pixel size (um)
% IS: image size
% bta: beta values, DR_NId: dynamic range calcualted for each beta
% value
% eg: [bta,DR_NId] = RF_Strehl(1.4, 0.95, 1.515, 1, 180000, 150000:1000:210000, 60, 40,...
% -100:0.25:100, 0.515, 5.3, 128)

bta = ones(numel(F2_range),1);
DR_NId = ones(numel(F2_range),1);

for kk = 1:numel(F2_range)

F2 = F2_range(kk); 
Mag1_Eff = Mag1*(F1/180000); Mag2_Eff = Mag2*(F2/180000);
f1 = F1/Mag1_Eff; 
f2 = F2/Mag2_Eff; 
M4f = F2/F1; M4f_Id = RI2*Mag1/(RI1*Mag2); 
bta(kk,1) = M4f/M4f_Id; % beta

sin_alpha1 = NA1/RI1; % maximum acceptance angle of O1
k0 = 2*pi/lambda; % wavenumber

% coordinate system
xc = (-IS/2):1:(IS/2)-1; % lateral coordinates
[xp, yp] = meshgrid(xc,xc);
rxy = sqrt(xp.^2+yp.^2); % radial coordinates

SF1 = Mag1_Eff / (PS * IS); % sampling frequency, pupil plane of O1

% sine, cos of the ray angles 
sine1 = (rxy .* SF1 .* lambda) ./ RI1; 
sine2 = sine1 * bta(kk); 
cos1 = real(sqrt(1-(sine1.^2)));
cos2 = real(sqrt(1-(sine2.^2)));

% find the limiting aperture of the system
% ideal system
rho1 = sine1./sin_alpha1;
pupil_mask1 = rho1<1; 

% non-ideal system
if (NA1/RI1)==(NA2/RI2)
if bta(kk)>1
    rho2 = sine2./sin_alpha1;
    pupil_mask2 = rho2<1; 
elseif bta(kk)<=1
    pupil_mask2 = pupil_mask1;
end
end

if (NA1/RI1)<(NA2/RI2)
if  bta(kk)>1 && sin_alpha1*bta(kk)<(NA2/RI2)
    clear rho2 pupil_mask2
    rho2 = rho1;
    pupil_mask2 = rho2<1; 
elseif bta(kk)>1 && sin_alpha1*bta(kk)>(NA2/RI2)
    clear rho2 pupil_mask2
    rho2 = sine2./(NA2/RI2);
    pupil_mask2 = rho2<1;
elseif bta(kk)<=1 
    pupil_mask2 = pupil_mask1;
end
end

Indx1 = pupil_mask1 ~= 0;  
Indx2 = pupil_mask2 ~= 0; 

SR_Id = ones(numel(Depz),1);
SR_NId = ones(numel(Depz),1);
maxInt_Id = ones(numel(Depz),1);
maxInt_NId = ones(numel(Depz),1);

for ii=1:numel(Depz)

    % Ideal Mapping phase (not required except to check dynamic range for ideal case)
    OPD_Id =  -((sine1.^2.*Depz(ii)^2)./f1); % OPD of RF system
    
    Dt = pupil_mask1.*((cos1)-(mean(mean(cos1(pupil_mask1==1))))); % Defocus calculation
    Df_CId = sum(sum(OPD_Id.*Dt))/sum(sum(Dt.^2)); 
    Df_Id = Df_CId.*Dt; 
    
    OPD_Id = OPD_Id - Df_Id; % removal of defocus
    
    Me1 = mean(OPD_Id(Indx1)); % mean
    Var1 = mean((OPD_Id(Indx1)-Me1).^2); % variance
    SR_Id(ii,1) = exp(-Var1.*k0^2); % strehl ratio
    
    clear Dt
    
    % Non-Ideal Mapping Phase

    OPD_NId1 = -(cos1.*Depz(ii))-((sine1.^2.*Depz(ii)^2)./(2.*f1)); % OPD, O1
    OPD_NId2 = (cos2.*Depz(ii))-((sine2.^2.*Depz(ii)^2)./(2.*f2)); %OPD, O2
    OPD_NId = OPD_NId1 + OPD_NId2;
    
    % defocus function for under and overmagnified cases
    if bta(kk)<=1
    Dt = pupil_mask2.*(cos2-(mean(mean(cos2(pupil_mask2==1)))));
    Df_CNId = sum(sum ((OPD_NId2+OPD_NId1).*Dt))/sum(sum(Dt.^2)); 
    Df_NId = Df_CNId.*Dt; 
    elseif bta(kk)>1
    Dt = pupil_mask2.*(cos1-(mean(mean(cos1(pupil_mask2==1)))));
    Df_CNId = sum(sum ((OPD_NId2+OPD_NId1).*Dt))/sum(sum(Dt.^2)); 
    Df_NId = Df_CNId.*Dt; 
    end
    
    OPD_NId = OPD_NId - Df_NId; % removal of defocus
    
    Me2 = mean(OPD_NId(Indx2)); % mean
    Var2 = mean((OPD_NId(Indx2)-Me2).^2); % variance
    SR_NId(ii,1) = exp(-Var2*k0^2); % strehl ratio

            pupilfn_Id = pupil_mask1.*exp(1i.*k0.*RI1.*OPD_Id); % pupil function   
            psf_Id = fftshift(fft2(ifftshift(pupilfn_Id))); % point spread function
            psf_Id = abs(psf_Id).^2; % PSF intensity

            pupilfn_NId = pupil_mask2.*exp(1i.*k0.*RI1.*OPD_NId); 
            psf_NId = fftshift(fft2(pupilfn_NId)); 
            psf_NId = abs(psf_NId).^2; 
                         
            maxInt_Id(ii,1) = max(psf_Id(:));
            maxInt_NId(ii,1) = max(psf_NId(:));   
  
        
   clear  psf_Ideal psf_NIdeal

end

% Strehl_Ratio_Id = maxInt_Id./max(maxInt_Id);
% DR_Id(kk,1) = numel(find(Strehl_Ratio_Id>0.8)).*abs(Depz(end)-Depz(end-1));
Strehl_Ratio_NId = maxInt_NId./max(maxInt_NId);
DR_NId(kk,1) = numel(find(Strehl_Ratio_NId>0.8)).*abs(Depz(end)-Depz(end-1));

% FDR_Id(kk,1) = numel(find(SR_Id>0.8)).*abs(Depz(end)-Depz(end-1));
% FDR_NId(kk,1) = numel(find(SR_NId>0.8)).*abs(Depz(end)-Depz(end-1));

% figure(1),plot(Depz,Strehl_Ratio_Id,'r',Depz,Strehl_Ratio_NId,'r*'),title('Strehl Ratio'),legend('Ideal','Non-Ideal')
% figure(2),plot(Depz,SR_Id,'r',Depz,SR_NId,'r*'),title('Strehl Ratio'),legend('Ideal','Non-Ideal')

disp(bta(kk))

end

figure(3),plot(bta,DR_NId,'r*'),title('Strehl Ratio'),legend('Ideal','Non-Ideal')

