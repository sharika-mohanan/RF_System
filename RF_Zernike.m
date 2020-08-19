function [SAc_z] = RF_Zernike(NA1, NA2, RI1, RI2, F1, F2, Mag1, Mag2, Depz, lambda, PS, IS)
%
% NA1, NA2: numerical aperture of imaging objective O1 and reference 
% objective O2
% RI1, RI2: Refractive index of objective O1 and O2
% F1, F2: focal length of tube lens (um)
% Mag1, Mag2: Magnification of objectives
% Depz: distance from focal plane (um)
% lambda: emission wavelength (um)
% PS: pixel size (um)
% IS: image size
% SAc_z: first order spherical aberration coeff, 1st column: ideal, 2nd
% column: non-ideal
% eg: SAc_z = RF_Zernike(1.15, 0.95, 1.33, 1, 180000, 135000, 40, 40,...
% -200:5:200, 0.532, 5.3, 1280)
% to be used along with zernike_coeffs.m and zernfun.m

Mag1_Eff = Mag1*(F1/180000); Mag2_Eff = Mag2*(F2/180000); % Effective magnification of the objectives
f1 = F1/Mag1_Eff; f2 = F2/Mag2_Eff;% focal length of objectives
M4f = F2/F1; M4f_Id = RI2*Mag1/(RI1*Mag2); bta = M4f/M4f_Id; % non-ideal, ideal and beta factor


sin_alpha1 = NA1/RI1; % maximum acceptance angle of O1

% coordinate system
xc = (-IS/2):1:(IS/2)-1; % lateral coordinates
[xp, yp] = meshgrid(xc,xc);
rxy = sqrt(xp.^2+yp.^2); % radial coordinates

SF1 = Mag1_Eff / (PS * IS); % sampling frequency, pupil plane of O1

% sine, cos of the ray angles 
sine1 = (rxy .* SF1 .* lambda) ./ RI1; 
sine2 = sine1 * bta; 
cos1 = sqrt(1-(sine1.^2));
cos2 = sqrt(1-(sine2.^2));

%% calculate the limiting aperture of the RF system
% ideal system
rho1 = sine1./sin_alpha1;
pupil_mask_Id = rho1<1; 

% non-ideal system
% for objective with same angular aperture

if (NA1/RI1)==(NA2/RI2)
if bta>1 
    rho2 = sine2./sin_alpha1;
    pupil_mask_NId = rho2<1; 
elseif bta<=1
    pupil_mask_NId = pupil_mask_Id;
end
end

% for objectives where sine(alpha2)>sine(alpha1)

if (NA2/RI2)>(NA1/RI1) 
    if bta>1 && sin_alpha1*bta>(NA2/RI2)
        rho2 = sine2./(NA2/RI2);
        pupil_mask_NId = rho2<1;
    elseif  bta>1 && sin_alpha1*bta<(NA2/RI2) 
        pupil_mask_NId = pupil_mask_Id; 
    elseif bta<=1 
        pupil_mask_NId = pupil_mask_Id;
    end
end

[row1,col1] = find(pupil_mask_Id);
[row2,col2] = find(pupil_mask_NId);

%%
ZC_Ideal = zeros(25,numel(Depz));
ZC_NIdeal = zeros(25,numel(Depz));

for ii=1:numel(Depz)

    % Ideal Mapping wavefront
    
    OPD_Id =  -((sine1.^2.*Depz(ii)^2)./f1); % OPD of wavefront for an ideal RF system
    
    Df = pupil_mask_Id.*((cos1)-(mean(mean(cos1(pupil_mask_Id==1))))); % Defocus mode
    Df_CId = sum(sum(OPD_Id.*Df))/sum(sum(Df.^2)); % Defocus coefficient
    Df_Id = Df_CId.*Df; 
    
    OPD_Id = OPD_Id - Df_Id; % removal of defocus term

    ZC_Ideal(:,ii) = zernike_coeffs(OPD_Id(min(row1):max(row1),min(col1):max(col1)), 25, 256); % calculation of zernike coefficients 
    
    % Non-Ideal Mapping wavefront

    OPD_NId1 = -(cos1.*Depz(ii))-((sine1.^2.*Depz(ii)^2)./(2.*f1)); % OPD, O1
    OPD_NId2 = (cos2.*Depz(ii))-((sine2.^2.*Depz(ii)^2)./(2.*f2)); %OPD, O2
    OPD_NId = OPD_NId1 + OPD_NId2;
    
    % defocus function for under and overmagnified cases
    if bta<=1
    Df = pupil_mask_NId.*(cos2-(mean(mean(cos2(pupil_mask_NId==1)))));
    Df_CNId = sum(sum ((OPD_NId2+OPD_NId1).*Df))/sum(sum(Df.^2)); 
    Df_NId = Df_CNId.*Df; 
    elseif bta>1
    Df = pupil_mask_NId.*(cos1-(mean(mean(cos1(pupil_mask_NId==1)))));
    Df_CNId = sum(sum ((OPD_NId2+OPD_NId1).*Df))/sum(sum(Df.^2)); 
    Df_NId = Df_CNId.*Df; 
    end
    
    OPD_NId = OPD_NId - Df_NId; 
    
    ZC_NIdeal(:,ii) = zernike_coeffs(OPD_NId(min(row2):max(row2),min(col2):max(col2)), 25, 256); % calculation of zernike coefficient
 
   disp(Depz(ii))

end

SAc_z(:,1) =  real(ZC_Ideal(11,:)); % ideal config. first order spherical aberration
SAc_z(:,2) =  real(ZC_NIdeal(11,:)); % non-ideal config. first order spherical aberration

% plot final results
figure(1),plot(Depz,(ZC_Ideal(11,:)),'r',Depz,(ZC_NIdeal(11,:)),'r*'),title('First Order Spherical'),legend('Ideal','Non-Ideal')



