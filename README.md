RF_Zernike.m calculates the first order spherical aberration (SA) at each 'z'.
The ideal case is automatically calculated, the non-ideal case is calculated depending on the magnification of the relay lenses.
The code requires the lens parameters for the imaging and reference objective (O1 and O2) and the relay lenses (L1 and L2).
The final output figure compares the SA at each z for the Ideal and Non-Ideal cases.

zernike_coeffs.m and zernfun.m are third-party MATLAB codes.
They are used for decomposing the wavefront into individual zernike modes. 
https://www.mathworks.com/matlabcentral/fileexchange/7687
https://www.mathworks.com/matlabcentral/fileexchange/27072

RF_Strehl calulates the decrease in dynamic range for each beta value. Final output figure plots dynamic range vs beta

#RF_System
