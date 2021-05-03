# MRI-SMS-EPI-Reconstruction-with-motion-correction
Movements change the location of the subject in a coil and alter coil sensitivity encoding for imaging. So the correction of coil sensitivity maps is required for SMS reconstruction. 
This code is to recon SMS with coil sensitivity maps interpolation. 
Please let me know if you have any questions. 

Code consists of the following parts:   

1. Load SMS acquisition k-space data and imaging protocol parameters. Motion parameters are also loaded. 
The main file of the project is SMS_recon_with_motion.m. 

2. Coil sensitivity maps (CSMs) are generated with phase alignment. 
CSMs were calculated by eigenvector-based parallel imaging technique (ESPIRiT). We used the set of ESPIRiT estimated coil sensitivity maps with the largest eigenvalues as CSMs. We modified the kernelEig.m to get phase aligned for each channel. 
CalCoilSens.m is to calculate the coil sensitivity maps. 

3. SENSE recon with original and updated CSMs.
DoSENSE_SMS.m is to reconstruct SMS aliased images by SENSE method. 

4. Split slice-GRAPPA recon with original and updated CSMs. 
DoSplitSliceGrappa.m is to reconstruct SMS aliased images by split slice-GRAPPA method
