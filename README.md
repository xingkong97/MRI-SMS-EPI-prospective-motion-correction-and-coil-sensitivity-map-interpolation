# MRI SMS-EPI Prospective motion correction by real-time receiver phase correction and coil sensitivity map interpolation
Movements change the location of the subject and alter coil sensitivity encoding for imaging. The corrections of receiver phase and coil sensitivity maps are required for SMS reconstruction. 

Please read the paper for more details. "Simultaneous multislice EPI prospective motion correction by realtime receiver phase correction and coil sensitivity map interpolation. Magn Reson Med. 2023;1-17. doi:10.1002/mrm.29789"

Since the receiver phase correction has been done during acquisition, the demo shows that the correction of coil sensitivity maps by interpolation can substantially attenuate image artifacts that caused by coil sensitivity mismatch between the reference and updated slice locations. 

Run demo.m script in MATLAB.

Code consists of the following parts:   

1. Coil sensitivity maps (CSMs) interpolation. 
  getCoilSen.m is to obtail coil sensitivity maps by ESPIRIT method. 
  getUpdatedCoilSen.m is to interpolate CSMs. 
  
3. SENSE recon with original and updated CSMs.
DoSENSE_SMS.m is to reconstruct SMS aliased images by SENSE method. 

4. Split slice-GRAPPA recon with original and updated CSMs. 
DoSplitSliceGrappa.m is to reconstruct SMS aliased images by split slice-GRAPPA method
