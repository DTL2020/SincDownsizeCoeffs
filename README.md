# SincDownsizeCoeffs
Supplementary tool for UserDefined2Resize (b,c,s) coeffs and param computing for different use cases

Usage: Run in command line and enter mode 1 to 4 and enter downsize ratio and estimated or better measured with waveform monitor or some 'pixelprofile' software overshoot amplitude in the input source in % (do not forgot to upsample it with sincresize to 2x size or more to see real overshoot at display scaler or other interpolation software).

Example: The downscale of some sharp content to 0.1 of initial size with UserDefined2Resize(width/10, height/10, b=80, c=-20) creates about 12% overshoot.
