# NeuralProcessingTools
My Matlab code for various neural data processing/decoding functions in my PhD work

Contains:

* LEDSyncBoard: Code used to process video data recordings using Radu's IR LED board.
  * getLEDExtractionRegions.m : Get rectangles deliminating the ROI in the frames for each LED on the board.
  * extractLEDs.m : Given the ROIs for each LED, determine which LEDs are turned on for the given frame.
  * getFrameInds.m : Given which LED's are on in a frame, detemine how many strobes have elapsed since the beginning of the strobing.
  
* Decoding: Code use to do neural decoding
  * divideBlocks.m : Used to split a vector of numbers (trial indices) into even (or as even as possible) blocks for cross validation.
  * addHistory.m : Used to add a history component to matrices of decoder inputs/outputs data.
  * calcPerformanceMetrics.m : calculate decoding performance metrics such as R2, MSE, CC, ect.
  
* BlackrockProcessing: Various tools for analyszing/processing data from Blackrock Microsystems recording files.
  * artifactRejection.m : code used to preprocess .ns5 neural recordings from Utah arrays via the wireless Cereplex W system. Finds and removes antenna switching artifacts, rail jumps, and signal dropouts.
  * sortElectrodesNEV.m : Sorts the list of all the spike events in .nev files into electrode numbers.
  * getEncoderCounts.m : converts quadrature A/B/Z encoder signals into an encoder count.
  * getToneEvents.m : extract the times different tones are played in my treadmill rig.
  
* Misc: Other random stuff
  * rasterplot.m : make raster plot of neural spikes
