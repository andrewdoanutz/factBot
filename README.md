# factBot

A bot allowing the user to say a voice command (currently only "learn 1", "learn 2", "learn 3") and get a random fact back. The program is implemented in MATLAB using spectral flux and Mel Scale Cepstral Coefficients(MFCC) as the features of the audio input. These features are then classified by a k-Nearest Neighbors (10 neighbors) and a Gaussian classifer in different versions of the program. Two different classifiers were used to test for accuracy. The datasets used to the train the program were the Free Spoken Digit Dataset and Speech Commands Dataset.

## Instructions

1. Add the training folder as a subfolder in MATLAB
2. Run the program and hit 1 to record a command, or 2 to use a prerecorded command file
3. The program will either load the training sets, or attempt to grab the .wav files in the training folder to create a model for spectral flux and a model for MFCC
4. If successful recognized, the program will output a fact in the MATLAB console

