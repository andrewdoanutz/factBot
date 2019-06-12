
function [results,acc] = classifyWord(x,fs)
    if isfile('training_data1.mat') && isfile('training_data2.mat')
          mfccModel = load('training_data.mat');
           specFluxModel = load('training_data2.mat');
           [mfccRes,mfccAcc] =mfccPredict(mfccModel,x,fs);
           [fluxRes,fluxAcc] =fluxPredict(specFluxModel,x,fs);
           mfccAcc=max(mfccAcc)/sum(mfccAcc);
           fluxAcc=max(fluxAcc)/sum(fluxAcc);
           if(mfccAcc>=fluxAcc)
               results = mfccRes;
               acc=mfccAcc*100;
           else
               acc=fluxAcc*100;
               results=fluxRes;
           end
    else
        %do the training
         mfccModel = train();
         specFluxModel = train2();
         [mfccRes,mfccAcc] =mfccPredict(mfccModel,x,fs);
           [fluxRes,fluxAcc] =fluxPredict(specFluxModel,x,fs);
           %mfccAcc=max(mfccAcc)/sum(mfccAcc);
           %fluxAcc=max(fluxAcc)/sum(fluxAcc);
           if(mfccAcc>=fluxAcc)
               results = mfccRes;
               acc=mfccAcc*100;
           else
               acc=fluxAcc*100;
               results=fluxRes;
           end
    end
end

function mfccModel = train()
    training = [];
    label = [];
    type = '';
    filenames=dir('training/*.wav');
    first=1;;
    for i=1:length(filenames)
        [x,fs]=audioread(filenames(i).name);
        if(strcmp(filenames(i).name(1:1),'1'))
            type='1';
        elseif(strcmp(filenames(i).name(1:1),'2'))
            type='2';
        elseif(strcmp(filenames(i).name(1:1),'3'))
            type='3';
        elseif(strcmp(filenames(i).name(1:1),'4'))
            type='4';
        elseif(strcmp(filenames(i).name(1:1),'5'))
            type='5';
        elseif(strcmp(filenames(i).name(1:1),'6'))
            type='6';
        elseif(strcmp(filenames(i).name(1:1),'7'))
            type='7';
        elseif(strcmp(filenames(i).name(1:1),'8'))
            type='8';
        elseif(strcmp(filenames(i).name(1:1),'9'))
            type='9';
        elseif(strcmp(filenames(i).name(1:1),'l'))
            type='0';
        else
            continue;
        end
        mfcc=genmfcc(x,fs,20,7,20,0.020);
        mfcc = reshape(mfcc,1,[]);
        label = [label ;type];
        training=[training ; mfcc];
    end
    mfccModel = fitcknn(training,label,'NumNeighbors',10,'Distance','hamming');
    save('training_data.mat','mfccModel');
end

function specFluxModel = train2()
    training2 = [];
    label2 = [];
    type = '';
    filenames=dir('training/*.wav');
    first=1;
    for i=1:length(filenames)
        [x,fs]=audioread(filenames(i).name);
        audio = x;
        if(strcmp(filenames(i).name(1:1),'1'))
            type='1';
        elseif(strcmp(filenames(i).name(1:1),'2'))
            type='2';
        elseif(strcmp(filenames(i).name(1:1),'3'))
            type='3';
        elseif(strcmp(filenames(i).name(1:1),'4'))
            type='4';
        elseif(strcmp(filenames(i).name(1:1),'5'))
            type='5';
        elseif(strcmp(filenames(i).name(1:1),'6'))
            type='6';
        elseif(strcmp(filenames(i).name(1:1),'7'))
            type='7';
        elseif(strcmp(filenames(i).name(1:1),'8'))
            type='8';
        elseif(strcmp(filenames(i).name(1:1),'9'))
            type='9';
        elseif(strcmp(filenames(i).name(1:1),'l'))
            type='0';
        else
            continue;
        end
        Len=24*fs/1000;   %Get Length for Finding Spectral Flux
        hopLen=Len/(3/2); %How much hops each time
        numBlocks=floor(length(audio)/hopLen)-3; %find number of blocks required
        Tblock=((0:numBlocks-1)*hopLen+(Len/2))/fs; 
        sFlux=zeros(numBlocks,1); 
        OldDFT=0; %to hold the previous dft block
        iHopLength = 256;  %get hop length for finding energy
        L = 4096; %get length of each block
        iNumBlocks = floor(length(audio)/iHopLength)-3; %find number of blocks required
        Tblock2 = ((0:iNumBlocks-1)*iHopLength + (L/2))/fs;
        energy = zeros(iNumBlocks,1);
        for(n=1:numBlocks)
            start=(n-1)*hopLen+1;
            stop=min(length(audio),start+Len-1);
            templen=stop-start+1;
            window=hamming(templen);
            fftsize=2*templen;
            NewDFT=fft(audio(start:stop).*window,fftsize);
            rms_specflux(n)=sqrt(sum((abs(NewDFT)-abs(OldDFT)).^2)); %calculate the spectral flux
            OldDFT=NewDFT;
        end
        traincov = cov(rms_specflux);
        
        training2=[training2;traincov];
        label2 = [label2;type];
    end
     specFluxModel = fitcknn(training2,label2,'NumNeighbors',10,'Distance','hamming');
     save('training_data2.mat','specFluxModel');
end

function [fluxGenre,score] = fluxPredict(model,s,fs)
    audio = s;
    Len=24*fs/1000;   %Get Length for Finding Spectral Flux
    hopLen=Len/(3/2); %How much hops each time
    numBlocks=floor(length(audio)/hopLen)-3; %find number of blocks required
    Tblock=((0:numBlocks-1)*hopLen+(Len/2))/fs; 
    sFlux=zeros(numBlocks,1); 
    OldDFT=0; %to hold the previous dft block
     rms_specflux = [];
    iHopLength = 256;  %get hop length for finding energy
    L = 4096; %get length of each block
    iNumBlocks = floor(length(audio)/iHopLength)-3; %find number of blocks required
    Tblock2 = ((0:iNumBlocks-1)*iHopLength + (L/2))/fs;
    energy = zeros(iNumBlocks,1);
    
    for(n=1:numBlocks)
        start=(n-1)*hopLen+1;
        stop=min(length(audio),start+Len-1);
        templen=stop-start+1;
        window=hamming(templen);
        fftsize=2*templen;
        NewDFT=fft(audio(start:stop).*window,fftsize);
        rms_specflux(n)=sqrt(sum((abs(NewDFT)-abs(OldDFT)).^2)); %calculate the spectral flux
        OldDFT=NewDFT;
    end
   testVal = cov(rms_specflux);
   [class,score] = predict(model,testVal);
   display(["SpecFlux Predicts" , class]);
   fluxGenre = class;
end



function [mfccGenre,score] = mfccPredict(model,s,fs)
    testVal = genmfcc(s,fs,20,7,20,0.020);
    testVal = reshape(testVal,1,[]);
    [class,score] = predict(model,testVal);
    display(["MFCC Predicts", class]);
    mfccGenre = class;
end

function MFCC_coeff = genmfcc(s,fs, binCount, frameCount, coeffCount, stepTime)

    % song_file_name is a string specifying a .mp3 or .au audio file.
    % MFCC_MATRIX is a matrix of order of No. of Frames X Number of Mel
    % Cofficients.
    % COUNT_FRAMES is how many frames (20 ms samples) of the song to use.
    % COUNT_BINS is how many mel frequency bins to map each frame spectrum to.
    % COUNT_COEFF is how many mel coefficients to keep (COUNT_COEFF <= COUNT_BINS).
    % STEP_TIME is the frequency after which frame is captured.
    %
           seconds = 1;  % take a full 3-second sample (then split into 20ms frames)
    %offset = floor(length(s)/1);
    %s=s(offset:offset+fs*seconds);
    fftlen = 128;              % length of fft
    timeFrame = 0.020;         % # seconds/frame
    frameLength = floor(fs*timeFrame);  % # samples/frame
    % STEP_TIME = 0.010;          % # seconds/frame step
    stepLength = floor(fs*stepTime);    % # samples/frame step

    MFCC_coeff = zeros(frameCount, coeffCount);

    window = hamming(frameLength);
    stop = 1+floor(fftlen/2);  

    %[x, ~] = size(s);
    x=length(s);
    begin = 0;
 
    for i = 1:frameCount
        firstIndex = begin + (i-1)*stepLength + 1;
        lastIndex = begin + (i-1)*stepLength + frameLength;
        if(lastIndex>x)
            s=[s;zeros(lastIndex-x,1)];
        end
        %fprintf("first: %d \n last: %d\n",firstIndex,lastIndex);
        f = s(firstIndex:lastIndex);
        heightF = size(f, 1);
        widthF = size(f, 2);
        if (heightF > widthF)
            f = f.*window;
        else
            f = f'.*window;
        end
        f = fft(f, frameLength);
        mel_bins = mel_scale(binCount, fftlen, fs);
        f = abs(real(f(1:stop))).^2; 
        f=mel_bins*f;
        f(f==0)=1;
        m = log10(f);
        m = dct(m);          % dim(m): COUNT_BINS x 1
        m = m(1:coeffCount);
        MFCC_coeff(i, :) = m'; 
    end
end

function [xmatrix,lowestBin,highestBin] = mel_scale(numOfFilters,fftlen,fs,fMin,fMax,w)
    %  
    %  Determine matrix for a mel-spaced filterbank 
    %
    %  Inputs:   p   number of filters in filterbank
    %            n   length of fft
    %            fs  sample rate in Hz
    %            fl  low end of the lowest filter as a fraction of fs (default = 0)
    %            fh  high end of highest filter as a fraction of fs (default = 0.5)
    %            w   any sensible combination of the following:
    %                't'  triangular shaped filters in mel domain (default)
    %                'n'  hanning shaped filters in mel domain
    %                'm'  hamming shaped filters in mel domain
    %
    %                'z'  highest and lowest filters taper down to zero (default)
    %                'y'  lowest filter remains at 1 down to 0 frequency and
    %                     highest m remains at 1 up to nyquist freqency
    %
    %                If 'ty' or 'ny' is specified, the total power in the FFT 
    %                is preserved.
    %
    %  Outputs:  x   sparse matrix containing the filterbank amplitudes
    %                If x is the only output argument then 
    %                  size(x)=[p,1+floor(n/2)]
    %                otherwise 
    %                  size(x)=[p,mx-mn+1]
    %            mn  the lowest fft bin with a non-zero coefficient
    %            mx  the highest fft bin with a non-zero coefficient

    if nargin < 6
      w='tz';
      if nargin < 5
        fMax=0.5;
        if nargin < 4
          fMin=0;
        end
      end
    end
    f0=700/fs;
    %f0=1000/fs;
    fn2=floor(fftlen/2);
    lr=log((f0+fMax)/(f0+fMin))/(numOfFilters+1);
    % convert to fft bin numbers with 0 for DC term
    binl=fftlen*((f0+fMin)*exp([0 1 numOfFilters numOfFilters+1]*lr)-f0);
    bin2=ceil(binl(2));
    bin3=floor(binl(3));
    if any(w=='y')
      pf=log((f0+(bin2:bin3)/fftlen)/(f0+fMin))/lr;
      fp=floor(pf);
      r=[ones(1,bin2) fp fp+1 numOfFilters*ones(1,fn2-bin3)];
      c=[1:bin3+1 bin2+1:fn2+1];
      v=2*[0.5 ones(1,bin2-1) 1-pf+fp pf-fp ones(1,fn2-bin3-1) 0.5];
      lowestBin=1;
      highestBin=fn2+1;
    else
      b1=floor(binl(1))+1;
      b4=min(fn2,ceil(binl(4)))-1;
      pf=log((f0+(b1:b4)/fftlen)/(f0+fMin))/lr;
      fp=floor(pf);
      pm=pf-fp;
      k2=bin2-b1+1;
      k3=bin3-b1+1;
      k4=b4-b1+1;
      r=[fp(k2:k4) 1+fp(1:k3)];
      c=[k2:k4 1:k3];
      v=2*[1-pm(k2:k4) pm(1:k3)];
      lowestBin=b1+1;
      highestBin=b4+1;
    end
    if any(w=='n')
      v=1-cos(v*pi/2);  %hanning shaped filter
    elseif any(w=='m')
      v=1-0.92/1.08*cos(v*pi/2); %hamming shaped filter
    end
    if nargout > 1
      xmatrix=sparse(r,c,v);
    else
      xmatrix=sparse(r,c+lowestBin-1,v,numOfFilters,1+fn2);
    end

end




