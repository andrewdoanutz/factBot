function wordSignals=splitWords(audio,fs)


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

        %if(((1/Len)*(sum(aud(start:stop))))>-.0001)
            %sFlux(n)=(1/Len)*(sum(aud(start:stop)));
        %end
    end
    for(n=1:iNumBlocks)
            i_start = (n-1)*iHopLength+1;
            i_stop = min(length(audio),i_start+L-1);
            energy(n) = mean((audio(i_start:i_stop).^2.*hann(i_stop-i_start+1))); %get energy 
    end
    
    threshold = energy > 0.03;  %if below threshold, index entry is 0, if above entry is 1
    segStarted = false;
    segEnded = false;
    starter = []; %array that holds index for where the word starts
    ender = []; %array that holds index for where word stops
    for a = 1:length(threshold) 
        if (threshold(a) == 1) && segStarted == false
            starter = [starter,Tblock2(a)];
            segStarted = true;
            segEnded = false;
        end
        if(threshold(a) == 0) &&segEnded == false
            if(length(starter)~=0)
                ender = [ender,Tblock2(a)];
            end
            segStarted = false;
            segEnded = true;
        end
    end
    factor = length(audio)/Tblock2(length(Tblock2));
    ender = ceil(ender*factor)+1000;
    starter = floor(starter*factor)-1000;
    zero=zeros(1000,1);

    for i=1:length(starter)
        if(i==1)
%             res=audio(starter(i):ender(i));
            wordSignals={audio(starter(i):ender(i))};
        else
            
%             res=[zero;res];
%             res=[res;audio(starter(i):ender(i))];
            wordSignals=[wordSignals;{audio(starter(i):ender(i))}];
        end
        

    end

    %soundsc(res,fs);

%     sFlux=smoothdata(rms_specflux);
%     baseFlux = sFlux > 1; %If spectral flux is smaller than 0.
%     factor = length(audio)/length(baseFlux);  %find the factor for the flux to scale back into original size
%     baseFlux = imresize(baseFlux,factor);
%     baseFlux = transpose(baseFlux(1,:));
%     audio = audio.*baseFlux(1:length(audio),:); %remove areas of samples to improve coherency.
    %subplot(2,1,1);plot(audio);
     %subplot(2,1,2);plot(res);


end
