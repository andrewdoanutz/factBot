  x=2;
    if((x==1) || (x==2))
        filename="sample.wav";
        [aud,fs]=audioread(filename);
        filterBlock=ones(1,50)/50;
        aud=filter(filterBlock,1,aud);
        aud=aud*100;
        wordSignals=splitWords(aud,fs);
        temp=[];
        for i=1:length(wordSignals)
            if(i==1)
               learn=[wordSignals{i}];
            elseif(i==2)
                temp=[wordSignals{i}];
            else
                temp=[temp; wordSignals{i}];
            end
        end
        [digitResult,digitAcc]=classifyGaussianWord(temp,fs);
        [learnResult,learnAcc]=classifyGaussianWord(learn,fs);
        if(learnResult=='0')
            learnResult='learn ';
        end
        words=[learnResult];
        words=[words digitResult];
        %if(length(words)>1)
        %command=strjoin(words);
        %else
            command=words;
        %end
        if(command=="learn 1")
            %do command 1
            output="The Roman People claim their heritage from the Trojans, as described in the Iliad by Virgil";
            fprintf("doing %s\n",command);
        elseif(command=="learn 2")
            %do command 2
            output="The mitochondria is the powerhouse of the cell";
            fprintf("doing %s\n",command);
        elseif(command=="learn 3")
            %do command 3
            output="The blue whale is the largest animal to ever exist on the planet";
            fprintf("doing %s\n",command);
        else
            fprintf("%s is not a valid command\n",command);
            output="error";
        end
    else
        fprintf("Entered %d, only 1 or 2 is allowed\n",x);
        output="error";
    end
    