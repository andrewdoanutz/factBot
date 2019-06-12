  x=2;
    if((x==1) || (x==2))
        filename="output.wav";
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
        [digitResult,digitAcc]=classifyWord(temp,fs);
        [learnResult,learnAcc]=classifyWord(learn,fs);
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
            output=command;
            fprintf("doing %s\n",command);
        elseif(command=="learn 2")
            %do command 2
            output=command;
            fprintf("doing %s\n",command);
        elseif(command=="learn 3")
            %do command 3
            output=command;
            fprintf("doing %s\n",command);
        else
            fprintf("%s is not a valid command\n",command);
            output="error\n";
        end
    else
        fprintf("Entered %d, only 1 or 2 is allowed\n",x);
        output="error\n";
    end
