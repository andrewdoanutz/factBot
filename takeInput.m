
function filename=takeInput(option)
    if(option==1)
       fprintf("Recording\n");
       recObj = audiorecorder(44100,16,1);
       recordblocking(recObj, 3);
       myRecording = getaudiodata(recObj);
       audiowrite('input.wav', myRecording,44100);
       filename="input.wav";
       fprintf("Done recording\n");
    elseif(option==2)
       filename=input("Enter filename\n",'s');
    end
end
