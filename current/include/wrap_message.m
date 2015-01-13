function wrapped = wrap_message(msg,ch)
%wrap a message (string) into character ch
msgSize = length(msg);
wrapped = [repmat(ch,1,msgSize+3*2) sprintf('\n%s%s %s %s%s\n',ch,ch,msg,ch,ch) repmat(ch,1,msgSize+3*2)];
end