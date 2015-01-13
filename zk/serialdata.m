classdef serialdata < handle
    %a channel with data received in serial format
    %it is a subclass of virtual serial
   
    properties(SetAccess = protected)
        %serial = (a virtual serial object)
        stream       = []; % the data vector
        streamT      = [];
        bytes        = [];
        words        = [];
        wordLength   = 2 ; % the number of bytes transmitted every time
        bitDur       = 0 ; % the length of a bit in ms 
        samplingRate = 19385; %the sampling frequency
        serial       = virtualserial;
        
    end
    
    methods
        %constructor
        function obj = serialdata(varargin)
            % variable arguments in:
            % the parameters for the serialdata (all manda
            %tory, except
            % optional wordLength that has a default)
            %and a vector (5 element vector) for the serial port properties.
            %the serial port properties are assigned altogether and in order
            %(baudRate,nBits,parity,stopBits)
            
            inPar=inputParser;
            inPar.addParamValue('serial',{300,8,0,1,1000},@(x) iscell(x) && numel(x)==5);
            inPar.addParamValue('wordLength',2,@(x) isscalar(x));
            inPar.addParamValue('samplingFrequency',19385,@(x) isscalar(x))
            
            inPar.parse(varargin{:});
            serialConfig = inPar.Results.serial;
            
            %initialize the object
            fprintf('Initializing serialdata object...\n')
            obj.wordLength   = inPar.Results.wordLength;
            obj.samplingRate = inPar.Results.samplingFrequency;
            obj.serial       = virtualserial(serialConfig{:});
            obj.bitDur       = 1./obj.serial.baudRate;
            
            %check nyquist condition for the baudRate of the serial
            if obj.samplingRate < 2 * obj.serial.baudRate
                error('Sorry, I cant handle such a high baudRate with such a low sampling frequency')
            end
            disp(obj);
        end %constructor
        
        %get the stream
        function obj = read_stream_bin(obj,fileName,chan,NChan,type)
            %read a stream from a binary file into the object
            fprintf('Reading serial stream from binary file %s (chan %d/%d, type %s) ...',fileName,chan,NChan,type);
            fileID=fopen(fileName,'r');
            fseek(fileID, 2*(chan-1), 'bof');
            obj.stream  = fread(fileID, inf, type, (NChan-1)*2);
            obj.streamT = (0:numel(obj.stream)-1)*1./obj.samplingRate;
            fclose(fileID);
            fprintf('done!\n');
        end %read_stream_bin
        
        %get a stream from a vector
        function obj = put_stream(obj,streamIn)
            %read a stream from a vector into the object
            fprintf('Getting serial stream from an input vector...');
            obj.stream  = streamIn;
            obj.streamT = (0:numel(obj.stream)-1)*1./obj.samplingRate;
            fprintf('done!\n');
        end %put_stream
        
        %make it a lo/hi vector
        function obj = clip_stream(obj)
            fprintf('Clipping stream into logical states...');
            midPoint = min(obj.stream) + range(obj.stream)/2;
            hi = find([obj.stream]>midPoint);
            lo = find([obj.stream]<=midPoint);
            obj.stream(hi)=1;
            obj.stream(lo)=0;
            fprintf('done!\n');
        end %clip_stream(obj)
        
        %get the beginnings of words
        function obj = get_bytes(obj)
            fprintf('Parsing bytes out of the stream...');
            s=obj.serial;
            
            bytePars = struct('dur',[], 'len',[], 'bitDur',[],'bitLen', [],'totalBits',[], 'bitMidPoints',[],...
                'start',[],'end',[],'error',[],'bitValues',[]); 
            
            auxByte = struct('bitMidPoints',[],'start',[],'end',[],'error',[],'bitValues',[]); 
            
            bytes(1) = auxByte;
            
            down = [];
            up   = [];
            loSegments   = [];
            
            %check that there is a stream
            if isempty(obj.stream)
                error('trying to parse bytes in an empty stream')
            end
            
            %check if the stream has been clipped (1 0)
            if range(obj.stream)>1
                fprintf('\n\t');
                obj.clip_stream;
            end
            
            %get the ups and downs
            down = find(diff(obj.stream)==-1);
            up   = find(diff(obj.stream)==1);
            
            %make sure that it starts high and finishes high
            if up(1)<down(1)
                up(1) = [];
            end
            if down(end)>up(end)
                down(end)=[];
            end
          
            %get the low segments
            loSegments = [down+1 up(1:length(down))]';
                
            %search lo segments and make bytes.
            bytePars.totalBits = (1+s.nBits+s.stopBits+s.parity);
            bytePars.dur     = (1+s.nBits+s.stopBits+s.parity) * obj.bitDur; %duration of a byte in ms
            bytePars.len     = round(bytePars.dur*obj.samplingRate);
            bytePars.bitDur  = obj.bitDur;
            bytePars.bitLen  = round(obj.bitDur*obj.samplingRate);
            
            
            auxByte.end  = 0;
            moreBytes = 1;
            while moreBytes;
                %start of next byte is first lo after end of last byte
                nextStartLeftBound = auxByte.end - ceil(bytePars.bitLen*(obj.serial.stopBits-1/2));
                nextStartSample    = down(find(down>nextStartLeftBound,1));
                auxByte = obj.get_a_byte(nextStartSample,bytePars);
                %if it was a good byte, add a new byte to the array,
                %otherwise replace it with the next one
                if ~auxByte.error
                    bytes =[bytes auxByte];
                    if auxByte.end - ceil(bytePars.bitLen*(obj.serial.stopBits-1/2)) > down(end)
                        moreBytes = 0;
                    end
                elseif auxByte.error==4
                    moreBytes=0;
                end
            end
            bytes(1)=[];
            obj.bytes = bytes;
            fprintf(' done! (found %d bytes)\n',numel(bytes));

        end %get_bytes
  
        function gotByte = get_a_byte(obj,startSample,bytePars)
            % get a byte starting at a position
            % un-Synched: don't re-sinchronize every up/down
            % just use the midpoint of every bit position counting from start
                        
            % start is the start of first lo sample
            gotByte.start = startSample;
            % end is the end of the last stop bit
            gotByte.end   = startSample+bytePars.len;
            gotByte.error = 10;
            
            %mid position of every bit (in samples)
            if gotByte.start == 44945326
            disp(gotByte.start)
            end
            gotByte.bitMidPoints = gotByte.start + round((1:bytePars.totalBits)*bytePars.bitLen-bytePars.bitLen/2.);
            
            % check that it ends with n stopBits in high (first check that it
            % is well defined.
            if obj.stream(gotByte.bitMidPoints(1))>0
                return_error(1);
                %it doesnt start with a low start bit
            elseif sum(obj.stream(gotByte.bitMidPoints(end-obj.serial.stopBits:end)))<obj.serial.stopBits
                return_error(2);
                %it doesn't end with n high stop bits
            elseif obj.serial.parity && ~parity_check(gotByte)
                return_error(3)
                %it doesn't check parity
            elseif gotByte.end > length(obj.stream)
                return_error(4)
                % it got to the end of the stream before the end of the
                % byte
            end
            
            %if everything is ok read the byte
     
            gotByte = read_byte_asynch(obj,gotByte);
            
            
            %_____________________________________%
            function gotByte = read_byte_asynch(obj,gotByte)
                %read the byte just by the midpoint positions (no
                %resynching every up/down
                gotByte.bitValues = zeros(1,obj.serial.nBits);
                gotByte.bitValues = obj.stream(gotByte.bitMidPoints(1+(1:obj.serial.nBits)));
                gotByte.error     = 0;
            end
            %____%
            
            function parityCheck=parity_check()
                % I still don't know how to check parity
                parityCheck=0;
            end
            %____%
            
            function return_error(errorCode)
                gotByte.error=errorCode;
                warning('Bad byte sarting in sample %d (error code %d)',gotByte.start,gotByte.error);
            end
            %____%

        end % function get_a_byte
        
        function obj = get_words(obj)
            % goes through the bytes and makes the n-bytw words
            % returns an array of words:
            % word.bytes = [wordLength x nBits] : array of the bytes and the
            % bits values.
            % word.starts = [wordLenght x 1] : vector of starts of the bytes in samples
            % word.formatted = int16 : for now, convert the word into a 2-byte
            % int
            
            word = struct('bytes',[],'starts',[],'formatted',[]);
            fprintf('Getting all the %d-bytes words in the serial stream...',obj.wordLength);
            
            if isempty(obj.bytes)
                fprintf('\n\t');
                obj.get_bytes();
            end
            
            starts = [obj.bytes.start];
            ends   = [obj.bytes.end];
            
            %find groups of bytes within wordLength between consecutive ends
            %and starts.
            samplesTimeOut = obj.serial.timeOut*1e-3 * obj.samplingRate;
            lags = starts(2:end)-ends(1:end-1);
            consecutive = zeros(size(lags));
            
            consecutive(find(lags<samplesTimeOut))=1;
            %get blocks of at least wordLength consecutive:
            % gets starts of chunks between timeOuts
            chunkStarts = find(diff(consecutive)==1)-1;
            % if the last block is shorter, dont count it
            if chunkStarts(end) > numel(consecutive) - obj.wordLength
                chunkStarts(end) = [];
            end
            chunkEnds = find(diff(consecutive)==-1)+1;
            if chunkEnds(1) < chunkStarts(1)
                chunkEnds(1) = [];
            end
            %get blocks of at least wordLength consecutive:
            wordStarts = [];
            for iChunk = 1: numel(chunkStarts)
                wordStarts = [wordStarts parse_chunk(chunkStarts(iChunk):chunkEnds(iChunk),obj.wordLength)];
            end
            
            % got all the wordStarts
            % each wordStart is the timeStamp fo the beginning of a byte
            % now get those words.
            for iW = 1: numel(wordStarts)
                words(iW) = get_word(wordStarts(iW),obj.wordLength);
            end
            
            obj.words = words;
            fprintf(' done! (found %d words)\n',numel(words));
            %_____________________________________%
            function wStarts = parse_chunk(chunk, wLen)
                %gets a chunk of consecutive byte time stamps and parses it into words
                %of wLen bytes
                %get the number of words in the chunk
                wStarts = [];
                nWords  = floor(numel(chunk)/wLen);
                if nWords > 0
                    wStarts = chunk((0:nWords-1)*nWords + 1);
                end
            end
            %____%          
            
            function word = get_word(wStart, wLen)
                %get a word starting at sample wStart
                word = struct('bytes',[],'starts',[],'formatted',[]);
                
                startingByteIndex = wStart;
                word.bytes     = [obj.bytes((0:wLen-1) + startingByteIndex)];
                word.starts    = [word.bytes.start];
                word.formatted = bi2de(reshape([word.bytes.bitValues],1,wLen*obj.serial.nBits));
            end

            %____%
            
            
        end %get_words
        
        function [stamps, numbers] = get_numbers(obj)
            fprintf('Getting all the formatted words in the serial stream...');
            if isempty(obj.words)
                fprintf('\n\t');
                obj.get_words();
            end
            wStarts = reshape([obj.words.starts],obj.wordLength,numel(obj.words));
            stamps  = wStarts(1,:);
            numbers = [obj.words.formatted];
        end %get_numbers
        
    end %methods
        
end %classdef