%set of functions to read/write individual channels from binary files, in
%int16 format (arranged as in spikeGL: ch1(i) ch2(i) ...; ch1(n)
%ch2(n)..chm(n)
%zk, based on DR data_management_tools;
%Feb 2014

function ft=file_tools()

ft.read_analog_channel  = @read_analog_channel;
ft.write_analog_channel = @write_analog_channel;
end

function Y = read_analog_channel(fid, ich, nch)
fprintf('Reading channel %d/%d ...',ich,nch);
fseek(fid, 2*(ich-1), 'bof');
Y = fread(fid, inf, 'int16', (nch-1)*2);
fprintf(' read!\n');
end

function write_analog_channel(X, fid, offset, kch, nch)
fprintf('Writing channel %d/%d ...',kch,nch);
fseek(fid, (offset*nch+kch-1)*2, 'bof');
fwrite(fid, X(1), 'int16');
fwrite(fid, X(2:end), 'int16', (nch-1)*2);
fprintf(' written!\n');
end