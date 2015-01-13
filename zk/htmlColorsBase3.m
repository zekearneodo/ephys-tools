b3colors=[];
for i=1:26
    thisN=dec2base(i*2,3,3)-'0';
    b3colors=[b3colors;thisN(end-2:end)];
end
b3hColors=(ceil(b3colors/2*255))
% a=dec2hex(b3hColors',2)
a=reshape(transpose(dec2hex(b3hColors',2)),[],1)
colors=transpose(reshape(a,6,[]))