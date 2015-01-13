function [ color ] = number2color( n )
%gets a number and gives back an html
%hope the number is smaller than 27
%if it is larger than 27, return white for all of them
%make it base 3
colors=colors_list();
color='#';
if n>length(colors) || n<1
    color='#ffffff'
else
    color=[color colors(n,:)];
   
end
   
    function colors=colors_list()
        b3colors=[];
        for i=1:26
            thisN=dec2base(i+9,3,3)-'0';
            b3colors=[b3colors;thisN(end-2:end)];
        end
        b3hColors=(ceil(b3colors/2*255));
        % a=dec2hex(b3hColors',2)
        a=reshape(transpose(dec2hex(b3hColors',2)),[],1);
        colors=transpose(reshape(a,6,[]));
    end

end

