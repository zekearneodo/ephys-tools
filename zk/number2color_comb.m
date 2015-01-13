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
   colorWeights=(colors(n,:));
   for i=1:3
       color=[color dec2hex(colorWeights(i),2)];
   end
end
   
    
    
    function colors=colors_list()
    %make 3 by 27 matrix with permutations of 0,1,2, n distinct elements
    nel=3;
    dim=3;
    nrep=1;
    m=[];
    
    for elem=1:nel
        list=[1:elem];
        ma=zeros(1,dim);
        ma(1:length(list))=list;
        m=[m;unique(perms(ma),'rows')];
    end
    %now repeating 2 numbers
    biglist=combnk(1:nel,2);
    for elem=1:nel
        for ibl=1:numel(biglist(:,1));
            list=[elem biglist(ibl,:)];
            ma=zeros(1,dim);
            ma(1:length(list))=list;
            m=[m;unique(perms(ma),'rows')];
        end
    end
    m=unique(m,'rows');
    colors=ceil(m*255/2);
    end

end

