function out=plot_if(varargin)
%conditional plot:
%plot if optional variable (first one) is 'plot' (default), don't plot if it is
%'noplot'
%Example:
%plot_if(x,y,...)
%will plot x,y, and all the options and return the handle to the plot
%plot_if('noplot',x,y,...)
%will not plot anything and return -1.


doplot=true;

xplot=varargin{1};
if any(strcmpi(xplot,{'plot','noplot'}))
    if strcmpi(xplot,'noplot')
        doplot=false;
    end
end

if doplot
    out=plot(varargin{2:end});
else
    out=-1;
end

end