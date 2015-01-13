latencies=[  76 37 66 70 70;
            nan 27 nan 67 70;
            56 35 nan nan nan;
            78 40 nan 69 nan]
theLegend= {'Et tiglate','acetophenone', 'methyl salic' , '4-methyl'};
gfV=figure
plot(1:length(latencies),latencies,'.','MarkerSize',24);
xlim([0 length(latencies)+1]);
ylim([0 max(max(latencies))+10]);
legend(theLegend)

        
            