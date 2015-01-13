
function plot_solubility_affinity()
load('M72_odorProperties.mat')

lig_colors = {'r', 'g', 'b', 'k', 'm', 'r', 'g', 'b', 'k', 'm', 'r', 'g', 'b', 'k', 'm', 'r', 'g', 'b', 'k', 'm'};
lig_shapes = {'s', 'v', '^', 'o'};

figure(1); clf; hold on
for ilig = 1:numel(M72_odorProperties)
    
%     EC50 = M72_odorProperties(ilig).EC50;
%     if isfinite(EC50)
%         plot(EC50,M72_odorProperties(ilig).Solubility_water,...
%             lig_shapes{5-ceil((4*ilig)/20)},...
%             'Color', lig_colors{ilig},...
%             'MarkerSize',15)
%     end
    
    pA = M72_odorProperties(ilig).pA_10uM;
    if isfinite(pA)
        plot(pA,M72_odorProperties(ilig).Solubility_water,...
            lig_shapes{5-ceil((4*ilig)/20)},...
            'Color', lig_colors{ilig},...
            'MarkerSize',15)
    end
    xlim([0 180])
    ylim([0 25])
    xlabel('affinity')
    ylabel('solubility')
    
    sprintf('odor is %s, marker %s %s', M72_odorProperties(ilig).odorName, lig_colors{ilig}, lig_shapes{5-ceil((4*ilig)/20)})
%     pause()
end

end
