function [qCells] = find_lightrals()

fn=file_names;

nCells=0;
find_qCells=[];
cellBaseName='*.mat';
cellsList=dir(fullfile(fn.fold_unit_db,cellBaseName));
find_qCells=cellfun(@(x) getCell(fullfile(fn.fold_unit_db,x)),{cellsList.name},'UniformOutput',false);
find_qCells(cellfun('isempty',find_qCells))=[];
qCells(nCells+1:nCells+numel(find_qCells))=[find_qCells{:}];

%_______________________________________________________________________%
    function theCell=getCell(unit_filename)
        theCell=load(unit_filename);
        if theCell.quality==1 && theCell.light==1 && theCell.odor==1
            theCell.resp = get_resp_struct(unit_filename,'odor');
        else
            theCell='';
        end
    end
end