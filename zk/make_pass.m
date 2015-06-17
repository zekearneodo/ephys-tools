%quick and dirty
%from cellsArray get a unit meta and run the visualize response
cp = cell_passport_tools();
%open cellsArray

%select cellsarray
cellsArray(~([cellsArray.light]==1))=[]
for i=12:numel(cellsArray)
    a_unit=cellsArray(i);
    try
    pa = cp.make_passport(a_unit);
    catch me
        continue
    end
end
