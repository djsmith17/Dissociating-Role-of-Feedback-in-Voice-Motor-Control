function boxPos = dfSetDialBoxPos(wMon)

monitorDims = get(0, 'Monitor');
numMon      = size(monitorDims, 1);

boxW = 237;
boxH = 145;

boxCPos = [0.5 0.6];

if wMon == 2 && numMon > 1
    sizeMon = monitorDims(2, 3:4);
    
    boxWpos = boxCPos(1) - (boxW/2)/sizeMon(1) + 1;
    boxHpos = boxCPos(2) - (boxH/2)/sizeMon(2);
    boxPos = [boxWpos boxHpos];
else
    sizeMon = monitorDims(1, 3:4);
    
    boxWpos = boxCPos(1) - (boxW/2)/sizeMon(1);
    boxHpos = boxCPos(2) - (boxH/2)/sizeMon(2);
    boxPos = [boxWpos boxHpos];
end
end