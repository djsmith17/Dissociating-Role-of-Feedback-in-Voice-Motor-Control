function [expParam, p] = dfInitAudapater(dirs, expParam)

% Paradigm Configurations
expParam.sRate              = 48000;
expParam.frameLen           = 96;
expParam.downFact           = 3;
expParam.sRateAnal          = expParam.sRate/expParam.downFact;
expParam.frameLenDown       = expParam.frameLen/expParam.downFact;
expParam.audioInterfaceName = 'MOTU MicroBook'; %'ASIO4ALL' 'Komplete'

% Set up Audapter
Audapter('deviceName', expParam.audioInterfaceName);
Audapter('setParam', 'downFact', expParam.downFact, 0);
Audapter('setParam', 'sRate', expParam.sRateAnal, 0);
Audapter('setParam', 'frameLen', expParam.frameLenDown, 0);
p = getAudapterDefaultParams(expParam.gender);

if strcmp(expParam.expType(1:3), 'Aud')
    ostName = 'AFPerturbOST.ost';
    pcfName = 'AFPerturbPCF.pcf';
else
    ostName = 'SFPerturbOST.ost';
    pcfName = 'SFPerturbPCF.pcf';
end

% Set up OST and PCF Files
expParam.ostFN = fullfile(dirs.Prelim, ostName); check_file(expParam.ostFN);
expParam.pcfFN = fullfile(dirs.Prelim, pcfName); check_file(expParam.pcfFN);
end