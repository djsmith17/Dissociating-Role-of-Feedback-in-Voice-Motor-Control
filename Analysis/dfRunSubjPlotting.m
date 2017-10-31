function dfRunSubjPlotting()
%The function that draws the plots from analyzed results

%Plot Toggles. This could eventually become an input variable
PltTgl.ForceSensor     = 0; %Voltage trace of force sensor signal
PltTgl.IntraTrial_T    = 0; %SPL trace of individual trial
PltTgl.IntraTrial_f0   = 0; %f0 trace of individual trial
PltTgl.InterTrial_f0   = 1; %Average f0 trace over all trials of a run
PltTgl.InterRun_f0       = 1; %Average f0 trace over all runs analyzed
PltTgl.InterTrial_AudRes = 1; %Average f0 response trace to auditory pert trials of a run
PltTgl.InterRun_AudRes   = 1; %Average f0 response trace to auditory pert over all runs analyzed
PltTgl.InterTrial_Force  = 0;
PltTgl.InterRun_Force    = 0;
PltTgl.svInflaRespRoute  = 0;


drawDAQAll(niAn, 2, dirs.SavResultsDir, 1)
drawInterTrialf0(auRes.timeSec, auRes.meanTrialf0_St, auRes.meanTrialf0_Sp, auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
drawAllTrialf0(auRes.time, auRes.allTrialf0, auRes.runTrialOrder, auAn.trigsT, auRes.f0Limits, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
drawAudResp_AllTrial(auRes, auAn.curSess, DRF.expParam.curSess, dirs.SavResultsDir)

drawAudResp_InterTrial(auRes.timeSec, auRes.meanTrialf0_St(:,:,2), auRes.meanTrialf0_Sp(:,:,2), auRes.f0LimitsSec, auRes.trialCount, auRes.meanTrialf0b, auAn.curSess, '', dirs.SavResultsDir)
         
if PltTgl.svInflaRespRoute == 1
    InflaRespRoute = CalcInflationResponse(auAn, auRes.meanTrialf0b, auRes.meanTrialf0_St, auRes.InflaRespLimits, dirs.SavResultsDir);
    tStep = auAn.tStep;
    save(dirs.InflaRespFile, 'InflaRespRoute', 'tStep')
end