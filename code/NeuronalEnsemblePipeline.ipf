#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later



// -----------------------------------------------------------------------------
// Detects positive and negative calcium events using a local pre-stimulus
// baseline and standard deviation calculated independently for each ROI
// and stimulus interval.
//
// Output:
//   <DB>_events : signed event matrix (+1 positive, -1 negative, 0 no event)
// -----------------------------------------------------------------------------

Function GetCalciumEvents_PreStimWindows1(DB,thresholdfactor, minDistance, startTimes, endTimes)


    Wave DB
    Variable thresholdfactor,minDistance 
    Wave startTimes
    Wave endTimes

    Variable i, j, w
    Variable nwaves,numpoints2
    Variable baseline, sigma, threshold 
    Variable p1, p2
    Variable preStart, preEnd
    Variable detectStart, detectEnd
    Variable lastEventPos
    Variable nIntervals
    

    String arrayname = NameOfWave(DB)
    String eventname = NameOfWave(DB) + "_events"

    nwaves = DimSize($arrayname, 1)
    numpoints2 = DB[%npoints][0]

    Variable deltat = DB[%XDelta][0]
    nIntervals = numpnts(startTimes)

    // =========================
    // OUTPUT WAVES
    // =========================
    Make/O/N=(numpoints2, nwaves) events_exc, events_inh, events_all
    events_exc = 0
    events_inh = 0
    events_all = 0

    Make/O/N=(nwaves) events_per_cell_exc, events_per_cell_inh
    events_per_cell_exc = 0
    events_per_cell_inh = 0

    SetScale/P x 0, deltat, "", events_exc, events_inh, events_all

    Make/O/N=(numpoints2) tempwave

    // =========================
    // PRE-STIMULUS THRESHOLD WAVES
    // rows = stimuli, columns = ROIs
    // =========================
    Make/O/N=(nIntervals, nwaves) prestim_baseline_byStim
    Make/O/N=(nIntervals, nwaves) prestim_sigma_byStim
    Make/O/N=(nIntervals, nwaves) prestim_threshold_byStim
    Make/O/N=(nIntervals, nwaves) prestim_posLimit_byStim
    Make/O/N=(nIntervals, nwaves) prestim_negLimit_byStim

    
    prestim_baseline_byStim = NaN
    prestim_sigma_byStim = NaN
    prestim_threshold_byStim = NaN
    prestim_posLimit_byStim = NaN
    prestim_negLimit_byStim = NaN

  
    // =========================
    // MAIN LOOP ROIs
    // =========================
    for(i=0; i<nwaves; i+=1)

        Duplicate/O/R=[0,(numpoints2-1)][i] $arrayname, tempwave
        Redimension/N=(numpoints2), tempwave


        events_per_cell_exc[i] = 0
        events_per_cell_inh[i] = 0

        lastEventPos = -minDistance

        // =========================
        // LOOP STIM BY STIM
        // =========================
        for(w=0; w<nIntervals; w+=1)

            p1 = round(startTimes[w])
            p2 = round(endTimes[w])

            if (p1 < 0)
                p1 = 0
            endif

            if (p2 >= numpoints2)
                p2 = numpoints2 - 1
            endif

            if (p2 <= p1)
                continue
            endif

            // =========================
            // DEFINE CURRENT ANALYSIS SEGMENT
            // =========================
            if (w == 0)
                detectStart = 0
                preStart = 0
            else
                detectStart = round(endTimes[w-1]) + 1
                preStart = detectStart
            endif

            detectEnd = p2
            preEnd = p1 - 1

            if (detectStart < 0)
                detectStart = 0
            endif

            if (detectEnd >= numpoints2)
                detectEnd = numpoints2 - 1
            endif

            // Skip the stimulus if the pre-stimulus segment is not long enough
            if (preEnd <= preStart)
                continue
            endif

            // =========================
            // LOCAL PRE-STIMULUS BASELINE
            // =========================
            WaveStats/Q/R=[preStart, preEnd] tempwave  
           

            baseline = V_avg
            sigma = V_sdev
            threshold = thresholdfactor * sigma
          

            prestim_baseline_byStim[w][i] = baseline
            prestim_sigma_byStim[w][i] = sigma
            prestim_threshold_byStim[w][i] = threshold
            prestim_posLimit_byStim[w][i] = baseline + threshold
            prestim_negLimit_byStim[w][i] = baseline - threshold


            // =========================
            // EVENT DETECTION ACROSS THE CURRENT SEGMENT
            // detection is not restricted to the stimulus window
            // =========================
            for(j=detectStart; j<=detectEnd; j+=1)

                if ((j - lastEventPos) < minDistance)
                    continue
                endif

                // POSITIVE EVENT
                if (tempwave[j] > baseline + threshold)

                    Variable amp = tempwave[j] - baseline

                    if (amp > 0)

                        events_exc[j][i] = 1
                        events_all[j][i] = 1

                        events_per_cell_exc[i] += 1
                        lastEventPos = j

                    endif

                endif

                // NEGATIVE EVENT
                if (tempwave[j] < baseline - threshold)

                    Variable amp2 = baseline - tempwave[j]

                    if (amp2 > 0)

                        events_inh[j][i] = -1
                        events_all[j][i] = -1

                        events_per_cell_inh[i] += 1
                        lastEventPos = j

                    endif

                endif

            endfor

        endfor

        // =========================
        // MARK ACTIVE ROI
        // =========================
        if ((events_per_cell_exc[i] + events_per_cell_inh[i]) > 0)
            DB[%ONOFF][i] = 1
        else
            DB[%ONOFF][i] = 0
        endif

    endfor

    // =========================
    // MERGE
    // =========================
    events_all = events_exc + events_inh

    Duplicate/O events_all, $eventname

    Print "============================"
    Print "PRE-STIM THRESHOLD EVENT DETECTION"
    Print "Number of stim windows = " + num2str(nIntervals)
    Print "Output wave = " + eventnam
    Print "============================"

    KillWaves tempwave

End



// -----------------------------------------------------------------------------
// Classifies each stimulus-response window according to the balance of
// positive and negative events.
//
// Window codes:
//   1  = ON
//  -1  = OFF
//   0  = no response
//   2  = mixed response
//
// Output:
//   stim_class_by_window : rows = stimuli, columns = ROIs
// -----------------------------------------------------------------------------

Function ClassifyByStim2(DB, startTimes, endTimes)

    Wave DB
    Wave startTimes
    Wave endTimes

    Variable i, j, w
    Variable nwaves, numpoints2
    Variable p1, p2
    Variable nIntervals

    String spikearrayname = NameOfWave(DB) + "_events"
    Wave spikes = $spikearrayname

    nIntervals = numpnts(startTimes)

    nwaves     = DimSize(spikes, 1)
    numpoints2 = DimSize(spikes, 0)

    // =========================
    // ENSURE REQUIRED DATABASE FIELD
    // =========================
    Variable stimIndex = FindDimLabel(DB,0,"stims")
    if (stimIndex < 0)
        Variable oldRows = DimSize(DB,0)
        Redimension/N=(oldRows+1, DimSize(DB,1)) DB
        SetDimLabel 0, oldRows, stims, DB
    endif

    // =========================
    // WINDOW-LEVEL CLASSIFICATION OUTPUT
	// rows = stimuli, columns = ROIs
    // =========================
    Make/O/N=(nIntervals, nwaves) stim_class_by_window

   
    stim_class_by_window = 0

    // stim_class_by_window:
    //  1  = ON window
    // -1  = OFF window
    //  0  = no response
    //  2  = mixed/tie

    // =========================
    // LOOP ROIs
    // =========================
    for(i=0; i<nwaves; i+=1)

        Variable stimResponses = 0
        Variable Excitatory = 0
        Variable Inhibitory = 0
        Variable Mixed = 0

        for(w=0; w<nIntervals; w+=1)

            p1 = round(startTimes[w])
            p2 = round(endTimes[w])

            if (p1 < 0)
                p1 = 0
            endif

            if (p2 >= numpoints2)
                p2 = numpoints2 - 1
            endif

            if (p2 <= p1)
                continue
            endif

            Variable nPos = 0
            Variable nNeg = 0

            // Count all positive and negative events within the response window
            for(j=p1; j<=p2; j+=1)

                if(spikes[j][i] == 1)
                    nPos += 1
                elseif(spikes[j][i] == -1)
                    nNeg += 1
                endif

            endfor

          

            // Assign the window-level response class
            if(nPos > nNeg)

                Excitatory += 1
                stimResponses += 1
                stim_class_by_window[w][i] = 1

            elseif(nNeg > nPos)

                Inhibitory += 1
                stimResponses += 1
                stim_class_by_window[w][i] = -1

            elseif(nPos == 0 && nNeg == 0)

                stim_class_by_window[w][i] = 0

            else

                Mixed += 1
                stimResponses += 1
                stim_class_by_window[w][i] = 2

            endif

        endfor
		// Legacy database fields:
		// Baseline stores the number of ON windows
		// Tsus stores the number of OFF windows
		
        DB[%stim][i] = stimResponses
        DB[%Baseline][i] = Excitatory
        DB[%Tsus][i] = Inhibitory

    endfor

    Print "============================"
    Print "CLASSIFY BY STIM 2 FINISHED"
    Print "Classification is based on balance within each window"
    Print "stim_class_by_window: 1=ON, -1=OFF, 0=no response, 2=mixed/tie"
    Print "============================"

End

// -----------------------------------------------------------------------------
// Visualizes the calcium trace, detected positive and negative events,
// and the local pre-stimulus thresholds used for event detection.
//
// Requires prestim_posLimit_byStim and prestim_negLimit_byStim,
// generated by the event-detection functions.
//
// This function is intended for visual inspection of individual ROIs.
// It does not modify the final ROI classification.
// -----------------------------------------------------------------------------


Function OverlayCaPN_WithThresholds(DB, wavenum, startTimes, endTimes, plot)

    Wave DB
    Variable wavenum, plot
    Wave startTimes
    Wave endTimes

    Variable numpoints
    Variable deltat = DB[%XDelta][0]

    String arrayname = NameOfWave(DB)
    String spikearrayname = NameOfWave(DB) + "_events"

    Wave spikes = $spikearrayname
    Wave posLimit = prestim_posLimit_byStim
    Wave negLimit = prestim_negLimit_byStim

    
    numpoints = DB[%npoints][0]

    Make/O/N=(numpoints) Cawave, Spikewave
    Make/O/N=(numpoints) Spikewave_pos, Spikewave_neg
    Make/O/N=(numpoints) PosThresholdWave, NegThresholdWave

    Cawave = NaN
    Spikewave = 0
    Spikewave_pos = 0
    Spikewave_neg = 0
    PosThresholdWave = NaN
    NegThresholdWave = NaN

    Duplicate/O/R=[0, numpoints-1][wavenum] $arrayname, Cawave
    Duplicate/O/R=[0, numpoints-1][wavenum] $spikearrayname, Spikewave

    SetScale/P x 0, deltat, "", Cawave, Spikewave
    SetScale/P x 0, deltat, "", Spikewave_pos, Spikewave_neg
    SetScale/P x 0, deltat, "", PosThresholdWave, NegThresholdWave

    // Separate positive and negative events
    Spikewave_pos = (Spikewave == 1) ? 1 : 0
    Spikewave_neg = (Spikewave == -1) ? -1 : 0

    Variable w, j
    Variable nIntervals = numpnts(startTimes)
    Variable p1, p2
    Variable detectStart, detectEnd

    // =========================
    // BUILD THRESHOLD WAVES FOR EACH ANALYSIS SEGMENT
    // =========================
    for(w=0; w<nIntervals; w+=1)

        p1 = round(startTimes[w])
        p2 = round(endTimes[w])

        if (p1 < 0)
            p1 = 0
        endif

        if (p2 >= numpoints)
            p2 = numpoints - 1
        endif

        if (p2 <= p1)
            continue
        endif

        if (w == 0)
            detectStart = 0
        else
            detectStart = round(endTimes[w-1]) + 1
        endif

        detectEnd = p2

        if (detectStart < 0)
            detectStart = 0
        endif

        if (detectEnd >= numpoints)
            detectEnd = numpoints - 1
        endif

        for(j=detectStart; j<=detectEnd; j+=1)
            PosThresholdWave[j] = posLimit[w][wavenum]
            NegThresholdWave[j] = negLimit[w][wavenum]
        endfor

    endfor

    if (plot == 1)

        Display/K=1 Cawave

        // Positive events
        AppendToGraph/R Spikewave_pos
        ModifyGraph rgb(Spikewave_pos) = (0,65535,0)

        // Negative events
        AppendToGraph/R Spikewave_neg
        ModifyGraph rgb(Spikewave_neg) = (65535,0,0)

        SetAxis right -1,1

        // Threshold lines displayed on the calcium-signal axis
        AppendToGraph PosThresholdWave
        ModifyGraph rgb(PosThresholdWave) = (0,40000,0)
        ModifyGraph lsize(PosThresholdWave) = 2
        ModifyGraph lstyle(PosThresholdWave) = 3

        AppendToGraph NegThresholdWave
        ModifyGraph rgb(NegThresholdWave) = (50000,0,0)
        ModifyGraph lsize(NegThresholdWave) = 2
        ModifyGraph lstyle(NegThresholdWave) = 3

        ModifyGraph rgb(Cawave) = (0,0,0)
        ModifyGraph gFont = "Helvetica", gfSize = 14

        Label bottom "Time (s)"
        Label left "\\F'Symbol'Δ\\F'Helvetica'F/F"
        Label right "event (+1 / -1)"

    endif

End


// -----------------------------------------------------------------------------
// Detects positive and negative calcium events using independent positive and
// negative threshold factors.
//
// For each ROI and stimulus interval, a local pre-stimulus baseline and
// standard deviation are calculated. Separate threshold factors are then
// applied to positive and negative responses.
//
// Output:
//   <DB>_events : signed event matrix (+1 positive, -1 negative, 0 no event)
//
// This alternative implementation is intended for datasets requiring
// asymmetric positive and negative detection thresholds.
// -----------------------------------------------------------------------------


Function GetCalciumEvents_PreStimWindows_TwoThresholds(DB,thresholdfactor_pos, thresholdfactor_neg, minDistance, startTimes, endTimes)

    Wave DB
    Variable thresholdfactor_pos, thresholdfactor_neg, minDistance
    Wave startTimes
    Wave endTimes

    Variable i, j, w
    Variable nwaves, numpoints2
    Variable threshold_pos, threshold_neg
    Variable baseline, sigma
    Variable p1, p2
    Variable preStart, preEnd
    Variable detectStart, detectEnd
    Variable lastEventPos
    Variable nIntervals
    Variable amp, amp2

    String arrayname = NameOfWave(DB)
    String eventname = NameOfWave(DB) + "_events"

    nwaves = DimSize($arrayname, 1)
    numpoints2 = DB[%npoints][0]

    Variable deltat = DB[%XDelta][0]
    nIntervals = numpnts(startTimes)

    // =========================
    // OUTPUT WAVES
    // =========================
    Make/O/N=(numpoints2, nwaves) events_exc, events_inh, events_all
    events_exc = 0
    events_inh = 0
    events_all = 0

    Make/O/N=(nwaves) events_per_cell_exc, events_per_cell_inh
    events_per_cell_exc = 0
    events_per_cell_inh = 0

    SetScale/P x 0, deltat, "", events_exc, events_inh, events_all

    Make/O/N=(numpoints2) tempwave

    // =========================
   	// PRE-STIMULUS THRESHOLD WAVES
	// rows = stimuli, columns = ROIs
    // =========================
    Make/O/N=(nIntervals, nwaves) prestim_baseline_byStim
    Make/O/N=(nIntervals, nwaves) prestim_sigma_byStim

    Make/O/N=(nIntervals, nwaves) prestim_threshold_pos_byStim
    Make/O/N=(nIntervals, nwaves) prestim_threshold_neg_byStim

    Make/O/N=(nIntervals, nwaves) prestim_posLimit_byStim
    Make/O/N=(nIntervals, nwaves) prestim_negLimit_byStim

 

    prestim_baseline_byStim = NaN
    prestim_sigma_byStim = NaN
    prestim_threshold_pos_byStim = NaN
    prestim_threshold_neg_byStim = NaN
    prestim_posLimit_byStim = NaN
    prestim_negLimit_byStim = NaN

   

    // =========================
    // MAIN LOOP ROIs
    // =========================
    for(i=0; i<nwaves; i+=1)

        Duplicate/O/R=[0,(numpoints2-1)][i] $arrayname, tempwave
        Redimension/N=(numpoints2), tempwave

       

        events_per_cell_exc[i] = 0
        events_per_cell_inh[i] = 0

        lastEventPos = -minDistance

        // =========================
        // LOOP STIM BY STIM
        // =========================
        for(w=0; w<nIntervals; w+=1)

            p1 = round(startTimes[w])
            p2 = round(endTimes[w])

            if (p1 < 0)
                p1 = 0
            endif

            if (p2 >= numpoints2)
                p2 = numpoints2 - 1
            endif

            if (p2 <= p1)
                continue
            endif

            // =========================
            // DEFINE CURRENT ANALYSIS SEGMENT
            // =========================
            if (w == 0)
                detectStart = 0
                preStart = 0
            else
                detectStart = round(endTimes[w-1]) + 1
                preStart = detectStart
            endif

            detectEnd = p2
            preEnd = p1 - 1

            if (detectStart < 0)
                detectStart = 0
            endif

            if (detectEnd >= numpoints2)
                detectEnd = numpoints2 - 1
            endif

            // Skip the stimulus if the pre-stimulus segment is not long enough
            if (preEnd <= preStart)
                continue
            endif

            if ((preStart + 5) > preEnd)
                continue
            endif

            // =========================
            // LOCAL PRE-STIMULUS BASELINE
            // =========================
            WaveStats/Q/R=[preStart + 5, preEnd] tempwave

            baseline = V_avg
            sigma = V_sdev

            threshold_pos = thresholdfactor_pos * sigma
            threshold_neg = thresholdfactor_neg * sigma

            prestim_baseline_byStim[w][i] = baseline
            prestim_sigma_byStim[w][i] = sigma

            prestim_threshold_pos_byStim[w][i] = threshold_pos
            prestim_threshold_neg_byStim[w][i] = threshold_neg

            prestim_posLimit_byStim[w][i] = baseline + threshold_pos
            prestim_negLimit_byStim[w][i] = baseline - threshold_neg

           // =========================
			// EVENT DETECTION ACROSS THE CURRENT SEGMENT
			// detection is not restricted to the stimulus window
			// =========================
			for(j=detectStart; j<=detectEnd; j+=1)

   				 if ((j - lastEventPos) < minDistance)
       					 continue
    		endif

                // =========================
                // POSITIVE EVENT
                // =========================
                if (tempwave[j] > baseline + threshold_pos)

                    amp = tempwave[j] - baseline

                    if (amp > 0)

                        events_exc[j][i] = 1
                        events_all[j][i] = 1

                        events_per_cell_exc[i] += 1
                        lastEventPos = j

                    endif

                endif

                // =========================
                // NEGATIVE EVENT
                // =========================
                if (tempwave[j] < baseline - threshold_neg)

                    amp2 = baseline - tempwave[j]

                    if (amp2 > 0)

                        events_inh[j][i] = -1
                        events_all[j][i] = -1

                        events_per_cell_inh[i] += 1
                        lastEventPos = j

                    endif

                endif

            endfor

        endfor

        // =========================
        // MARK ACTIVE ROI
        // =========================
        if ((events_per_cell_exc[i] + events_per_cell_inh[i]) > 0)
            DB[%ONOFF][i] = 1
        else
            DB[%ONOFF][i] = 0
        endif

    endfor

    // =========================
    // MERGE
    // =========================
    events_all = events_exc + events_inh

    Duplicate/O events_all, $eventname


    Print "============================"
    Print "PRE-STIM TWO-THRESHOLD EVENT DETECTION"
    Print "Number of stim windows = " + num2str(nIntervals)
    Print "Positive threshold factor = " + num2str(thresholdfactor_pos)
    Print "Negative threshold factor = " + num2str(thresholdfactor_neg)
    Print "Output wave = " + eventname
    Print "============================"

    KillWaves tempwave

End


// -----------------------------------------------------------------------------
// Assigns each ROI to an ON, OFF, or NON ensemble according to the
// reproducibility of its stimulus-evoked responses across repeated stimuli.
//
// Input:
//   StimClassWindows : rows = stimuli, columns = ROIs
//
// Window codes:
//   1  = ON
//  -1  = OFF
//   2  = mixed response
//   0  = no response
//
// Final ensemble codes:
//   1  = ON
//  -1  = OFF
//   0  = NON
// -----------------------------------------------------------------------------
Function ClassifyROIEnsembles(StimClassWindows, nStim)

    Wave StimClassWindows
    Variable nStim

    Variable nROIs = DimSize(StimClassWindows, 1)

    Make/O/N=(nROIs) EnsembleCode
    Make/O/T/N=(nROIs) EnsembleLabel

    Make/O/N=(nROIs) Ensemble_ON_Raw
    Make/O/N=(nROIs) Ensemble_OFF_Raw
    Make/O/N=(nROIs) Ensemble_MIXED_Raw
    Make/O/N=(nROIs) Ensemble_RESP_Total
    Make/O/N=(nROIs) Ensemble_ON_Final
    Make/O/N=(nROIs) Ensemble_OFF_Final

    Variable roi, stimIdx
    Variable val
    Variable onCount, offCount, mixedCount, respCount
    Variable onFinal, offFinal

    for (roi = 0; roi < nROIs; roi += 1)

        onCount = 0
        offCount = 0
        mixedCount = 0
        respCount = 0

        // Count responses across stimuli
        for (stimIdx = 0; stimIdx < nStim; stimIdx += 1)

            val = StimClassWindows[stimIdx][roi]

            if (val == 1)

                onCount += 1
                respCount += 1

            elseif (val == -1)

                offCount += 1
                respCount += 1

            elseif (val == 2)

                mixedCount += 1
                respCount += 1

            endif

        endfor

        // Store original response counts
        Ensemble_ON_Raw[roi] = onCount
        Ensemble_OFF_Raw[roi] = offCount
        Ensemble_MIXED_Raw[roi] = mixedCount
        Ensemble_RESP_Total[roi] = respCount

        // Initial final counts correspond to the original ON/OFF counts
        onFinal = onCount
        offFinal = offCount

        // Mixed responses are assigned to the majority response class.
        // If ON and OFF counts are equal, mixed responses are not assigned
        // to either class in order to avoid forcing the classification.
        if (mixedCount > 0)

            if (onCount > offCount)

                onFinal += mixedCount

            elseif (offCount > onCount)

                offFinal += mixedCount

            endif

        endif

        Ensemble_ON_Final[roi] = onFinal
        Ensemble_OFF_Final[roi] = offFinal

        // Default classification
        EnsembleCode[roi] = 0
        EnsembleLabel[roi] = "NON"

        // ============================
        // CLASSIFICATION FOR 5 STIMULI
        // ============================
        if (nStim == 5)

            if (respCount >= 4)

                if (respCount == 4)

                    // Four responsive windows:
                    // at least 3/4 must belong to the same response class
                    if ((onFinal >= 3) && (onFinal > offFinal))

                        EnsembleCode[roi] = 1
                        EnsembleLabel[roi] = "ON"

                    elseif ((offFinal >= 3) && (offFinal > onFinal))

                        EnsembleCode[roi] = -1
                        EnsembleLabel[roi] = "OFF"

                    endif

                elseif (respCount == 5)

                    // Five responsive windows:
                    // at least 4/5 must belong to the same response class.
                    // A 3/2 distribution remains classified as NON.
                    if ((onFinal >= 4) && (onFinal > offFinal))

                        EnsembleCode[roi] = 1
                        EnsembleLabel[roi] = "ON"

                    elseif ((offFinal >= 4) && (offFinal > onFinal))

                        EnsembleCode[roi] = -1
                        EnsembleLabel[roi] = "OFF"

                    endif

                endif

            endif

        // ============================
        // CLASSIFICATION FOR 4 STIMULI
        // ============================
        elseif (nStim == 4)

            if (respCount >= 3)

                if (respCount == 3)

                    // Three responsive windows:
                    // all three must belong to the same response class
                    if ((onFinal == 3) && (onFinal > offFinal))

                        EnsembleCode[roi] = 1
                        EnsembleLabel[roi] = "ON"

                    elseif ((offFinal == 3) && (offFinal > onFinal))

                        EnsembleCode[roi] = -1
                        EnsembleLabel[roi] = "OFF"

                    endif

                elseif (respCount == 4)

                    // Four responsive windows:
                    // at least 3/4 must belong to the same response class
                    if ((onFinal >= 3) && (onFinal > offFinal))

                        EnsembleCode[roi] = 1
                        EnsembleLabel[roi] = "ON"

                    elseif ((offFinal >= 3) && (offFinal > onFinal))

                        EnsembleCode[roi] = -1
                        EnsembleLabel[roi] = "OFF"

                    endif

                endif

            endif

        endif

    endfor

    Print "================================"
    Print "ENSEMBLE CLASSIFICATION COMPLETED"
    Print "Generated outputs:"
    Print " - EnsembleLabel: ON / OFF / NON"
    Print " - EnsembleCode: 1 = ON, -1 = OFF, 0 = NON"
    Print " - Ensemble_ON_Raw, Ensemble_OFF_Raw, Ensemble_MIXED_Raw"
    Print " - Ensemble_RESP_Total"
    Print " - Ensemble_ON_Final, Ensemble_OFF_Final"
    Print "================================"

End


// -----------------------------------------------------------------------------
// Displays an individual calcium trace and overlays the stimulus-search windows.
//
// StartTimes and EndTimes are provided in frames and converted to time using
// the x-axis scaling of the input trace.
//
// This function is intended for visualization only and does not modify
// event detection or ROI classification.
// -----------------------------------------------------------------------------
Function PlotSearchWindows(traceWave, StartTimes, EndTimes)

    Wave traceWave
    Wave StartTimes
    Wave EndTimes

    Variable w
    Variable tStart, tEnd
    Variable dt = DimDelta(traceWave, 0)
    Variable x0 = DimOffset(traceWave, 0)

    // Display the selected trace
    Display/K=1 traceWave

    Label bottom "Time (s)"
    Label left "\F'Symbol'D\F'Helvetica'F/F"

    // Draw stimulus-search windows behind the calcium trace
    SetDrawLayer UserBack

    for(w = 0; w < numpnts(StartTimes); w += 1)

        // StartTimes and EndTimes are expressed in frames
        tStart = x0 + StartTimes[w] * dt
        tEnd   = x0 + EndTimes[w] * dt

        // Draw the shaded search window
        SetDrawEnv xcoord=bottom, ycoord=prel, fillpat=1, fillfgc=(50000,58000,65535), linefgc=(55000,45000,0)
        DrawRect tStart, 0, tEnd, 1

    endfor

End