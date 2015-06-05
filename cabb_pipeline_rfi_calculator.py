from __future__ import print_function
from cabb_pipeline_user_routines import *
import numpy as np
import matplotlib.pyplot as plt
import re
import json

"""
This routine is run as part of the CABB data reduction pipeline.
It is designed to use the flagging fractions calculated from the automatic
flagging to figure out how much RFI is present and where.

Author: Jamie Stevens
email: Jamie.Stevens@csiro.au

"""

def USER_rfi_calculator(progressDict, options):
    USER_name = 'rfi_calculator'
    USER_version = '0.1'
    USER_print_name(USER_name, USER_version, options)

    # Check that we have the dictionary elements we require.
    if ('loadFlagStats' not in progressDict or
        'autoFlagStats' not in progressDict):
        # For some reason we haven't been automatically flagged.
        print("Unable to run RFI calculator. Exiting.")
        return USER_routine_unsuccessful(USER_name, USER_version, options)

    specUse = [ { 'lowFreq': 960, 'highFreq': 1164,
                  'usage': [ 'Aeronautical Radionavigation',
                             'Aeronautical Mobile' ] },
                { 'lowFreq': 1164, 'highFreq': 1215,
                  'usage': [ 'Aeronautical Radionavigation',
                             'Radionavigation - Satellite' ] },
                { 'lowFreq': 1215, 'highFreq': 1240,
                  'usage': [ 'Earth Exploration - Satellite',
                             'Radiolocation',
                             'Radionavigation - Satellite',
                             'Space Research' ] },
                { 'lowFreq': 1240, 'highFreq': 1300,
                  'usage': [ 'Earth Exploration - Satellite',
                             'Radiolocation',
                             'Radionavigation - Satellite',
                             'Space Research' ] },
                { 'lowFreq': 1300, 'highFreq': 1350,
                  'usage': [ 'Aeronautical Radionavigation',
                             'Radiolocation',
                             'Radionavigation - Satellite' ] },
                { 'lowFreq': 1350, 'highFreq': 1400,
                  'usage': [ 'Radiolocation' ] },
                { 'lowFreq': 1400, 'highFreq': 1427,
                  'usage': [ 'Passive' ] },
                { 'lowFreq': 1427, 'highFreq': 1429,
                  'usage': [ 'Space Operation',
                             'Fixed', 'Mobile' ] },
                { 'lowFreq': 1429, 'highFreq': 1452,
                  'usage': [ 'Fixed', 'Mobile' ] },
                { 'lowFreq': 1452, 'highFreq': 1492,
                  'usage': [ 'Broadcasting',
                             'Broadcasting - Satellite',
                             'Fixed', 'Mobile' ] },
                { 'lowFreq': 1492, 'highFreq': 1518,
                  'usage': [ 'Fixed', 'Mobile' ] },
                { 'lowFreq': 1518, 'highFreq': 1525,
                  'usage': [ 'Fixed', 'Mobile',
                             'Mobile - Satellite' ] },
                { 'lowFreq': 1525, 'highFreq': 1530,
                  'usage': [ 'Space Operation', 'Fixed',
                             'Mobile - Satellite' ] },
                { 'lowFreq': 1530, 'highFreq': 1535,
                  'usage': [ 'Space Operation',
                             'Mobile - Satellite' ] },
                { 'lowFreq': 1535, 'highFreq': 1559,
                  'usage': [ 'Mobile - Satellite' ] },
                { 'lowFreq': 1559, 'highFreq': 1610,
                  'usage': [ 'Aeronautical Radionavigation',
                             'Radionavigation - Satellite' ] },
                { 'lowFreq': 1610, 'highFreq': 1626.5,
                  'usage': [ 'Mobile - Satellite',
                             'Aeronautical Navigation',
                             'Radiodetermination - Satellite' ] },
                { 'lowFreq': 1626.5, 'highFreq': 1660.5,
                  'usage': [ 'Mobile - Satellite' ] },
                { 'lowFreq': 1660.5, 'highFreq': 1668,
                  'usage': [ 'Passive' ] },
                { 'lowFreq': 1668, 'highFreq': 1670,
                  'usage': [ 'Meteorological Aids',
                             'Fixed', 'Mobile',
                             'Mobile - Satellite' ] },
                { 'lowFreq': 1670, 'highFreq': 1675,
                  'usage': [ 'Meteorological Aids',
                             'Fixed', 'Meteorological - Satellite',
                             'Mobile', 'Mobile - Satellite' ] },
                { 'lowFreq': 1675, 'highFreq': 1690,
                  'usage': [ 'Meteorological Aids',
                             'Fixed', 'Meteorological - Satellite' ] },
                { 'lowFreq': 1690, 'highFreq': 1700,
                  'usage': [ 'Meteorological Aids',
                             'Meteorological - Satellite' ] },
                { 'lowFreq': 1700, 'highFreq': 1710,
                  'usage': [ 'Fixed', 'Meteorological - Satellite',
                             'Mobile' ] },
                { 'lowFreq': 1710, 'highFreq': 1980,
                  'usage': [ 'Fixed', 'Mobile' ] },
                { 'lowFreq': 1980, 'highFreq': 2010,
                  'usage': [ 'Fixed', 'Mobile', 'Mobile - Satellite' ] },
                { 'lowFreq': 2010, 'highFreq': 2025,
                  'usage': [ 'Fixed', 'Mobile' ] },
                { 'lowFreq': 2025, 'highFreq': 2110,
                  'usage': [ 'Space Operation',
                             'Earth Exploration - Satellite', 'Fixed',
                             'Mobile', 'Space Research' ] },
                { 'lowFreq': 2110, 'highFreq': 2120,
                  'usage': [ 'Fixed', 'Mobile', 'Space Research' ] },
                { 'lowFreq': 2120, 'highFreq': 2170,
                  'usage': [ 'Fixed', 'Mobile' ] },
                { 'lowFreq': 2170, 'highFreq': 2200,
                  'usage': [ 'Fixed', 'Mobile', 'Mobile - Satellite' ] },
                { 'lowFreq': 2200, 'highFreq': 2290,
                  'usage': [ 'Space Operation',
                             'Earth Exploration - Satellite', 'Fixed',
                             'Mobile', 'Space Research' ] },
                { 'lowFreq': 2290, 'highFreq': 2300,
                  'usage': [ 'Fixed', 'Mobile', 'Space Research' ] },
                { 'lowFreq': 2300, 'highFreq': 2483.5,
                  'usage': [ 'Fixed', 'Mobile', 'Radiolocation' ] },
                { 'lowFreq': 2483.5, 'highFreq': 2500,
                  'usage': [ 'Fixed', 'Mobile', 'Mobile - Satellite',
                             'Radiolocation',
                             'Radiodetermination - Satellite' ] },
                { 'lowFreq': 2500, 'highFreq': 2520,
                  'usage': [ 'Fixed', 'Fixed - Satellite',
                             'Mobile', 'Mobile - Satellite' ] },
                { 'lowFreq': 2520, 'highFreq': 2535,
                  'usage': [ 'Fixed', 'Fixed - Satellite',
                             'Mobile', 'Broadcasting - Satellite' ] },
                { 'lowFreq': 2535, 'highFreq': 2655,
                  'usage': [ 'Fixed', 'Mobile',
                             'Broadcasting - Satellite' ] },
                { 'lowFreq': 2655, 'highFreq': 2670,
                  'usage': [ 'Fixed', 'Fixed - Satellite',
                             'Mobile', 'Broadcasting - Satellite' ] },
                { 'lowFreq': 2670, 'highFreq': 2690,
                  'usage': [ 'Fixed', 'Fixed - Satellite',
                             'Mobile', 'Mobile - Satellite' ] },
                { 'lowFreq': 2690, 'highFreq': 2700,
                  'usage': [ 'Passive' ] },
                { 'lowFreq': 2700, 'highFreq': 2900,
                  'usage': [ 'Aeronautical Radionavigation',
                             'Radiolocation' ] },
                { 'lowFreq': 2900, 'highFreq': 3100,
                  'usage': [ 'Radiolocation', 'Radionavigation' ] },
                { 'lowFreq': 3100, 'highFreq': 3400,
                  'usage': [ 'Radiolocation' ] } ]
    
    # We examine the flagging statistics per channel and try to make
    # a heat map.
    progressDict['rfiCalculator'] = {}
    for n in range(0, len(progressDict['datasets'])):
        p = progressDict['datasets'][n]
        if 'autoFlagStats' not in progressDict:
            print('Dataset', p, 'was not flagged automatically.')
            continue
        if options['verbose']:
            print('Dataset', p)
        progressDict['logs'][n].write('Summarising the automatic flagging.\n')
        # We want the amount of flagging done by the automatic flagger.
        diffFlags = []
        for i in range(1, (len(progressDict['loadFlagStats'][p]['channel']) + 1)):
            q = str(i)
#            autoPct = float(re.sub('%', '',
#                                   progressDict['autoFlagStats'][p]['channel'][q]))
            autoPct = progressDict['autoFlagStats'][p]['channel'][q] * 100.0
#            loadPct = float(re.sub('%', '',
#                                   progressDict['loadFlagStats'][p]['channel'][q]))
            loadPct = progressDict['loadFlagStats'][p]['channel'][q] * 100.0
            if (loadPct == 100):
                diffFlags.append(loadPct)
            else:
                diffFlags.append(autoPct - loadPct)
        # Now get the frequency of each channel.
        diffFreqs = []
        for f in progressDict['freqConfigs']:
            c = progressDict['freqConfigs'][f]
            if p in c['dataset']:
                # Now find the index.
                j = c['dataset'].index(p)
                diffFreqs = c['chanFreqs'][j]

        # Determine the amount of flagging for each spectrum-use range.

        ruin = {}
        for u in range(0, len(specUse)):
            # The channels covering this frequency range.
            diffFreqsNp = np.array(diffFreqs)
            lf = specUse[u]['lowFreq'] / 1000.0
            hf = specUse[u]['highFreq'] / 1000.0
            rn = np.where(np.logical_and(diffFreqsNp >= lf,
                                         diffFreqsNp < hf))
            if (len(rn) < 1):
                continue
            # The average amount of flagging in this range.
            diffFlagsNp = np.array(diffFlags)
            fp = np.mean(diffFlagsNp[rn])
            # The bandwidth of this region.
            bw = max(diffFreqsNp[rn]) - min(diffFreqsNp[rn])
            # Now distribute equally amongst the culprits.
            v = (bw / len(specUse[u]['usage'])) * (fp / 100)
            for c in range(0, len(specUse[u]['usage'])):
                usageName = specUse[u]['usage'][c]
                if usageName in ruin:
                    ruin[usageName].append(v)
                else:
                    ruin[usageName] = [ v ]
        
        # Now summarise the ruiners.
        totBw = (max(diffFreqs) - min(diffFreqs)) * 1000.0
        totCost = 0.0
        ruinJSON = {}
        for r in ruin:
            freqCost = sum(ruin[r]) * 1000.0
            costPct = freqCost / totBw
            costStr = 'RFI generator %(r)s cost us %(freqCost).1f MHz (%(costPct).1f%%)' % \
                      {'r': r, 'freqCost': freqCost,
                       'costPct': (costPct * 100.0)}
            if options['verbose']:
                print(costStr)
            progressDict['logs'][n].write(costStr + '\n')
            totCost += freqCost
            ruinJSON[r] = { 'bandwidthCost': { 'absolute': freqCost,
                                               'fraction': costPct } }
        costStr = 'Overall, RFI cost us %(totCost).1f MHz (%(totPct).1f%%)' % \
                  {'totCost': totCost, 'totPct': (totCost / totBw * 100.0)}
        if options['verbose']:
            print(costStr)
        progressDict['logs'][n].write(costStr + '\n')
        # Dump the JSON string to a file.
        outFile = 'rfi_calculator.' + p + '.json'
        if options['verbose']:
            print('Writing summary to file', outFile)
        with open(outFile, 'w') as f:
            f.write(json.dumps(ruinJSON))
        # Pass this information back to the rest of the pipeline too.
        progressDict['rfiCalculator'][p] = ruinJSON
        
        # Make a line plot.
#        plt.fill_between(diffFreqs, 0, diffFlags, facecolor='black')
#        # Make a histogram plot.
#        width = 100.0
#        pos = np.arange(0, len(diffFlags), width)
#        ax = plt.axes()
#        ax.set_xticks(pos + (width / 2))
#        ax.set_xticklabels(diffFreqs)
#        plt.bar(pos, diffFlags, width, color='r')
#        plt.show()

    # Simple example first: flagging comparison statistics.
#    for n in range(0, len(progressDict['datasets'])):
#        p = progressDict['datasets'][n]
#        for b in progressDict['loadFlagStats'][p]['baseline']:
#            print('Flagging on baseline', b, 'before:',
#                  progressDict['loadFlagStats'][p]['baseline'][b],
#                  'after:',
#                  progressDict['autoFlagStats'][p]['baseline'][b])

    return USER_routine_successful(USER_name, USER_version, options)

