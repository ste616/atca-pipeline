from __future__ import print_function
import sys
from optparse import OptionParser
sys.path.append('/n/ste616/atcacode/cabb_pipeline')
sys.path.append('/home/jstevens/usr/src/atcacode/cabb_pipeline')
from cabb_pipeline_main_modules import *
from cabb_pipeline_rfi_calculator import USER_rfi_calculator

version = '0.1'

#print('CABB pipeline v' + version)

# Go through the arguments.
usage = "usage: %prog [options] rpfits-file [rpfits-file ...]"
parser = OptionParser(usage=usage, version="CABB pipeline " + version)
parser.add_option('--no-atlod',
                  help='disable loading step, if data has already been loaded',
                  action='store_true')
parser.add_option('--no-split',
                  help='disable split step, if data has already been split' +
                  ' (cannot be used without --no-atlod)',
                  action='store_true')
parser.add_option('--no-midweek', action='store_true',
                  help='disable midweek RFI detection step')
parser.add_option('--no-flag',
                  help='disable automatic flagging step',
                  action='store_true')
parser.add_option('--keep-flags',
                  help='use the current flagging table to begin with',
                  action='store_true')
parser.add_option('--use-flags',
                  help='specify the flag table to use upon starting' +
                  ' (only useful when not loading and/or splitting)',
                  default='original')
parser.add_option('--no-user',
                  help='disable running of the user routines',
                  action='store_true')
parser.add_option('--no-casa', action='store_true',
                  help='disable all parts that depend on casapy')
parser.add_option('--keep-reduction', action='store_true',
                  help='keep the current reduction split directory, if it ' +
                  'exists, but delete the current calibration tables')
parser.add_option('-v', '--verbose', action='store_true', dest='verbose',
                  help='allow all output from routines')
parser.add_option('-q', '--quiet', action='store_true', dest='quiet',
                  help='supress all output from pipeline, except for errors')
(options, args) = parser.parse_args()

# Have to change options into ioptions so we can use the script interactively.
ioptions = {}
ioptions['no_atlod'] = options.no_atlod
ioptions['no_split'] = options.no_split
ioptions['no_midweek'] = options.no_midweek
ioptions['no_flag'] = options.no_flag
ioptions['keep_flags'] = options.keep_flags
ioptions['use_flags'] = options.use_flags
ioptions['no_user'] = options.no_user
ioptions['no_casa'] = options.no_casa
ioptions['keep_reduction'] = options.keep_reduction
ioptions['verbose'] = options.verbose
ioptions['quiet'] = options.quiet

# Check for incompatible arguments.
if ioptions['no_split'] and not ioptions['no_atlod']:
    # Can't not split a newly loaded dataset.
    ioptions['no_split'] = False
if ioptions['verbose'] and ioptions['quiet']:
    # Can't be quiet and verbose: we choose verbose.
    ioptions['quiet'] = False
if ((not ioptions['no_atlod'] or not ioptions['no_split']) and ioptions['keep_flags']):
    # Can't keep flags for a newly split or loaded set.
    ioptions['keep_flags'] = False
# And arguments that won't work.
if not cmd_exists('uv2ms'):
    # Won't be able to convert to MeasurementSet format.
    ioptions['no_casa'] = True

# The positional arguments to the script are the RPFITS files to load.
rpfitsFiles = args

if len(rpfitsFiles) < 1:
    print('Please supply the RPFITS files to load as command line arguments.')
    parser.print_usage()
    sys.exit(0)

# Load the files and start the results dictionary.
progressDict = cabbLoad(rpfitsFiles, ioptions)
if ioptions['verbose']:
    print('Miriad data file:', progressDict['miriadData'])

# Split the master file into its component IFs.
progressDict['datasets'] = splitIFs(progressDict['miriadData'],
                                    progressDict['freqConfigs'], ioptions)

# Determine the flagging fraction of these data sets.
progressDict['loadFlagStats'] = {}
origLoad = {}
for n in range(0, len(progressDict['datasets'])):
    p = progressDict['datasets'][n]
    # Load the original flags for this step.
    origLoad[p] = False
    if checkMiriadFlagTable(p, 'original', ioptions):
        origLoad[p] = True
        keepMiriadFlagTable(p, 'pipelineStart', ioptions)
        restoreMiriadFlagTable(p, 'original', ioptions)
    progressDict['loadFlagStats'][p] = flagStats(p, ioptions)
    if origLoad[p]:
        restoreMiriadFlagTable(p, 'pipelineStart', ioptions)

if ioptions['no_split'] and not ioptions['keep_flags']:
    # Restore the user-specified flag table if it is present.
    for n in range(0, len(progressDict['datasets'])):
        restoreMiriadFlagTable(progressDict['datasets'][n],
                               ioptions['use_flags'], ioptions)

# Output a summary of the dataset and determine future actions.
masterLog = open('log.' + progressDict['miriadData'], 'w')
masterLog.write('CABB pipeline v' + version + '\n')
masterLog.write('Files loaded:\n')
for i in range(0, len(progressDict['rpfitsFiles'])):
    masterLog.write(progressDict['rpfitsFiles'][i] + '\n')
nConfigs = 0
for fc in progressDict['freqConfigs']:
    nConfigs += 1
masterLog.write('There are ' + str(nConfigs) + ' frequency configs present.\n')
for fc in progressDict['freqConfigs']:
    masterLog.write('Frequency Configuration ' + fc + ':\n')
    # Output the IF configuration.
    sfc = progressDict['freqConfigs'][fc]
    for i in range(1, 3):
        masterLog.write('IF ' + str(i) + ':\n')
        someOutput = False
        for j in range(0, len(sfc['ifChain'])):
            if sfc['ifChain'][j] == i:
                someOutput = True
                chainStr = (sfc['classification'][j]['bandType'] + ' band: ' +
                            'BW=' + str(sfc['classification'][j]['bandwidth']) +
                            ' [GHz] NCHAN=' + str(sfc['channels'][j]) +
                            ' CFREQ=' + str(sfc['centreFreq'][j]) + ' [GHz]' +
                            ' CWIDTH=' + str(sfc['chanWidth'][j]) + ' [GHz]' +
                            ' SIDE=' + sfc['sideband'][j])
                if sfc['classification'][j]['bandType'] == 'zoom':
                    chainStr += (' ZTYPE=' + 
                                 sfc['classification'][j]['zoomType'])
                masterLog.write(chainStr + '\n')
        if not someOutput:
            masterLog.write('None')

# Open some logs for each dataset.
progressDict['logs'] = []
for l in range(0, len(progressDict['datasets'])):
    progressDict['logs'].append(open('log.' + progressDict['datasets'][l], 'w'))
                
"""
The progress dictionary layout:
progressDict
 -> 'rpfitsFiles': the list of RPFITS files fed into the pipeline
 -> 'miriadData': the name of the Miriad dataset output from atlod
 -> 'sources': each key is the name of a source
   -> 'rightAscension': the R.A. of the named source
   -> 'declination': the Dec of the named source
   -> 'freqConfig': an array of the frequency config each time this source
                    was observed
   -> 'startTime': an array of the start time each time this source was observed
   -> 'calCode': an array of the calcode each time this source was observed
 -> 'freqConfigs': each key is the frequency config from uvindex
   -> 'channels': an array (each index is an IF); number of channels
   -> 'firstFreq': an array (each index is an IF); frequency of channel 1 (GHz)
   -> 'chanWidth': an array (each index is an IF); the width of each channel (GHz)
   -> 'restFreq': an array (each index is an IF); the supplied rest freq (GHz)
   -> 'centreFreq': an array (each index is an IF); the band centre freq (GHz)
   -> 'sideband': an array (each index is an IF); the sideband (USB or LSB)
   -> 'dataset': an array (each index is an IF); the name of the split dataset
   -> 'ifChain': an array (each index is an IF); the IF chain for this IF; this
                 parameter will only usually be available if zooms are specified
   -> 'chanFreqs': an array (each index is an IF); an array of length 'channels',
                   where each value is the freq of the channel (GHz)
   -> 'classification': an array (each index is an IF); the type of IF, as a
                        dictionary
     -> 'bandwidth': the total bandwidth of this IF (GHz)
     -> 'bandType': the type of band ('wide', 'zoom')
     -> 'zoomType': if a zoom band, the type of zoom ('single', 'consolidated')
     -> 'zoomConfig': if a zoom band, the zoom config ('1', '4', '16', '64')
 -> 'datasets': an array listing all the datasets from the main file
 -> 'logs': an array, a log for each of the datasets
 -> 'measurementSets': an array listing all the MS from the main file (same length
                       and direct correspondence to 'datasets')
 -> 'calibrationSources': the calibration to use for each dataset, and the
                          dataset name is the key for this dictionary
   -> 'bandpass': the name of the bandpass calibrator to use
   -> 'flux': the name of the flux density calibrator to use
   -> 'leakages': the name of the leakage calibrator to use
   -> 'gain': an array of gain calibrator names
   -> 'warnings': an array of text warnings about the potential problems with
                  the calibration
 -> 'reductionDir': the directory being used for the Miriad reduction, and the
                    dataset name is the key for this dictionary
 -> 'loadFlagStats': the flag statistics immediately after atlod and uvsplit;
                     each key is the name of a dataset
   -> 'stokes': each key is the name of a Stokes parameter, each value
                is the percentage of visibilities flagged
   -> 'baseline': each key is the name of a baseline, each value is the
                  percentage of visibilities flagged
   -> 'antenna': each key is the name of an antenna, each value is the
                 percentage of visibilities flagged
   -> 'channel': each key is the number of a channel, each value is the
                 percentage of visibilities flagged
 -> 'startFlagStats': the flag statistics upon starting this run of the pipeline;
                      the layout is the same as for 'loadFlagStats'
 -> 'autoFlagStats': the flag statistics after automatic flagging; the layout
                     is the same as for 'loadFlagStats'
     
"""

# Keep a copy of this original flagging table.
progressDict['startFlagStats'] = {}
for n in range(0, len(progressDict['datasets'])):
    p = progressDict['datasets'][n]
    if ioptions['keep_flags']:
        keepMiriadFlagTable(p, 'user', ioptions)
    else:
        keepMiriadFlagTable(p, 'original', ioptions)
    # Get the flagging statistics again.
    if origLoad[p]:
        progressDict['startFlagStats'][p] = flagStats(p, ioptions)
    else:
        progressDict['startFlagStats'][p] = progressDict['loadFlagStats'][p]

# Check for mid-week RFI.
progressDict['midweekRFI'] = {}
progressDict['midweekStats'] = {}
for n in range(0, len(progressDict['datasets'])):
    p = progressDict['datasets'][n]
    # Check that we can actually detect midweek RFI in this IF.
    t = getIF(progressDict, p)
    if (t is None or t['classification']['bandType'] != 'wide' or
        frequencyBand(t['centerFreq']) != '16cm'):
        # Can't flag based on this IF.
        continue
    didMidweek = False
    d = chainDatasets(progressDict, t['freqConfig'], t['ifChain'])
    if not ioptions['no_midweek']:
        progressDict['logs'][n].write('Checking this dataset for midweek RFI.\n')
        progressDict['midweekRFI'][p] = midweekDetector(p, ioptions)
        if midweekFlagger(p, progressDict['midweekRFI'][p]['flagRegions'], ioptions):
            didMidweek = True
            progressDict['logs'][n].write('Detected and flagged midweek RFI in the ' +
                                          'time ranges:\n')
            for i in range(0, len(progressDict['midweekRFI'][p]['flagRegions'])):
                progressDict['logs'][n].write(
                    progressDict['midweekRFI'][p]['flagRegions'][i]['start'] + ' - ' +
                    progressDict['midweekRFI'][p]['flagRegions'][i]['stop'] + '\n')
            keepMiriadFlagTable(p, 'midweek', ioptions)
            for i in range(0, len(d)):
                if d[i] == p:
                    continue
                midweekFlagger(d[i], progressDict['midweekRFI'][p]['flagRegions'],
                               ioptions)
                keepMiriadFlagTable(d[i], 'midweek', ioptions)
        else:
            progressDict['logs'][n].write('No midweek RFI detected.\n')
    elif checkMiriadFlagTable(p, 'midweek', ioptions):
        for i in range(0, len(d)):
            restoreMiriadFlagTable(d[i], 'midweek', ioptions)
        didMidweek = True
    # New flagging statistics.
    for i in range(0, len(d)):
        if didMidweek:
            progressDict['midweekStats'][d[i]] = flagStats(d[i], ioptions)
        else:
            progressDict['midweekStats'][d[i]] = progressDict['startFlagStats'][d[i]]

progressDict['autoFlagStats'] = {}
for n in range(0, len(progressDict['datasets'])):
    p = progressDict['datasets'][n]
    autoRestored = False
    # Check for zoom bands.
    t = getIF(progressDict, p)
    if (t is None or t['classification']['bandType'] != 'wide'):
        # Don't flag zoom bands.
        continue
    if not ioptions['no_flag']:
        autoPgflag(p, ioptions)
        # Keep a copy of this flagging.
        keepMiriadFlagTable(p, 'auto', ioptions)
    else:
        # Check if we should restore a previous automatic flag set.
        if checkMiriadFlagTable(p, 'auto', ioptions):
            autoRestored = True
            keepMiriadFlagTable(p, 'pipelineAutoSwitch', ioptions)
            restoreMiriadFlagTable(p, 'auto', ioptions)

    # Determine the flagging fraction now.
    progressDict['autoFlagStats'][p] = flagStats(p, ioptions)
    if autoRestored:
        restoreMiriadFlagTable(p, 'pipelineAutoSwitch', ioptions)

# Turn the dataset into a measurement set.
progressDict['measurementSets'] = []
if not ioptions['no_casa']:
    for n in range(0, len(progressDict['datasets'])):
        progressDict['measurementSets'].append(
            uvToMeasurementSet(progressDict['datasets'][n], ioptions))

# USER routines go here.
if not ioptions['no_user']:
    # RFI detector: Jamie Stevens
    # Adds to progress dictionary:
    # -> 'rfiCalculator': a summary of the emitters that have been flagged,
    #                     both as an absolute amount of bandwidth, and as the
    #                     fraction of the total. Each key is the name of an
    #                     emitter, from the ACMA spectrum plan.
    #   -> 'bandwidthCost'
    #     -> 'absolute': the amount of bandwidth flagged (MHz)
    #     -> 'fraction': the fraction of the total bandwidth flagged
    USER_rfi_calculator(progressDict, ioptions)
    
# Identify the calibrators.
fcSources = {}
for c in progressDict['freqConfigs']:
    fcSources[c] = determineCalibrators(progressDict, c, ioptions)
#    if ioptions['verbose']:
#        print('Calibration sources for frequency config', c)
#        print(fcSources[c])

progressDict['calibrationSources'] = {}
for n in range(0, len(progressDict['datasets'])):
    p = progressDict['datasets'][n]
    for c in progressDict['freqConfigs']:
        if p in progressDict['freqConfigs'][c]['dataset']:
            progressDict['calibrationSources'][p] = fcSources[c]
            progressDict['logs'][n].write('Calibration sources determined:\n')
            progressDict['logs'][n].write('Bandpass: ' +
                                          fcSources[c]['bandpass'] + '\n')
            progressDict['logs'][n].write('Flux Density: ' +
                                          fcSources[c]['flux'] + '\n')
            progressDict['logs'][n].write('Leakages: ' +
                                          fcSources[c]['leakages'] + '\n')
            progressDict['logs'][n].write('Gain: ' +
                                          ", ".join(fcSources[c]['gain']) + '\n')
            break

# Perform calibration.
progressDict['reductionDir'] = {}
progressDict['refAnt'] = {}
for n in range(0, len(progressDict['datasets'])):
    p = progressDict['datasets'][n]
    ourLog = progressDict['logs'][n]
    reflagged = {}
    # Make a directory for the overall reduction and split this dataset
    # into it.
    progressDict['reductionDir'][p] = prepareReductionDir(p, ioptions)
    reductionDir = progressDict['reductionDir'][p]

    # Change into this directory now.
    if not chReductionDir(reductionDir, ioptions):
        sys.exit(0)

    # Determine the reference antenna we will use.
    if not ioptions['no_flag']:
        progressDict['refAnt'][p] = determineRefAnt(
            progressDict['autoFlagStats'][p], ioptions)
    else:
        progressDict['refAnt'][p] = determineRefAnt(
            progressDict['startFlagStats'][p], ioptions)
    refant = progressDict['refAnt'][p]

    # Begin with bandpass calibration.
    bpcal = progressDict['calibrationSources'][p]['bandpass']
    if bpcal is None:
        print('Unable to bandpass calibrate dataset', p)
        sys.exit(0)

    bpset = findDataset(bpcal)
    if bpset is None:
        print('Cannot find bandpass calibrator dataset.')
        sys.exit(0)

    bandpassState = calibrateBandpass(bpset, refant, ioptions)
    if not bandpassState['converged']:
        if not ioptions['quiet']:
            print('Bandpass calibration with ' + bpcal +
                  ' failed to converge!')
        ourLog.write('Bandpass calibration did not ' +
                                      'converge.\n')
    else:
        ourLog.write('Bandpass calibration converged ' +
                     'in ' + str(bandpassState['iterations']) +
                     ' iterations.\n')

    # Flag the bandpass calibrator again.
    reflagged[bpcal] = autoPgflag(bpset, ioptions)

    # And redo the bandpass calibration.
    bandpassState = calibrateBandpass(bpset, refant, ioptions)
    if not bandpassState['converged']:
        if not ioptions['quiet']:
            print('Bandpass calibration with ' + bpcal +
                  ' failed to converge!')
        ourLog.write('Bandpass calibration did not ' +
                                      'converge.\n')
    else:
        ourLog.write('Bandpass calibration converged ' +
                     'in ' + str(bandpassState['iterations']) +
                     ' iterations.\n')

    # Fix the bandpass solution if required.
    fluxcal = progressDict['calibrationSources'][p]['flux']
    if fluxcal is None:
      print('Unable to flux calibrate dataset', p)
    else:
        fluxset = findDataset(fluxcal)
        if fluxset is None:
            print('Cannot find flux calibrator dataset.')
            sys.exit(0)
        if bpcal != fluxcal:
            # Copy the bandpass table to the flux calibrator.
            copyCalTables(bpset, fluxset, {}, ioptions)
            # If required, flag the flux calibrator.
            if (fluxcal not in reflagged or
                reflagged[fluxcal] == False):
                reflagged[fluxcal] = autoPgflag(fluxset, ioptions)
            # Correct the bandpass table.
            correctBandpass(fluxcal, fluxset, ioptions)
            # Copy the bandpass table back to the bandpass calibrator.
            copyCalTables(fluxset, bpset, { 'pol': False, 'cal': False,
                                            'pass': True }, ioptions)

    # Calibrate the leakages.
    leakagecal = progressDict['calibrationSources'][p]['leakages']
    if leakagecal is None:
        print('Unable to calibrate leakages for dataset', p)
    else:
        leakageset = findDataset(leakagecal)
        if leakageset is None:
            print('Cannot find leakage calibrator dataset.')
            sys.exit(0)
        if leakagecal != bpcal:
            # Copy the bandpass table to the leakage calibrator.
            copyCalTables(bpset, leakageset, {}, ioptions)
        # If required, flag the leakage calibrator.
        if (leakagecal not in reflagged or
            reflagged[leakagecal] == False):
            reflagged[leakagecal] = autoPgflag(leakageset, ioptions)
        # Do the gain calibration.
        if leakagecal != '1934-638':
            leakageState = calibrateGains(leakageset, refant, {}, ioptions)
        else:
            leakageState = calibrateGains(leakageset, refant,
                                          { 'leakages': True }, ioptions)
        if leakagecal != bpcal:
            # Transfer the gains back to the bandpass calibrator.
            copyCalTables(leakageset, bpset, { 'cal': False, 'pass': False },
                          ioptions)

    # Calibrate the flux density calibrator.
    

    # Go back to the previous directory.
    if reductionDir != '.':
        if not chReductionDir('..', options):
            sys.exit(0)

# Wrap it all up.
masterLog.close()
for n in range(0, len(progressDict['datasets'])):
    progressDict['logs'][n].close()
if not ioptions['quiet']:
    print('Pipeline operation completed.')

