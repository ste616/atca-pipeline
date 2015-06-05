from __future__ import print_function
import sys
import os
import re
import math
from mirpy import miriad
import numpy as np
import shutil
import glob
from datetime import datetime, timedelta
import calendar
import time
import subprocess
import atca_calibrator_database as caldb

# A collection of routines used by the system part of the
# CABB pipeline. Users should not need to alter this file.

# The following routines are algorithmic only: they don't print
# output to the screen.

def lineFilter(output):
    # Filter output from a Miriad task line-by-line.
    return output.splitlines()

def num(s):
    # Try to turn some string into a number.
    try:
        return int(s)
    except ValueError:
        return float(s)

def frequencyArray(firstFreq, chanWidth, channels):
    # Turn an IF specification into an array of channel frequencies.
    retArr = [ float(firstFreq) ]
    for i in range(2, (int(channels) +1)):
        retArr.append(float(i - 1) * float(chanWidth) + float(firstFreq))
    return retArr

def characteriseIF(channels, firstFreq, chanWidth):
    # Figure out what kind of IF this is.
    # Easiest way is to use the channel width, but we need to adjust
    # the channel width for 1 MHz zooms first.
    fChanWidth = float(chanWidth)
    fChannels = int(channels)
    fFirstFreq = float(firstFreq)
    if math.fabs(fChanWidth) < 5e-7:
        # This is a 1 MHz zoom width, which is actually 488 Hz.
        t = 1e-6 / 2048.0
        if fChanWidth < 0:
            t *= -1.0
        fChanWidth = t
    rDict = {}
    rDict['bandwidth'] = math.fabs(fChanWidth * fChannels)
    if rDict['bandwidth'] > 2:
        # A wide band.
        rDict['bandType'] = 'wide'
    elif fChannels >= 2049:
        # Must be a zoom band.
        rDict['bandType'] = 'zoom'
        if fChannels > 2049:
            # A consolidated zoom.
            rDict['zoomType'] = 'consolidated'
        else:
            rDict['zoomType'] = 'single'
        if math.fabs(fChanWidth) < 5e-7:
            rDict['zoomConfig'] = '1'
        elif math.fabs(fChanWidth) < 2e-6:
            rDict['zoomConfig'] = '4'
        elif math.fabs(fChanWidth) < 8e-6:
            rDict['zoomConfig'] = '16'
        else:
            rDict['zoomConfig'] = '64'
    return rDict

def getIF(progressDict, dataset):
    # Return the IF this dataset is from.
    for f in progressDict['freqConfigs']:
        for i in range(0, len(progressDict['freqConfigs'][f]['dataset'])):
            if dataset == progressDict['freqConfigs'][f]['dataset'][i]:
                rDict = {
                    'freqConfig': f,
                    'channels': progressDict['freqConfigs'][f]['channels'][i],
                    'firstFreq': progressDict['freqConfigs'][f]['firstFreq'][i],
                    'chanWidth': progressDict['freqConfigs'][f]['chanWidth'][i],
                    'restFreq': progressDict['freqConfigs'][f]['restFreq'][i],
                    'centerFreq': progressDict['freqConfigs'][f]['centreFreq'][i],
                    'chanFreqs': progressDict['freqConfigs'][f]['chanFreqs'][i],
                    'classification': progressDict['freqConfigs'][f]['classification'][i],
                    'ifChain': progressDict['freqConfigs'][f]['ifChain'][i] }
                return rDict
                
    # We should never get here.
    return None

def frequencyBand(frequency):
    # Return the band name that this frequency is in.
    if frequency < 3.5:
        return '16cm'
    elif frequency < 12:
        return '4cm'
    elif frequency < 28:
        return '15mm'
    elif frequency < 60:
        return '7mm'
    else:
        return '3mm'

def chainDatasets(progressDict, freqConfig, ifChain):
    # Return an array of all the datasets in this IF chain.
    rArr = []
    for i in range(0, len(progressDict['freqConfigs'][freqConfig]['ifChain'])):
        if progressDict['freqConfigs'][freqConfig]['ifChain'][i] == ifChain:
            rArr.append(progressDict['freqConfigs'][freqConfig]['dataset'][i])
    return rArr
    
def interactiveMode():
    # Make an options object with defaults set for use in interactive mode.
    options = { 'no_atlod': True, 'no_split': True, 'no_flag': True,
                'use_flags': 'original', 'no_user': False, 'keep_flags': False,
                'no_casa': False, 'keep_reduction': False,
                'verbose': True, 'quiet': False }
    return options

def cmd_exists(cmd):
    # Check for the command in the user's PATH.
    return subprocess.call("type " + cmd, shell=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

def findDataset(prefix):
    # Find a dataset that starts with the prefix.
    po = glob.glob(prefix + '.*')
    if len(po) == 1:
        return po[0]
    # We shouldn't get in this situation if we've done things right.
    return None

# The following routines are for process control and may print
# output to the screen.

def cabbLoad(rpfitsFiles, options):
    # Load a set of RPFITS files into a Miriad data set.
    if not options['quiet']:
        print('Loading RPFITS files...', ", ".join(rpfitsFiles))
    
    # Check first that each RPFITS file is accessible to us.
    for f in rpfitsFiles:
        if not os.path.isfile(f):
            print('File', f, 'is not accessible.')
            sys.exit(0)
    # Form the name of the output file from the first input.
    outEls = re.split('^(....-..-..).*\.(.*)$', rpfitsFiles[0])
    rDict = {}
    rDict['miriadData'] = outEls[2] + '_' + outEls[1] + '.uv'
    # Store the RPFITS files as a list as well.
    rDict['rpfitsFiles'] = rpfitsFiles
    # Check if this dataset exists, and remove it if it does.
    if options['no_atlod']:
        if options['verbose']:
            print("Not deleting exisiting data set.")
    else:
        if os.path.isdir(rDict['miriadData']):
            try:
                shutil.rmtree(rDict['miriadData'])
            except:
                print('Cannot delete existing data file.')
                sys.exit(0)
            
    # Run atlod.
    miriad.set_filter('atlod', lineFilter)
    loadOptions = ['birdie', 'rfiflag', 'xycorr', 'noauto']
    if options['no_atlod']:
        if options['verbose']:
            print("Not performing load operation.")
    else:
        r = miriad.atlod(In=",".join(rpfitsFiles), out=rDict['miriadData'],
                         options=",".join(loadOptions))
    rDict['sources'] = {}
    if not options['quiet']:
        print('Loading complete.')

    # Get the sources and frequencies from running uvindex.
    if not options['quiet']:
        print('Compiling sources and frequencies...')
    miriad.set_filter('uvindex', lineFilter)
    ur = miriad.uvindex(vis=rDict['miriadData'])
    cc = 0
    cs = 0
    rDict['freqConfigs'] = {}
    for l in ur:
        lsp = re.split('\s+', l)
        if (len(lsp) == 8 and
            re.search('^\d\d\D\D\D\d\d\:\d\d\:\d\d\:\d\d\.\d$', lsp[0]) is not None):
            if lsp[1] not in rDict['sources']:
                rDict['sources'][lsp[1]] = { 'freqConfig': [ lsp[6] ],
                                             'startTime': [ lsp[0] ],
                                             'calCode': [ lsp[2] ] }
            else:
                rDict['sources'][lsp[1]]['freqConfig'].append(lsp[6])
                rDict['sources'][lsp[1]]['startTime'].append(lsp[0])
                rDict['sources'][lsp[1]]['calCode'].append(lsp[2])
        elif len(lsp) == 6:
            if lsp[0] in rDict['sources']:
                rDict['sources'][lsp[0]]['rightAscension'] = lsp[2]
                rDict['sources'][lsp[0]]['declination'] = lsp[3]
        elif lsp[0] == 'Frequency' and lsp[1] == 'Configuration':
            cc = lsp[2]
            cs = -1
        elif cs == 1:
            cs = -1
        elif cs == -1:
            centreFreq = (num(lsp[1]) - 1) / 2 * num(lsp[3]) + num(lsp[2])
            sideband = 'USB'
            if num(lsp[2]) > centreFreq:
                sideband = 'LSB'
            rDict['freqConfigs'][cc] = {
                'channels': [ num(lsp[1]) ],
                'firstFreq': [ num(lsp[2]) ],
                'chanWidth': [ num(lsp[3]) ],
                'restFreq': [ num(lsp[4]) ],
                'centreFreq': [ centreFreq ],
                'sideband': [ sideband ],
                'dataset': [ '' ],
                'chanFreqs': [ frequencyArray(lsp[2], lsp[3], lsp[1]) ],
                'classification': [ characteriseIF(lsp[1], lsp[2], lsp[3]) ] }
            if len(lsp) == 7:
                # We have zooms here probably.
                rDict['freqConfigs'][cc]['ifChain'] = [ num(lsp[6]) ]
            else:
                rDict['freqConfigs'][cc]['ifChain'] = [ 1 ]
            cs = -2
        elif cs == -2:
            if len(lsp) == 1:
                cs = 0
            else:
                centreFreq = (num(lsp[1]) - 1) / 2 * num(lsp[3]) + num(lsp[2])
                sideband = 'USB'
                if num(lsp[2]) > centreFreq:
                    sideband = 'LSB'
                rDict['freqConfigs'][cc]['channels'].append(num(lsp[1]))
                rDict['freqConfigs'][cc]['firstFreq'].append(num(lsp[2]))
                rDict['freqConfigs'][cc]['chanWidth'].append(num(lsp[3]))
                rDict['freqConfigs'][cc]['restFreq'].append(num(lsp[4]))
                rDict['freqConfigs'][cc]['centreFreq'].append(centreFreq)
                rDict['freqConfigs'][cc]['sideband'].append(sideband)
                rDict['freqConfigs'][cc]['dataset'].append('')
                rDict['freqConfigs'][cc]['chanFreqs'].append(
                    frequencyArray(lsp[2], lsp[3], lsp[1]))
                rDict['freqConfigs'][cc]['classification'].append(
                    characteriseIF(lsp[1], lsp[2], lsp[3]))
                if len(lsp) == 7:
                    rDict['freqConfigs'][cc]['ifChain'].append(num(lsp[6]))
                else:
                    nn = len(rDict['freqConfigs'][cc]['ifChain']) + 1
                    rDict['freqConfigs'][cc]['ifChain'].append(nn)
    if not options['quiet']:
        print('Source and frequency compilation complete.')
        print('Load process complete.')
    return rDict

def uvPresent(uvFile, options):
    # Check that the dataset directory is accessible to us.
    if not os.path.isdir(uvFile):
        if options['verbose']:
            print('Dataset', uvFile, 'is not accessible.')
        return False
    return True

def chReductionDir(redDir, options):
    # Try to change directory.
    try:
        os.chdir(redDir)
    except:
        if options['verbose']:
            print('Unable to change directory to', redDir)
        return False
    return True


def splitIFs(uvFile, freqConfigs, options):
    # Split out the component IFs from the master Miriad uv file.
    if not options['quiet']:
        print('Splitting out IFs from dataset', uvFile)

    # Check that the file is accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)

    # Check that there are no uvsplit directories already,
    # excluding MeasurementSets.
    uvsplits = [f for f in glob.glob('uvsplit.*') if '.ms' not in f and
                '.def' not in f]
    if len(uvsplits) > 0 and not options['no_split']:
        for u in uvsplits:
            # Delete this tree.
            try:
                shutil.rmtree(u)
            except:
                print('Cannot delete uvsplit set', u)
                sys.exit(0)

    # Split the set now.
    ur = []
    if options['no_split']:
        ur = uvsplits
    else:
        miriad.set_filter('uvsplit', lineFilter)
        ur = miriad.uvsplit(vis=uvFile, options='nosource')
    dsets = []
    for l in ur:
        lsp = re.split('\s+', l)
        if options['no_split']:
            lsp.insert(0, 'Creating')
        if lsp[0] == 'Creating':
            dsets.append(lsp[1])
            p = re.split('^uvsplit.([^\.]+)', lsp[1])
            # Search for the matching frequency config.
            for c in freqConfigs:
                for f in range(0, len(freqConfigs[c]['centreFreq'])):
                    m = int(round(freqConfigs[c]['centreFreq'][f] * 1000))
                    if m == int(p[1]) and freqConfigs[c]['dataset'][f] == '':
                        freqConfigs[c]['dataset'][f] = lsp[1]
                        break

    if not options['quiet']:
        print('Dataset splitting complete.')
    return dsets

def flagStats(uvFile, options):
    # Determine the amount of data flagged in a uv dataset.
    if not options['quiet']:
        print('Checking flagging in dataset', uvFile)

    # Check that the file is accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)

    # Run uvfstats.
    rDict = { 'stokes': {},
              'baseline': {},
              'antenna': {},
              'channel': {} }
    # Flagging stats by the parameters.
    miriad.set_filter('uvfstats', lineFilter)
    for mode in rDict:
        fr = miriad.uvfstats(vis=uvFile, mode=mode)
        cc = 0
        for l in fr:
            lsp = re.split('\s+', l)
            if len(lsp) > 1 and re.match('---', lsp[1]) is not None:
                cc = 1
            elif len(lsp) == 3 and cc == 1:
                rDict[mode][lsp[1]] = float(
                    re.sub('%', '', lsp[2])) / 100.0

    if not options['quiet']:
        print('Flagging check complete.')
    return rDict

def copySet(uvFile, suffix, options):
    # Make a new uvFile that is a copy of the other file.
    if not options['quiet']:
        print('Making a copy of dataset', uvFile)
    prefix = re.split('^(.*)\.uv$', uvFile)
    destFile = prefix[1] + '.' + suffix + '.uv'

    # Check if the file is accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)
        
    # Check if the new directory is already present.
    if os.path.isdir(destFile):
        try:
            shutil.rmtree(destFile)
        except:
            print('Cannot delete exisiting directory', destFile)
            sys.exit(0)
    # Do the copy.
    shutil.copytree(uvFile, destFile)
    if not options['quiet']:
        print('Copy to', destFile, 'complete.')
    return destFile

def keepMiriadFlagTable(uvFile, suffix, options):
    # Copy the flag table in the Miriad dataset to a new file.
    newTable = 'flags.' + suffix
    if options['verbose']:
        print('Adding new flag version table', newTable,
              'to the dataset', uvFile)

    # Check if the dataset and file are accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)
    currTable = uvFile + '/flags'
    outTable = uvFile + '/' + newTable
    if not os.path.isfile(currTable):
        print('Flag table', currTable, 'is not accessible.')
        sys.exit(0)

    # Do the copy.
    try:
        shutil.copyfile(currTable, outTable)
    except:
        print('Unable to make copy of flag table.')
        sys.exit(0)
    if options['verbose']:
        print('Flag version table created.')
    return True

def restoreMiriadFlagTable(uvFile, suffix, options):
    # Restore a stored flag table in the Miriad dataset.
    newTable = 'flags.' + suffix
    if options['verbose']:
        print('Restoring the flag version table', newTable,
              'in the dataset', uvFile)

    # Check if the dataset and file are accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)
    restTable = uvFile + '/' + newTable
    defaultTable = uvFile + '/flags'
    if os.path.isfile(restTable):
        try:
            shutil.copyfile(restTable, defaultTable)
        except:
            print('Unable to restor flag table.')
            sys.exit(0)
    if options['verbose']:
        print('Flag version table restored.')
    return True

def checkMiriadFlagTable(uvFile, suffix, options):
    # Check if a named flag table exists in the Miriad dataset.
    checkTable = 'flags.' + suffix
    if options['verbose']:
        print('Checking for the flag version table', checkTable,
              'in the dataset', uvFile)

    # Check if the dataset and file are accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)
    cTable = uvFile + '/' + checkTable
    if os.path.isfile(cTable):
        if options['verbose']:
            print('Flag version table exists.')
        return True
    else:
        if options['verbose']:
            print('Flag version table does not exist.')
        return False

def autoPgflag(uvFile, options):
    # Use the AOFlagger function in pgflag.
    if not options['quiet']:
        print('Automatically flagging dataset', uvFile)

    # Check if the file is accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)
    
    # We use Craig Anderson's flagging strategy here.
    stokes = [ 'i,q,u,v', 'v,q,u,i', 'v,q,u', 'u,v,q' ]
    flagpar = [ '8,5,5,3,6,3', '10,2,2,3,7,3',
                '8,2,2,3,6,3', '8,2,2,3,6,3' ]
    miriad.set_filter('pgflag', lineFilter)
    for n in range(0, len(stokes)):
        if options['verbose']:
            print(' Stage', (n+1))
        pr = miriad.pgflag(vis=uvFile, stokes=stokes[n], flagpar=flagpar[n],
                           options='nodisp', command='<b')

    if not options['quiet']:
        print('Flagging complete.')
    return True

def uvToMeasurementSet(uvFile, options):
    # Use the uv2ms Miriad task.
    msFile = uvFile + '.ms'
    if not options['quiet']:
        print('Converting Miriad dataset', uvFile,
              'into MeasurementSet', msFile)

    # Check if the file is accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)

    # See if the MeasurementSet already exists.
    if uvPresent(msFile, options):
        # Delete the tree.
        try:
            shutil.rmtree(msFile)
        except:
            print('Cannot delete MeasurementSet', msFile)
            sys.exit(0)

    miriad.set_filter('uv2ms', lineFilter)
    ur = miriad.uv2ms(vis=uvFile, ms=msFile)

    if not options['quiet']:
        print('Conversion complete.')
    return msFile

def summariseFile(progressDict, options):
    # Make a standard summary for a CABB data set, and try to
    # make some decisions about which reductions to do later on.
    if options['verbose']:
        print('Producing observation summary.')
    # We support multiple frequency configurations.
    for f in progressDict['freqConfigs']:
        fc = progressDict['freqConfigs'][f]
        fc['details'] = []
        for i in range(0, len(fc['channels'])):
            bw = fc['channels'][i] * fc['chanWidth'][i]
            if bw == 2.048:
                # This is a wideband IF.
                fc['details'].append({
                    'type': 'wideband',
                    'bandwidth': bw
                    })
            else:
                # This must be a zoom band IF.
                
                fc['details'].append({
                    'type': 'zoom',
                    'bandwidth': bw
                    })

    if options['verbose']:
        print('Observation summary complete.')
    return True
        
def midweekDetector(uvFile, options):
    # We try to detect if we have been affected by midweek RFI
    # and if so, when it started and stopped.
    if not options['quiet']:
        print('Checking for midweek RFI.')

    # We identify time ranges where we might have midweek RFI
    # using the uv variable 'xyamp'. This amplitude should be
    # fairly constant across the observation, but will become
    # very large and have a large spread when midweek RFI is
    # present.
    # Check that the dataset is accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)
        
    miriad.set_filter('varplt', lineFilter)
    # We output the data from varplt to a file.
    outLog = 'varplt.xyamp.out'
    if os.path.isfile(outLog):
        # Try to delete it.
        try:
            os.remove(outLog)
        except:
            print('Cannot delete existing varplt log.')
            sys.exit(0)
    rv = miriad.varplt(vis=uvFile, device='/null', xaxis='time',
                       yaxis='xyamp', log=outLog)
    # Check for the output file.
    if not os.path.isfile(outLog):
        print('Failed to make output for midweek RFI detection.')
        sys.exit(0)
    logData = []
    with open(outLog) as ol:
        for line in ol:
            els = line.rstrip('\n').split()
            if els[0] == '#':
                if els[1] == 'Base':
                    # This is the base time.
                    baseTime = datetime.strptime(els[4], '%y%b%d:%H:%M:%S.%f')
            else:
                startN = 0
                if len(els) == 6:
                    # Includes the offset time.
                    tels = els[1].split(':')
                    offTime = timedelta(days=int(els[0]), hours=int(tels[0]),
                                        minutes=int(tels[1]),
                                        seconds=int(tels[2]))
                    atime = baseTime + offTime
                    etime = calendar.timegm(atime.timetuple())
                    logData.append({ 'time': atime, 'epochTime': etime,
                                     'xyamp': [] })
                    startN = 2
                for i in range(startN, len(els)):
                    logData[-1]['xyamp'].append(float(els[i]))

    # Bin the data into 1 minute bins.
    binSizeSeconds = 60
    etimes = [logData[i]['epochTime'] for i, b in enumerate(logData)]
    timeBins = np.linspace(min(etimes), max(etimes),
                           int((max(etimes) - min(etimes)) / binSizeSeconds))
    digTimes = np.digitize(etimes, timeBins)
    # The quantity we want is the variance of the xy amplitudes.
    amps = np.array([np.ptp(logData[i]['xyamp']) for i, b in enumerate(logData)])
    # We classify as possible midweek RFI if the variance is greater than 200.
    tinc = np.array(
        [(amps[digTimes == i].var() > 200) for i in range(1, len(timeBins) + 1)])
    # And now extend the areas where midweek RFI may have been detected, such
    # that if we have detected midweek RFI within the previous 5 minutes,
    # and there is more midweek RFI within 5 minutes, we say it is still on.
    jinc = np.array(
        [True if (tinc[i - 5:i].any() and tinc[i:i + 5].any())
         else False for i, b in enumerate(tinc)])
    # Make a binary truth table.
    dinc = np.array([1 if b else 0 for i, b in enumerate(jinc)])
    # Map the transitions.
    transOn = np.array([1 if (dinc[i] == 0 and dinc[i + 1] == 1) else 0
                        for i in range(0, len(dinc) - 1)])
    transOff = np.array([1 if (dinc[i] == 1 and dinc[i + 1] == 0) else 0
                         for i in range(0, len(dinc) - 1)])
    # Get the time of each transition.
    onTimes = timeBins[np.where(transOn == 1)]
    offTimes = timeBins[np.where(transOff == 1)]
    # Extend to the edge of the previous/next time bin to be sure to catch it all.
    extendLength = 3 * binSizeSeconds / 2
    for i in range(0, len(onTimes)):
        onTimes[i] -= extendLength
    for i in range(0, len(offTimes)):
        offTimes[i] += extendLength
    # Account for special cases.
    if (len(onTimes) > 0):
        if min(offTimes) < min(onTimes):
            # Must have been on at the start.
            onTimes.insert(0, min(etimes))
        if max(transOn) > max(transOff):
            # Must still be on at the end.
            offTimes.append(max(etimes))
    # Make some start and stop dates compatible with Miriad's flagger.
    flagRegions = []
    for i in range(0, len(onTimes)):
        o = { 'start': time.strftime('%y%b%d:%H:%M:%S', time.gmtime(onTimes[i])).upper(),
              'stop': time.strftime('%y%b%d:%H:%M:%S', time.gmtime(offTimes[i])).upper() }
        flagRegions.append(o)

    # Our return dictionary.
    rDict = {
        'onTimes': onTimes, 'offTimes': offTimes, 'flagRegions': flagRegions
        }

    if not options['quiet']:
        print('Midweek RFI check complete.')
    return rDict

def midweekFlagger(uvFile, timeRegions, options):
    # Flag times when midweek RFI was detected.
    if not options['quiet']:
        print('Flagging times due to midweek RFI.')

    # Check that we actually have some regions to flag.
    if len(timeRegions) == 0:
        if not options['quiet']:
            print('No midweek RFI detected.')
            return False

    # Check we have access to the data.
    if not uvPresent(uvFile, options):
        sys.exit(0)

    # Run the uvflag task as required.
    miriad.set_filter('uvflag', lineFilter)
    for i in range(0, len(timeRegions)):
        rf = miriad.uvflag(vis=uvFile, flagval='flag',
                           select='time(' + timeRegions[i]['start'] + ',' +
                           timeRegions[i]['stop'] + ')')

    if not options['quiet']:
        print('Midweek RFI flagging complete.')
    return True

def determineCalibrators(progressDict, freqConfig, options):
    # Work out which calibrators are present and how to use them.
    if options['verbose']:
        print('Finding calibrators.')

    rDict = { 'flux': None,
              'bandpass': None,
              'leakages': None,
              'gain': [ ],
              'warnings': [] }
    # Get all the sources for this frequency config.
    sources = {}
    for s in progressDict['sources']:
        if freqConfig in progressDict['sources'][s]['freqConfig']:
            sources[s] = progressDict['sources'][s]

    # Determine the band we're in.
    band = frequencyBand(progressDict['freqConfigs'][freqConfig]['centreFreq'][0])

    # First, look for the obvious preferred flux calibrators.
    if band == '16cm' or band == '4cm' or band == '15mm':
        if '1934-638' in sources:
            rDict['flux'] = '1934-638'
        elif '0823-500' in sources:
            rDict['flux'] = '0823-500'
            rDict['warnings'].append(
                'Non-preferred flux calibration source 0823-500 used. ' +
                'Flux density scale is unlikely to be correct.')
        else:
            rDict['warnings'].append(
                'No flux calibration source found. Flux density ' +
                'calibration will not be attempted.')
    elif band == '7mm' or band == '3mm':
        fluxPriority = [ 'uranus', 'neptune', 'mars' ]
        for i in range(0, len(fluxPriority)):
            if fluxPriority[i] in sources:
                rDict['flux'] = fluxPriority[i]
                break
    elif ((band == '7mm') and ('1934-638' in sources)):
        rDict['flux'] = '1934-638'
        rDict['warnings'].append(
            'Non-preferred flux calibration source 1934-638 used. ' +
            'Flux density scale may not be correct.')
    if rDict['flux'] is None:
        rDict['warnings'].append(
            'No flux calibration source found. Flux density ' +
            'calibration will not be attempted.')

    # And bandpass calibrators.
    if band == '16cm' or band == '4cm':
        bandpassPriority = [ '1934-638', '0823-500', '1921-293',
                             '1253-055', '0537-441' ]
    elif band == '15mm':
        bandpassPriority = [ '1921-293', '1253-055', '0537-441',
                             '0420-014', '1934-638' ]
    else:
        bandpassPriority = [ '1921-293', '1253-055', '0537-441',
                             '0420-014' ]

    for i in range(0, len(bandpassPriority)):
        if bandpassPriority[i] in sources:
            rDict['bandpass'] = bandpassPriority[i]
            break
    if rDict['bandpass'] is None:
        rDict['warnings'].append(
            'No preferred bandpass calibration source has been ' +
            'found.')

    # Look for a possible leakage calibrator.
    if band == '16cm' or band == '4cm':
        if '1934-638' in sources:
            rDict['leakages'] = '1934-638'
    if rDict['leakages'] is None:
        # Look for the source with the most number of cuts.
        mn = 0
        ms = None
        for s in sources:
            nf = sources[s]['freqConfig'].count(freqConfig)
            if nf > mn:
                mn = nf
                ms = s
        if mn >= 3:
            # Three cuts will usually give us success with
            # leakages.
            rDict['leakages'] = ms
    if rDict['leakages'] is None:
        # Still no very good candidate for leakage calibration.
        rDict['warnings'].append(
            'No good calibrator was found for leakage calibration. ' +
            'Polarisation measurements may be inaccurate.')
        # And use 1934-638 as a last resort at all bands.
        if '1934-638' in sources:
            rDict['leakages'] = '1934-638'

    # Now we find the gain calibrator(s), beginning with the
    # sources marked with the calibrator code.
    for s in sources:
        if 'C' in sources[s]['calCode']:
            rDict['gain'].append(s)
        else:
            # Check whether this is a source in the ATCA calibrator
            # database.
            if s in caldb.calibrators:
                rDict['gain'].append(s)

    # We're finished.
    if options['verbose']:
        print('Calibrators found.')
    return rDict

def deleteCalTables(sets, options):
    # In the current directory, delete the calibration tables
    # from sets matching the glob-compatible 'sets' argument.
    tables = [ 'bandpass', 'gains', 'leakage', 'gainsf', 'leakagef' ]
    miriad.set_filter('delhd', lineFilter)
    success = True
    for t in tables:
        dfiles = glob.glob(sets + '/' + t)
        for d in dfiles:
            if options['verbose']:
                print('Deleting calibration table', d)
            dhr = miriad.delhd(In=d)
            # Check it did actually get deleted.
            if os.path.isfile(d):
                if options['verbose']:
                    print(' Delete failed!')
                success = False
    return success

def copyCalTables(origSet, destSet, copyOptions, options):
    # Copy the calibration tables from one dataset to another.
    if options['verbose']:
        print('Copying calibration tables from ' + origSet +
              ' to ' + destSet)

    # Check we have access to both the datasets.
    if not uvPresent(origSet, options):
        sys.exit(0)
    if not uvPresent(destSet, options):
        sys.exit(0)

    # Determine the options required.
    optionList = []
    if 'pol' in copyOptions and copyOptions['pol'] == False:
        optionList.append('nopol')
    if 'cal' in copyOptions and copyOptions['cal'] == False:
        optionList.append('nocal')
    if 'pass' in copyOptions and copyOptions['pass'] == False:
        optionList.append('nopass')
    optionString = ','.join(optionList)

    # Run gpcopy.
    miriad.set_filter('gpcopy', lineFilter)
    if len(optionList) > 0:
        gcr = miriad.gpcopy(vis=origSet, out=destSet, options=optionString)
    else:
        gcr = miriad.gpcopy(vis=origSet, out=destSet)

    if options['verbose']:
        print('Calibration tables copied.')
    return True

def correctBandpass(fluxCal, fluxSet, options):
    # Correct a bandpass solution that was obtained on a non-flux
    # calibration source.
    if options['verbose']:
        print('Correcting bandpass table using flux calibrator',
              fluxSet)

    # Check we have access to the dataset.
    if not uvPresent(fluxSet, options):
        sys.exit(0)

    # Run mfboot.
    miriad.set_filter('mfboot', lineFilter)
    specCorr = 1e6
    selString = 'source(' + fluxCal + ')'
    corrCount = 0
    while (math.fabs(specCorr) > 0.005 and corrCount < 5):
        mbr = miriad.mfboot(vis=fluxSet, select=selString, device='/null')
        corrCount += 1
        for l in mbr:
            lsp = re.split('\s+', l)
            if lsp[0] == 'Adjusting':
                specCorr = float(lsp[4])

    if options['verbose']:
        print('Bandpass correction complete.')
        
    if corrCount == 5:
        return False
    return True

def prepareReductionDir(uvFile, options):
    # Make a directory that can be used for the Miriad reduction
    # of the specified uv dataset.
    redDir = 'reduction.' + uvFile
    if options['verbose']:
        print('Preparing reduction directory', redDir)

    # Check that the file is accessible to us.
    if not uvPresent(uvFile, options):
        sys.exit(0)

    # Check if the reduction directory exists already.
    if os.path.isdir(redDir):
        # Are we allowed to delete it?
        if not options['keep_reduction']:
            # It does exist, let's delete it.
            try:
                shutil.rmtree(redDir)
            except:
                print('Cannot delete existing reduction directory',
                      redDir)
                sys.exit(0)

    # Make the reduction directory.
    if not options['keep_reduction']:
        try:
            os.mkdir(redDir)
        except:
            print('Cannot make reduction directory', redDir)
            sys.exit(0)

    # Change to the reduction directory and split out the
    # uv dataset.
    if not chReductionDir(redDir, options):
        sys.exit(0)

    if not options['keep_reduction']:
        # Split out the data.
        bFile = '../' + uvFile
        ur = miriad.uvsplit(vis=bFile)
    else:
        # Delete the calibration tables of any datasets here.
        deleteCalTables('*', options)

    # Go back to the lower directory.
    if not chReductionDir('..', options):
        sys.exit(0)

    if options['verbose']:
        print('Reduction directory prepared.')
    return redDir

def determineRefAnt(flagStats, options):
    # Try to work out the best antenna to use as calibration
    # reference. This is done by selecting the antenna with
    # the least amount of flagging, or CA03 by preference if
    # more than one has the same amount of flagging.
    mf = min(flagStats['antenna'].values())
    if (flagStats['antenna']['3'] == mf):
        if options['verbose']:
            print('Will use CA03 as reference antenna.')
        return '3'
    else:
        for a in flagStats['antenna']:
            if flagStats['antenna'][a] == mf:
                if options['verbose']:
                    print('Will use CA0' + a +
                          ' as reference antenna.')
                return a
    if options['verbose']:
        print('Reference antenna is undetermined. Defaulting ' +
              'to CA03.')
    return '3'

def calibrateBandpass(calSet, refAnt, options):
    # Perform a bandpass calibration.
    if options['verbose']:
        print('Performing bandpass calibration with dataset', calSet)

    # Check we have access to the dataset.
    if not uvPresent(calSet, options):
        sys.exit(0)

    # Run mfcal.
    miriad.set_filter('mfcal', lineFilter)
    mfr = miriad.mfcal(vis=calSet, refant=refAnt, interval=0.1)
    rDict = { 'converged': True,
              'iterations': 0,
              'fluxDensity': None }
    # Check for conditions.
    for l in mfr:
        lsp = re.split('\s+', l)
        if len(lsp) < 1:
            continue
        # Some possible matches.
        a = re.match('Iter=(\d+)\,', lsp[0])
        if lsp[0] == 'I' and lsp[1] == 'flux':
            rDict['fluxDensity'] = float(lsp[3])
        elif lsp[0] == '###' and lsp[2] == 'Failed' and lsp[4] == 'converge':
            rDict['converged'] = False
        elif a is not None:
            rDict['iterations'] = int(a.group(1))
            
    if options['verbose']:
        print('Bandpass calibration complete.')
    return rDict

def calibrateGains(calSet, refAnt, calOptions, options):
    # Perform gain calibration.
    if options['verbose']:
        print('Performing gain calibration with dataset', calSet)

    # Check we have access to the dataset.
    if not uvPresent(calSet, options):
        sys.exit(0)

    # Set up the options.
    optList = [ 'xyvary' ]
    if 'leakages' in calOptions:
        optList.append('qusolve')
        if calOptions['leakages'] == False:
            optList.append('nopol')

    optStr = ','.join(optList)

    # Run gpcal.
    miriad.set_filter('gpcal', lineFilter)
    gcr = miriad.gpcal(vis=calSet, refant=refAnt, interval=0.1,
                       options=optStr)
    rDict = { 'converged': True,
              'iterations': 0,
              'fluxDensity': None,
              'leakageSolved': False }
    for l in gcr:
        lsp = re.split('\s+', l)
        if len(lsp) < 1:
            continue
        # Some possible matches.
        a = re.match('Iter=\s*(\d+),', lsp[0])
        if lsp[0] == 'I' and lsp[1] == 'flux':
            rDict['fluxDensity'] = float(lsp[3])
        elif lsp[0] == 'Leakage' and lsp[1] == 'terms:':
            rDict['leakageSolved'] = True
        elif a is not None:
            rDict['iterations'] = int(a.group(1))

    if options['verbose']:
        print('Gain calibration complete.')
    return rDict
    
    
