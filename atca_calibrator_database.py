import json
import os

# A module for stuff that we need the ATCA calibrator database for.

# Read in the list of calibrators.
def readCalibratorList():
    # Locate the file we can use.
    checkPaths = [ '/n/ste616/atcacode/cabb_pipeline',
                   '/home/jstevens/usr/src/atcacode/cabb_pipeline' ]
    a = None
    for i in range(0, len(checkPaths)):
        cc = checkPaths[i] + '/atca-caldb.json'
        if os.path.isfile(cc):
            with open(cc) as f:
                a = f.read()
            break

    if a is None:
        return None

    b = json.loads(a)
    calibrators = b['caldb']
    return calibrators

# Do this by default.
calibrators = readCalibratorList()


