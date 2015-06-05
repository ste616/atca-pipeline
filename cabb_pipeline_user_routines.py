from __future__ import print_function

"""
This is a set of routines that any user-made routines should import
to assist them in making a compliant part of the CABB data reduction
pipeline.
"""

def USER_print_name(name, version, options):
    # USER routines must always print out their name and version
    # when they are entered and exited.
    if not options['quiet']:
        print('USER', name, 'v' + version)
    return True

def USER_routine_successful(name, version, options):
    USER_print_name(name, version, options)
    # USER routines must return True if they exit successfully.
    return True

def USER_routine_unsuccessful(name, version, options):
    USER_print_name(name, version, options)
    # User routines must return False if they exit unsuccessfully.
    return False
