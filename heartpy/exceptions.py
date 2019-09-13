'''
Custom exceptions and warnings for HeartPy
'''

__all__ = ['BadSignalWarning']

class BadSignalWarning(UserWarning):
    '''
    warning class to raise when no heart rate is detectable
    in supplied signal.
    
    This warning notifies the user that the supplied signal is
    of insufficient quality and/or does not contain enough information
    to properly process.
    '''

class IncorrectFileType(UserWarning):
    '''
    warning class to raise when incorrect file type or incorrectly
    formatted file is provided to data loader
    '''
