class PangenomeHeritabilityError(Exception):
    """Base exception for the package"""
    pass

class InputError(PangenomeHeritabilityError):
    """Raised when input files are invalid"""
    pass

class AlignmentError(PangenomeHeritabilityError):
    """Raised when MAFFT alignment fails"""
    pass
