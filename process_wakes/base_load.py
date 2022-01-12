

class LoadWakes:
    ANALYSIS_TYPES = {
        'dx',  # horizontal impedance
        'dy',  # vertical impedance
        'db',  # both planes are symmetric
        'll'   # longitudinal and transverse quadrupolar impedances
        }

    def __init__(self, plane, simul_data):
        self.plane = plane
        self.simul_data = simul_data

    def load_and_process_raw_data(self):
        return NotImplementedError('Not Implemented.')
