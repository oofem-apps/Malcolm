import math
from enum import Enum


class Task_status(Enum):
    UNSELECTED = 1
    SELECTED = 2
    PROGRESS = 3
    COMPLETED = 4
    PREDEFINED = 5

class task:
    def __init__(self, status, eccentricity_normalized):

        self.status = status

        self.file_path = []
        
        self.eccentricity_normalized =  eccentricity_normalized
        self.eccentricity_actual = []
        self.engng_definition = []

        self.max_load = math.nan
        self.max_MN = []
        
        self.load_level = []
        self.N = []
        self.M = []
        self.eps = []
        self.kappa = []

    def reset(self):

        self.max_load = math.nan
        self.max_MN = []
        
        self.load_level = []
        self.N = []
        self.M = []
        self.eps = []
        self.kappa = []

        if(self.status == Task_status.COMPLETED):
            self.status = Task_status.SELECTED


        
