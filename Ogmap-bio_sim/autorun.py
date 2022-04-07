# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 22:56:49 2021

@author: jpw
"""


#run -i gte_utils.py

#run -i ogmap.py

#run -i bottom_survey.py

#make_arena_directory('O:\OneDrive\Masters 2020-2022\Main Project\Ogmap\local.arena', 'SFA6 test', SFA=6)

#preprocess_data("O:\OneDrive\Masters 2020-2022\Main Project\Ogmap\PB_fall.dat", 35)

#use_ogmap()

from random import randint

from gte_utils import *

from ogmap import *

from bottom_survey import *

make_arena_directory('C:\Users\jpw\OneDrive\Masters 2020-2022\Main Project\Ogmap2 - Demo\local.arena', 'SFA6 test', SFA=6)

preprocess_data("C:\Users\jpw\OneDrive\Masters 2020-2022\Main Project\Ogmap2 - Demo\PB_fall.dat", 2)

use_ogmap()