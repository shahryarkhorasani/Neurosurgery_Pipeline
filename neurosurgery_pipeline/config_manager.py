from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import datetime
import time


class Flags:
    def __init__(self):
        self.template0 = '/mnt/projects/charite-brains/derived/MICCAI2012_Multi_Atlas_Challenge_Data/T_template0.nii.gz'
        self.template_cerebellum = '/mnt/projects/charite-brains/derived/MICCAI2012_Multi_Atlas_Challenge_Data/T_template0_BrainCerebellumProbabilityMask.nii.gz'
        self.template_cerebellum_mask = '/mnt/projects/charite-brains/derived/MICCAI2012_Multi_Atlas_Challenge_Data/T_template0_BrainCerebellumRegistrationMask.nii.gz'

FLAGS = Flags()