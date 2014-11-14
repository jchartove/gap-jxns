#!/bin/bash

#$ -j y

matlab -nodisplay -r "cd /projectnb/crc-nak/chartove, matlabpool open, inhib_gj(5000,2,1,5,1), matlabpool close force local, quit"