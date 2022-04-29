#!/usr/bin/env python

import sys
import numpy as np
import math
import os
from subprocess import call
from shutil import copyfile
from datetime import datetime

sys.path.insert(0, '../')
import pc

dir_name = "/home/nblinov/HEP/MG5_aMC_v2_6_4/"+sys.argv[1]+'/'
mg_path = dir_name+"bin/generate_events"
f_param_card_name = dir_name+"Cards/param_card.dat"
f_run_card_name = dir_name+"Cards/run_card.dat"

if not os.path.exists(dir_name):
    print "Run directory not found!"
    sys.exit(0)

# Back up the original param_card
print "Backed up original param_card and run_card..."
copyfile(f_param_card_name, f_param_card_name + '.backup')
copyfile(f_run_card_name, f_run_card_name + '.backup')

param_card = pc.ParamCard(f_param_card_name)

mass_block = "mass"
Ap_code = 622
Nuc_code = 623
Electron_code = 11
Muon_code = 13

model_param_block = "frblock"
Anuc_code = 6
Znuc_code = 7
inelastic1_code = 8
inelastic2_code = 9
aval_code = 10
apval_code = 11
dval_code = 12

decay_block = "decay"

mAp_list = np.logspace(-2.,np.log10(2.),10) 
mAp_list = np.logspace(-3.,np.log10(0.03),10) 
#mAp_list = np.logspace(-2.,-1.,10) 
#mAp_list = np.array([mAp_list[0]])


# Run settings
scale_factor = 10. # multiplies all relevant energy scales
N_event = 10000
ebeam1 = 0.05
#ebeam1 = 4.0
#ebeam1 = 8.0
ebeam2 = 171.246 

# default values
me = 0.51099895/1000.
mmu = 105.6583745/1000.
mnuc = ebeam2

Anuc = 184.0
Znuc = 74.0
inelastic1 = 1.9276
inelastic2 = 1.40845
aval = 111.0/(0.0005111*np.power(Znuc,1./3.))
apval = 773.0/(0.0005111*np.power(Znuc,2./3.))
dval = 0.164/np.power(Anuc,2./3.)


# Modify run card with new number of events and possibly rescale the bea energies
with open(f_run_card_name + '.backup') as f_in, open(f_run_card_name,'w') as f_out:
    for line in f_in:
        new_line = line
        # Adjust the number of events to generate
        if 'nevents' in line:
            s = str.split(line,' = ')
            new_line = '  ' + str(N_event) + ' = ' + ''.join(s[1:]) 
        # Adjust the beam energies
        if 'ebeam' in line:
            s = str.split(line,' = ')

            if 'ebeam1' in line:
                ebeam = scale_factor*ebeam1 #float(s[0].strip())
            if 'ebeam2' in line:
                ebeam = scale_factor*ebeam2 #float(s[0].strip())

            new_line = '     ' + str(ebeam) + ' = ' + ''.join(s[1:])

        f_out.write(new_line)

f_readme = open(dir_name+"Events/README",'w', 0)
f_readme.write('This folder contains events generated with the following parameters:\n')
f_readme.write('\n\n')
f_readme.write('Runs:\n')
run_name_prev = ""

def rescale_param(block, param_code, r):
    param_val = r*param_card[block].get(param_code).value
    param_card.mod_param(block,param_code,value=param_val)

for mAp in mAp_list:
    start_time = datetime.now()
    run_name = "mAp_"+str(int(math.floor(mAp*1000.)))# + "_lh"
    if run_name == run_name_prev:
        run_name = run_name + "b"

    run_name_prev = run_name

    if param_card.has_param("mass", Ap_code):

        # Rescale all relevant dimensionful scales
        # Update masses 
        param_card.mod_param(mass_block,Nuc_code,value=scale_factor*mnuc)
        param_card.mod_param(mass_block,Electron_code,value=scale_factor*me)
        param_card.mod_param(mass_block,Muon_code,value=scale_factor*mmu)

        # Update form-factor numbers
        # the rescaling factors come from the dimension of the FF parameters
        param_card.mod_param(model_param_block, inelastic1_code, value=inelastic1/(scale_factor**2))
        param_card.mod_param(model_param_block, inelastic2_code, value=inelastic2/(scale_factor**2))
        param_card.mod_param(model_param_block, aval_code, value=aval/scale_factor)
        param_card.mod_param(model_param_block, apval_code, value=apval/scale_factor)
        param_card.mod_param(model_param_block, dval_code, value=dval*scale_factor**2)

        param_card.mod_param(mass_block,Ap_code,value=scale_factor*mAp)
        #param_card.mod_param(decay_block,Ap_code,value='Auto')
        param_card.mod_param(decay_block,Ap_code,value='1e-6')

        f_readme.write('mAp = ' + str(mAp) + '; start: ' + str(start_time))
        print "Changed mAp to ", param_card[mass_block].get(Ap_code).value

        if scale_factor > 1:
            print "Changed mNuc to ", param_card[mass_block].get(Nuc_code).value
            print "Changed me to ", param_card[mass_block].get(Electron_code).value
            print "Changed mmu to ", param_card[mass_block].get(Muon_code).value

        print "Writing new param_card.dat"
        param_card.write(f_param_card_name)

        print "OK!"

    print "Running MG5..."
    call([mg_path,"-f", run_name])
    print "Done!"
    end_time = datetime.now()
    elapsed_time = end_time - start_time
    f_readme.write('; finish: ' + str(end_time) + '; elapsed: ' + str(elapsed_time) + '\n')

f_readme.close()

