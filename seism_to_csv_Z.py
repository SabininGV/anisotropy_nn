# Z-component
# column -- timesample, row -- trace
import obspy
import numpy as np
import pandas as pd

st = obspy.read('model_1\\model+GathNP-Z.sgy')
N_traces = len(st)
N_timesamples = len(st[0])


for i in range(1,1101):
    print('\t',i,'\t/\t1100')
    st = obspy.read('model_{}/model+GathNP-Z.sgy'.format(i))
    data = np.empty([N_timesamples])
    for n in range(N_traces):
        data = np.vstack((data,st[n]))
    df = pd.DataFrame(data[1:])
    df.to_csv('csv_models_2frac_Thomsen_saturated_full_formulae_Z\\model_{}.csv'.format(i),index=None)
