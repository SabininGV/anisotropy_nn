from shutil import copyfile
import os
from random import random as rand
import numpy as np
import pandas as pd
from subprocess import call
import time
import datetime

Vs = 2750
Vp = 5000
Ro = 2550
g = (Vs**2)/(Vp**2)

def dndt_to_stiffness(dn1,dt1, dn2,dt2):
    # так представляется в Tesseral Pro
    Vp_ = Vp/1000
    Vs_ = Vs/1000
    Ro_ = Ro/1000
    mu = Ro_*(Vs_**2)
    lambd = Ro_*(Vp_**2)-2*mu
    r = 1 - 2*g
    d = 1 - (r**2)*dn1*dn2
    L1 = 1 - dn1
    L2 = 1 - r*dn1
    L3 = 1 - (r**2)*dn1
    L4 = 4*(r**2)*(g**2)*dn1*dn2
    M1 = 1 - dn2
    M2 = 1 - r*dn2
    M3 = 1 - (r**2)*dn2
    
    c11 = ((lambd+2*mu)*L1*M3)/d
    c12 = lambd*L1*M1/d
    c13 = lambd*L1*M2/d
    c22 = ((lambd+2*mu)*L3*M1)/d
    c23 = lambd*L2*M1/d
    c33 = (lambd+2*mu)*(L3*M3-L4)/d
    c44 = mu*(1-dt2)
    c55 = mu*(1-dt1)
    c66 = mu*(1-dt1)*(1-dt2)/(1-dt1*dt2)
    
    return "\"0\" {0:.6f} {1:.6f} {2:.6f} 0 0 0 {3:.6f} {4:.6f} 0 0 0 {5:.6f} 0 0 0 {6:.6f} 0 0 {7:.6f} 0 {8:.6f}".format(c11,c12,c13,c22,c23,c33,c44,c55,c66)

def copy_files(dir):
    # tgr, tgr.cr0, tgr.cr2 files
    for word in ['Gath']:
        for form in ['tgr','tgr.cr0','tgr.cr2']:
            str = 'model+{0}NP.{1}'.format(word,form)
            file = os.path.join(dir,str)
            temp_file = os.path.join(cwd,str)
            copyfile(temp_file, file)
            os.remove(temp_file)

    # SEGY files
    for ch in ['C','X','Z']:
        str = 'model+GathNP-{}.sgy'.format(ch)
        file = os.path.join(dir,str)
        temp_file = os.path.join(cwd,str)
        copyfile(temp_file, file)
        os.remove(temp_file)

    # the rest of the files
    #file = os.path.join(dir,'model+GathNP-Z.sgy.cr0')
    #copyfile(os.path.join(cwd,'model+GathNP-Z.sgy.cr0'), file)
    file = os.path.join(dir,'model+WaveNP-1.tgr')
    copyfile(os.path.join(cwd,'model+WaveNP-1.tgr'), file)
    os.remove(os.path.join(cwd,'model+WaveNP-1.tgr'))

print('Enter the number of generated models:', end=' ')
N = int(input())
print('Enter the first model number:',end=' ')
K = int(input())

start_time = datetime.datetime.now()

cwd = os.getcwd()
pattern = os.path.join(cwd,'pattern\\pattern_orthorhombic.tam')
runtask_original = os.path.join(cwd,'runtask_orthorhombic.ini')
runtask = os.path.join(cwd,'runtask.ini')
copyfile(runtask_original, runtask)
data = np.empty([9]) # массив, в к-м будут храниться данные по всем моделям

for i in range(K,N+1):
    dir = os.path.join(cwd,"model_"+str(i))
    if not os.path.exists(dir):
        os.mkdir(dir)
    #model = os.path.join(dir,'model1.tam')
    model = os.path.join(cwd,'model1.tam')
    form = os.path.join(cwd,'form2.txt')
    copyfile(pattern, model)

    # dn, dt generation from Bakulin Grechka Tsvankin article "Estimation of fracture parameters from reflection seismic data -- Part I"
    e1 = 0.1*rand()
    e2 = 0.1*rand()

    # saturated
    dn1 = 0
    dt1 = 16*e1/(3*(3-2*g))
    dn2 = 0
    dt2 = 16*e2/(3*(3-2*g))

    # записываем элементы тензора в файл
    f = open(form,'tw')
    f.write(Thomsen_to_stiffness(eps1,eps2, delta1,delta2,delta3, gamma1,gamma2))
    f.close()
    
    copyfile(form,os.path.join(dir,'form2.txt'))
    copyfile(model, os.path.join(dir,'model1.tam'))

    frac_params = np.array([eps1,eps2, delta1,delta2,delta3, gamma1,gamma2, e1, e2])
    data = np.vstack((data,frac_params)) # записываем осреднённые параметры в массив

    # запускаем вычислительный модуль
    log_file = open(os.path.join(dir,'log.txt'),'w')
    call([os.path.join(cwd,"Tesseral2D_Win64.exe")], universal_newlines=True, stdout=log_file)
    log_file.close()

    # копируем выходные файлы в папку соответствующей модели
    copy_files(dir)
    # удаляем model1.tam и form2.txt (на всякий пожарный)
    os.remove(os.path.join(cwd,'model1.tam'))
    os.remove(os.path.join(cwd,'form2.txt'))


# записываем параметры в csv-файл
df = pd.DataFrame(data[1:])
df.to_csv(str(cwd) + '\\frac_params_2frac_Thomsen_saturated.csv',index=None, header=['eps1','eps2', 'delta1','delta2','delta3', 'gamma1','gamma2','e1','e2'])

# прошедшее время
sec = (datetime.datetime.now() - start_time).total_seconds()
print(sec,'sec')
print('Time elapsed:',int(sec/60),'min',sec%60,'sec')


