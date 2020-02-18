from shutil import copyfile
import os
from random import random as rand
import numpy as np
import pandas as pd
from subprocess import call
import time
import datetime
from math import sqrt

Vs = 2750
Vp = 5000
Ro = 2550
    
def Thomsen_to_stiffness(eps1,eps2, delta1,delta2,delta3, gamma1,gamma2):
    # так представляется в Tesseral Pro
    Vp_ = Vp/1000
    Vs_ = Vs/1000
    Ro_ = Ro/1000
    
    c33 = Ro_ * (Vp_**2)
    c55 = Ro_ * (Vs_**2)
    c11 = 2*c33*eps2 + c33
    c66 = 2*c55*gamma1 + c55
    c22 = 2*c33*eps1 + c33
    c44 = c66/(2*gamma2 + 1)
    c13 = sqrt(2*c33*(c33-c55)*delta2 + (c33-c55)**2) - c55
    c23 = sqrt(2*c33*(c33-c44)*delta1 + (c33-c44)**2) - c44
    #c12 = (c23*c11 - c13*c22)/(c13 - c23)
    c12 = sqrt(2*c11*(c11-c66)*delta3+(c11-c66)**2) - c66
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
    Vs = 2750
    Vp = 5000
    Ro = 2550
    g = (Vs**2)/(Vp**2)
    e1 = 0.1*rand()
    e2 = 0.1*rand()

    # saturated
    dn1 = 0
    dt1 = 16*e1/(3*(3-2*g))
    dn2 = 0
    dt2 = 16*e2/(3*(3-2*g))
    eps1 = -2*g*(1-g)*dn2/(1-dn2*(1-2*g)**2)
    delta1 = -2*g*((1-2*g)*dn2 + dt2)*(1-(1-2*g)*dn2)/( (1-dn2*(1-2*g)**2)*(1+(g*dt2 - dn2*(1-2*g)**2)/(1-g)) )
    gamma1 = -1*dt2/2
    eps2 = -2*g*(1-g)*dn1/(1-dn1*(1-2*g)**2)
    delta2 = -2*g*((1-2*g)*dn1 + dt1)*(1-(1-2*g)*dn1)/( (1-dn1*(1-2*g)**2)*(1+(g*dt1 - dn1*(1-2*g)**2)/(1-g)) )
    gamma2 = -1*dt1/2
    delta3 = delta1 + delta2 - 2*eps2

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


