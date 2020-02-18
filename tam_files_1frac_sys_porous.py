from shutil import copyfile
import os
from random import random as rand
import numpy as np
import pandas as pd
from subprocess import call
import time
import datetime
from math import pi

def anisotropic_parameters(dn,dt,alpha=5.0):
    return "lh= '0 (base)' ls= 65536 lac= 0.000000 las= 0.000000 laa= 0.000000 lab=0.000000 tf=0.000000 tg=0.000000 tf1={0:.6f} tf2=0.000000 tf3=0.000000 tdn1={1:.6f} tdn2=0.000000 tdn3=0.000000 tdt1={2:.6f} tdt2=0.000000 tdt3=0.000000 tfb=1\n".format(alpha,dn,dt)

def copy_files(dir):
    # tgr, tgr.cr0, tgr.cr2 files
    for word in ['Gath','Snap']:
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
pattern = os.path.join(cwd,'pattern\\pattern.tam')
runtask_original = os.path.join(cwd,'runtask_original.ini')
runtask = os.path.join(cwd,'runtask.ini')
copyfile(runtask_original, runtask)
data = np.empty([5]) # массив, в к-м будут храниться данные по всем моделям

######################################
# добавить чтение data из файла .csv #
######################################

for i in range(K,N+1):
    dir = os.path.join(cwd,"model_"+str(i))
    if not os.path.exists(dir):
        os.mkdir(dir)
    #model = os.path.join(dir,'model1.tam')
    model = os.path.join(cwd,'model1.tam')
    copyfile(pattern, model)

    # записываем случайные параметры в файл
    f = open(model, 'r')    # pass an appropriate path of the required file
    lines = f.readlines()
    n = 46
    # dn, dt generation from Bakulin Grechka Tsvankin article "Estimation of fracture parameters from reflection seismic data -- Part I"
    Vs = 2750
    Vp = 5000
    Ro = 2550
    g = (Vs**2)/(Vp**2)
    mu = Ro*(Vs**2)
    lambd = Ro*(Vp**2)-2*mu
    
    e = 0.1*rand()

    # Petroleum!
    k_1 = 1.25*(10**9)
    aspect_ratio = 1/10000 + rand()*(1/1000 - 1/10000) # mesofractures -- 1/10000 through 1/1000
    phi_c = 4*pi*e/3 * aspect_ratio # from Bakulin Grechka Tsvankin
    phi_p = 1/1000 # 0.1%
    A_p = 3/(4*g) # "Numerical Modeling of P-wave AVOA in Media Containing Vertical Fracures" Chichinina, Sabinin, Ronquillo Jarillo
    A_c = 4/9 * (3-4*g)/(g*(1-g))

    Dcp = 1/( 1 - k_1/(lambd + 2*mu/3) + k_1/((lambd + 2*mu/3)*(phi_c + phi_p))*(A_p*phi_p + A_c*e) )
    dn = 4*e/(3*g*(1-g)) * (1 - k_1/(lambd + 2*mu/3)) * Dcp
    
    dt = 16*e/(3*(3-2*g))
    alpha = 5.0
    lines[n-1] = anisotropic_parameters(dn,dt,alpha)    # n is the line number you want to edit; subtract 1 as indexing of list starts from 0
    f.close()   # close the file and reopen in write mode to enable writing to file; you can also open in append mode and use "seek", but you will have some unwanted old data if the new data is shorter in length.
    f = open(model, 'w')
    f.writelines(lines)
    # do the remaining operations on the file
    f.close()

    copyfile(model, os.path.join(dir,'model1.tam'))

    frac_params = np.array([dn,dt,e,aspect_ratio,alpha])
    data = np.vstack((data,frac_params)) # записываем осреднённые параметры в массив

    """
    # указываем путь к tam-файлу в runtask.ini:
    f = open(runtask, 'r')    # pass an appropriate path of the required file
    lines_runtask = f.readlines()
    #for line in lines_runtask:
    #    print(line)
    n = 4
    #"Model Name=model_" + str(i) + "\\model.tam\n"
    lines_runtask[n-1] = "Model Name=model.tam\n"   # n is the line number you want to edit; subtract 1 as indexing of list starts from 0
    f.close()   # close the file and reopen in write mode to enable writing to file; you can also open in append mode and use "seek", but you will have some unwanted old data if the new data is shorter in length.
    #for line in lines_runtask:
    #    print(line)
    f = open(runtask, 'w')
    f.writelines(lines_runtask)
    # do the remaining operations on the file
    f.close()
    """

    # запускаем вычислительный модуль
    log_file = open(os.path.join(dir,'log.txt'),'w')
    call([os.path.join(cwd,"Tesseral2D_Win64.exe")], universal_newlines=True, stdout=log_file)
    log_file.close()

    # копируем выходные файлы в папку соответствующей модели
    copy_files(dir)
    # удаляем model1.tam (на всякий пожарный)
    os.remove(os.path.join(cwd,'model1.tam'))


# записываем параметры в csv-файл
df = pd.DataFrame(data[1:])
df.to_csv(str(cwd) + '\\frac_params_porous.csv',index=None, header=['dn','dt','e','aspect_ratio','alpha'])

# прошедшее время
sec = (datetime.datetime.now() - start_time).total_seconds()
print(sec,'sec')
print('Time elapsed:',int(sec/60),'min',sec%60,'sec')


