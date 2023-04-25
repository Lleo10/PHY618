# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 14:07:13 2023

@author: Nishant10
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


tlrf = r"D:\Study Material\College\MS Thesis\My server\nd50 10m\timesteps.csv"
tlr = open(tlrf)
time_loopavg_arr = np.loadtxt(tlr, delimiter=",")

e25 = r"D:\Study Material\College\MS Thesis\My server\nd25 10m\e2e_nd0_ts10.0M.csv"
e50 = r"D:\Study Material\College\MS Thesis\My server\nd50 10m\e2e_nd0_ts10.0M.csv"
e75 = r"D:\Study Material\College\MS Thesis\My server\nd75 10m\e2e_nd0_ts10.0M.csv"
e100 = r"D:\Study Material\College\MS Thesis\My server\nd100 10m\e2e_nd0_ts10.0M.csv"
e150 = r"D:\Study Material\College\MS Thesis\My server\nd150 10m\e2e_nd0_ts10.0M.csv"
e200 = r"D:\Study Material\College\MS Thesis\My server\nd200 10m\e2e_nd0_ts10.0M.csv"
e225 = r"D:\Study Material\College\MS Thesis\My server\nd225 10m\e2e_nd0_ts10.0M.csv"
e245 = r"D:\Study Material\College\MS Thesis\My server\nd245 10m\e2e_nd0_ts10.0M.csv"
filex1 = open(e25)
filex = open(e50)
filey = open(e75)
filez = open(e150)
filez1 = open(e100)
filez2 = open(e200)
filez3 = open(e225)
filez4 = open(e245)

x1 = np.loadtxt(filex1, delimiter=",")
x1 = np.delete(x1,0)
x = np.loadtxt(filex, delimiter=",")
y = np.loadtxt(filey, delimiter=",")
z = np.loadtxt(filez, delimiter=",")
z1 = np.loadtxt(filez1, delimiter=",")
z2 = np.loadtxt(filez2, delimiter=",")
z2 = np.delete(z2,9999)
z3 = np.loadtxt(filez3, delimiter=",")
z3 = np.delete(z3,9999)
z4 = np.loadtxt(filez4, delimiter=",")
z4 = np.delete(z4,9999)




crop_index = int(len(x) * 0.65)






plt.plot(time_loopavg_arr,x1,label='25')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,x,label='50')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,y,label='75')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,z1,label='100')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,z,label='150')

plt.plot(time_loopavg_arr,z2,label='200')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,z3,label='225')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,z4,label='245')
plt.legend(title = 'No. of Dimers')


plt.xlabel('Timesteps')
plt.ylabel('End to End distance')
plt.legend(title = 'No. of Dimers')
plt.show()



r25 = r"D:\Study Material\College\MS Thesis\My server\nd25 10m\rog_nd0_ts10.0M.csv"
r50 = r"D:\Study Material\College\MS Thesis\My server\nd50 10m\rog_nd0_ts10.0M.csv"
r75 = r"D:\Study Material\College\MS Thesis\My server\nd75 10m\rog_nd0_ts10.0M.csv"
r100 = r"D:\Study Material\College\MS Thesis\My server\nd100 10m\rog_nd0_ts10.0M.csv"
r150 = r"D:\Study Material\College\MS Thesis\My server\nd150 10m\rog_nd0_ts10.0M.csv"
r200 = r"D:\Study Material\College\MS Thesis\My server\nd200 10m\rog_nd0_ts10.0M.csv"
r225 = r"D:\Study Material\College\MS Thesis\My server\nd225 10m\rog_nd0_ts10.0M.csv"
r245 = r"D:\Study Material\College\MS Thesis\My server\nd245 10m\rog_nd0_ts10.0M.csv"
filexr1 = open(r25)
filexr = open(r50)
fileyr = open(r75)
filezr = open(r150)
filezr1 = open(r100)
filezr2 = open(r200)
filezr3 = open(r225)
filezr4 = open(r245)
  
xr1 = np.loadtxt(filexr1, delimiter=",")
xr1 = np.delete(xr1,0)
xr = np.loadtxt(filexr, delimiter=",")
yr = np.loadtxt(fileyr, delimiter=",")
zr = np.loadtxt(filezr, delimiter=",")
zr1 = np.loadtxt(filezr1, delimiter=",")
zr2 = np.loadtxt(filezr2, delimiter=",")
zr2 = np.delete(zr2,9999)
zr3 = np.loadtxt(filezr3, delimiter=",")
zr3 = np.delete(zr3,9999)
zr4 = np.loadtxt(filezr4, delimiter=",")
zr4 = np.delete(zr4,9999)

'''plt.plot(time_loopavg_arr,xr1,label='25')
plt.legend(title='No. of Dimers')

plt.plot(time_loopavg_arr,xr,label='50')
plt.legend(title='No. of Dimers')

plt.plot(time_loopavg_arr,yr,label='75')
plt.legend(title='No. of Dimers')

plt.plot(time_loopavg_arr,zr1,label='100')
plt.legend(title='No. of Dimers')

plt.plot(time_loopavg_arr,zr2,label='200')
plt.legend(title='No. of Dimers')

plt.plot(time_loopavg_arr,zr,label='150')
plt.xlabel('Timesteps')
plt.ylabel('ROG')
plt.legend(title='No. of Dimers')
plt.show()



cropped_25 = xr1[crop_index:]
cropped_50 = xr[crop_index:]
cropped_75 = yr[crop_index:]
cropped_100 = zr1[crop_index:]
cropped_150 = zr[crop_index:]
cropped_200 = zr2[crop_index:]

time_cropped = time_loopavg_arr[crop_index:]


plt.plot(time_cropped,cropped_25,label='25')
plt.legend(title='No. of Dimers')

plt.plot(time_cropped,cropped_50,label='50')
plt.legend(title='No. of Dimers')

plt.plot(time_cropped,cropped_75,label='75')
plt.legend(title='No. of Dimers')

plt.plot(time_cropped,cropped_100,label='100')
plt.legend(title='No. of Dimers')

plt.plot(time_cropped,cropped_150,label='150')
plt.legend(title='No. of Dimers')

plt.plot(time_cropped,cropped_200,label='200')
plt.xlabel('Timesteps')
plt.ylabel('ROG zoomed in')
plt.legend(title='No. of Dimers')
plt.show()'''






lp25 = r"D:\Study Material\College\MS Thesis\My server\nd25 10m\loop_nd0_ts10.0M.csv"
lp50 = r"D:\Study Material\College\MS Thesis\My server\nd50 10m\loop_nd0_ts10.0M.csv"
lp75 = r"D:\Study Material\College\MS Thesis\My server\nd75 10m\loop_nd0_ts10.0M.csv"
lp100 = r"D:\Study Material\College\MS Thesis\My server\nd100 10m\loop_nd0_ts10.0M.csv"
lp150 = r"D:\Study Material\College\MS Thesis\My server\nd150 10m\loop_nd0_ts10.0M.csv"
lp200 = r"D:\Study Material\College\MS Thesis\My server\nd200 10m\loop_nd0_ts10.0M.csv"
lp225 = r"D:\Study Material\College\MS Thesis\My server\nd225 10m\loop_nd0_ts10.0M.csv"
lp245 = r"D:\Study Material\College\MS Thesis\My server\nd245 10m\loop_nd0_ts10.0M.csv"
filexl1 = open(lp25)
filexl = open(lp50)
fileyl = open(lp75)
filezl = open(lp150)
filezl1 = open(lp100)
filezl2 = open(lp200)
filezl3 = open(lp225)
filezl4 = open(lp245)

xl1 = np.loadtxt(filexl1, delimiter=",")  
xl1 = np.delete(xl1,0)
xl = np.loadtxt(filexl, delimiter=",")
yl = np.loadtxt(fileyl, delimiter=",")
zl = np.loadtxt(filezl, delimiter=",")
zl1 = np.loadtxt(filezl1, delimiter=",")
zl2 = np.loadtxt(filezl2, delimiter=",")
zl2 = np.delete(zl2,9999)
zl3 = np.loadtxt(filezl3, delimiter=",")
zl3 = np.delete(zl3,9999)
zl4 = np.loadtxt(filezl4, delimiter=",")
zl4 = np.delete(zl4,9999)

plt.plot(time_loopavg_arr,xl1,label='25')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,xl,label='50')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,yl,label='75')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,zl1,label='100')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,zl,label='150')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,zl2,label='200')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,zl3,label='225')
plt.legend(title = 'No. of Dimers')

plt.plot(time_loopavg_arr,zl4,label='245')
plt.legend(title = 'No. of Dimers')


plt.xlabel('Timesteps')
plt.ylabel('Average Loop Length')
plt.legend(title = 'No. of Dimers')
plt.show()







cropped_25 = xr1[crop_index:]
cropped_50 = xr[crop_index:]
cropped_75 = yr[crop_index:]
cropped_100 = zr1[crop_index:]
cropped_150 = zr[crop_index:]
cropped_200 = zr2[crop_index:]
cropped_225 = zr3[crop_index:]
cropped_245 = zr4[crop_index:]


avg_cropped_25 = np.average(xr1[crop_index:])
avg_cropped_50 = np.average(xr[crop_index:])
avg_cropped_75 = np.average(yr[crop_index:])
avg_cropped_100 = np.average(zr1[crop_index:])
avg_cropped_150 = np.average(zr[crop_index:])
avg_cropped_200 = np.average(zr2[crop_index:])
avg_cropped_225 = np.average(zr3[crop_index:])
avg_cropped_245 = np.average(zr4[crop_index:])

compaction_25 = xr1[0]/np.average(xr1[crop_index:])
compaction_50 = xr[0]/np.average(xr[crop_index:])
compaction_75 = yr[0]/np.average(yr[crop_index:])
compaction_100 = zr1[0]/np.average(zr1[crop_index:])
compaction_150 = zr[0]/np.average(zr[crop_index:])
compaction_200 = zr2[0]/np.average(zr2[crop_index:])
compaction_225 = zr3[0]/np.average(zr3[crop_index:])
compaction_245 = zr4[0]/np.average(zr4[crop_index:])


initial25 = xr1[0]
initial50 = xr[0]
initial75 = yr[0]
initial100 = zr1[0]
initial150 = zr[0]
initial200 = zr2[0]
initial225 = zr3[0]
initial245 = zr4[0]

'''c50 = initial50/xr
c75 = initial75/yr
c100 = initial100/zr1
c150 = initial150/zr


plt.plot(time_loopavg_arr,c50,label='50')
plt.legend(title='No. of Dimers')

plt.plot(time_loopavg_arr,c75,label='75')
plt.legend(title='No. of Dimers')

plt.plot(time_loopavg_arr,c100,label='100')
plt.legend(title='No. of Dimers')

plt.plot(time_loopavg_arr,c150,label='150')
plt.xlabel('Timesteps')
plt.ylabel('Compaction')
plt.legend(title='No. of Dimers')
plt.show()'''

N_d_arr = [25,50,75,100,150,200,225,245]





'''print("final compaction_25 rog = ", compaction_25)
print("final compaction_50 rog = ", compaction_50)
print("final compaction_75 rog = ", compaction_75)
print("final compaction_100 rog = ", compaction_100)
print("final compaction_150 rog = ", compaction_150)
print("final compaction_200 rog = ", compaction_200)

compaction_r = [compaction_25,compaction_50,compaction_75,compaction_100,compaction_150,compaction_200]

plt.plot(N_d_arr,compaction_r,marker='*')
plt.xlabel('Number of dimers')
plt.ylabel('Compaction ROG')

plt.show()'''


cropped_25e = x1[crop_index:]
cropped_50e = x[crop_index:]
cropped_75e = y[crop_index:]
cropped_100e = z1[crop_index:]
cropped_150e = z[crop_index:]
cropped_200e = z2[crop_index:]
cropped_225e = z3[crop_index:]
cropped_245e = z4[crop_index:]


avg_cropped_25e = np.average(x1[crop_index:])
avg_cropped_50e = np.average(x[crop_index:])
avg_cropped_75e = np.average(y[crop_index:])
avg_cropped_100e = np.average(z1[crop_index:])
avg_cropped_150e = np.average(z[crop_index:])
avg_cropped_200e = np.average(z2[crop_index:])
avg_cropped_225e = np.average(z3[crop_index:])
avg_cropped_245e = np.average(z4[crop_index:])

compaction_25e = x1[0]/np.average(x1[crop_index:])
compaction_50e = x[0]/np.average(x[crop_index:])
compaction_75e = y[0]/np.average(y[crop_index:])
compaction_100e = z1[0]/np.average(z1[crop_index:])
compaction_150e = z[0]/np.average(z[crop_index:])
compaction_200e = z2[0]/np.average(z2[crop_index:])
compaction_225e = z3[0]/np.average(z3[crop_index:])
compaction_245e = z4[0]/np.average(z4[crop_index:])

print("final compaction_25 e2e = ", compaction_25e)
print("final compaction_50 e2e = ", compaction_50e)
print("final compaction_75 e2e = ", compaction_75e)
print("final compaction_100 e2e = ", compaction_100e)
print("final compaction_150 e2e = ", compaction_150e)
print("final compaction_200 e2e = ", compaction_200e)
print("final compaction_225 e2e = ", compaction_225e)
print("final compaction_245 e2e = ", compaction_245e)


compaction_e = [compaction_25e,compaction_50e,compaction_75e,compaction_100e,compaction_150e,compaction_200e,compaction_225e,compaction_245e]

plt.plot(N_d_arr,compaction_e,marker='*')
plt.xlabel('Number of dimers')
plt.ylabel('Compaction E2E')

plt.show()



croppedl_25 = xl1[crop_index:]
croppedl_50 = xl[crop_index:]
croppedl_75 = yl[crop_index:]
croppedl_100 = zl1[crop_index:]
croppedl_150 = zl[crop_index:]
croppedl_200 = zl2[crop_index:]
croppedl_225 = zl3[crop_index:]
croppedl_245 = zl4[crop_index:]


avg_croppedl_25 = np.average(croppedl_25)
avg_croppedl_50 = np.average(croppedl_50)
avg_croppedl_75 = np.average(croppedl_75)
avg_croppedl_100 = np.average(croppedl_100)
avg_croppedl_150 = np.average(croppedl_150)
avg_croppedl_200 = np.average(croppedl_200)
avg_croppedl_225 = np.average(croppedl_225)
avg_croppedl_245 = np.average(croppedl_245)

print("average loop length 25 = ", avg_croppedl_25)
print("average loop length 50 = ", avg_croppedl_50)
print("average loop length 75 = ", avg_croppedl_75)
print("average loop length 100 = ", avg_croppedl_100)
print("average loop length 150 = ", avg_croppedl_150)
print("average loop length 200 = ", avg_croppedl_200)
print("average loop length 225 = ", avg_croppedl_225)
print("average loop length 245 = ", avg_croppedl_245)



avg_ll = [avg_croppedl_25,avg_croppedl_50,avg_croppedl_75,avg_croppedl_100,avg_croppedl_150,avg_croppedl_200,avg_croppedl_225,avg_croppedl_245]

plt.plot(N_d_arr,avg_ll,marker='*')
plt.xlabel('Number of dimers')
plt.ylabel('Average loop length')

plt.show()


avg_e2e = [avg_cropped_25e,avg_cropped_50e,avg_cropped_75e,avg_cropped_100e,avg_cropped_150e,avg_cropped_200e,avg_cropped_225e,avg_cropped_245e]

# AVG LOOP LENGTH + E2E

combined_compaction_25 = 1/(np.average(x1[crop_index:])+np.average(xl1[crop_index:]))
combined_compaction_50 = 1/(np.average(x[crop_index:])+np.average(xl[crop_index:]))
combined_compaction_75 = 1/(np.average(y[crop_index:])+np.average(yl[crop_index:]))
combined_compaction_100 = 1/(np.average(z1[crop_index:])+np.average(zl1[crop_index:]))
combined_compaction_150 = 1/(np.average(z[crop_index:])+np.average(zl[crop_index:]))
combined_compaction_200 = 1/(np.average(z2[crop_index:])+np.average(zl2[crop_index:]))
combined_compaction_225 = 1/(np.average(z3[crop_index:])+np.average(zl3[crop_index:]))
combined_compaction_245 = 1/(np.average(z4[crop_index:])+np.average(zl4[crop_index:]))

print("final combined_compaction_25 = ", combined_compaction_25)
print("final combined_compaction_50 = ", combined_compaction_50)
print("final combined_compaction_75 = ", combined_compaction_75)
print("final combined_compaction_100 = ", combined_compaction_100)
print("final combined_compaction_150 = ", combined_compaction_150)
print("final combined_compaction_200 = ", combined_compaction_200)
print("final combined_compaction_225 = ", combined_compaction_225)
print("final combined_compaction_225 = ", combined_compaction_245)



combined_compaction = [combined_compaction_25,combined_compaction_50,combined_compaction_75,combined_compaction_100,combined_compaction_150,combined_compaction_200,combined_compaction_225,combined_compaction_245]

plt.plot(N_d_arr,combined_compaction,marker='*')
plt.xlabel('Number of dimers')
plt.ylabel('combined_compaction')
plt.title('Combined compaction (E2E+Avg Loop Length)')

plt.show()





# define the function to fit
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

def func_log(x, a, b):
    return a * np.log(x) + b

def langmuir_eq(x, K, qm):
    return (qm * K * x) / (1 + K * x)

def linear(x, a, b):
    return a*x + b





popt2, pcov2 = curve_fit(func_log, N_d_arr, avg_ll)

# plot the results
plt.scatter(N_d_arr, avg_ll)
plt.plot(N_d_arr, func_log(N_d_arr, *popt2), 'b-', label='Avg Loop Length')
plt.xlabel('Number of dimers')
plt.ylabel('avg loop length')
plt.legend()
plt.show()


popt3, pcov3 = curve_fit(func_log, N_d_arr, compaction_e)

# plot the results
plt.scatter(N_d_arr, compaction_e)
plt.plot(N_d_arr, func_log(N_d_arr, *popt3), 'b-', label='Compaction in E2E distance')
plt.xlabel('Number of dimers')
plt.ylabel('e2e compaction')
plt.legend()
plt.show()


# fit the curve
popt, pcov = curve_fit(func_log, N_d_arr, combined_compaction)

# plot the results
plt.scatter(N_d_arr, combined_compaction)
plt.plot(N_d_arr, func_log(N_d_arr, *popt), 'b-', label='Compaction (E2E+loop)')
plt.xlabel('Number of dimers')
plt.ylabel('combined_compaction')
plt.legend()
plt.show()

popt4, pcov4 = curve_fit(func_log, N_d_arr, avg_e2e)

# plot the results
plt.scatter(N_d_arr, avg_e2e)
plt.plot(N_d_arr, func_log(N_d_arr, *popt4), 'b-', label='Average end to end distance')
plt.xlabel('Number of dimers')
plt.ylabel('end-to-end distance')
plt.legend()
plt.show()



