import random
import numpy as np
#import matplotlib.pyplot as plt 
import math
from tqdm import tqdm as tqdm
from numba import njit
import warnings
warnings.filterwarnings("ignore")



#Parameters

N   = 5000
Nd_arr = []
n = int(input("Enter number of elements : "))
for i in range(0, n):
    ele = int(input("{}/{})No. of dimers: ".format(i+1,n)))
    Nd_arr.append(ele)
#N_d = int(input("No. of dimers: "))
Nd_g = 1
#rho = 0.1     #density
kB_T   = 1
dt = 0.001
zeta = 1
kappa = 1000
kappa_dp = 50
kappa_dd = 100
lb = 1 
lb_d = 1
rc2 = np.cbrt(2)
#frequency = 5000
force_strength = int(input("Force (def. take 350): "))
timesteps = int(input("Timesteps: "))
start = 1
dissociation = 'y'    #dissociate on 'y'
rate_d = 10000
rate_a = 10000
loopavg_time = 1000


for q in Nd_arr:
    q = int(q)
    new_path = r'D:\Study Material\College\MS Thesis\Simulations\n500_nd'+str(q)+'_ndg'+str(Nd_g)+'_ts'+str(timesteps/1000000)+'M_rd'+str(rate_d/1000)+'k_ra'+str(rate_a/1000)+'k.xyz'
    file = open(new_path, 'w')
    file.close()
    
    
    def intialize_pos(x,y,z,N):
        for i in range(N):
            x = i*lb
            y = 0
            z = 0
        return x,y,z
    
    
    def intialize_vel(Vx,Vy,Vz,N):
      KE_temp=0.0
    
      for i in range(N):
        Vx[i]=random.random()-0.5
        Vy[i]=random.random()-0.5
        Vz[i]=random.random()-0.5
        KE_temp = KE_temp +  (1/2)*(Vx[i]**2+Vy[i]**2+Vz[i]**2)
    
      Vx_sum=sum(Vx)
      Vy_sum=sum(Vy)
      Vz_sum=sum(Vz)
    
      for i in range(N):
        Vx[i]=Vx[i]-(Vx_sum/N)
        Vy[i]=Vy[i]-(Vy_sum/N)
        Vz[i]=Vz[i]-(Vz_sum/N)
    
      for i in range(N):
        Vx[i]=(math.sqrt((3*N*kB_T)/(2*KE_temp)))*Vx[i]
        Vy[i]=(math.sqrt((3*N*kB_T)/(2*KE_temp)))*Vy[i]
        Vz[i]=(math.sqrt((3*N*kB_T)/(2*KE_temp)))*Vz[i]
      return Vx, Vy, Vz
    
    def Kinetic_E(Vx,Vy,Vz):
      KE=0
      for i in range(N):
        KE+=(1/2)*((Vx[i]**2)+(Vy[i]**2)+(Vz[i]**2))
      return KE
    
    @njit
    def force(fx,fy,fz,x,y,z,rc2): 
      ulj=0.0
      for i in range(N):
        fx[i]=0.0
        fy[i]=0.0
        fz[i]=0.0
    
      for i in range(N-1):
        j = (i+1)
        dx = x[i]-x[j]
        dy = y[i]-y[j]
        dz = z[i]-z[j]
    
        dr2 = dx**2 + dy**2 + dz**2
        dr  = math.sqrt(dr2)
        ulj   = ulj+(dr-lb)**2
        f   = kappa*(1-lb/dr)   #F_stretch
        
        Fx = -f*dx
        Fy = -f*dy 
        Fz = -f*dz
    
        fx[i] += Fx
        fy[i] += Fy
        fz[i] += Fz
    
        fx[j] -= Fx
        fy[j] -= Fy
        fz[j] -= Fz
    
      ulj = ulj*kappa/2   
    
      #Ulj
      for i in range(N-1):
        for j in range(i+1,N):
            x_diff_n=x[i]-x[j]
            y_diff_n=y[i]-y[j]
            z_diff_n=z[i]-z[j]
    
            r2=((x_diff_n**2)+(y_diff_n**2)+(z_diff_n**2))
    
            if(r2 < rc2 and r2!=0 and j != i+1 and j != i-1):
              force=48.0*(((1.0/r2)**7)-((1.0/2.0)*((1.0/r2)**4)))          #F_LJ
              fox= force*x_diff_n 
              foy= force*y_diff_n 
              foz= force*z_diff_n 
    
              r6 = r2**3
              rc6=rc2**3
              ulj += 4.0*(1.0/r6)*(1/r6 - 1.0)-(4.0*(1/rc6)*(1.0/rc6 - 1.0))
              #print("passed")
            else:
              fox=0.0
              foy=0.0
              foz=0.0
    
    
            fx[i] = fx[i] + fox
            fx[j] = fx[j] - fox
                
                
            fy[i] = fy[i] + foy
            fy[j] = fy[j] - foy
                
                
            fz[i] = fz[i] + foz
            fz[j] = fz[j] - foz
            
    
      return fx,fy,fz,ulj
    
    def force_pull_left(fx,fy,fz,x,y,z,rc2,selected_mono1,selected_mono2):
        
        #Parameters inside function
        x_diff_tan1 = [0.0]*N
        y_diff_tan1 = [0.0]*N
        z_diff_tan1 = [0.0]*N
        
        r2_1 = [0.0]*N
        r_1 = [0.0]*N
    
        
        force1 = [0.0]*N
        fox1 = [0.0]*N
        foy1 = [0.0]*N
        foz1 = [0.0]*N
        
        for j in range(len(selected_mono1)):
            if ((selected_mono1[j]-1) in selected_mono2 or selected_mono1[j]==0):
                continue
            else:
                x_diff_tan1[j]=x[selected_mono1[j]]-x[selected_mono1[j]-1]
                y_diff_tan1[j]=y[selected_mono1[j]]-y[selected_mono1[j]-1]
                z_diff_tan1[j]=z[selected_mono1[j]]-z[selected_mono1[j]-1]
                
                
                r2_1[j]=((x_diff_tan1[j]**2)+(y_diff_tan1[j]**2)+(z_diff_tan1[j]**2))
                r_1[j] = np.sqrt(r2_1[j])
                
                
                force1[j] = (force_strength/r_1[j])
                fox1[j] = force1[j]*x_diff_tan1[j] 
                foy1[j] = force1[j]*y_diff_tan1[j] 
                foz1[j] = force1[j]*z_diff_tan1[j] 
                  
                
                fx[selected_mono1[j]-1] = fx[selected_mono1[j]-1] + fox1[j]
                fy[selected_mono1[j]-1] = fy[selected_mono1[j]-1] + foy1[j]
                fz[selected_mono1[j]-1] = fz[selected_mono1[j]-1] + foz1[j]
    
        
        return fx,fy,fz
    
    def force_pull_right(fx,fy,fz,x,y,z,rc2,selected_mono1,selected_mono2):
        
        #Parameters inside function
    
        x_diff_tan2 = [0.0]*N
        y_diff_tan2 = [0.0]*N
        z_diff_tan2 = [0.0]*N
        
        r2_2 = [0.0]*N
        r_2 = [0.0]*N
        
        
        force2 = [0.0]*N
        fox2 = [0.0]*N
        foy2 = [0.0]*N
        foz2 = [0.0]*N   
        
        for j in range(len(selected_mono1)):
            if ((selected_mono2[j]+1) in selected_mono1 or selected_mono2[j]==499):
                continue
            else:
                x_diff_tan2[j]=x[selected_mono2[j]]-x[selected_mono2[j]+1]
                y_diff_tan2[j]=y[selected_mono2[j]]-y[selected_mono2[j]+1]
                z_diff_tan2[j]=z[selected_mono2[j]]-z[selected_mono2[j]+1]
                
                r2_2[j]=((x_diff_tan2[j]**2)+(y_diff_tan2[j]**2)+(z_diff_tan2[j]**2))
                r_2[j] = np.sqrt(r2_2[j])
                
                force2[j] = (force_strength/r_2[j])
                fox2[j] = force2[j]*x_diff_tan2[j] 
                foy2[j] = force2[j]*y_diff_tan2[j] 
                foz2[j] = force2[j]*z_diff_tan2[j]            
            
                fx[selected_mono2[j]+1] = fx[selected_mono2[j]+1] + fox2[j]
                fy[selected_mono2[j]+1] = fy[selected_mono2[j]+1] + foy2[j]
                fz[selected_mono2[j]+1] = fz[selected_mono2[j]+1] + foz2[j]
        
        return fx,fy,fz
        
    def distance_check_1(x,y,z,xd,yd,zd,selected_mono1,i):
        
        #Parameters
        x_dist_1 = [0.0]*len(selected_mono1)
        y_dist_1 = [0.0]*len(selected_mono1)
        z_dist_1 = [0.0]*len(selected_mono1)
        
        r2_dist_1 = [0.0]*len(selected_mono1)
        r_dist_1 = [0.0]*len(selected_mono1)
        
        
        x_dist_2 = [0.0]*len(selected_mono1)
        y_dist_2 = [0.0]*len(selected_mono1)
        z_dist_2 = [0.0]*len(selected_mono1)
        
        r2_dist_2 = [0.0]*len(selected_mono1)
        r_dist_2 = [0.0]*len(selected_mono1)
        
        for j in range(len(selected_mono1)):
            if ((selected_mono1[j]-1) in selected_mono2 or selected_mono1[j]==0):
                continue
            else:
                x_dist_1[j] = x[selected_mono1[j]] - xd[j][0]
                y_dist_1[j] = y[selected_mono1[j]] - yd[j][0]
                z_dist_1[j] = z[selected_mono1[j]] - zd[j][0]
                
                r2_dist_1[j] = x_dist_1[j]**2 + y_dist_1[j]**2 + z_dist_1[j]**2
                r_dist_1[j] = np.sqrt(r2_dist_1[j])
                
                
                x_dist_2[j] = x[selected_mono1[j]-1] - xd[j][0]
                y_dist_2[j] = y[selected_mono1[j]-1] - yd[j][0]
                z_dist_2[j] = z[selected_mono1[j]-1] - zd[j][0]
                
                r2_dist_2[j] = x_dist_2[j]**2 + y_dist_2[j]**2 + z_dist_2[j]**2
                r_dist_2[j] = np.sqrt(r2_dist_2[j])
                
                if(r_dist_1[j]>=r_dist_2[j]):
                    selected_mono1[j] = selected_mono1[j]-1
        
        return selected_mono1
    
    def distance_check_2(x,y,z,xd,yd,zd,selected_mono2,i):
    
        #Parameters
    
        x_dist_1 = [0.0]*len(selected_mono1)
        y_dist_1 = [0.0]*len(selected_mono1)
        z_dist_1 = [0.0]*len(selected_mono1)
        
        r2_dist_1 = [0.0]*len(selected_mono1)
        r_dist_1 = [0.0]*len(selected_mono1)
        
        
        x_dist_2 = [0.0]*len(selected_mono1)
        y_dist_2 = [0.0]*len(selected_mono1)
        z_dist_2 = [0.0]*len(selected_mono1)
        
        r2_dist_2 = [0.0]*len(selected_mono1)
        r_dist_2 = [0.0]*len(selected_mono1)
        
        for j in range(len(selected_mono1)):
            if ((selected_mono2[j]+1) in selected_mono1 or selected_mono2[j]==499):
                continue
            else:
                x_dist_1[j] = x[selected_mono2[j]] - xd[j][1]
                y_dist_1[j] = y[selected_mono2[j]] - yd[j][1]
                z_dist_1[j] = z[selected_mono2[j]] - zd[j][1]
                
                r2_dist_1[j] = x_dist_1[j]**2 + y_dist_1[j]**2 + z_dist_1[j]**2
                r_dist_1[j] = np.sqrt(r2_dist_1[j])
                
                
                x_dist_2[j] = x[selected_mono2[j]+1] - xd[j][1]
                y_dist_2[j] = y[selected_mono2[j]+1] - yd[j][1]
                z_dist_2[j] = z[selected_mono2[j]+1] - zd[j][1]
                
                r2_dist_2[j] = x_dist_2[j]**2 + y_dist_2[j]**2 + z_dist_2[j]**2
                r_dist_2[j] = np.sqrt(r2_dist_2[j])
                
                if(r_dist_1[j]>=r_dist_2[j]):
                    selected_mono2[j] = selected_mono2[j]+1
        
        return selected_mono2
        
        
    
    def update_position(x,y,z,N,dt):
      for i in range(N):
        x[i]=x[i]+Vx[i]*C
        y[i]=y[i]+Vy[i]*C
        z[i]=z[i]+Vz[i]*C
    
      return x,y,z
    
    def update_velocities(Vx,Vy,Vz,N,dt,A,zeta,eta1,eta2,eta3):
      for i in range(N):
        Vx[i]=A*Vx[i]+fx[i]*dt*(1/2)+B*eta1[i]
        Vy[i]=A*Vy[i]+fy[i]*dt*(1/2)+B*eta2[i]
        Vz[i]=A*Vz[i]+fz[i]*dt*(1/2)+B*eta3[i]
      return Vx,Vy,Vz
    
    def distance(N):
      r2 = 0
      dx = x[0]-x[N-1]
      dy = y[0]-y[N-1]
      dz = z[0]-z[N-1]
      r2 = dx**2 + dy**2 + dz**2
      dist = np.sqrt(r2)
      return dist

    
    """DIMER CODE"""
    
    #Initialise dimer (check where to place (on polymer directly) or (far from polymer, snaps onto polymer))
    def intialize_dimer(xd,yd,zd):
        for j in range(len(selected_mono1)):
            for i in range(2):
                xd[j][i] = x[selected_mono1[j]+i]
                yd[j][i] = y[selected_mono1[j]+i]
                zd[j][i] = z[selected_mono1[j]+i]
        return xd,yd,zd
    
    def force_dimer(fxd,fyd,fzd,xd,yd,zd,rc2): 
       #Parameters inside function
    
        ulj_total=0
        
        dxd = [0.0]*len(selected_mono1)
        dyd = [0.0]*len(selected_mono1)
        dzd = [0.0]*len(selected_mono1)
        
        dr2_d = [0.0]*len(selected_mono1)
        dr_d = [0.0]*len(selected_mono1)
        
        uljd = [0.0]*len(selected_mono1)
        f_d = [0.0]*len(selected_mono1)
        
        Fx_d = [0.0]*len(selected_mono1)
        Fy_d = [0.0]*len(selected_mono1)
        Fz_d = [0.0]*len(selected_mono1)
        
        dx_dp1 = [0.0]*len(selected_mono1)
        dy_dp1 = [0.0]*len(selected_mono1)
        dz_dp1 = [0.0]*len(selected_mono1)
        
        dr2_dp1 = [0.0]*len(selected_mono1)
        dr_dp1 = [0.0]*len(selected_mono1)
        
        f_dp = [0.0]*len(selected_mono1)
        
        Fx_dp1 = [0.0]*len(selected_mono1)
        Fy_dp1 = [0.0]*len(selected_mono1)
        Fz_dp1 = [0.0]*len(selected_mono1)
        
        dx_dp2 = [0.0]*len(selected_mono1)
        dy_dp2 = [0.0]*len(selected_mono1)
        dz_dp2 = [0.0]*len(selected_mono1)
        
        dr2_dp2 = [0.0]*len(selected_mono1)
        dr_dp2 = [0.0]*len(selected_mono1)
        
        f_dp2 = [0.0]*len(selected_mono1)
        
        Fx_dp2 = [0.0]*len(selected_mono1)
        Fy_dp2 = [0.0]*len(selected_mono1)
        Fz_dp2 = [0.0]*len(selected_mono1)
        
        for j in range(len(selected_mono1)):
            
            for i in range(2):
                fxd[j][i]=0
                fyd[j][i]=0
                fzd[j][i]=0
        
          #force between 2 monomers of dimer(cohesin)
          
            dxd[j] = xd[j][0]-xd[j][1]
            dyd[j] = yd[j][0]-yd[j][1]
            dzd[j] = zd[j][0]-zd[j][1]
          
            dr2_d[j] = dxd[j]**2 + dyd[j]**2 + dzd[j]**2
            dr_d[j]  = math.sqrt(dr2_d[j])
            uljd[j]  = (kappa_dd/2)*(dr_d[j]-lb_d)**2
            f_d[j]   = kappa_dd*(1-lb_d/dr_d[j])   #F_stretch
            
            Fx_d[j] = -f_d[j]*dxd[j]
            Fy_d[j] = -f_d[j]*dyd[j] 
            Fz_d[j] = -f_d[j]*dzd[j]
          
            fxd[j][0] += Fx_d[j]
            fyd[j][0] += Fy_d[j]
            fzd[j][0] += Fz_d[j]
          
            fxd[j][1] -= Fx_d[j]
            fyd[j][1] -= Fy_d[j]
            fzd[j][1] -= Fz_d[j]
          
          
          #force between dimer and polymer
          #FOR 1
            dx_dp1[j] = x[selected_mono1[j]]-xd[j][0]
            dy_dp1[j] = y[selected_mono1[j]]-yd[j][0]
            dz_dp1[j] = z[selected_mono1[j]]-zd[j][0]
          
            dr2_dp1[j] = dx_dp1[j]**2 + dy_dp1[j]**2 + dz_dp1[j]**2
            dr_dp1[j]  = math.sqrt(dr2_dp1[j])
            uljd[j]   = (dr_dp1[j]-(0))**2
            f_dp[j]   = kappa_dp*(1-(0))   
            #(dr_dp would be zero, sort that out)
              
            Fx_dp1[j] = -f_dp[j]*dx_dp1[j]
            Fy_dp1[j] = -f_dp[j]*dy_dp1[j]
            Fz_dp1[j] = -f_dp[j]*dz_dp1[j]
              
            fxd[j][0] -= Fx_dp1[j]
            fyd[j][0] -= Fy_dp1[j]
            fzd[j][0] -= Fz_dp1[j]
              
            fx[selected_mono1[j]] += Fx_dp1[j]
            fy[selected_mono1[j]] += Fy_dp1[j]
            fz[selected_mono1[j]] += Fz_dp1[j]
            
            #dr_dp1_arr.append(dr_dp1)
            
            
            #FOR 2
            dx_dp2[j] = x[selected_mono2[j]]-xd[j][1]
            dy_dp2[j] = y[selected_mono2[j]]-yd[j][1]
            dz_dp2[j] = z[selected_mono2[j]]-zd[j][1]
          
            dr2_dp2[j] = dx_dp2[j]**2 + dy_dp2[j]**2 + dz_dp2[j]**2
            dr_dp2[j]  = math.sqrt(dr2_dp2[j])
            uljd[j]   = (dr_dp2[j]-(0))**2
            f_dp2[j]   = kappa_dp*(1-(0))   
            #(dr_dp would be zero, sort that out)
              
            Fx_dp2[j] = -f_dp2[j]*dx_dp2[j]
            Fy_dp2[j] = -f_dp2[j]*dy_dp2[j]
            Fz_dp2[j] = -f_dp2[j]*dz_dp2[j]
              
            fxd[j][1] -= Fx_dp2[j]
            fyd[j][1] -= Fy_dp2[j]
            fzd[j][1] -= Fz_dp2[j]
              
            fx[selected_mono2[j]] += Fx_dp2[j]
            fy[selected_mono2[j]] += Fy_dp2[j]
            fz[selected_mono2[j]] += Fz_dp2[j]
            
            #dr_dp2_arr.append(dr_dp2)
              
          
        for i in range(len(selected_mono1)):
            ulj_total = ulj_total + uljd[i]
        ulj_total = kappa_dp*ulj_total
        
        return fxd,fyd,fzd,ulj_total
    
    def update_position_dimer(xd,yd,zd,dt):
        for j in range(len(selected_mono1)):  
            for i in range(2):
                xd[j][i]=xd[j][i]+Vxd[j][i]*dt+fxd[j][i]*dt*dt*(1/2)
                yd[j][i]=yd[j][i]+Vyd[j][i]*dt+fyd[j][i]*dt*dt*(1/2)
                zd[j][i]=zd[j][i]+Vzd[j][i]*dt+fzd[j][i]*dt*dt*(1/2)
        return xd,yd,zd
    
    def update_velocities_dimer(Vxd,Vyd,Vzd,fxd,fyd,fzd,dt):
        for j in range(len(selected_mono1)): 
            for i in range(2):
                Vxd[j][i]=Vxd[j][i]+fxd[j][i]*dt*(1/2)
                Vyd[j][i]=Vyd[j][i]+fyd[j][i]*dt*(1/2)
                Vzd[j][i]=Vzd[j][i]+fzd[j][i]*dt*(1/2)
        return Vxd,Vyd,Vzd
    
    """END DIMER CODE"""
    
    def ROG():
        Rcm_x=0
        Rcm_y=0
        Rcm_z=0
        Rg2=0
        for i in range(N):
          Rcm_x = Rcm_x + x[i]
          Rcm_y = Rcm_y + y[i]
          Rcm_z = Rcm_z + z[i]
        Rcm_x = Rcm_x/N
        Rcm_y = Rcm_y/N
        Rcm_z = Rcm_z/N
        for i in range(N):
          dx = x[i]-Rcm_x
          dy = y[i]-Rcm_y
          dz = z[i]-Rcm_z 
      
          dr2 = dx**2 + dy**2 + dz**2
          Rg2 +=dr2
        Rg2 = Rg2/N  
        return Rg2
    
    def E2E():
        x_e = 0
        y_e = 0
        z_e = 0
        r_e = 0
        
        x_e = x[N-1] - x[0]
        y_e = y[N-1] - y[0]
        z_e = z[N-1] - z[0]
        
        r_e = x_e**2 + y_e**2 + z_e**2
        e2e = np.sqrt(r_e)
    
        return e2e
    
    
      
    '''ARRAYS'''
    
    dist_arr1 = []
    KE_arr1 = []
    PE_arr1 = []
    PED_arr1 = []
    PE_total_arr1 = []
    RG_arr1 = []
    dr_dp1_arr = []
    dr_dp2_arr = []
    
    position_arr = []
    positiond_arr = []
    frames = []
    
    avgs = []
    Rg_arr = []
    e2e_arr = []
    
    
    time=np.linspace(0,20,20000)
    time_loopavg_arr = []
    
    x   = [0.0]*N
    y   = [0.0]*N
    z   = [0.0]*N
    
    Vx  = [0.0]*N
    Vy  = [0.0]*N
    Vz  = [0.0]*N
    
    fx  = [0.0]*N
    fy  = [0.0]*N
    fz  = [0.0]*N
    
    x2  = [0.0]*N
    y2  = [0.0]*N
    z2  = [0.0]*N
    
    #For Dimer
    
    selected_mono1 = [0]*q
    selected_mono2 = [0]*q
    
    xd   = [ [0.0]*2 for i in range(len(selected_mono1))]
    yd   = [ [0.0]*2 for i in range(len(selected_mono1))]
    zd   = [ [0.0]*2 for i in range(len(selected_mono1))]
    
    Vxd  = [ [0.0]*2 for i in range(len(selected_mono1))]
    Vyd  = [ [0.0]*2 for i in range(len(selected_mono1))]
    Vzd  = [ [0.0]*2 for i in range(len(selected_mono1))]
    
    fxd  = [ [0.0]*2 for i in range(len(selected_mono1))]
    fyd  = [ [0.0]*2 for i in range(len(selected_mono1))]
    fzd  = [ [0.0]*2 for i in range(len(selected_mono1))]
    
    
      #End
    
    for k in range(1000):
        for i in range(len(selected_mono1)):
            selected_mono1[i] = random.randint(1,N-2)
            selected_mono2[i] = selected_mono1[i] + 1
            
        check = [x-1 for x in selected_mono1]
        check1 = [x+1 for x in selected_mono1]
        check2 = [x for x in selected_mono1]
        
        if any(x in selected_mono1 for x in check):
            continue
        elif any(x in selected_mono1 for x in check1):
            continue
        elif any(x in selected_mono1 for x in check2):
            continue
        else:
            break
                
            
        
        
    A=(2-zeta*dt)/(2+zeta*dt)
    B=math.sqrt(kB_T*zeta*dt/2)
    C=2*dt/(2+zeta*dt)
     
    x,y,z=intialize_pos(x,y,z,N)
     #Saving initial positions 
    x2=np.copy(x)
    y2=np.copy(y)
    z2=np.copy(z)
    Vx,Vy,Vz=intialize_vel(Vx,Vy,Vz,N)
    fx,fy,fz,en= force(fx,fy,fz,x,y,z,rc2)
    
    xd,yd,zd = intialize_dimer(xd,yd,zd)
    fxd,fyd,fzd,en_d = force_dimer(fxd,fyd,fzd,xd,yd,zd,rc2)
    
      
    
    
    for i in tqdm(range(timesteps)):
        
          
        #For graphs on python  
        ''' if(i==0):
            fig = plt.figure(figsize=(20,10))  
            ax = fig.add_subplot(111, projection='3d')
            ax.set_xlim3d(-1, 501)
            ax.set_ylim3d(-5, 5)
            ax.set_zlim3d(-5, 5)
                
            ax.plot(x,y,z, color='black')
            ax.scatter(x[selected_mono1],y[selected_mono1],z[selected_mono1], color='blue')
            ax.scatter(x[selected_mono2],y[selected_mono2],z[selected_mono2], color='blue')
            ax.plot(xd,yd,zd, color='red')
            ax.scatter(xd[0],yd[0],zd[0], color='red')
            ax.scatter(xd[1],yd[1],zd[1], color='red')
          #  ax.plot(loop_x,loop_y,loop_z, color='green')
            print('dimer 1 = ',xd[0],yd[0],zd[0])
            print('mono 1 = ',x[selected_mono1],y[selected_mono1],z[selected_mono1])
            print('selected mono -1 = ',x[selected_mono1-1],y[selected_mono1-1],z[selected_mono1-1])
            print('dimer 2 = ',xd[1],yd[1],zd[1])
            print('mono 2 = ',x[selected_mono1+1],y[selected_mono1+1],z[selected_mono1+1])
            print('selected mono +2 = ',x[selected_mono1+2],y[selected_mono1+2],z[selected_mono1+2])
            plt.show()'''
                  
                
        
        
        if(i>=start):
    
            if(len(selected_mono1)>0):
                fx,fy,fz = force_pull_left(fx,fy,fz,x,y,z,rc2,selected_mono1,selected_mono2)
                selected_mono1 = distance_check_1(x,y,z,xd,yd,zd,selected_mono1,i)
                fx,fy,fz = force_pull_right(fx,fy,fz,x,y,z,rc2,selected_mono1,selected_mono2)
                selected_mono2 = distance_check_2(x,y,z,xd,yd,zd,selected_mono2,i)
                
                if(i%loopavg_time==0 and i != 0):
                    loop_lengths = []
                    RG = 0
                    for v in range(len(selected_mono1)):
                        for j in range(N):
                            if(j==selected_mono1[v]):
                                a = j+1
                            if(j==selected_mono2[v]):
                                b = j
                        c = b-a
                        loop_lengths.append(c)
                        #print statement for loops at that time step
                    avg_lp_len = np.average(loop_lengths)
                    avgs.append(avg_lp_len)
                    step_no = i
                    time_loopavg_arr.append(step_no)
                    #try for std dev array
                    
                    RG = ROG()
                    Rg_arr.append(RG)
                    e2e = E2E()
                    e2e_arr.append(e2e)
            
                if(dissociation == 'y'):
                    if(i%rate_d == 0):
                        for p in range(Nd_g):
                            rand_v = random.randint(0,len(selected_mono1)-1)
                            selected_mono1.pop(rand_v)
                            selected_mono2.pop(rand_v)
                            
                            xd.pop(rand_v)
                            yd.pop(rand_v)
                            zd.pop(rand_v)
                            
                            Vxd.pop(rand_v)
                            Vyd.pop(rand_v)
                            Vzd.pop(rand_v)
                            
                            fxd.pop(rand_v)
                            fyd.pop(rand_v)
                            fzd.pop(rand_v)
            
            
            if(i%rate_a == 0):
                for q in range(Nd_g):
                    all_mono = np.linspace(1,N-2,N-2).astype(int)
                    possible_mono1 = [x for x in all_mono if x not in selected_mono1]
                    possible_mono2 = [x for x in possible_mono1 if x+1 not in selected_mono1]
                    possible_mono3 = [x for x in possible_mono2 if x not in selected_mono2]
                    possible_mono_final = [x for x in possible_mono3 if x+1 not in selected_mono2]
                    selected_mono1_new = random.choice(possible_mono_final)
                    selected_mono2_new = selected_mono1_new+1
                    
                    selected_mono1.append(selected_mono1_new)
                    selected_mono2.append(selected_mono2_new)
                    
                    xd.append([x[selected_mono1_new],x[selected_mono2_new]])
                    yd.append([y[selected_mono1_new],y[selected_mono2_new]])
                    zd.append([z[selected_mono1_new],z[selected_mono2_new]])
                    
                    Vxd.append([0.0]*2)
                    Vyd.append([0.0]*2)
                    Vzd.append([0.0]*2)
                    
                    fxd.append([0.0]*2)
                    fyd.append([0.0]*2)
                    fzd.append([0.0]*2)
                    
                   # fxd,fyd,fzd,en_d = force_dimer(fxd,fyd,fzd,xd,yd,zd,rc2)
                   # Vxd,Vyd,Vzd = update_velocities_dimer(Vxd,Vyd,Vzd,fxd,fyd,fzd,dt)
        
            
            
                        
        '''if(i%rate_a==1):
            for j in range(len(selected_mono1)):
                print(selected_mono1[j],selected_mono2[j])'''
                
                    
        eta1=np.random.normal(0,1,N)  #for random force (brownian motion)
        eta2=np.random.normal(0,1,N)
        eta3=np.random.normal(0,1,N)
        Vx,Vy,Vz = update_velocities(Vx,Vy,Vz,N,dt,1,zeta,eta1,eta2,eta3)
        x,y,z = update_position(x,y,z,N,dt)
        fx,fy,fz,en= force(fx,fy,fz,x,y,z,rc2)
        Vx,Vy,Vz = update_velocities(Vx,Vy,Vz,N,dt,A,zeta,eta1,eta2,eta3)
                
        xd,yd,zd = update_position_dimer(xd,yd,zd,dt)
        Vxd,Vyd,Vzd = update_velocities_dimer(Vxd,Vyd,Vzd,fxd,fyd,fzd,dt)
        fxd,fyd,fzd,en_d= force_dimer(fxd,fyd,fzd,xd,yd,zd,rc2)
        Vxd,Vyd,Vzd = update_velocities_dimer(Vxd,Vyd,Vzd,fxd,fyd,fzd,dt)
        
        
        '''dimer1_distance_x = x[selected_mono1] - xd[0]
        dimer1_distance_y = y[selected_mono1] - yd[0]
        dimer1_distance_z = z[selected_mono1] - zd[0]
          
        dimer2_distance_x = x[selected_mono2] - xd[1]
        dimer2_distance_y = y[selected_mono2] - yd[1]
        dimer2_distance_z = z[selected_mono2] - zd[1]
        
        dimer1_distance = np.sqrt(dimer1_distance_x**2 + dimer1_distance_y**2 + dimer1_distance_z**2)
        dimer2_distance = np.sqrt(dimer2_distance_x**2 + dimer2_distance_y**2 + dimer2_distance_z**2)
        
        dimer1_distance_arr.append(dimer1_distance)
        dimer2_distance_arr.append(dimer2_distance)'''
            
            
        
        '''if(selected_mono1 == 0 or selected_mono2 == 499):
            selected_mono1 = 999998
            selected_mono2 = 999999
            for i in range(2):
                xd[i] = 250+i
                yd[i] = -250+i
                zd[i] = -250+i'''
                
                
    
        
        ''' Try to make loop green for multiple dimers case'''
        if(2*i%(20*(timesteps/10000))==0):
            combined_position_arr = []
            for j in range(len(x)):
                '''for m in range(len(selected_mono1)):
                    if(j>=selected_mono1[m] and j<=selected_mono2[m]):
                        combined_position_arr.append(['G',x[j],y[j],z[j]])'''
                if(1>0):
                    combined_position_arr.append(['R',x[j],y[j],z[j]])
            for k in range(len(selected_mono1)):
                for l in range(2):
                    combined_position_arr.append(['B',xd[k][l],yd[k][l],zd[k][l]])
                    
            #Write File
            file = open(new_path, 'a')
            tot_par = N + 2*(len(selected_mono1))
            
            file.write(str(tot_par)+"\n")
            file.write(str(i+1)+"\n")
            for i in range(len(combined_position_arr)):
                file.write("{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                    combined_position_arr[i][0], combined_position_arr[i][1],
                    combined_position_arr[i][2], combined_position_arr[i][3]))
    
    
                
                
    print(new_path)
     
    '''plt.plot(time_loopavg_arr,avgs,label='Average Loop Length')
    plt.xlabel('Timesteps')
    plt.ylabel('Average Loop Length')
    plt.legend()
    plt.show()
    
    plt.plot(time_loopavg_arr,Rg_arr,label='Radius of Gyration')
    plt.xlabel('Timesteps')
    plt.ylabel('Radius of gyration')
    plt.legend()
    plt.show()
    
    plt.plot(time_loopavg_arr,e2e_arr,label='End to End distance')
    plt.xlabel('Timesteps')
    plt.ylabel('End to End distance')
    plt.legend()
    plt.show()'''
    
            
    
    np.savetxt('e2e_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', e2e_arr, delimiter=",")
    np.savetxt('rog_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', Rg_arr, delimiter=",")
    np.savetxt('loop_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', avgs, delimiter=",")
    
    
    
    np.savetxt('fx_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', fx, delimiter=",")
    np.savetxt('fy_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', fy, delimiter=",")
    np.savetxt('fz_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', fz, delimiter=",")
    
    np.savetxt('x_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', x, delimiter=",")
    np.savetxt('y_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', y, delimiter=",")
    np.savetxt('z_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', z, delimiter=",")
    
    np.savetxt('vx_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', Vx, delimiter=",")
    np.savetxt('vy_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', Vy, delimiter=",")
    np.savetxt('vz_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', Vz, delimiter=",")
    
    np.savetxt('fxd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', fxd, delimiter=",")
    np.savetxt('fyd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', fyd, delimiter=",")
    np.savetxt('fzd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', fzd, delimiter=",")
    
    np.savetxt('xd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', xd, delimiter=",")
    np.savetxt('yd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', yd, delimiter=",")
    np.savetxt('zd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', zd, delimiter=",")
    
    np.savetxt('vxd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', Vxd, delimiter=",")
    np.savetxt('vyd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', Vyd, delimiter=",")
    np.savetxt('vzd_nd'+str(q*1)+'_ts'+str(timesteps/1000000)+'M.csv', Vzd, delimiter=",")

np.savetxt('timesteps.csv', time_loopavg_arr, delimiter=",")











'''plt.plot(selected_mono1_arr, color='red', label='selected mono 1')
plt.xlabel('steps')
plt.xlim(0,timesteps)
plt.legend()
#plt.show()

plt.plot(selected_mono2_arr, color='magenta', label='selected mono 2')
plt.xlabel('steps')
plt.xlim(0,timesteps)
plt.legend()
plt.show()

plt.plot(dimer1_distance_arr, label='distance b/w dimer1 and selected mono1')
plt.xlabel('steps')
plt.xlim(0,timesteps)
plt.legend()
plt.show()

plt.plot(dimer2_distance_arr, label='distance b/w dimer2 and selected mono2')
plt.xlabel('steps')
plt.xlim(0,timesteps)
plt.legend()
plt.show()
'''
#print(selected_mono2_arr)
#print(selected_mono1_arr)
'''plt.plot(dist_arr1, color='blue', label='MSD')
plt.xlabel('steps')
plt.ylabel('distance b/w first and last molecule')'''

'''fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z)
ax.scatter(x[selected_mono1],y[selected_mono1],z[selected_mono1])
ax.scatter(x[selected_mono1+1],y[selected_mono1+1],z[selected_mono1+1])
ax.scatter(xd[0],yd[0],zd[0], color='red')
ax.scatter(xd[1],yd[1],zd[1], color='red')
print('dimer 1 = ',xd[0],yd[0],zd[0])
print('mono 1 = ',x[selected_mono1],y[selected_mono1],z[selected_mono1])
print('selected mono -1 = ',x[selected_mono1-1],y[selected_mono1-1],z[selected_mono1-1])
print('dimer 2 = ',xd[1],yd[1],zd[1])
print('mono 2 = ',x[selected_mono1+1],y[selected_mono1+1],z[selected_mono1+1])
print('selected mono +2 = ',x[selected_mono1+2],y[selected_mono1+2],z[selected_mono1+2])
plt.show()'''