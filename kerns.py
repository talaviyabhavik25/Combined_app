import streamlit as st
import numpy as np
import pandas as pd
from scipy.optimize import fsolve,root
import ht
@st.cache_data
def load_data_table():
    url ='http://raw.githubusercontent.com/Ahmedhassan676/htcalc/main/data_tables.csv'

    return pd.read_csv(url)

tube_table = load_data_table().iloc[1:11,1:5]
def get_index(series, value):
            n = 0 
            series = series.reset_index()
            del series['index']
            try:
                n = int(series[series.iloc[:,0]==str(value)].index[0])
            except IndexError:
                n = int(series[series.iloc[:,0]==str(int(value*1000))].index[0]) 
            except NameError: pass
            return n  
def main_kern(U_assumed, Tube_list, Shell_list,HB_data,j_const,Do,thick,L,geo_input_list,dp_sin,dp_tin,s3,assumptions):
    m_t,t1_t,t2_t,rho_t,Cp_t,mu_t,k_t,fouling_t = Tube_list[0], Tube_list[1], Tube_list[2], Tube_list[3], Tube_list[4], Tube_list[5]*0.001, Tube_list[6], Tube_list[7]
    m_s,t1_s,t2_s,rho_s,Cp_s,mu_s,k_s,fouling_s = Shell_list[0], Shell_list[1], Shell_list[2], Shell_list[3], Shell_list[4], Shell_list[5]*0.001, Shell_list[6], Shell_list[7]
    Cp_t = Cp_t*4184
    Cp_s = Cp_s*4184
    dp_sin,dp_tin = dp_sin/10**8,dp_tin/10**8
    L = L/1000
    
    pn = 2 # assumed
    tpitch = assumptions[2]*Do # assumed
    Q, dTlm, ft = HB_data[0], HB_data[1], HB_data[2]
    #U_assumed = 350
    corrected_LMTD = dTlm*ft
    
    velocity_t = 0
    percentage_diff = -1
    iterv = 0
    iteru,iteru2  =0,0
    iterdp = 0
    error_dp_t,error_dp_s = 1,1
    Do_ind=Do_ind_init = get_index(tube_table.iloc[1:10,0],Do)
    while  (percentage_diff < 10 ) and iteru2 <= 20:
        A_required = Q/(corrected_LMTD*U_assumed)
        tn = int(A_required/(np.pi*L*Do*0.001*s3))
        while (error_dp_s > 0 or error_dp_t > 0.2 ) and iterdp <= 50:
            Di = (Do - 2*thick)*0.001
            tpitch = 1.25*Do
            while (percentage_diff < 10 ) and iteru <= 20: #or percentage_diff > 30
                A_required = Q/(corrected_LMTD*U_assumed)
                tn = int(A_required/(np.pi*L*Do*0.001*s3))
                #st.write(U_assumed,A_required)
                selected_velocity = assumptions[0]
                while velocity_t < selected_velocity and iterv <= 10:
                    
                    
                    
                    cross_A=(np.pi*0.25*(Di**2))*(tn/pn)
                   
                    velocity_t = m_t/(rho_t*3600*cross_A)
                    
                    bundle = ht.hx.DBundle_for_Ntubes_Phadkeb(tn, Do/1000, tpitch/1000, pn, angle=30)
                    m,c=0.027,0.0446
                    shell_D = int(bundle*1000+(0.027*bundle+0.0446)*1000)
                    Ret=(rho_t*velocity_t*Di)/mu_t
                    #st.write(Ret)
                    if Ret >= 2300:
                        f_t =1/(1.58*np.log(Ret)-3.28)**2 # valid for Re 2300 to 5,000,000 and Pr 0.5 to 2000
                    else: f_t = 64/(4*Ret)
                    #port_1 = f_t*L*pn/Di
                    #port_2 = rho_t*(velocity_t**2)/2
                    #dp_t = (4*(port_1)+4*(pn))*port_2*0.000010197
                    Pr = Cp_t*mu_t/k_t
                    Nu = ((0.5*f_t*(Ret-1000)*Pr))/(1+12.7*((0.5*f_t)**0.5)*((Pr**(2/3))-1)) # valid for Re 2300 to 5,000,000 (Gnielinski)
                    h_t = Nu *k_t/Di
                    if velocity_t < selected_velocity and pn <8:
                        pn +=2
                        tn = int(A_required/(np.pi*L*Do*0.001*s3))
                    
                        cross_A=(np.pi*0.25*(Di**2))*(tn/pn)
                        velocity_t = m_t/(rho_t*3600*cross_A)
                    #st.write(pn,velocity_t,shell_D,h_t,Ret,Nu)
                    #elif pn == 8 and velocity_t < selected_velocity:
                        #selected_velocity = 1
                        #pn = 2
                        #iterv = 0
                    
                    if iterv ==9:
                        st.write('velocity iteration failed')
                    iterv +=1 
                bundle = ht.hx.DBundle_for_Ntubes_Phadkeb(tn, Do/1000, tpitch/1000, pn, angle=30)
                m,c=0.027,0.0446
                shell_D = int(bundle*1000+(0.027*bundle+0.0446)*1000)
                Ret=(rho_t*velocity_t*Di)/mu_t
                f_t =1/(1.58*np.log(Ret)-3.28)**2 # valid for Re 2300 to 5,000,000 and Pr 0.5 to 2000
                
                Pr = Cp_t*mu_t/k_t
                Nu = ((0.5*f_t*(Ret-1000)*Pr))/(1+12.7*((0.5*f_t)**0.5)*((Pr**(2/3))-1)) # valid for Re 2300 to 5,000,000 (Gnielinski)
                h_t = Nu *k_t/Di
                b_space = shell_D/5 # assumed
                C = tpitch-Do
                As = (shell_D*b_space*C)/(tpitch*1000000)
                #st.warning(As)
                #st.warning(shell_D)
                Gs = m_s/(As*3600)
                velocity_s = Gs/rho_s 
                pitch = assumptions[1]
                if pitch == 'square' or 'rotated square 45':
                    De = 4*(((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.25))/(3.14*Do*0.001)
                else:
                    De = 8*(0.43301*((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.125))/(3.14*Do*0.001)  
                
                h_shell = (0.36*((De*Gs/mu_s)**0.55)*((Cp_s*mu_s/k_s)**(1/3)))*k_s/De
                d_ratio = Do/(Di*1000)
                Uc = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell))
                Ud = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell)+fouling_s+(d_ratio*fouling_t))
                
                percentage_diff = ((Ud-U_assumed)/U_assumed)*100
                #st.write(percentage_diff)
                if percentage_diff < 10 : # or percentage_diff > 30:
                    
                    U_assumed = U_assumed*0.9
                    #st.write(A_required,U_assumed)
                    #percentage_diff = ((Ud-U_assumed)/U_assumed)*100
                    
                    #st.write(U_assumed)
                    
                if iteru ==20:
                        st.write('U iteration failed')
                iteru +=1 
                #st.write(percentage_diff,U_assumed)
            #st.write(U_assumed,iteru)
            #st.write(A_required)
            As = (shell_D*b_space*C)/(tpitch*1000000)
            Gs = m_s/(As*3600)
            tn = int(A_required/(np.pi*L*Do*0.001*s3)) 
            cross_A=(np.pi*0.25*(Di**2))*(tn/pn)
            velocity_t = m_t/(rho_t*3600*cross_A)
            Res = (De*Gs)/mu_s
            f = np.exp(0.576-(0.19*np.log(Res)))
            Nb = int((L*1000/b_space)-1)
            #st.warning(Nb)
            dp_s = ((f*(Gs**2)*(Nb+1)*shell_D)/(2*rho_s*De*1000))*0.000010197
            
            Ret=(rho_t*velocity_t*Di)/mu_t
            if Ret >= 2300:
                f_t =1/(1.58*np.log(Ret)-3.28)**2 # valid for Re 2300 to 5,000,000 and Pr 0.5 to 2000
            else: f_t = 64/(4*Ret)
            
            port_1 = f_t*L*pn/Di
            port_2 = rho_t*(velocity_t**2)/2
            #st.write(f_t,L,pn,Di,rho_t,velocity_t,rho_s)
            dp_t = (4*(port_1)+4*(pn))*port_2*0.000010197
            
            error_dp_s = dp_s-(dp_sin)
            error_dp_t = dp_t-(dp_tin)
            #st.write(error_dp_t)
            #sst.write(dp_s)
            #st.write(dp_sin)
            #if error_dp_s > 0 and error_dp_t > 0:
                #b_space +=(shell_D/5)*0.1
                #Do_ind += 1
                #Do = float(tube_table.iloc[Do_ind,0])
            if error_dp_s > 0 and iterdp < 50:
                b_space +=(shell_D/5)*0.1
            if error_dp_t > 0.2 and iterdp < 50:
                Do_ind += 1
                Do = float(tube_table.iloc[Do_ind,0])
                
                #Di = (Do - 2*thick)*0.001
               
            if iterdp ==50:
                st.write('dp iteration failed')
            #st.write(dp_s,dp_t)
            #st.write(dp_s,dp_t,Ud,percentage_diff,velocity_t,Do_ind,Di)
            #float(tube_table.iloc[2,0])
            iterdp +=1
            #st.warning('errt is {} while Do is {}'.format(error_dp_t,Do))
            #st.warning(dp_tin)
            #st.write(dp_s,Res,Gs,velocity_s,De)
        Ret=(rho_t*velocity_t*Di)/mu_t
        if Ret >= 2300:
            f_t =1/(1.58*np.log(Ret)-3.28)**2 # valid for Re 2300 to 5,000,000 and Pr 0.5 to 2000
        else: f_t = 64/(4*Ret)
        
        Pr = Cp_t*mu_t/k_t
        Nu = ((0.5*f_t*(Ret-1000)*Pr))/(1+12.7*((0.5*f_t)**0.5)*((Pr**(2/3))-1)) # valid for Re 2300 to 5,000,000 (Gnielinski)
        h_t = Nu *k_t/Di
        #b_space = shell_D/5 # assumed
        C = tpitch-Do
        As = (shell_D*b_space*C)/(tpitch*1000000)
        #st.warning(As)
        #st.warning(shell_D)
        Gs = m_s/(As*3600)
        velocity_s = Gs/rho_s 
        pitch = assumptions[1]
        if pitch == 'square' or 'rotated square 45':
            De = 4*(((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.25))/(3.14*Do*0.001)
        else:
            De = 8*(0.43301*((tpitch*0.001)**2)-(3.14*((Do*0.001)**2)*0.125))/(3.14*Do*0.001)  
        
        h_shell = (0.36*((De*Gs/mu_s)**0.55)*((Cp_s*mu_s/k_s)**(1/3)))*k_s/De
        d_ratio = Do/(Di*1000)
        Uc = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell))
        Ud = 1/((d_ratio/h_t)+(Do*0.001*np.log(d_ratio)/(2*60))+(1/h_shell)+fouling_s+(d_ratio*fouling_t))
        percentage_diff = ((Ud-U_assumed)/U_assumed)*100
        if percentage_diff < 10:
            U_assumed = U_assumed*0.9
            velocity_t = 0
            percentage_diff = -1
            iterv = 0
            iteru,iteru2  =0,0
            iterdp = 0
            error_dp_t,error_dp_s = 1,1
             
        iteru2 +=1  
        #st.write(percentage_diff,U_assumed,Ud,h_shell,h_t)
        #st.write(shell_D,b_space,tn,A_required)
    tn = int(A_required/(np.pi*L*Do*0.001*s3))        
    pitch = assumptions[1]
    b_cut = 25
    #st.write(percentage_diff,U_assumed,Ud)
    #st.write(dp_s,dp_t,Ud,percentage_diff,velocity_t,Do)
    geo_input_df = pd.DataFrame(index=geo_input_list)
    geo_input_df.loc[['Shell D','Baffle Spacing','Number of baffles','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut'],'Kern_summary']=geo_input_df.loc[['Shell D','Baffle Spacing','Number of baffles','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut'],'Bell_summary'] =  [int(shell_D),b_space,int(Nb),Do,Di*1000,L,int(tn),pn,tpitch,pitch,b_cut]
    geo_list =  [int(tn),pn,Do,Di*1000,pitch,tpitch,L*1000,b_space,b_cut,int(shell_D)/1000]
    #st.write(geo_list)
    return geo_list, geo_input_df