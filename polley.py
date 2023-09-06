import streamlit as st
import numpy as np
import pandas as pd
from scipy.optimize import fsolve,root
import ht
print('------------------begin------------------')
def main_polley(Tube_list, Shell_list,HB_data,j_const,Do,thick,geo_input_list,pitch_ratio,pitch,dp_s,dp_t,s3):
    m_t,t1_t,t2_t,rho_t,Cp_t,mu_t,k_t,fouling_t = Tube_list[0], Tube_list[1], Tube_list[2], Tube_list[3], Tube_list[4], Tube_list[5], Tube_list[6]/1.163, Tube_list[7]*1.163
    m_s,t1_s,t2_s,rho_s,Cp_s,mu_s,k_s,fouling_s = Shell_list[0], Shell_list[1], Shell_list[2], Shell_list[3], Shell_list[4], Shell_list[5], Shell_list[6]/1.163, Shell_list[7]*1.163
    
    Q, dTlm, ft = HB_data[0], HB_data[1], HB_data[2]
    #initialise baffle cut
    b_cut = 25
    #initialise shell diameter
    if dp_s < 0.5*(10**8) and dp_t < 0.5*(10**8):
        error_selected = 0.01
    else: error_selected = 0.1
    pn = 2
    error =2
    Di = (Do - 2*thick)
    # assumed pitch ratio
    tpitch = pitch_ratio*Do
    if pitch == 'square':
      t_p_angle = 90
    elif pitch == 'rotated square 45':
      t_p_angle = 45
    else:
      t_p_angle = 30
    def initialise(b_cut,shell_D, pn, L,t_p_angle,Do,pitch_ratio,pitch,s3):
            A_ratio = 1.5
            if pitch == 'square':
                t_p_angle = 90
            elif pitch == 'rotated square 45':
                t_p_angle = 45
            else:
                t_p_angle = 30
            while A_ratio > 1:
                #Do = 19.05
                Di = Do- 2*thick
                #print('Di in init is {}'.format(Di))
                #t_p_angle = 30
                shell_D = shell_D
                b_cut= b_cut
                L= L/1000
                L_tb = 0.4 # Diametral Tube-Baffle Clearance
                D_sb = 3.1+0.004*shell_D
                L_s = 0.1*shell_D
                L_eff = L-2*L_s/1000 # Effective tube length
                #Estimate: tube count and baffle spacing
                tpitch = pitch_ratio*Do /1000
                t_p = pitch_ratio*Do
                if t_p_angle == 45:
                    t_p_effective = t_p/(2**0.5)
                else: t_p_effective = t_p
                if t_p_angle == 30:
                    t_arrg = np.sqrt(3)/2
                elif t_p_angle == 45:
                    t_arrg = 1/np.sqrt(2)
                else : t_arrg = 1
                
                DBundle = shell_D/1000 - 2*ht.shell_clearance(DShell=shell_D/1000)
                tn = ht.hx.Ntubes(DBundle, Do/1000, tpitch, Ntp=pn, angle=t_p_angle, Method=None)
                print('number of tube in init is {}'.format(tn))
                D_otl = shell_D - (12.5+(shell_D/200))

                theta_Ds = 2*np.arccos(1-(2*b_cut)/100) # Baffle window angle
                theta_CTL = 2*np.arccos(shell_D*(1-2*b_cut/100)/(D_otl-Do)) #b_cut baffle cut, shell_D isnide shell diameter mm
                F_w = (theta_CTL-np.sin(theta_CTL))/(2*np.pi)
                S_wg = (((shell_D/1000)**2)/8)*(theta_Ds-np.sin(theta_Ds)) #Gross window area
                S_wt = tn*F_w*np.pi*((Do/1000)**2)/4 #Window area occupied with tubes
                S_w = S_wg - S_wt # Net Cross flow area
                S_m = S_w
                G_w = (m_s/3600)/np.sqrt(S_m*S_w)
                v_w = G_w/rho_s

                #S_m = (Lb_cut/1000)*((shell_D-D_otl)+(D_otl-Do)*(t_p-Do)/t_p_effective)/1000
                Lb_cut = LB_in = LB_out = ((S_m*1000)/((shell_D-D_otl)+(D_otl-Do)*(t_p-Do)/t_p_effective))*1000
                Lb_cut

                # Calculate shell side correction factors
                Gs = (m_s/3600)/S_m
                Re_s = (Do/1000)*Gs/(mu_s*0.001)
                #print('dp for Re_s '+str(Re_s))
                Pr_s = (mu_s/1000)*Cp_s/k_s*3600
                if Re_s > 100000:
                    Re_temp = Re_s 
                    Re_s = 99999

                # j- ideal factor coefficients
                a1 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['a_{1}'],0).sum()
                a2 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['a_{2}'],0).sum()
                a3 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['a_{3}'],0).sum()
                a4 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['a_{4}'],0).sum()
                a=a3/(1+0.14*((Re_s)**a4))
                #print(a)
                #print(Re_s)
                if Re_s == 99999:
                    Re_s = Re_temp
                j=a1*((1.33/(t_p/Do))**a)*(Re_s**a2)
                #print('dp for a3 '+str(a3))
                #print('dp for a4 '+str(a4))
                #print('value for j '+str(j))
                #print(j)
                mu_s_w = mu_s
                h_s_i = Cp_s*(m_s/S_m)*j*(Pr_s**(-2/3))*((mu_s/mu_s_w)**0.14)
                #print('value of hs in init is {}'.format(h_s_i))

                # 1. correction factor for baffle window flow
                N_b = 1 +int((L-(2*L_s*0.001)-(LB_in+LB_out)*0.001)/(Lb_cut*0.001)) # number of baffles
                theta_CTL = 2*np.arccos(shell_D*(1-2*b_cut/100)/(D_otl-Do)) #b_cut baffle cut, shell_D isnide shell diameter mm
                F_w = (theta_CTL-np.sin(theta_CTL))/(2*np.pi)

                F_c = 1-2*F_w # Fraction of tubes in cross flow
                j_c = 0.55+0.72*F_c

                # 2. correction factor for baffle leakage
                theta_Ds = 2*np.arccos(1-(2*b_cut)/100) # Baffle window angle
                S_sb = (shell_D/1000)*(D_sb/1000)*(np.pi-0.5*theta_Ds)    # Shell to baffle leakage area
                S_tb = (np.pi/4)*((((Do/1000)+(L_tb/1000))**2)-((Do/1000)**2))*tn*(1-F_w)  # Tube to baffle leakage area

                r_L = (S_sb+S_tb)/S_m # ratio of leakage to cross flow
                r_S = S_sb/(S_sb+S_tb) # ratio of shell-baffle to total area

                j_l = 0.44*(1-r_S)+(1-0.44*(1-r_S))*np.exp(-2.2*r_L) # Baffle leakage correction factor

                # 3. correction factor for bundle bypass
                    # number_of_sealing_strips   
                P_p = t_p * t_arrg # Tube row distance in flow direction
                N_TCC = (shell_D/P_p)*(1-(2*b_cut)/100) #N_TCC number of tube rows between baffle
                N_ss = int(N_TCC/6)
                r_ss = N_ss/N_TCC # Nss Number of sealing strips
                #print('N_ss is '+str(N_ss)+' While N_Tcc ia '+str(N_TCC))
                S_b = (Lb_cut/1000)*(shell_D-D_otl-(Do/2))/1000 #bundle pybass area
                if Re_s < 100:
                    C_j = 1.35
                else: C_j = 1.25 # Correlation constant

                if r_ss >= 0.5:
                    j_b = 1
                else: j_b = np.exp(-1*C_j*(S_b/S_m)*(1-((2*r_ss)**(1/3))))

                # 4. Correction factor for adverse temperature gradient
                N_tcw = (0.8/P_p)*(shell_D*b_cut/100-(shell_D-(D_otl-Do))/2)
                N_b = 1 +int((L-(2*L_s*0.001)-(LB_in+LB_out)*0.001)/(Lb_cut*0.001)) # number of baffles
                #print('value for N_b '+str(N_b))
                N_c = (N_tcw +N_TCC)*(1+N_b) # tube rows crossed in entire exchanger
                j_RL = (10/N_c)**0.18
                if Re_s <= 20:
                    j_R = j_RL
                elif Re_s <100:
                    j_R = j_RL+((20-Re_s)/80)*(j_RL-1)
                else: j_R = 1

                # 5. correction factor for unequal baffle spacing
                if Re_s < 100:
                    n1 = 1/3
                else: n1 = 0.6
                j_s =((N_b-1)+(LB_in/Lb_cut)**(1-n1)+(LB_out/Lb_cut)**(1-n1))/((N_b-1)+(LB_in/Lb_cut)+(LB_out/Lb_cut))


                h_shell = j_s * j_R *j_b * j_l * j_c * h_s_i
                #print(h_shell)
                if Re_s > 100000:
                    Re_temp = Re_s 
                    Re_s = 99999
                ### Shell side pressure drop
                b1 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{1}'],0).sum()
                b2 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{2}'],0).sum()
                b3 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{3}'],0).sum()
                b4 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{4}'],0).sum()
                if Re_s == 99999:
                    Re_s = Re_temp
                b = b3/(1+0.14*(Re_s**b4))
                f_s = b1*((1.33/(t_p/Do))**b)*(Re_s**b2)
                dp_shell_ideal = 2*f_s*(Gs**2)*N_TCC*(mu_s_w/mu_s)**0.14/(rho_s)/100000 #N_TCC number of tube rows between baffle
                #print(dp_shell_ideal)

                # 1. correction factor for baffle leakage
                f_p = 0.8-0.15*(1+r_S) # r_s ratio of shell-baffle to total area
                R_L = np.exp(-1.33*(1+r_S)*(r_L**f_p))# r_l ratio of  leakage area to cross flow

                # pressure drop for an ideal window section
                S_wg = (((shell_D/1000)**2)/8)*(theta_Ds-np.sin(theta_Ds)) #Gross window area
                S_wt = tn*F_w*np.pi*((Do/1000)**2)/4 #Window area occupied with tubes
                S_w = S_wg - S_wt # Net Cross flow area
                G_w = (m_s/3600)/np.sqrt(S_m*S_w)
                v_w = G_w/rho_s
                #print('dp for G_w '+str(G_w))
                D_w = 4*S_w/(np.pi*(Do/1000)*tn*F_w+(shell_D/1000)*theta_Ds)

                #pressure drop for turbluent flow in ideal window section
                if Re_s >= 100:
                    dp_window = N_b*R_L*(2+0.6*N_tcw)*(G_w**2)/(2*rho_s)/100000

                    
                else:
                    dp_window = N_b*R_L*((26*G_w*Cp_s/rho_s)*(N_tcw/(tpitch-Do)/1000 + (Lb_cut/1000)/(D_w**2)) + (G_w**2)/rho_s)/100000
                    

                # 2. Correction factor for bundle bypass effect
                if Re_s < 100:
                    C_r = 4.5
                else: C_r = 3.7

                if r_ss >= 0.5:
                    R_b = 1
                else:
                    R_b = np.exp(-1*C_r*(S_b/S_m)*(1-(2*r_ss)**(1/3)))

                # 3. Correction for unequal baffle spacing inlet/ outlet
                if Re_s < 100:
                    n = 1 # Slope of friction factor curve
                else: n = 0.2
                R_s = (Lb_cut/LB_in)**(2-n) + (Lb_cut/LB_out)**(2-n) # correction factor

                # Shell side pressure drop (Excluding nozzles)
                dp_cent_baff = dp_shell_ideal*(N_b-1)*R_L*R_b # pressure drop in baffle windows
                dp_baff_window = dp_window # pressure drop in baffle windows
                dp_entrance_exit = dp_shell_ideal*R_s*R_b*(1+N_tcw/N_TCC) # pressure drop in entrance / exit baffles
                #print('dp for cent baffle '+str(dp_cent_baff))
                #print('dp for baffle wind '+str(dp_baff_window))
                #print('dp for entrance '+str(dp_entrance_exit))

                total_dp_shell = s3*(dp_cent_baff + dp_baff_window + dp_entrance_exit)
                total_dp_shell
                A_available = s3*L_eff * tn *np.pi*(Do/1000)
                ### Tube Side Heat transfer coeficient
                a_tube = (np.pi*((Di/1000)**2)*tn)/(4*pn) # Flow area
                v_t = (m_t/3600)/(rho_t*a_tube) #velocity through tube
                #print('value for v_t '+str(v_t))
                Re_t = (Di/1000)*rho_t*v_t/(mu_t/1000)
                #print('value for Re_t '+str(Re_t))
                Pr_t = Cp_t*(mu_t/1000)/k_t*3600
                #print('value for Pr_t '+str(Pr_t))
                L_eff = L-2*L_s/1000 # Effective tube length

                # Nusselt Number Calculation
                Nu_laminar = 1.86*(Re_t*Pr_t*(Di/1000)/L_eff)**(1/3)
                # Turbulent flow Petukhov-Kirillov
                f_turbulent = (1.58*np.log(Re_t)-3.28)**-2
                #print('value for f_turbulent '+str(f_turbulent))
                Nu_turb = (f_turbulent/2)*Re_t*Pr_t/(1.07+12.7*((f_turbulent/2)**0.5)*((Pr_t**(2/3))-1))
                # Transition flow Nu
                Nu_2300 = 1.86*(2300*Pr_t*(Di/1000)/L_eff)**(1/3)

                f_Re_10000 = (1.58*np.log(10000)-3.28)**-2
                Nu_10000 = (f_Re_10000/2)*10000*Pr_t/(1.07+12.7*((f_Re_10000/2)**0.5)*((Pr_t**(2/3))-1))

                Nu_trans = Nu_2300+(Nu_10000-Nu_2300)*(Re_t-2300)/(10000-2300)


                if Re_t <= 2300:
                    Nu_tube = Nu_laminar
                elif Re_t < 10000:
                    Nu_tube = Nu_trans
                else: Nu_tube = Nu_turb
                #print('value for Nu_tube '+str(Nu_tube))
                mu_t_w = mu_t
                h_t_i=Nu_tube*k_t/(Di/1000)*(mu_t/mu_t_w)**0.14
                #print('value for h_t_i '+str(h_t_i))
                total_dp_tube =s3*((4*f_turbulent*L*pn/(Di/1000))+4*pn)*rho_t*(v_t**2)/2/100000
                #print(total_dp_tube,L)
                dict_of_conductivity = {'Carbon Steel':38.69,'Copper':324.42,'Inconel':12.95,'Monel':21.28,'Nickel':52.09,'Stainless Steel':13.54}
                k_w_t = dict_of_conductivity['Carbon Steel']*1.163
                wall_resistance = (Do/2000)*np.log(Do/Di)/k_w_t
                U_clean = (1/((1/h_shell)+(Do/(h_t_i*Di))+wall_resistance))*1.163
                #print('dp for U_clean '+str(U_clean))
                U_dirty = 1/((1/U_clean)+(Do*fouling_t/(Di))+fouling_s)
                
                #print('dp for U_dirty '+str(U_dirty))
                
                Q, dTlm, ft = HB_data[0], HB_data[1], HB_data[2]
                A_available = s3*L_eff * tn *np.pi*(Do/1000)
                corrected_LMTD = dTlm*ft
                A_required = Q/(corrected_LMTD*U_dirty*1.163)
                A_ratio = A_required/A_available
                #print('value for A_required '+str(A_required))
                if A_ratio >1:
                    L = L *1500
            print('intialization data')
            print([total_dp_shell*10**8, h_shell , A_required,h_t_i,total_dp_tube*10**8,A_ratio,shell_D,A_available])
            return [total_dp_shell*10**8, h_shell , A_required,h_t_i,total_dp_tube*10**8,A_ratio,shell_D,A_available]
    shell_D = 387
    b_cut=25
    err_s=2
    
    #dps =91271301.16605562
    dict_of_conductivity = {'Carbon Steel':38.69,'Copper':324.42,'Inconel':12.95,'Monel':21.28,'Nickel':52.09,'Stainless Steel':13.54}
    k_w_t = dict_of_conductivity['Carbon Steel']*1.163
    wall_resistance = (Do/2000)*np.log(Do/Di)/k_w_t
    iteration = 0
    shell_d_iter = 0
    shell_no_iteration = 0
    #while shell_no_iteration < 8:
    while abs(err_s)>error_selected and pn <= 8 and iteration <= 100 and b_cut < 50:
        while abs(error) > 0.1 and shell_d_iter <= 50 :
                sol1 = initialise(b_cut,shell_D,pn,6000,t_p_angle,Do,pitch_ratio,pitch,s3)
                #sol2 = initialise(b_cut,shell_D,pn,3000,t_p_angle,Do,pitch_ratio)
                print(sol1) #,sol2)
                dps1_hs1 = sol1[0]/(sol1[1]**4.412)
                #dps2_hs2 = sol2[0]/(sol2[1]**4.412)
                #A1,A2 = sol1[7],sol2[7]
                A1 = sol1[7]
                #k2 = dps1_hs1-(k1*A1)
                #print('k1 and k2'+str(k1)+' '+str(k2))
                #print('balance shell'+str((k1*A1+k2)*(sol1[1]**2))+' '+str(sol1[0]))
                #print((k1*A2+k2)*(sol2[1]**2),sol2[0])
                dpt1_ht1 = sol1[4]/(sol1[3]**3.5)
                #dpt2_ht2 = sol2[4]/(sol2[3]**3.5)
                #A1,A2 = sol1[7],sol2[7]
                #k1 = (dps2_hs2-dps1_hs1)/(A2-A1)
                #k3 = (dpt2_ht2-dpt1_ht1)/(A2-A1)
                #k4 = dpt1_ht1-(k3*A1)
                
                #k5 = sol1[0]-((k3* A1+k4)* np.sign(w[1]) * (np.abs(w[1]) ** 3)*(np.sqrt(w[1])))
                #print('k3 k4 balance '+str((k3* sol2[7]+k4)* np.sign(sol2[3]) * (np.abs(sol2[3]) ** 3)*(np.sqrt(sol2[3]))  - sol2[4]) )
                #print('k3 and k4'+str(k3)+','+str(k4))
                
                Q, dTlm, ft = HB_data[0]/1.163, HB_data[1], HB_data[2]
                print(HB_data)
                c1 =Q/(dTlm*ft)
                c2=Q/(dTlm*ft)*(fouling_s+(fouling_t*((Do)/Di)))

                c3 = c1*((Do)/Di)
                c4 =wall_resistance
                
                # for a given variable w, this function returns F(w)
                # if w is the solution of the nonlinear system, then 
                # F(w)=0
                # F can be interpreted as the residual
                #print('c1 and c2='+str(c1)+','+str(c2))
                #print('c3 and c4='+str(c3)+','+str(c4))
                #def nonlinearEquation(w):
                    # A = w[2], h_s is w[0] and h_t is w[1]
                #   F=np.zeros(3)
                #  F[0]=(((c1/w[0])+(c3/w[1])+c2 + (c1*c4)) -w[2])
                # F[1]=((k3* w[2]+k4)* np.sign(w[1]) * (np.abs(w[1]) ** 3.5)  - dp_t) 
                    #F[2]=(k1*w[2]*(w[0]**4.412)+k2*(w[0]**4.412)-dp_s)
                    #F[1]=((k3* w[2])* np.sign(w[1]) * (np.abs(w[1]) ** 3.5)  - dp_t) 
                    #F[2]=(k1*w[2]*(w[0]**4.412)-dp_s)
                    #return F
                # generate an initial guess
                #initialGuess=np.array([sol2[1],sol2[3],sol2[7]])   
                
                # solve the problem    
                #solutionInfo=fsolve(nonlinearEquation,initialGuess,maxfev = 10000,full_output=1)
                
                import sympy as sy
                k1 = (dps1_hs1)/(A1)
                k3 = (dpt1_ht1)/(A1)
                print('kshell and ktube='+str(k1)+','+str(k3))
                x, y,z = sy.symbols("x y z", real=True, positive = True)
                # solve the problem    
                #solutionInfo=fsolve(nonlinearEquation,initialGuess,maxfev = 10000,full_output=1)
                exp1=sy.nsimplify((((c1/x)+(c3/y)+c2 + (c1*c4)) -(z)), rational=1)
                #exp2=sy.nsimplify((((k3* z+k4)* (y ** 3.5))  - dp_t), rational=1)
                #exp3=sy.nsimplify((k1*z*(x**4.412)+k2*(x**4.412)-dp_s), rational=1)
                exp2=sy.nsimplify((((k3* z)* (y ** 3.5))  - dp_t), rational=1)
                exp3=sy.nsimplify((k1*z*(x**4.412)-dp_s), rational=1)
                solutionInfo = sy.nsolve((exp1,exp2,exp3),(x,y,z),(sol1[1],sol1[3],sol1[7]),verify=False)
                print('area solved is ' + str(solutionInfo))
                #print(nonlinearEquation(solutionInfo))
                #print((nonlinearEquation(solutionInfo)))
                h_t_i = solutionInfo[1] #[1]
                Di = Do-2*thick
                Pr_t = Cp_t*(mu_t/1000)/k_t*3600
                mu_t_w = mu_t
                Nu_turb = h_t_i/(k_t/(Di/1000)*(mu_t/mu_t_w)**0.14)
                def nonlinearEquation2(w):
                    # A = w[2], h_s is w[0] and h_t is w[1]
                    F=np.zeros(2)
                
                    F[0]=np.exp((np.sqrt(1/w[0])+3.28)/1.58) - w[1]
                    F[1]= (w[0]/2)*w[1]*Pr_t/(1.07+12.7*((w[0]/2)**0.5)*((Pr_t**(2/3))-1))-Nu_turb
                    
                    return F
                # generate an initial guess
                initialGuess=np.array([0.0046,86385])   
                
                # solve the problem    
                solutionInfo2=fsolve(nonlinearEquation2,initialGuess,full_output=1)


                #f_turbulent = (1.58*np.log(Re_t)-3.28)**-2
                #Nu_turb = (solutionInfo[0][0]/2)*solutionInfo[0][1]*Pr_t/(1.07+12.7*((solutionInfo[0][0]/2)**0.5)*((Pr_t**(2/3))-1))
                print(solutionInfo2)
                itert_ht ,err_ht= 0,2
                Re_t = solutionInfo2[0][1] 
                itert_ht = 0
                #st.write(s3)
                while itert_ht <10 and abs(err_ht) > 0.1:
                    
                    v_t = Re_t/((Di/1000)*rho_t/(mu_t/1000))
                    #print(' tubes v_t is'+str(v_t))
                    #   a_tube = (np.pi*((Di/1000)**2)*tn)/(4*pn)
                    a_tube = (m_t/3600)/(v_t*rho_t)
                    tn = int(a_tube*4*pn/(np.pi*((Di/1000)**2)))
                    #print('Number of tubes is'+str(tn))
                    L_s = 0.1*shell_D
                    A_available = solutionInfo[2] #[2]
                    L_eff = A_available/(s3*tn *np.pi*(Do/1000))
                    # Nusselt Number Calculation
                    Nu_laminar = 1.86*(Re_t*Pr_t*(Di/1000)/L_eff)**(1/3)
                    # Turbulent flow Petukhov-Kirillov
                    f_turbulent = (1.58*np.log(int(Re_t))-3.28)**-2
                    #print('value for f_turbulent '+str(f_turbulent))
                    Nu_turb = (f_turbulent/2)*Re_t*Pr_t/(1.07+12.7*((f_turbulent/2)**0.5)*((Pr_t**(2/3))-1))
                    # Transition flow Nu
                    Nu_2300 = 1.86*(2300*Pr_t*(Di/1000)/L_eff)**(1/3)

                    f_Re_10000 = (1.58*np.log(10000)-3.28)**-2
                    Nu_10000 = (f_Re_10000/2)*10000*Pr_t/(1.07+12.7*((f_Re_10000/2)**0.5)*((Pr_t**(2/3))-1))

                    Nu_trans = Nu_2300+(Nu_10000-Nu_2300)*(Re_t-2300)/(10000-2300)


                    if Re_t <= 2300:
                        Nu_tube = Nu_laminar
                    elif Re_t < 10000:
                        Nu_tube = Nu_trans
                    else: Nu_tube = Nu_turb
                    #print('value for Nu_turb '+str(Nu_turb))
                   # print('value for Nu_tube '+str(Nu_tube))
                    h_t_i_calc=Nu_tube*k_t/(Di/1000)*(mu_t/mu_t_w)**0.14
                   # print('value for h_t_i_calc '+str(h_t_i_calc))
                    #print('value for Di '+str(Di)+' value for k_t '+str(k_t))
                    err_ht = (h_t_i-h_t_i_calc)/h_t_i
                    if err_ht > 0.1:
                        Re_t = Re_t*(1.1)
                    if err_ht <-0.1:
                        Re_t= Re_t*0.9
                    #st.write(h_t_i_calc,h_t_i,err_ht)
                    h_t_i=h_t_i_calc
                    itert_ht+=1
                
                
                print('Re ={} , ht ={} , vt ={} , at= {},Leff ={} , Nt ={}'.format(Re_t,h_t_i,v_t,a_tube,L_eff,tn))
                L_s = 0.1*shell_D
                A_available = solutionInfo[2] #[2]
                L_eff = A_available/(s3*tn *np.pi*(Do/1000))
                
                #st.write('shell no is  '+ str(s3))
                L = L_eff+2*L_s/1000
                #print('L required is '+ str(L))
                t_p = pitch_ratio*Do
                
                if t_p_angle == 45:
                    t_p_effective = t_p/(2**0.5)
                else: t_p_effective = t_p
                D_otl = shell_D - (12.5+(shell_D/200))

                theta_Ds = 2*np.arccos(1-(2*b_cut)/100) # Baffle window angle
                theta_CTL = 2*np.arccos(shell_D*(1-2*b_cut/100)/(D_otl-Do)) #b_cut baffle cut, shell_D isnide shell diameter mm
                F_w = (theta_CTL-np.sin(theta_CTL))/(2*np.pi)
                S_wg = (((shell_D/1000)**2)/8)*(theta_Ds-np.sin(theta_Ds)) #Gross window area
                S_wt = tn*F_w*np.pi*((Do/1000)**2)/4 #Window area occupied with tubes
                S_w = S_wg - S_wt # Net Cross flow area
                S_m = S_w
                G_w = (m_s/3600)/np.sqrt(S_m*S_w)
                v_w = G_w/rho_s

                #S_m = (Lb_cut/1000)*((shell_D-D_otl)+(D_otl-Do)*(t_p-Do)/t_p_effective)/1000
                Lb_cut = LB_in = LB_out = ((S_m*1000)/((shell_D-D_otl)+(D_otl-Do)*(t_p-Do)/t_p_effective))*1000
                #print(Lb_cut)
                N_b = 1 +int((L-(2*L_s*0.001)-(LB_in+LB_out)*0.001)/(Lb_cut*0.001)) # number of baffles
                #print(N_b,L_s)
                new_shell_D = ht.hx.DBundle_for_Ntubes_Phadkeb(tn, Do/1000, t_p/1000, pn, angle=30)*1000+3.1+0.004*shell_D
                error = (shell_D/new_shell_D)-1
                #print('shell diameter was'+str(shell_D))
                if error > 0.1:
                    shell_D=shell_D*0.9
                if error < -0.1:
                    shell_D=shell_D*1.1
                #print('shell diameter error is'+str(error))
                #print('new value for shell  '+str(shell_D)+' and tn is '+str(tn))
                #print('Area required is '+ str(A_available))
                shell_d_iter += 1
        print('new value for shell  '+str(shell_D)+' and tn is '+str(tn))
        mu_s_w = mu_s
        D_sb = 3.1+0.004*shell_D
        t_p = tpitch #p_ratio *Do
        L_tb = 0.4 # Diametral Tube-Baffle Clearance
        #t_p_effective 
        
        if t_p_angle == 45:
            t_p_effective = t_p/(2**0.5)
        else: t_p_effective = t_p
        D_otl = shell_D - (12.5+(shell_D/200))

        theta_Ds = 2*np.arccos(1-(2*b_cut)/100) # Baffle window angle
        theta_CTL = 2*np.arccos(shell_D*(1-2*b_cut/100)/(D_otl-Do)) #b_cut baffle cut, shell_D isnide shell diameter mm
        F_w = (theta_CTL-np.sin(theta_CTL))/(2*np.pi)
        S_wg = (((shell_D/1000)**2)/8)*(theta_Ds-np.sin(theta_Ds)) #Gross window area
        S_wt = tn*F_w*np.pi*((Do/1000)**2)/4 #Window area occupied with tubes
        S_w = S_wg - S_wt # Net Cross flow area
        S_m = S_w
        G_w = (m_s/3600)/np.sqrt(S_m*S_w)
        v_w = G_w/rho_s
        theta_Ds = 2*np.arccos(1-(2*b_cut)/100) # Baffle window angle
        S_sb = (shell_D/1000)*(D_sb/1000)*(np.pi-0.5*theta_Ds)    # Shell to baffle leakage area
        S_tb = (np.pi/4)*((((Do/1000)+(L_tb/1000))**2)-((Do/1000)**2))*tn*(1-F_w)  # Tube to baffle leakage area

        r_L = (S_sb+S_tb)/S_m # ratio of leakage to cross flow
        r_S = S_sb/(S_sb+S_tb) # ratio of shell-baffle to total area
        Lb_cut = LB_in = LB_out = ((S_m*1000)/((shell_D-D_otl)+(D_otl-Do)*(t_p-Do)/t_p_effective))*1000
        #print('dp for Lb_cut '+str(Lb_cut))
        #print('dp for b_cut '+str(b_cut))
        Gs = (m_s/3600)/S_m
        Re_s = (Do/1000)*Gs/(mu_s*0.001)
        #print('dp for Re_s '+str(Re_s))
        # 3. correction factor for bundle bypass
        # tube arrangement
        if t_p_angle == 30:
            t_arrg = np.sqrt(3)/2
        elif t_p_angle == 45:
            t_arrg = 1/np.sqrt(2)
        else : t_arrg = 1
        P_p = t_p * t_arrg # Tube row distance in flow direction
        N_TCC = (shell_D/P_p)*(1-(2*b_cut)/100) #N_TCC number of tube rows between baffle
        N_ss = int(N_TCC/6)
        r_ss = N_ss/N_TCC # Nss Number of sealing strips
        S_b = (Lb_cut/1000)*(shell_D-D_otl-(Do/2))/1000 #bundle pybass area
        if Re_s < 100:
            C_j = 1.35
        else: C_j = 1.25 # Correlation constant

        if r_ss >= 0.5:
            j_b = 1
        else: j_b = np.exp(-1*C_j*(S_b/S_m)*(1-((2*r_ss)**(1/3))))

        # 4. Correction factor for adverse temperature gradient
        N_tcw = (0.8/P_p)*(shell_D*b_cut/100-(shell_D-(D_otl-Do))/2)
        N_b = 1 +int((L-(2*L_s*0.001)-(LB_in+LB_out)*0.001)/(Lb_cut*0.001)) # number of baffles
        #print('value for N_b '+str(N_b))
        
        if Re_s > 100000:
                Re_temp = Re_s 
                Re_s = 99999
        ### Shell side pressure drop
        b1 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{1}'],0).sum()
        b2 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{2}'],0).sum()
        b3 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{3}'],0).sum()
        b4 = np.where((Re_s < j_const['Reynolds_max']) & (Re_s > j_const['Reynolds_min']) & (j_const['Layout'] == t_p_angle) ,j_const['b_{4}'],0).sum()
        if Re_s == 99999:
            Re_s = Re_temp
        b = b3/(1+0.14*(Re_s**b4))
        f_s = b1*((1.33/(t_p/Do))**b)*(Re_s**b2)
        dp_shell_ideal = 2*f_s*(Gs**2)*N_TCC*(mu_s_w/mu_s)**0.14/(rho_s)/100000 #N_TCC number of tube rows between baffle
        #print(dp_shell_ideal)

        # 1. correction factor for baffle leakage
        f_p = 0.8-0.15*(1+r_S) # r_s ratio of shell-baffle to total area
        R_L = np.exp(-1.33*(1+r_S)*(r_L**f_p))# r_l ratio of  leakage area to cross flow

        # pressure drop for an ideal window section
        S_wg = (((shell_D/1000)**2)/8)*(theta_Ds-np.sin(theta_Ds)) #Gross window area
        S_wt = tn*F_w*np.pi*((Do/1000)**2)/4 #Window area occupied with tubes
        S_w = S_wg - S_wt # Net Cross flow area
        G_w = (m_s/3600)/np.sqrt(S_m*S_w)
        v_w = G_w/rho_s
        #print('dp for G_w '+str(G_w))
        D_w = 4*S_w/(np.pi*(Do/1000)*tn*F_w+(shell_D/1000)*theta_Ds)

        #pressure drop for turbluent flow in ideal window section
        if Re_s >= 100:
            dp_window = N_b*R_L*(2+0.6*N_tcw)*(G_w**2)/(2*rho_s)/100000

            
        else:
            dp_window = N_b*R_L*((26*G_w*Cp_s/rho_s)*(N_tcw/(tpitch-Do)/1000 + (Lb_cut/1000)/(D_w**2)) + (G_w**2)/rho_s)/100000
            

        # 2. Correction factor for bundle bypass effect
        if Re_s < 100:
            C_r = 4.5
        else: C_r = 3.7

        if r_ss >= 0.5:
            R_b = 1
        else:
            R_b = np.exp(-1*C_r*(S_b/S_m)*(1-(2*r_ss)**(1/3)))

        # 3. Correction for unequal baffle spacing inlet/ outlet
        if Re_s < 100:
            n = 1 # Slope of friction factor curve
        else: n = 0.2
        R_s = (Lb_cut/LB_in)**(2-n) + (Lb_cut/LB_out)**(2-n) # correction factor

        # Shell side pressure drop (Excluding nozzles)
        dp_cent_baff = dp_shell_ideal*(N_b-1)*R_L*R_b # pressure drop in baffle windows
        dp_baff_window = dp_window # pressure drop in baffle windows
        dp_entrance_exit = dp_shell_ideal*R_s*R_b*(1+N_tcw/N_TCC) # pressure drop in entrance / exit baffles
        #print('dp for cent baffle '+str(dp_cent_baff))
        #print('dp for baffle wind '+str(dp_baff_window))
        #print('dp for entrance '+str(dp_entrance_exit))

        total_dp_shell = s3*(dp_cent_baff + dp_baff_window + dp_entrance_exit)
        #print('length '+str(L*1000))
        #print('total_dp_shell '+str(total_dp_shell))
        f_b_cut = b_cut
        err_s=total_dp_shell-(dp_s/10**8)
        
        #print('b_cut was'+str(b_cut) )
        if abs(err_s) > 0.1:
            if err_s >0.1:
                b_cut+=(1)#+err_s)*b_cut
            elif err_s <=-0.1:
                b_cut-=(1)#+err_s)*b_cut
            if (f_b_cut <= 10  ) and abs(err_s) >0.1: 
                pn +=2
                b_cut =25
                
            if b_cut == 49 and abs(err_s) >0.1:
                b_cut =25
                shell_D = 590
                error = 2
            if b_cut == 49 or (pn ==10 and abs(err_s) > 0.1):
                st.write('Couldnt fully converge')
                break
        else:
            if err_s >0.01:
                b_cut+=(0.01)#+err_s)*b_cut
            elif err_s <=-0.01:
                b_cut-=(0.1)#+err_s)*b_cut
            if (f_b_cut <= 10  ) and abs(err_s) >0.1: 
                pn +=2
                b_cut =25
                
            if b_cut == 49 and abs(err_s) >0.1:
                b_cut =25
                shell_D = 590
                error = 2
            if b_cut == 49 or (pn ==10 and abs(err_s) > 0.1):
                st.write('Couldnt fully converge')
                break
        
        #print(pn)
        iteration += 1
        #print(n)
        #print('shell dp error is'+str(err_s) )  
        
        #if iteration == 99:
           # s3 +=1
           # shell_no_iteration +=1
       # elif shell_no_iteration == 9:
          #  st.write('Couldnt fully converge')
          #  break
        #else: break

    
    print('final b_cut is'+str(b_cut) )
    print('final dps is'+str(total_dp_shell) )
    # 'Number of tubes','Number of passes','Do','Di','pitch type','Tube pitch','Length','Baffle Spacing','baffle cut','Shell D'
    #pitch = pitch
    geo_input_df = pd.DataFrame(index=geo_input_list)
    geo_input_df.loc[['Shell D','Baffle Spacing','Number of baffles','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut'],'Kern_summary']=geo_input_df.loc[['Shell D','Baffle Spacing','Number of baffles','Do','Di','Length','Number of tubes','Number of passes','Tube pitch','pitch type','baffle cut'],'Bell_summary'] =  [shell_D,Lb_cut,int(N_b),Do,Di,L,int(tn),pn,tpitch,pitch,f_b_cut]
    geo_list =  [int(tn),pn,Do,Di,pitch,tpitch,L*1000,Lb_cut,f_b_cut,shell_D/1000]
    return geo_list, geo_input_df,s3