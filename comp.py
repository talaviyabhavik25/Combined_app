import pandas as pd
import numpy as np
import streamlit as st

url_1 = 'http://raw.githubusercontent.com/Ahmedhassan676/compressor_evaluation/main/comp.csv'
df_comp_table = pd.read_csv(url_1)

def Z_calculations(df,t_suc,p_suc):
        pc = (df['mol%']*df['Pc']).sum() * 0.01
        tc = (df['mol%']*(df['Tc']+460)).sum() * 0.01  
        Tr = (t_suc*1.8 + 32+460)/tc
        Pr = (p_suc +1.03323)*14.2233/pc
        A = [1,0.31506237,-1.04670990,-0.57832729,0.53530771,-0.61232032,-0.10488813,0.68157001,0.68446549]
        rho = (0.27*Pr)/Tr
        Z = 1
        error = 10
        y = 1/0.84
        while error > 0.001:
            part_1 = (A[1]+(A[2]/Tr)+(A[3]/(Tr**3)))*rho*y
            part_2 = (A[4]+(A[5]/Tr))*(rho**2)*(y**2)
            part_3 = ((A[5]*A[6]*(rho**5)*(y**5))/Tr)
            part_4 = (A[7]*(rho**3)*(y**3))/((Tr**3)*(1+(A[8]*(rho**2)*(y**2)))*np.exp(-A[8]*(rho**2)*(y**2)))
            Z = A[0]+part_1+part_2+ part_3 + part_4 
            error = abs(Z-(1/y))
            y = 1/Z 
        return Z
def calculations(Q_nm3, suc_p,suc_t, disch_p,disch_t,m_wt,z,k):
        m_kg_hr = Q_nm3*m_wt*0.0446098580359918482953956685434
        suc_p = suc_p + 1.03323
        disch_p = disch_p + 1.03323
        suc_t = suc_t + 273.15
        disch_t = disch_t + 273.15
        k_1_k = (k-1)/k
        k_k_1 = k/(k-1)
        r = 8.3143
        td_ts = disch_t/suc_t
        dd_ds = disch_p/suc_p
        p_port = (dd_ds**k_1_k) - 1
        r_mwt = r / m_wt
        poly_coef = (1-(np.log(td_ts)/np.log(dd_ds)))**-1
        n_1_n = (poly_coef - 1)/ poly_coef
        poly_eff = (k_1_k*np.log(dd_ds)/np.log(td_ts)) * 100
        td_adiab = suc_t*(((disch_p/suc_p)**k_1_k)-1)+suc_t-273.15
        adiab_eff = (suc_t*(((disch_p/suc_p)**k_1_k)-1)/(disch_t-suc_t)) * 100
        power_kw = ((1/(adiab_eff*36)))*r_mwt*suc_t*z*m_kg_hr*(p_port)*k_k_1
        return poly_eff, adiab_eff, td_adiab, power_kw,m_kg_hr,dd_ds,poly_coef,td_ts

def k_calculations(df,df_comp_table,Q_nm3,suc_p,suc_t, disch_p,disch_t):
        
        temperatures = np.array([suc_t,disch_t])*1.8+ 32 
        df['y_MCp_suc']=[np.interp(temperatures[0],df_comp_table['Gas'][3:],df_comp_table['{}'.format(compound)][3:]) for compound in df.index]
        suc_MCp = (df['mol%']*df['y_MCp_suc']).sum()*0.01
        k_suc = suc_MCp/(suc_MCp-1.986)
        df['y_MCp_disch']=[np.interp(temperatures[1],df_comp_table['Gas'][3:],df_comp_table['{}'.format(compound)][3:]) for compound in df.index]
        disch_MCp = (df['mol%']*df['y_MCp_disch']).sum()*0.01
        k_disch = disch_MCp/(disch_MCp-1.986)
        k = (k_suc + k_disch)/2

        m_wt = np.sum(df['mol%']*df['m.wt'])*0.01
        df1= pd.DataFrame({'mol%':100,'m.wt':m_wt,'cp/cv':k}, index=['Total'])
        df = df[df['mol%'] != 0].sort_values(by='mol%', ascending=False).append(df1)
        
        p = [suc_p,disch_p]
        t = [suc_t,disch_t]
        
        
        p = (np.array(p)+1) * 14.2233
        t = np.array(t)*1.8 + 491.67
        z1 = Z_calculations(df,suc_t,suc_p)
        z2 = Z_calculations(df,disch_t,disch_p)
        z = (z1+z2)*0.5
        return df, z, m_wt, k
def choose_composition():
            
            df = pd.DataFrame({'Composition/property':df_comp_table.columns[1:],'mol%':np.zeros(len(df_comp_table.columns)-1), 'm.wt':df_comp_table.iloc[0,1:],'Pc':df_comp_table.iloc[1,1:],'Tc':df_comp_table.iloc[2,1:]})
            try:
                sum_of_comp = 0 
                c1,c2,c3,c15,c4,c5,c6,c7,c8,c9,c16,c10,c11,c13,c14,nh3,h2o = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
                options_list = [df_comp_table.columns[i] for i in [22,1,4,3,6,11,10,13,15,20,21,24,25,23,18,17,19]]
                while sum_of_comp != 100:
                    options = st.multiselect(
                    'Select your components', options_list
                    )
                    if df_comp_table.columns[22] in options:
                        c1 = st.number_input('hydrogen%', key = 'c1')
                    if df_comp_table.columns[1] in options:
                        c2 = st.number_input('methane%', key = 'c2')
                    if df_comp_table.columns[4] in options:
                        c3 = st.number_input('ethane%', key = 'c3')
                    if df_comp_table.columns[3] in options:
                        c15 = st.number_input('ethylene%', key = 'c15')
                    if df_comp_table.columns[6] in options:
                        c4 = st.number_input('propane%', key = 'c4')
                    if df_comp_table.columns[11] in options:
                        c5 = st.number_input('nbutane%', key = 'c5')
                    if df_comp_table.columns[10] in options:
                        c6 = st.number_input('isobutane%', key = 'c6')
                    if df_comp_table.columns[13] in options:
                        c7 = st.number_input('pentane%', key = 'c7')
                    if df_comp_table.columns[15] in options:
                        c8 = st.number_input('hexane%', key = 'c8')
                    if df_comp_table.columns[20] in options:
                        c9 = st.number_input('Oxygen%', key = 'c9')
                    if df_comp_table.columns[21] in options:
                        c16 = st.number_input('nitrogen%', key = 'c16')
                    if df_comp_table.columns[24] in options:
                        c10 = st.number_input('carbon monoxide%', key = 'c10')
                    if df_comp_table.columns[25] in options:
                        c11 = st.number_input('carbon dioxide%', key = 'c11')
                    if df_comp_table.columns[23] in options:
                        c13 = st.number_input('hydrogen sulphide%', key = 'c13')
                    if df_comp_table.columns[18] in options:
                        c14 = st.number_input('air%', key = 'c14')
                    if df_comp_table.columns[17] in options:
                        nh3 = st.number_input('Ammonia%', key = 'nh3')
                    if df_comp_table.columns[19] in options:
                        h2o = st.number_input('Water vapor%', key = 'h2o')
                    if c1 or c2 or c3 or c15 or c4 or c5 or c6 or c7 or c8 or c9 or c16 or c10 or c11 or c13 or c14 or nh3 or h2o:
                        c = []
                        for i in (c1,c2,c3,c15,c4,c5,c6,c7,c8,c9,c16,c10,c11,c13,c14,nh3,h2o):
                            c.append(i)
                        
                        for (i,j) in zip(options_list,c):
                            if j != None:
                                    df.loc[i,'mol%'] = j
                        
                        sum_of_comp = np.sum(df['mol%'])
                        
                st.success('Composition in Mol. percent completed!', icon="✅")
                
                return df[df['mol%'] != 0]

            except (ValueError, st.errors.DuplicateWidgetID): pass
            except (TypeError, KeyError, ZeroDivisionError):st.write('Please Check your data')
def main():
    html_temp="""
        <div style="background-color:lightblue;padding:16px">
        <h2 style="color:black"; text-align:center> Compressor Performance Evaluation </h2>
        </div>
        <style>
        table {
        font-family: arial, sans-serif;
        border-collapse: collapse;
        width: 100%;
        }

        td, th {
        border: 1px solid #dddddd;
        text-align: left;
        padding: 8px;
        }
        </style>
            """
    st.markdown(html_temp, unsafe_allow_html=True)
    stage_no = st.selectbox('Number of compression stages?',(1,2,3,4), key = 'test_type')
    if stage_no == 1:
        
    
        Q_nm3= st.number_input('Capacity (nm3/hr)', key = 'cap_1')
        suc_p= st.number_input('Suction pressure (kg/cm2.g)', key = 'sucp')
        suc_t= st.number_input('Suction temperature (C)', key = 'suct')
        disch_p= st.number_input('Discharge pressure (kg/cm2.g)', key = 'disp')
        disch_t= st.number_input('Discharge temperature (C)', key = 'dist')
        s1 = st.selectbox('Estimate M.wt, Cp/Cv and Z?',('I already have these values','Yes'), key = 'k_calculations')
        if s1 == 'I already have these values':
            m_wt= st.number_input('Molecular weight' , key = 'mwt')
            z= st.number_input('Compressibility factor', key = 'z')
            k= st.number_input('Cp/Cv', key = 'k')
        else:
            try:
                df_comp = choose_composition()
                df_comp, z, m_wt, k = k_calculations(df_comp,df_comp_table,Q_nm3,suc_p,suc_t, disch_p,disch_t)
            except (ValueError,TypeError, KeyError, ZeroDivisionError):st.write('your total mol. percent should add up to 100')
            except AttributeError: pass
    
        if st.button("Reveal Calculations", key = 'calculations_table'):
                try:
                    poly_eff, adiab_eff, td_adiab, power_kw,m_kg_hr,dd_ds,poly_coef,td_ts =   calculations(Q_nm3, suc_p,suc_t, disch_p,disch_t,m_wt,z,k)
                    url = 'http://raw.githubusercontent.com/Ahmedhassan676/compressor_evaluation/main/compressor_table.csv'
                    df = pd.read_csv(url, index_col=0)
                    df['Calculations'] = [Q_nm3,m_kg_hr,suc_p, suc_t, disch_p,disch_t,td_adiab, power_kw,m_wt,z,k,dd_ds,td_ts,poly_coef,poly_eff, adiab_eff]
                    st.dataframe(df)
                except (ValueError, TypeError, KeyError, ZeroDivisionError, UnboundLocalError): st.write('Please Check your data')
    else: # Multipule stages
        url = 'http://raw.githubusercontent.com/Ahmedhassan676/compressor_evaluation/main/multi_stage_table.csv'
        df = pd.read_csv(url, index_col=0)
        for i in range(stage_no):
            df['Stage No. '+'{}'.format(i+1)] = 0.00
        edited_df = st.experimental_data_editor(df)
        Z_multi = list(range(stage_no))
        m_wt_multi =  list(range(stage_no))
        k_multi = list(range(stage_no))
        poly_eff, adiab_eff, td_adiab, power_kw,m_kg_hr,dd_ds,poly_coef,td_ts = list(range(stage_no)),list(range(stage_no)),list(range(stage_no)),list(range(stage_no)),list(range(stage_no)),list(range(stage_no)),list(range(stage_no)),list(range(stage_no))
        Q_nm3 = np.array(df.iloc[0,:])
        suc_p = np.array(df.iloc[1,:])
        suc_t = np.array(df.iloc[2,:])
        disch_p = np.array(df.iloc[3,:])
        disch_t = np.array(df.iloc[4,:])
        try:
            s2 = st.selectbox('Estimate M.wt, Cp/Cv and Z?',('I already have these values','Yes'), key = 'k_calculations_multi')
            if s2 == 'I already have these values':
                m_wt= st.number_input('Molecular weight' , key = 'mwt_multi')
                z= st.number_input('Compressibility factor', key = 'z_multi')
                k= st.number_input('Cp/Cv', key = 'k_multi')
                
                if st.button("Reveal Calculations", key = 'calculations_table124'):
                    edited_df = pd.DataFrame(edited_df)
                    url = 'http://raw.githubusercontent.com/Ahmedhassan676/compressor_evaluation/main/compressor_table.csv'
                    df = pd.read_csv(url, index_col=0)
                    for i in range(stage_no):
                        poly_eff[i], adiab_eff[i], td_adiab[i], power_kw[i],m_kg_hr[i],dd_ds[i],poly_coef[i],td_ts[i] =   calculations(edited_df.iloc[0,i],edited_df.iloc[1,i],edited_df.iloc[2,i], edited_df.iloc[3,i],edited_df.iloc[4,i],m_wt,z,k)
                        df['Stage No. '+'{}'.format(i+1)] = [edited_df.iloc[0,i],m_kg_hr[i],edited_df.iloc[1,i],edited_df.iloc[2,i], edited_df.iloc[3,i],edited_df.iloc[4,i],td_adiab[i], power_kw[i],m_wt,z,k,dd_ds[i],td_ts[i],poly_coef[i],poly_eff[i], adiab_eff[i]]
                    st.dataframe(df)
            else:
                df_comp = pd.DataFrame(choose_composition())
                if st.button("Reveal Calculations", key = 'calculations_table22'):
                    edited_df = pd.DataFrame(edited_df)
                    url = 'http://raw.githubusercontent.com/Ahmedhassan676/compressor_evaluation/main/compressor_table.csv'
                    df = pd.read_csv(url, index_col=0)
                    for i in range(stage_no):
                        df_comp_new, Z_i, M_wt_i,K_i = k_calculations(df_comp,df_comp_table,edited_df.iloc[0,i],edited_df.iloc[1,i],edited_df.iloc[2,i], edited_df.iloc[3,i],edited_df.iloc[4,i])
                        Z_multi[i] = Z_i
                        m_wt_multi[i] = M_wt_i
                        k_multi[i] = K_i
                        poly_eff[i], adiab_eff[i], td_adiab[i], power_kw[i],m_kg_hr[i],dd_ds[i],poly_coef[i],td_ts[i] =   calculations(edited_df.iloc[0,i],edited_df.iloc[1,i],edited_df.iloc[2,i], edited_df.iloc[3,i],edited_df.iloc[4,i],m_wt_multi[i],Z_multi[i],k_multi[i])
                        
                        df['Stage No. '+'{}'.format(i+1)] = [edited_df.iloc[0,i],m_kg_hr[i],edited_df.iloc[1,i],edited_df.iloc[2,i], edited_df.iloc[3,i],edited_df.iloc[4,i],td_adiab[i], power_kw[i],m_wt_multi[i],Z_multi[i],k_multi[i],dd_ds[i],td_ts[i],poly_coef[i],poly_eff[i], adiab_eff[i]]
                    st.dataframe(df)
                
        except (ValueError,TypeError, KeyError, UnboundLocalError):st.write('Please Check your input data!')
        except ZeroDivisionError: pass

if __name__ == '__main__':
    main()
