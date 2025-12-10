import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

import sys
import os
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from ldvm_back_up import ldvm
from matplotlib.animation import FuncAnimation

import yaml
from matplotlib.animation import PillowWriter



class PostProcessing:
    def __init__(self, path_data,config_data):
        self.path_data = path_data
        try:
            self.data=np.load(self.path_data, allow_pickle=True).item()
        except Exception as e:
            print(f"Error loading data from {self.path_data}: {e}")
            
        print(config_data)
       
        with open(config_data, 'r') as file:
            try:
                self.config = yaml.safe_load(file)
            except yaml.YAMLError as e:
                print(f"Error loading config data from {config_data}: {e}")
        self.ldvm_instance = ldvm(self.config)

    
    def make_animmation(self,k,alpha0,h0,phi,data_save):
        self.ldvm_instance.make_parameterized_motions(k=k,alpha0=alpha0,h0=h0,phi=phi,save=False)


        #ldvm_instance.load_motion()
        print("Making animation...")
        print("k:",k," alpha0:",alpha0," h0:",h0," phi:",phi)
        print(self.ldvm_instance.n_period,self.ldvm_instance.ppp)
        self.ldvm_instance.initialize_computation()
        ani=self.ldvm_instance.make_ldvm_animation(n_frame=self.ldvm_instance.n_period*self.ldvm_instance.ppp-1, add_reference=False,colorscale=True,fixed_airfoil=True)
    
        # Move the generated animation to the data path
        destination_path = os.path.join(data_save, f'ldvm_animation_k_{k}_alpha0_{alpha0*180/np.pi}_h0_{h0}_Phi_{phi}.gif')
        #ani.save(destination_path, writer='ffmpeg', fps=20)
        ani.save(destination_path, writer='pillow', fps=20)
        # animation_path = '../ldvm_animation.gif'
        # if os.path.exists(animation_path):
            
        #     os.rename(animation_path, destination_path)
        #     print(f"Animation moved to {destination_path}")
        # else:
        #     print("Animation file not found.")

    def plot_colored_pareto(self, data_save,x_min=-1, y_min=-2, x_max=1, y_max=2.7):
        pop_size=self.config["pop_size"]
        n_gens=self.config["n_generations"]
        print(len(self.data["All"]),type(self.data["All"]))
        all_f = -np.array(self.data["All"]).reshape(-1,pop_size, 2)
        print(all_f.shape)
        #all_f = -all_f[~np.any(all_f == 1e6, axis=2)]
        print(all_f.shape)
        
        fig, ax = plt.subplots(figsize=(4, 3), tight_layout=True)
        n_gens, pop_size, _ = all_f.shape
        c = np.tile(np.arange(1, n_gens + 1), pop_size)
        scatter = ax.scatter(all_f[:,:, 0], all_f[:,:, 1], c=c, cmap='Blues', marker='o', label='Pareto Front',s=10)
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Generation number')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel(r'$\eta$')
        ax.set_ylabel(r'$C_T$')
        plt.savefig(f"{data_save}/colored_pareto_front.png")
    
    def plot_lift(self, k,alpha0,h0,phi,data_save,hysteresis):
        self.ldvm_instance.initialize_computation()
        self.ldvm_instance.make_parameterized_motions(k=k, alpha0=alpha0, h0=h0, phi=phi, save=False)
    
        cl_history = [0]
        cd_history = [0]
        cm_history = [0]
        cn_history = [0]
        cnc_history = [0]
        cnnc_history = [0]
        cs_history = [0]
        non_l_history = [0]
        for i in range(self.ldvm_instance.n_period*self.ldvm_instance.ppp - 1):
            cl, cd, cm, lesp, re_le, cn = self.ldvm_instance.step()
            cl_history.append(cl)
            cd_history.append(cd)
            cm_history.append(cm)
            cn_history.append(cn)
            # cnc_history.append(cnc)
            # cnnc_history.append(cnnc)
            # cs_history.append(cs)
            # non_l_history.append(non_l)
        cd_history = np.array(cd_history)
    
        cl_history = np.array(cl_history)
        cm_history = np.array(cm_history)
        cn_history = np.array(cn_history)
        cnc_history = np.array(cnc_history)
        cnnc_history = np.array(cnnc_history)
        cs_history = np.array(cs_history)
        non_l_history = np.array(non_l_history)
        eta,ct,cp = self.ldvm_instance.compute_thrust_efficiency(cl_history, cd_history, cm_history,k,self.ldvm_instance.ppp)
        print("Thrust contribution", ct)
        print("efficiency",eta)
        
        fig,ax= plt.subplots(3, 2, figsize=(8, 6), tight_layout=True)
        ax[0,0].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cl_history[self.ldvm_instance.ppp:],color='red', label='Lift Coefficient $C_L$')
        ax[0,0].set_xlabel('Time/Period')
        ax[0,0].set_ylabel('$C_L$')
        
        ax[0,1].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cd_history[self.ldvm_instance.ppp:],color='red', label='Drag Coefficient $C_D$')
        ax[0,1].set_xlabel('Time/Period')
        ax[0,1].set_ylabel('$C_D$')
        ax[1,0].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cm_history[self.ldvm_instance.ppp:],color='red', label='Moment Coefficient $C_M$')
        ax[1,0].set_xlabel('Time/Period')
        ax[1,0].set_ylabel('$C_M$')
        
        
        ax[1,1].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cn_history[self.ldvm_instance.ppp:], color='red',label='Normal Force Coefficient $C_N$')
        
        #ax[1,1].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cnc_history[self.ldvm_instance.ppp:], color='blue',label='Normal Force Coefficient $C_{N_c}$')
        #ax[1,1].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cnnc_history[self.ldvm_instance.ppp:], color='green',label='Normal Force Coefficient $C_{NN_c}$')
        #ax[1,1].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,non_l_history[self.ldvm_instance.ppp:], color='black',label='Non Lift Force Coefficient $C_{NL}$')
        ax[1,1].set_xlabel('Time/Period')
        ax[1,1].set_ylabel('$C_N$')
        
        #ax[2,0].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cs_history[self.ldvm_instance.ppp:],color='orange')
        ax[2,0].set_xlabel('Time/Period')
        ax[2,0].set_ylabel('$C_S$')
        
        ax[2,1].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cl_history[self.ldvm_instance.ppp:],color='purple')
        ax[2,1].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cd_history[self.ldvm_instance.ppp:],color='blue')
        ax[2,1].plot(self.ldvm_instance.time[self.ldvm_instance.ppp:]/self.ldvm_instance.period,cn_history[self.ldvm_instance.ppp:],color='green')
        ax[2,1].set_xlabel('Time/Period')
        ax[2,1].set_ylabel('Force Coefficients')

        if hysteresis:
            fig,ax= plt.subplots(1, 2, figsize=(6, 4), tight_layout=True)
            alpha=self.ldvm_instance.alpha
            hdot=self.ldvm_instance.hdot

            Uinf=self.ldvm_instance.u_ref

            alpha_eff=alpha-np.arctan(hdot/Uinf)
            print(alpha_eff.shape,cl_history.shape)
            ax[0].plot(alpha_eff[self.ldvm_instance.ppp:],cl_history[self.ldvm_instance.ppp:],color='red', label='Lift Coefficient $C_L$')
            ax[1].plot(alpha_eff[self.ldvm_instance.ppp:],cd_history[self.ldvm_instance.ppp:],color='red', label='Lift Coefficient $C_d$')
            ax[0].set_xlabel('alpha eff')
            ax[1].set_xlabel('alpha eff')
            ax[0].set_ylabel('$C_L$')
            ax[1].set_ylabel('$C_D$')
            fig.savefig(f"{data_save}/lift_drag_coefficients_hysteresis_k_{k}_alpha0_{alpha0*180/np.pi}_h0_{h0}_Phi_{phi}.png")


  

        
        
        fig.savefig(f"{data_save}/lift_drag_moment_coefficients_k_{k}_alpha0_{alpha0*180/np.pi}_h0_{h0}_Phi_{phi}.png")
    
    
    def plot_LESP_alpha_eff(self, k,alpha0,h0,phi,data_save):
        self.ldvm_instance.initialize_computation()
        self.ldvm_instance.make_parameterized_motions(k=k, alpha0=alpha0, h0=h0, phi=phi, save=False)

        alpha=self.ldvm_instance.alpha
        hdot=self.ldvm_instance.hdot

        Uinf=self.ldvm_instance.u_ref

        alpha_eff=alpha-np.arctan(hdot/Uinf)
    
        lesp_history = []
        # for i in range(self.ldvm_instance.n_period*self.ldvm_instance.ppp - 1):
        #     cl, cd, cm, lesp, re_le,cn = self.ldvm_instance.step()

        #     lesp_history.append(lesp)
        # lesp_history = np.array(lesp_history)

        fig,ax= plt.subplots(figsize=(6, 4), tight_layout=True)
        ax.plot(self.ldvm_instance.time/self.ldvm_instance.period,hdot,color='green', label='hdot')
        ax.plot(self.ldvm_instance.time/self.ldvm_instance.period,alpha,color='blue', label='alpha')
        ax.plot(self.ldvm_instance.time/self.ldvm_instance.period,alpha_eff,color='red', label='alpha_reff')
        ax.set_xlabel('Time/Period')
        ax.set_ylabel('LESP')
        plt.title(f'Alpha eff for k={k}, alpha0={alpha0*180/np.pi}, h0={h0}, phi={phi}')
        fig.savefig(f"{data_save}/Alpha_eff_k_{k}_alpha0_{alpha0*180/np.pi}_h0_{h0}_Phi_{phi}.png")
        # fig,ax= plt.subplots(figsize=(6, 4), tight_layout=True)
        # ax.plot(self.ldvm_instance.time[self.ldvm_instance.ppp+1:]/self.ldvm_instance.period,lesp_history[self.ldvm_instance.ppp:],color='red', label='Leading Edge Suction Parameter (LESP)')
        # ax.set_xlabel('Time/Period')
        # ax.set_ylabel('LESP')
        # plt.title(f'Leading Edge Suction Parameter for k={k}, alpha0={alpha0*180/np.pi}, h0={h0}, phi={phi}')
        # fig.savefig(f"{data_save}/LESP_k_{k}_alpha0_{alpha0*180/np.pi}_h0_{h0}_Phi_{phi}.png")
    
    
    def plot_i_generation_pareto(self,i_gen,data_save):
        if i_gen < 0 or i_gen >= len(self.data['All']):
            print(f"Generation {i_gen} is out of bounds. Valid range: 0 to {len(self.data['All'])-1}")
            return
        
        all_f = np.array(self.data["All"][i_gen]).reshape(-1, 2)
        all_f = -all_f[~np.any(all_f == 1e6, axis=1)]
        
        
        
        fig, ax = plt.subplots(figsize=(4, 3), tight_layout=True)
        ax.scatter(all_f[:, 0], all_f[:, 1], c='blue', marker='o', label='Pareto Front')
        ax.set_xlabel(r'$\eta$')
        ax.set_ylabel(r'$C_T$')
        plt.title(f'Pareto front for generation {i_gen}')
        fig.savefig(f"{data_save}/pareto_front_generation_{i_gen}.png")
        plt.show()
    def plot_i_generation_design(self,i_gen,data_save):
        if i_gen < 0 or i_gen >= len(self.data['All']):
            print(f"Generation {i_gen} is out of bounds. Valid range: 0 to {len(self.data['All'])-1}")
            return
        
        all_f = np.array(self.data["All"][i_gen]).reshape(-1, 2)
        mask= np.any(all_f == 1e6, axis=1)
        all_f = -all_f[~mask]
        designs = np.array(self.data["All_designs"][i_gen]).reshape(-1, 4)
        designs = designs[~mask]
        
        omega = 2 * designs[:, 0] 
        
    
        f= omega/2/ np.pi # Frequency in Hz
        St=f *designs[:, 2]*2/1 # Assuming velocity is 1  
        
        fig, ax = plt.subplots(3, 2, figsize=(12, 9), tight_layout=True)
        ax[0, 0].scatter(St, all_f[:, 0], c='red', marker='x', label=r'$\eta$')
        ax[0, 0].set_xlabel(r'$S_t$')
        ax[0, 0].set_ylabel(r'$\eta$')
        ax_ = ax[0, 0].twinx()
        ax_.scatter(St, all_f[:, 1], c='blue', marker='o', label=r'$C_T$')
        ax_.set_ylabel(r'$C_T$')
        ax[0, 0].legend()
        #ax_.legend(loc='upper left')

        ax[0, 1].scatter(designs[:, 1]*180/np.pi, all_f[:, 0], c='red', marker='x', label=r'$\eta$')
        ax[0, 1].set_xlabel(r'$\alpha_0$')
        ax[0, 1].set_ylabel(r'$\eta$')
        ax_ = ax[0, 1].twinx()
        ax_.scatter(designs[:, 1]*180/np.pi, all_f[:, 1], c='blue', marker='o', label=r'$C_T$')
        ax_.legend()
        ax_.set_ylabel(r'$C_T$')
        ax[1, 0].scatter(designs[:, 2], all_f[:, 0], c='red', marker='x', label=r'$\eta$')
        ax[1, 0].set_xlabel(r'$h_0/c$')
        ax[1, 0].set_ylabel(r'$\eta$')
        ax_ = ax[1, 0].twinx()
        ax_.scatter(designs[:, 2], all_f[:, 1], c='blue', marker='o', label=r'$C_T$')
        ax_.set_ylabel(r'$C_T$')
        ax[1, 1].scatter(designs[:, 3], all_f[:, 0], c='red', marker='x', label=r'$\eta$')
        ax[1, 1].set_xlabel(r'$\Psi$')
        ax[1, 1].set_ylabel(r'$\eta$')
        ax_ = ax[1, 1].twinx()
        ax_.scatter(designs[:, 3], all_f[:, 1], c='blue', marker='o', label=r'$C_T$')
        ax_.set_ylabel(r'$C_T$')
        
        
        ax[2, 0].scatter(designs[:,0], all_f[:, 0], c='red', marker='x', label=r'$\eta$')
        ax[2, 0].set_xlabel(r'$k$')
        ax[2, 0].set_ylabel(r'$\eta$')
        ax_ = ax[2, 0].twinx()
        ax_.scatter(designs[:,0], all_f[:, 1], c='blue', marker='o', label=r'$C_T$')
        ax_.set_ylabel(r'$C_T$')
        ax[2, 0].legend()
        fig.savefig(f"{data_save}/designs_generation_{i_gen}.png")
        plt.show()
folder='/scratch/disc/b.martin/Documents/energy_harvesting/LDVM/NSGAB/paretto_results/Results_NSGA2_2025_07_28_10_29_33'
folder='/scratch/disc/b.martin/Documents/energy_harvesting/LDVM/NSGAB/paretto_results/Results_NSGA2_2025_07_28_11_39_22'
folder='/scratch/disc/b.martin/Documents/energy_harvesting/LDVM/NSGAB/paretto_results/Results_NSGA2_2025_07_28_16_52_30'
folder='/scratch/disc/b.martin/Documents/energy_harvesting/LDVM/NSGAB/paretto_results/Results_NSGA2_2025_07_31_11_44_24'
folder='/scratch/disc/b.martin/Documents/energy_harvesting/LDVM/NSGAB/paretto_results/Results_NSGA2_2025_07_31_11_44_44'
folder='/scratch/disc/b.martin/Documents/energy_harvesting/LDVM/NSGAB/paretto_results/Results_NSGA2_2025_09_01_14_12_26'

data_save=os.path.join(folder,'callback_data.npy')
config= os.path.join(folder,'config.yaml')
pp= PostProcessing(data_save,config)
#pp.make_animmation(0.99993335,0.70839055,0.99949244,1.87209392,folder)
#pp.make_animmation(0.3231061 , 0.46522308, 0.99998693, 1.36552392,folder)
#pp.make_animmation(0.80997988, 0.82344625, 0.99998003, 1.45602329,folder)
# pp.make_animmation(0.3055602 , 0.74079869, 1.99735096, 1.42931189,folder)


# pp.plot_colored_pareto(folder,x_min=-1, y_min=-5, x_max=1, y_max=12)
# pp.plot_i_generation_pareto(399,folder)
# pp.plot_i_generation_design(399,folder)

data_plot=pp.data['All_designs'][399]


for i in range(len(data_plot[0])):
   
    print(i,type(i))
    
    
    print(data_plot[0][i,:])
    
    print('wesh')
    pp.plot_lift(data_plot[0][i,0],data_plot[0][i,1],data_plot[0][i,2],data_plot[0][i,3],folder,hysteresis=True)
    #pp.plot_LESP_alpha_eff(data_plot[0][i,0],data_plot[0][i,1],data_plot[0][i,2],data_plot[0][i,3],folder)
    #pp.make_animmation(data_plot[0][i,0],data_plot[0][i,1],data_plot[0][i,2],data_plot[0][i,3],folder)
dd
# print(data_plot)
# print(pp.data['All'][399])
# pp.make_animmation(0.99693632, 0.59986571, 0.99957501, 1.68807532,folder)
pp.make_animmation(0.37602386, 0.54671483, 0.88980057, 1.29050389,folder)
ss
pp.plot_lift(0.51089754, 0.62599875, 0.99964206, 1.39849274,folder)


# for i in [24]:
#     pp.make_animmation(data_plot[i][0],data_plot[i][1],data_plot[i][2],data_plot[i][3],folder)
#     print(data_plot[i][0],data_plot[i][1],data_plot[i][2],data_plot[i][3])
# print(data_plot,len(data_plot[0]))
