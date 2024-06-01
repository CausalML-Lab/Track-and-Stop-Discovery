import math
import numpy as np

def plot_shd(Data_save,color_line,plt,T,plt_label):
        dd  = Data_save
        if (plt == 0):
                        import matplotlib.pyplot as plt
                       
        m  = np.mean(dd,axis = 0)
        sd = np.std(dd,axis =0)/(math.sqrt(dd.shape[0]))
        mup = m+sd
        mlp = m-sd
        m = m[0:T]
        mup = mup[0:T]
        
        mlp = mlp[0:T]
        
       
        color_area = color_line + (1-color_line)*2.3/4 ;
        
        plt.plot(range((T)),m,color= color_line,label = plt_label)
        plt.fill_between(range((T)),mup,mlp,color=color_area)
        plt.xlabel('Interventional Samples')
        plt.ylabel('SHD')
    
        plt.grid(True)
        plt.legend()
        return  plt
    
def Run_CausalDiscovery_Algorithms_and_plot_results(Num_Grphs,nodes,degree,Max_samples,Gap):

        import math
        import numpy as np
        
        print( 'Running Our Track and Stop Discovery Algorithm...' )    
        import TsP
        Data_save = TsP.Track_and_stop(Num_Grphs,nodes,degree,Max_samples)
        color_line = np.array([179, 63, 64])/255 ;
        plt_label= 'Track and Stop Algorithm'
        plt = plot_shd(np.array(Data_save),color_line,0,Max_samples-Gap,plt_label)

        print( 'Running Active Learning using DCT Algorithm...' )    
        import DCT
        Data_save = DCT.DCT_discovery(Num_Grphs,nodes,degree,Max_samples,Gap)
        color_line = np.array([72, 161, 77])/255 ;
        plt_label= 'Active Learning using DCT'
        plt = plot_shd(np.array(Data_save),color_line,plt,Max_samples-Gap,plt_label)

        print( 'Running Radomized Traget Discovery Algorithm...' )    
        import Rnd
        Data_save = Rnd.Random_Interventions(Num_Grphs,nodes,degree,Max_samples,Gap)
        color_line = np.array([1, 119, 179])/255 ;
        plt_label= 'Random'
        plt =plot_shd(np.array(Data_save),color_line,plt,Max_samples-Gap,plt_label)
        
        
        
        
       
        print( 'Running Adaptivity-sensitive search Algorithm...' )     
        import Radpt
        Data_save= Radpt.r_apative(Num_Grphs,nodes,degree,Max_samples,Gap)
        color_line = np.array([0.83,0.09,0.99]) ;
        plt_label= 'Adaptivity-sensitive search'
        plot_shd(np.array(Data_save),color_line,plt,Max_samples-Gap,plt_label)
        
        print( 'Running GIES Algorithm...' )     
        import Radpt
        import GIES
        Data_save= GIES.GIES_discovery(Num_Grphs,nodes,degree,Max_samples,5000)
        color_line = np.array([0.13,0.19,0.45]) ;
        plt_label= 'GIES'
        plt = plot_shd(np.array(Data_save),color_line,plt,Max_samples-Gap,plt_label)
        plt.show()



       