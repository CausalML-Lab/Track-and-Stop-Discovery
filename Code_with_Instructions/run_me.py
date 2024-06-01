import results_generator

Max_samples = 100000 #(set upto 100,000)
Gap = 2000 #(Gives the resolution for the SHD plot versus samples)

# degree represents the how denese the DAGs are: density = 0 -> empty graph and  density = 1 -> complete graph

# To reproduce fig 3(a) of the paper
Num_Grphs= 50
nodes = 5
degree = 1 
results_generator.Run_CausalDiscovery_Algorithms_and_plot_results(Num_Grphs,nodes,degree,Max_samples,Gap)



# To reproduce fig 3(b) of the paper
Num_Grphs= 50
nodes = 6
degree = 1 
results_generator.Run_CausalDiscovery_Algorithms_and_plot_results(Num_Grphs,nodes,degree,Max_samples,Gap)



# To reproduce fig 3(c) of the paper
Num_Grphs= 50
nodes = 7
degree = 1 
results_generator.Run_CausalDiscovery_Algorithms_and_plot_results(Num_Grphs,nodes,degree,Max_samples,Gap)





# To reproduce fig 4(a) of the paper
Num_Grphs= 50
nodes = 10
degree = 0.1 
results_generator.Run_CausalDiscovery_Algorithms_and_plot_results(Num_Grphs,nodes,degree,Max_samples,Gap)



# To reproduce fig 4(b) of the paper
Num_Grphs= 50
nodes = 10
degree = 0.15
results_generator.Run_CausalDiscovery_Algorithms_and_plot_results(Num_Grphs,nodes,degree,Max_samples,Gap)



# To reproduce fig 4(c) of the paper
Num_Grphs= 50
nodes = 10
degree = 0.20
results_generator.Run_CausalDiscovery_Algorithms_and_plot_results(Num_Grphs,nodes,degree,Max_samples,Gap)
















