import numpy
import h5py

tau=0.1 ## tau value in ns



############## single mab #############

#f1 = h5py.File("../reweight_unbinding_multimab/west.h5", 'r')
#t1_400 = f1['summary']['n_particles'][:-1][:400]
#t1_500 = f1['summary']['n_particles'][:-1][:500]
#sum1_400=sum(t1_400)
#sum1_500=sum(t1_500)



#wallclock1_400 = f1['summary']['walltime'][:-1][:400]
#wallclock1_500 = f1['summary']['walltime'][:-1][:500]
#sum_wc_1_400=numpy.sum(wallclock1_400)
#sum_wc_1_500=numpy.sum(wallclock1_500)

############## multi mab #############

f2 = h5py.File("tmp.h5", 'r') 
t2 = f2['summary']['n_particles'][:-1][:62]
sum2=sum(t2)

wallclock2 = f2['summary']['walltime'][:-1][:62]
sum_wc_2=numpy.sum(wallclock2)



print("\n ### Comparing aggregate times ###")

#print(" Aggregate time (ns) for singlemab (until it62):", sum1_400*tau)
#print(" Aggregate time (ns) for singlemab (unti. it500):", sum1_500*tau)
print(" Aggregate time (ns) for multimab:", sum2*tau)



#print("\n ### Comparing wallclock times ###")


#print(" Wallclock time (days) for singlemab: (until it 400)", sum_wc_1_400*1.15E-5)
#print(" Wallclock time (days) for singlemab: (until it 500)", sum_wc_1_500*1.15E-5)
print(" Wallclock time (days) for multimab:", sum_wc_2*1.15E-5)


