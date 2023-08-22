import numpy
import h5py

h5file = h5py.File("west.h5", 'r')
fi = 502 
li = 600
arr=numpy.zeros((2544,4))

count=0
for i in range(fi,li+1):

  i = str(i)
  iteration = "iter_" + str(numpy.char.zfill(i,8))
  pc0 = h5file['iterations'][iteration]['pcoord'][:,-1,0]
  pc1 = h5file['iterations'][iteration]['pcoord'][:,-1,1]
  pc2 = h5file['iterations'][iteration]['pcoord'][:,-1,2]

  weights = h5file['iterations'][iteration]['seg_index']['weight']
  #  for val in pc:
  for val in pc2:
    #    for val2 in pc2:
    if val > 3.0:
      #        if val2 > 10.0:
      #print(val,val2)
      nw = numpy.where(pc2==val) #and numpy.where(pc2==val2)
      for seg in nw:
        #    seg_weight = weights[seg]
        pcoord0 = pc0[seg]
        pcoord1 = pc1[seg]
        pcoord2 = pc2[seg]
        count=count+1        

  j=int(i)-502
  arr[j,0]=seg[0]
  arr[j,1]=pcoord0[0]
  arr[j,2]=pcoord1[0]
  arr[j,3]=pcoord2[0]

print(count)
print(arr.shape)
print(arr[arr[:, 1].argsort()])
