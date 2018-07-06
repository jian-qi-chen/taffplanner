# resource file, use python syntax
# set resource types, e.g. LAB (resource name): [1,1](width, height)
resource = {'LAB':(1,1), 'RAM':(1,1), 'DSP':(1,2)}

# set FPGA width and height
WIDTH = 12
HEIGHT = 8

# set each resource block, e.g. ('LAB',0,0) - resource type is LAB, left bottom coordinate is (0,0)
res_cells = []
# LAB
for i in list(range(1,3))+list(range(4,6))+list(range(7,9))+list(range(10,12)):
    for j in range(8):
        res_cells.append( ('LAB',i,j) )
    
# RAM    
for i in [0,6]:
    for j in range(8):
        res_cells.append( ('RAM',i,j) )
        
# DSP
for i in [3,9]:
    for j in [x*2 for x in range(4)]:
        res_cells.append( ('DSP',i,j) )

