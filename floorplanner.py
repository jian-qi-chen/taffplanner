#!/usr/bin/env python3
""" author: Jianqi Chen """
import sys, svgwrite, random, copy, numpy
from SlicingTree import STree

resource_file = "./test.res"
module_file = "./test.module"
exec(open(resource_file).read())
exec(open(module_file).read())

s_floorplan = None # the slicing tree, just remember it's a global variable

colormap = {'LAB':'blue','RAM':'orange','DSP':'red'}
asp_max = 4 # maximum aspect ratio for a single module, minimum is 1/asp_max
alpha = 1.5 # for intermediate floorplans, maximum area is (alpha*WIDTH) * (alpha*HEIGHT), 1<alpha<=2
leaf_IRL = {} #IRL of leaves (functional modules), usually only generate once for given resource/module files
history_record = {} #stores slicing trees and their cost that are evaluated

def main(argv):
    FloorplaningSA()

# check whether the cell is inside the rectangle
# (cx,cy): coordinate of the cell, every cell's (x,y) actually represents 4 cells
# e.g. cx=1, cy=3, the four cells' location: (1,3),(WIDTH+1,3),(1,HEIGHT+3),(WIDTH+1,HEIGHT+3)
# cw: cell width, ch: cell height, (rx,ry):left bottom coordinate of rectangle, rw: rect width, rh: rect height
def CellInRect(cx,cy,cw,ch,rx,ry,rw,rh):
    #upper right coordinate
    cx1 = cx+cw-1
    cy1 = cy+ch-1
    rx1 = rx+rw-1
    ry1 = ry+rh-1
    
    count = 0
    if cx >= rx and cy >= ry and cx1 <= rx1 and cy1 <= ry1:
        count += 1
    if (cx+WIDTH) >= rx and cy >= ry and (cx1+WIDTH) <= rx1 and cy1 <=ry1:
        count += 1
    if cx >= rx and (cy+HEIGHT) >= ry and cx1 <= rx1 and (cy1+HEIGHT) <= ry1:
        count += 1
    if (cx+WIDTH) >= rx and (cy+HEIGHT) >= ry and (cx1+WIDTH) <= rx1 and (cy1+HEIGHT) <= ry1:
        count += 1
        
    return count

# irreducible realization list(IRL) generation for a single module, at a single location
# x,y is the coordinate of the left bottom corner, mod_res is a dict of resource required e.g. {'lAB':12,'RAM':1}
# return a list of rectangle r = (x,y,w,h), w=width, h=height
def SingleIRLGen(x,y,mod_res):
    v_WIDTH = int(alpha*WIDTH) #virtual width of the whole floorplanning area
    v_HEIGHT = int(alpha*HEIGHT)
    IRL = [] #irreducible realization list
    
    #start searching from a square
    s_width = 1
    
    while True:
        check_list = [False for i in range(len(mod_res))]
        i = 0
        for res in mod_res:
            cell_w = resource[res][0] #cell width
            cell_h = resource[res][1] #cell height
            
            count = 0
            for cell in res_cells:
                # check resource type
                if cell[0] == res:
                    count += CellInRect(cell[1],cell[2],cell_w,cell_h,x,y,s_width,s_width)
            
            if count >= mod_res[res]:
                check_list[i] = True
                i += 1
            else:
                break
            
        if False not in check_list:
            IRL.append( (x,y,s_width,s_width) )
            break
        elif (x + s_width) <= v_WIDTH and (y + s_width) <= v_HEIGHT:
            s_width += 1
        else:
            return IRL

    #searching for a long rectangle (width>height)
    l_height = s_width - 1
    l_width = s_width + 1
    
    while l_height >= 1:
        while l_width/l_height < asp_max:
            check_list = [False for i in range(len(mod_res))]
            i = 0
            for res in mod_res:
                cell_w = resource[res][0] #cell width
                cell_h = resource[res][1] #cell height
                
                count = 0
                for cell in res_cells:
                    # check resource type
                    if cell[0] == res:
                        count += CellInRect(cell[1],cell[2],cell_w,cell_h,x,y,l_width,l_height)
                
                if count >= mod_res[res]:
                    check_list[i] = True
                else:
                    break
                
            if False not in check_list:
                IRL.append( (x,y,l_width,l_height) )
                break
            else:
                l_width += 1
                
        l_height -= 1
        
    #searching for a tall rectangle (width<height)
    t_height = s_width + 1
    t_width = s_width - 1
    
    while t_width >= 1:
        while t_height/t_width < asp_max:
            check_list = [False for i in range(len(mod_res))]
            i = 0
            for res in mod_res:
                cell_w = resource[res][0] #cell width
                cell_h = resource[res][1] #cell height
                
                count = 0
                for cell in res_cells:
                    # check resource type
                    if cell[0] == res:
                        count += CellInRect(cell[1],cell[2],cell_w,cell_h,x,y,t_width,t_height)
                
                if count >= mod_res[res]:
                    check_list[i] = True
                else:
                    break
                
            if False not in check_list:
                IRL.append( (x,y,t_width,t_height) )
                break
            else:
                t_height += 1
                
        t_width -= 1
        
    return IRL

# IRL generation for a single module at all possible locations
# IRL belongs to {r=(x,y,w,h)|x>0,y>0,x+w<alpha*WIDTH,y+h<alpha*HEIGHT}
# mod_res is a dict of resource required e.g. {'lAB':12,'RAM':1}
# store the IRL in a list
def ModuleIRLGen(mod_res):
    IRL = []
    for x in range(int(alpha*WIDTH)):
        for y in range(int(alpha*HEIGHT)):
            IRL += SingleIRLGen(x,y,mod_res)
    
    return IRL

# get IRL for 'V' node(vertical slicing), given (x,y), IRL of left/right child
# returned list:  [(IRL of the v node, rectangles of each module which is a descendant of the v node)]
# e.g. [( (1,2,4,4),[('jpeg',(1,2,2,2)),('fir',(3,4,2,2))] ),(...)]
def VGetIRL(x,y,l_IRL, r_IRL):
    v_IRL = [] #return list
    
    ll_IRL = list(filter(lambda p: p[0][0]==x and p[0][1]==y, l_IRL))
    ll_IRL.sort(key=lambda p: p[0][2]) #sort according to width
    for i in range(len(ll_IRL)):
        l_w = ll_IRL[i][0][2]
        l_h = ll_IRL[i][0][3]
        l_mod_use = ll_IRL[i][1]
        if i == 0:
            l_h1 = alpha*HEIGHT
        else:
            l_h1 = ll_IRL[i-1][0][3]#upperbound for right child's height
        
        rr_IRL = list(filter(lambda p: p[0][0]==(x+l_w) and p[0][1]==y and p[0][3]<l_h1, r_IRL))
        rr_IRL.sort(key=lambda p: p[0][2])
        for j in range(len(rr_IRL)):
            r_w = rr_IRL[j][0][2]
            r_h = rr_IRL[j][0][3]
            r_mod_use = rr_IRL[i][1]
            
            w_new = l_w + r_w
            h_new = max(l_h,r_h)
            v_IRL.append( ((x,y,w_new,h_new), l_mod_use + r_mod_use) )
            if r_h <= l_h:
                break
            
    return v_IRL

# get IRL for 'H' node(horizontal slicing), given (x,y), left child is at bottom, right child top
def HGetIRL(x,y,l_IRL, r_IRL):
    h_IRL = [] #return list
    
    ll_IRL = list(filter(lambda p: p[0][0]==x and p[0][1]==y, l_IRL))
    ll_IRL.sort(key=lambda p: p[0][2]) #sort according to width
    for i in range(len(ll_IRL)):
        l_w = ll_IRL[i][0][2]
        l_h = ll_IRL[i][0][3]
        l_mod_use = ll_IRL[i][1]
        if i == 0:
            l_w1 = alpha*WIDTH
        else:
            l_w1 = ll_IRL[i+1][0][2]#upperbound for right child's width
        
        rr_IRL = list(filter(lambda p: p[0][0]==x and p[0][1]==(y+l_h) and p[0][2]<l_w1, r_IRL))
        rr_IRL.sort(key=lambda p: p[0][2])
        for j in reversed(range(len(rr_IRL))):
            r_w = rr_IRL[j][0][2]
            r_h = rr_IRL[j][0][3]
            r_mod_use = rr_IRL[i][1]
            
            w_new = max(l_w,r_w)
            h_new = l_h + r_h
            h_IRL.append( ((x,y,w_new,h_new), l_mod_use + r_mod_use) )
            if r_w <= l_w:
                break
            
    return h_IRL

# usually input the index of root node when call this function, which is always 0.
# Then get IRL for all nodes in the slicing tree
def EvaluateNode(index):
    global leaf_IRL
    # for leaves
    if s_floorplan.slicing_tree[index].l == None and s_floorplan.slicing_tree[index].r == None:
        mod_name = s_floorplan.slicing_tree[index].t
        mod_res = module_list[ mod_name ][0]
        
        if mod_name in leaf_IRL: #IRL already calculated
            IRL_tmp = leaf_IRL[mod_name]
        else: #leaf_IRL does not have the IRL
            IRL_tmp = ModuleIRLGen(mod_res)
            leaf_IRL[mod_name] = IRL_tmp

        s_floorplan.slicing_tree[index].IRL = [ (x,[(mod_name,x)]) for x in IRL_tmp]
        return
    
    EvaluateNode(s_floorplan.slicing_tree[index].l - 1)
    EvaluateNode(s_floorplan.slicing_tree[index].r - 1)
    
    IRL = []
    left_child_IRL = s_floorplan.slicing_tree[ s_floorplan.slicing_tree[index].l - 1 ].IRL
    right_child_IRL = s_floorplan.slicing_tree[ s_floorplan.slicing_tree[index].r - 1 ].IRL
    
    if s_floorplan.slicing_tree[index].t == 'V':
        for x in range(int(alpha*WIDTH)):
            for y in range(int(alpha*HEIGHT)):
                IRL += VGetIRL(x,y,left_child_IRL,right_child_IRL)
    elif s_floorplan.slicing_tree[index].t == 'H':
        for x in range(int(alpha*WIDTH)):
            for y in range(int(alpha*HEIGHT)):
                IRL += HGetIRL(x,y,left_child_IRL,right_child_IRL)
                
    s_floorplan.slicing_tree[index].IRL = IRL

# floorplanning algorithm based on Simulated Annealing
def FloorplaningSA():
    global s_floorplan
    global leaf_IRL
    global history_record
    # acceptance probability
    def accept_prob(old_cost, new_cost, T):
        return numpy.exp( (old_cost-new_cost)/T )
        
    # cost function
    def cost_func(area):
        cost = area/(WIDTH*HEIGHT)
        return cost
    
    # calculate IRL of root node
    def new_floorplan():
        global s_floorplan
        global history_record
        nonlocal best_cost
        nonlocal best_fp
        
        # if the slicing tree is evaluated before
        cur_polish = ' '.join( s_floorplan.Polish() )
        if cur_polish in history_record:
            return history_record[cur_polish]
        
        # slicing tree not evaluated before
        EvaluateNode(0)
        useful_floorplan = list(filter(lambda p: (p[0][0]+p[0][2])<WIDTH and (p[0][1]+p[0][3])<HEIGHT, s_floorplan.slicing_tree[0].IRL))
        root_area = [ a[0][2]*a[0][3] for a in s_floorplan.slicing_tree[0].IRL]
        
        if root_area:
            root_area_min = min(root_area)
            cost = cost_func(root_area_min)
        else:
            cost = 99999
            
        if useful_floorplan: #if there is floorplan fits the real maximum area
            real_root_area = [ a[0][2]*a[0][3] for a in useful_floorplan ]
            real_root_area_min = min(real_root_area)
            index_best = real_root_area.index(real_root_area_min)
            real_cost = cost_func(real_root_area_min)
            if real_cost < best_cost:
                best_cost = real_cost
                best_fp = useful_floorplan[index_best]
        
        history_record[cur_polish] = cost
                
        return cost

    mod_names = list(module_list)    
    # generate IRL for leaves (modules)
    for mod_n in mod_names:
        leaf_IRL[mod_n] = ModuleIRLGen(module_list[mod_n][0])
    
    s_floorplan = STree(mod_names) #ramdom slicing tree, nodes have (t,p,l,r), no (x,y,w,h)
    
    print(s_floorplan)#for debug
    polish_exp = ' '.join( s_floorplan.Polish() )
    print('0.Polish Expression: '+polish_exp)
    
 #   s_floorplan.M3()
    best_cost = 99999
    best_fp = ()  
 
    T = 1.0 # temperature
    T_min = 0.05
    coeff = 0.88
    
    print('Simulated Anealing started. Starting with T = '+str(T)+', will be stopped when T < '+str(T_min))
    
    old_cost = new_floorplan()
    
    while T > T_min:
        for i in range(25):
            s_floorplan.ClearIRL()
            tmp_tree = copy.deepcopy( s_floorplan.slicing_tree )
            
            select = random.randint(0,2)
            if select == 0:
                s_floorplan.M1()
            elif select == 1:
                s_floorplan.M2()
            else:
                if not s_floorplan.M3():
                    s_floorplan.M1()
                    
            new_cost = new_floorplan()
            ap = accept_prob(old_cost, new_cost, T)
            if ap > random.random():
                old_cost = new_cost
            else:
                s_floorplan.slicing_tree = copy.deepcopy( tmp_tree ) #recover the original tree
                
        T = T*coeff
        print('25 floorplan evaluated, T = '+str(T))
        print('So far, best cost found = '+str(best_cost)+' , best floorplan: '+str(best_fp))
        
    # draw result
    if best_cost < 0:
        print('No feasible floorplan found')
    else:
        print('best floorplan: '+str(best_fp))
        modules = best_fp[1]
        DrawFloorplan('floorplan.svg',modules)


def DrawFloorplan(svg_name, modules):
    global colormap
    cell_width = 10
    cell_spacing = 4
    label_spacing = 60
    label_spacing_to_floorplan = 30
    
    # draw resource labels
    dr = svgwrite.Drawing(svg_name,profile='tiny')
    i = 0
    for cur_r in list(resource):
        if cur_r not in colormap:
            colormap[cur_r] = 'svgwrite.rgb('+str(random.randint(128,255))+','+str(random.randint(128,255))+','+str(random.randint(128,255))+',"%")'
            
        dr.add(dr.rect(insert=(i*label_spacing+10, HEIGHT*(cell_width+cell_spacing)+label_spacing_to_floorplan-(resource[cur_r][1]-1)*cell_width),size=(resource[cur_r][0]*cell_width,resource[cur_r][1]*cell_width),stroke='gray',fill=colormap[cur_r]))
        dr.add(dr.text(cur_r,insert=(i*label_spacing+15+resource[cur_r][0]*cell_width, HEIGHT*(cell_width+cell_spacing)+label_spacing_to_floorplan+10)))
        
        i = i+1
    
    # draw empty floorplan
    dr.add(dr.rect(insert=(0,0),size=(cell_width*WIDTH+cell_spacing*(WIDTH+1),cell_width*HEIGHT+cell_spacing*(HEIGHT+1)),stroke='orange',fill='white'))
    for cell in res_cells:
        if cell[0] not in resource:
            print('Error: unknown resource in resource file')
            sys.exit(1)
            
        cwidth = resource[cell[0]][0] #width in number of width of a LAB
        cheight = resource[cell[0]][1]
        initial_x = cell[1]*cell_width + (cell[1]+1)*cell_spacing #coordinate of top right corner
        initial_y = (HEIGHT-cell[2]-cheight)*cell_width + (HEIGHT-cell[2]-cheight+1)*cell_spacing
        rwidth = cwidth*cell_width + (cwidth-1)*cell_spacing #width on the graph drawn (in pixel)
        rheight = cheight*cell_width + (cheight-1)*cell_spacing
        ccolor = colormap[cell[0]]
        
        dr.add(dr.rect(insert=(initial_x,initial_y),size=(rwidth,rheight),stroke='gray',fill=ccolor))
        
    # draw modules
    for mod in modules:
        mod_name = mod[0]
        mx = mod[1][0]
        my = mod[1][1]
        mw = mod[1][2]
        mh = mod[1][3]
        initial_x = mx*cell_width +(mx+0.5)*cell_spacing
        initial_y = (HEIGHT-my-mh)*cell_width + (HEIGHT-my-mh+0.5)*cell_spacing
        rwidth = mw*cell_width + mw*cell_spacing
        rheight = mh*cell_width + mh*cell_spacing
        
        dr.add(dr.rect(insert=(initial_x,initial_y),size=(rwidth,rheight),stroke='black',fill='white'))
        dr.add(dr.text(mod_name,insert=(initial_x,initial_y+20)))
    
    dr.save()
        

if __name__ == "__main__":
    main(sys.argv[1:])