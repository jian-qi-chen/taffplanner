#!/usr/bin/env python3
""" author: Jianqi Chen """
import sys, os, svgwrite, random, copy, numpy, getopt, re, threading, multiprocessing, json, datetime, itertools
from queue import Queue
from SlicingTree import STree

s_floorplan = None # the slicing tree, just remember it's a global variable

PROCESS_MAX = 5 # Maximal number of processes for RunHotSpot() and ModuleIRLGen()
colormap = {'LAB':'blue','RAM':'orange','DSP':'red'}
asp_max = 3.5 # maximum aspect ratio for a single module, minimum is 1/asp_max
alpha = 1.3 # for intermediate floorplans, maximum area is (alpha*WIDTH) * (alpha*HEIGHT), 1<alpha<=2
temp_am = 310 #ambient temperature for HotSpot
leaf_IRL = {} #IRL of leaves (functional modules), usually only generate once for given resource/module files
history_record = {} #stores slicing trees and their cost that are evaluated
final_result = None #[max_temperature,[(mod_name1,[x,y,w,h]),(mod_name2,[x,y,w,h])...]]

# coefficients in SA's cost function: cost = ( aaa*temp_max + bbb*tot_ap + ccc*tot_wl ) * ( 0.2*area_ex/(WIDTH*HEIGHT) + 1)
aaa = 2#0.2
bbb = 0.1
ccc = 0.0005#0.002

# coefficients in MCG's cost functional: cost = mcg_alpha*BoundingArea + mcg_beta*sum{sqrt(mod_area_self*mod_area_other)/(dist*(power_density difference + mcg_const))}
mcg_alpha = 0.97
mcg_beta = 0.03
mcg_const = 0.05

def usage():
    print("To run the program for [design_name], [design_name].module and [design_name].res should be provided. Then enter:")
    print("\t./floorplanner.py [design_name]")
    print("design_name = 'test' by default")

def main():
    print('Starting at '+str( datetime.datetime.now() ))
    
    # Draw empty floorplan
#    DrawFloorplan('./output_files/empty_floorplan.svg',[])
    
    # set alpha value
    PreEstimate()
    
    # floorplanning using simulated annealing
#    mod_loc_list = FloorplaningSA()

    # floorplanning using cluster growth
    FloorplaningMCG()
    
    print('Finishing at '+str( datetime.datetime.now() ))


# handle command line arguments by setting global variables    
def SetGlobalVar(argv):
    global design_name, beta, gamma
    os.system('mkdir -p output_files')
    os.system('chmod 744 ./hotspot/hotspot ./hotspot/grid_thermal_map.pl')
    # Output file
    if not os.path.isdir('./output_files/json'):
        os.system('mkdir ./output_files/json')
        
    with open('./output_files/json/final_result','w') as json_file:
            json.dump(final_result, json_file) 
    
    try:
        opts, args = getopt.getopt(argv,'h',['help','aaa','bbb','ccc'])
    except getopt.GetopError:
        usage()
        sys.exit(1)
        
    for opt, arg in opts:
        if opt in ('-h','--help'):
            usage()
            sys.exit(0)
        elif opt in ('--aaa'):
            aaa = float(arg)
        elif opt in ('--bbb'):
            bbb = float(arg)
        elif opt in ('--ccc'):
            ccc = float(arg)
        else:
            usage()
            sys.exit(2)
    
    if len(args) > 1:
        usage()
        sys.exit(2)
    elif len(args) == 1:
        design_name = args[0]
        
def PreEstimate():
    global alpha
    cell_type = list( map(lambda p: p[0],res_cells) )
    usage_ratio = []
    
    for res in resource:
        count_a = cell_type.count(res) #number of available resources
        count_n = 0 #number of resources needed
        
        for mod in module_list:
            count_n += module_list[mod][0][res]
            
        ratio_u = count_n/count_a
        usage_ratio.append(ratio_u)
        
    key_ratio = max(usage_ratio)
    if key_ratio > 0.86:
        print('Not enough FPGA resources')
        sys.exit(0)
    elif key_ratio > 0.8:
        alpha = 1.5
    elif key_ratio > 0.5:
        alpha = 1.4
    elif key_ratio > 0.4:
        alpha = 1.25
    elif key_ratio > 0.3:
        alpha = 1.1
    else:
        alpha = 1
        
        
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
            res_cells_useful = filter(lambda p: p[0] == res and p[1]>=x and p[2]>=y, res_cells)
            for cell in res_cells_useful:
                count += CellInRect(cell[1],cell[2],cell_w,cell_h,x,y,s_width,s_width)
            
            if count >= mod_res[res]:
                check_list[i] = True
                i += 1
            else:
                break
  
        if False not in check_list:
            IRL.append( (x,y,s_width,s_width) )
            rm_flg = True
            break
        elif (x + s_width) <= v_WIDTH and (y + s_width) <= v_HEIGHT:
            s_width += 1
        else:
            break

    #searching for a long rectangle (width>height)
    l_height = s_width - 1
    l_width = s_width
    
    while l_height >= 1:
        while l_width/l_height < asp_max and (x + l_width) <= v_WIDTH:
            check_list = [False for i in range(len(mod_res))]
            i = 0
            for res in mod_res:
                cell_w = resource[res][0] #cell width
                cell_h = resource[res][1] #cell height
                
                count = 0
                res_cells_useful = filter(lambda p: p[0] == res and p[1]>=x and p[2]>=y, res_cells)
                for cell in res_cells_useful:
                    count += CellInRect(cell[1],cell[2],cell_w,cell_h,x,y,l_width,l_height)

                if count >= mod_res[res]:
                    check_list[i] = True
                    i += 1
                else:
                    break
   
            if False not in check_list:
                if rm_flg:
                    IRL.pop()
                IRL.append( (x,y,l_width,l_height) )
                rm_flg = True
                break
            else:
                l_width += 1
                rm_flg = False
        
        if (x + l_width) > v_WIDTH:
            break
                
        l_height -= 1
        
    #searching for a tall rectangle (width<height)
    t_height = s_width
    t_width = s_width - 1
    
    while t_width >= 1:
        while t_height/t_width < asp_max and (y + t_height) <= v_HEIGHT:
            check_list = [False for i in range(len(mod_res))]
            i = 0
            for res in mod_res:
                cell_w = resource[res][0] #cell width
                cell_h = resource[res][1] #cell height
                
                count = 0
                res_cells_useful = filter(lambda p: p[0] == res and p[1]>=x and p[2]>=y, res_cells)
                for cell in res_cells_useful:
                    count += CellInRect(cell[1],cell[2],cell_w,cell_h,x,y,t_width,t_height)
                
                if count >= mod_res[res]:
                    check_list[i] = True
                    i += 1
                else:
                    break
                
            if False not in check_list:
                if rm_flg:
                    IRL.pop()
                IRL.append( (x,y,t_width,t_height) )
                rm_flg = True
                break
            else:
                t_height += 1
                rm_flg = False
        
        if (y + t_height) > v_HEIGHT:
            break
        
        t_width -= 1
        
    return IRL

# IRL generation for a single module at all possible locations
# IRL belongs to {r=(x,y,w,h)|x>0,y>0,x+w<alpha*WIDTH,y+h<alpha*HEIGHT}
# mod_res is a dict of resource required e.g. {'lAB':12,'RAM':1}
# store the IRL in a list
def ModuleIRLGen(mod_res, mod_name):
    sema_irlgen.acquire()
    IRL = []

    for x in range(int(alpha*WIDTH)):
        for y in range(int(alpha*HEIGHT)):
            IRL += SingleIRLGen(x,y,mod_res)
    
    if not os.path.isdir('./output_files/json'):
        os.system('mkdir ./output_files/json')
             
    irl_dict = {mod_name: IRL}
    with open('./output_files/json/'+mod_name+'_IRL','w') as json_file:
        json.dump(irl_dict, json_file) 
        
    sema_irlgen.release()
    return IRL

# Generate all IRL for leaves of the slicing tree, using multiprocessing
def AllLeavesIRLGen():
    global leaf_IRL
    mod_names = list(module_list) 
    os.system('mkdir -p ./output_files/json')
    
    print('Generating irreducible realization list for all modules')
    process_list = []
    for mod_n in mod_names:
        t = multiprocessing.Process( target = ModuleIRLGen, args = (module_list[mod_n][0], mod_n) )
        process_list.append(t)
        
    for process in process_list:
        process.start()
                
    for process in process_list:
        process.join()
        
    for mod_n in mod_names:
        with open('./output_files/json/'+mod_n+'_IRL','r') as f:
            new_IRL = json.load(f)
        leaf_IRL.update(new_IRL)

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
            r_mod_use = rr_IRL[j][1]
            
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
        if i == len(ll_IRL)-1:
            l_w1 = alpha*WIDTH
        else:
            l_w1 = ll_IRL[i+1][0][2]#upperbound for right child's width
        
        rr_IRL = list(filter(lambda p: p[0][0]==x and p[0][1]==(y+l_h) and p[0][2]<l_w1, r_IRL))
        rr_IRL.sort(key=lambda p: p[0][2])
        for j in reversed(range(len(rr_IRL))):
            r_w = rr_IRL[j][0][2]
            r_h = rr_IRL[j][0][3]
            r_mod_use = rr_IRL[j][1]
            
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
            leaf_IRL[mod_name] = ModuleIRLGen(mod_res, mod_name)
            IRL_tmp = leaf_IRL[mod_name]
            

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

# generate floorplan file for thermal simulator HotSpot
# module_loc_list = [(mod_name1,(x,y,w,h)),(mod_name2..)..], folder is the name of ./hotspot/folder, can be empty string
def FLPgen(module_loc_list,folder):
    # how long is 1 unit width in IRLs
    cell_width = 0.0001
    cell_height = 0.0001
    
    flp_file = open('./hotspot/'+folder+'/'+design_name+'.flp','w')
    flp_file.write('root\t'+str( WIDTH*cell_width )+'\t'+str( HEIGHT*cell_height )+'\t0\t0\n')
    
    for mod in module_loc_list:
        width = str(mod[1][2]*cell_width)
        height = str(mod[1][3]*cell_height)
        left_x = str(mod[1][0]*cell_width)
        bottom_y = str(mod[1][1]*cell_height)
        flp_file.write(mod[0]+'\t'+width+'\t'+height+'\t'+left_x+'\t'+bottom_y+'\n')
        
    flp_file.close()

# generate power trace file for thermal simulator HotSpot
# mod_list should have the same format of the module_list in file test.module
def PTRACEgen(mod_list, root_power,folder):
    ptrace_file = open('./hotspot/'+folder+'/'+design_name+'.ptrace','w')
    
    name = 'root' #first line
    power = str(root_power/1000) #second line
    for mod in mod_list:
        name += '\t'
        name += mod
        power += '\t'
        power += str( mod_list[mod][1]/1000 ) #unit: Watt
        
    ptrace_file.write(name+'\n'+power+'\n')
    ptrace_file.close()
    
# generate flp and ptrace file from module_list and IRL. root_power is static power of the whole floorplanning area(unit: mW)
# run HotSpot and get thermal map file [design_name].grid.steady as output
def RunHotSpot(module_loc_list, root_power, folder):
    global sema_hotspot
    sema_hotspot.acquire()
    FLPgen(module_loc_list, folder)
    PTRACEgen(module_list, root_power, folder)
    
    os.system('cd hotspot && ./hotspot -c hotspot.config -f ./'+folder+'/'+design_name+'.flp -p ./'+folder+'/'+design_name+'.ptrace -steady_file ./'+folder+'/'+design_name+'.steady -model_type grid -grid_steady_file ./'+folder+'/'+design_name+'.grid.steady > /dev/null 2>%1')
    sema_hotspot.release()

# read [design_name].grid.steady and return the highest temperature(unit: Kelvin)    
def ReadTempMax(folder):
    thermal_map = open('./hotspot/'+folder+'/'+design_name+'.grid.steady','r')
    temperature_str = re.findall(r'\d{3}\.\d{2}', thermal_map.read() )
    thermal_map.close()
    temperature = map(float, temperature_str)
    max_temp = max(temperature)
    
    return max_temp
    
# given the location of every module[(mod_name,(x,y,w,h)),(..).], return total interconnection wire length = sum( Manhattan distance between block centers)
def TotalWireLen(module_loc_list):
    mod_loc_dict = {}
    # read the list into a dict, and only store the center location
    for mod_loc in module_loc_list:
        mod_loc_dict[ mod_loc[0] ] = (mod_loc[1][0]+mod_loc[1][2]/2, mod_loc[1][2]+mod_loc[1][3]/2)
        
    wl_sum = 0
    for mod in module_list: #module_list is the global variable in .module file_name
        center1 = mod_loc_dict[mod]
        for other_mod in module_list[mod][2]:
            center2 = mod_loc_dict[other_mod]
            wire_num = module_list[mod][2][other_mod]
            wl_sum += wire_num * ( abs(center1[0]-center2[0]) + abs(center1[1]-center2[1]) )
            
    wl_sum = wl_sum/2 #actually calculated twice above
    return wl_sum
    
# given the location of every module[(mod_name,(x,y,w,h)),(..).], return aspect ratio sum = sum( longer side length/shorter side length of all blocks), this a metric for internal wire length
def RatioSum(module_loc_list):
    r_sum = 0
    for mod_loc in module_loc_list:
        ap_ratio = max( [ mod_loc[1][2]/mod_loc[1][3], mod_loc[1][3]/mod_loc[1][2] ] )
        r_sum += ap_ratio
        
    return r_sum

# floorplanning algorithm based on Simulated Annealing
def FloorplaningSA():
    global s_floorplan
    global leaf_IRL
    global history_record
    global final_result
    # acceptance probability
    def accept_prob(old_cost, new_cost, T):
        return numpy.exp( (old_cost-new_cost)/T )
        
    # cost function, area_ex - area exceeds WIDTH*HEIGHT, temp_max - max temperature(in Kelvin), tot_wl - total external interconnection wire length, tot_ap - aspect ratio sum of all modules
    def cost_func(area_ex, temp_max, tot_ap, tot_wl):
        cost = cost = ( aaa*temp_max + bbb*tot_ap + ccc*tot_wl ) * ( 0.2*area_ex/(WIDTH*HEIGHT) + 1)
        return cost
    
    # calculate IRL of root node
    def new_floorplan():
        global s_floorplan
        global history_record
        nonlocal max_tp
        nonlocal best_cost
        nonlocal best_fp
        nonlocal best_temp_max
        
        # if the slicing tree is evaluated before
        cur_polish = ' '.join( s_floorplan.Polish() )
        if cur_polish in history_record:
            return history_record[cur_polish]
        
        # slicing tree not evaluated before
        EvaluateNode(0)           
        
        useful_floorplan = list(map(lambda p: (p[0][0]+p[0][2])<=WIDTH and (p[0][1]+p[0][3])<=HEIGHT, s_floorplan.slicing_tree[0].IRL))
        
        if s_floorplan.slicing_tree[0].IRL:
            ex_root_area = [] # area exceeds WIDTH*HEIGHT 
            thread_list = [] # mutithreading for RunHotSpot
            temp_max_list = [] # maximal temperature from HotSpot
            total_ratio = [] # aspect ratio rum of all modules
            total_wirelen = [] # total length of external interconnection wire between modules
            cost_list = []
            for IR in s_floorplan.slicing_tree[0].IRL:
                if (IR[0][0]+IR[0][2] > WIDTH) and (IR[0][1]+IR[0][3] <= HEIGHT):
                    area_exceed = IR[0][3] * (IR[0][0]+IR[0][2]-WIDTH)
                elif (IR[0][0]+IR[0][2] <= WIDTH) and (IR[0][1]+IR[0][3] > HEIGHT):
                    area_exceed = IR[0][2] * (IR[0][1]+IR[0][3]-HEIGHT)
                elif (IR[0][0]+IR[0][2] > WIDTH) and (IR[0][1]+IR[0][3] > HEIGHT):
                    area_exceed = IR[0][2]*IR[0][3] - (WIDTH-IR[0][0])*(HEIGHT-IR[0][1])
                else:
                    area_exceed = 0
                    
                ex_root_area.append(area_exceed)
                total_wirelen.append( TotalWireLen(IR[1]) )
                total_ratio.append( RatioSum(IR[1]) )
                    
                cur_index = s_floorplan.slicing_tree[0].IRL.index(IR)
                os.system('mkdir -p ./hotspot/'+str(cur_index))
                if area_exceed == 0:
                    t = threading.Thread( target = RunHotSpot, args = (IR[1],0,str(cur_index)) )
                    thread_list.append(t)
              
            for thread in thread_list:
                thread.start()
                
            for thread in thread_list:
                thread.join()
            
            for i in range(len(useful_floorplan)):
                if useful_floorplan[i] == True:
                    t_max = ReadTempMax(str(i))
                    if t_max > max_tp:
                        max_tp = t_max
                    temp_max_list.append(t_max)
                else:
                    if max_tp == 0:
                        temp_max_list.append(315)
                    else:
                        temp_max_list.append(max_tp)
            
            for i in range(len(temp_max_list)):
                cur_max_temp = temp_max_list[i]
                area_exceed = ex_root_area[i]
                ratio_sum = total_ratio[i]
                wire_length = total_wirelen[i]
                cost_list.append( cost_func(area_exceed, cur_max_temp, ratio_sum, wire_length) )
                    
            cost = min(cost_list)
        else:
            cost = 999
            
        if 1 in useful_floorplan: #if there is floorplan fits the real maximum area
            for i, useful in enumerate(useful_floorplan):
                if useful == True:
                    real_cost = cost_list[i]
                    if real_cost < best_cost:
                        best_cost = real_cost
                        best_temp_max = temp_max_list[i]
                        best_fp = s_floorplan.slicing_tree[0].IRL[i]

        history_record[cur_polish] = cost
                
        return cost


    mod_names = list(module_list)    
    # generate IRL for leaves (modules)
    AllLeavesIRLGen()
        
    s_floorplan = STree(mod_names) #ramdom slicing tree, nodes have (t,p,l,r), no (x,y,w,h)
    
    print(s_floorplan)#for debug
    polish_exp = ' '.join( s_floorplan.Polish() )
    print('0.Polish Expression: '+polish_exp)
    
    max_tp = 0
 #   s_floorplan.M3()
    best_cost = 999
    best_temp_max = temp_am
    best_fp = ()  
 
    T = 1.0 # temperature
    T_min = 0.005
    coeff = 0.8
    
    print('Simulated Anealing started. Starting with T = '+str(T)+', will be stopped when T < '+str(T_min)+' or 2 consecutive temperature without cost improvement')
    
    old_cost = new_floorplan()
    
    max_inner_iter = 20
    nobetter_t_count = 0
    while T > T_min:
        old_best_cost = best_cost
        inner_iter_count = 0
        nojump_count = 0
        while inner_iter_count < max_inner_iter:
            inner_iter_count += 1
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
                nojump_count = 0
            else:
                s_floorplan.slicing_tree = copy.deepcopy( tmp_tree ) #recover the original tree
                nojump_count += 1
                if nojump_count >= 3:
                    break
                
        T = T*coeff
        print(inner_iter_count,'floorplan evaluated, T =', T)
        print('So far, best cost found = '+str(best_cost)+' , best floorplan: '+str(best_fp))
        
        if best_cost == old_best_cost:
            nobetter_t_count += 1
        else:
            nobetter_t_count =0
            
        if nobetter_t_count >= 2 and best_cost != 999:
            break
                 
    # draw result
    if best_cost == 999:
        print('No feasible floorplan found')
        with open('./output_files/json/final_result','w') as json_file:
            json.dump(final_result, json_file) 
        return
    else:
        print('best floorplan: '+str(best_fp))
        modules = best_fp[1]
        final_result = [best_temp_max, modules]
        with open('./output_files/json/final_result','w') as json_file:
            json.dump(final_result, json_file) 
        DrawFloorplan('./output_files/'+design_name+'_floorplan.svg',modules)
        DrawThermalMap('./output_files/'+design_name+'_thermal_map.svg',modules, 0)
        
    return modules

# check if module1 [x1,y1,w1,h1] overlap with module2 [x2,y2,w2,h2]
def CheckOverlapMod(module1, module2):
    x1,y1,w1,h1 = module1
    x2,y2,w2,h2 = module2
    
    if x1 < x2+w2 and x1+w1 > x2 and y1 < y2+h2 and y1+h1 > y2:
        return True
    else:
        return False

# check if new_module [x,y,w,h] overlap with modules in floorplan [['name',[x,y,w,h]],[..],..]
def CheckOverlap(floorplan, new_module):
    for mod in floorplan:
        if CheckOverlapMod(mod[1], new_module):
            return True
            
    return False
    
# floorplan format: [['name',[x,y,w,h]],[..],..]
def BoundingArea(floorplan):
    max_x = 0
    max_y = 0
    for mod in floorplan:
        right_x = mod[1][0]+mod[1][2]
        upper_y = mod[1][1]+mod[1][3]
        if right_x > max_x:
            max_x = right_x
        if upper_y > max_y:
            max_y = upper_y
            
    return max_x*max_y
    
# Euclidean distance of the centers of two modules, module format:[x,y,w,h]
def EuclideanDist(module1,module2):
    center1_x = module1[0] + module1[2]/2
    center1_y = module1[1] + module1[3]/2
    center2_x = module2[0] + module2[2]/2
    center2_y = module2[1] + module2[3]/2
    
    Edist = ((center1_x - center2_x)**2 + (center1_y - center2_y)**2 )**0.5
    return Edist
    
# module format: ['name',[x,y,w,h]] 
def PowerDensity(module):
    power = module_list[module[0]][1]
    area = module[1][2]*module[1][3]
    
    return power/area

# floorplanning algorithm based on (modified) cluster growth
def FloorplaningMCG():
    global alpha, final_result
    alpha = 1
    
    # Cost function when select module locations, module: [x,y,w,h]
    def cost_func(module):
        if CheckOverlap(floorplan_mcg, module):
            cost = float('inf')
        else:
            new_floorplan = floorplan_mcg + [['new_module',module]]
            bound = BoundingArea(new_floorplan)
            b_area_norm = bound/(WIDTH*HEIGHT) #normalized bounding area
            power_density = PowerDensity([current_mod_name,module])
            mod_area_self = module[2]*module[3]
            
            second_term_sum = 0
            for mod in floorplan_mcg:
                mod_area = mod[1][2] * mod[1][3]
                dist = EuclideanDist(module,mod[1])
                other_power_density = PowerDensity(mod)
                pd_diff = abs(power_density - other_power_density)
                second_term_sum += (mod_area*mod_area_self)**0.5 /(dist * (pd_diff+mcg_const))
                
            cost = mcg_alpha*b_area_norm + mcg_beta*second_term_sum             
            
        return cost
    
    # first module to place is given as seed
    def LinearOrdering(first_mod_name):
        module_order = [first_mod_name]
        mod_list = list(module_list)
        gain_list = [0]*len(mod_list)
        unselected_list = [True]*len(mod_list)
        
        selected_index = mod_list.index(module_order[-1])
        unselected_list[selected_index] = False

        for i,mod in enumerate(mod_list):
            for other_mod in module_list[mod][2]:
                gain_list[i] -= module_list[mod][2][other_mod] #new nets
                       
        while True in unselected_list:
            net_dict = module_list[ module_order[-1] ][2]
            for mod in net_dict:
                mod_index = mod_list.index(mod)
                gain_list[mod_index] += net_dict[mod] #nets going to be terminated
                
            remaining_mods = list( itertools.compress(mod_list,unselected_list) )
            remaining_gains = list( itertools.compress(gain_list,unselected_list) )
            
            module_order.append( remaining_mods[ remaining_gains.index( max(remaining_gains) ) ] )
            selected_index = mod_list.index(module_order[-1])
            unselected_list[selected_index] = False
        
        return module_order

    os.system('mkdir -p ./hotspot/mcg')

    # generate IRL for leaves (modules), generate leaf_IRL:{'fir':[[x,y,w,h]..],...}
    AllLeavesIRLGen()
    
    floorplan_list = []
    max_temp_list = []

    for module_n in list(module_list):
        module_order = LinearOrdering(module_n)
        print('order:',module_order)
        
        floorplan_mcg = [] # in formant [['name',[x,y,w,h]],[..],..]
        mcg_success = True
        
        for mod_name in module_order:
            mod_location_list = leaf_IRL[mod_name]
            current_mod_name = mod_name
            mod_cost_list = list( map(cost_func,mod_location_list) )
            min_cost_ind = mod_cost_list.index( min(mod_cost_list) )
            
            if min_cost_ind != float('inf'):
                min_cost_location = mod_location_list[min_cost_ind]
                floorplan_mcg.append([mod_name,min_cost_location])
                print(mod_name,'placed')
            else:
                mcg_success = False
                break
                
        if mcg_success == True:
            print('succeed')
            floorplan_list.append(copy.deepcopy(floorplan_mcg))
            RunHotSpot(floorplan_mcg, 0, 'mcg')
            cur_max_temp = ReadTempMax('mcg')
            max_temp_list.append(cur_max_temp)
    
    # draw result
    if not mcg_success:
        print('No feasible floorplan found')
        with open('./output_files/json/final_result','w') as json_file:
            json.dump(final_result, json_file) 
        return
    else:
        floorplan_mcg = copy.deepcopy( floorplan_list[max_temp_list.index( min(max_temp_list) )] )
        print('best floorplan:',floorplan_mcg)
        os.system('mkdir -p ./hotspot/mcg')
        DrawThermalMap('./output_files/'+design_name+'_thermal_map_mcg.svg',floorplan_mcg, 0)
        temp_max = ReadTempMax('')
        final_result = [temp_max, floorplan_mcg]
        with open('./output_files/json/final_result','w') as json_file:
            json.dump(final_result, json_file) 
        DrawFloorplan('./output_files/'+design_name+'_floorplan_mcg.svg',floorplan_mcg)
        DrawThermalMap('./output_files/'+design_name+'_thermal_map_mcg.svg',floorplan_mcg, 0)
        
    return floorplan_mcg
    

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
    
def DrawThermalMap(file_name,module_loc_list, root_power):
    RunHotSpot(module_loc_list, root_power, '')
    os.system('cd hotspot && ./grid_thermal_map.pl '+design_name+'.flp '+design_name+'.grid.steady > ../'+file_name)


design_name = 'test' #default

if __name__ == "__main__":
    SetGlobalVar(sys.argv[1:])
    
resource_file = "./"+ design_name +".res"
module_file = "./"+ design_name +".module"
exec(open(resource_file).read())
exec(open(module_file).read())

sema_hotspot = threading.Semaphore(value = PROCESS_MAX )
sema_irlgen = multiprocessing.Semaphore(PROCESS_MAX)

if __name__ == "__main__":
    main()
