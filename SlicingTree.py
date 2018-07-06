#!/usr/bin/env python3
""" author: Jianqi Chen """
import random

# parameter of Node: (t,p,l,r,x,y,w,h)
# t = type(veritcal or horizontal), p = parent, l = left child, r = right child
# x = x coordinate of left bottom corner, y = y coordinate, w = width, h = height
# all indices(p,l,r) start from 1 instead of 0
class Node:
    #the list 'IRL' will be defined for every node during evaluation
    def __init__(self, t, p, l, r):
        self.t = t # 'V'(verical slicing),'H'(horizontal) or module name(string)
        self.p = p # integer number(index+1 of slicing_tree in STree), None for root
        self.l = l # integer number, None for leaves
        self.r = r # integer number, None for leaves

# modules : list of all module names
# n : module number, the length of tree(the list) is 2*n - 1
class STree:
    def __init__(self, modules):
        self.mod_num = len(modules)
        self.slicing_tree = []
        for i in range(1,2*self.mod_num):
            if i < self.mod_num:
                node_type = random.choice (['V','H'])
            else:
                node_type = modules[i - self.mod_num]
                
            if i == 1:
                parent = None
            else:
                parent = int(i/2)
                
            if i < self.mod_num:
                left_child = i*2
                right_child = i*2 + 1
            else:
                left_child = None
                right_child = None
                
            self.slicing_tree.append( Node(node_type, parent, left_child, right_child) )
            
        self.MakeNormalized()
    
    # return stardard array representation of the slicing tree, using dict instead of list
    # (p,l,r) are meaningless in this case        
    def TreeArray(self):
        bin_tree_array = {} # the binary tree in standard array representation
        visit = [False for i in range(len(self.slicing_tree))]
        bin_tree_array[1] = self.slicing_tree[0]
        
        i = 0
        while False in visit:
            if i in bin_tree_array:
                visit[ self.slicing_tree.index(bin_tree_array[i]) ] = True
                if bin_tree_array[i].l != None and bin_tree_array[i].r != None:
                    bin_tree_array[i*2] = self.slicing_tree[ bin_tree_array[i].l - 1 ]
                    bin_tree_array[i*2+1] = self.slicing_tree[ bin_tree_array[i].r - 1 ]
                    
            i = i+1
                    
        return bin_tree_array
      
    # show the tree, put a space when there is no node        
    def __str__(self):
        bin_tree_array = self.TreeArray()
        return_str = ''
        
        line_count = 1
        for i in range( 1, max(bin_tree_array)+1 ):
            if i in bin_tree_array:
                return_str += bin_tree_array[i].t
            else:
                return_str += ' '
                
            if (i+1) == (2**line_count):
                return_str += '\n'
                line_count += 1
            else:
                return_str +=','
                
        return return_str
    
    # return a list in which nodes are put in the order of Polish expression of the slicing tree (DFS)
    def DFS(self):
        polish_list = []
        def DFSTraverse(node):
            nonlocal polish_list
            if node.l == None and node.r == None:
                polish_list.append(node)
            else:
                DFSTraverse( self.slicing_tree[node.l-1] )
                DFSTraverse( self.slicing_tree[node.r-1] )
                polish_list.append(node)
        
        DFSTraverse( self.slicing_tree[0] )
        return polish_list
    
    # return Polish expression in a list
    def Polish(self):
        polish_list = self.DFS()
        return [ node.t for node in polish_list ]
    
    # a Polish expression is normalized iff there is no consecutive H or V
    # if the Polish expression of the slicing tree is normalized, return True
    def CheckNormalized(self):
        polish_exp = self.Polish()
        for i in range( 1, len(polish_exp) ):
            pre_node = polish_exp[i-1]
            cur_node = polish_exp[i]
            if pre_node == cur_node:
                return False
        
        return True
    
    # if the Polish expression of the slicing tree is not normalized, this
    # function will make it normalized by inverting H or V's
    def MakeNormalized(self):
        if self.CheckNormalized():
            return
        
        # only when it's not normalized
        polish_list = self.DFS()
        for i in range( 1, len(polish_list) ):
            pre_node = polish_list[i-1]
            cur_node = polish_list[i]
            if pre_node.t == 'V' and cur_node.t == 'V':
                cur_node.t = 'H'
            elif pre_node.t == 'H' and cur_node.t == 'H':
                cur_node.t = 'V'
                
    # remove IRLs from all nodes
    def ClearIRL(self):
        for node in self.slicing_tree:
            node.IRL = None
            del node.IRL
    
    # Slicing Tree movement M1: Swap two adjacent operands (modules), randomly
    def M1(self):
        self.ClearIRL()
        op_list = list(filter(lambda p: p.t!='H' and p.t!='V', self.DFS() ))
        index = random.randint(0, len(op_list)-2)
        
        tmp = op_list[index].t
        op_list[index].t = op_list[index+1].t
        op_list[index+1].t = tmp
        
    # Slicing Tree movement M2: complement chain of nonzero length, randomly. see definitions in D.F.Wong 1986 DAC paper 
    def M2(self):
        self.ClearIRL()
        polish_list = self.DFS()
        chain_list = [] #list of lists of nodes (chains are lists )
        
        i = 0
        while i<len(polish_list):
            chain = []
            while polish_list[i].t == 'H' or polish_list[i].t == 'V':
                chain.append(polish_list[i])
                i += 1
                if i>=len(polish_list):
                    break
            
            if chain != []:
                chain_list.append(chain)
            else:
                i += 1
            
        index = random.randint(0,len(chain_list)-1)
        for node in chain_list[index]:
            if node.t == 'H':
                node.t = 'V'
            elif node.t == 'V':
                node.t = 'H'
     
    # rebuild the tree from Polish Expression
    def BuildTree(self, polish_list):        
        def BuildNode(cur_index):
            nonlocal index
            if polish_list[cur_index].t != 'H' and polish_list[cur_index].t != 'V':
                polish_list[cur_index].l = None
                polish_list[cur_index].r = None
                return
                
            else:
                polish_list[cur_index].r = len(polish_list) - (index-1) # reversed order for the tree array
                index = index - 1
                BuildNode(index)
                polish_list[cur_index].l = len(polish_list) - (index-1)
                index = index - 1
                BuildNode(index)
                    
        # assign left and right child index for all nodes
        index = len(polish_list)-1
        
#        print('index: '+str(index)) #for debug
#        for i in range(len(polish_list)):
#            print( [ polish_list[i].t, polish_list[i].l, polish_list[i].r, polish_list[i].p] )
            
        BuildNode(index)
        
        # assign parent index for all nodes
        for i in range(len(polish_list)):
            if polish_list[i].l != None:
                polish_list[ len(polish_list) - polish_list[i].l ].p = len(polish_list) - i
                polish_list[ len(polish_list) - polish_list[i].r ].p = len(polish_list) - i
                
        # rebuild the tree
        self.slicing_tree = []
        for i in reversed(range(len(polish_list))):
            self.slicing_tree.append(polish_list[i])
            
#        print('index: '+str(index)) #for debug
#        for i in range(len(polish_list)):
#            print( [ polish_list[i].t, polish_list[i].l, polish_list[i].r, polish_list[i].p] )

                
           
    # Slicing Tree movement M3: Swap two adjacent operand and operator randomly.
    # if return true, succeed. if return false, failed.
    def M3(self):
        # exchange 2 nodes in the Polish expression, without adjusting the tree accordingly
        def Exchange(polish_list, index1, index2):
            tmp_node = polish_list[index1]
            polish_list[index1] = polish_list[index2]
            polish_list[index2] = tmp_node
            
        # check if the Polish expression does NOT have consecutive 'H' or 'V's
        def CheckNoCons(polish_list):
            for i in range( 1, len(polish_list) ):
                pre_node = polish_list[i-1]
                cur_node = polish_list[i]
                if pre_node.t == cur_node.t:
                    return False
        
            return True
        
        # check if the Polish expression does NOT violate the balloting property, for M3 at index i
        def CheckBalloting(polish_list, i):
            count = 0 #number of operators
            for j in range(i+2):
                if polish_list[j].t == 'H' or polish_list[j].t == 'V':
                    count += 1
            
            if 2*count < i+1:
                return True
            else:
                return False
        
        self.ClearIRL()
        polish_list = self.DFS()
        ret_val = False
        
        upper_ind_limit = len(polish_list)-3
        while upper_ind_limit >= 0:
            start_index = random.randint(0,upper_ind_limit)
            i = start_index
            while i<=len(polish_list)-3:
                if ((polish_list[i].t == 'H' or polish_list[i].t == 'V') and (polish_list[i+1].t != 'H' and polish_list[i+1].t != 'V')) or ((polish_list[i].t != 'H' and polish_list[i].t != 'V') and (polish_list[i+1].t == 'H' or polish_list[i+1].t == 'V')):
                    if CheckBalloting(polish_list,i):
                        Exchange(polish_list, i, i+1)
                        if CheckNoCons(polish_list):
                            self.BuildTree(polish_list)
                            ret_val = True
                            break
                        
                i += 1
                
            if ret_val == True:
                break
            else:
                upper_ind_limit = start_index-1
                
        return ret_val
            
                    
            
        
    


        

