from __main__ import *

import sys, os, time

import numpy as np
from itertools import permutations
import multiprocessing
from functools import partial

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import imshow, show, loadtxt

###################
# POLARIS BACKEND #
###################

################################################################################
# user inputs from GUI:

def user():
    global user_df, user_pts, user_R, user_N, user_proc, user_border, user_rate
    global pyDir, outDir
    user_df = [] #2d landscape user data
    user_pts = [] #user coordinates
    user_R = [] #user permutational orders
    user_N = [] #user segmentation depths
    user_proc = 1 #user processor count
    user_rate = False #transition rate weighting
    pyDir = os.path.dirname(os.path.realpath(__file__)) #python file location
    outDir = os.path.join(pyDir, 'data_output') #output directory
    if not os.path.exists(outDir):
        os.makedirs(outDir)

# multiprocessing worker:
def perm_energy(perms, start, end): 
    rx = []
    ry = []
    lineSegs = []
    for r in range(len(perms)): 
        rx.append(r)
        ry.append(r)
    for s in range(len(perms) + 1):
        lineSegs.append(s)
    
    tempLine = []
    tempLine.append([start[0], start[1]])
    for i in range(len(perms)): #grab coords for items in each combination
        rx[i],ry[i] = perms[i]#r[i]
        tempLine.append([rx[i],ry[i]])
    tempLine.append([end[0],end[1]]) #S->permA->permB->...->pernN->E
                       #...next loop: S->permB->permA->...->pernN->E, etc.
    # start -> p1:
    lineSegs[0] = np.vstack(line_coords(start[0],start[1],rx[0],ry[0]))
    # p1 -> p2; p2 -> p3... p_(n-1) -> p_n
    for l in range(len(perms) - 1):
        lineSegs[l+1] = np.vstack(line_coords(rx[l],ry[l],rx[l+1],ry[l+1]))
        lineSegs[l+1] = np.delete(lineSegs[l+1], (0), axis=0) #don't overcount junction
    # p_n -> end:
    lineSegs[len(perms)] = np.vstack(line_coords(rx[-1],ry[-1],end[0],end[1]))
    lineSegs[len(perms)] = np.delete(lineSegs[len(perms)], (0), axis=0) #don't overcount junction

    linePtsAll = []
    for i in range(len(perms) + 1):
        linePtsAll.append(lineSegs[i])
    linePtsAll = np.vstack(linePtsAll)
    lineErg = line_energy(linePtsAll)
    lineMean = np.mean(linePtsAll)

    return tempLine, lineErg, lineMean

################################################################################
# main algorithm:

def init():
    print('\n\
    ##########################################################################\n\
    #    POLARIS (Path of Least Action Recusive Survey)                      #\n\
    #    v.1.0 Python 3.x                                                    #\n\
    ##########################################################################\n\
    #    Copyright (C) 2018-2019 Evan Elliott Seitz                               #\n\
    ##########################################################################\n\
    #                                                                        #\n\
    #                               LICENSE                                  #\n\
    #                                                                        #\n\
    ##########################################################################\n\
    #                                                                        #\n\
    #    POLARIS is a free software: you can redistribute it and/or modify    #\n\
    #    it under the terms of the GNU General Public License as published   #\n\
    #    by the Free Software Foundation, either version 2 of the License,   #\n\
    #    or (at your option) any later version.                              #\n\
    #                                                                        #\n\
    #    This program is distributed in the hope that it will be useful,     #\n\
    #    but WITHOUT ANY WARRANTY; without even the implied warranty of      #\n\
    #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the        #\n\
    #    GNU General Public License for more details.                        #\n\
    #                                                                        #\n\
    #    You should have received a copy of the GNU General Public License   #\n\
    #    along with this program, found at <http://www.gnu.org/licenses>.    #\n\
    #                                                                        #\n\
    #                                                                        #\n\
    ##########################################################################\n\
    #   Contact: evan.e.seitz@gmail.com                                      #\n\
    ##########################################################################\n\
    #   DOI: 10.1021/acs.jcim.9b01108                                        #\n\
    ##########################################################################\n\
    ')
    
################################################################################
# add temporary frame to landscape:

    def add_border(df): #create padding around manifold for subdivisions

        if np.shape(df)[0] % 2 != 0: #trim data to even if x is odd
            df = np.delete(df, np.shape(df)[0]-1, 0)
        if np.shape(df)[1] % 2 != 0: #trim data to even if y is odd
            df = np.delete(df, np.shape(df)[1]-1, 1)

        df_dim = (int(np.shape(df)[0]), int(np.shape(df)[1]))

        ii = 0
        while True:
            if 2**ii >= max(df_dim):
                upper_bound = 2**ii
                break
            else:
                ii+=1

        sideY = int((upper_bound - int(np.shape(df)[0])) / 2)
        sideX = int((upper_bound - int(np.shape(df)[1])) / 2)
        x_box = np.full((upper_bound, sideX), np.amax(df))
        y_box = np.full((sideY, int(np.shape(df)[1])), np.amax(df))
        df1 = np.concatenate((y_box,df), axis=0)
        df2 = np.concatenate((df1,y_box), axis=0)
        df3 = np.concatenate((x_box,df2), axis=1)
        df_framed = np.concatenate((df3,x_box), axis=1)

        lenX, lenY = np.shape(df_framed)
        maxE = np.amax(df_framed)
        
        return df_framed, lenX, sideX, sideY
    
################################################################################
# prepare data:

    try:
        LS2d = loadtxt(user_df,float,delimiter=',') #load datafile
    except ValueError:
        LS2d = loadtxt(user_df,float)

    if user_rate is True: #transition state weighting enabled
        idx = 0
        for i in LS2d:
#            LS2d[idx] = LS2d[idx]**2
            LS2d[idx] = 2**(LS2d[idx]) - 1
            idx += 1

    df_framed, lenX, sideX, sideY = add_border(LS2d) #reframe data

    ii = 0
    for j,k in user_pts: #shift coords to match reframing
        ii += 1
        user_pts[ii-1] = [j+sideX, k+sideY]

    # needed for multiprocessing:
    global line_coords, line_energy
               
################################################################################
# least action recursions:

    def least_action(start, end, master):

        print('\n\n')
        print('#########################################')
        print('# NEW RECURSION:')
        dist = np.sqrt((start[0]-end[0])**2+(start[1]-end[1])**2)
        print('# Euclidean Distance:',dist)
        
        R_minLine = []
        R_minLineEnergy = []
        N_index = []

        # 0th order permutation (energy of start->end line only):
        R_minLine.append([[start[0],start[1]],[end[0],end[1]]])
        R_minLineEnergy.append(line_energy(np.vstack(line_coords(start[0],start[1],end[0],end[1]))))
        N_index.append(0)
        
        print('#\n# RUNNING R0, N0...')
        print('#  N0 minimum line:',R_minLine[-1])
        print('#  N0 minimum line energy:',np.round(R_minLineEnergy[-1],2))

        for u in range(len(user_R)):
            R = user_R[u]
            N = user_N[u]

            print('#\n# RUNNING R%s, N%s...' % (R,N))
            if dist < (np.shape(df_framed)[0]/(2**N)):
                print('#  PASS')
                pass #skip for distances shorter than segmentation width
            
            else:
                ti = time.time() #start timer

                minLine = []
                blocks_min = image_segment(start,end,N)
                r_minLine = []
                r_minLineEnergy = []
                
                minLine, minLineEnergy = line_permute(start,end,blocks_min,N,R)
                r_minLine.append(minLine)
                r_minLineEnergy.append(minLineEnergy)
                n_min_index = np.argmin(r_minLineEnergy) #only keep lowest
                R_minLine.append(r_minLine[n_min_index])
                R_minLineEnergy.append(r_minLineEnergy[n_min_index])
                N_index.append(N)

                tf = time.time() #end timer
                timer = round((tf-ti)/60,2) #function duration
                print('#  R%s, N%s duration: %s min' % (R,N, timer))
                print('#  R%s, N%s minimum line:' % (R,N),R_minLine[-1])
                print('#  R%s, N%s minimum line energy:' % (R,N),np.round(R_minLineEnergy[-1],2))

        # loop to grab out first instance of minima:
        minima = []
        N_min_idx = []
        idx = 0
        for i in R_minLineEnergy:
            if not minima:
                minima.append(i)
                N_min_idx.append(idx)
            elif minima:
                if i < np.amin(minima): #grab lowest index (via n)
                    minima.append(i)
                    N_min_idx.append(idx)
            idx += 1

        N_minLine = R_minLine[N_min_idx[-1]]
        N_minLineEnergy = R_minLineEnergy[N_min_idx[-1]]
        N_index = N_index[N_min_idx[-1]]

        print('#\n# BEST = N%s:' % (N_index))
        print('#  minimum line:', N_minLine)
        print('#  minimum line energy:', np.round(N_minLineEnergy,2))
        print('#########################################')
        
        # branching recursion:
        for i in range(len(N_minLine)-1):
            x0,y0 = N_minLine[i]
            x1,y1 = N_minLine[i+1]

            if N_index == 0:
                master.append(line_coords(x0,y0,x1,y1))
            
            else: #recurse:
                least_action(N_minLine[i],N_minLine[i+1],master)
                            
################################################################################
# image segmentation:

    def image_segment(start,end,N):
        shift_x = 0
        shift_y = 0
        blocks_min = [] #all local minima across all quadrants

        for b in range(4**N): #number of segmented blocks
            block_all = [] #all coords in block, no filters
            block_coords = [] #initial coordinates per block
            block_ergs = [] #energy per coordinate in block_coord
            for x in range(0 + shift_x, int(lenX/(2**N)) + shift_x): #dim of each block
                for y in range(0 + shift_y, int(lenX/(2**N)) + shift_y):
                    if ([x,y] != [start[1],start[0]]) and\
                       ([x,y] != [end[1],end[0]]): #don't record input start and end pts
                            block_all.append([x,y])
                            block_coords.append([x,y]) #coords for pixels in one quadrant
                            block_ergs.append(df_framed[x,y]) #energy of pixels in one quadrant               
            ii = 0
            for i in block_ergs:
                ii += 1 #count upward in block_coord to match block_erg index
                if i != np.amax(LS2d) and i == np.amin(block_ergs): #ignore global max, include local min
                    blocks_min.append(block_coords[ii-1]) #keeping track of index between coord/erg lists
                    break #only one minima from each segment, others captured later by larger N
                        
            shift_x += int(lenX/(2**N)) #permute across all blocks to find min in each
            if shift_x == int(lenX):
                shift_x = 0
                shift_y += int(lenX/(2**N))

        ii = 0
        for y,x in blocks_min:
            ii += 1
            blocks_min[ii-1] = [x,y]

        return blocks_min

################################################################################
# obtain line coordinates btw any 2 points:

    def line_coords(x0,y0,x1,y1):
        def line_maker(x0,y0,x1,y1):
            "Bresenham's Line Algorithm"
            linePts = []
            dx = abs(x1 - x0)
            dy = abs(y1 - y0)
            x, y = x0, y0
            sx = -1 if x0 > x1 else 1
            sy = -1 if y0 > y1 else 1
            if dx > dy:
                err = dx / 2.0
                while x != x1:
                    linePts.append((x,y))
                    err -= dy
                    if err < 0:
                        y += sy
                        err += dx
                    x += sx
            else:
                err = dy / 2.0
                while y != y1:
                    linePts.append((x,y))
                    err -= dx
                    if err < 0:
                        x += sx
                        err += dy
                    y += sy       
            linePts.append((x,y))
            return linePts, line_energy(linePts)

        # fix for slight discrepancy btw forward and reverse lines:
        F_line, F_energy = line_maker(x0,y0,x1,y1)
        R_line, R_energy = line_maker(x1,y1,x0,y0)
        if F_energy <= R_energy:
            return F_line
        else:
            return R_line[::-1]
        
################################################################################
# integrated energy from line_coords:      
        
    def line_energy(linePts):
        lineEnergy = []
        for i in range(len(linePts)):
            j, k = linePts[i]
            E = df_framed[k,j]
            lineEnergy.append(E)

        return sum(lineEnergy)

################################################################################
# find all line permutations for a given order:

    def line_permute(start, end, blocks_min, N, R):
        # first check to see that factorials can even be done...
        # ...(error if too many max_energy regions were ignored above):
        if (len(blocks_min) - R) < 0:
            # abort function:
            minLine = [[0,0],[0,0],[0,0]]
            minLineEnergy = 999999999
            print('#  INSUFFICIENT POINTS')
            return minLine, minLineEnergy

        else: #run through permutations
            tempLines = []
            lineErgs = []
            lineMeans = []
    
            pool = multiprocessing.Pool(processes=user_proc)
            perms = permutations(blocks_min, R)
            for r in pool.imap_unordered(partial(perm_energy, start=start, end=end), perms):
                tempLine, lineErg, lineMean = r
                tempLines.append(tempLine)
                lineErgs.append(lineErg)
                lineMeans.append(lineMean)

            min_1st = np.argmin(lineErgs) #returns first instance of minimum
            min_idxs = np.where(lineErgs == np.amin(lineErgs)) #count number of equal-energy minima
            degen = np.shape(min_idxs)[1] - 1

            if degen > 0: #if more than one equal-minimum-energy option found
                lowest_mean = np.amax(lineMeans)
                lowest_mean_idx = 0
                for idx in min_idxs[0]: #e.g., idx = [781 822 909]
                    if lineMeans[idx] < lowest_mean:
                        lowest_mean = lineMeans[idx]
                        lowest_mean_idx = idx
                minLine = tempLines[lowest_mean_idx]
                minLineEnergy = lineErgs[lowest_mean_idx]
                
            else:
                minLine = tempLines[min_1st]
                minLineEnergy = lineErgs[min_1st]
                    
            return minLine, minLineEnergy


################################################################################
# slice out duplicate points (if ever any):

    def slicer(master):
        master_unique = []
        for i in np.vstack(master): 
            if [i[0],i[1]] not in master_unique:
                master_unique.append([i[0],i[1]])

        return np.vstack(master_unique)
    
################################################################################
# apply discrete perturbations across globally defined path
# to optimize energies locally:

    def wiggler(master):
        master = np.vstack(master)
        masterEnergy = line_energy(master)
        lineCoords = []
        lineEnergy = []
        minLineCoords = []
        minLineEnergy = []
        temp = []
        temp = list(master)
        temp = np.vstack(temp)
        index = 0
        
        for i in range(1,len(master)-1): #don't perturb start or end point
            j,k = master[i]
            Wblock, Eblock, Sblock, Nblock = False, False, False, False
            # delete excess pts (if ever any):
            temp1 = temp[:i]
            temp2 = temp[i+1:]
            temp = np.concatenate((temp1,temp2),axis=0)
            gap = gap_finder(temp)
            if (line_energy(temp) < masterEnergy) and (gap == False): # <= if deleting E=0 pts
                lineCoords.append(temp)
                lineEnergy.append(line_energy(temp))
                index += 1
            temp = list(master)
            temp = np.vstack(temp)
            # west occupation:
            for [t,u] in master:
                if [j-1,k] == [t,u]:
                    Wblock = True
                    break
            # west perturbation:
            if Wblock is False:
                temp[i][0] = j-1
                gap = gap_finder(temp)
                if (line_energy(temp) < masterEnergy) and (gap == False):
                    lineCoords.append(temp)
                    lineEnergy.append(line_energy(temp))
                    index += 1
                temp = list(master)
                temp = np.vstack(temp)
            # east occupation:
            for [t,u] in master:
                if [j+1,k] == [t,u]:
                    Eblock = True
                    break
            # east perturbation:
            if Eblock is False:
                temp[i][0] = j+1
                gap = gap_finder(temp)
                if (line_energy(temp) < masterEnergy) and (gap == False):
                    lineCoords.append(temp)
                    lineEnergy.append(line_energy(temp))
                    index += 1
                temp = list(master)
                temp = np.vstack(temp)           
            # south occupation:
            for [t,u] in master:
                if [j,k-1] == [t,u]: 
                    Sblock = True
                    break
            # south perturbation:
            if Sblock is False:
                temp[i][1] = k-1
                gap = gap_finder(temp)
                if (line_energy(temp) < masterEnergy) and (gap == False):
                    lineCoords.append(temp)
                    lineEnergy.append(line_energy(temp))
                    index += 1
                temp = list(master)
                temp = np.vstack(temp)
            # north occupation:
            for [t,u] in master:
                if [j,k+1] == [t,u]: 
                    Nblock = True
                    break
            # north perturbation:
            if Nblock is False:
                temp[i][1] = k+1
                gap = gap_finder(temp)
                if (line_energy(temp) < masterEnergy) and (gap == False):
                    lineCoords.append(temp)
                    lineEnergy.append(line_energy(temp))
                    index += 1
                temp = list(master)
                temp = np.vstack(temp)
                    
        ii = 0
        for i in lineEnergy:
            ii += 1
            if i == np.amin(lineEnergy):
                minLineCoords = lineCoords[ii-1] 
                minLineEnergy = np.amin(lineEnergy)
                break
            
        if minLineEnergy == []:
            return master
        else:
            return wiggler(minLineCoords)

    # check for gaps in line due to perturbations:
    def gap_finder(master):
        gap = False
        for i in range(0,len(master)):
            if i == len(master)-1:
                break
            else:
                j,k = master[i]
                m,n = master[i+1]
                dist = np.sqrt((j - m)**2 + (k - n)**2)
                if dist > np.sqrt(2):
                    gap = True
        return gap

################################################################################
# 2d plot of energy vs. coupled reaction coordinates:

    def curves_plot(path, timestr):
        fig2, ax2 = plt.subplots()
        pathLen = len(path)
        X = np.empty([pathLen,1])
        Y = np.empty([pathLen,1])
        XY = np.empty([pathLen,1])
        Z = np.empty([pathLen,1])

        ii = 0
        for i in range(pathLen):
            ii+=1
            j,k,l = path[i]
            X[i] = [int(j)]
            Y[i] = [int(k)]
            XY[i] = ii-1
            Z[i] = [float(l)]

        ax2.plot(XY,Z)
        ax2.set_xlabel('Coupled Reaction Coordinates (RC1, RC2)')
        ax2.set_ylabel('Energy')
        fig2.tight_layout()
        outname_EP = os.path.join(outDir, 'POLARIS_EP_%s.png' % (timestr))
        fig2.savefig(fname=outname_EP, bbox_inches='tight', dpi=200)
        ax2.clear()
        fig2.clf()
        plt.close(fig2)
        
################################################################################
# call functions and plot:

    ti = time.time() #start timer
    master = []
    master_all = []

    if len(user_pts) > 2:
        for i in range(0,len(user_pts)-1):
            start = user_pts[i]
            end = user_pts[i+1]
            least_action(start, end, master)
            master_all.append(np.vstack(master))
            master = []
        master_all = np.vstack(master_all)

    elif len(user_pts) == 2:
        start = user_pts[0]
        end = user_pts[1]
        least_action(start, end, master)
        master_all = master
        master_all.append(user_pts[-1]) #add end point to list
        master_all = np.vstack(master_all)

    master_all = slicer(master_all) #remove duplicate points (if ever any)
    master_all = wiggler(master_all) #optimize line via local energy perturbation

    if user_rate is True:
        idx = 0
        for i in LS2d:
            #LS2d[idx] = LS2d[idx]**(1/2)
            LS2d[idx] = np.log(LS2d[idx]+1)/np.log(2)
            idx += 1

    final_x = []
    final_y = []
    final_e = []

    fig1, ax1 = plt.subplots()
    # rescale coordinates to original framing:
    for i in np.vstack(master_all):
        j,k = i
        j -= sideX
        k -= sideY
        final_x.append(j)
        final_y.append(k)
        final_e.append(LS2d[k,j])
        ax1.scatter([j],[k], c='k', s=10)
        ax1.scatter([j],[k], c='w', s=5)
    user_pts_final = []
    for i in user_pts:
        j,k = i
        j -= sideX
        k -= sideY
        user_pts_final.append([j,k])
        ax1.scatter([j],[k], c='k', s=.5)
        
    final_all = np.column_stack((final_x, final_y, final_e))
    final_e_all = sum(final_e)
    final_e_all = round(final_e_all,3)

    tf = time.time() #end timer
    timer = round((tf-ti)/60,2) #function duration
    timestr = time.strftime('%Y%m%d-%H%M%S')

    outname_DF = os.path.join(outDir, 'POLARIS_DF_%s.txt' % (timestr))
    np.savetxt(outname_DF,
               final_all,
               fmt='%i %i %10.4f',
               header=('x  y \t  energy'),
               comments=('POLARIS output\n\ninput points: %s\ninput settings: R: %s, N: %s\npath energy: %s\npath length: %s\nelapsed time: %s min\n\n'\
                         % (user_pts_final,user_R,user_N,final_e_all,len(np.vstack(final_e)),timer)))

    im = ax1.imshow(LS2d, cmap='jet', origin='lower')
    ax1.set_xlabel('Path Energy = %s, Path Length = %s' % (final_e_all, len(np.vstack(final_e))))
    fig1.colorbar(im, orientation='vertical')

    # save final energy landscape path:
    outname_LS = os.path.join(outDir, 'POLARIS_LS_%s.png' % (timestr))
    fig1.savefig(fname=outname_LS, bbox_inches='tight', dpi=200)
    ax1.clear()
    fig1.clf()
    plt.close(fig1)

    # save energy vs. coords:
    curves_plot(final_all, timestr)

    print('\n\nPOLARIS COMPLETE')
    print('Computation time: %s minutes\n\n' % (timer))


if __name__ == '__main__':
    init()
