#!/usr/bin/python 

#Ibai scripts
#02/07/2018


from subprocess import Popen, PIPE
import random
from collections import defaultdict
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
import sys
from matplotlib.cbook import flatten

#check also bedtools: multiIntersectTest for many datasets
#we could do a parser for multiIntersect



# function that counts number of occurrences from the insertBed tool output
def count_occurrences_from_output(output):
    string = ""
    lista = []
    for line in output:
        string = string + line
        if line == "\n":
            lista.append(string)
            string = ""
    counter = 0

    for line in lista:
        counter+=1
    return counter

# the same as above but for many datasets
#MultiIntersect output
#chr1	6	8	1	1	1	0	0
#chr1	8	12	2	1,3	1	0	1
#chr1	12	15	3	1,2,3	1	1	1
#chr1	15	20	2	1,2	1	1	0

#output for 3 datasets(A,B,C):000,001,010,011,100,101,110,111 -> 0,1,2,3,4,5,6,7
#output for 4 datasets(A,B,C,D):0000,0001,0010,0011,... -> 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
#the output will be delivered in a list, starting from 0

#we can also if needed use  int("10011",2) or "{0:b}".format(42)

def count_multiple_occurrences_from_output(output,n_datasets):
    string = ""
    lista = []
    n = 2**n_datasets
    output_result = [0]*n
    popped_result = [0]*n
    for line in output:
        string = string + line
        if line == "\n":
            lista.append(string)
            string = ""

    #now, merge the contiguous regions
    #for cases like this.
    #chrX	162009842	162009953	1	2	0	1	0
    #chrX	162009953	162010406	2	2,3	0	1	1
    #chrX	162010406	162010563	1	3	0	0	1

    #for each line, populate our final list
    previous_buffer = [0,0,0,0,0,0,0,0]
    max_num_overlap = 0
    merged_bed = "merged_occurrences.bed"
    pop_next = False
    for line in lista:
        values = line.split("\t")
        occurrences = "".join(line.split("\t")[5:])
        binary = int(occurrences,2)
        output_result[binary] += 1 
        if int(previous_buffer[2]) == int(values[1]):
            print max_num_overlap
            if int(values[3]) > max_num_overlap:  #we want to count like one enhancers in these format: ####    ###    these would be 1 enhancer, not 2
                                                                                                 #############
                                            # but if the ones above are from different datasets? fix that.
                if not pop_next:
                    print "X B ({})".format(last_binary)
                    output_result[last_binary] -= 1
                    popped_result[last_binary] += 1
                    saved = binary
                else:
                    print "X B previous({})".format(saved)
                    output_result[saved] -= 1
                    popped_result[saved] += 1

                max_num_overlap = int(values[3])
            else:
                pop_next = True
                output_result[binary] -= 1 
                print "X ({})".format(binary),
                popped_result[binary] += 1
        else:
            if pop_next:
                pop_next = False
                max_num_overlap = 0
        last_binary = binary 
        print values,output_result,popped_result
        previous_buffer = values
    multiplying_factor = [0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1]
    try:
        sum_ = 0
        for i in range(len(popped_result)):
            sum_ = sum_ + popped_result[i]*multiplying_factor[i]
        print "TOTAL 'DELETED' REGIONS: {}".format(sum_)
    except:
        pass
    print "POPPED LIST: {}".format(popped_result)
    try:
        sum_ = 0
        for i in range(len(output_result)):
            sum_ = sum_ + output_result[i]*multiplying_factor[i]
        print "TOTAL REGIONS: {}".format(sum_)
    except:
        pass
    print "OUTPUT LIST: {}".format(output_result)
    return output_result


 
#DOES NOT WORK STILL
# determines p-value for two groups over 10000 iterations
def pvalue2(N,a,b,times,sampling):
    counter = 0
    dic = defaultdict(int)
    for i in range(sampling):
        A = random.sample(range(N),a)
        B = random.sample(range(N),b)
        cooccur = len(set(A).intersection(B))
        if cooccur in dic:
            dic[cooccur] += 1
        else:
            dic[cooccur] = 1        
    for x in dic:
        ratio = float(dic[x])/sampling
        print str(x)+":"+str(ratio)
    return float(dic[times])/sampling


# function that uses bedtools_intersect to get the number of overlapping between two datasets
def get_intersect (comp1,comp2):
    process = Popen([bedtools+"intersectBed","-u","-a",comp1,"-b",comp2], stdout=PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()
    occurrences = count_occurrences_from_output(output)

    process = Popen(["wc","-l",comp1],stdout=PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()
    total1 = int(output.split()[0])
    print comp1+": "+str(total1)

    process = Popen(["wc","-l",comp2],stdout=PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()
    total2 = int(output.split()[0])
    print comp2+": "+str(total2)
    print "intersection: "+str(occurrences)
    print "#----------------#" 
    return (total1-occurrences,total2-occurrences,occurrences)

    #print "Pvalue: "+str(pvalue2(total1+total2,total1,total2,occurrences,10000))
    #print pvalue2(10,9,9,1,1000)

# Intersect for more than 2 datasets. Function only support 3 and 4 now.
def get_multi_intersect(comp_list):

    if len(comp_list) == 3:
        process = Popen([bedtools+"multiIntersectBed", "-i",comp_list[0],comp_list[1],comp_list[2]], stdout=PIPE)
    elif len(comp_list) == 4:
        process = Popen([bedtools+"multiIntersectBed", "-i",comp_list[0],comp_list[1],comp_list[2],comp_list[3]], stdout=PIPE)
    elif len(comp_list) == 5:
        process = Popen([bedtools+"multiIntersectBed", "-i",comp_list[0],comp_list[1],comp_list[2],comp_list[3],comp_list[4]], stdout=PIPE)
    elif len(comp_list) == 6:
        process = Popen([bedtools+"multiIntersectBed", "-i",comp_list[0],comp_list[1],comp_list[2],comp_list[3],comp_list[4],comp_list[5]], stdout=PIPE)
    else:
        print "Check comparison list size."
        sys.exit(0)
    (output, err) = process.communicate()
    exit_code = process.wait()
    occurrences = count_multiple_occurrences_from_output(output,len(comp_list))
    return occurrences



# venn scaler for generate_venn2
def set_venn_scale(ax, true_area, reference_area):
    s = np.sqrt(float(reference_area)/true_area)
    ax.set_xlim(-s, s)
    ax.set_ylim(-s, s)

# function to generate venn diagrams for 2 datasets
def generate_venn2(nA_list,nB_list,nBoth_list,nameA,nameB_list,layout_line,layout_col):
    data = []
    colors = []
    for i in range(len(nB_list)):
        tuple_ = (nA_list[i],nB_list[i],nBoth_list[i])
        data.append(tuple_)

    max_area = max(map(sum, data))
    colors=["green","green","green","green","red","red","red","red","red"]
    
    figure, axes = plt.subplots(layout_line, layout_col, figsize=(30,30))
    #figure, axes = plt.subplots(len(nB_list), 1, figsize=(12,12))
    #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.01, hspace=0.01)
    line = 0
    column = 0
    for i in range(len(nB_list)):
        v = venn2(subsets=data[i], set_labels = (nameA.split("/")[-1][:20], nameB_list[i].split("/")[-1][:20]), ax=axes[line][column])
#        v = venn2(subsets=data[i], ax=axes[line][column])
        v.get_patch_by_id('100').set_color('lime')
        v.get_patch_by_id('010').set_color(colors[i])
        #plt.title(nameA+" vs "+nameB_list[i].split("/")[-1])
        column += 1
        print line,column
        if column == layout_col:
            line += 1
            column = 0
    for a, d in zip(flatten(axes), data):
        set_venn_scale(a, sum(d),max_area)
    #plt.show()
    plt.savefig("venn2.png")
# function to generate venn diagrams
def generate_venn3(binary_values,name1,name2,name3):

    binary_values=binary_values[:-1] #take 0 and last value
    values = tuple(binary_values)
    plt.figure(figsize=(12,12))
    v = venn3(subsets=binary_values, set_labels = (name1.split("/")[-1],name2.split("/")[-1],name3.split("/")[-1]))
    #v.get_patch_by_id('100').set_alpha(1.0)
    #v.get_patch_by_id('100').set_color('white')
    #v.get_label_by_id('100').set_text('Unknown')
    #v.get_label_by_id('A').set_text('Set "A"')
    #c = venn3_circles(subsets=binary_values, linestyle='dashed')
    #c[0].set_lw(1.0)
    #c[0].set_ls('dotted')
    #plt.title("Sample Venn diagram")
    #plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
    #             ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    #             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    plt.show()

#### MAIN BODY ####

# PATHS


total1_list = []
total2_list = []
occur_list = []
for i in range(len(comp2)): 
    total1,total2,occurrences = get_intersect(comp1,comp2[i])
    total1_list.append(total1)
    total2_list.append(total2)
    occur_list.append(occurrences)

#generate_venn2(total1_list,total2_list,occur_list,comp1,comp2,2,2)



# 46C
comp_list = [extragenic_data+"46c_31_day0_4h8.day0_dig.BCP.extragenic_mESC_complete_ucsc.mm10.BCP_FP5.classes_specified.no_header.colors_specified.bed",
            extragenic_data+"46c_33_day1.4h8.day0_dig.BCP.extragenic_mESC_complete_ucsc.mm10.4h8.46c_33_day1.threshold.bed",
            extragenic_data+"46c_33_day3.4h8.day0_dig.BCP.extragenic_mESC_complete_ucsc.mm10.4h8.46c_33_day3.threshold.bed",
            extragenic_data+"46c_59_day16.4h8_AB.BCP.extragenic_mESC_complete_ucsc.mm10.4h8_AB.46c_59_day16.threshold.bed",
            extragenic_data+"46c_59_day30.4h8.BCP.extragenic_mESC_complete_ucsc.mm10.4h8.46c_59_day30.threshold.bed"]
counter = 0
for dataset in comp_list:
    counter += 1
    print str(counter)+": "+dataset
value_list = get_multi_intersect(comp_list)
print comp_list
#print value_list
binary_list = []
n = len(comp_list)
for i in range(n**2):
    binary_list.append("{0:b}".format(i))
print binary_list

#this is wrong
#generate_venn3(value_list,comp_list[2],comp_list[1],comp_list[0]) #CHECK ORDER, FIRST 2, THEN 1 THEN 0
















