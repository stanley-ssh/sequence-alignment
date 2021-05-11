
"""
Author : Anamelechi Stanley 
Assignemnet # : 1 
Date : 23/ 02/ 2021 


"""


import numpy as np 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys 
import time 

def bruteForce(pattern ,text ):
    ##This function does the brute  force method of pattern matching
    pat_array = list(pattern)
    text_array = list(text)

    counter = len(text_array) - len (pat_array)
    # return value ()
    ret_value =0 
    counting_value = 0 
    while (counting_value<= counter):
        # we slice the list based on k counter 
        
        if (text_array[counting_value :counting_value + len(pat_array)] == pat_array):
            ret_value +=1
        
        counting_value +=1
    

    return ret_value


def KMP (pattern , sequence):
    pat_array = list(pattern)
    seq_array = list(sequence)
    table = preprocessing(pattern)

    #loop_counter = len(seq_array) - len(pat_array)
    N = len(pat_array)
    M = len (seq_array)
    counter = M - N 
    number_occurence = 0 
    j =0 
    i=0 
    
    while (j < counter):

        match = False # check if we have gotten a full match 

        value_check = (seq_array[j+i] == pat_array[i])
        

        if (value_check  ):
            
            
            i +=1
            if (i == N):
                match = True
                #print ("we have a match")
                number_occurence+=1
                #print(number_occurence)
            # check if we have found a match 
           

        if  ( not value_check or  match==True):
            # update according ly
            j = j + i - table[i]
            if (table[i] < 0 ):
                i=0 
            else :
                i = table[i]

    return number_occurence
    

#uses tow indexes one for the pattern and one for text '


def preprocessing(pattern):
    pat_array = list(pattern)
    n = len(pat_array)
    table = np.zeros(n +1,dtype=int)
    i=0 
    j = table[0] =-1 
    while(i < n ):
        
        while (j> -1 and (pat_array[i]!= pat_array[j])):
            j= table[j]
        i+=1
        j+=1
        
        if ( i < n and pat_array[i] == pat_array[j]):
            table[i] = table[j]
        else :
            table[i] = j 
    
    return table

def run (sequence,pattern,algo, direct):

    # this function runs an algorithm based on the input 
    if (algo=="KMP"):
        # run the kmp algorithm
        ret =  KMP(pattern,sequence)
        output_result(ret,"KMP",direct)
    elif (algo =="BF"):
        ret=bruteForce(pattern,sequence)
        output_result(ret,"BF",direct)




def output_result(ret_value, algo_type,direct):
    #print("Using the {} algorithm this are the following results:".format(algo_type))
    if direct == "F":

       
        print("The number of occurences in the forward strand = {}".format(ret_value))
    elif direct == "R":
        print("The total number of occurences in the reverse strand = {}".format(ret_value))

def processFile(fileName,pattern,options):
    # The function will proces sthe file based oin the algorithm specified \
    algo ="KMP"

    if options is not None :
        if options =="-b":

            algo = "BF"
    
    for fastaRec in SeqIO.parse(fileName,'fasta'):
        rev = fastaRec.seq.reverse_complement() # changing the strand to recverse strand 
        
        run(fastaRec.seq,pattern, algo,"F")
        run(fastaRec.seq.reverse_complement(),pattern,algo,"R")

        

def main():
    if ( len(sys.argv) == 3 ):
        # the we use the defaulkt kmp algo 
        
        fileName = sys.argv[1]
        
        pattern = sys.argv[2]

        start_time = time.time()
        processFile(fileName,pattern,None)
        print("---Completed the search in  %s seconds ---" % (time.time() - start_time))

    elif(len(sys.argv) == 4) :
        
        fileName = sys.argv[1]
        
        pattern = sys.argv[2]
        
        algo = sys.argv[3]

        start_time = time.time()
        processFile(fileName,pattern,algo)
        print("---Completed the search in  %s seconds ---" % (time.time() - start_time))

    else :
        # if there is any thing else 
        print("Please input the right arguments.")
        
        
        
if __name__ == '__main__':
    #start_time = time.time()

    main()
    #print("The number of occurences is {}".format(KMP("ATG","ATGACGATTCCGAAACACGAGCCGCGCGAGGTCTTCGATCGGGCGATCGAGCATACGCGGGCGCTTTCCCATG")))
    #print("The preprocesisng table is {}".format(preprocessing("ABABCABAB")))

    

    #print("---Completed the search in  %s seconds ---" % (time.time() - start_time))






