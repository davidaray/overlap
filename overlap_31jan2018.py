#######################################################################
#
# overlap.py
#
# version 0.1
# 2017-09-13
#
#  This script reads in RepeatMasker output and then identifies
#   overlapping insertions.  This version of the script prints
#   out all non-overlapping annotations and when two annotations
#   overlap, the longest annotation is printed.
#
#  input: <repeat_masker_file>
#  output: <repeat_masker_file_no_overlaps>
#
#  input and output files are currently designated in the script
#
#
########################################################################
#TO DO LIST
#
#1. Create category outputs
	#1A 	Resolution idea: When we put R's and L's in the SEQ-ID column, maybe put an indicator for which category this element falls under... from there you can search for how many "C1"s (for example... don't use C1!) there are in the output file.
#1B		Make this an option????!!!
#2. Create a "get rid of small hits" argument
#3. Summary statistics :-O
#4. percent divergence decisions after iteration 1

#load modules
import os
import re
import sys
import argparse
import random


def modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence, q_begin, q_end, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID, CATEGORY):
	"""read in RM input and convert to strings"""
	if CATEGORY:
		output_string=	str(SW_score)+"\t"+str(perc_div)+"\t"+str(perc_del)+"\t"+str(perc_ins)+"\t"+str(query_sequence)+"\t"+str(q_begin)+"\t"+str(q_end)+"\t"+str(q_left)+"\t"+str(orient)+"\t"+str(matching_repeat)+"\t"+str(class_family)+"\t"+str(r_begin)+"\t"+str(r_end)+"\t"+str(r_left)+"\t"+ str(ID)+"_"+CATEGORY
	else:
		output_string=	str(SW_score)+"\t"+str(perc_div)+"\t"+str(perc_del)+"\t"+str(perc_ins)+"\t"+str(query_sequence)+"\t"+str(q_begin)+"\t"+str(q_end)+"\t"+str(q_left)+"\t"+str(orient)+"\t"+str(matching_repeat)+"\t"+str(class_family)+"\t"+str(r_begin)+"\t"+str(r_end)+"\t"+str(r_left)+"\t"+ str(ID)
	return output_string

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 0 - Setting Arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Set arguments for command line
def get_args():
	parser=argparse.ArgumentParser()
	parser.add_argument('-i', '--repeatmasker', type=str, 
						help='Repeat Masker output file', required=True)
	parser.add_argument('-l', '--location', type=str, help='location to start', 
						required=True)
	parser.add_argument('-o', '--output', type=str, help='Output file name')
	args=parser.parse_args()
	RM=args.repeatmasker
	LOC=args.location
	OUT=args.output

	return RM, LOC, OUT
	
#Set categories
def get_category(q_begin,q_begin_i,q_end,q_end_i,query_sequence,query_sequence_i):
	
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Step 4 Identify overlaps to Category
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Overlapping repeat annotations can occur in 5 different ways
		#  defined as "categories" here

		# CATEGORY 1 OVERLAP
		# chrX   |q_begin----------- q -----------q_end|
		# chrX	          |q_begin_i------------- qi --------q_end_i|
		
		# CATEGORY 2 OVERLAP
		# chrX                 |q_begin------- q ------q_end|
		# chrX            |q_begin_i----------  qi --------q_end_i|

		# CATEGORY 3 OVERLAP
		# chrX   |q_begin----------- q -------------------q_end|
		# chrX            |q_begin_i----qi ------q_end_i|

		# CATEGORY 4 OVERLAP
		# chrX          |q_begin----------- q -----------q_end|
		# chrX   |q_begin_i----- qi -----q_end_i|
		
		# CATEGORY 5 OVERLAP
		# chrX   |q_begin---------- q -----------q_end|
		# chrX   |q_begin_i---------qi --------q_end_i|


	if (q_begin <= q_begin_i) and (q_end > q_begin_i) and (q_end < q_end_i) and (query_sequence == query_sequence_i):
		return str("CATEG1")
	elif (q_begin > q_begin_i) and (q_end < q_end_i) and ( query_sequence == query_sequence_i):
		return str("CATEG2")
	elif (q_begin < q_begin_i) and (q_end >= q_end_i) and (query_sequence == query_sequence_i):
		return str("CATEG3")
	elif (q_begin >= q_begin_i) and (q_end > q_end_i) and (q_begin < q_end_i) and (query_sequence == query_sequence_i):
		return str("CATEG4")
	elif (q_begin == q_begin_i) and (q_end == q_end_i) and (query_sequence == query_sequence_i):
		return str("CATEG5")
	else:
		return None

		
def get_hit_array(rm_input,Iteration=0):
	if Iteration==0:
		lines2skip=3
	else:
		lines2skip=0
	
	HIT_ENTRY=1
	HIT_ARRAY=[]

        for HIT in rm_input:

                #RM records start on line 3 after the header info
                if HIT_ENTRY > lines2skip:

					#format record so that its tab delimited and no leading/tr$
					HIT = HIT.replace("*","")
					HIT = re.sub(' +', '\t', HIT).rstrip().lstrip()

					#appending record to list of all other RM records
					HIT_ARRAY.append(HIT)

                HIT_ENTRY = HIT_ENTRY + 1

	print("Input file read into memory")
	return HIT_ARRAY

def sort_hit_array(HIT_ARRAY):
	#convert hit array from list of strings into list of lists
	#for every line in hit array make tab delimited
	HIT_ARRAY=[x.split("\t") for x in HIT_ARRAY]
	
	#sort list of lists
	HIT_ARRAY.sort(key=lambda x: (x[4], int(x[5])))
	
	#convert back to tab-delmited list of strings
	HIT_ARRAY=['\t'.join(x) for x in HIT_ARRAY]
	
	return HIT_ARRAY

def main():

	RM, LOC, OUT = get_args()

	#Sanity Checks!!!!
	print('The repeatMasker output file is', RM)
	print('Directory to begin is', LOC)
	print('The output file name is', OUT)

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Step 1 - initializing
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	#set to the workingdirectory and open the output file
	os.chdir(LOC)

	#Generate Initial Output File
	#OUTPUT=open(LOC + '/' + OUT,'w')

	#Read in input file
	rm_input = open(LOC + '/' + RM, 'r')

	#Generate category specific output files:
	#CAT1OUT=open(LOC + '/' + 'cat1.out', 'w')
	#CAT2OUT=open(LOC + '/' + 'cat2.out', 'w')
	#CAT3OUT=open(LOC + '/' + 'cat3.out', 'w')
	#CAT4OUT=open(LOC + '/' + 'cat4.out', 'w')
	#CAT5OUT=open(LOC + '/' + 'cat5.out', 'w')

	#Create a Lengths file that prints how many basepairs are in each overlap.
	#LENGTHS=open(LOC + '/' + 'lengths.csv', 'w')

	#Set file iterator
	arraynum=0
	IterateStop=False
	
	NEW_HIT_ARRAY=[]	
	HIT_ARRAY=get_hit_array(rm_input,Iteration=0)
	
	#sort original RM output (HIT_ARRAY) here
	print("Sorting input file")
	HIT_ARRAY=sort_hit_array(HIT_ARRAY)
	
	while not IterateStop:
		#Call super bigfunction here. make sure all the variables are defined.
		print("Beginning Iteration {}".format(arraynum))
		NEW_HIT_ARRAY=superBIGfunction(HIT_ARRAY)	
		
		#sort NEW_HIT_ARRAY here
		print("Sorting output of iteration {}".format(arraynum))
		NEW_HIT_ARRAY=sort_hit_array(NEW_HIT_ARRAY)

		#Determine whether new hit array (start,end) locations of all hits is the same as last iteration
		if [(i.split("\t")[5],i.split("\t")[6]) for i in HIT_ARRAY] == [(i.split("\t")[5],i.split("\t")[6]) for i in NEW_HIT_ARRAY]:
			IterateStop=True
	
		else:
			HIT_ARRAY=list(NEW_HIT_ARRAY)
			
			OUTPUT=open(LOC + '/' + OUT + '.' + str(arraynum) + '.out','w')
			OUTPUT.write("\n".join(NEW_HIT_ARRAY)+'\n')
			OUTPUT.close()
			arraynum+=1
			print("Iterations Complete:{}".format(arraynum))


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Step 2. Storing RM data in memory
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# cycle through repeat masker input and load each record into a list 
	#   (stored in mem)
#	for HIT in rm_input:
#
#		#RM records start on line 3 after the header info
#		if HIT_ENTRY > 3:
#
#			#format record so that its tab delimited and no leading/trailing spaces
#			HIT = re.sub(' +', '\t', HIT).rstrip().lstrip()
#
#			#appending record to list of all other RM records
#			HIT_ARRAY.append(HIT)
#
#		HIT_ENTRY = HIT_ENTRY + 1

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Step 3 Cycling through RM to identify overlaps
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Start new function here..........................

def superBIGfunction(HIT_ARRAY):
	#iterator
	i = 0

	NEW_HIT_ARRAY=[]

	#examine each "HIT" (RM insertion) in the "HIT_ARRAY" (repeat_masker file)
	for HIT in HIT_ARRAY:
		
		#cleaning HIT of all white spaces and making tab delimited (to split on)
		HIT = HIT.lstrip()
		HIT=re.sub('\s+','\t',HIT)

		#splitting the current "HIT" on tabs to get the RM data
		(SW_score, 		perc_div,		perc_del,
		 perc_ins, 		query_sequence,		q_begin,
		 q_end, 		q_left,			orient,
		 matching_repeat,	class_family, 		r_begin,
		 r_end, 		r_left, 		ID)	=HIT.split("\t")
		#splitting the next record in the RM data (these will be compared to e/o)
		(SW_score_i,		perc_div_i,		perc_del_i,
		 perc_ins_i,		query_sequence_i,	q_begin_i,
		 q_end_i,		q_left_i,		orient_i,
		 matching_repeat_i,	class_family_i,		r_begin_i,
		 r_end_i,		r_left_i,		ID_i) =HIT_ARRAY[i+1].split("\t")

		#forcing the data into proper types
		q_begin, q_end, q_begin_i, q_end_i=int(q_begin), int(q_end), int(q_begin_i), int(q_end_i)
		perc_div, perc_div_i=float(perc_div), float(perc_div_i)

		# q_begin, q_end, etc.. variables store the start and stop coordinates of the repeat
		# query_sequence is chromosome

		#calulate the lengths of each repeat (and other metrics) to determine which
		# insertion will be saved or how they will be treated
		length_q=q_end-q_begin
		length_q_i=q_end_i-q_begin_i



		CATEGORY = get_category(q_begin,q_begin_i,q_end,q_end_i,query_sequence,query_sequence_i)
#		if not CATEGORY:
#			print "No category detected."
#			return
			
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# CATEGORY 1 OVERLAP
		#
		# chrX   |q_begin----------- q -----------q_end|
		# chrX	          |q_begin_i------------- qi --------q_end_i|

		if CATEGORY=="CATEG1":

			#BASIC STATS
			#Print the number of base pairs in overlap to LENGTHS file
			CAT1_LENGTH=q_end-q_begin_i
			#LENGTHS.write(str(CAT1_LENGTH) + "\n")

			#Print overlap to category 1 file
			#CAT1OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

			if perc_div < perc_div_i:
				q_begin_i=q_end+1

			elif perc_div_i < perc_div:
				q_end=q_begin_i-1

			else:
				if length_q > length_q_i:
					q_begin_i=q_end+1

				elif length_q_i > length_q:
					q_end=q_begin_i-1

				else:
					k=random.randint(0,1)

					if k == 0:
						q_begin_i=q_end+1
					else:
						q_end=q_begin_i-1
			length_q=q_end-q_begin+1
			length_q_i=q_end_i-q_begin_i+1				

			if length_q>1:
				output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence, q_begin, q_end, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID,CATEGORY)
				
#Replace the string below with "NEW_HIT_ARRAY.append(output_string)"
				NEW_HIT_ARRAY.append(output_string)

			if length_q_i>1:
				output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i, q_begin_i, q_end_i, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i,CATEGORY)
				NEW_HIT_ARRAY.append(output_string)


			#delete q_i from future comparisons
			del HIT_ARRAY[i+1]
			

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# CATEGORY 2 OVERLAP
			#
			# chrX                 |q_begin------- q ------q_end|
			# chrX            |q_begin_i----------  qi --------q_end_i|


		elif CATEGORY=="CATEG2":

			CAT2_LENGTH=q_end-q_begin
			#LENGTHS.write(str(CAT2_LENGTH) + "\n")

			#CAT2OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

			if perc_div <= perc_div_i:
				#Write top hit
				output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence, q_begin, q_end, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID,CATEGORY)
				NEW_HIT_ARRAY.append(output_string)
			
				#Write bottom left portion 
				q_begin_i1 = q_begin_i
				q_end_i1 = q_begin-1
			
				output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i,
					q_begin_i1, q_end_i1, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i+"L",CATEGORY)
				NEW_HIT_ARRAY.append(output_string)
			
				#Write bottom right portion
				q_begin_i2 = q_end+1
				q_end_i2 = q_end_i

				output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i,
					q_begin_i2, q_end_i2, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i+"R",CATEGORY)
#				OUTPUT.write(output_string)
				NEW_HIT_ARRAY.append(output_string)

				

			else:
				output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i, q_begin_i, q_end_i, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i,CATEGORY)
#				OUTPUT.write(output_string)
				NEW_HIT_ARRAY.append(output_string)

			del HIT_ARRAY[i+1]
				
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# CATEGORY 3 OVERLAP
			#
			# chrX   |q_begin----------- q -------------------q_end|
			# chrX            |q_begin_i----qi ------q_end_i|

		elif CATEGORY=="CATEG3":

			#Print the number of base pairs in overlap to LENGTHS file	
			CAT3_LENGTH=q_end_i-q_begin_i
			#LENGTHS.write(str(CAT3_LENGTH) + "\n")

			#Print overlap to category 2 file
			#CAT3OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

			if perc_div <= perc_div_i:		
				output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence, q_begin, q_end, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID,CATEGORY)
#				OUTPUT.write(output_string)
				NEW_HIT_ARRAY.append(output_string)
			
			else:
			#Write top left portion	
				q_begin_1=q_begin
				q_end_1=q_begin_i-1
				output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence,
					q_begin_1, q_end_1, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID+"L",CATEGORY)
#				OUTPUT.write(output_string)
				NEW_HIT_ARRAY.append(output_string)
			
			#Write top right portion
				q_begin_2=q_end_i+1
				q_end_2=q_end
				if not q_end == q_end_i:
					output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence,
						q_begin_2, q_end_2, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID+"R",CATEGORY)
#				OUTPUT.write(output_string)			
					NEW_HIT_ARRAY.append(output_string)
			
			#Write bottom line
				output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i, q_begin_i, q_end_i, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i,CATEGORY)
#				OUTPUT.write(output_string)		
				NEW_HIT_ARRAY.append(output_string)
			
			#remove overlap from additional comparisons
			del HIT_ARRAY[i+1]

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# CATEGORY 4 OVERLAP
		#
		# chrX          |q_begin----------- q -----------q_end|
		# chrX   |q_begin_i----- qi -----q_end_i|

		elif CATEGORY=="CATEG4":

			#We need to see if this overlap is significant
			# ex. 2bp vs 2,000 bp overlap
			#OVERLAP=q_end_i-q_begin
			CAT4_LENGTH=q_end_i-q_begin
			#LENGTHS.write(str(CAT4_LENGTH) + "\n")

			#CAT4OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")
		
			#if the overlap is greater than 10 bp (it is significant)
			# and we need identify the longest one
			if perc_div < perc_div_i:
				q_end_i=q_begin-1

			elif perc_div_i < perc_div:
				q_begin=q_end_i+1

			else:
				if length_q > length_q_i:
					q_end_i=q_begin-1

				elif length_q_i > length_q:
					q_begin=q_end_i+1
				
				else:
					k=random.randint(0,1)
					
					if k == 0:
						q_end_i=q_begin-1
					else:
						q_begin=q_end_i+1


			length_q=q_end-q_begin+1
			length_q_i=q_end_i-q_begin_i+1

			if length_q>1:
	#			write_hit()
				output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence, q_begin, q_end, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID,CATEGORY)
#				OUTPUT.write(output_string)
				NEW_HIT_ARRAY.append(output_string)


			if length_q_i>1:
	#			write_hit_i()
				output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i, q_begin_i, q_end_i, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i,CATEGORY)
#				OUTPUT.write(output_string)
				NEW_HIT_ARRAY.append(output_string)


			#delete q_i from future comparisons
			del HIT_ARRAY[i+1]

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# CATEGORY 5 OVERLAP
			#
			# chrX   |q_begin---------- q -----------q_end|
			# chrX   |q_begin_i---------qi --------q_end_i|

		elif CATEGORY=="CATEG5":
	#		print ("Comparing", SW_score, "to", SW_score_i, ".", "cat_5")

			#Print the number of base pairs in overlap to LENGTHS file	
			CAT5_LENGTH=q_end-q_begin
			#LENGTHS.write(str(CAT5_LENGTH) + "\n")

			#Print overlap to category 2 file
			#CAT5OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

			#if cat 5: choose the element most similar to consensus
			if perc_div < perc_div_i:
				output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence, q_begin, q_end, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID,CATEGORY)
#				OUTPUT.write(output_string)
				NEW_HIT_ARRAY.append(output_string)

			
			elif perc_div > perc_div_i:
				output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i, q_begin_i, q_end_i, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i,CATEGORY)
#				OUTPUT.write(output_string) 
				NEW_HIT_ARRAY.append(output_string)

			else:
				k=random.randint(0,1)

				if k == 0:
					output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence, q_begin, q_end, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID,CATEGORY)
#					OUTPUT.write(output_string)
					NEW_HIT_ARRAY.append(output_string)

				else:
					output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i, q_begin_i, q_end_i, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i,CATEGORY)
#					OUTPUT.write(output_string)
					NEW_HIT_ARRAY.append(output_string)

			#remove overlap from additional comparisons
			del HIT_ARRAY[i+1]

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# DO NOT OVERLAP
		#
		else:
		# We are assuming that any insertions not in CAT1-5 are non-overlaps
		#  and printing them out	
	#		print ("Comparing", SW_score, "to", SW_score_i, ".", "no overlap")
	#			OUTPUT.write(HIT)
				output_string=modRMhit(SW_score, perc_div, perc_del, perc_ins, query_sequence, q_begin, q_end, q_left, orient, matching_repeat, class_family, r_begin, r_end, r_left, ID,CATEGORY)
	#			OUTPUT.write(output_string)
				NEW_HIT_ARRAY.append(output_string)


		#We have printed the RM hit, but need to add new line to output file
		# or all would be on the same line
		#OUTPUT.write("\n")

		#Before continuing with the loop, make sure not to continue all the way so
		# that the last element in the RM file is "q" since there will be no "qi"
		#this is minus two because we cant do a comparison with the last entry, and array starts at entry 0
		if i < len(HIT_ARRAY)-2:
			i = i+ 1 
		else:
			if not CATEGORY:
				NEW_HIT_ARRAY.append(HIT_ARRAY[i+1])
#				output_string=modRMhit(SW_score_i, perc_div_i, perc_del_i, perc_ins_i, query_sequence_i, q_begin_i, q_end_i, q_left_i, orient_i, matching_repeat_i, class_family_i, r_begin_i, r_end_i, r_left_i, ID_i)
			break


	#NEW_HIT_ARRAY
#	OUTPUT.write("\n".join(NEW_HIT_ARRAY)+'\n')
	
	return NEW_HIT_ARRAY
	
	
#.........................end of function


if __name__ =="__main__":main()
