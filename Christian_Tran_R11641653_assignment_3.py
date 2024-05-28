#!/usr/bin/env python3
# =====================================================================
# Title             : Christian_Tran_R11641653_assignment_3.py
# Description       : This program reads command line arguments for input FASTA files, aligns 2 sequences from score matrix and outputs a FASTA style alignment using Smith-Waterman
# Author            : Christian Tran (R11641653)
# Date              : 10/13/2023
# Usage             : python3 Christian_Tran_R11641653_assignment_3.py -i <input file path for .fna> -o <output file path for .fna> -s <path for score matrix>
# Notes             : This program requires compatibility for Python 3.11.3 and has to accept command line arguments -i and -o respectively for HPCC grading script
# Python Version    : Python 3.11
# =====================================================================

import argparse #import argparse for command line arguments 
from operator import itemgetter #import the itemgetter to specify which element in the list of lists to base the sorting off of

# Function Name : readFastA
# Parameters    : 
#                   inputFile (File type) - Accepts a input file (in this case .fna) to read and sort
# Returns       : list_of_lists[] - mainly used to write the sorted sequences back to another .fna file
# Purpose       : Reads the given inputFile (from the -i (inputpath) command), sorts the sequences based on the length from longest to shortest, then returns the sorted list to main
def readFastA(inputFile):
    list_of_lists = [] #Combined list of all the elements together
    title_sequence = [] #list for the titles of each sequence starting with >
    full_sequence_list = [] #list for the full combined sequence
    sequence = "" #temp string for reading full sequences
    j = 0 #Counter for the sequence element

    #for loop to read the inputFile given and appends them to their respective list depending on if its a sequence title, or part of the actual sequence
    for line in inputFile: #For each line in the file given
        if line[0] == ">": #if the first character in the line starts with a > then:
            title_sequence.append(line.strip()) #Add the title of the sequence to the title_sequence list
            if sequence: #if the sequence HAS something, then add it to the full_sequence_list before getting to the next > symbol
                full_sequence_list.append(sequence)
            sequence = "" #empty the list whenever we reach a new sequence that starts with the character >
        else: #if the first character in a line doesn't start with >, then we can assume that its part of a sequence
            sequence = sequence + line.strip() #until we get to the next sequence title, add lines to make the full sequence
    if sequence: #at the end of the file, if there still is something inside this sequence variable empty it into the full_sequence_list
        full_sequence_list.append(sequence)
        sequence = "" #empty the list
    
    #for loop to create the list of lists to compare and sort the elements
    for element in title_sequence: 
        list_of_lists.append([element, full_sequence_list[j], len(full_sequence_list[j])])#gets the length of the sequence
        j = j + 1 #adds 1 to the counter to iterate through the full_sequence_list
    
    #Sorts the list of lists based off of the length of the sequence string value in list index 2
    list_of_lists.sort(key = itemgetter(2), reverse = True) 

    #Removes the length value from the list of lists
    for h in list_of_lists: 
        h.pop(2) #pops element 2 out of list

    return list_of_lists

#   Function Name: readScoreMatrix
#   Parameteres: 
#               ScoreMatrixFile - file to read score matrix from
#   Returns: A set of tuples for keys (letter pairs) assigned to their corresponding score
#   Purpose: Reads the score matrix file using nested for loops and yield calls to generate a set of tuples to be turned into a dictionary
def readScoreMatrix(ScoreMatrixFile):
  #note: Use .rstrip() (gets rid of all trailing characters) and .split() (makes a list based on the separated elements)
  #List comprehension: A way to apply the given condition to all of the elements in a list. In this case, we want all of the lines in the file
  #                    to have no trailing characters, but also in it's own list, so a list for each row with each score as elements of the list.
  inputFile_lines = (line.rstrip().split() for line in ScoreMatrixFile) #Creates a list of the rows/lines in the file with list comprehension
  
  Score_header = next(inputFile_lines) #iterate the rows once to get the next row, which in this case is the header letters (A, B, C, etc.)

  #for each row in the list of lines from the file that are .rstripped and .split
  for row in inputFile_lines:
    row_letter = row[0] #get the index of the letter
    #Set zip(Score_header, to current whole row) 
    #The zip takes the individual score header (ex: A) and "zips" it in a tuple pairs, so in format (A, 1)
    #Goes through whole list line/row by row and creates tuples of tuples based on the row_letter and col_letter, and goes through all possibilities and assigns them at the current score position
    for col_letter, Score_num in zip(Score_header, row[1:]):
      yield ((row_letter, col_letter), Score_num) #yield a set of tuples as the two letter key values and their score to be used in dict() to turn into dictionary


#   Function Name: smithWaterman
#   Parameteres: 
#               seq1 - The first sequence to align
#               seq2 - The second sequence to align with the first sequence
#               score_dict - The score dictionary from the score matrix file; used to look up matching score values
#   Returns: aligned sequences and the optimal alignment score (score to the bottom right corner) as a tuple
#   Purpose: Implements the smith-waterman algorithm to calculate the optimal score and alignment for 2 given sequences
def smithWaterman(seq1, seq2, score_dict):

    #Note for assignment 3 - Modify this code to do smith waterman implementation instead of needleman wunsch
    
    #initialize empty matrix
    Initial_Matrix = [] 

    #Initialize Traceback matrix
    Traceback_Matrix = []

    #Alignment score variable
    optimal_alignment_score = 0

    #finished alignment variable strings
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    #Get gap penalty from score_dict
    #Note: originally I thought that the gap score came from element key tuple ("-", "-"), but that seems wrong so I am testing it with the letter tuple for the new gap value which might be correct (letter, "-") or ("-", letter)
    #note: getting a dict element from the key is a string, make sure to convert it to an integer for calculation purposes
    #gap_pen = -int(score_dict.get(("-","-"))) #Remember to turn the gap score negative per the algorithm requirment
    #gap_pen = -2 #Test gap score
    #gap_pen = int(score_dict.get()) #Unfinished test gap score
    #Note: changed gap score to (letter and "-") or vice versa, and commented out the ("-", "-") gap score to see if this works.

    #print("CURRENT GAP PENALTY: "+str(gap_pen)) #Test to display current gap penalty

    #Fill empty matrix with zeroes to start off, also fill the traceback matrix with lists and empty strings
    for row in range(len(seq1) + 1): #each row made from the amount of len in seq 1
        Initial_Matrix.append([]) #add empty list inside list to make a list of list (matrix)
        Traceback_Matrix.append([]) #add empty lists to traceback matrix also
        # print(Initial_Matrix)  #used to test and see where the added lists were going
        for col in range(len(seq2) + 1): #each column made from the amount of len in seq 2
          #  print(str(Initial_Matrix[-1])) #used to test and track where the zeroes were going
            Initial_Matrix[-1].append(0) #add zero to each column in each row
            Traceback_Matrix[-1].append("?") #add a character random/question mark into the traceback
            #note: we append to [-1] because we want to add the zeroes to the list we JUST added in the row for loop

    #Fill out the first row and first column with default needleman_wunsch values by taking multiples of the gap score
    #From lecture 05 needleman-wunsch algorithm slides:
    for row_iterator in range(len(seq1) + 1): #Fill out score matrix rows with Zeroes (smith-waterman)
        Initial_Matrix[row_iterator][0] = 0

    for row_iterator in range(len(seq1) + 1): #Fill out traceback matrix rows with the default question mark variable to mark DONE for smith-waterman
        Traceback_Matrix[row_iterator][0] = "?"

    for col_iterator in range(len(seq2) + 1): #Fill out score matrix columns with Zeroes (smith-waterman)
        Initial_Matrix[0][col_iterator] = 0

    for col_iterator in range(len(seq2) + 1): #Fill out traceback matrix rows with the default question mark variable to mark DONE for smith-waterman
        Traceback_Matrix[0][col_iterator] = "?"

    #Note: fill out the Traceback[0,0] with question mark to keep accuracy and to mark where the traceback ends.
    #Note for smith-waterman: Traceback ends whenever the traceback matrix hits an X (or ? in this case) and doesn't include the rest of the letters.
    Traceback_Matrix[0][0] = "?"

    #Test display for initial creation of the matrix
    #print("INITIAL MATRIX SCORE AND TRACEBACK")
    #displayMatrix(Initial_Matrix)
    #print("\n\n\n")
    #displayMatrix(Traceback_Matrix)


    #use the needleman wunsch rules to find the greatest value of the score matrix and fill out all cells
    for row_iterator in range(1, len(seq1) + 1): #Start at range 1 to the end of seq1 length because the first row and column do NOT need to be recalculated due to having default values
        for col_iterator in range(1, len(seq2) + 1): #Start at range 1 because at position (0,0) there is a 0, so theres a small shift in the table by 1 for the default values of the starting score matrix
            #Algorithm pseudocode is on Slides 05 - pairwise alignment slide 26 for Needleman-wunsch rules

            diag_score = Initial_Matrix[row_iterator - 1][col_iterator-1] + int(score_dict.get((seq1[row_iterator-1],seq2[col_iterator-1]))) #add diagonal score from score matrix
            left_score = Initial_Matrix[row_iterator][col_iterator-1] - abs(int(score_dict.get(("-", seq2[col_iterator-1])))) #its initially negative, but to follow the subtract rule, we absolute value it back to positive
            up_score = Initial_Matrix[row_iterator - 1][col_iterator] - abs(int(score_dict.get((seq1[row_iterator-1], "-")))) #same thing with up calculation
            zero_score = 0 #In smith-waterman, you always compare the normal diag, left, and up score with zero. If the greatest score is zero, then put an X (or ? in this case) as it's direction
            #note: I mixed up the left and up scores because I wasn't paying attention. This should be fixed after I changed it (10/10/23)

            max_value = max(diag_score, up_score, left_score, zero_score) #finds the max/largest value out of the 3 scores, compares scores for max one 

            #Note: for smith-waterman I added the comparison to the zero_score as well.

            Initial_Matrix[row_iterator][col_iterator] = max_value #adds the max/largest value out of the 3 and stores it as the main score in this current position

            #Input Traceback directions into traceback matrix
            if(zero_score == max_value):  #Added a zero comparison first to ensure that all zeroes are turned into X (or ? in this case) before making other comparisons to the directions
                Traceback_Matrix[row_iterator][col_iterator] = "?"
            elif(diag_score == max_value):
                Traceback_Matrix[row_iterator][col_iterator] = "diag"
            elif(up_score == max_value):
                Traceback_Matrix[row_iterator][col_iterator] = "up"
            elif(left_score == max_value):
                Traceback_Matrix[row_iterator][col_iterator] = "left"
            

            #Test display for the score tracking
            #print("DIAG SCORE: "+str(diag_score))
            #print("LEFT SCORE: "+str(left_score))
            #print("UP SCORE:   "+str(up_score))
            #print("MAX SCORE FOUND: "+str(max_value))
            #print("ROW: "+str(row_iterator))
            #print("COL: "+str(col_iterator))
            #print("GAP SCORE UP: "+str(score_dict.get(seq1[row_iterator-1], "-")))
            #print("GAP SCORE LEFT: "+str(score_dict.get("-", seq2[col_iterator-1])))
            #print("\n\n\n")

            

    #Test display for both score and traceback matrix:
    #displayMatrix(Initial_Matrix)
    #print("\n\n")
    #displayMatrix(Traceback_Matrix)

    #Note: for smith-waterman, the traceback matrix does NOT start at the bottom right like needleman wunsch, but it starts at the HIGHEST SCORING MATRIX CELL
    #Then it goes through the normal traceback process until it hits an X or in this case, ?, so we need to add code to find the highest scoring matrix cell and
    #start at that position.
    default_max_num = 0 #variable to represent the current max value of the score matrix

    #Goes through each row of the 2d matrix and finds the max of that individual row array, then compares it with the current max value found
    #until it goes through all rows and finds the largest score value.
    for row_iterator in Initial_Matrix:
            default_max_num = max(default_max_num, max(row_iterator)) 

    #Test to see what was the largest score number out of all the score matrix was
    #print("THIS IS THE MAX VALUE FOUND OUT OF ALL THE SCORE MATRIX: "+str(default_max_num))

    #Another nested for loop to find what position the max value is and assign row_iterator and col_iterator to THAT current position
    for row_iterator in range(1, len(seq1) + 1): #Start at range 1 to the end of seq1 length because the first row and column do NOT need to be recalculated due to having default values
        for col_iterator in range(1, len(seq2) + 1): #Start at range 1 because at position (0,0) there is a 0, so theres a small shift in the table by 1 for the default values of the starting score matrix
            if (Initial_Matrix[row_iterator][col_iterator] == default_max_num):
                    break
        else:
            continue #If the inner loop doesn't break, keep going through the rows of the 2d matrix
        break #If the inner loop DOES break, then break out of the nested for loops

    #Test print statement to see where the largest score number position was at
    #print("\n\n")
    #print("THE LARGEST SCORE NUMBER WAS FOUND AT:")
    #print("ROW: "+str(row_iterator))
    #print("COL: "+str(col_iterator))
    #print("\n\n")

    #Changed the optimal alignment score from needleman-wunsch to smith-waterman.
    #Unlike needleman-wunsch, smith-waterman's alignment score is whatever the highest score was found in the 2d score matrix.
    #Since the alignment starts at the highest score in the 2d score matrix, THAT is the optimal alignment score for local alignment.
    #because of this face, smith-waterman's optimal alignment score doesn't always appear in the bottom right like in needleman-wunsch.
    optimal_alignment_score = default_max_num

    #Iterate through the traceback matrix and combine sequence depending on traceback direction
    #Note: Changelog (10/12/23): 
    #I originally used a nested for loop to go through the traceback matrix, which did not work.
    #I went through and saw that I got the score and the traceback matrix right, its just that
    #I wasn't traversing through the traceback matrix properly. Whenever there was a diag, up, or left
    #operation, the for loop would decrement (since I started in the bottom right) along with the 
    #movement operations, adding an extra -1 to the row and column iterator which made the traversal
    #incorrect. I changed it to a while loop to cause the loop to not decrement the row and column iterator
    #like a for loop, and to ONLY CHANGE POSITION WHEN IT CHECKS SPECIFIC DIRECTION. No other thing can
    #cause the row and column iterators to change value except those.

    #Another note: I do NOT have to redefine/reassign col_iterator and row_iterator because they both
    #have a starting value of len(row_iterator) and len(col_iterator) (they both start at the end, bottom right) 
    #due to the previous nested for loop from the score matrix. I am simply just reusing those length values 
    #and going backwards from bottom right to top left. I verified this method works on (10/12/23)

    #Changed for Assignment 3 for smith-waterman. Replaced the DONE check with the ? check to tell the program to stop whenever it encounters a ? (just like how the smith-waterman stops whenever it hits a 0 score or an X)
    while col_iterator > -1 and row_iterator > -1:
        #print("INNER LOOP: Curr row: "+str(row_iterator))
        #print("INNER LOOP: Curr col: "+str(col_iterator))
        #Start at bottom right and go to DONE 
        if(Traceback_Matrix[row_iterator][col_iterator] == "?"): #Check if reached any X (in this case, ? symbol), then exit nested for loops and stop aligning.
           # print("DONE DETECTED")
            break #Exit inner for loop
                
        elif(Traceback_Matrix[row_iterator][col_iterator] == "diag"): #if the current traceback is diagonal, then align both sequence letters together and move diagonal
            aligned_seq1 += seq1[row_iterator-1] #align seq1 letter with sequence 2
            aligned_seq2 += seq2[col_iterator-1]
            row_iterator -= 1 #move diagonally backwards since we started at the bottom right and moving to the upper left
            col_iterator -= 1
          #  print("DIAG CHANGE: Curr row: "+str(row_iterator))
          #  print("DIAG CHANGE: Curr col: "+str(col_iterator))
          #  print("diag DETECTED")
            continue

        elif(Traceback_Matrix[row_iterator][col_iterator] == "up"): #if the current traceback is up, then align just seq1 and gap the up letter from seq2      
            aligned_seq1 += seq1[row_iterator-1] #align just seq1 current letter
            aligned_seq2 += "-" #gap the up from seq2
            row_iterator -= 1 #move up on the traceback matrix
            #col_iterator += 1 #stay on the same column
          #  print("UP CHANGE: Curr row: "+str(row_iterator))
          #  print("UP CHANGE: Curr col: "+str(col_iterator))
          #  print("up DETECTED")
            continue
               
        elif(Traceback_Matrix[row_iterator][col_iterator] =="left"): #if the current traceback is left, then align just seq2 and gap the left letter from seq1    
            aligned_seq1 += "-" #gap the left from seq1
            aligned_seq2 += seq2[col_iterator-1] #align just seq2 current letter
            col_iterator -= 1 #move left on the traceback matrix
            #row_iterator += 1 #stay on same row
           # print("LEFT CHANGE: Curr row: "+str(row_iterator))
           # print("LEFT CHANGE: Curr col: "+str(col_iterator))
           # print("left DETECTED")
            continue
              
    #reverse alignment strings to get the correct alignment since it aligns starting at the bottom right like we do when we're drawing it out, from right to left
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return (aligned_seq1, aligned_seq2, optimal_alignment_score) #return finished aligned sequences as a tuple



#   Function Name: displayMatrix
#   Parameteres: 
#               matrix - The input matrix to be displayed
#   Returns: nothing
#   Purpose: Used mainly for test purposes, this operation displays a given list of lists in matrix style
def displayMatrix(matrix):
    for i in matrix:
        print(str(i) )


def main():
    parser = argparse.ArgumentParser() #Establish the parser
    sorted_sequences = [] #List to represent the finished, sorted list_of_lists from the readFastA function
    score_dict = {} #List to represent the scores from each letter in the dictionary from the score matrix BLOSUM50 or BLOSUM62
    # Adding an argument/option to the command line
    # -i is used to specify the file path to the input FASTA file
    # called in this format:
    # python3 <python file name> -i <input file path>
    parser.add_argument("-i", "--inputpath", help = "This option retrieves the file path to the starting cellular matrix", type = argparse.FileType('r'))
    
    # -o is used to specify the file path to the output file
    # I think that he wants us to have the capability to write out the output to an output file
    # called in similar format to our -i:
    # python3 <python file name> -o <output file path> 
    parser.add_argument("-o", "--outputpath", help = "This option retrieves the file path for the final output file", type = argparse.FileType('w'))
    #Note: I verified that option -o specifies the and writes to a new file with the output, but you can specify the filepath to where the output goes. Requires this code to be ran in the HPCC grading script (verified on 9/18/2023)
    

    # -s is used to specify the file path to the score matrix
    # calledin format:
    # python3 <python file name> -s <file path to the score matrix>
    parser.add_argument("-s", "--scorepath", help = "This operation retrieves the file path for the score matrix", type = argparse.FileType('r'))

    args = parser.parse_args() #runs parser and parses the arguments from the command line input

    #Calls the readFastA function to read all sequences, sort them, then return the sorted list of lists
    sorted_sequences = readFastA(args.inputpath.readlines()) 

    
    #Reminder:
    #The first line of the output printed to the screen during execution MUST be in the format:
    # "Assignmetn 1 :: R#"
    #Where R# is the TTU R#. 
    #This is REQUIRED for the grading script
    #Displays the required 1st line printed output fromt he program output specifications:
    print("Assignment 1 :: R11641653")

    SCORE_FILE = open(args.scorepath.name) #Opens the Score Matrix file to be used in readScoreMatrix

    score_dict = dict(readScoreMatrix(SCORE_FILE)) #Calls the readScoreMatrix and generates a dictionary based on the yield tuples of letter pairs WITH their given matching score
    
    
    sequence_1 = sorted_sequences[0][1] #Gets the second sequence from FASTA 

    sequence_2 = sorted_sequences[1][1] #Gets the first sequence from FASTA


    #Test display to check for sequences
    #print("INPUT SEQUENCES: \n")
    #print("INPUT 1: "+sequence_1+"\n")
    #print("INPUT 2: "+sequence_2+"\n")
    #print(sorted_sequences)
    

    #run smith-waterman algorithm and align the two given sequences with the given score matrix, then return a tuple of the aligned matricies
    #note: I changed from needleman-wunsch to smith-waterman for assignment 3. It was exactly the same as needleman-wunsch except for a couple of differences/changes:
    # 1. Smith waterman's optimal alignment score isn't in the bottom right corner like needleman-wunsch, but at whatever position the MAXIMUM score found was in the entire
    #   2d score matrix.
    # 2. Smith waterman's alignment doesn't start at the lower right corner like needleman-wunsch, but it starts at whereever the MAXIMUM score found was in the score matrix.
    # 3. Instead of the filling the initial matrix with up, left, DONE, and multiples of the gap score, smith-waterman fills the initial score matrix with 0 and X (I used ? symbols in this case)
    # 4. Instead of going through alignment until the code finds DONE, the alignment will STOP whenever it hits a ? symbol (like the X symbol in the hand-written examples)
    aligned_sequences = smithWaterman(sequence_1, sequence_2, score_dict)

    #note: sequence_2 goes first because the "readFastA" sorts the 

    #Test display for tuple output check
    #print(aligned_sequences)
    #print(aligned_sequences[0])
    #print(aligned_sequences[1])

    #Take the output tuple of "aligned_sequences" and reassign them back into the list with appended score values and new aligned sequences
    sorted_sequences[0][1] = aligned_sequences[0] #Reassign the sequence 1 array with the aligned version of it back to position 0 since the algorithm always puts sequence/longest sequence at position 0
    sorted_sequences[1][1] = aligned_sequences[1] #Reassign the sequence 2 array with the aligned version of it back to position 1 since the algorithm has the shortest sequence at position 1
    sorted_sequences[1][0] = sorted_sequences[1][0]+"; score="+str(aligned_sequences[2]) #append the optimal alignment score to the end of the their sequence headers
    sorted_sequences[0][0] = sorted_sequences[0][0]+"; score="+str(aligned_sequences[2]) #append the optimal alignment score to the end of the their sequence headers
    
    #Test display to see sorted_sequences
    #print(sorted_sequences)

    #Writes out the output to the specified output location and file name from -o argument from command line, adding a new line for each sequence
    #In this case, it takes the "sorted_sequences" with the newly assigned
    for elements in sorted_sequences:
        for h in elements:
            args.outputpath.write(h + "\n")

#ASSIGNMENT 3 FINISHED ON 10/15/23, make sure to test it on the HPCC before turning it in 
#Runs the Main() function
if __name__ == '__main__': 
    main()