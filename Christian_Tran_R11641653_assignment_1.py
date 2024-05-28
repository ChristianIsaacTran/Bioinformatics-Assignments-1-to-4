import argparse #import argparse for command line arguments 
from operator import itemgetter

#Reads the given FASTA file and then puts elements into a list of lists for each sequence
def readFastA(inputFile):
    list_of_lists = [] #Combined list of all the elements together
    title_sequence = [] #list for the titles of each sequence starting with >
    sequence = "" #temp string for reading full sequences
    full_sequence_list = [] #list for the full combined sequence
    j = 0 #Counter for the sequence element

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

    for element in title_sequence: #for loop to create the list of lists to compare and sort the elements
        list_of_lists.append([element, full_sequence_list[j], len(full_sequence_list[j])])#gets the length of the sequence
        j = j + 1 #adds 1 to the counter to iterate through the full_sequence_list
           
    list_of_lists.sort(key = itemgetter(2), reverse = True) #Sorts the list of lists based off of the length of the sequence string value in list index 2

    for h in list_of_lists: #Removes the length value from the list of lists
        h.pop(2) #pops element 2 out of list

   # for i in list_of_lists: #display sorted list
    #    print(str(i))

    return list_of_lists



        
        
    



def Main():
    parser = argparse.ArgumentParser() #Establish the parser
    sorted_sequences = []
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
    #Note: I verified that option -o specifies the and writes to a new file with the output, but you can specify the filepath to where the output goes
    args = parser.parse_args() #runs parser and parses the arguments from the command line input

    #Calls the readFastA function to read all sequences, sort them, then return the sorted list of lists
    sorted_sequences = readFastA(args.inputpath.readlines()) 

    for elements in sorted_sequences:
        for h in elements:
            args.outputpath.write(h + "\n")

    #print(args.inputpath.readlines()) #basic test to read in all elements in the FASTA file
    #args.outputpath.write("Bruh moment!") #basic test to write the output to a output .txt file 

    #Note: Go to class today and verify if the -o option is supposed to write to a output file or should specify the filepath


#Runs the Main() function
if __name__ == '__main__': 
    Main()