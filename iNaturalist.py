import csv
import os
import argparse
import math

#Parse arguments from the terminal
def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help = "The name of the given file")
    args = parser.parse_args()

    return args.name

#Finds the file in a specified directory
def finder(name):
    for root, dir, files in os.walk('/mnt/c/Users/AGATA/Downloads'):
        for f in files:
            if name==f:
                return os.path.join(root,name)

#Opens csv file and returns the labels of columns and columns without the labels
def opener(path):
    with open(path) as op:
        csv_file = csv.reader(op, delimiter=",")
        csv_lines = [row for row in csv_file]
        return csv_lines[0], csv_lines[1:]

#Selects unrepeated taxons
def unrepeated(file):
    taxons = []
    for row in lines:
        if row[5] not in taxons:
            taxons.append(row[5])
    return taxons

#Calculates the Shannon-Wiener and Simpson indexes. Requires the list of the unrepeated names of the species and the original list with information
def shannon_wiener_simpson(unrep_taxon, orig_lines):
    track = unrep_taxon.copy() #Copy of the unrepeated taxon names
    all_specimens = [row[5] for row in orig_lines] #Repeated names of the organisms
    sum_shw = 0 #The sum for Shannon-Wiener index (eventually negative)
    sum_simpson = 0 #The sum for Simpson index
    for species in range(len(unrep_taxon)): #Sums up to the number of the unique names of the organisms
        for row in orig_lines: #Checks every observation in the orginal file containing possible duplicates
            if row[5] in track: #If the species has not been counted
                num_tax_one = all_specimens.count(row[5]) #The number of the occurences of the specimens of the particular species among all specimens
                fraction = float(num_tax_one)/float(len(orig_lines)) #The proportion of the number of specimen of the specific species among all specimens
                sum_el = fraction*math.log(fraction) #The element of the sum
                sum_shw = sum_shw-sum_el #Substraction gives the negative value of the overall expression
                sum_simpson = sum_simpson+(num_tax_one*(num_tax_one-1))
                track.remove(row[5]) #Prevents from repeating the taxon names
    return sum_shw, float(sum_simpson)/float(len(orig_lines)*(len(orig_lines)-1))






if __name__ == "__main__":
    id  = arguments()
    path = finder(id)
    columns, lines = opener(path)
    tax = unrepeated(lines)
    shw, simpson = shannon_wiener_simpson(tax, lines)
    print("Liczba obserwacji: ", len(lines), "\nBogactwo gatunkowe: ", len(tax), "\nWspółczynnik Shannona-Wienera: ", shw,
          "\nWspółczynnik Simpsona: ", simpson)

