import csv
from ftplib import FTP
#Download sequences from NCBI
from Bio import Entrez
from Bio import SeqIO
#Using records
from Bio.Blast import NCBIWWW
import tkFileDialog
import tkMessageBox
import winsound
from Bio.Blast import NCBIXML
from StringIO import StringIO
import numpy as np
import panda as pd
from Tkinter import *
import string
import os

#update data, create globa table with name, id, type, of sequences
def update_files():
    staticRows=4  #static rows in global table
    numLines=0  #lines in global table
    ftp = FTP('ftp.ncbi.nlm.nih.gov')  #FTP site
    ftp.login(user='anonymous', passwd = 'anonymous')
    ftp.cwd('/genomes/GENOME_REPORTS/')  #specific file in site
    save_path=os.path.dirname(os.path.abspath(__file__))+r"\databases"  #save the ftp files to databases folder
    if not os.path.exists(save_path):  #if folder not exist, create it
        os.makedirs(save_path)
    ID_data=[]  #list with all id in that appear in glaobal table, help to know which sequences is not in table yet
    len_row=0  #len of coulomn in row, help to know how many coulomns added for blast
    if os.path.isfile(name_table) :  #if the table exist, fill list with id that exist in table
        f=open(name_table, 'rb')
        _table = csv.reader(f)
        first_raw=0
        for line in _table:
            if first_raw==0:  #name of coulomn
                len_row=len(line)
                first_raw=1
            ID_data+=[line[1]]
        f.close()
        numLines=len(ID_data)
    else:  #table not exist
        with open(name_table, 'wb') as myfile:  #create file and add titles
            global_table = csv.writer(myfile)
            global_table.writerow(["NUMBER", "NAME", "ID", "TYPE"])
    with open(name_table, 'ab') as myfile:  #open in mode 'ab' =add in the end of the file
        global_table = csv.writer(myfile)
        for name in names_of_DB:  #for each in list of DB
            filename = name+'.txt'
            print "Updating ",filename
            completeName = os.path.join(save_path, filename)  #path to save in the file
            localfile = open(completeName, 'wb')
            ftp.retrbinary('RETR ' + filename, localfile.write, 1024)  #download
          #Conver to csv
            txt_file = completeName
            filename=name+".csv"
            csv_file =  os.path.join(save_path, filename)
            in_txt = csv.reader(open(txt_file, "rb"), delimiter = '\t')
            out_csv = csv.writer(open(csv_file, 'wb'))
            out_csv.writerows(in_txt)
            with open(csv_file, 'rb') as handle:  #open XML
                my_csv = csv.reader(handle)
                first_line=0
                for line in my_csv:  #add to table if not there yet
                    if first_line!=0:
                        if len(line)>coloumn[name] and not line[coloumn[name]] in ID_data :
                            numLines+=1
                            ID_data+=[line[coloumn[name]]]
                            global_table.writerow([numLines, line[0], line[coloumn[name]], name]+ ["?"]*(len_row-staticRows))
                    first_line=1  #dont add first line==titles
            print "updated!"
    print("All Updated, global table updated")

#download sequence from ncbi
def download_seq():
    Entrez.email = "A.N.Other@example.com"  #change to real address to get messages (about problems, warnning)
    
    root = Tk()
    root.attributes("-topmost", True)  #pop in front other windows
    root.wm_title("DNA or Protein?")  #title
    Label(root, text="Blast DNA or protein?  Choose and close box", justify=LEFT, padx=20).pack()
    v = IntVar()
    v.set(2)
    Radiobutton(root, text="Nucleotids", variable=v, value=1).pack(anchor=W)  #radio button, allow user choose nucleotide or protein
    Radiobutton(root, text="Protein", variable=v, value=2).pack(anchor=W)
    root.mainloop()

    p=v.get()
    if p==1:
        DataBase="nucleotide"
    else:
        DataBase="protein"

    def close_window():  #internal function, activate by button "OK"
        global SequenceID, SequenceName
        SequenceID = e1.get()  #input to variables the user inputs
        SequenceName = e2.get()
        master.destroy()

    master = Tk()  #allow user input sequence id and file to save in the file
    master.wm_title("sequence information")  #title
    Label(master, text="sequence ID").grid(row=0)
    Label(master, text="sequence file name").grid(row=1)
    e1 = Entry(master)
    e2 = Entry(master)
    e1.grid(row=0, column=1)
    e2.grid(row=1, column=1)
    B = Button(master, text = "OK", command = close_window).grid(row=2)  #button to close the input window
    master.mainloop()

    handle = Entrez.efetch(db=DataBase, rettype="gb", retmode="text", id=SequenceID)  #download through entrz the sequence
    Sequence = SeqIO.read(handle, "gb") #using "gb" as an alias for "genbank"
    handle.close()

    save_path=os.path.dirname(os.path.abspath(__file__))+r"\results\\"+SequenceName  #path to save the file in results
    if not os.path.exists(save_path):  #if folder not exist, create it
        os.makedirs(save_path)
    completeName = os.path.join(save_path, SequenceName+".gb") 
    
    save_file = open(completeName, "w")  #open for write
    SeqIO.write(Sequence,completeName,"gb")
    save_file.close()

    print SequenceName," was downloaded"

#Open a gb file with DialogBox
def open_gb_fasta(BlastFileType,ftype):
    root=Tk()
    save_path=os.path.dirname(os.path.abspath(__file__))+r"\results"
    root.Fullname = tkFileDialog.askopenfilename(title = "choose your file",initialdir =
                (save_path),filetypes =
                (ftype,("all","*.*")))
    Fullname=root.Fullname
    root.destroy()
    print Fullname


    handle = open(Fullname, "rU")
    temp=(Fullname.split("/")[-1]).split(".")
    SequenceName=temp[0]
    Ftype="."+temp[1]

    save_path=os.path.dirname(os.path.abspath(__file__))+r"\results\\"+SequenceName
    if not os.path.exists(save_path):
        handle.close()
        os.makedirs(save_path)
        completeName = os.path.join(save_path, SequenceName+Ftype)
        save_file = open(completeName, "w")
        save_file.write(open(Fullname).read())
        save_file.close()
        Fullname=completeName
        
        
    global Save_path  #define the folder we work in this run
    Save_path=os.path.dirname(Fullname)
    handle = open(Fullname, "rU")
    seq1= list(SeqIO.parse(handle, BlastFileType))
    handle.close()

    # take nucleotide seq
    Seq_nuc=seq1[0].seq
    return Seq_nuc

#blast, allow dowmload/ open exist sequence and choose mode nuc/protein
def blast():
    root = Tk()
    root.attributes("-topmost", True)  #pop in front all windows
    root.wm_title("Files")  #title
    Label(root, text="Download file", justify=LEFT, padx=20).pack()
    v = IntVar()
    v.set(2)
    #radiobutton to allow user chioce for souece to blast
    Radiobutton(root, text="download gb file", variable=v, value=1).pack(anchor=W) 
    Radiobutton(root, text="Open local gb file", variable=v, value=2).pack(anchor=W)
    Radiobutton(root, text="Open local fasta file", variable=v, value=3).pack(anchor=W)
    root.mainloop()

    query = ''  #get query according to user choice
    p = v.get()
    if p==1:
        download_seq()
    if p==2 or p==1:
        BlastFileType = "gb"
        ftype = "GeneBank","*.gb"
        query=open_gb_fasta(BlastFileType,ftype)
    if p==3:
        BlastFileType = "fasta"
        ftype = "Fasta", "*.fasta"
        query = open_gb_fasta(BlastFileType,ftype)

    print "Now BLASTing"

    #BLAST
        # DNA or protein
    root = Tk()
    root.attributes("-topmost", True)  #pop in front other windows
    root.wm_title("DNA or protein?")  #title
    Label(root, text="Blast DNA or protein?", justify=LEFT, padx=20).pack()
    v = IntVar()
    v.set(2)
    #choice nuc/protein
    Radiobutton(root, text="Nucleotids", variable=v, value=1).pack(anchor=W)
    Radiobutton(root, text="Protein", variable=v, value=2).pack(anchor=W)
    root.mainloop()

    p=v.get()  #according the choice select the program for blast
    if p==1:
        Program="blastn"
    else:
       Program="tblastn"

    #get the info from file- definitions
    info=open("Definitions.txt", "r").read()
    info=info.split()
    Dbase = info[info.index('Dbase')+1]
    Eval = float(info[info.index('Eval')+1])
    HitList = int(info[info.index('HitList')+1])
    Align = int(info[info.index('Align')+1])
    Description = int(info[info.index('Description')+1])
    Ftype = info[info.index('Ftype')+1]

    #running blast
    BlastResults = NCBIWWW.qblast(Program,Dbase,query,expect=Eval,hitlist_size=HitList,alignments=Align,
                                  descriptions=Description,format_type=Ftype)

    #Save the results
    SaveFileName = "my_blast." + Ftype
    print "Results in ",SaveFileName
    completeName = os.path.join(Save_path, SaveFileName)  #path to save in the results
    save_file = open(completeName, "w")
    save_file.write(BlastResults.read()) #save results to file "my_blast"
    save_file.close()
    BlastResults.close()

    winsound.Beep(300,2000)  #Alert Sound

#add blast results to global table, for each sequence mark yes/no 
def Add_global_table(positive, negative):
    line_num=0
    with open(name_table, 'rb') as myfile:  #open file of global table for read
        global_table = csv.reader(myfile)
        with open("temp_table.csv", 'wb') as handle:  #open temp file for write
            temp_table = csv.writer(handle)
            for line in global_table: #each line in table
                if line_num==0:  #if first line , add title= name of the sequence of blast
                    line_num=1
                    temp_list=Save_path.split('/')
                    name=temp_list[temp_list.index('results')+1]
                    temp_table.writerow(line+[name])
                    id_coulomn=line.index("ID")
                elif line[id_coulomn] in positive:  #if id in positive add yes
                    temp_table.writerow(line+["yes"])
                elif line[id_coulomn] in negative:  #if id in negative add no
                    temp_table.writerow(line+["no"])

                else:
                    temp_table.writerow(line+["?"])  #if not in positive or negative
    os.remove('global_table.csv')  #delete old table
    os.rename('temp_table.csv', 'global_table.csv')  #rename the table

#analyze blast results
def analyze():

    Fullname="my_blast.XML"
    print "Now Uploading Blast results in XML file: ",Fullname
    completeName = os.path.join(Save_path, Fullname)
    result_handle=open(completeName)  #open results
    blast_results=NCBIXML.read(result_handle)

    #Convert from XML to List
    BlastList=[]  #list hold information about the match seq
    BlastArray=[]  #list hold ID of match seq
    for alignment in blast_results.alignments:
        BlastList.append(alignment.title)
        BlastArray+=[alignment.accession]
    BlastList2 = BlastList
    BlastList = [BlastList.replace('>','\n') for BlastList in BlastList2]

    #Write the List to a file
    Fullname ="BLASToutput.txt"
    completeName = os.path.join(Save_path, Fullname)
    f = open(completeName,"w") #W for writing
    for item in BlastList:
      f.write("%s\n" % item)
    f.close()
    print "Blast results are in ",Fullname

    print "There are ", len(BlastArray), " BLAST results"

    #checkbox for user
    my_var={}
    root=Tk()
    root.attributes("-topmost", True) #window in front
    Label(root, text="Choose which DBs to compare with and close the box").grid(row=0, sticky=W)
    root.wm_title("Choose Dbs")  #user choose DBs to analyze with
    for i in range(len(names_of_DB)):
        my_var[i] = IntVar()
        my_var[i].set(1)
        Checkbutton(root, text=names_of_DB[i], variable=my_var[i]).grid(row=i+1, sticky=W)
    root.mainloop()
    
    IDblast=set(BlastArray)  #make set of IDs

    #hold id of positive / negative to help in fill of global table 
    All_positiveID=[]
    All_negativeID=[]
    count_positive=0
    count_negative=0

    root=Tk()  #ask if user want to add results to global table
    root.withdraw()
    root.AddOrNot=tkMessageBox.askquestion("Gloabal table", "Do you want to add the results to global table?")
    AddOrNot=root.AddOrNot
    root.destroy()
    #compare the result of blast to DB
    for i in range(len(names_of_DB)):
        save_path = os.path.dirname(os.path.abspath(__file__)) + r"\databases"  #path to DB
        Positives=[]  #hold the info in pos/ neg
        Negatives =[]
        if my_var[i].get():  #if user chose it
            name=names_of_DB[i]
            Fullname=name+".csv"
            completeName = os.path.join(save_path, Fullname)  #path of the DB file
            print "Working on ",Fullname
            mone=0
            with open(completeName, "r") as infile:
                FileLength = sum(1 for line in infile)
                infile.seek(0) #Return to begining of file
                CSVreader = csv.reader(infile)
                for line in CSVreader:
                    mone+=1
                    if mone==1:  #add title
                        Positives.append(line)
                        Negatives.append(line)
                    else:
                        if not mone % 1000:  #help to see the progress
                            print mone, " of ", FileLength, " of ", name
                        #if line is full(not part of line) and if id in set of blast results=>add to positive list
                        if len(line)>coloumn[name] and [Id for Id in IDblast if Id in line[coloumn[name]]]:
                            Positives.append(line)
                            if AddOrNot == "yes":  #if needed to add to global table
                                All_positiveID+=[line[coloumn[name]]]
                        elif len(line)>coloumn[name]:  #add to negative
                            Negatives.append(line)
                            if AddOrNot == "yes":  #if needed to add to global table
                                All_negativeID+=[line[coloumn[name]]]
            print "So far there are in "+name+":"
            print "Positives ",len(Positives)-1
            print "Negatives ",len(Negatives)-1
            count_positive+=len(Positives)-1
            count_negative+=len(Negatives)-1
            print "So far there are in general:"
            print "Positives ",count_positive
            print "Negatives ",count_negative  

            #write positive and negative to files
            completeName = os.path.join(Save_path, "Positives_"+name+".csv") 
            with open(completeName, 'wb') as PositiveFile:
                wr = csv.writer(PositiveFile, dialect='excel')
                wr.writerows(Positives)
            completeName = os.path.join(Save_path, "Negative_"+name+".csv") 
            with open(completeName, 'wb') as NegativeFile:
                wr = csv.writer(NegativeFile, dialect='excel')
                wr.writerows(Negatives)
        
    if root.AddOrNot == "yes":  #update the results in global table
        Add_global_table(All_positiveID, All_negativeID)

    winsound.Beep(300,2000)  #Alert Sound


#-------------------------------------main----------------------------------------#

names_of_DB=['viruses','plasmids', 'prokaryotes', 'eukaryotes']  #list for use in some functions
coloumn={'plasmids':5, 'eukaryotes':9, 'prokaryotes':8, 'viruses':9}  #dictionary, the value is the coulomn in ftp file there is the id (genebank)
global Save_path  #represent path to file we work on
global name_table  #table that save global information
name_table="global_table.csv"  
my_root=Tk()  #open root
my_root.withdraw()  #hide the black (Tk) window 
my_root.UpdateOrNot=tkMessageBox.askquestion("Updat db", "Do you want to updated the databases?")
if my_root.UpdateOrNot == "yes":
    print "Now Updating db"
    print "Please wait..."
    update_files()  #update files and global table
my_root.BlastOrNot=tkMessageBox.askquestion("Blast", "Do you want to Blast a sequence?")
BlastOrNot=my_root.BlastOrNot
my_root.destroy()  #destroy root
if BlastOrNot == "yes":
    blast()  #blast
else:
    root=Tk()
    root.withdraw()
    root.attributes("-topmost", True)  #pop window (Tk) in front all other windows
    #in case need to do analyze on exist file, ask for the spesific file
    save_path=os.path.dirname(os.path.abspath(__file__))+r"\results"  
    Save_path=tkFileDialog.askdirectory(parent=root, initialdir=save_path,
                                    title='Please select the folder to analyze')
    print Save_path
    root.destroy()
print "Okay, lets analyze.."
analyze()  #analyze the blast results
