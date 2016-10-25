#!/usr/bin/python

#xipe_comparison.py 0.1
# This program uses FOCUS output (http://edwards.sdsu.edu/FOCUS) into XIPE (http://edwards.sdsu.edu/cgi-bin/xipe.cgi)
# For a a non-parametric statistical analysis of the distribution of samples to determine which samples are statistically
# significantly different.
# The script parses FOCUS output and compare the each pair of samples and reports which targeted level (Genus for example)
# were statistically different in a confidence level

import numpy as np, random,os,sys

parameters={"-q":"","-l":"Genus","-c":"95","-o":"my_xipe_output.xls"}

def setParameters():
    r=0
    if "/" in sys.argv[0]:
        parameters["dir"]="/".join(sys.argv[0].split("/")[:-1])+"/"
    userParameters=sys.argv[1:]
    if "-h" in userParameters:
        parameters["-h"]=1
        print "xipe_comparison.py: Uses XIPE for a non-parametric statistical comparision of the samples in the FOCUS ouput"
        print "     -q FOCUS output (*__STAMP_tabular.spf)"
        print "     -c Minimum Confidence Level (Default: 95)"
        print "     -l Comparison Level [kingdom,phylum,class,order,family,genus,species,strain,all](Default: genus)"
        print "     -o Output Name (Default: my_xipe_output.xls)\n"
        print "     example> python xipe_comparison.py -q input__STAMP_tabular.spf.xls -c 95 -l genus -o xipe__genus___sharks_ouput"
    else:
        for i in range(0,len(userParameters),2):
            try:
                parameters[userParameters[i]]=userParameters[i+1]
            except:
                if userParameters[i] in parameters:
                    print "Please inform a value for "+userParameters[i]
                else:
                    print userParameters[i]+" is not a valid parameter"
                    
setParameters()#store user parameters

INPUT=parameters["-q"]
LEVEL=parameters["-l"].lower()
confidence_level=int(parameters["-c"])
OUPUTNAME=parameters["-o"]

def focus2xipe():
    global INPUT,LEVEL,confidence_level,OUPUTNAME
    f=open(INPUT)

    #head with all the levels + file names
    head=f.readline().replace("\n","").split("\t")
    #only the levels
    #lower case all the levels names to help in the identification for users who messes up with the level name
    levels=[x.lower() for x in head[:8]]
    #only yhe filenames
    filenames=head[8:]
    #Index for the level chosen by the user
    I=levels.index(LEVEL.lower())

    #reduces the big table up to the level showed by the user
    h={}
    for line in f:
        line=line.replace("\n","").split("\t")

        target=line[I]#level
        #Abundanc for comparison. Xipe only accepts interget, so i round it to 4 and multiple the number by 10000
        info=[int(round(float(xx),4)*10000) for xx in line[8:]]

        #add to the hash to write it later
        if target not in h:
            h[target]=[]
        h[target]+=[info]
    f.close()

    #info of each file into a matrix
    m=[]
    for i in h:
        info=[str(x) for x in np.sum(h[i], axis=0)]
        info=[i]+info
        m.append(info)
    m=np.array(m).T



    #parse the xipe result and return the list of targets that passed by the confidence_level
    def parse_xipe():
        f=open("xipe_out.txt")

        result=[]
        for line in f:
            for j in ["(",")","\n"," "]:
                line=line.replace(j,"")
            line=line.replace(",","\t").split("\t")

            if int(line[1])>=confidence_level:
                result.append(line)
        f.close()
        return result


    temp=[]
    k=range(1,len(m))

    #Does all the combination of files
    xipeout=open(LEVEL+"__"+OUPUTNAME,"w+")
    xipeout.write("Parameters:\n")
    xipeout.write("Input: "+INPUT+"\n")
    xipeout.write("Level: "+LEVEL+"\n")
    xipeout.write("Min. Confidence Level: "+str(confidence_level)+"\n\n")
    xipeout.write(LEVEL+"\tCategory 1 over-represented\tCompared to Category 2\tConfidence Level\n")

    print "\n--> Level Comparison: "+LEVEL
    for i in k:
        for j in k:
            t=[i,j];t.sort()

            #found a pair that was not check before
            if i!=j and t not in temp:
                temp.append(t)
                filekeys={'1':filenames[i-1],'2':filenames[j-1]}

                print "Comparing "+filekeys["1"]+" vs "+filekeys["2"]

                #write file 1
                o=open("file1.txt","w+")
                for ii in np.array([m[0],m[i]]).T:
                    o.write("\t".join(ii)+"\n")
                o.close()

                #write file 2 
                o=open("file2.txt","w+")
                for ii in np.array([m[0],m[j]]).T:
                    o.write("\t".join(ii)+"\n")
                o.close()
                #run xipe for pair of files
                os.system("python xipe.py --file1 file1.txt --file2 file2.txt --outfile xipe_out.txt --nrepetitions 1000 --samplesize 1000")
                xipeouput=parse_xipe()

                #write final output
                for xipe_target in xipeouput:
                    over=filekeys[xipe_target[2]]
                    under=[xx for xx in filekeys.values() if xx != over][0]
                    xipeout.write(xipe_target[0]+"\t"+over+"\t"+under+"\t"+xipe_target[1]+"\n")
                xipeout.write("\n")
    xipeout.close()

    print "\nDone :)"
    print "Please check "+LEVEL+"__"+OUPUTNAME
    os.system("rm xipe_out.txt file1.txt file2.txt")



def main():
    global INPUT,LEVEL,confidence_level,OUPUTNAME
    #check if focus output exists

    if os.path.isfile(parameters["-q"]):
        l=["kingdom","phylum","class","order","family","genus","species","strain","all"]
        if LEVEL in l:
            if LEVEL == "all":
                for LEVEL in l[:-1]:
                    focus2xipe()
            else:
                focus2xipe()
        else:
            print LEVEL+"  is invalid! Please enter on of these:"
            for i in l:
                print i
            print "\n"
    else:
        print "file in [-q] doesn't exist!!!"

if "-h" not in parameters:
    print "Running xipe_comparison.py"
    main()
    
