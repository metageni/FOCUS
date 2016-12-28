#!/usr/bin/python
import random,os,sys,numpy as np
parameters={"-q":"","-p":"50","-r":"1000"}

#######################################################
#Read and save the parameters given by the user       #
#######################################################
def setParameters():
    if "/" in sys.argv[0]:
        parameters["dir"]="/".join(sys.argv[0].split("/")[:-1])+"/"
    userParameters=sys.argv[1:]
    for i in range(0,len(userParameters),2):
        try:
            parameters[userParameters[i]]=userParameters[i+1]
        except:
            if userParameters[i] in parameters:
                print "Please inform a value for "+userParameters[i]
            else:
                if "-h" not in userParameters:
                    print userParameters[i]+" is not a valid parameter"
setParameters()#store user parameters

INPUT=parameters["-q"]
PERC=int(parameters["-p"])
number_resamples=int(parameters["-r"])


def fasta2hash(fasta):
    c=1;my_hash={}
    for line in "".join(file(fasta).readlines()).split(">")[1:]:
        line=line.split("\n")
        my_hash[c]="".join(line[1:])
        c+=1
    return my_hash

def writeRandsample(ID,keys,sequences):
    o=open("resample/"+str(ID+1)+".fasta","w+")
    for k in keys:
        o.write(">"+str(k)+"\n"+sequences[k]+"\n")
    o.close()
        
#Creates path for bins output
if not os.path.exists("resample"):
    os.makedirs("resample")
    os.makedirs("resample_result")

def main():
    sequences=fasta2hash(INPUT)
    number_randsample=int(len(sequences)*(float(PERC)/100))
    k=range(1,len(sequences)+1)

    print "1) Writing "+str(number_resamples)+" samples with "+str(PERC)+"(%) of sequences"
    for sample in xrange(number_resamples):
        random_sample=random.sample(k,number_randsample)
        writeRandsample(sample,random_sample,sequences)

    print "2) Running FOCUS in the resample sequences"
    os.system("python focus.py -q resample/")

    print "3) Calculating confidence interval"

    #parses the FOCUS output for all the resamples by taking the AVG and STD of all the samples
    for myfile in [i for i in os.listdir(".") if "resample__STAMP_tabular." in i]:
        f=open(myfile)
        o=open("resample_result/"+INPUT+"__"+myfile,"w+")
        head=f.readline().split("\t")
        head="\t".join(head[:len(head)-number_resamples])+"Resample Average\tResample Standard Deviation\n"
        o.write(head)
        for line in f:
            line=line.split("\t")
            
            abundance=[float(x) for x in line[len(line)-number_resamples:]]
            taxa="\t".join(line[:len(line)-number_resamples])
            avg=str(np.average(abundance));std=str(np.std(abundance))
            o.write(taxa+"\t"+avg+"\t"+std+"\n") 
        f.close()
        o.close()
        
if parameters["-q"]!="":
    main()
    os.system("rm *resample__STAMP_tabular.* -r")
    os.system("rm resample/ -r")
    print "4) Done! Please check the folder 'resample_result' for results"
else:
    print "Please set a query [for example: python focus_confidence_interval.py -q sequences.fasta]"
