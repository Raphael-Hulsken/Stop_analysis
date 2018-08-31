import threading
import Queue

import time, random, os
from os import listdir
from os.path import isfile, join
import subprocess

WORKERS = 16


############################################################################
# The goal is to run combine tool on cards which are stored in a directory
# It uses a //ation and the nof worker can be tuned
# As a reference  - running on 350 points with 12 nodes takes 7'30 
# with 15 boxes and only one systematics
# The best would be to used the largest nof workers: 16 ?
############################################################################



executable = "combine -M ProfileLikelihood --significance -t -1 --expectSignal=1"
directory = "cards/"
extension = ".tab"
ofilename = "significance.txt"
counter = 0

class Worker(threading.Thread):

    def __init__(self, queue):
        self.__queue = queue
        threading.Thread.__init__(self)

    def run(self):
        while 1:
            item = self.__queue.get()
            if item is None:
                break # reached end of queue

	    
	    # here we launch our job
	    commandline = executable + " " + directory + "/" + item
	    #p = subprocess.Popen([executable,directory+"/"+item], shell=True, stdout=subprocess.PIPE)
	    p = subprocess.Popen([commandline], shell=True, stdout=subprocess.PIPE)
	    output, err = p.communicate()
	    #print output
	    sig = 0
	    for line in output.split('\n'):
		#if line.find("Significance:")>0: 
		if "Significance:" in line.split():
			elts = line.split()
    			sig =  elts.pop()
	    masses = item.split(".")[0].split("_")
	    #print "The significance for ", masses[0], masses[1], " is ", sig
	    result = masses[0] + " " + masses[1] + " " + sig + "\n"
	    with open(ofilename, "a") as myfile:
	        myfile.write(result)


	    global counter
	    counter+=1
	    if (counter/10) == 0:
	      print counter, " jobs done"
            #print "task", item, "finished"

#
# run with limited queue

queue = Queue.Queue(0)

for i in range(WORKERS):
    Worker(queue).start() # start a worker


files =  [f for f in listdir(directory) if f.find(".tab")>0]
#print files

# output results
# remove the file
command = "rm " + ofilename
os.system(command)


#for item in range(njobs):
njobs = len(files)
indice = 0
for item in files:
    #print "push", item
    queue.put(item)
    if indice/100 == 0: 
      print indice/100, "%"
    indice+=1 

for i in range(WORKERS):
    queue.put(None) # add end-of-queue markers
