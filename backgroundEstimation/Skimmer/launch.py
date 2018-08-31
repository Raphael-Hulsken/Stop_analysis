import threading
import Queue

import time, random, os

WORKERS = 8

executable = "./a.out "
ifilename = "/opt/sbg/scratch1/cms/echabert/store/tmp/tmp/JetHT/PBSjob"
ofilename = "output_"
treename  = "babyTuple"
#selection = "lep1_pt>20&&pfmet>150"
selection = "lep1_pt>20"
njobs = 550

	#let number=$number+1 ;./a.out ${ifilename}_${number}_${number}.root ${treename} ${ofilename}_${number}.root ${selection} &

class Worker(threading.Thread):

    def __init__(self, queue):
        self.__queue = queue
        threading.Thread.__init__(self)

    def run(self):
        while 1:
            item = self.__queue.get()
            if item is None:
                break # reached end of queue

            # pretend we're doing something that takes 10-100 ms
            time.sleep(random.randint(10, 100) / 10.0)
	    
	    # here we launch our job
	    commandline = executable + " " + ifilename+"_"+str(item)+"_"+str(item)+".root" + " " + treename + " " + ofilename+str(item)+".root" + " " + selection
	    os.system(commandline)

            print "task", item, "finished"

#
# run with limited queue

queue = Queue.Queue(0)

for i in range(WORKERS):
    Worker(queue).start() # start a worker

for item in range(njobs):
    print "push", item
    queue.put(item)

for i in range(WORKERS):
    queue.put(None) # add end-of-queue markers
