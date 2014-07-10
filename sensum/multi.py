'''
.. module:: multi
   :platform: Unix, Windows
   :synopsis: This module includes functions related to the high-level classification of multi-spectral satellite images.

.. moduleauthor:: Mostapha Harb <mostapha.harb@eucentre.it>
.. moduleauthor:: Daniele De Vecchi <daniele.devecchi03@universitadipavia.it>
.. moduleauthor:: Daniel Aurelio Galeazzo <dgaleazzo@gmail.com>
   :organization: EUCENTRE Foundation / University of Pavia
'''
'''
---------------------------------------------------------------------------------
Created on April 28, 2014
Last modified on May 13, 2014

---------------------------------------------------------------------------------
Project: Framework to integrate Space-based and in-situ sENSing for dynamic 
         vUlnerability and recovery Monitoring (SENSUM)

Co-funded by the European Commission under FP7 (Seventh Framework Programme)
THEME [SPA.2012.1.1-04] Support to emergency response management
Grant agreement no: 312972

---------------------------------------------------------------------------------
License: This file is part of SensumTools.

    SensumTools is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SensumTools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SensumTools.  If not, see <http://www.gnu.org/licenses/>.
---------------------------------------------------------------------------------
'''

import multiprocessing
import sys
import time
import os,sys
'''
sys.path.append("C:\\OSGeo4W64\\apps\\Python27\\Lib\\site-packages")
sys.path.append("C:\\OSGeo4W64\\apps\\orfeotoolbox\\python")
os.environ["PATH"] = os.environ["PATH"] + "C:\\OSGeo4W64\\bin"
'''
class Multi(object):
    
    
    def __init__(self, num_consumers=0):
        # Establish communication queues
        self.tasks = multiprocessing.JoinableQueue()
        self.results = multiprocessing.Queue()
        # Start consumers
        if num_consumers:
            self.num_consumers = num_consumers
        else:
            self.num_consumers = multiprocessing.cpu_count() * 2
        print 'Creating %d consumers' % self.num_consumers
        consumers = [ Consumer(self.tasks, self.results)
                      for i in xrange(self.num_consumers) ]
        for w in consumers:
            w.start()

    def put(self,Task):
        self.tasks.put(Task)

    def kill(self):
        # Poison killer
        for i in xrange(self.num_consumers):
            self.tasks.put(None)
        # Wait for all of the tasks to finish
        self.tasks.join()

    def result(self):
        result = self.results.get()
        return result


class Consumer(multiprocessing.Process):
    
    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                print '%s: Exiting' % proc_name
                self.task_queue.task_done()
                break
            print '%s: %s' % (proc_name, next_task)
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return