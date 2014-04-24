import multiprocessing
import sys
import time

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