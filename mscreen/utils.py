# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 20:17:51 2021

@author: o_o
"""
import multiprocessing

def count_processors(num_inputs, num_processors):
    """
    Checks processors available and returns a safe number of them to
    utilize.

    :param int num_inputs: The number of inputs.
    :param int num_processors: The number of desired processors.

    :returns: The number of processors to use.
    """
    # first, if num_processors<= 0, determine the number of processors to
    # use programatically
    if num_processors<= 0:
        num_processors = multiprocessing.cpu_count()

    # reduce the number of processors if too many have been specified
    if num_inputs < num_processors:
        num_processors = num_inputs

    return num_processors

def start_processes(inputs, num_processors):
    """
    Creates a queue of inputs and outputs
    """

    # Create queues
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()

    # Submit tasks
    for item in inputs:
        task_queue.put(item)

    # Start worker processes
    for i in range(num_processors):
        multiprocessing.Process(target = worker, args = (task_queue, done_queue)).start()

    # Get and print results
    results = []
    for i in range(len(inputs)):
        results.append(done_queue.get())

    # Tell child processes to stop
    for i in range(num_processors):
        task_queue.put('STOP')

    results.sort(key = lambda tup: tup[0])

    return  [item[1] for item in map(list, results)]
