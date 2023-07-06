# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 14:17:51 2022

@author: maste
"""

def on_generation(ga_instance):
    #global last_fitness
    if(ga_instance.generations_completed%100 == 0):
        print("Generation = {generation},\tFitness = {fitness}".format(generation=ga_instance.generations_completed,fitness=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]))
    
    #print("Generation = {generation},\tFitness = {fitness},\tChange = {change}".format(generation=ga_instance.generations_completed,
    #                                                                                    fitness=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1],
    #                                                                                    change=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1] - last_fitness_old))
    #print("Generation = {generation}".format(generation=ga_instance.generations_completed))
    #print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]))
    #print("Change     = {change}".format(change=ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1] - last_fitness))
    last_fitness = ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]
    return last_fitness