# test
import sys
import random
prompts = [">","#","END"]
sys.stdout.write(prompts[random.randint(0,2)])
input = sys.stdin.read()
sys.stdout.write('Received: %s\n'%input)
