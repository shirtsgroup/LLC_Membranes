#!/usr/bin/env python

with open('jobs_running_old', 'r') as f:

    old = []
    for line in f:
        old.append(line)

with open('jobs_running', 'r') as f:

    current = []
    for line in f:
        current.append(line)

old = [float(str.strip(old[i])) for i in range(len(old))]
current = [float(str.strip(current[i])) for i in range(len(current))]

restarts = []
for i in old:
    if i not in current:
        restarts.append(i)

if restarts:
    for i in restarts:
        print(int(i))
else:
    print(None)