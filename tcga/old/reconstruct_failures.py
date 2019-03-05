import os
import sys

ROOT = sys.argv[1]
failures = set()

for d in os.listdir(ROOT):
    patient = d
    d = os.path.join(ROOT, d)
    if os.path.isdir(d):
       look_for = os.path.join(d, 'optimal.graphml')
       if not os.path.isfile(look_for):
           failures.add(patient)

with open(os.path.join(ROOT, 'failed.txt'), 'a') as fp:
    for failure in failures:
        print(failure)
        fp.write(failure+'\n')
